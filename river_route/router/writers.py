from __future__ import annotations

from concurrent.futures import Executor
from typing import TYPE_CHECKING

import netCDF4 as nc
import numpy as np
import pandas as pd
import zarr
from zarr.codecs import BloscCodec

from ..types import DatetimeArray, FloatArray, PathInput

if TYPE_CHECKING:
    from .Router import Router

__all__ = ['null_writer', 'netcdf_writer', 'zarr_writer', 'bitround']

ZARR_RIVERS_PER_CHUNK = 500
ZARR_CHUNKS_PER_SHARD: int | None = None
ZARR_KEEPBITS = 12
ZARR_COMPRESSOR = BloscCodec(cname='lz4', clevel=5, shuffle='bitshuffle')


def bitround(values: FloatArray, keepbits: int) -> FloatArray:
    """
    A copy of a float32 array with all but the leading ``keepbits`` mantissa bits zeroed, rounded to nearest even.
    The relative error is bounded by ``2 ** -(keepbits + 1)`` for every normal float. With 23 or more keepbits, every
    float32 mantissa bit, or an array that is not float32, the array itself is returned, not a copy.
    """
    if keepbits >= 23 or values.dtype != np.float32:
        return values
    shift = np.uint32(23 - keepbits)
    bits = values.view(np.uint32)
    tie = np.right_shift(bits, shift)  # the first array this allocates is reused for the whole calculation
    tie &= np.uint32(1)  # 1 when the kept part is odd, which is what rounds a tie to even
    tie += bits
    tie += np.uint32((1 << int(shift)) - 1) >> np.uint32(1)
    tie &= ~np.uint32((1 << int(shift)) - 1)
    return tie.view(np.float32)


def _check_discharge_shape(router: Router, dates: DatetimeArray, discharge_array: FloatArray, path: PathInput) -> None:
    """Check the (river, time) discharge against the original rivers and the dates it is written with."""
    n_rivers = router.network.original_river_ids.shape[0]
    if discharge_array.shape[0] != n_rivers:
        raise ValueError(f'Cannot write {path}: {discharge_array.shape[0]} rows of discharge for {n_rivers} rivers')
    if dates.shape[0] != discharge_array.shape[1]:
        raise ValueError(
            f'Cannot write {path}: {dates.shape[0]} dates for {discharge_array.shape[1]} columns of discharge'
        )
    return


def null_writer(*args, **kwargs) -> None:
    """
    Return without writing anything. Useful for benchmarking and auditing other processes. If you actually need to
    write to the null device, your custom writer may need to do something like this:

    with open(os.devnull, 'wb') as sink:
        discharge_array.tofile(sink)
    """
    return


def netcdf_writer(
    router: Router,
    dates: DatetimeArray,
    discharge_array: FloatArray,
    discharge_file: PathInput,
    runoff_file: PathInput = '',
    *,
    thread_pool: Executor | None = None,
    threads: int = 1,
) -> None:
    """
    Write routed discharge, Q, to an uncompressed netCDF file with dimensions (riverId, time).

    The variable is dimensioned as the routed array is laid out, so the whole array is written in one call with
    nothing transposed. Writing (time, riverId) instead costs about twice the kernel time on a large network,
    because each river's series then has to be transposed out in blocks.

    Args:
        router: the Router that routed the discharge
        dates: datetime array corresponding to the discharge columns
        discharge_array: routed discharge values, C-order with shape (river, time)
        discharge_file: path to write the discharge data to
        runoff_file: path to the runoff file the discharge was routed from, or '' for channel routing
        thread_pool: not used, since the whole array is written in one call
        threads: not used, since the whole array is written in one call
    """
    _check_discharge_shape(router, dates, discharge_array, discharge_file)
    n_rivers, n_steps = discharge_array.shape
    with nc.Dataset(str(discharge_file), mode='w', format='NETCDF4') as ds:
        ds.createDimension('riverId', size=n_rivers)
        ds.createDimension('time', size=n_steps)
        ds.runoff_file = str(runoff_file)
        time_var = ds.createVariable('time', 'f8', ('time',))
        time_var.units = f'seconds since {pd.Timestamp(dates[0]).strftime("%Y-%m-%d %H:%M:%S")}'
        time_var[:] = (dates - dates[0]).astype('timedelta64[s]').astype(np.int64)
        id_var = ds.createVariable('riverId', 'i8', 'riverId')
        id_var[:] = router.network.original_river_ids
        discharge_var = ds.createVariable('Q', 'f4', ('riverId', 'time'))
        discharge_var[:] = discharge_array  # one call: blocking the write only slows it down
        discharge_var.long_name = 'Discharge at catchment outlet'
        discharge_var.standard_name = 'discharge'
        discharge_var.aggregation_method = 'mean'
        discharge_var.units = 'm3 s-1'
    return


def zarr_writer(
    router: Router,
    dates: DatetimeArray,
    discharge_array: FloatArray,
    discharge_file: PathInput,
    runoff_file: PathInput = '',
    *,
    thread_pool: Executor | None = None,
    threads: int = 1,
) -> None:
    """
    Write routed discharge, Q, to a zarr store with dimensions (riverId, time), replacing anything already at the path.
    This is the default writer.

    Args:
        router: the Router that routed the discharge
        dates: datetime array corresponding to the discharge columns
        discharge_array: routed discharge values, C-order with shape (river, time)
        discharge_file: path of the zarr store to write
        runoff_file: path to the runoff file the discharge was routed from, or '' for channel routing
        thread_pool: optional pool the chunks are rounded and written on, used as given and never shut down here
        threads: the concurrency of zarr's writes, and of the chunk writes when ``thread_pool`` is given
    """
    _check_discharge_shape(router, dates, discharge_array, discharge_file)
    n_rivers, n_steps = discharge_array.shape
    # zarr sizes its worker pool once per process, so the per-operation concurrency limit is what follows threads
    with zarr.config.set({'async.concurrency': threads}):
        group = zarr.create_group(str(discharge_file), overwrite=True, attributes={'runoff_file': str(runoff_file)})
        discharge = group.create_array(
            'Q',
            shape=(n_rivers, n_steps),
            chunks=(ZARR_RIVERS_PER_CHUNK, -1),
            shards=(ZARR_RIVERS_PER_CHUNK * ZARR_CHUNKS_PER_SHARD, -1) if ZARR_CHUNKS_PER_SHARD else None,
            dtype=discharge_array.dtype,
            filters=None,
            compressors=ZARR_COMPRESSOR,
            config={'write_empty_chunks': True},  # skips comparing every chunk to the fill value before writing it
            dimension_names=('riverId', 'time'),
            attributes={
                'long_name': 'Discharge at catchment outlet',
                'standard_name': 'discharge',
                'aggregation_method': 'mean',
                'units': 'm3 s-1',
            },
        )
        if ZARR_KEEPBITS >= 23:
            discharge[:] = discharge_array
        else:  # round a chunk of rivers at a time, so the copy this makes stays the size of one chunk

            def write_chunk(r0: int) -> None:
                r1 = min(r0 + ZARR_RIVERS_PER_CHUNK, n_rivers)
                discharge[r0:r1] = bitround(discharge_array[r0:r1], ZARR_KEEPBITS)

            starts = range(0, n_rivers, ZARR_RIVERS_PER_CHUNK)
            if thread_pool is not None and threads > 1:
                # one write call per chunk would otherwise serialize what zarr does concurrently for a whole-array
                # assignment: 10.0 s against 1.4 s on a year of the Amazon
                list(thread_pool.map(write_chunk, starts))  # list() so a worker exception propagates
            else:
                for chunk_start in starts:
                    write_chunk(chunk_start)
        time_var = group.create_array(
            'time',
            shape=(n_steps,),
            dtype='int64',
            dimension_names=('time',),
            attributes={
                'units': f'seconds since {pd.Timestamp(dates[0]).strftime("%Y-%m-%d %H:%M:%S")}',
                'calendar': 'standard',
            },
        )
        time_var[:] = (dates - dates[0]).astype('timedelta64[s]').astype(np.int64)
        id_var = group.create_array('riverId', shape=(n_rivers,), dtype='int64', dimension_names=('riverId',))
        id_var[:] = router.network.original_river_ids
        zarr.consolidate_metadata(str(discharge_file))
    return
