from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING

import netCDF4 as nc
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import zarr
from zarr.codecs import BloscCodec

from ..types import DatetimeArray, FloatArray, PathInput
from . import _river_kernels as river_kernels
from ._river_kernels import to_time_major

if TYPE_CHECKING:
    from .Router import Router

__all__ = ['null_writer', 'netcdf_writer', 'zarr_writer', 'parquet_writer', 'to_time_major', 'bitround']

ZARR_RIVERS_PER_CHUNK = 500
ZARR_CHUNKS_PER_SHARD: int | None = None
ZARR_KEEPBITS = 12
ZARR_COMPRESSOR = BloscCodec(cname='lz4', clevel=5, shuffle='bitshuffle')
PARQUET_WRITE_OPTIONS = {'compression': 'none', 'use_dictionary': False, 'write_statistics': False}
NETCDF_RIVERS_PER_WRITE = 1024


def _decoded(discharge_array: FloatArray) -> FloatArray:
    """The discharge as the float dtype it holds: float32 passes through, a narrowed uint16 buffer is float16."""
    return river_kernels.decode_discharge(discharge_array)


def bitround(values: FloatArray, keepbits: int) -> FloatArray:
    """
    A copy of a float32 array with all but the leading ``keepbits`` mantissa bits zeroed, rounded to nearest even.
    The relative error is bounded by ``2 ** -(keepbits + 1)`` for every normal float.
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
    """Check the (river, time) discharge against the network and the dates it is written with."""
    if discharge_array.shape[0] != router.network.river_ids.shape[0]:
        raise ValueError(
            f'Cannot write {path}: {discharge_array.shape[0]} rows of discharge for '
            f'{router.network.river_ids.shape[0]} rivers'
        )
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
) -> None:
    """
    Write routed discharge to an uncompressed netCDF file with dimensions (river_id, time).

    The variable is dimensioned as the routed array is laid out, so the whole array is written in one call with
    nothing transposed. Writing (time, river_id) instead costs about twice the kernel time on a large network,
    because each river's series then has to be transposed out in blocks.

    The variable is always float32, because netCDF has no half precision type. A ``discharge_dtype='float16'`` run
    is widened a block of ``NETCDF_RIVERS_PER_WRITE`` rivers at a time, so only a bounded float32 block is allocated,
    but the file is the same size as a float32 run and the widening costs time. Use ``zarr_writer`` or
    ``parquet_writer`` to keep float16 on disk.

    Args:
        router: the Router that routed the discharge
        dates: datetime array corresponding to the discharge columns
        discharge_array: routed discharge values, C-order with shape (river, time)
        discharge_file: path to write the discharge data to
        runoff_file: path to the lateral inflow used to generate the discharge values, if applicable
    """
    _check_discharge_shape(router, dates, discharge_array, discharge_file)
    rid = router.configs.var_river_id
    n_rivers, n_steps = discharge_array.shape
    with nc.Dataset(str(discharge_file), mode='w', format='NETCDF4') as ds:
        ds.createDimension(rid, size=n_rivers)
        ds.createDimension('time', size=n_steps)
        ds.runoff_file = str(runoff_file)
        time_var = ds.createVariable('time', 'f8', ('time',))
        time_var.units = f'seconds since {pd.Timestamp(dates[0]).strftime("%Y-%m-%d %H:%M:%S")}'
        time_var[:] = (dates - dates[0]).astype('timedelta64[s]').astype(np.int64)
        id_var = ds.createVariable(rid, 'i4', rid)
        id_var[:] = router.network.river_ids
        flow_var = ds.createVariable(router.configs.var_discharge, 'f4', (rid, 'time'))
        if discharge_array.dtype == np.float32:
            flow_var[:] = discharge_array  # one call: nothing to widen, and blocking the write only slows it down
        else:  # netCDF has no half type, so a float16 buffer is widened a block of rivers at a time
            decoded = _decoded(discharge_array)
            for r0 in range(0, n_rivers, NETCDF_RIVERS_PER_WRITE):
                r1 = min(r0 + NETCDF_RIVERS_PER_WRITE, n_rivers)
                flow_var[r0:r1] = decoded[r0:r1].astype(np.float32)
        flow_var.long_name = 'Discharge at catchment outlet'
        flow_var.standard_name = 'discharge'
        flow_var.aggregation_method = 'mean'
        flow_var.units = 'm3 s-1'
    return


def zarr_writer(
    router: Router,
    dates: DatetimeArray,
    discharge_array: FloatArray,
    discharge_file: PathInput,
    runoff_file: PathInput = '',
) -> None:
    """
    Write routed discharge to a zarr store with dimensions (river_id, time), replacing anything already at the path.
    This is the default writer.

    The river dimension comes first so that each chunk is a contiguous span of the routed buffer, copied out without
    a transpose. Chunks hold ``ZARR_RIVERS_PER_CHUNK`` rivers and span every time step of the file.

    Values are rounded to ``ZARR_KEEPBITS`` mantissa bits and compressed with ``ZARR_COMPRESSOR``, which on a year of
    the Amazon is 2.71x smaller than the raw array and faster to write than storing it uncompressed. The rounding is
    done here in numpy rather than by a zarr filter, a block of rivers at a time so nothing the size of the whole
    array is allocated and the array handed in is never modified. A narrowed ``discharge_dtype`` is stored as it is,
    since float16 has fewer mantissa bits than the rounding would keep.

    Setting ``ZARR_CHUNKS_PER_SHARD`` packs that many chunks into each shard file instead of writing one file per
    chunk. The store opens with ``xarray.open_zarr``.

    Args:
        router: the Router that routed the discharge
        dates: datetime array corresponding to the discharge columns
        discharge_array: routed discharge values, C-order with shape (river, time)
        discharge_file: path of the zarr store to write
        runoff_file: path to the lateral inflow used to generate the discharge values, if applicable
    """
    _check_discharge_shape(router, dates, discharge_array, discharge_file)
    rid = router.configs.var_river_id
    n_rivers, n_steps = discharge_array.shape
    # zarr sizes its worker pool once per process, so the per-operation concurrency limit is what follows threads
    with zarr.config.set({'async.concurrency': router.threads}):
        group = zarr.create_group(str(discharge_file), overwrite=True, attributes={'runoff_file': str(runoff_file)})
        flow = group.create_array(
            router.configs.var_discharge,
            shape=(n_rivers, n_steps),
            chunks=(ZARR_RIVERS_PER_CHUNK, -1),
            shards=(ZARR_RIVERS_PER_CHUNK * ZARR_CHUNKS_PER_SHARD, -1) if ZARR_CHUNKS_PER_SHARD else None,
            dtype=_decoded(discharge_array).dtype,
            filters=None,
            compressors=ZARR_COMPRESSOR,
            config={'write_empty_chunks': True},  # skips comparing every chunk to the fill value before writing it
            dimension_names=(rid, 'time'),
            attributes={
                'long_name': 'Discharge at catchment outlet',
                'standard_name': 'discharge',
                'aggregation_method': 'mean',
                'units': 'm3 s-1',
            },
        )
        by_river = _decoded(discharge_array)  # (river, time), C-order as the kernels wrote it
        if by_river.dtype != np.float32 or ZARR_KEEPBITS >= 23:
            flow[:] = by_river
        else:  # round a chunk of rivers at a time, so the copy this makes stays the size of one chunk

            def write_chunk(r0: int) -> None:
                r1 = min(r0 + ZARR_RIVERS_PER_CHUNK, n_rivers)
                flow[r0:r1] = bitround(by_river[r0:r1], ZARR_KEEPBITS)

            starts = range(0, n_rivers, ZARR_RIVERS_PER_CHUNK)
            if router.threads > 1:
                # one write call per chunk would otherwise serialize what zarr does concurrently for a whole-array
                # assignment: 10.0 s against 1.4 s on a year of the Amazon
                with ThreadPoolExecutor(router.threads) as pool:
                    list(pool.map(write_chunk, starts))  # list() so a worker exception propagates
            else:
                for r0 in starts:
                    write_chunk(r0)
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
        id_var = group.create_array(rid, shape=(n_rivers,), dtype='int64', dimension_names=(rid,))
        id_var[:] = router.network.river_ids
        zarr.consolidate_metadata(str(discharge_file))
    return


def parquet_writer(
    router: Router,
    dates: DatetimeArray,
    discharge_array: FloatArray,
    discharge_file: PathInput,
    runoff_file: PathInput = '',
) -> None:
    """
    Write routed discharge to a parquet file with one row per river and one column per time step.

    The first column holds the river ids and every other column is named by its time step as
    ``YYYY-MM-DDTHH:MM:SS``. Parquet is columnar, so this is the one writer that needs the transpose of what the
    kernels produce: ``to_time_major`` turns the routed ``(river, time)`` buffer into a C-order ``(time, river)``
    one, and each time step column is then a row of that array, which pyarrow wraps without a copy.

    ``PARQUET_WRITE_OPTIONS`` is passed to ``pyarrow.parquet.write_table`` and sets the compression. pyarrow writes a
    parquet file on one thread, so the ``threads`` given to ``Router.route`` does not apply.

    Args:
        router: the Router that routed the discharge
        dates: datetime array corresponding to the discharge columns
        discharge_array: routed discharge values, C-order with shape (river, time)
        discharge_file: path of the parquet file to write
        runoff_file: path to the lateral inflow used to generate the discharge values, if applicable
    """
    _check_discharge_shape(router, dates, discharge_array, discharge_file)
    by_step = _decoded(to_time_major(discharge_array))  # parquet needs each step's rivers contiguous
    names = [router.configs.var_river_id, *pd.DatetimeIndex(dates).strftime('%Y-%m-%dT%H:%M:%S')]
    columns = [pa.array(router.network.river_ids), *(pa.array(row) for row in by_step)]
    metadata = {'runoff_file': str(runoff_file), 'variable': router.configs.var_discharge, 'units': 'm3 s-1'}
    table = pa.table(columns, names=names).replace_schema_metadata(metadata)
    pq.write_table(table, str(discharge_file), **PARQUET_WRITE_OPTIONS)
    return
