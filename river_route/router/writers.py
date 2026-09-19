from __future__ import annotations

from typing import TYPE_CHECKING

import netCDF4 as nc
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import zarr

from ..types import DatetimeArray, FloatArray, PathInput

if TYPE_CHECKING:
    from .Router import Router

__all__ = ['null_writer', 'netcdf_writer', 'zarr_writer', 'parquet_writer']

ZARR_RIVERS_PER_CHUNK = 500
ZARR_CHUNKS_PER_SHARD: int | None = None
PARQUET_WRITE_OPTIONS = {'compression': 'none', 'use_dictionary': False, 'write_statistics': False}


def _check_discharge_shape(router: Router, dates: DatetimeArray, discharge_array: FloatArray, path: PathInput) -> None:
    if dates.shape[0] != discharge_array.shape[0]:
        raise ValueError(
            f'Cannot write {path}: {dates.shape[0]} dates for {discharge_array.shape[0]} rows of discharge'
        )
    if discharge_array.shape[1] != router.network.river_ids.shape[0]:
        raise ValueError(
            f'Cannot write {path}: {discharge_array.shape[1]} columns of discharge for '
            f'{router.network.river_ids.shape[0]} rivers'
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
    Write routed discharge to an uncompressed netCDF file with dimensions (time, river_id).

    Args:
        router: the Router that routed the discharge
        dates: datetime array corresponding to the discharge rows
        discharge_array: routed discharge values with shape (time, river)
        discharge_file: path to write the discharge data to
        runoff_file: path to the lateral inflow used to generate the discharge values, if applicable
    """
    _check_discharge_shape(router, dates, discharge_array, discharge_file)
    with nc.Dataset(str(discharge_file), mode='w', format='NETCDF4') as ds:
        ds.createDimension('time', size=discharge_array.shape[0])
        ds.createDimension(router.configs.var_river_id, size=discharge_array.shape[1])
        ds.runoff_file = str(runoff_file)
        time_var = ds.createVariable('time', 'f8', ('time',))
        time_var.units = f'seconds since {pd.Timestamp(dates[0]).strftime("%Y-%m-%d %H:%M:%S")}'
        time_var[:] = (dates - dates[0]).astype('timedelta64[s]').astype(np.int64)
        id_var = ds.createVariable(router.configs.var_river_id, 'i4', router.configs.var_river_id)
        id_var[:] = router.network.river_ids
        flow_var = ds.createVariable(router.configs.var_discharge, 'f4', ('time', router.configs.var_river_id))
        flow_var[:] = discharge_array
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
    Write routed discharge to an uncompressed zarr store with dimensions (time, river_id), replacing anything already
    at the path. Built to write as fast as possible.

    The array is stored exactly as routed, with no filters or compression, in chunks of ``ZARR_RIVERS_PER_CHUNK``
    rivers that each span every time step of the file. Chunks are written one at a time with ``threads=1`` and up to
    the ``threads`` given to ``Router.route`` at once otherwise. Setting ``ZARR_CHUNKS_PER_SHARD`` packs that many
    chunks into each shard file instead of writing one file per chunk. The store opens with ``xarray.open_zarr``.

    Args:
        router: the Router that routed the discharge
        dates: datetime array corresponding to the discharge rows
        discharge_array: routed discharge values with shape (time, river)
        discharge_file: path of the zarr store to write
        runoff_file: path to the lateral inflow used to generate the discharge values, if applicable
    """
    _check_discharge_shape(router, dates, discharge_array, discharge_file)
    rid = router.configs.var_river_id
    n_steps, n_rivers = discharge_array.shape
    # zarr sizes its worker pool once per process, so the per-operation concurrency limit is what follows threads
    with zarr.config.set({'async.concurrency': router.threads}):
        group = zarr.create_group(str(discharge_file), overwrite=True, attributes={'runoff_file': str(runoff_file)})
        flow = group.create_array(
            router.configs.var_discharge,
            shape=(n_steps, n_rivers),
            chunks=(-1, ZARR_RIVERS_PER_CHUNK),
            shards=(-1, ZARR_RIVERS_PER_CHUNK * ZARR_CHUNKS_PER_SHARD) if ZARR_CHUNKS_PER_SHARD else None,
            dtype='float32',
            filters=None,
            compressors=None,
            config={'write_empty_chunks': True},  # skips comparing every chunk to the fill value before writing it
            dimension_names=('time', rid),
            attributes={
                'long_name': 'Discharge at catchment outlet',
                'standard_name': 'discharge',
                'aggregation_method': 'mean',
                'units': 'm3 s-1',
            },
        )
        flow[:] = discharge_array
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
    ``YYYY-MM-DDTHH:MM:SS``. Each time step column is a row of the routed array, which pyarrow wraps without a
    copy. ``PARQUET_WRITE_OPTIONS`` is passed to ``pyarrow.parquet.write_table`` and sets the compression. pyarrow
    writes a parquet file on one thread, so the ``threads`` given to ``Router.route`` does not apply.

    Args:
        router: the Router that routed the discharge
        dates: datetime array corresponding to the discharge rows
        discharge_array: routed discharge values with shape (time, river)
        discharge_file: path of the parquet file to write
        runoff_file: path to the lateral inflow used to generate the discharge values, if applicable
    """
    _check_discharge_shape(router, dates, discharge_array, discharge_file)
    discharge_array = np.ascontiguousarray(discharge_array, dtype=np.float32)
    names = [router.configs.var_river_id, *pd.DatetimeIndex(dates).strftime('%Y-%m-%dT%H:%M:%S')]
    columns = [pa.array(router.network.river_ids), *(pa.array(row) for row in discharge_array)]
    metadata = {'runoff_file': str(runoff_file), 'variable': router.configs.var_discharge, 'units': 'm3 s-1'}
    table = pa.table(columns, names=names).replace_schema_metadata(metadata)
    pq.write_table(table, str(discharge_file), **PARQUET_WRITE_OPTIONS)
    return
