"""
Tests for the premade discharge writers in river_route.writers, run on the synthetic network from conftest.
"""

import os
import warnings

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
import xarray as xr
import zarr
from conftest import SyntheticNetwork, build_network

import river_route as rr


def route_with_writer(network: SyntheticNetwork, out_name: str, writer=None, **kwargs) -> str:
    out = network.path(out_name)
    router = rr.Router(
        forcing='vlateral',
        params_file=str(network.params_file),
        vlateral_files=[str(network.vlateral_file)],
        channel_state_init_file=str(network.state_file),
        dt_routing=network.dt_runoff,
        discharge_files=[out],
        log=False,
        progress_bar=False,
        **kwargs,
    )
    if writer is not None:
        router.set_write_discharges(writer)
    router.route()
    return out


def test_netcdf_writer_is_the_default(network: SyntheticNetwork):
    router = rr.Router(params_file=str(network.params_file), discharge_files=[network.path('q.nc')])
    assert router._write_discharges is rr.writers.netcdf_writer


def test_zarr_writer_matches_netcdf_exactly(network: SyntheticNetwork):
    """The zarr store holds exactly the netCDF discharge, in the same (time, river_id) layout and coordinates."""
    netcdf_out = route_with_writer(network, 'q.nc')
    zarr_out = route_with_writer(network, 'q.zarr', writer=rr.writers.zarr_writer)
    with xr.open_dataset(netcdf_out) as ds_nc, xr.open_zarr(zarr_out) as ds_zarr:
        assert ds_zarr['Q'].dims == ('time', 'river_id')
        assert ds_zarr['Q'].dtype == np.float32
        assert ds_nc['Q'].values.any()
        np.testing.assert_array_equal(ds_zarr['Q'].values, ds_nc['Q'].values)
        np.testing.assert_array_equal(ds_zarr['river_id'].values, ds_nc['river_id'].values)
        np.testing.assert_array_equal(ds_zarr['time'].values, ds_nc['time'].values)
        assert ds_zarr.attrs['runoff_file'] == ds_nc.attrs['runoff_file']
        assert ds_zarr['Q'].attrs['units'] == 'm3 s-1'


def test_zarr_writer_layout_is_uncompressed_with_whole_time_chunks(network: SyntheticNetwork):
    out = route_with_writer(network, 'q.zarr', writer=rr.writers.zarr_writer)
    array = zarr.open_group(out, mode='r')['Q']
    assert array.chunks == (array.shape[0], rr.writers.ZARR_RIVERS_PER_CHUNK)
    assert array.filters == ()
    assert array.compressors == ()


def test_zarr_writer_splits_rivers_into_chunks(tmp_path):
    """A network wider than one chunk is split only on the river dimension and still reads back exactly."""
    n_rivers = rr.writers.ZARR_RIVERS_PER_CHUNK * 2 + 7
    router = rr.Router(
        params_file=str(build_network(tmp_path / 'net').params_file), discharge_files=[str(tmp_path / 'q.nc')]
    )
    router.river_ids = np.arange(1, n_rivers + 1, dtype=np.int64)
    dates = pd.date_range('2000-01-01', periods=6, freq='h').to_numpy()
    discharge = np.random.default_rng(0).random((6, n_rivers), dtype=np.float32)
    discharge[:, : rr.writers.ZARR_RIVERS_PER_CHUNK] = 0  # an all-zero chunk is still written, not skipped as empty
    out = str(tmp_path / 'wide.zarr')

    rr.writers.zarr_writer(router, dates, discharge, out)

    array = zarr.open_group(out, mode='r')['Q']
    assert array.nchunks == 3
    assert array.nchunks_initialized == 3
    with xr.open_zarr(out) as ds:
        np.testing.assert_array_equal(ds['Q'].values, discharge)
        np.testing.assert_array_equal(ds['river_id'].values, router.river_ids)
        np.testing.assert_array_equal(ds['time'].values, dates)


def test_zarr_writer_packs_chunks_into_shards(tmp_path, monkeypatch):
    """With sharding on, chunks are packed into shard files and a partial last shard still reads back exactly."""
    monkeypatch.setattr(rr.writers, 'ZARR_CHUNKS_PER_SHARD', 2)
    n_rivers = rr.writers.ZARR_RIVERS_PER_CHUNK * 3 + 7
    router = rr.Router(
        params_file=str(build_network(tmp_path / 'net').params_file), discharge_files=[str(tmp_path / 'q.nc')]
    )
    router.river_ids = np.arange(1, n_rivers + 1, dtype=np.int64)
    dates = pd.date_range('2000-01-01', periods=6, freq='h').to_numpy()
    discharge = np.random.default_rng(0).random((6, n_rivers), dtype=np.float32)
    out = str(tmp_path / 'sharded.zarr')

    rr.writers.zarr_writer(router, dates, discharge, out)

    array = zarr.open_group(out, mode='r')['Q']
    assert array.chunks == (6, rr.writers.ZARR_RIVERS_PER_CHUNK)
    assert array.shards == (6, rr.writers.ZARR_RIVERS_PER_CHUNK * 2)
    assert sum(len(files) for _, _, files in os.walk(os.path.join(out, 'Q', 'c'))) == 2
    with xr.open_zarr(out) as ds:
        np.testing.assert_array_equal(ds['Q'].values, discharge)
        np.testing.assert_array_equal(ds['river_id'].values, router.river_ids)


def test_zarr_writer_threads_do_not_change_output(network: SyntheticNetwork):
    single = route_with_writer(network, 'q_single.zarr', writer=rr.writers.zarr_writer)
    threaded = route_with_writer(network, 'q_threaded.zarr', writer=rr.writers.zarr_writer, threads=4)
    with xr.open_zarr(single) as ds_single, xr.open_zarr(threaded) as ds_threaded:
        np.testing.assert_array_equal(ds_single['Q'].values, ds_threaded['Q'].values)


def test_zarr_writer_replaces_an_existing_store(network: SyntheticNetwork):
    out = route_with_writer(network, 'q.zarr', writer=rr.writers.zarr_writer)
    os.makedirs(os.path.join(out, 'stale'))
    route_with_writer(network, 'q.zarr', writer=rr.writers.zarr_writer)
    assert not os.path.exists(os.path.join(out, 'stale'))
    with xr.open_zarr(out) as ds:
        assert ds['Q'].shape == (ds.sizes['time'], network.n_rivers)


def test_zarr_store_opens_without_warnings(network: SyntheticNetwork):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        out = route_with_writer(network, 'q.zarr', writer=rr.writers.zarr_writer)
        with xr.open_zarr(out) as ds:
            assert ds['Q'].values.shape == (ds.sizes['time'], network.n_rivers)
    assert not [str(w.message) for w in caught if 'zarr' in str(w.message).lower()]


@pytest.mark.parametrize('writer', [rr.writers.netcdf_writer, rr.writers.zarr_writer, rr.writers.parquet_writer])
def test_writers_reject_mismatched_shapes(network: SyntheticNetwork, writer):
    router = rr.Router(params_file=str(network.params_file), discharge_files=[network.path('q.nc')])
    router.river_ids = np.arange(1, network.n_rivers + 1, dtype=np.int64)
    dates = pd.date_range('2000-01-01', periods=4, freq='h').to_numpy()
    with pytest.raises(ValueError, match='dates for'):
        writer(router, dates[:3], np.zeros((4, network.n_rivers), np.float32), network.path('a.out'))
    with pytest.raises(ValueError, match='columns of discharge'):
        writer(router, dates, np.zeros((4, network.n_rivers + 1), np.float32), network.path('b.out'))


def read_parquet_discharge(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (discharge as (time, river), river ids, datetime64 times) from a parquet_writer file."""
    df = pd.read_parquet(path).set_index('river_id')
    return df.to_numpy().T, df.index.to_numpy(), pd.to_datetime(df.columns).to_numpy()


def test_parquet_writer_matches_netcdf_exactly(network: SyntheticNetwork):
    """The parquet file holds exactly the netCDF discharge with one row per river and one column per time step."""
    netcdf_out = route_with_writer(network, 'q.nc')
    parquet_out = route_with_writer(network, 'q.parquet', writer=rr.writers.parquet_writer)
    discharge, river_ids, times = read_parquet_discharge(parquet_out)
    with xr.open_dataset(netcdf_out) as ds_nc:
        assert ds_nc['Q'].values.any()
        assert discharge.dtype == np.float32
        np.testing.assert_array_equal(discharge, ds_nc['Q'].values)
        np.testing.assert_array_equal(river_ids, ds_nc['river_id'].values)
        np.testing.assert_array_equal(times, ds_nc['time'].values)
        metadata = pq.read_schema(parquet_out).metadata
        assert metadata[b'runoff_file'].decode() == ds_nc.attrs['runoff_file']
        assert metadata[b'units'] == b'm3 s-1'


def test_parquet_writer_uses_the_configured_write_options(network: SyntheticNetwork, monkeypatch):
    monkeypatch.setattr(rr.writers, 'PARQUET_WRITE_OPTIONS', {'compression': 'zstd', 'use_dictionary': False})
    out = route_with_writer(network, 'q.parquet', writer=rr.writers.parquet_writer)
    column = pq.ParquetFile(out).metadata.row_group(0).column(1)
    assert column.compression == 'ZSTD'
    assert 'RLE_DICTIONARY' not in column.encodings


def test_parquet_writer_replaces_an_existing_file(network: SyntheticNetwork):
    out = network.path('q.parquet')
    with open(out, 'w') as f:
        f.write('not parquet')
    route_with_writer(network, 'q.parquet', writer=rr.writers.parquet_writer)
    discharge, river_ids, _ = read_parquet_discharge(out)
    assert discharge.shape[1] == river_ids.shape[0] == network.n_rivers
