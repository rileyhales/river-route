"""
Tests for the gridded runoff path on a small synthetic grid and catchment set: building a weight table,
aggregating grid depths to catchment volumes, and routing straight from grid files. Needs no downloaded data.
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import geopandas as gpd
import netCDF4 as nc
import numpy as np
import pandas as pd
import pytest
import scipy.sparse
import xarray as xr
from shapely.geometry import box

import river_route as rr

N_RIVERS = 4
NX, NY, NT = 6, 5, 12
RUNOFF_DEPTH = 0.001  # metres per timestep
PULSE_STEPS = 3


def write_runoff_grid(path: Path, cumulative: bool = False) -> None:
    """A regular lon/lat grid holding a uniform runoff depth for the first PULSE_STEPS timesteps."""
    lons = np.linspace(0.25, 2.75, NX)
    lats = np.linspace(0.25, 2.25, NY)
    values = np.zeros((NT, NY, NX), dtype=np.float32)
    values[:PULSE_STEPS] = RUNOFF_DEPTH
    if cumulative:
        values = np.cumsum(values, axis=0)
    with nc.Dataset(str(path), 'w') as ds:
        ds.createDimension('time', NT)
        ds.createDimension('lon', NX)
        ds.createDimension('lat', NY)
        time_var = ds.createVariable('time', 'f8', ('time',))
        time_var.units = 'seconds since 2000-01-01 00:00:00'
        time_var[:] = np.arange(NT) * 3600
        ds.createVariable('lon', 'f8', ('lon',))[:] = lons
        ds.createVariable('lat', 'f8', ('lat',))[:] = lats
        runoff = ds.createVariable('ro', 'f4', ('time', 'lat', 'lon'))
        runoff[:] = values
        runoff.units = 'm'
    return


def write_catchments(path: Path) -> None:
    """Four rectangular catchments tiling the grid, chained headwater to outlet."""
    gpd.GeoDataFrame(
        {'river_id': np.arange(1, N_RIVERS + 1, dtype=np.int64)},
        geometry=[box(0, 0, 1.5, 1.25), box(1.5, 0, 3, 1.25), box(0, 1.25, 1.5, 2.5), box(1.5, 1.25, 3, 2.5)],
        crs=4326,
    ).to_parquet(path)
    return


def write_chain_params(path: Path) -> None:
    pd.DataFrame(
        {
            'river_id': np.arange(1, N_RIVERS + 1, dtype=np.int64),
            'next_river_id': np.array([2, 3, 4, -1], dtype=np.int64),
            'k': np.full(N_RIVERS, 3600.0),
            'x': np.full(N_RIVERS, 0.2),
        }
    ).to_parquet(path, index=False)
    return


@pytest.fixture
def grid_case(tmp_path: Path) -> dict:
    """A grid, catchments, params, and the weight table that connects them."""
    grid_file = tmp_path / 'runoff.nc'
    catchments_file = tmp_path / 'catchments.parquet'
    params_file = tmp_path / 'params.parquet'
    weights_file = tmp_path / 'weights.nc'
    write_runoff_grid(grid_file)
    write_catchments(catchments_file)
    write_chain_params(params_file)
    weights = rr.runoff.grid_weights(
        grid_file,
        catchments_file,
        var_x='lon',
        var_y='lat',
        save_weights_path=weights_file,
        routing_params_path=params_file,
    )
    return {
        'directory': tmp_path,
        'grid_file': grid_file,
        'catchments_file': catchments_file,
        'params_file': params_file,
        'weights_file': weights_file,
        'weights': weights,
        'total_area': float(weights.groupby('river_id')['area_sqm'].sum().sum()),
    }


def prepare(grid_case: dict, **kwargs) -> rr.Runoff:
    """A Runoff for the synthetic grid's variable names and weight table."""
    return rr.Runoff(grid_case['weights_file'], var_runoff='ro', var_x='lon', var_y='lat', **kwargs)


def test_grid_weights_proportions_sum_to_one(grid_case):
    weights = grid_case['weights']
    expected_columns = {'river_id', 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion'}
    assert expected_columns.issubset(weights.columns)
    sums = weights.groupby('river_id')['proportion'].sum()
    np.testing.assert_allclose(sums.to_numpy(), 1.0, rtol=1e-9)
    assert set(sums.index) == set(range(1, N_RIVERS + 1))


def test_grid_weights_pass_deep_validation(grid_case):
    """The weight table this writes must satisfy the checks route() runs on it."""
    rr.Configs(
        forcing='vlateral',
        params_file=str(grid_case['params_file']),
        grid_weights_file=str(grid_case['weights_file']),
        grid_runoff_files=[str(grid_case['grid_file'])],
        discharge_files=[str(grid_case['directory'] / 'q.nc')],
    ).deep_validate()


def test_to_dataset_conserves_volume(grid_case):
    """Aggregating a uniform depth over the catchments must give depth times total catchment area."""
    ds = prepare(grid_case, as_volumes=True).to_dataset(grid_case['grid_file'])
    expected = grid_case['total_area'] * RUNOFF_DEPTH * PULSE_STEPS
    assert float(ds['vlateral'].sum()) == pytest.approx(expected, rel=1e-6)
    assert ds['vlateral'].shape == (NT, N_RIVERS)
    assert ds['vlateral'].attrs['units'] == 'm3'


def test_to_dataset_returns_depths_by_default(grid_case):
    ds = prepare(grid_case).to_dataset(grid_case['grid_file'])
    assert ds['vlateral'].attrs['units'] == 'm'
    np.testing.assert_allclose(ds['vlateral'].values[:PULSE_STEPS], RUNOFF_DEPTH, rtol=1e-6)
    np.testing.assert_allclose(ds['vlateral'].values[PULSE_STEPS:], 0.0, atol=1e-12)


def test_to_dataset_river_order_matches_params(grid_case):
    """The router indexes vlateral by position, so the column order has to follow the params file."""
    ds = prepare(grid_case).to_dataset(grid_case['grid_file'])
    params = pd.read_parquet(grid_case['params_file'])
    np.testing.assert_array_equal(ds['river_id'].values, params['river_id'].to_numpy())


def test_cumulative_runoff_matches_incremental(grid_case, tmp_path):
    """Cumulative grids must de-accumulate to the same volumes as the incremental grid."""
    cumulative_file = tmp_path / 'runoff_cumulative.nc'
    write_runoff_grid(cumulative_file, cumulative=True)

    incremental = prepare(grid_case, as_volumes=True).to_dataset(grid_case['grid_file'])
    converted = prepare(grid_case, cumulative=True, as_volumes=True).to_dataset(cumulative_file)
    np.testing.assert_allclose(converted['vlateral'].values, incremental['vlateral'].values, rtol=1e-5, atol=1e-6)


def test_route_from_grid_files(grid_case):
    """Routing straight from grid files must conserve the water the grid delivered."""
    out = grid_case['directory'] / 'q_grid.nc'
    rr.Router(
        forcing='vlateral',
        params_file=str(grid_case['params_file']),
        grid_weights_file=str(grid_case['weights_file']),
        grid_runoff_files=[str(grid_case['grid_file'])],
        discharge_files=[str(out)],
        var_grid_runoff='ro',
        var_x='lon',
        var_y='lat',
        log=False,
        progress_bar=False,
    ).route()

    with xr.open_dataset(out) as ds:
        q = ds['Q'].values
        np.testing.assert_array_equal(ds['river_id'].values, np.arange(1, N_RIVERS + 1))
    assert q.shape == (NT, N_RIVERS)
    assert not np.isnan(q).any()
    assert np.all(q >= 0)
    # the outlet drains the whole grid, but 12 hours is not long enough for all of it to arrive
    delivered = q[:, -1].sum() * 3600
    assert 0 < delivered < grid_case['total_area'] * RUNOFF_DEPTH * PULSE_STEPS


def test_route_from_grid_matches_route_from_vlateral(grid_case):
    """Precomputing vlateral and routing it must equal routing the grid files directly."""
    vlateral_file = grid_case['directory'] / 'vlateral.nc'
    ds = prepare(grid_case, as_volumes=True).to_dataset(grid_case['grid_file'])
    ds.to_netcdf(vlateral_file)

    shared = dict(
        params_file=str(grid_case['params_file']),
        var_grid_runoff='ro',
        var_x='lon',
        var_y='lat',
        log=False,
        progress_bar=False,
    )
    from_grid = grid_case['directory'] / 'q_from_grid.nc'
    from_vlateral = grid_case['directory'] / 'q_from_vlateral.nc'
    rr.Router(
        forcing='vlateral',
        grid_weights_file=str(grid_case['weights_file']),
        grid_runoff_files=[str(grid_case['grid_file'])],
        discharge_files=[str(from_grid)],
        **shared,
    ).route()
    rr.Router(
        forcing='vlateral', vlateral_files=[str(vlateral_file)], discharge_files=[str(from_vlateral)], **shared
    ).route()

    with xr.open_dataset(from_grid) as a, xr.open_dataset(from_vlateral) as b:
        np.testing.assert_allclose(a['Q'].values, b['Q'].values, rtol=1e-5, atol=1e-6)


def test_unknown_runoff_units_raise(grid_case, tmp_path):
    bad_units = tmp_path / 'bad_units.nc'
    write_runoff_grid(bad_units)
    with nc.Dataset(str(bad_units), 'a') as ds:
        ds['ro'].units = 'furlongs'
    with pytest.raises(ValueError, match='Unknown units'):
        prepare(grid_case).to_dataset(bad_units)


def test_millimeter_runoff_is_converted(grid_case, tmp_path):
    """A grid in mm must produce depths 1000x smaller than the same numbers read as metres."""
    mm_file = tmp_path / 'runoff_mm.nc'
    write_runoff_grid(mm_file)
    with nc.Dataset(str(mm_file), 'a') as ds:
        ds['ro'].units = 'mm'
    in_mm = prepare(grid_case).to_dataset(mm_file)
    in_m = prepare(grid_case).to_dataset(grid_case['grid_file'])
    np.testing.assert_allclose(in_mm['vlateral'].values * 1000, in_m['vlateral'].values, rtol=1e-5)


def write_irregular_runoff_grid(path: Path) -> np.ndarray:
    """A grid whose time axis skips steps, which triggers the resample branch in Runoff.aggregate."""
    offsets = np.array([0, 3600, 7200, 14400, 21600], dtype='f8')  # 0, 1, 2, 4, 6 hours
    with nc.Dataset(str(path), 'w') as ds:
        ds.createDimension('time', offsets.shape[0])
        ds.createDimension('lon', NX)
        ds.createDimension('lat', NY)
        time_var = ds.createVariable('time', 'f8', ('time',))
        time_var.units = 'seconds since 2000-01-01 00:00:00'
        time_var[:] = offsets
        ds.createVariable('lon', 'f8', ('lon',))[:] = np.linspace(0.25, 2.75, NX)
        ds.createVariable('lat', 'f8', ('lat',))[:] = np.linspace(0.25, 2.25, NY)
        runoff = ds.createVariable('ro', 'f4', ('time', 'lat', 'lon'))
        runoff[:] = np.full((offsets.shape[0], NY, NX), RUNOFF_DEPTH, dtype=np.float32)
        runoff.units = 'm'
    return offsets


def test_irregular_timesteps_are_resampled(grid_case, tmp_path):
    """Irregular input is resampled onto the first timestep. This path converts to numpy through pandas,
    which returns a read-only array, and the volume conversion then writes into it in place."""
    irregular = tmp_path / 'runoff_irregular.nc'
    write_irregular_runoff_grid(irregular)

    ds = prepare(grid_case, as_volumes=True).to_dataset(irregular)
    times = pd.to_datetime(ds['time'].values)
    assert ds['vlateral'].shape == (7, N_RIVERS)  # 0 through 6 hours, filled in hourly
    np.testing.assert_array_equal(np.diff(times).astype('timedelta64[s]').astype(int), 3600)
    assert not np.isnan(ds['vlateral'].values).any()
    assert float(ds['vlateral'].sum()) > 0


def hourly_times(n_steps: int) -> np.ndarray:
    return np.datetime64('2000-01-01T00:00', 's') + np.arange(n_steps) * np.timedelta64(3600, 's')


@pytest.mark.parametrize('conversion_factor', [1, 0.001])
@pytest.mark.parametrize('as_volumes', [False, True])
@pytest.mark.parametrize('force_positive', [False, True])
@pytest.mark.parametrize('cumulative', [False, True])
def test_one_pass_aggregation_matches_array_passes(
    grid_case, cumulative, force_positive, as_volumes, conversion_factor
):
    """The one-pass kernel must equal the same conversion done as separate whole-array passes, NaN included."""
    preparer = prepare(grid_case, cumulative=cumulative, force_positive_runoff=force_positive, as_volumes=as_volumes)
    n_cells = preparer.x_index.shape[0]
    runoff = np.random.default_rng(42).normal(0.0, RUNOFF_DEPTH, size=(NT, n_cells)).astype(np.float32)
    runoff[4, 0] = np.nan

    matrix = scipy.sparse.csr_matrix(
        (preparer.proportion * conversion_factor, preparer.cell, preparer.indptr), shape=(N_RIVERS, n_cells)
    )
    expected = np.asarray(matrix @ runoff.T).T.copy()
    if cumulative:
        expected[1:] = np.diff(expected, axis=0)
    if force_positive:
        np.clip(expected, 0, None, out=expected)
    expected[np.isnan(expected)] = 0.0
    if as_volumes:
        expected *= preparer.catchment_area[np.newaxis, :]

    times = hourly_times(NT)
    vlateral, vlateral_times = preparer.aggregate(runoff, times, conversion_factor)
    assert vlateral.flags.c_contiguous
    np.testing.assert_array_equal(vlateral_times, times)
    np.testing.assert_allclose(vlateral, expected, rtol=1e-5, atol=1e-6 * float(np.abs(expected).max()))


def test_aggregation_writes_into_a_reused_buffer(grid_case):
    """The router hands the kernels a view of one reused buffer, so aggregation must write into it, not copy."""
    preparer = prepare(grid_case)
    runoff = np.full((NT, preparer.x_index.shape[0]), RUNOFF_DEPTH, dtype=np.float32)
    buffer = np.full((NT + 5, N_RIVERS), -1.0, dtype=np.float32)
    vlateral, _ = preparer.aggregate(runoff, hourly_times(NT), out=buffer)
    assert vlateral.shape == (NT, N_RIVERS)
    assert np.shares_memory(vlateral, buffer)
    assert vlateral.flags.c_contiguous
    np.testing.assert_allclose(vlateral, RUNOFF_DEPTH, rtol=1e-6)
    with pytest.raises(ValueError, match='out must be'):
        preparer.aggregate(runoff, hourly_times(NT), out=buffer[: NT - 1])


@pytest.mark.parametrize('cumulative', [False, True])
def test_threaded_aggregation_matches_single_threaded(grid_case, cumulative):
    """Aggregating river ranges concurrently must give exactly what one range over every river gives."""
    shared = dict(cumulative=cumulative, force_positive_runoff=True, as_volumes=True)
    single_preparer = prepare(grid_case, **shared)
    runoff = np.random.default_rng(7).normal(0.0, RUNOFF_DEPTH, size=(NT, single_preparer.x_index.shape[0]))
    runoff = runoff.astype(np.float32)
    single, _ = single_preparer.aggregate(runoff, hourly_times(NT))
    with ThreadPoolExecutor(3) as pool:
        threaded, _ = prepare(grid_case, thread_pool=pool, threads=3, **shared).aggregate(runoff, hourly_times(NT))
    np.testing.assert_array_equal(threaded, single)


def test_route_with_thread_pool_matches_single_threaded(grid_case):
    """A thread_pool opened in a with block routes the same discharge as a single thread."""
    shared = dict(
        forcing='vlateral',
        params_file=str(grid_case['params_file']),
        grid_weights_file=str(grid_case['weights_file']),
        grid_runoff_files=[str(grid_case['grid_file'])],
        var_grid_runoff='ro',
        var_x='lon',
        var_y='lat',
        log=False,
        progress_bar=False,
    )
    single = grid_case['directory'] / 'q_single.nc'
    threaded = grid_case['directory'] / 'q_threaded.nc'
    rr.Router(discharge_files=[str(single)], **shared).route()
    with ThreadPoolExecutor(2) as pool:
        rr.Router(discharge_files=[str(threaded)], threads=2, **shared).route(thread_pool=pool)
    with xr.open_dataset(single) as a, xr.open_dataset(threaded) as b:
        np.testing.assert_array_equal(a['Q'].values, b['Q'].values)


@pytest.mark.parametrize('n_ranges', [1, 2, 3, 7, 50])
def test_river_ranges_cover_every_river_once(n_ranges):
    indptr = np.array([0, 5, 6, 6, 20, 21, 22, 40], dtype=np.int32)
    bounds = rr.Runoff._river_ranges(indptr, n_ranges)
    assert bounds[0] == 0
    assert bounds[-1] == indptr.shape[0] - 1
    assert np.all(np.diff(bounds) > 0)
    assert bounds.shape[0] - 1 <= n_ranges


def test_router_grid_forcing_needs_no_copy(grid_case):
    """The kernels need C-order float32 vlateral; the grid path must produce it so route() never copies it."""
    router = rr.Router(
        forcing='vlateral',
        params_file=str(grid_case['params_file']),
        grid_weights_file=str(grid_case['weights_file']),
        grid_runoff_files=[str(grid_case['grid_file'])],
        discharge_files=[str(grid_case['directory'] / 'q.nc')],
        var_grid_runoff='ro',
        var_x='lon',
        var_y='lat',
        log=False,
        progress_bar=False,
    )
    router._set_vectors_from_params()
    _, vlateral, _, _ = next(router._vlateral_generator())
    assert vlateral.shape == (NT, N_RIVERS)
    assert vlateral.dtype == np.float32
    assert vlateral.flags.c_contiguous
