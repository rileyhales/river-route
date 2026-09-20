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
from river_route.router import writers

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


def prepare(grid_case: dict, **kwargs) -> rr.RunoffGaussianGrid:
    """A RunoffGaussianGrid for the synthetic grid's variable names and weight table.

    RunoffGaussianGrid takes a Configs and nothing else, so the per-test options are folded into the Configs here.
    """
    configs = rr.Configs(
        grid_weights_file=grid_case['weights_file'], var_grid_runoff='ro', var_x='lon', var_y='lat', **kwargs
    )
    return rr.RunoffGaussianGrid.from_configs(configs)


def test_runoff_builds_the_same_thing_directly_and_from_configs(grid_case):
    """from_configs only reads the values off a Configs; the direct constructor takes those same values."""
    configs = rr.Configs(
        grid_weights_file=grid_case['weights_file'], var_grid_runoff='ro', var_x='lon', var_y='lat', as_volumes=True
    )
    direct = rr.RunoffGaussianGrid(
        grid_case['weights_file'], var_grid_runoff='ro', var_x='lon', var_y='lat', as_volumes=True
    )
    from_configs = rr.RunoffGaussianGrid.from_configs(configs)
    np.testing.assert_array_equal(direct.river_ids, from_configs.river_ids)
    np.testing.assert_array_equal(direct.proportion, from_configs.proportion)
    assert (direct.var_runoff, direct.var_x, direct.var_y, direct.as_volumes) == (
        from_configs.var_runoff,
        from_configs.var_x,
        from_configs.var_y,
        from_configs.as_volumes,
    )
    with pytest.raises(TypeError, match='from_configs takes a Configs'):
        rr.RunoffGaussianGrid.from_configs(grid_case['weights_file'])
    with pytest.raises(ValueError, match='grid_weights_file is required'):
        rr.RunoffGaussianGrid(None)


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
    converted = prepare(grid_case, grid_accumulation_type='cumulative', as_volumes=True).to_dataset(cumulative_file)
    np.testing.assert_allclose(converted['vlateral'].values, incremental['vlateral'].values, rtol=1e-5, atol=1e-6)


def test_route_from_grid_files(grid_case):
    """Routing straight from grid files must conserve the water the grid delivered."""
    out = grid_case['directory'] / 'q_grid.nc'
    rr.Router(
        rr.Configs(
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
        )
    ).set_discharge_writer(writers.netcdf_writer).route()

    with xr.open_dataset(out) as ds:
        q = ds['Q'].transpose('time', 'river_id').values
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
        rr.Configs(
            forcing='vlateral',
            grid_weights_file=str(grid_case['weights_file']),
            grid_runoff_files=[str(grid_case['grid_file'])],
            discharge_files=[str(from_grid)],
            **shared,
        )
    ).set_discharge_writer(writers.netcdf_writer).route()
    rr.Router(
        rr.Configs(
            forcing='vlateral', vlateral_files=[str(vlateral_file)], discharge_files=[str(from_vlateral)], **shared
        )
    ).set_discharge_writer(writers.netcdf_writer).route()

    with xr.open_dataset(from_grid) as a, xr.open_dataset(from_vlateral) as b:
        np.testing.assert_allclose(
            a['Q'].transpose('time', 'river_id').values,
            b['Q'].transpose('time', 'river_id').values,
            rtol=1e-5,
            atol=1e-6,
        )


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
    """A grid whose time axis skips steps, which triggers the resample branch in RunoffGaussianGrid.aggregate."""
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
    """The one-pass kernel must equal the same conversion done as separate whole-array passes, NaN cells as zero."""
    preparer = prepare(
        grid_case,
        grid_accumulation_type='cumulative' if cumulative else 'incremental',
        force_positive_runoff=force_positive,
        as_volumes=as_volumes,
    )
    n_cells = preparer.x_index.shape[0]
    runoff = np.random.default_rng(42).normal(0.0, RUNOFF_DEPTH, size=(NT, n_cells)).astype(np.float32)
    runoff[4, 0] = np.nan

    matrix = scipy.sparse.csr_matrix(
        (preparer.proportion * conversion_factor, preparer.cell, preparer.indptr), shape=(N_RIVERS, n_cells)
    )
    # a NaN cell contributes nothing, while the other cells of its catchments still count
    expected = np.asarray(matrix @ np.nan_to_num(runoff, nan=0.0).T).T.copy()
    if cumulative:
        expected[1:] = np.diff(expected, axis=0)
    if force_positive:
        np.clip(expected, 0, None, out=expected)
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
    shared = dict(
        grid_accumulation_type='cumulative' if cumulative else 'incremental',
        force_positive_runoff=True,
        as_volumes=True,
    )
    single_preparer = prepare(grid_case, **shared)
    runoff = np.random.default_rng(7).normal(0.0, RUNOFF_DEPTH, size=(NT, single_preparer.x_index.shape[0]))
    runoff = runoff.astype(np.float32)
    single, _ = single_preparer.aggregate(runoff, hourly_times(NT))
    with ThreadPoolExecutor(3) as pool:
        threaded, _ = prepare(grid_case, **shared).aggregate(runoff, hourly_times(NT), thread_pool=pool, threads=3)
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
    rr.Router(rr.Configs(discharge_files=[str(single)], **shared)).set_discharge_writer(writers.netcdf_writer).route()
    with ThreadPoolExecutor(2) as pool:
        rr.Router(rr.Configs(discharge_files=[str(threaded)], **shared)).set_discharge_writer(
            writers.netcdf_writer
        ).route(thread_pool=pool, threads=2)
    with xr.open_dataset(single) as a, xr.open_dataset(threaded) as b:
        np.testing.assert_array_equal(
            a['Q'].transpose('time', 'river_id').values, b['Q'].transpose('time', 'river_id').values
        )


@pytest.mark.parametrize('n_ranges', [1, 2, 3, 7, 50])
def test_river_ranges_cover_every_river_once(n_ranges):
    indptr = np.array([0, 5, 6, 6, 20, 21, 22, 40], dtype=np.int32)
    bounds = rr.RunoffGaussianGrid._river_ranges(indptr, n_ranges)
    assert bounds[0] == 0
    assert bounds[-1] == indptr.shape[0] - 1
    assert np.all(np.diff(bounds) > 0)
    assert bounds.shape[0] - 1 <= n_ranges


def test_runoff_reader_needs_no_copy(grid_case):
    """The kernels need C-order float32 vlateral; the grid path must produce it so route() never copies it."""
    configs = rr.Configs(
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
    runoff = rr.RunoffGaussianGrid.from_configs(configs)
    reader = runoff.reader(configs.grid_runoff_files)
    _, vlateral, _ = next(reader)
    assert vlateral.shape == (NT, N_RIVERS)
    assert vlateral.dtype == np.float32
    assert vlateral.flags.c_contiguous


def test_to_netcdf_is_read_by_runoff_vlateral(grid_case):
    """A file written by to_netcdf must read back through RunoffVlateral as the same dates and volumes."""
    runoff = prepare(grid_case)
    dates, vlateral, _ = next(runoff.reader([grid_case['grid_file']]))
    path = grid_case['directory'] / 'vlateral_written.nc'
    runoff.to_netcdf(path, dates, vlateral, runoff.river_ids)
    read_dates, read_vlateral, _ = next(rr.RunoffVlateral().reader([path]))
    np.testing.assert_array_equal(read_dates, dates)
    np.testing.assert_array_equal(read_vlateral, vlateral)


def test_runoff_is_abstract():
    with pytest.raises(TypeError):
        rr.runoff.Runoff()


def grid_configs(grid_case: dict, out_name: str) -> rr.Configs:
    return rr.Configs(
        forcing='vlateral',
        params_file=str(grid_case['params_file']),
        grid_weights_file=str(grid_case['weights_file']),
        grid_runoff_files=[str(grid_case['grid_file'])],
        discharge_files=[str(grid_case['directory'] / out_name)],
        var_grid_runoff='ro',
        var_x='lon',
        var_y='lat',
        log=False,
        progress_bar=False,
    )


def test_router_builds_and_reuses_one_runoff(grid_case):
    """A Router given no RunoffGaussianGrid has none until it routes, then keeps the one the reader built."""
    router = rr.Router(grid_configs(grid_case, 'q_lazy.nc'))
    assert router.runoff is None
    router.set_discharge_writer(writers.netcdf_writer).route()
    assert isinstance(router.runoff, rr.RunoffGaussianGrid)


def test_router_takes_a_runoff_at_construction(grid_case):
    """A RunoffGaussianGrid built by hand, or subclassed, is routed as given instead of one built from the weights."""
    configs = grid_configs(grid_case, 'q_given.nc')
    runoff = rr.RunoffGaussianGrid(
        grid_case['weights_file'], var_grid_runoff='ro', var_x='lon', var_y='lat', as_volumes=True
    )
    router = rr.Router(configs, runoff=runoff)
    assert router.runoff is runoff
    router.set_discharge_writer(writers.netcdf_writer).route()

    rr.Router(grid_configs(grid_case, 'q_built.nc')).set_discharge_writer(writers.netcdf_writer).route()
    with (
        xr.open_dataset(grid_case['directory'] / 'q_given.nc') as given,
        xr.open_dataset(grid_case['directory'] / 'q_built.nc') as built,
    ):
        np.testing.assert_array_equal(
            given['Q'].transpose('time', 'river_id').values, built['Q'].transpose('time', 'river_id').values
        )


def test_router_sets_as_volumes_on_a_given_runoff(grid_case):
    """Routing consumes volumes, so a RunoffGaussianGrid handed over as depths is switched, not silently wrong."""
    runoff = rr.RunoffGaussianGrid(grid_case['weights_file'], var_grid_runoff='ro', var_x='lon', var_y='lat')
    assert not runoff.as_volumes
    rr.Router(grid_configs(grid_case, 'q_depths.nc'), runoff=runoff).set_discharge_writer(writers.netcdf_writer).route()
    assert runoff.as_volumes


def test_router_rejects_a_runoff_that_is_not_one(grid_case):
    with pytest.raises(TypeError, match='must be a RunoffGaussianGrid'):
        rr.Router(grid_configs(grid_case, 'q.nc'), runoff='not a runoff')


# ── fused aggregate and route ───────────────────────────────────────────────


def dfs_ordered_tree(rng: np.random.Generator, n_rivers: int, n_outlets: int) -> np.ndarray:
    """Random tree whose downstream index is always above the river's own, mostly the next river as in DFS order."""
    downstream = np.full(n_rivers, -1, dtype=np.int32)
    for i in range(n_rivers - n_outlets):
        downstream[i] = min(n_rivers - 1, i + int(rng.geometric(0.6)))
    return downstream


@pytest.mark.parametrize('n_substeps', [1, 3])
@pytest.mark.parametrize('cumulative', [False, True])
@pytest.mark.parametrize('force_positive', [False, True])
def test_fused_grid_kernel_matches_aggregate_then_route(n_substeps, cumulative, force_positive):
    """Aggregating and routing in one pass must equal aggregating to a vlateral array and routing that."""
    from river_route.router import _river_kernels
    from river_route.runoff import _numba_kernels as runoff_kernels

    rng = np.random.default_rng(3)
    n_rivers, n_steps, n_cells, dt_runoff = 700, 48, 60, 3600
    dt_routing = dt_runoff // n_substeps
    downstream = dfs_ordered_tree(rng, n_rivers, n_outlets=4)
    k = rng.uniform(dt_routing * 0.6, dt_routing * 10, n_rivers)
    x = rng.uniform(0.0, 0.4, n_rivers)
    dt_div_k = dt_routing / k
    denominator = dt_div_k + 2 * (1 - x)
    c1 = ((dt_div_k - 2 * x) / denominator).astype(np.float32)
    c2 = ((dt_div_k + 2 * x) / denominator).astype(np.float32)
    c3 = ((2 * (1 - x) - dt_div_k) / denominator).astype(np.float32)
    c4_dt = ((c1 + c2) / dt_runoff).astype(np.float32)
    runoff_by_cell = rng.normal(RUNOFF_DEPTH, RUNOFF_DEPTH, (n_cells, n_steps)).astype(np.float32)
    if cumulative:
        runoff_by_cell = np.ascontiguousarray(np.cumsum(np.abs(runoff_by_cell), axis=1, dtype=np.float32))
    indptr = np.concatenate(([0], np.cumsum(rng.integers(1, 4, n_rivers)))).astype(np.int32)
    cell = rng.integers(0, n_cells, indptr[-1]).astype(np.int32)
    weight = rng.uniform(0.1, 1.0, indptr[-1]).astype(np.float32)
    scale = rng.uniform(1e5, 1e7, n_rivers).astype(np.float32)
    q_init = rng.uniform(0, 50, n_rivers).astype(np.float32)
    aggregation = dict(
        indptr=indptr, cell=cell, weight=weight, scale=scale, cumulative=cumulative, force_positive=force_positive
    )

    vlateral = np.empty((n_steps, n_rivers), dtype=np.float32)
    runoff_kernels.aggregate_to_rivers(
        runoff_by_cell=runoff_by_cell,
        r_start=0,
        r_stop=n_rivers,
        scratch=np.empty((64, n_steps), np.float32),
        out=vlateral,
        zero=np.float32(0),
        **aggregation,
    )
    q_expected = q_init.copy()
    expected = np.empty((n_rivers, n_steps), dtype=np.float32)
    _river_kernels.static_vlateral(
        q_t=q_expected,
        discharge_array=expected,
        downstream_indices=downstream,
        c1=c1,
        c2=c2,
        c3=c3,
        c4_dt=c4_dt,
        n_substeps=n_substeps,
        vlateral=vlateral,
        by_river=False,
        block=_river_kernels.BLOCK,
    )

    q_fused = q_init.copy()
    fused = np.empty((n_rivers, n_steps), dtype=np.float32)
    _river_kernels.static_grid(
        q_t=q_fused,
        discharge_array=fused,
        downstream_indices=downstream,
        c1=c1,
        c2=c2,
        c3=c3,
        c4_dt=c4_dt,
        n_substeps=n_substeps,
        runoff_by_cell=runoff_by_cell,
        block=_river_kernels.BLOCK,
        **aggregation,
    )
    scale_q = float(np.abs(expected).max())
    np.testing.assert_allclose(fused, expected, rtol=1e-5, atol=1e-6 * scale_q)
    np.testing.assert_allclose(q_fused, q_expected, rtol=1e-5, atol=1e-6 * scale_q)


def test_fused_inflow_rows_follow_the_open_confluences():
    """The fused kernel keeps one inflow row per river whose upstreams are partly routed, never one per river."""
    from river_route.router._river_kernels import plan_inflow_rows

    def whole(downstream):
        n = downstream.shape[0]
        return plan_inflow_rows(downstream, np.array([0]), np.array([n]), np.array([-1]), np.zeros(0, np.int32))

    # while a river is routed its own row and its downstream's row are both live
    assert whole(np.array([1, 2, 3, -1], dtype=np.int32)) == 2
    # 0 -> 1 -> 5 and 2 -> 3 -> 4 -> 5: river 5's row stays open while the second branch is routed
    assert whole(np.array([1, 5, 3, 4, 5, -1], dtype=np.int32)) == 3


def test_route_irregular_grid_falls_back_to_resampled_vlateral(grid_case, tmp_path):
    """A file that must be resampled cannot be fused, so it is aggregated first and routed like a vlateral file."""
    irregular = tmp_path / 'runoff_irregular.nc'
    write_irregular_runoff_grid(irregular)
    vlateral_file = tmp_path / 'vlateral_irregular.nc'
    prepare(grid_case, as_volumes=True).to_dataset(irregular).to_netcdf(vlateral_file)

    shared = dict(params_file=str(grid_case['params_file']), log=False, progress_bar=False)
    from_grid = tmp_path / 'q_irregular_grid.nc'
    from_vlateral = tmp_path / 'q_irregular_vlateral.nc'
    rr.Router(
        rr.Configs(
            forcing='vlateral',
            grid_weights_file=str(grid_case['weights_file']),
            grid_runoff_files=[str(irregular)],
            discharge_files=[str(from_grid)],
            var_grid_runoff='ro',
            var_x='lon',
            var_y='lat',
            **shared,
        )
    ).set_discharge_writer(writers.netcdf_writer).route()
    rr.Router(
        rr.Configs(
            forcing='vlateral', vlateral_files=[str(vlateral_file)], discharge_files=[str(from_vlateral)], **shared
        )
    ).set_discharge_writer(writers.netcdf_writer).route()
    with xr.open_dataset(from_grid) as a, xr.open_dataset(from_vlateral) as b:
        assert a['Q'].transpose('time', 'river_id').shape == (7, N_RIVERS)
        np.testing.assert_array_equal(
            a['Q'].transpose('time', 'river_id').values, b['Q'].transpose('time', 'river_id').values
        )


def test_nan_runoff_is_replaced_once_when_prepared(grid_case, tmp_path):
    """NaN runoff becomes zero in the cell series before any kernel runs, on the fused path and the aggregate path."""
    preparer = prepare(grid_case, as_volumes=True)
    runoff = np.full((NT, preparer.x_index.shape[0]), RUNOFF_DEPTH, dtype=np.float32)
    runoff[2, 0] = np.nan
    by_cell = preparer._by_cell(runoff)
    assert by_cell.flags.c_contiguous and not np.isnan(by_cell).any() and by_cell[0, 2] == 0
    vlateral, _ = preparer.aggregate(runoff, hourly_times(NT))
    assert not np.isnan(vlateral).any()
    assert (vlateral[2] > 0).all(), 'a NaN cell must not zero the other cells of its catchments'


def test_nan_vlateral_is_replaced_when_read(tmp_path):
    from conftest import write_vlateral

    volumes = np.ones((4, 3), dtype=np.float32)
    volumes[1, 2] = np.nan
    path = tmp_path / 'vlateral_nan.nc'
    write_vlateral(path, volumes, np.arange(1, 4))
    ((_, vlateral, _),) = list(rr.RunoffVlateral().reader([path]))
    assert not np.isnan(vlateral).any() and vlateral[1, 2] == 0
