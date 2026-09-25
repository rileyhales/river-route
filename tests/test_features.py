"""
Each feature does what it says, and refuses what it should: the configs, the network, the runoff readers, the router,
the discharge writers, basin subsetting, the metrics, and the command line, followed by checks of the package as a
whole. Every test uses the Willamette River cut from the reference solution, or a copy of its real inputs with one
thing deliberately broken.
"""

import importlib
import json
import os
import pkgutil
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
import zarr
from reference_solution import GRID_NAMES, WILLAMETTE, Basin, assert_same, route

import river_route as rr
from river_route import metrics
from river_route.router import writers

DT = 3600
ROOT = Path(__file__).parent.parent
TOOLS = Path(sys.executable).parent  # the environment's ruff, mypy, and zensical


def routing_options(basin: Basin, directory: Path) -> dict:
    """Configs options that are valid for routing the basin's March runoff."""
    return basin.gridded(months=1) | {'discharge_files': [directory / 'discharge.zarr']}


################################################
# Configs
################################################

REFUSED_CONFIGS = {
    'runoff forcing needs a runoff type': (
        lambda basin, d: {'runoff_type': None},
        ValueError,
        'runoff_type is required',
    ),
    'selectors take only their values': (
        lambda basin, d: {'coefficients': 'sideways'},
        ValueError,
        'coefficients must be one of',
    ),
    'catchment runoff takes no grid weights': (
        lambda basin, d: {'runoff_type': 'catchment'},
        ValueError,
        'grid_weights_file is not used',
    ),
    'gridded runoff needs grid weights': (
        lambda basin, d: {'grid_weights_file': None},
        ValueError,
        'grid_weights_file is required',
    ),
    'one output per runoff file': (
        lambda basin, d: {'runoff_files': basin.runoff_files[:2]},
        ValueError,
        'must match number of input files',
    ),
    'outputs must be distinct': (
        lambda basin, d: {'runoff_files': basin.runoff_files[:2], 'discharge_files': [d / 'q.zarr', d / 'q.zarr']},
        ValueError,
        'duplicate paths',
    ),
    'channel routing needs an initial state': (
        lambda basin, d: {'forcing': 'channel', 'runoff_files': [], 'dt_routing': DT, 'dt_total': 86400},
        ValueError,
        'channel_state_init_file is required',
    ),
    'unit hydrographs need a kernel file': (
        lambda basin, d: {'transform': 'unit_hydrograph'},
        ValueError,
        'uh_kernel_file is required',
    ),
    'input files must exist': (lambda basin, d: {'params_file': d / 'missing.parquet'}, FileNotFoundError, 'not found'),
    'output folders must exist': (
        lambda basin, d: {'discharge_files': [d / 'missing' / 'q.zarr']},
        NotADirectoryError,
        'Directory not found',
    ),
}


@pytest.mark.parametrize(('changes', 'error', 'message'), REFUSED_CONFIGS.values(), ids=REFUSED_CONFIGS.keys())
def test_configs_refuse_inconsistent_options(changes, error, message, willamette: Basin, tmp_path: Path) -> None:
    with pytest.raises(error, match=message):
        rr.Configs(**(routing_options(willamette, tmp_path) | changes(willamette, tmp_path))).validate_routing()


def test_configs_round_trip_through_json(willamette: Basin, tmp_path: Path) -> None:
    configs = rr.Configs(**routing_options(willamette, tmp_path))
    configs.to_json(tmp_path / 'configs.json')
    assert rr.Configs.from_json(tmp_path / 'configs.json') == configs


def test_unknown_config_keys_are_refused(tmp_path: Path) -> None:
    (tmp_path / 'configs.json').write_text(json.dumps({'params_file': 'routing.parquet', 'not_an_option': 1}))
    with pytest.raises(ValueError, match='Unrecognized config key'):
        rr.Configs.from_json(tmp_path / 'configs.json')


def test_deep_validate_accepts_the_real_inputs(willamette: Basin) -> None:
    rr.Configs(params_file=willamette.params_file, grid_weights_file=willamette.weights_file).deep_validate()


BROKEN_PARAMS = {
    'a downstream river is missing': (
        lambda table: table.drop(index=table.index[table['river_id'].isin(table['next_river_id'])][0]),
        'next_river_id values must exist',
    ),
    'rows out of upstream to downstream order': (lambda table: table.iloc[::-1], 'not topologically sorted'),
    'negative travel time': (lambda table: table.assign(k=-table['k']), 'k column must be positive'),
}


@pytest.mark.parametrize(('breaks', 'message'), BROKEN_PARAMS.values(), ids=BROKEN_PARAMS.keys())
def test_deep_validate_finds_broken_params(breaks, message, willamette: Basin, tmp_path: Path) -> None:
    breaks(pd.read_parquet(willamette.params_file)).to_parquet(tmp_path / 'routing.parquet')
    with pytest.raises(ValueError, match=message):
        rr.Configs(params_file=tmp_path / 'routing.parquet').deep_validate()


def test_deep_validate_finds_weights_that_do_not_sum_to_one(willamette: Basin, tmp_path: Path) -> None:
    with xr.open_dataset(willamette.weights_file) as weights:
        broken = weights.load()
    broken['proportion'].values[0] *= 0.5
    broken.to_netcdf(tmp_path / 'weights.nc')
    with pytest.raises(ValueError, match='must sum to 1'):
        rr.Configs(params_file=willamette.params_file, grid_weights_file=tmp_path / 'weights.nc').deep_validate()


################################################
# Network
################################################

BROKEN_TABLES = {
    'duplicate river ids': (lambda table: pd.concat([table, table.iloc[:1]]), 'duplicate river_id'),
    'an unknown downstream id': (
        lambda table: table.assign(
            next_river_id=np.where(np.arange(len(table)) == 0, table['river_id'].max() + 1, table['next_river_id'])
        ),
        'is not in the river_id column',
    ),
    'rows out of upstream to downstream order': (lambda table: table.iloc[::-1], 'topologically sorted'),
}


@pytest.mark.parametrize(('breaks', 'message'), BROKEN_TABLES.values(), ids=BROKEN_TABLES.keys())
def test_networks_refuse_broken_tables(breaks, message, willamette: Basin) -> None:
    with pytest.raises(ValueError, match=message):
        rr.Network(breaks(pd.read_parquet(willamette.params_file)))


@pytest.mark.parametrize('threads', [1, 4])
def test_the_routing_schedule_routes_every_river_once(threads: int, package: Path, manifest: dict) -> None:
    network = rr.Network(package / manifest['runs']['static_standard']['configs']['params_file'])
    jobs, cut_target = network.routing_schedule(threads=threads, concurrent=True)
    rivers = np.concatenate(
        [np.arange(start, stop) for job in jobs for start, stop in zip(job[0], job[1], strict=True)]
    )
    np.testing.assert_array_equal(np.sort(rivers), np.arange(network.size))
    downstream = network.downstream_indices
    for starts, stops, outlet, region in jobs[:-1]:
        assert starts.shape == (1,), 'a region is one contiguous block of rivers'
        members = np.arange(starts[0], stops[0])
        leaving = members[(downstream[members] < starts[0]) | (downstream[members] >= stops[0])]
        np.testing.assert_array_equal(leaving, [outlet])  # only a region's outlet drains out of it
        assert cut_target[region] == downstream[outlet]
    assert (len(jobs) > 1) == (threads > 1)


@pytest.mark.parametrize('dt', [3600, 900])
def test_conditioning_makes_every_river_stable(dt: int, package: Path, manifest: dict) -> None:
    network = rr.Network(package / manifest['runs']['static_standard']['configs']['params_file'])
    subdivisions, substeps, resolvable = network.conditioning(dt)
    assert resolvable.all()
    too_long, too_short = network.unstable_mask(dt, subdivisions, substeps)
    assert not too_long.any() and not too_short.any()


def test_the_largest_stable_dt_divides_the_runoff_step(willamette: Basin) -> None:
    dt = willamette.network.largest_stable_dt(DT)
    assert DT % dt == 0 and willamette.network.conditioning(dt)[2].all()


def test_write_stabilized_writes_only_its_output(willamette: Basin, tmp_path: Path) -> None:
    source = tmp_path / 'routing.parquet'
    shutil.copy2(willamette.params_file, source)
    original = source.read_bytes()
    written = pd.read_parquet(rr.Network(source).write_stabilized(DT, tmp_path / 'stabilized.parquet'))
    assert source.read_bytes() == original
    stabilized = rr.Network(source).stabilize(DT)
    np.testing.assert_array_equal(written['river_id'], stabilized.river_ids)
    np.testing.assert_array_equal(written['synthetic'], stabilized.synthetic)
    with pytest.raises(ValueError, match='needs a path'):
        rr.Network(pd.read_parquet(source)).write_stabilized(DT)


def test_networks_report_unstable_rivers(willamette: Basin) -> None:
    network = rr.Network(willamette.params_file)
    assert f'n_rivers={network.size:,}' in repr(network)
    with pytest.warns(UserWarning, match='not Muskingum-stable'):
        network.check_stability(DT, action='warn')
    with pytest.raises(ValueError, match='not Muskingum-stable'):
        network.check_stability(DT, action='raise')


################################################
# Runoff
################################################


@pytest.fixture(scope='module')
def march(willamette: Basin) -> tuple[np.ndarray, np.ndarray, float]:
    """The basin's March runoff depths at the grid cells, NaN free, with their times and unit conversion factor."""
    cells, time_index, factor = grid(willamette).read_runoff(willamette.runoff_files[0])
    return np.nan_to_num(cells), time_index, factor


def grid(basin: Basin, **options) -> rr.GaussianGridRunoff:
    return rr.GaussianGridRunoff(basin.weights_file, **GRID_NAMES, **options)


def test_depth_units_are_converted(willamette: Basin) -> None:
    meters, _ = grid(willamette).catchment_runoff(willamette.runoff_files[0])
    millimeters, _ = grid(willamette, runoff_depth_unit='mm').catchment_runoff(willamette.runoff_files[0])
    assert_same(millimeters, 0.001 * meters)


def test_cumulative_runoff_becomes_runoff_per_step(willamette: Basin, march) -> None:
    cells, time_index, factor = march
    per_step, _ = grid(willamette).aggregate(cells, time_index, factor)
    totals = np.cumsum(cells, axis=0, dtype=np.float32)
    from_totals, _ = grid(willamette, grid_accumulation_type='cumulative').aggregate(totals, time_index, factor)
    assert_same(from_totals, per_step)


def test_missing_cells_contribute_nothing(willamette: Basin, march) -> None:
    cells, time_index, factor = march
    missing, zero = cells.copy(), cells.copy()
    missing[:, 0], zero[:, 0] = np.nan, 0
    np.testing.assert_array_equal(
        grid(willamette).aggregate(missing, time_index, factor)[0],
        grid(willamette).aggregate(zero, time_index, factor)[0],
    )


def test_force_positive_runoff_clips_negative_catchment_runoff(willamette: Basin, march) -> None:
    cells, time_index, factor = march
    shifted = cells - cells.mean(dtype=np.float32)  # plenty of negative runoff
    unclipped, _ = grid(willamette).aggregate(shifted, time_index, factor)
    clipped, _ = grid(willamette, force_positive_runoff=True).aggregate(shifted, time_index, factor)
    np.testing.assert_array_equal(clipped, np.maximum(unclipped, 0))


def test_volumes_are_depths_times_catchment_area(willamette: Basin, march) -> None:
    cells, time_index, factor = march
    depths, _ = grid(willamette).aggregate(cells, time_index, factor)
    volume_grid = grid(willamette, as_volumes=True)
    volumes, _ = volume_grid.aggregate(cells, time_index, factor)
    assert_same(volumes, depths * volume_grid.catchment_area[:, np.newaxis])


def test_irregular_time_steps_are_resampled_keeping_the_total(willamette: Basin, march) -> None:
    cells, time_index, factor = march
    keep = np.ones(time_index.size, dtype=bool)
    keep[5::7] = False  # drop every seventh step after the fifth
    every_step, _ = grid(willamette).aggregate(cells, time_index, factor)
    resampled, resampled_times = grid(willamette).aggregate(cells[keep], time_index[keep], factor)
    assert np.all(np.diff(resampled_times) == np.timedelta64(DT, 's'))
    assert (resampled_times[0], resampled_times[-1]) == (time_index[0], time_index[-1])
    assert_same(resampled.sum(axis=1), every_step[:, keep].sum(axis=1))


def test_depth_files_are_read_as_volumes(willamette: Basin, march, tmp_path: Path) -> None:
    cells, time_index, factor = march
    depth_grid = grid(willamette)
    depths, dates = depth_grid.aggregate(cells, time_index, factor)
    depth_grid.to_netcdf(tmp_path / 'depths.nc', dates, depths, depth_grid.river_ids, depth_grid.catchment_area)
    volumes, _ = grid(willamette, as_volumes=True).aggregate(cells, time_index, factor)
    ((_, read, _),) = rr.CatchmentRunoff().generator([tmp_path / 'depths.nc'])
    assert_same(read.runoff, volumes)


def test_catchment_runoff_volumes_check_their_layout() -> None:
    volumes = rr.runoff.CatchmentRunoffVolumes(np.zeros((3, 10), dtype=np.float32))
    volumes.check(3, 10)
    with pytest.raises(ValueError, match='has shape'):
        volumes.check(4, 10)
    with pytest.raises(ValueError, match='row contiguous'):
        rr.runoff.CatchmentRunoffVolumes(np.zeros((10, 3), dtype=np.float32).T).check(3, 10)
    first = volumes.first_steps(4)
    assert first.runoff.shape == (3, 4) and first.runoff.strides[1] == first.runoff.itemsize


def test_grid_cell_runoff_checks_its_weight_table(willamette: Basin) -> None:
    dates, cell_runoff, _ = next(grid(willamette, as_volumes=True).generator(willamette.runoff_files[:1]))
    n_rivers = willamette.network.size
    cell_runoff.check(n_rivers, dates.size)
    with pytest.raises(ValueError, match='does not describe'):
        cell_runoff.check(n_rivers + 1, dates.size)
    assert cell_runoff.first_steps(10) is cell_runoff


################################################
# Router
################################################

UNSUPPORTED_OPTIONS = {
    'dynamic coefficients on a stabilized network': (
        lambda basin: {'coefficients': 'dynamic', 'network_type': 'stabilized'},
        'cannot route a stabilized network',
    ),
    'the unit hydrograph transform': (
        lambda basin: {'transform': 'unit_hydrograph', 'uh_kernel_file': basin.params_file},
        'transform is not implemented',
    ),
    'reduced gaussian grid runoff': (lambda basin: {'runoff_type': 'reduced_gaussian_grid'}, 'not implemented yet'),
}


@pytest.mark.parametrize(('changes', 'message'), UNSUPPORTED_OPTIONS.values(), ids=UNSUPPORTED_OPTIONS.keys())
def test_unsupported_options_raise_before_routing(changes, message, willamette: Basin, tmp_path: Path) -> None:
    with pytest.raises(NotImplementedError, match=message):
        route(tmp_path, **(willamette.gridded(months=1) | changes(willamette)))


@pytest.mark.parametrize('threads', [0, True, 2.0])
def test_threads_must_be_a_positive_integer(threads, willamette: Basin, tmp_path: Path) -> None:
    router = rr.Router(rr.Configs(**routing_options(willamette, tmp_path)))
    with pytest.raises(ValueError, match='threads must be an integer'):
        router.route(threads=threads)


def test_a_runoff_of_the_wrong_class_is_refused(willamette: Basin, tmp_path: Path) -> None:
    options = routing_options(willamette, tmp_path) | {'runoff_type': 'catchment', 'grid_weights_file': None}
    with pytest.raises(TypeError, match='is read by CatchmentRunoff'):
        rr.Router(rr.Configs(**options), runoff=grid(willamette))


def test_an_initial_state_of_the_wrong_size_is_refused(willamette: Basin, tmp_path: Path) -> None:
    pd.DataFrame({'Q': np.zeros(willamette.network.size - 1, dtype=np.float32)}).to_parquet(tmp_path / 'state.parquet')
    with pytest.raises(ValueError, match='channel_state_init_file has'):
        route(
            tmp_path,
            willamette.months[:1],
            params_file=willamette.params_file,
            channel_state_init_file=tmp_path / 'state.parquet',
        )


def test_split_rivers_write_only_the_original_rivers(willamette: Basin, tmp_path: Path) -> None:
    network = rr.Network(willamette.params_file).stabilize(DT)
    routed = route(tmp_path, network=network, writer=writers.zarr_writer, **willamette.gridded(months=1))
    assert routed.discharge[0].shape[0] == willamette.network.size
    with xr.open_zarr(tmp_path / routed.discharge_files[0]) as written:
        np.testing.assert_array_equal(written['river_id'].to_numpy(), willamette.network.river_ids)


def test_the_writer_is_called_once_per_runoff_file(willamette: Basin, tmp_path: Path) -> None:
    routed = route(tmp_path, willamette.months, params_file=willamette.params_file)
    assert routed.discharge_files == [f'discharge_{i}.zarr' for i in range(len(willamette.months))]
    assert [dates.size for dates in routed.dates] == [dates.size for dates, _ in willamette.months]


################################################
# Discharge writers
################################################


@pytest.fixture(scope='module')
def march_discharge(willamette: Basin, tmp_path_factory: pytest.TempPathFactory):
    """The basin's routed March discharge, as the router hands it to a writer."""
    return route(tmp_path_factory.mktemp('march'), willamette.months[:1], params_file=willamette.params_file)


def test_the_zarr_writer_layout(willamette: Basin, march_discharge, tmp_path: Path) -> None:
    written = tmp_path / 'discharge.zarr'
    writers.zarr_writer(march_discharge.router, march_discharge.dates[0], march_discharge.discharge[0], written)
    zarr.open_consolidated(str(written))
    with xr.open_zarr(written) as discharge:
        assert discharge['Q'].dims == ('river_id', 'time')
        np.testing.assert_array_equal(discharge['time'].to_numpy(), march_discharge.dates[0].astype('datetime64[ns]'))
        np.testing.assert_array_equal(discharge['river_id'].to_numpy(), willamette.network.river_ids)
        rounded = writers.bitround(march_discharge.discharge[0], writers.ZARR_KEEPBITS)
        np.testing.assert_array_equal(discharge['Q'].to_numpy(), rounded)


def test_bitround_error_is_bounded(march_discharge) -> None:
    values = march_discharge.discharge[0]
    rounded = writers.bitround(values, 15)
    assert np.all(np.abs(rounded.astype(np.float64) - values) <= np.abs(values.astype(np.float64)) * 2.0**-16)


def test_the_netcdf_writer_round_trips(willamette: Basin, march_discharge, tmp_path: Path) -> None:
    writers.netcdf_writer(
        march_discharge.router, march_discharge.dates[0], march_discharge.discharge[0], tmp_path / 'q.nc'
    )
    with xr.open_dataset(tmp_path / 'q.nc') as discharge:
        np.testing.assert_array_equal(discharge['Q'].to_numpy(), march_discharge.discharge[0])
        np.testing.assert_array_equal(discharge['river_id'].to_numpy(), willamette.network.river_ids)


def test_the_parquet_writer_round_trips(willamette: Basin, march_discharge, tmp_path: Path) -> None:
    writers.parquet_writer(
        march_discharge.router, march_discharge.dates[0], march_discharge.discharge[0], tmp_path / 'q.parquet'
    )
    table = pd.read_parquet(tmp_path / 'q.parquet')
    np.testing.assert_array_equal(table['river_id'].to_numpy(), willamette.network.river_ids)
    np.testing.assert_array_equal(table.drop(columns='river_id').to_numpy(), march_discharge.discharge[0])


def test_to_time_major_is_the_transpose(march_discharge) -> None:
    np.testing.assert_array_equal(writers.to_time_major(march_discharge.discharge[0]), march_discharge.discharge[0].T)


################################################
# Subsetting, metrics, and the command line
################################################


def test_a_subset_is_a_closed_upstream_basin(willamette: Basin, package: Path, manifest: dict) -> None:
    table = pd.read_parquet(willamette.params_file)
    ids = set(table['river_id'])
    assert set(table['next_river_id']) <= ids | {-1}
    np.testing.assert_array_equal(table.loc[table['next_river_id'] == -1, 'river_id'], [WILLAMETTE])
    with xr.open_dataset(willamette.weights_file) as weights:
        assert set(np.unique(weights['river_id'].to_numpy())) == ids
    columbia = rr.Network(package / manifest['runs']['static_standard']['configs']['params_file'])
    in_basin = columbia.river_ids == WILLAMETTE
    for river in range(
        columbia.size - 1, -1, -1
    ):  # downstream rivers come later, so each is decided before its upstreams
        downstream = columbia.downstream_indices[river]
        in_basin[river] |= downstream >= 0 and in_basin[downstream]
    assert in_basin.sum() == len(table), 'every river upstream of the outlet is in the subset'


def test_metrics_of_identical_series(march_discharge) -> None:
    series = march_discharge.discharge[0][-1].astype(np.float64)
    assert metrics.kge2012(series, series) == pytest.approx(1)
    assert metrics.me(series, series) == metrics.mae(series, series) == metrics.mse(series, series) == 0


def test_the_command_line_routes_a_config_file(willamette: Basin, tmp_path: Path) -> None:
    configs = rr.Configs(**routing_options(willamette, tmp_path), log=False, progress_bar=False)
    configs.to_json(tmp_path / 'configs.json')
    command = [sys.executable, '-m', 'river_route._cli', 'route', str(tmp_path / 'configs.json')]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    expected = route(tmp_path, **willamette.gridded(months=1)).discharge[0]
    with xr.open_zarr(tmp_path / 'discharge.zarr') as written:
        assert_same(written['Q'].to_numpy(), expected)


################################################
# The package as a whole
################################################

MODULES = [module.name for module in pkgutil.walk_packages(rr.__path__, 'river_route.')]


@pytest.mark.parametrize('module', MODULES)
def test_every_module_imports_in_a_fresh_process(module: str) -> None:
    result = subprocess.run([sys.executable, '-c', f'import {module}'], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize('module', ['river_route', *MODULES])
def test_every_name_in_all_exists(module: str) -> None:
    imported = importlib.import_module(module)
    missing = [name for name in getattr(imported, '__all__', []) if not hasattr(imported, name)]
    assert not missing, f'{module}.__all__ lists names it does not define: {missing}'


STATIC_CHECKS = {
    'ruff check': ['ruff', 'check', 'river_route', 'tests', 'examples'],
    'ruff format': ['ruff', 'format', '--check', 'river_route', 'tests', 'examples'],
    'mypy': ['mypy', 'river_route', '--warn-unreachable', f'--cache-dir={TOOLS.parent / ".mypy_cache"}'],
    'docs build': ['zensical', 'build', '--strict'],
}


@pytest.mark.parametrize('command', STATIC_CHECKS.values(), ids=STATIC_CHECKS.keys())
def test_static_checks_pass(command: list[str]) -> None:
    result = subprocess.run([str(TOOLS / command[0]), *command[1:]], cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.slow
def test_freshly_compiled_kernels_match_the_reference(tmp_path: Path) -> None:
    """Compile every kernel into an empty numba cache and route the reference solution with them. Freshly compiled and
    cached kernels have differed in the last bit before, for sub-cycled rivers."""
    command = [
        sys.executable,
        '-m',
        'pytest',
        str(ROOT / 'tests' / 'test_reference_solution.py'),
        '-q',
        '-p',
        'no:cacheprovider',
    ]
    result = subprocess.run(
        command, cwd=ROOT, env=os.environ | {'NUMBA_CACHE_DIR': str(tmp_path)}, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout[-4000:]
