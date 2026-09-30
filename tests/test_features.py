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
from river_route.network import streams
from river_route.router import writers

DT = 3600
ROOT = Path(__file__).parent.parent
TOOLS = Path(sys.executable).parent  # the environment's ruff, mypy, and zensical

REFUSED_CONFIGS = {
    'forcing names the form of its runoff': (
        lambda basin, d: {'forcing': 'runoff'},
        ValueError,
        'forcing must be one of',
    ),
    'selectors take only their values': (
        lambda basin, d: {'coefficients': 'sideways'},
        ValueError,
        'coefficients must be one of',
    ),
    'catchment runoff takes no grid weights': (
        lambda basin, d: {'forcing': 'catchment'},
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
    'input files must exist': (
        lambda basin, d: {'network_file': d / 'missing.parquet'},
        FileNotFoundError,
        'not found',
    ),
    'output folders must exist': (
        lambda basin, d: {'discharge_files': [d / 'missing' / 'q.zarr']},
        NotADirectoryError,
        'Directory not found',
    ),
    'outputs are all discarded or all written': (
        lambda basin, d: {'runoff_files': basin.runoff_files[:2], 'discharge_files': [os.devnull, d / 'q.zarr']},
        ValueError,
        'mixes the null device',
    ),
}

BROKEN_NETWORK_FILES = {
    'a downstream river is missing': (
        lambda table: table.drop(index=table.index[table['riverId'].isin(table['nextRiverId'])][0]),
        'nextRiverId values must exist',
    ),
    'rows out of upstream to downstream order': (lambda table: table.iloc[::-1], 'not topologically sorted'),
    'a river with the outlet id': (
        lambda table: table.assign(riverId=np.where(np.arange(len(table)) == 0, -1, table['riverId'])),
        'riverId must not be -1',
    ),
    'negative travel time': (
        lambda table: table.assign(muskingumK=-table['muskingumK']),
        'muskingumK column must be positive',
    ),
    'no river index': (lambda table: table.drop(columns='riverIndex'), 'missing riverIndex column'),
    'river indexes that skip': (
        lambda table: table.assign(riverIndex=2 * table['riverIndex']),
        'riverIndex must increase by one',
    ),
    'an upstream count that is off by one': (
        lambda table: table.assign(upstreamCount=table['upstreamCount'] + (np.arange(len(table)) == len(table) - 1)),
        'upstreamCount must count',
    ),
}

BROKEN_TABLES = {
    'duplicate river ids': (lambda table: pd.concat([table, table.iloc[:1]]), 'duplicate riverId'),
    'an unknown downstream id': (
        lambda table: table.assign(
            nextRiverId=np.where(np.arange(len(table)) == 0, table['riverId'].max() + 1, table['nextRiverId'])
        ),
        'is not in the riverId column',
    ),
    'rows out of upstream to downstream order': (lambda table: table.iloc[::-1], 'topologically sorted'),
    'no upstream count': (lambda table: table.drop(columns='upstreamCount'), 'missing the required column'),
    'an upstream count that is off by one': (
        lambda table: table.assign(upstreamCount=table['upstreamCount'] + (np.arange(len(table)) == len(table) - 1)),
        'upstreamCount must count',
    ),
}

UNSUPPORTED_OPTIONS = {
    'dynamic coefficients on a stabilized network': (
        lambda basin: {'coefficients': 'dynamic', 'network_type': 'stabilized'},
        'cannot route a stabilized network',
    )
}

MODULES = [module.name for module in pkgutil.walk_packages(rr.__path__, 'river_route.')]

STATIC_CHECKS = {
    'ruff check': ['ruff', 'check', 'river_route', 'tests', 'examples'],
    'ruff format': ['ruff', 'format', '--check', 'river_route', 'tests', 'examples'],
    'mypy': ['mypy', 'river_route', '--warn-unreachable', f'--cache-dir={TOOLS.parent / ".mypy_cache"}'],
    'docs build': ['zensical', 'build', '--strict'],
}


def routing_options(basin: Basin, directory: Path) -> dict:
    """Configs options that are valid for routing the basin's March runoff."""
    return basin.gridded(months=1) | {'discharge_files': [directory / 'discharge.zarr']}


@pytest.mark.parametrize(('changes', 'error', 'message'), REFUSED_CONFIGS.values(), ids=REFUSED_CONFIGS.keys())
def test_configs_refuse_inconsistent_options(changes, error, message, willamette: Basin, tmp_path: Path) -> None:
    with pytest.raises(error, match=message):
        rr.Configs(**(routing_options(willamette, tmp_path) | changes(willamette, tmp_path))).validate_routing()


def test_discharge_dir_names_outputs_for_the_default_writer(tmp_path: Path) -> None:
    configs = rr.Configs(discharge_dir=tmp_path, runoff_files=['runoff/2020.nc', 'runoff/2021.zarr'])
    assert configs.discharge_files == [str(tmp_path / 'discharge_2020.zarr'), str(tmp_path / 'discharge_2021.zarr')]
    assert rr.Configs(discharge_dir=tmp_path).discharge_files == [str(tmp_path / 'discharge.zarr')]
    with pytest.raises(ValueError, match='same output file'):
        rr.Configs(discharge_dir=tmp_path, runoff_files=['a/2020.nc', 'b/2020.zarr'])


@pytest.mark.parametrize('outputs', ['files', 'null device'])
def test_configs_round_trip_through_json(outputs: str, willamette: Basin, tmp_path: Path) -> None:
    options = routing_options(willamette, tmp_path)
    if outputs == 'null device':
        options |= {'discharge_files': [], 'discharge_dir': Path(os.devnull)}
    configs = rr.Configs(**options)
    configs.to_json(tmp_path / 'configs.json')
    assert rr.Configs.from_json(tmp_path / 'configs.json') == configs


def test_preparing_runoff_does_not_skip_routing_validation(willamette: Basin, tmp_path: Path) -> None:
    configs = rr.Configs(**routing_options(willamette, tmp_path), coefficients='dynamic', network_type='stabilized')
    runoff = rr.GridRunoff.from_configs(configs)  # validates the configs for preparing runoff, not for routing
    with pytest.raises(NotImplementedError, match='cannot route a stabilized network'):
        rr.Router(configs, runoff=runoff).route()


def test_unknown_config_keys_are_refused(tmp_path: Path) -> None:
    (tmp_path / 'configs.json').write_text(json.dumps({'network_file': 'network.parquet', 'not_an_option': 1}))
    with pytest.raises(ValueError, match='Unrecognized config key'):
        rr.Configs.from_json(tmp_path / 'configs.json')


def test_deep_validate_accepts_the_real_inputs(willamette: Basin) -> None:
    rr.Configs(network_file=willamette.network_file, grid_weights_file=willamette.weights_file).deep_validate()


@pytest.mark.parametrize(('breaks', 'message'), BROKEN_NETWORK_FILES.values(), ids=BROKEN_NETWORK_FILES.keys())
def test_deep_validate_finds_broken_network_files(breaks, message, willamette: Basin, tmp_path: Path) -> None:
    breaks(pd.read_parquet(willamette.network_file)).to_parquet(tmp_path / 'network.parquet')
    with pytest.raises(ValueError, match=message):
        rr.Configs(network_file=tmp_path / 'network.parquet').deep_validate()


def test_deep_validate_checks_a_state_against_the_network_type(willamette: Basin, tmp_path: Path) -> None:
    final_state = tmp_path / 'final_state.parquet'
    options = {'network_file': willamette.network_file, 'channel_state_final_file': final_state}
    route(tmp_path, willamette.months[:1], network_type='stabilized', **options)
    stabilized_state = {'network_file': willamette.network_file, 'channel_state_init_file': final_state}
    rr.Configs(**stabilized_state, network_type='stabilized').deep_validate()  # one row per sub-reach
    with pytest.raises(ValueError, match='same number of rows'):
        rr.Configs(**stabilized_state).deep_validate()


def test_deep_validate_finds_weights_that_do_not_sum_to_one(willamette: Basin, tmp_path: Path) -> None:
    with xr.open_dataset(willamette.weights_file) as weights:
        broken = weights.load()
    broken['proportion'].values[0] *= 0.5
    broken.to_netcdf(tmp_path / 'weights.nc')
    with pytest.raises(ValueError, match='must sum to 1'):
        rr.Configs(network_file=willamette.network_file, grid_weights_file=tmp_path / 'weights.nc').deep_validate()


@pytest.mark.parametrize(('breaks', 'message'), BROKEN_TABLES.values(), ids=BROKEN_TABLES.keys())
def test_networks_refuse_broken_tables(breaks, message, willamette: Basin) -> None:
    with pytest.raises(ValueError, match=message):
        rr.Network(breaks(pd.read_parquet(willamette.network_file)))


@pytest.mark.parametrize('threads', [1, 4])
def test_the_routing_blocks_route_every_river_once(threads: int, package: Path, manifest: dict) -> None:
    network = rr.Network(package / manifest['runs']['static_standard']['configs']['network_file'])
    job_blocks, cut_target = network.routing_blocks(threads)
    *sub_watershed_jobs, main_stem = job_blocks

    def rivers_of(starts: np.ndarray, stops: np.ndarray) -> np.ndarray:
        return np.concatenate([np.arange(start, stop) for start, stop in zip(starts, stops, strict=True)])

    rivers = np.concatenate([rivers_of(starts, stops) for starts, stops, _, _ in job_blocks])
    np.testing.assert_array_equal(np.sort(rivers), np.arange(network.size))
    downstream = network.downstream_indices
    main_stem_rivers = rivers_of(main_stem[0], main_stem[1])
    for job in sub_watershed_jobs:
        for start, stop, outlet, block_number in zip(*job, strict=True):
            members = np.arange(start, stop)
            leaving = members[(downstream[members] < start) | (downstream[members] >= stop)]
            np.testing.assert_array_equal(leaving, [outlet])  # only a block's outlet drains out of it
            entering = np.flatnonzero((downstream >= start) & (downstream < stop))
            assert np.all((entering >= start) & (entering < stop)), 'a block holds every river upstream of it'
            assert cut_target[block_number] == downstream[outlet]
            assert cut_target[block_number] < 0 or cut_target[block_number] in main_stem_rivers
    assert (len(job_blocks) > 1) == (threads > 1)


@pytest.mark.parametrize('dt', [3600, 900])
def test_conditioning_makes_every_river_stable(dt: int, package: Path, manifest: dict) -> None:
    network = rr.Network(package / manifest['runs']['static_standard']['configs']['network_file'])
    substeps, subcycles, resolvable = network.conditioning(dt)
    assert resolvable.all()
    too_long, too_short = network.unstable_mask(dt, substeps, subcycles)
    assert not too_long.any() and not too_short.any()


def test_the_largest_stable_dt_divides_the_runoff_step(willamette: Basin) -> None:
    dt = willamette.network.largest_stable_dt(DT)
    assert DT % dt == 0 and willamette.network.conditioning(dt)[2].all()


def test_write_stabilized_writes_only_its_output(willamette: Basin, tmp_path: Path) -> None:
    source = tmp_path / 'network.parquet'
    shutil.copy2(willamette.network_file, source)
    original = source.read_bytes()
    written = pd.read_parquet(rr.Network(source).write_stabilized(DT, tmp_path / 'stabilized.parquet'))
    assert source.read_bytes() == original
    stabilized = rr.Network(source).stabilize(DT)
    np.testing.assert_array_equal(written['riverId'], stabilized.river_ids)
    np.testing.assert_array_equal(written['synthetic'], stabilized.synthetic)
    assert np.all(np.diff(written['riverIndex']) == 1), 'the added sub-reaches are numbered in with the rivers'
    rr.Configs(network_file=tmp_path / 'stabilized.parquet', network_type='stabilized').deep_validate()
    with pytest.raises(ValueError, match='needs a path'):
        rr.Network(pd.read_parquet(source)).write_stabilized(DT)


def test_stabilize_splits_each_river_into_stable_pieces_of_its_travel_time(willamette: Basin) -> None:
    river_ids, k, x = willamette.network.river_ids, willamette.network.k.astype(np.float64), willamette.network.x
    pieces = {}
    for mode in ('uniform', 'nonuniform'):
        network = rr.Network(willamette.network_file).stabilize(DT, mode=mode)
        river = pd.Index(river_ids).get_indexer(network.parent_river_ids)  # the original river of each piece
        pieces[mode] = np.bincount(river, minlength=river_ids.shape[0])
        assert_same(np.bincount(river, weights=network.k.astype(np.float64)), k)
        split = pieces[mode][river] > 1  # only a river whose every piece is stable is split
        piece_k, piece_x = network.k[split].astype(np.float64), x[river[split]].astype(np.float64)
        assert np.all(piece_k >= DT / (2 * (1 - piece_x)) * (1 - 1e-6))
        assert np.all(piece_k * 2 * piece_x <= DT * (1 + 1e-6))
    assert np.any(pieces['uniform'] > 1), 'no river of the basin is too long to route at dt'
    np.testing.assert_array_equal(pieces['nonuniform'], pieces['uniform'])


def test_stabilize_follows_explicit_weights(willamette: Basin) -> None:
    n_rivers, k = willamette.network.size, willamette.network.k.astype(np.float64)
    network = rr.Network(willamette.network_file).stabilize(DT, weights=[[1, 3]] * n_rivers)
    np.testing.assert_allclose(network.k[0::2], k / 4, rtol=1e-6)  # the added piece upstream of each river
    np.testing.assert_allclose(network.k[1::2], 3 * k / 4, rtol=1e-6)  # the river itself, now its outlet piece
    np.testing.assert_array_equal(network.synthetic, np.tile([True, False], n_rivers))
    np.testing.assert_array_equal(network.river_ids[1::2], willamette.network.river_ids)
    np.testing.assert_array_equal(network.next_river_ids[0::2], network.river_ids[1::2])


def test_networks_report_unstable_rivers(willamette: Basin) -> None:
    network = rr.Network(willamette.network_file)
    assert f'n_rivers={network.size:,}' in repr(network)
    with pytest.warns(UserWarning, match='not Muskingum-stable'):
        network.check_stability(DT, action='warn')
    with pytest.raises(ValueError, match='not Muskingum-stable'):
        network.check_stability(DT, action='raise')


@pytest.fixture(scope='module')
def march(willamette: Basin) -> tuple[np.ndarray, np.ndarray, float]:
    """The basin's March runoff depths at the grid cells, NaN free, with their times and unit conversion factor."""
    cells, time_index, factor = grid(willamette).read_runoff(willamette.runoff_files[0])
    return np.nan_to_num(cells), time_index, factor


def grid(basin: Basin, **options) -> rr.GridRunoff:
    return rr.GridRunoff(basin.weights_file, **GRID_NAMES, **options)


@pytest.mark.parametrize('unit', ['mm', 'kg m-2'])  # a kilogram of water on a square meter is a millimeter deep
def test_depth_units_are_converted(unit: str, willamette: Basin) -> None:
    meters, _ = grid(willamette).catchment_runoff(willamette.runoff_files[0])
    millimeters, _ = grid(willamette, runoff_depth_unit=unit).catchment_runoff(willamette.runoff_files[0])
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


def test_reading_runoff_for_routing_leaves_depths_as_depths(willamette: Basin, march) -> None:
    cells, time_index, factor = march
    depth_grid = grid(willamette)
    depths, _ = depth_grid.aggregate(cells, time_index, factor)
    next(depth_grid.generator(willamette.runoff_files[:1]))  # both read volumes, which routing uses
    next(depth_grid.catchment_reader(willamette.runoff_files[:1]))
    assert not depth_grid.as_volumes
    np.testing.assert_array_equal(depth_grid.aggregate(cells, time_index, factor)[0], depths)


def test_catchment_runoff_in_another_river_order_is_refused(willamette: Basin, tmp_path: Path) -> None:
    dates, volumes = willamette.months[0]
    volume_grid = grid(willamette, as_volumes=True)
    reversed_file = tmp_path / 'reversed.nc'
    river_ids, areas = volume_grid.river_ids[::-1], volume_grid.catchment_area[::-1]
    volume_grid.to_netcdf(reversed_file, dates, volumes[::-1], river_ids, areas)
    changes = {'forcing': 'catchment', 'grid_weights_file': None, 'runoff_files': [reversed_file]}
    with pytest.raises(ValueError, match='in the same order'):
        route(tmp_path, **(willamette.gridded(months=1) | changes))


def test_a_weight_table_in_another_river_order_is_refused(willamette: Basin, tmp_path: Path) -> None:
    with xr.open_dataset(willamette.weights_file) as weights:
        position = pd.Index(willamette.network.river_ids).get_indexer(weights['river_id'].to_numpy())
        reordered = weights.isel(index=np.argsort(-position, kind='stable')).load()
    reordered.assign_coords(index=np.arange(reordered.sizes['index'])).to_netcdf(tmp_path / 'reversed.nc')
    with pytest.raises(ValueError, match='in the same order'):
        route(tmp_path, **(willamette.gridded(months=1) | {'grid_weights_file': tmp_path / 'reversed.nc'}))


def test_catchment_runoff_volumes_check_their_layout(willamette: Basin) -> None:
    dates, volumes = willamette.months[0]
    river_ids = willamette.network.river_ids
    runoff = rr.runoff.CatchmentRunoffVolumes(volumes, river_ids)
    runoff.check(river_ids, dates.size)
    with pytest.raises(ValueError, match='has shape'):
        runoff.check(river_ids, dates.size + 1)
    with pytest.raises(ValueError, match='in the same order'):
        runoff.check(river_ids[::-1], dates.size)
    with pytest.raises(ValueError, match='row contiguous'):
        rr.runoff.CatchmentRunoffVolumes(np.asfortranarray(volumes), river_ids).check(river_ids, dates.size)
    first = runoff.first_steps(4)
    assert first.runoff.shape == (river_ids.shape[0], 4) and first.runoff.strides[1] == first.runoff.itemsize
    assert first.river_ids is river_ids


def test_grid_cell_runoff_checks_its_weight_table(willamette: Basin) -> None:
    dates, cell_runoff, _ = next(grid(willamette, as_volumes=True).generator(willamette.runoff_files[:1]))
    river_ids = willamette.network.river_ids
    cell_runoff.check(river_ids, dates.size)
    with pytest.raises(ValueError, match='does not describe'):
        cell_runoff.check(river_ids[:-1], dates.size)
    with pytest.raises(ValueError, match='in the same order'):
        cell_runoff.check(river_ids[::-1], dates.size)
    assert cell_runoff.first_steps(10) is cell_runoff


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
    options = routing_options(willamette, tmp_path) | {'forcing': 'catchment', 'grid_weights_file': None}
    with pytest.raises(TypeError, match='is read by CatchmentRunoff'):
        rr.Router(rr.Configs(**options), runoff=grid(willamette))


def test_an_initial_state_of_the_wrong_size_is_refused(willamette: Basin, tmp_path: Path) -> None:
    pd.DataFrame({'Q': np.zeros(willamette.network.size - 1, dtype=np.float32)}).to_parquet(tmp_path / 'state.parquet')
    with pytest.raises(ValueError, match='channel_state_init_file has'):
        route(
            tmp_path,
            willamette.months[:1],
            network_file=willamette.network_file,
            channel_state_init_file=tmp_path / 'state.parquet',
        )


def test_a_stabilized_run_refuses_sub_reaches_that_change_between_files(willamette: Basin, tmp_path: Path) -> None:
    (march_dates, march), (april_dates, april) = willamette.months[:2]
    n_steps = april.shape[1] // 3 * 3
    three_hourly = april[:, :n_steps].reshape(april.shape[0], -1, 3).sum(axis=2)
    network = willamette.network
    hourly_reaches, three_hourly_reaches = (
        network.conditioning(network.largest_stable_dt(dt))[0].sum() for dt in (DT, 3 * DT)
    )
    assert hourly_reaches != three_hourly_reaches, 'both time steps route the same sub-reaches'
    months = [(march_dates, march), (april_dates[:n_steps:3], three_hourly)]
    with pytest.raises(ValueError, match='sub-reaches'):
        route(tmp_path, months, network_file=willamette.network_file, network_type='stabilized')


def test_null_device_outputs_are_routed_and_discarded(willamette: Basin, tmp_path: Path) -> None:
    options = willamette.gridded(months=1) | {'discharge_dir': os.devnull, 'log': False, 'progress_bar': False}
    written = []
    router = rr.Router(rr.Configs(**options)).set_discharge_writer(lambda *arguments: written.append(arguments))
    router.route()
    assert not written and not any(tmp_path.iterdir())


def test_split_rivers_write_only_the_original_rivers(willamette: Basin, tmp_path: Path) -> None:
    network = rr.Network(willamette.network_file).stabilize(DT)
    routed = route(tmp_path, network=network, writer=writers.zarr_writer, **willamette.gridded(months=1))
    assert routed.discharge[0].shape[0] == willamette.network.size
    with xr.open_zarr(tmp_path / routed.discharge_files[0]) as written:
        np.testing.assert_array_equal(written['river_id'].to_numpy(), willamette.network.river_ids)


def test_the_writer_is_called_once_per_runoff_file(willamette: Basin, tmp_path: Path) -> None:
    routed = route(tmp_path, willamette.months, network_file=willamette.network_file)
    assert routed.discharge_files == [f'discharge_{i}.zarr' for i in range(len(willamette.months))]
    assert [dates.size for dates in routed.dates] == [dates.size for dates, _ in willamette.months]


@pytest.fixture(scope='module')
def march_discharge(willamette: Basin, tmp_path_factory: pytest.TempPathFactory):
    """The basin's routed March discharge, as the router hands it to a writer."""
    return route(tmp_path_factory.mktemp('march'), willamette.months[:1], network_file=willamette.network_file)


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


def test_a_subset_is_a_closed_upstream_basin(willamette: Basin, package: Path, manifest: dict) -> None:
    table = pd.read_parquet(willamette.network_file)
    ids = set(table['riverId'])
    assert set(table['nextRiverId']) <= ids | {-1}
    np.testing.assert_array_equal(table.loc[table['nextRiverId'] == -1, 'riverId'], [WILLAMETTE])
    with xr.open_dataset(willamette.weights_file) as weights:
        assert set(np.unique(weights['river_id'].to_numpy())) == ids
    columbia = rr.Network(package / manifest['runs']['static_standard']['configs']['network_file'])
    in_basin = columbia.river_ids == WILLAMETTE
    for river in range(
        columbia.size - 1, -1, -1
    ):  # downstream rivers come later, so each is decided before its upstreams
        downstream = columbia.downstream_indices[river]
        in_basin[river] |= downstream >= 0 and in_basin[downstream]
    assert in_basin.sum() == len(table), 'every river upstream of the outlet is in the subset'


def test_metrics_of_known_series(march_discharge) -> None:
    series = march_discharge.discharge[0][-1].astype(np.float64)
    assert metrics.kge2012(series, series) == pytest.approx(1)
    assert metrics.me(series, series) == metrics.mae(series, series) == metrics.mse(series, series) == 0
    # adding its mean keeps the correlation at 1, doubles the mean, and halves the coefficient of variation
    assert metrics.kge2012(series, series + series.mean()) == pytest.approx(1 - np.sqrt(1.25))
    assert metrics.me(series, series + 1) == pytest.approx(-1)
    assert metrics.mae(series, series + 1) == metrics.mse(series, series + 1) == pytest.approx(1)


def test_the_command_line_subsets_a_basin(willamette: Basin, package: Path, manifest: dict, tmp_path: Path) -> None:
    configs = manifest['runs']['static_standard']['configs']
    network_file, weights = package / configs['network_file'], package / configs['grid_weights_file']
    out_network_file, out_weights = tmp_path / 'network.parquet', tmp_path / 'gridweights.nc'
    subset = [sys.executable, '-m', 'river_route._cli', 'subset', str(WILLAMETTE)]
    network_options = ['--network', str(network_file), '--out-network', str(out_network_file)]
    weights_options = ['--weights', str(weights), '--out-weights', str(out_weights)]
    result = subprocess.run(subset + network_options + weights_options, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    pd.testing.assert_frame_equal(pd.read_parquet(out_network_file), pd.read_parquet(willamette.network_file))
    with xr.open_dataset(out_weights) as written, xr.open_dataset(willamette.weights_file) as expected:
        xr.testing.assert_identical(written, expected)
    result = subprocess.run(subset + weights_options, capture_output=True, text=True)
    assert result.returncode == 2 and '--network' in result.stderr, result.stderr  # it holds the basin's rows
    with pytest.raises(ValueError, match='together'):
        streams.subset_network_to_river(WILLAMETTE, network_file, out_network_file, weights)


def test_the_command_line_routes_a_config_file(willamette: Basin, tmp_path: Path) -> None:
    configs = rr.Configs(**routing_options(willamette, tmp_path), log=False, progress_bar=False)
    configs.to_json(tmp_path / 'configs.json')
    command = [sys.executable, '-m', 'river_route._cli', 'route', str(tmp_path / 'configs.json')]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    expected = route(tmp_path, **willamette.gridded(months=1)).discharge[0]
    with xr.open_zarr(tmp_path / 'discharge.zarr') as written:
        assert_same(written['Q'].to_numpy(), expected)


@pytest.mark.parametrize('module', MODULES)
def test_every_module_imports_in_a_fresh_process(module: str) -> None:
    result = subprocess.run([sys.executable, '-c', f'import {module}'], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize('module', ['river_route', *MODULES])
def test_every_name_in_all_exists(module: str) -> None:
    imported = importlib.import_module(module)
    missing = [name for name in getattr(imported, '__all__', []) if not hasattr(imported, name)]
    assert not missing, f'{module}.__all__ lists names it does not define: {missing}'


@pytest.mark.parametrize('command', STATIC_CHECKS.values(), ids=STATIC_CHECKS.keys())
def test_static_checks_pass(command: list[str]) -> None:
    result = subprocess.run([str(TOOLS / command[0]), *command[1:]], cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.slow
def test_freshly_compiled_kernels_match_the_reference(tmp_path: Path) -> None:
    """Compile every kernel into an empty numba cache and route the reference solution with them. Freshly compiled and
    cached kernels have differed in the last bit before, for subcycled rivers."""
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
