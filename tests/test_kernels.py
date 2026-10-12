"""
The two routing kernels, routing.c and its numba twin, route the arrays of the reference solution's standard run to its
discharge, all at once and a window of steps per call, and refuse the same broken arrays with the same problem. Build
the C kernel first, with rr compile. The reference runs route with both kernels in test_reference_solution.py.
"""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from reference_solution import GRID_NAMES, assert_same

import river_route as rr
from river_route import kernels
from river_route.router import static_muskingum

KERNELS = ['numba', 'c']
STEPS_PER_WINDOW = 64  # the windows an interactive app routes in


@pytest.fixture(scope='module')
def standard_month(package: Path, manifest: dict) -> dict:
    """The first month of the standard run in the arrays the kernels route, and the discharge the package holds."""
    run = manifest['runs']['static_standard']
    options = run['configs']
    configs = rr.Configs(network_file=package / options['network_file'], invalid_coefficients='ignore')
    network = rr.Network.from_configs(configs)
    grid = rr.GridRunoff(package / options['grid_weights_file'], **GRID_NAMES)
    dates, runoff, _ = next(grid.generator([package / options['forcing_files'][0]]))
    dt_runoff = int((dates[1] - dates[0]) / np.timedelta64(1, 's'))
    coefficients = static_muskingum.prepare_routing(network, configs, dt_runoff, None)
    kernel_network, depth = kernels.kernel_arrays(network, coefficients, runoff, dates.shape[0])
    with xr.open_zarr(package / run['discharge'][0]) as reference:
        expected = reference['Q'].to_numpy()
    return {'network': kernel_network, 'depth': depth[:, : dates.shape[0]], 'expected': expected}


def route_every_river(
    kernel: str, network: kernels.KernelNetwork, depth: np.ndarray, discharge: np.ndarray, **changes
) -> None:
    """Route every river of the network in one job over every step of depth, from a channel state of zero unless
    one is given."""
    n_rivers = network.sub_reach_starts.shape[0] - 1
    arguments = {
        'n_steps': depth.shape[1],
        'block_starts': np.array([0], dtype=np.int32),
        'block_stops': np.array([n_rivers], dtype=np.int32),
        'channel_state': np.zeros(network.c1.shape[0], dtype=np.float32),
    } | changes
    kernels.route(
        kernel,
        network,
        depth,
        arguments['n_steps'],
        arguments['block_starts'],
        arguments['block_stops'],
        arguments['channel_state'],
        discharge,
    )


@pytest.mark.parametrize('kernel', KERNELS)
def test_both_kernels_reproduce_the_standard_run(kernel: str, standard_month: dict) -> None:
    network, depth = standard_month['network'], standard_month['depth']
    discharge = np.empty((network.sub_reach_starts.shape[0] - 1, depth.shape[1] + 1), dtype=np.float32)
    route_every_river(kernel, network, depth, discharge)
    assert_same(discharge[:, 1:], standard_month['expected'])


@pytest.mark.parametrize('kernel', KERNELS)
def test_both_kernels_route_a_run_a_window_at_a_time(kernel: str, standard_month: dict) -> None:
    """
    A caller that holds only a window of steps at a time, as an interactive app does, calls once per window with the
    runoff of the window and the channel state the window before left, and routes the discharge of the whole run.
    """
    network, depth = standard_month['network'], standard_month['depth']
    n_rivers = network.sub_reach_starts.shape[0] - 1
    channel_state = np.zeros(network.c1.shape[0], dtype=np.float32)
    windows = []
    for first_step in range(0, depth.shape[1], STEPS_PER_WINDOW):
        window = depth[:, first_step : first_step + STEPS_PER_WINDOW]
        discharge = np.empty((n_rivers, window.shape[1] + 1), dtype=np.float32)
        route_every_river(kernel, network, window, discharge, channel_state=channel_state)
        windows.append(discharge[:, 1:])
    assert_same(np.concatenate(windows, axis=1), standard_month['expected'])


@pytest.mark.parametrize('kernel', KERNELS)
def test_both_kernels_refuse_steps_that_are_not_one_after_another(kernel: str, standard_month: dict) -> None:
    network, depth = standard_month['network'], standard_month['depth']
    discharge = np.empty((depth.shape[1] + 1, network.sub_reach_starts.shape[0] - 1), dtype=np.float32).T
    with pytest.raises(ValueError, match='must hold each series of steps one after another'):
        route_every_river(kernel, network, depth, discharge)


def drain_into_the_last_river(network: kernels.KernelNetwork, n_rivers: int) -> dict:
    return {'network': replace(network, inflow_rivers=np.full_like(network.inflow_rivers, n_rivers - 1))}


def read_cells_beyond_the_runoff(network: kernels.KernelNetwork, n_rivers: int) -> dict:
    return {'network': replace(network, weight_cells=network.weight_cells + np.int32(1_000_000))}


def give_the_first_river_no_sub_reaches(network: kernels.KernelNetwork, n_rivers: int) -> dict:
    sub_reach_starts = network.sub_reach_starts.copy()
    sub_reach_starts[1] = sub_reach_starts[0]
    return {'network': replace(network, sub_reach_starts=sub_reach_starts)}


def route_the_blocks_out_of_order(network: kernels.KernelNetwork, n_rivers: int) -> dict:
    middle = n_rivers // 2
    return {'block_starts': np.array([middle, 0], np.int32), 'block_stops': np.array([n_rivers, middle], np.int32)}


def route_more_steps_than_the_runoff_holds(network: kernels.KernelNetwork, n_rivers: int) -> dict:
    return {'n_steps': 10**9}


@pytest.mark.parametrize(
    ('breaks', 'problem'),
    [
        (drain_into_the_last_river, 'a river drains into a river before it in routing order'),
        (read_cells_beyond_the_runoff, 'a river has weights out of order or a weight reads a cell outside'),
        (give_the_first_river_no_sub_reaches, 'a river has no sub-reaches'),
        (route_the_blocks_out_of_order, 'the blocks of rivers are not in routing order'),
        (route_more_steps_than_the_runoff_holds, 'a count of rivers, sub-reaches, cells, steps, or blocks'),
    ],
)
@pytest.mark.parametrize('kernel', KERNELS)
def test_both_kernels_refuse_broken_arrays(kernel: str, breaks, problem: str, standard_month: dict) -> None:
    network, depth = standard_month['network'], standard_month['depth']
    n_rivers = network.sub_reach_starts.shape[0] - 1
    changes = breaks(network, n_rivers)
    discharge = np.empty((n_rivers, depth.shape[1] + 1), dtype=np.float32)
    with pytest.raises(ValueError, match=f'the {kernel} kernel refused to route: {problem}'):
        route_every_river(kernel, changes.pop('network', network), depth, discharge, **changes)
