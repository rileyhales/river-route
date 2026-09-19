"""
Configs(network_conditioning='stabilized') routes every river too long for dt_routing as the fewest equal sub-reaches in
series that are each Muskingum-stable, and sub-cycles every river too short for it in the fewest equal steps that are.
Splitting must be exactly a standard network in which each sub-reach is its own river with k/N and 1/N of the river's
lateral inflow, and sub-cycling must be exactly routing that river at dt/m, which is what these tests compare against.
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from conftest import write_vlateral

import river_route as rr

N_RIVERS = 300
N_STEPS = 48
DT = 3600


def dfs_tree(rng: np.random.Generator, n_rivers: int, n_outlets: int) -> np.ndarray:
    """Random branching network numbered in DFS post-order, so every subtree is a contiguous block."""
    parent = np.full(n_rivers, -1)
    for node in range(n_outlets, n_rivers):
        parent[node] = rng.integers(0, node)
    children = [[] for _ in range(n_rivers)]
    for node in range(n_outlets, n_rivers):
        children[parent[node]].append(node)
    position = np.empty(n_rivers, dtype=np.int64)
    counter = 0
    for root in range(n_outlets):
        stack = [(root, False)]
        while stack:
            node, done = stack.pop()
            if done:
                position[node] = counter
                counter += 1
                continue
            stack.append((node, True))
            stack.extend((child, False) for child in reversed(children[node]))
    downstream = np.full(n_rivers, -1, dtype=np.int64)
    has_parent = parent >= 0
    downstream[position[has_parent]] = position[parent[has_parent]]
    return downstream


@pytest.fixture
def case(tmp_path: Path) -> dict:
    """A network where most rivers are too long for dt, a few are too short, with lateral inflow and a state."""
    rng = np.random.default_rng(5)
    downstream = dfs_tree(rng, N_RIVERS, n_outlets=2)
    river_ids = np.arange(1, N_RIVERS + 1, dtype=np.int64)
    k = rng.uniform(0.5, 12, N_RIVERS) * DT
    x = np.full(N_RIVERS, 0.2)
    params = pd.DataFrame(
        dict(river_id=river_ids, next_river_id=np.where(downstream >= 0, river_ids[downstream], -1), k=k, x=x)
    )
    params_file = tmp_path / 'params.parquet'
    params.to_parquet(params_file, index=False)
    state = rng.uniform(1, 20, N_RIVERS)
    state_file = tmp_path / 'state.parquet'
    pd.DataFrame({'Q': state}).to_parquet(state_file, index=False)
    vlateral = rng.uniform(0, 5e4, (N_STEPS, N_RIVERS)).astype(np.float32)
    vlateral_file = tmp_path / 'vlateral.nc'
    write_vlateral(vlateral_file, vlateral, river_ids, dt=DT)
    return dict(
        tmp_path=tmp_path,
        params_file=params_file,
        state_file=state_file,
        vlateral_file=vlateral_file,
        vlateral=vlateral,
        state=state,
        river_ids=river_ids,
    )


def configs(case: dict, out: Path, **kwargs) -> rr.Configs:
    options = dict(
        params_file=case['params_file'],
        vlateral_files=[case['vlateral_file']],
        channel_state_init_file=case['state_file'],
        discharge_files=[out],
        forcing='vlateral',
        routing_order='river',
        dt_routing=DT,
        unstable_coefficients='ignore',
        log=False,
        progress_bar=False,
    )
    return rr.Configs(**{**options, **kwargs})


def sub_reach_case(case: dict) -> dict:
    """The stabilized network written as a standard one: each sub-reach a river, with its share of the inflow."""
    stable = rr.Network(case['params_file']).stabilize(DT, mode='uniform')
    directory = case['tmp_path'] / 'sub_reaches'
    directory.mkdir()
    reach_ids = np.arange(1, stable.n_reaches + 1, dtype=np.int64)
    params_file = directory / 'params.parquet'
    pd.DataFrame(
        dict(
            river_id=reach_ids,
            next_river_id=np.where(stable.downstream_index >= 0, reach_ids[stable.downstream_index], -1),
            k=stable.k,
            x=stable.x,
        )
    ).to_parquet(params_file, index=False)
    vlateral_file = directory / 'vlateral.nc'
    write_vlateral(vlateral_file, case['vlateral'][:, stable.parent_index] * stable.lateral_scale, reach_ids, dt=DT)
    state_file = directory / 'state.parquet'
    pd.DataFrame({'Q': stable.broadcast_state(case['state'])}).to_parquet(state_file, index=False)
    return dict(stable=stable, params_file=params_file, vlateral_file=vlateral_file, state_file=state_file)


def discharge(path: Path) -> np.ndarray:
    with xr.open_dataset(path) as ds:
        return ds['Q'].values


def test_stabilized_matches_sub_reaches_routed_as_rivers(case):
    params = pd.read_parquet(case['params_file'])
    params['k'] = params['k'].clip(lower=0.7 * DT)  # no river too short, so none is sub-cycled
    params.to_parquet(case['params_file'], index=False)
    reference = sub_reach_case(case)
    stable = reference['stable']
    assert stable.reaches_added > N_RIVERS, 'most rivers should be split'

    final_expected = case['tmp_path'] / 'final_expected.parquet'
    rr.Router(
        configs({**case, **reference}, case['tmp_path'] / 'q_expected.nc', channel_state_final_file=final_expected)
    ).route()
    final = case['tmp_path'] / 'final.parquet'
    router = rr.Router(
        configs(case, case['tmp_path'] / 'q.nc', network_conditioning='stabilized', channel_state_final_file=final)
    )
    router.route()

    np.testing.assert_array_equal(router.subdivisions, stable.subdivisions)
    expected = discharge(case['tmp_path'] / 'q_expected.nc')[:, stable.outlet_index]
    got = discharge(case['tmp_path'] / 'q.nc')
    np.testing.assert_allclose(got, expected, rtol=1e-5, atol=1e-5 * float(np.abs(expected).max()))
    # the final state has one row per sub-reach, in the stabilized network's order
    np.testing.assert_allclose(
        pd.read_parquet(final)['Q'].to_numpy(), pd.read_parquet(final_expected)['Q'].to_numpy(), rtol=1e-5, atol=1e-3
    )


def test_stabilized_restarts_from_its_own_final_state(case):
    """Routing two halves through a per-sub-reach state file gives what one run over both gives."""
    whole = rr.Router(configs(case, case['tmp_path'] / 'q_whole.nc', network_conditioning='stabilized'))
    whole.route()

    half = N_STEPS // 2
    first_file = case['tmp_path'] / 'first.nc'
    second_file = case['tmp_path'] / 'second.nc'
    write_vlateral(first_file, case['vlateral'][:half], case['river_ids'], dt=DT)
    write_vlateral(second_file, case['vlateral'][half:], case['river_ids'], dt=DT)
    middle = case['tmp_path'] / 'middle.parquet'
    rr.Router(
        configs(
            case,
            case['tmp_path'] / 'q_first.nc',
            vlateral_files=[first_file],
            network_conditioning='stabilized',
            channel_state_final_file=middle,
        )
    ).route()
    assert pd.read_parquet(middle).shape[0] == int(whole.reach_indptr[-1])
    rr.Router(
        configs(
            case,
            case['tmp_path'] / 'q_second.nc',
            vlateral_files=[second_file],
            channel_state_init_file=middle,
            network_conditioning='stabilized',
        )
    ).route()
    joined = np.concatenate([discharge(case['tmp_path'] / 'q_first.nc'), discharge(case['tmp_path'] / 'q_second.nc')])
    expected = discharge(case['tmp_path'] / 'q_whole.nc')
    np.testing.assert_allclose(joined, expected, rtol=1e-5, atol=1e-5 * float(np.abs(expected).max()))


def test_stabilized_with_a_thread_pool_matches_single_threaded(case):
    single = rr.Router(configs(case, case['tmp_path'] / 'q_single.nc', network_conditioning='stabilized'))
    single.route()
    router = rr.Router(configs(case, case['tmp_path'] / 'q_threaded.nc', network_conditioning='stabilized'))
    with ThreadPoolExecutor(4) as pool:
        router.route(thread_pool=pool, threads=4)
    assert len(router.routing_jobs) > 2, 'the network should split into concurrent regions'
    expected = discharge(case['tmp_path'] / 'q_single.nc')
    np.testing.assert_allclose(
        discharge(case['tmp_path'] / 'q_threaded.nc'), expected, rtol=1e-5, atol=1e-6 * float(np.abs(expected).max())
    )


def test_standard_network_is_unchanged_by_the_option(case):
    """A network where every river is already stable routes the same either way."""
    params = pd.read_parquet(case['params_file'])
    params['k'] = DT * 2.0  # 2*k*x = 0.8 dt and 2*k*(1-x) = 3.2 dt, stable without splitting
    params.to_parquet(case['params_file'], index=False)
    standard = discharge_of(case, 'q_standard.nc')
    stabilized = discharge_of(case, 'q_stabilized.nc', network_conditioning='stabilized')
    np.testing.assert_array_equal(stabilized, standard)


def discharge_of(case: dict, name: str, **kwargs) -> np.ndarray:
    rr.Router(configs(case, case['tmp_path'] / name, **kwargs)).route()
    return discharge(case['tmp_path'] / name)


def test_stabilized_picks_the_largest_stable_dt(case):
    network = rr.Network(case['params_file'])
    assert network.largest_stable_dt(DT) == DT  # x = 0.2, so every river that is not too short can be split
    router = rr.Router(configs(case, case['tmp_path'] / 'q_auto.nc', network_conditioning='stabilized', dt_routing=0))
    router.route()
    assert router.dt_routing == DT


def test_largest_stable_dt_steps_down_when_a_window_holds_no_whole_count(tmp_path):
    """At x = 0.45 the window [0.9 k/dt, 1.1 k/dt] of sub-reach counts can miss every integer."""
    params_file = tmp_path / 'params.parquet'
    pd.DataFrame(dict(river_id=[1], next_river_id=[-1], k=[1.5 * DT], x=[0.45])).to_parquet(params_file, index=False)
    network = rr.Network(params_file)
    dt = network.largest_stable_dt(DT)
    assert dt < DT and DT % dt == 0
    subdivisions, resolvable = network.subdivisions(dt)
    assert resolvable.all()
    too_long, too_short = network.unstable_mask(dt, subdivisions)
    assert not too_long.any() and not too_short.any()


@pytest.mark.parametrize('routing_order,coeff', [('time', 'static'), ('river', 'dynamic')])
def test_stabilized_needs_river_order_and_static_coefficients(case, routing_order, coeff):
    params = pd.read_parquet(case['params_file'])
    params['alpha'] = 1.0
    params['beta'] = 0.0
    params.to_parquet(case['params_file'], index=False)
    router = rr.Router(
        configs(
            case, case['tmp_path'] / 'q.nc', network_conditioning='stabilized', routing_order=routing_order, coeff=coeff
        )
    )
    with pytest.raises(NotImplementedError, match='network_conditioning=stabilized'):
        router.route()


def test_network_conditioning_is_validated():
    with pytest.raises(ValueError, match='network_conditioning must be one of'):
        rr.Configs(network_conditioning='expanded')


def write_params(path: Path, river_ids, next_river_ids, k, x) -> Path:
    pd.DataFrame(dict(river_id=river_ids, next_river_id=next_river_ids, k=k, x=x)).to_parquet(path, index=False)
    return path


def test_every_river_is_stable_once_conditioned(case):
    router = rr.Router(configs(case, case['tmp_path'] / 'q.nc', network_conditioning='stabilized'))
    router.route()
    assert router.substeps.shape[0] and router.substeps.max() > 1, 'some rivers should be sub-cycled'
    for name in ('c1', 'c2', 'c3'):
        assert getattr(router, name).min() >= -1e-6, f'{name} is negative for a conditioned river'
    network = router.network
    subdivisions, substeps, resolvable = network.conditioning(DT)
    assert resolvable.all()
    assert not np.any((subdivisions > 1) & (substeps > 1)), 'no river is both split and sub-cycled'
    too_long, too_short = network.unstable_mask(DT, subdivisions, substeps)
    assert not too_long.any() and not too_short.any()


def test_sub_cycling_isolated_rivers_matches_routing_them_at_the_shorter_step(tmp_path):
    """Rivers with no upstream that each need 4 steps route exactly as they would at dt_routing = dt / 4."""
    rng = np.random.default_rng(3)
    n = 20
    river_ids = np.arange(1, n + 1)
    k = rng.uniform(570, 740, n)  # ceil(3600 / (1.6 k)) = 4
    params_file = write_params(tmp_path / 'params.parquet', river_ids, np.full(n, -1), k, np.full(n, 0.2))
    vlateral_file = tmp_path / 'vlateral.nc'
    write_vlateral(vlateral_file, rng.uniform(0, 5e4, (N_STEPS, n)), river_ids, dt=DT)
    state_file = tmp_path / 'state.parquet'
    pd.DataFrame({'Q': rng.uniform(1, 20, n)}).to_parquet(state_file, index=False)
    case = dict(params_file=params_file, vlateral_file=vlateral_file, state_file=state_file)

    router = rr.Router(configs(case, tmp_path / 'q_cycled.nc', network_conditioning='stabilized'))
    router.route()
    np.testing.assert_array_equal(router.substeps, 4)
    rr.Router(configs(case, tmp_path / 'q_fine.nc', dt_routing=DT // 4)).route()
    expected = discharge(tmp_path / 'q_fine.nc')
    np.testing.assert_allclose(
        discharge(tmp_path / 'q_cycled.nc'), expected, rtol=1e-5, atol=1e-6 * float(np.abs(expected).max())
    )


def test_sub_cycled_river_downstream_of_a_smooth_inflow_tracks_the_fine_step(tmp_path):
    """
    Upstream inflow reaches a sub-cycled river once per routing step and is interpolated between those levels. With
    a smooth inflow that is close to routing everything at the shorter step.
    """
    river_ids = np.array([1, 2])
    params_file = write_params(
        tmp_path / 'params.parquet', river_ids, np.array([2, -1]), np.array([2 * DT, 650.0]), np.array([0.0, 0.2])
    )
    hours = np.arange(N_STEPS * 4)
    volumes = np.zeros((hours.size, 2))
    volumes[:, 0] = 3e5 * (1 + np.sin(hours / 10.0))  # a smooth inflow at the upstream river only
    vlateral_file = tmp_path / 'vlateral.nc'
    write_vlateral(vlateral_file, volumes, river_ids, dt=DT)
    case = dict(params_file=params_file, vlateral_file=vlateral_file, state_file=None)

    rr.Router(
        configs(case, tmp_path / 'q_cycled.nc', network_conditioning='stabilized', channel_state_init_file=None)
    ).route()
    rr.Router(configs(case, tmp_path / 'q_fine.nc', dt_routing=DT // 4, channel_state_init_file=None)).route()
    cycled = discharge(tmp_path / 'q_cycled.nc')[:, 1]
    fine = discharge(tmp_path / 'q_fine.nc')[:, 1]
    # the upstream river itself is routed at one hour in the first run and fifteen minutes in the second
    assert np.abs(cycled - fine).max() < 0.01 * fine.max()
