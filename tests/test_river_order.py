"""
River order kernels (Configs(routing_order='river')) route each river's whole series before the next river. They must
give the time order kernels' results to float32 round-off on every combination they implement.
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from conftest import write_vlateral

import river_route as rr
from river_route.router import _kernel_registry, _numba_kernels, _river_kernels, writers

N_RIVERS = 400
N_STEPS = 60
DT_RUNOFF = 3600


def dfs_ordered_tree(rng: np.random.Generator, n_rivers: int, n_outlets: int) -> np.ndarray:
    """Random tree whose downstream index is always above the river's own, mostly the next river as in DFS order."""
    downstream = np.full(n_rivers, -1, dtype=np.int32)
    for i in range(n_rivers - n_outlets):
        downstream[i] = min(n_rivers - 1, i + int(rng.geometric(0.6)))
    return downstream


def dfs_tree(rng: np.random.Generator, n_rivers: int, n_outlets: int) -> np.ndarray:
    """
    Random branching network numbered in DFS post-order, so every subtree is a contiguous block ending at its outlet,
    which is what concurrent routing requires.
    """
    parent = np.full(n_rivers, -1)
    for node in range(n_outlets, n_rivers):
        parent[node] = rng.integers(0, node)  # grown from the outlets toward the headwaters
    children = [[] for _ in range(n_rivers)]
    for node in range(n_outlets, n_rivers):
        children[parent[node]].append(node)
    position = np.empty(n_rivers, dtype=np.int64)
    counter = 0
    for root in range(n_outlets):
        stack = [(root, False)]
        while stack:
            node, expanded = stack.pop()
            if expanded:
                position[node] = counter
                counter += 1
                continue
            stack.append((node, True))
            stack.extend((child, False) for child in reversed(children[node]))
    downstream = np.full(n_rivers, -1, dtype=np.int32)
    has_parent = parent >= 0
    downstream[position[has_parent]] = position[parent[has_parent]]
    return downstream


def coefficients(rng: np.random.Generator, n_rivers: int, dt_routing: int):
    k = rng.uniform(dt_routing * 0.6, dt_routing * 10, n_rivers)
    x = rng.uniform(0.0, 0.4, n_rivers)
    dt_div_k = dt_routing / k
    denominator = dt_div_k + 2 * (1 - x)
    c1 = ((dt_div_k - 2 * x) / denominator).astype(np.float32)
    c2 = ((dt_div_k + 2 * x) / denominator).astype(np.float32)
    c3 = ((2 * (1 - x) - dt_div_k) / denominator).astype(np.float32)
    return k, x, c1, c2, c3


def one_pass(n_rivers: int, n_routing_steps: int, buffers: int = 1) -> dict:
    """Keywords for one unthreaded pass of a time order kernel: one block, no outlet, no cuts."""
    boundary = [np.zeros((1, n_routing_steps), np.float32) for _ in range(buffers)]
    names = ['boundary'] if buffers == 1 else ['boundary_old', 'boundary_new']
    return dict(
        block_starts=np.array([0], np.int32),
        block_stops=np.array([n_rivers], np.int32),
        outlet=-1,
        region=0,
        cut_target=np.zeros(0, np.int32),
        **dict(zip(names, boundary, strict=True)),
    )


def assert_same_routing(q_river, d_river, q_time, d_time):
    scale = float(np.abs(d_time).max())
    np.testing.assert_allclose(d_river, d_time, rtol=1e-5, atol=1e-6 * scale)
    np.testing.assert_allclose(q_river, q_time, rtol=1e-5, atol=1e-6 * scale)


# ── kernels ─────────────────────────────────────────────────────────────────


@pytest.mark.parametrize('n_substeps', [1, 3])
def test_static_channel_matches_time_order(n_substeps):
    rng = np.random.default_rng(1)
    downstream = dfs_ordered_tree(rng, N_RIVERS, n_outlets=3)
    _, _, c1, c2, c3 = coefficients(rng, N_RIVERS, DT_RUNOFF // n_substeps)
    has_downstream = downstream >= 0
    q_init = rng.uniform(1, 100, N_RIVERS).astype(np.float32)

    q_time, d_time = q_init.copy(), np.empty((N_STEPS, N_RIVERS), np.float32)
    _numba_kernels.static_channel(
        q_t=q_time,
        discharge_array=d_time,
        downstream_indices=downstream,
        downstream_c1=np.where(has_downstream, c1[downstream], 0).astype(np.float32),
        downstream_c2=np.where(has_downstream, c2[downstream], 0).astype(np.float32),
        c3=c3,
        n_rivers=N_RIVERS,
        n_steps=N_STEPS,
        n_substeps=n_substeps,
        **one_pass(N_RIVERS, N_STEPS * n_substeps),
    )
    q_river, d_river = q_init.copy(), np.empty((N_STEPS, N_RIVERS), np.float32)
    _river_kernels.static_channel(
        q_t=q_river,
        discharge_array=d_river,
        downstream_indices=downstream,
        c1=c1,
        c2=c2,
        c3=c3,
        n_substeps=n_substeps,
        block=_river_kernels.BLOCK,
    )
    assert_same_routing(q_river, d_river, q_time, d_time)


@pytest.mark.parametrize('by_river', [False, True])
@pytest.mark.parametrize('n_substeps', [1, 3])
def test_static_vlateral_matches_time_order(n_substeps, by_river):
    rng = np.random.default_rng(2)
    downstream = dfs_ordered_tree(rng, N_RIVERS, n_outlets=3)
    _, _, c1, c2, c3 = coefficients(rng, N_RIVERS, DT_RUNOFF // n_substeps)
    c4_dt = ((c1 + c2) / DT_RUNOFF).astype(np.float32)
    has_downstream = downstream >= 0
    q_init = rng.uniform(0, 50, N_RIVERS).astype(np.float32)
    vlateral = rng.uniform(0, 5e4, (N_STEPS, N_RIVERS)).astype(np.float32)

    q_time, d_time = q_init.copy(), np.empty((N_STEPS, N_RIVERS), np.float32)
    _numba_kernels.static_vlateral(
        q_t=q_time,
        discharge_array=d_time,
        downstream_indices=downstream,
        downstream_c1=np.where(has_downstream, c1[downstream], 0).astype(np.float32),
        downstream_c2=np.where(has_downstream, c2[downstream], 0).astype(np.float32),
        c3=c3,
        n_rivers=N_RIVERS,
        n_steps=N_STEPS,
        n_substeps=n_substeps,
        vlateral=vlateral,
        c4_dt=c4_dt,
        **one_pass(N_RIVERS, N_STEPS * n_substeps),
    )
    q_river, d_river = q_init.copy(), np.empty((N_STEPS, N_RIVERS), np.float32)
    _river_kernels.static_vlateral(
        q_t=q_river,
        discharge_array=d_river,
        downstream_indices=downstream,
        c1=c1,
        c2=c2,
        c3=c3,
        c4_dt=c4_dt,
        n_substeps=n_substeps,
        vlateral=np.ascontiguousarray(vlateral.T) if by_river else vlateral,
        by_river=by_river,
        block=_river_kernels.BLOCK,
    )
    assert_same_routing(q_river, d_river, q_time, d_time)


@pytest.mark.parametrize('by_river', [False, True])
@pytest.mark.parametrize('n_substeps', [1, 3])
def test_dynamic_vlateral_matches_time_order(n_substeps, by_river):
    rng = np.random.default_rng(3)
    dt_routing = DT_RUNOFF // n_substeps
    downstream = dfs_ordered_tree(rng, N_RIVERS, n_outlets=3)
    _, x, _, _, _ = coefficients(rng, N_RIVERS, dt_routing)
    x = x.astype(np.float32)
    # K = alpha q^beta stays between about 1 and 10 dt over the discharges this network reaches, so it is stable
    beta = rng.uniform(-0.1, -0.05, N_RIVERS).astype(np.float32)
    alpha = (rng.uniform(2, 8, N_RIVERS) * dt_routing * 1e3**-beta).astype(np.float32)
    q_init = rng.uniform(0, 50, N_RIVERS).astype(np.float32)
    vlateral = rng.uniform(0, 5e4, (N_STEPS, N_RIVERS)).astype(np.float32)
    shared = dict(
        downstream_indices=downstream,
        alpha=alpha,
        beta=beta,
        x=x,
        dt_routing=np.float32(dt_routing),
        dt_runoff=np.float32(DT_RUNOFF),
        n_substeps=n_substeps,
    )

    q_time, d_time = q_init.copy(), np.empty((N_STEPS, N_RIVERS), np.float32)
    _numba_kernels.dynamic_vlateral(
        q_t=q_time,
        discharge_array=d_time,
        n_rivers=N_RIVERS,
        n_steps=N_STEPS,
        vlateral=vlateral,
        **shared,
        **one_pass(N_RIVERS, N_STEPS * n_substeps, buffers=2),
    )
    q_river, d_river = q_init.copy(), np.empty((N_STEPS, N_RIVERS), np.float32)
    _river_kernels.dynamic_vlateral(
        q_t=q_river,
        discharge_array=d_river,
        vlateral=np.ascontiguousarray(vlateral.T) if by_river else vlateral,
        by_river=by_river,
        block=_river_kernels.BLOCK,
        **shared,
    )
    # the time order kernel builds its coefficients partly in float64, so the nonlinear result drifts a little more
    scale = float(np.abs(d_time).max())
    np.testing.assert_allclose(d_river, d_time, rtol=1e-4, atol=1e-5 * scale)
    np.testing.assert_allclose(q_river, q_time, rtol=1e-4, atol=1e-5 * scale)


def test_layout_is_read_from_the_array():
    """A reader that builds (river, time) arrays hands over their transpose, which is read one river row at a time."""
    by_river = np.zeros((7, 5), np.float32)  # (river, time)
    array, read_by_river = _kernel_registry._river_layout(by_river.T)
    assert read_by_river and np.shares_memory(array, by_river)
    array, read_by_river = _kernel_registry._river_layout(np.zeros((5, 7), np.float32))
    assert not read_by_river and array.flags.c_contiguous
    # a truncated transpose is neither layout, so it is copied to (time, river)
    array, read_by_river = _kernel_registry._river_layout(by_river.T[:3])
    assert not read_by_river and array.shape == (3, 7) and array.flags.c_contiguous


def test_registry_keys_include_routing_order():
    for order in ('time', 'river'):
        for coeff, forcing in (('static', 'channel'), ('static', 'vlateral'), ('dynamic', 'vlateral')):
            assert _kernel_registry.resolve_kernel(coeff, forcing, 'uniform', 'standard', order).routing_order == order
    with pytest.raises(NotImplementedError, match='routing_order=river'):
        _kernel_registry.resolve_kernel('static', 'vlateral', 'uniform', 'expanded', 'river')


# ── through the Router ──────────────────────────────────────────────────────


def write_tree_case(directory: Path, dynamic: bool = False) -> dict:
    """A branching DFS ordered network, a random lateral inflow file, and a nonzero initial state."""
    rng = np.random.default_rng(11)
    downstream = dfs_tree(rng, N_RIVERS, n_outlets=3)
    river_ids = np.arange(1, N_RIVERS + 1, dtype=np.int64)
    k, x, _, _, _ = coefficients(rng, N_RIVERS, DT_RUNOFF)
    columns = dict(river_id=river_ids, next_river_id=np.where(downstream >= 0, river_ids[downstream], -1), k=k, x=x)
    if dynamic:
        columns['beta'] = rng.uniform(-0.1, -0.05, N_RIVERS)
        columns['alpha'] = rng.uniform(2, 8, N_RIVERS) * DT_RUNOFF * 1e3 ** -columns['beta']
    directory.mkdir(parents=True, exist_ok=True)
    params_file = directory / 'params.parquet'
    pd.DataFrame(columns).to_parquet(params_file, index=False)
    state_file = directory / 'state.parquet'
    pd.DataFrame({'Q': rng.uniform(1, 20, N_RIVERS)}).to_parquet(state_file, index=False)
    vlateral_file = directory / 'vlateral.nc'
    write_vlateral(vlateral_file, rng.uniform(0, 5e4, (N_STEPS, N_RIVERS)), river_ids, dt=DT_RUNOFF)
    return dict(params_file=str(params_file), channel_state_init_file=str(state_file), vlateral_file=vlateral_file)


def route(case: dict, out: Path, routing_order: str, forcing: str = 'vlateral', **kwargs) -> np.ndarray:
    options = dict(params_file=case['params_file'], channel_state_init_file=case['channel_state_init_file'])
    if forcing == 'vlateral':
        options['vlateral_files'] = [str(case['vlateral_file'])]
    else:
        options.update(dt_routing=DT_RUNOFF, dt_total=DT_RUNOFF * N_STEPS)
    configs = rr.Configs(
        forcing=forcing,
        routing_order=routing_order,
        discharge_files=[str(out)],
        log=False,
        progress_bar=False,
        **options,
        **kwargs,
    )
    rr.Router(configs).route()
    with xr.open_dataset(out) as ds:
        return ds['Q'].values


@pytest.mark.parametrize('coeff,forcing', [('static', 'vlateral'), ('dynamic', 'vlateral'), ('static', 'channel')])
def test_router_river_order_matches_time_order(tmp_path, coeff, forcing):
    case = write_tree_case(tmp_path, dynamic=coeff == 'dynamic')
    time_order = route(case, tmp_path / 'q_time.nc', 'time', forcing=forcing, coeff=coeff)
    river_order = route(case, tmp_path / 'q_river.nc', 'river', forcing=forcing, coeff=coeff)
    rtol = 1e-4 if coeff == 'dynamic' else 1e-5
    np.testing.assert_allclose(river_order, time_order, rtol=rtol, atol=rtol * float(np.abs(time_order).max()))


@pytest.mark.parametrize('coeff,forcing', [('static', 'vlateral'), ('dynamic', 'vlateral'), ('static', 'channel')])
def test_river_order_with_a_thread_pool_matches_single_threaded(tmp_path, coeff, forcing):
    """Regions routed concurrently hand their outlet series to the main stem, which must reproduce one pass."""
    case = write_tree_case(tmp_path, dynamic=coeff == 'dynamic')
    single = route(case, tmp_path / 'q_single.nc', 'river', forcing=forcing, coeff=coeff)
    options = dict(params_file=case['params_file'], channel_state_init_file=case['channel_state_init_file'])
    if forcing == 'vlateral':
        options['vlateral_files'] = [str(case['vlateral_file'])]
    else:
        options.update(dt_routing=DT_RUNOFF, dt_total=DT_RUNOFF * N_STEPS)
    out = tmp_path / 'q_threaded.nc'
    configs = rr.Configs(
        forcing=forcing,
        coeff=coeff,
        routing_order='river',
        discharge_files=[str(out)],
        log=False,
        progress_bar=False,
        **options,
    )
    router = rr.Router(configs)
    with ThreadPoolExecutor(4) as pool:
        router.route(thread_pool=pool, threads=4)
    assert len(router.routing_jobs) > 2, 'the network should split into concurrent regions'
    # a confluence adds a region's series after its main stem upstreams, which differs from one pass only in rounding
    with xr.open_dataset(out) as ds:
        np.testing.assert_allclose(ds['Q'].values, single, rtol=1e-5, atol=1e-6 * float(np.abs(single).max()))


def test_routing_order_is_validated():
    with pytest.raises(ValueError, match='routing_order must be one of'):
        rr.Configs(routing_order='basin')


def test_river_order_writes_in_the_layout_the_writer_declares(tmp_path):
    """A writer that declares discharge_layout = 'river' gets the transpose of a (river, time) buffer, same values."""
    case = write_tree_case(tmp_path)
    received = {}

    def capture(name):
        def writer(router, dates, discharge_array, discharge_file, runoff_file=''):
            received[name] = discharge_array.copy(order='K')

        return writer

    by_river = capture('river')
    by_river.discharge_layout = 'river'
    configs = rr.Configs(
        forcing='vlateral',
        routing_order='river',
        params_file=case['params_file'],
        channel_state_init_file=case['channel_state_init_file'],
        vlateral_files=[str(case['vlateral_file'])],
        discharge_files=[str(tmp_path / 'unused.nc')],
        log=False,
        progress_bar=False,
    )
    rr.Router(configs).set_discharge_writer(capture('time')).route()
    rr.Router(configs).set_discharge_writer(by_river).route()
    assert received['time'].flags.c_contiguous
    assert received['river'].T.flags.c_contiguous and not received['river'].flags.c_contiguous
    np.testing.assert_array_equal(received['river'], received['time'])
    assert writers.zarr_writer.discharge_layout == 'river'
    assert getattr(writers.netcdf_writer, 'discharge_layout', 'time') == 'time'
