"""
The routing kernels solve each river's whole series before moving to the next river, eight steps per serial
operation. They must agree with the plainest reading of the Muskingum math: the whole network stepped one routing
step at a time, written here as an independent reference so the kernels are checked against something other than
themselves.
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from conftest import write_vlateral

import river_route as rr
from river_route.router import _kernel_registry, _river_kernels, writers

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


def reference_static(downstream, c1, c2, c3, c4_dt, q_init, vlateral, n_steps, n_substeps):
    """
    Muskingum stepped one routing step at a time over the whole network, in float64.

    Q_i(g) = c3_i Q_i(g-1) + c4dt_i vlateral_i(t) + c1_i U_i(g) + c2_i U_i(g-1), where U_i sums the upstream
    discharges. Every upstream index is below its downstream index, so U_i is complete by the time river i is
    solved. Returns the final state and the per runoff step mean, clamped at zero, in (river, time) as the kernels
    write it.
    """
    n_rivers = downstream.shape[0]
    q = q_init.astype(np.float64)
    upstream_prev = np.zeros(n_rivers)
    for i in range(n_rivers):
        if downstream[i] >= 0:
            upstream_prev[downstream[i]] += q[i]
    out = np.zeros((n_steps, n_rivers))
    for t in range(n_steps):
        interval_sum = np.zeros(n_rivers)
        for _ in range(n_substeps):
            upstream = np.zeros(n_rivers)
            q_next = np.empty(n_rivers)
            for i in range(n_rivers):
                lateral = 0.0 if vlateral is None else c4_dt[i] * vlateral[t, i]
                q_next[i] = c3[i] * q[i] + lateral + c1[i] * upstream[i] + c2[i] * upstream_prev[i]
                if downstream[i] >= 0:
                    upstream[downstream[i]] += q_next[i]
            q, upstream_prev = q_next, upstream
            interval_sum += q
        out[t] = np.maximum(interval_sum / n_substeps, 0.0)
    return q, out.T


def reference_dynamic(downstream, alpha, beta, x, dt_routing, dt_runoff, q_init, vlateral, n_steps, n_substeps):
    """The same walk with the coefficients rebuilt from each river's own discharge every substep, K = alpha q^beta.
    Returns the per runoff step mean in (river, time), as the kernels write it."""
    n_rivers = downstream.shape[0]
    q = q_init.astype(np.float64)
    upstream_prev = np.zeros(n_rivers)
    for i in range(n_rivers):
        if downstream[i] >= 0:
            upstream_prev[downstream[i]] += q[i]
    out = np.zeros((n_steps, n_rivers))
    for t in range(n_steps):
        interval_sum = np.zeros(n_rivers)
        for _ in range(n_substeps):
            upstream = np.zeros(n_rivers)
            q_next = np.empty(n_rivers)
            for i in range(n_rivers):
                k = alpha[i] * max(q[i], 1e-6) ** beta[i]
                dt_div_k = dt_routing / k
                denominator = dt_div_k + 2 * (1 - x[i])
                c1 = (dt_div_k - 2 * x[i]) / denominator
                c2 = (dt_div_k + 2 * x[i]) / denominator
                c3 = (2 * (1 - x[i]) - dt_div_k) / denominator
                external = vlateral[t, i] / dt_runoff
                q_next[i] = c3 * q[i] + ((c1 + c2) * external + c1 * upstream[i] + c2 * upstream_prev[i])
                if downstream[i] >= 0:
                    upstream[downstream[i]] += q_next[i]
            q, upstream_prev = q_next, upstream
            interval_sum += q
        out[t] = np.maximum(interval_sum / n_substeps, 0.0)
    return q, out.T


def river_major(n_steps: int, n_rivers: int) -> np.ndarray:
    """An empty discharge array in the only layout the kernels accept: C-order (river, time)."""
    return np.empty((n_rivers, n_steps), np.float32)


def assert_matches_reference(q_kernel, d_kernel, q_reference, d_reference, rtol=1e-4):
    """The kernels work in float32 with fused multiply-adds; the reference is float64, so they agree to float32."""
    scale = float(np.abs(d_reference).max())
    np.testing.assert_allclose(d_kernel, d_reference, rtol=rtol, atol=rtol * 0.1 * scale)
    np.testing.assert_allclose(q_kernel, q_reference, rtol=rtol, atol=rtol * 0.1 * scale)


# ── kernels ─────────────────────────────────────────────────────────────────


@pytest.mark.parametrize('n_substeps', [1, 3])
def test_static_channel_matches_the_reference(n_substeps):
    rng = np.random.default_rng(1)
    downstream = dfs_ordered_tree(rng, N_RIVERS, n_outlets=3)
    _, _, c1, c2, c3 = coefficients(rng, N_RIVERS, DT_RUNOFF // n_substeps)
    q_init = rng.uniform(1, 100, N_RIVERS).astype(np.float32)

    q_expected, d_expected = reference_static(downstream, c1, c2, c3, None, q_init, None, N_STEPS, n_substeps)
    q_kernel, d_kernel = q_init.copy(), river_major(N_STEPS, N_RIVERS)
    _river_kernels.static_channel(
        q_t=q_kernel,
        discharge_array=d_kernel,
        downstream_indices=downstream,
        c1=c1,
        c2=c2,
        c3=c3,
        n_substeps=n_substeps,
        block=_river_kernels.BLOCK,
    )
    assert_matches_reference(q_kernel, d_kernel, q_expected, d_expected)


@pytest.mark.parametrize('by_river', [False, True])
@pytest.mark.parametrize('n_substeps', [1, 3])
def test_static_vlateral_matches_the_reference(n_substeps, by_river):
    rng = np.random.default_rng(2)
    downstream = dfs_ordered_tree(rng, N_RIVERS, n_outlets=3)
    _, _, c1, c2, c3 = coefficients(rng, N_RIVERS, DT_RUNOFF // n_substeps)
    c4_dt = ((c1 + c2) / DT_RUNOFF).astype(np.float32)
    q_init = rng.uniform(0, 50, N_RIVERS).astype(np.float32)
    vlateral = rng.uniform(0, 5e4, (N_STEPS, N_RIVERS)).astype(np.float32)

    q_expected, d_expected = reference_static(downstream, c1, c2, c3, c4_dt, q_init, vlateral, N_STEPS, n_substeps)
    q_kernel, d_kernel = q_init.copy(), river_major(N_STEPS, N_RIVERS)
    _river_kernels.static_vlateral(
        q_t=q_kernel,
        discharge_array=d_kernel,
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
    assert_matches_reference(q_kernel, d_kernel, q_expected, d_expected)


@pytest.mark.parametrize('by_river', [False, True])
@pytest.mark.parametrize('n_substeps', [1, 3])
def test_dynamic_vlateral_matches_the_reference(n_substeps, by_river):
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

    q_expected, d_expected = reference_dynamic(
        downstream, alpha, beta, x, dt_routing, DT_RUNOFF, q_init, vlateral, N_STEPS, n_substeps
    )
    q_kernel, d_kernel = q_init.copy(), river_major(N_STEPS, N_RIVERS)
    _river_kernels.dynamic_vlateral(
        q_t=q_kernel,
        discharge_array=d_kernel,
        downstream_indices=downstream,
        alpha=alpha,
        beta=beta,
        x=x,
        dt_routing=np.float32(dt_routing),
        dt_runoff=np.float32(DT_RUNOFF),
        n_substeps=n_substeps,
        vlateral=np.ascontiguousarray(vlateral.T) if by_river else vlateral,
        by_river=by_river,
        block=_river_kernels.BLOCK,
    )
    # rebuilding the coefficients every substep compounds float32 rounding, so this one is held to a looser bound
    assert_matches_reference(q_kernel, d_kernel, q_expected, d_expected, rtol=1e-3)


def test_discharge_must_be_river_major():
    """A (river, time) array that is not C-order is refused: the kernels never transpose their output."""
    rng = np.random.default_rng(4)
    downstream = dfs_ordered_tree(rng, 8, n_outlets=1)
    _, _, c1, c2, c3 = coefficients(rng, 8, DT_RUNOFF)
    with pytest.raises(ValueError, match='C-order'):
        _river_kernels.static_channel(
            q_t=np.zeros(8, np.float32),
            discharge_array=np.zeros((4, 8), np.float32).T,
            downstream_indices=downstream,
            c1=c1,
            c2=c2,
            c3=c3,
            n_substeps=1,
            block=_river_kernels.BLOCK,
        )


@pytest.mark.parametrize('dtype', [np.float32, np.uint16])
def test_to_time_major_matches_numpy(dtype):
    """The tiled transpose is the one way out of (river, time), so it must equal numpy's copy, dtype for dtype."""
    rng = np.random.default_rng(5)
    by_river = rng.integers(0, 4096, (37, 53)).astype(dtype)  # (river, time), as a writer is handed it

    transposed = _river_kernels.to_time_major(by_river)

    assert transposed.shape == (53, 37)
    assert transposed.flags.c_contiguous and transposed.dtype == dtype
    np.testing.assert_array_equal(transposed, by_river.T)


def test_to_time_major_views_an_already_time_major_buffer():
    """A (river, time) array that is the transpose of a C-order (time, river) buffer is that buffer, not a copy."""
    by_step = np.zeros((4, 6), np.float32)
    viewed = _river_kernels.to_time_major(by_step.T)
    assert viewed.shape == (4, 6) and np.shares_memory(viewed, by_step)


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


def test_registry_covers_the_implemented_combinations():
    for coeff, forcing in (('static', 'channel'), ('static', 'vlateral'), ('dynamic', 'vlateral')):
        assert _kernel_registry.resolve_kernel(coeff, forcing, 'uniform', 'standard') is not None
    for forcing in ('channel', 'vlateral'):
        assert _kernel_registry.resolve_kernel('static', forcing, 'uniform', 'stabilized') is not None
    with pytest.raises(NotImplementedError, match='network_conditioning=expanded'):
        _kernel_registry.resolve_kernel('static', 'vlateral', 'uniform', 'expanded')


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


def route(case: dict, out: Path, forcing: str = 'vlateral', **kwargs) -> np.ndarray:
    options = dict(params_file=case['params_file'], channel_state_init_file=case['channel_state_init_file'])
    if forcing == 'vlateral':
        options['vlateral_files'] = [str(case['vlateral_file'])]
    else:
        options.update(dt_routing=DT_RUNOFF, dt_total=DT_RUNOFF * N_STEPS)
    configs = rr.Configs(
        forcing=forcing, discharge_files=[str(out)], log=False, progress_bar=False, **options, **kwargs
    )
    rr.Router(configs).set_discharge_writer(writers.netcdf_writer).route()
    with xr.open_dataset(out) as ds:
        return ds['Q'].values


@pytest.mark.parametrize('coeff,forcing', [('static', 'vlateral'), ('dynamic', 'vlateral'), ('static', 'channel')])
def test_routing_with_a_thread_pool_matches_single_threaded(tmp_path, coeff, forcing):
    """Regions routed concurrently hand their outlet series to the main stem, which must reproduce one pass."""
    case = write_tree_case(tmp_path, dynamic=coeff == 'dynamic')
    single = route(case, tmp_path / 'q_single.nc', forcing=forcing, coeff=coeff)
    options = dict(params_file=case['params_file'], channel_state_init_file=case['channel_state_init_file'])
    if forcing == 'vlateral':
        options['vlateral_files'] = [str(case['vlateral_file'])]
    else:
        options.update(dt_routing=DT_RUNOFF, dt_total=DT_RUNOFF * N_STEPS)
    out = tmp_path / 'q_threaded.nc'
    configs = rr.Configs(
        forcing=forcing, coeff=coeff, discharge_files=[str(out)], log=False, progress_bar=False, **options
    )
    router = rr.Router(configs)
    with ThreadPoolExecutor(4) as pool:
        router.set_discharge_writer(writers.netcdf_writer).route(thread_pool=pool, threads=4)
    assert len(router.routing_jobs) > 2, 'the network should split into concurrent regions'
    # a confluence adds a region's series after its main stem upstreams, which differs from one pass only in rounding
    with xr.open_dataset(out) as ds:
        np.testing.assert_allclose(ds['Q'].values, single, rtol=1e-5, atol=1e-6 * float(np.abs(single).max()))


def test_writers_receive_a_river_major_array(tmp_path):
    """Every writer is handed the same array: C-order (river, time), the layout the kernels route in."""
    case = write_tree_case(tmp_path)
    received = {}

    def capture(router, dates, discharge_array, discharge_file, runoff_file=''):
        received['array'] = discharge_array.copy(order='K')

    configs = rr.Configs(
        forcing='vlateral',
        params_file=case['params_file'],
        channel_state_init_file=case['channel_state_init_file'],
        vlateral_files=[str(case['vlateral_file'])],
        discharge_files=[str(tmp_path / 'unused.nc')],
        log=False,
        progress_bar=False,
    )
    rr.Router(configs).set_discharge_writer(capture).route()
    array = received['array']
    assert array.shape[0] == N_RIVERS
    assert array.flags.c_contiguous
