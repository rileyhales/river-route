import numba
import numpy as np

from ..runoff._numba_kernels import aggregate_river

__all__ = ['BLOCK', 'plan_inflow_rows', 'static_channel', 'static_vlateral', 'static_grid', 'dynamic_vlateral']

BLOCK = 64  # rivers per block of the (time, river) transposes, the block size aggregate_to_rivers was tuned with

# River order kernels: each river's whole time series is routed before the next river, instead of sweeping every river
# once per time step as the time order kernels in _numba_kernels do. Selected with Configs(routing_order='river').
#
# Why it is faster. In DFS order a river's downstream is usually the next index, so the time order sweep is one chain
# of dependent loads and stores: each iteration reads the value the previous one just added to. Muskingum couples a
# river only to its upstreams at the same step and to itself at the previous step, so with every upstream index below
# its downstream index the result is the same when rivers are routed one at a time:
#
#   Q_i(g) = c3_i Q_i(g-1) + c4dt_i vlateral_i(t) + c1_i U_i(g) + c2_i U_i(g-1),   U_i(g) = sum of upstream Q_u(g)
#
# where g counts routing steps and t is the runoff step that g falls in. Everything but c3_i Q_i(g-1) is known before
# the step starts, so only one multiply-add per step is serial. The trade is that every step of a river's forcing must
# be available when that river is routed, and the kernels are single threaded: one call routes the whole network.
#
# Inflow rows. U_i is accumulated in a row of `inflow` that river i's upstreams add their unclamped series into as they
# are routed. Element 0 holds U_i(-1), the sum of the upstreams' initial states, and elements 1.. hold U_i(g). A row
# is only live from when a river's first upstream is routed until the river itself is, so a small pool of rows sized
# by plan_inflow_rows is reused rather than holding an (n_rivers, n_routing_steps) array. Two rows follow the pool: a
# row of zeros that headwaters read as their inflow, and a sink that basin outlets write into.
#
# Layout. Rivers are handled in blocks so the (time, river) C-order arrays are read and written one short contiguous
# run per time step rather than one strided element per river per step. Lateral inflow arrives one of three ways:
#
#   static_vlateral / dynamic_vlateral with by_river=False   (time, river) C-order, copied into a scratch block
#   static_vlateral / dynamic_vlateral with by_river=True    (river, time) C-order, each river's row read in place
#   static_grid                                               gridded runoff aggregated per block, no vlateral array
#
# No fastmath in this module. static_grid calls runoff._numba_kernels.aggregate_river, whose NaN replacement fastmath
# would remove (see the note there), and numba compiles every helper with its caller's flags. The routing loop keeps
# only one multiply-add on the serial chain without fastmath.


@numba.njit(cache=True, nogil=True)
def plan_inflow_rows(downstream_indices):
    """Peak number of inflow rows live at once when routing in index order."""
    n_rivers = downstream_indices.shape[0]
    is_open = np.zeros(n_rivers, dtype=np.bool_)
    live = 0
    peak = 0
    for i in range(n_rivers):
        d = downstream_indices[i]
        if d >= 0 and not is_open[d]:
            is_open[d] = True
            live += 1
            if live > peak:
                peak = live
        if is_open[i]:
            live -= 1
    return peak


@numba.njit(cache=True, nogil=True)
def _inflow_pool(downstream_indices, n_routing_steps):
    """The inflow rows, the row each river holds (-1 for none), and the stack of free rows with its top index."""
    n_rows = plan_inflow_rows(downstream_indices)
    inflow = np.empty((n_rows + 2, n_routing_steps + 1), dtype=np.float32)
    inflow[n_rows, :] = 0.0  # headwaters read this row as their inflow
    slot_of = np.full(downstream_indices.shape[0], -1, dtype=np.int64)
    free = np.arange(n_rows)
    top = np.full(1, n_rows, dtype=np.int64)
    return inflow, slot_of, free, top


@numba.njit(cache=True, nogil=True)
def _upstream_row(r, inflow, slot_of):
    """The inflow row river r reads: the one its upstreams filled, or the zero row for a headwater."""
    s = slot_of[r]
    return inflow[s] if s >= 0 else inflow[inflow.shape[0] - 2]


@numba.njit(cache=True, nogil=True)
def _downstream_row(d, inflow, slot_of, free, top):
    """The inflow row river d accumulates into, taken from the pool and zeroed on first use; the sink if d < 0."""
    if d < 0:
        return inflow[inflow.shape[0] - 1]
    s = slot_of[d]
    if s < 0:
        top[0] -= 1
        s = free[top[0]]
        slot_of[d] = s
        inflow[s, :] = 0.0
    return inflow[s]


@numba.njit(cache=True, nogil=True)
def _release_row(r, slot_of, free, top):
    """River r has been routed, so the inflow row it read goes back to the pool."""
    s = slot_of[r]
    if s >= 0:
        free[top[0]] = s
        top[0] += 1
        slot_of[r] = -1


@numba.njit(cache=True, nogil=True)
def _read_block(vlateral, r0, r1, n_steps, scratch):
    """Copy columns r0:r1 of a (time, river) array into scratch rows, one contiguous run per time step."""
    for t in range(n_steps):
        source = vlateral[t]
        for b in range(r1 - r0):
            scratch[b, t] = source[r0 + b]


@numba.njit(cache=True, nogil=True)
def _write_block(out, r0, r1, discharge_array):
    """Copy out rows into columns r0:r1 of the (time, river) discharge, one contiguous run per time step."""
    for t in range(discharge_array.shape[0]):
        destination = discharge_array[t]
        for b in range(r1 - r0):
            destination[r0 + b] = out[b, t]


@numba.njit(cache=True, nogil=True)
def _route_static_river(q, c1, c2, c3, c4_dt, lateral, up, down, out, n_steps, n_substeps):
    """
    Route one river's whole series. Its unclamped series is added into ``down`` for the downstream river, its
    clamped per-step mean is written into ``out``, and its final state is returned. ``lateral`` is None for channel
    routing.
    """
    zero = np.float32(0.0)
    u_prev = up[0]
    down[0] += q
    inv_substeps = np.float32(1.0 / n_substeps)
    g = 1
    for t in range(n_steps):
        external = zero if lateral is None else c4_dt * np.float32(lateral[t])
        if n_substeps == 1:
            u = up[g]
            q = c3 * q + (external + c1 * u + c2 * u_prev)
            u_prev = u
            down[g] += q
            g += 1
            # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
            out[t] = q if q > zero else zero
            continue
        interval_sum = zero
        for _ in range(n_substeps):
            u = up[g]
            q = c3 * q + (external + c1 * u + c2 * u_prev)
            u_prev = u
            down[g] += q
            interval_sum += q
            g += 1
        val = interval_sum * inv_substeps
        # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
        out[t] = val if val > zero else zero
    return q


@numba.njit(cache=True, nogil=True)
def _route_dynamic_river(q, alpha, beta, x, dt_routing, inv_dt_runoff, lateral, up, down, out, n_steps, n_substeps):
    """
    Nonlinear Muskingum for one river's whole series. The coefficients are rebuilt from the river's own discharge
    every substep, as in _numba_kernels.dynamic_vlateral, and the upstream series is weighted by them.
    """
    zero = np.float32(0.0)
    qmin = np.float32(1e-6)
    two_x = np.float32(2.0) * x
    two_one_minus_x = np.float32(2.0) * (np.float32(1.0) - x)
    inv_substeps = np.float32(1.0 / n_substeps)
    u_prev = up[0]
    down[0] += q
    g = 1
    for t in range(n_steps):
        external = np.float32(lateral[t]) * inv_dt_runoff
        interval_sum = zero
        for _ in range(n_substeps):
            k = alpha * (q if q > qmin else qmin) ** beta
            dt_div_k = dt_routing / k
            denominator = dt_div_k + two_one_minus_x
            c1 = (dt_div_k - two_x) / denominator
            c2 = (dt_div_k + two_x) / denominator
            c3 = (two_one_minus_x - dt_div_k) / denominator
            u = up[g]
            q = c3 * q + ((c1 + c2) * external + c1 * u + c2 * u_prev)
            u_prev = u
            down[g] += q
            interval_sum += q
            g += 1
        val = interval_sum * inv_substeps
        # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
        out[t] = val if val > zero else zero
    return q


@numba.njit(cache=True, nogil=True)
def static_channel(
    *,
    q_t,  # Array (n_rivers,) of discharge state, updated in-place to the final state
    discharge_array,  # Array (n_steps, n_rivers) to write discharge time series into
    downstream_indices,  # Array (n_rivers,) of downstream river indices, -1 for none; upstream always below downstream
    c1,  # Array (n_rivers,) of c1 for river at index i
    c2,  # Array (n_rivers,) of c2 for river at index i
    c3,  # Array (n_rivers,) of c3 for river at index i
    n_substeps,  # integer number of routing substeps per output step
    block,  # integer rivers per block of the discharge transpose
):
    """Route the initial state through the network river by river with no lateral inflow; see the module header."""
    n_steps, n_rivers = discharge_array.shape
    inflow, slot_of, free, top = _inflow_pool(downstream_indices, n_steps * n_substeps)
    out = np.empty((block, n_steps), dtype=np.float32)
    for r0 in range(0, n_rivers, block):
        r1 = min(r0 + block, n_rivers)
        for r in range(r0, r1):
            up = _upstream_row(r, inflow, slot_of)
            down = _downstream_row(downstream_indices[r], inflow, slot_of, free, top)
            q_t[r] = _route_static_river(
                q_t[r], c1[r], c2[r], c3[r], np.float32(0.0), None, up, down, out[r - r0], n_steps, n_substeps
            )
            _release_row(r, slot_of, free, top)
        _write_block(out, r0, r1, discharge_array)
    return


@numba.njit(cache=True, nogil=True)
def static_vlateral(
    *,
    q_t,  # Array (n_rivers,) of discharge state, updated in-place to the final state
    discharge_array,  # Array (n_steps, n_rivers) to write discharge time series into
    downstream_indices,  # Array (n_rivers,) of downstream river indices, -1 for none; upstream always below downstream
    c1,  # Array (n_rivers,) of c1 for river at index i
    c2,  # Array (n_rivers,) of c2 for river at index i
    c3,  # Array (n_rivers,) of c3 for river at index i
    c4_dt,  # Array (n_rivers,) of c4 / dt_runoff for river at index i, scales lateral volumes to a rate
    n_substeps,  # integer number of routing substeps per runoff value
    vlateral,  # Array of lateral volumes, C-order (river, >= n_steps) if by_river else (>= n_steps, river)
    by_river,  # bool, vlateral is (river, time) and each river's row is read in place
    block,  # integer rivers per block of the transposes
):
    """Route the whole network river by river with lateral inflow from a vlateral array; see the module header."""
    n_steps, n_rivers = discharge_array.shape
    inflow, slot_of, free, top = _inflow_pool(downstream_indices, n_steps * n_substeps)
    scratch = np.empty((0 if by_river else block, n_steps), dtype=np.float32)
    out = np.empty((block, n_steps), dtype=np.float32)
    for r0 in range(0, n_rivers, block):
        r1 = min(r0 + block, n_rivers)
        if not by_river:
            _read_block(vlateral, r0, r1, n_steps, scratch)
        for r in range(r0, r1):
            lateral = vlateral[r] if by_river else scratch[r - r0]
            up = _upstream_row(r, inflow, slot_of)
            down = _downstream_row(downstream_indices[r], inflow, slot_of, free, top)
            q_t[r] = _route_static_river(
                q_t[r], c1[r], c2[r], c3[r], c4_dt[r], lateral, up, down, out[r - r0], n_steps, n_substeps
            )
            _release_row(r, slot_of, free, top)
        _write_block(out, r0, r1, discharge_array)
    return


@numba.njit(cache=True, nogil=True)
def static_grid(
    *,
    q_t,  # Array (n_rivers,) of discharge state, updated in-place to the final state
    discharge_array,  # Array (n_steps, n_rivers) to write discharge time series into
    downstream_indices,  # Array (n_rivers,) of downstream river indices, -1 for none; upstream always below downstream
    c1,  # Array (n_rivers,) of c1 for river at index i
    c2,  # Array (n_rivers,) of c2 for river at index i
    c3,  # Array (n_rivers,) of c3 for river at index i
    c4_dt,  # Array (n_rivers,) of c4 / dt_runoff for river at index i, scales lateral volumes to a rate
    n_substeps,  # integer number of routing substeps per runoff value
    runoff_by_cell,  # Array (n_cells, >= n_steps) of runoff depths, each cell's time series contiguous
    indptr,  # Array (n_rivers + 1,) of sparse row pointers, the weights of river r are indptr[r]:indptr[r + 1]
    cell,  # Array (n_weights,) of the row of runoff_by_cell each weight applies to
    weight,  # Array (n_weights,) of area proportions, already multiplied by the depth unit conversion factor
    scale,  # Array (n_rivers,) of per-river multipliers (catchment area for volumes); EMPTY to skip scaling
    zero,  # scalar zero in the dtype of scratch
    cumulative,  # bool, de-accumulate each river's series from cumulative to incremental
    force_positive,  # bool, clip negative runoff to zero
    replace_nan,  # bool, replace NaN runoff with zero
    scratch,  # Array (rivers_per_block, n_steps) of working space for one block's lateral inflow series
):
    """
    Aggregate gridded runoff onto rivers and route it river by river in one pass. No vlateral array is built: each
    block of rivers is aggregated into scratch and routed while it is still in cache.
    """
    n_steps, n_rivers = discharge_array.shape
    block = scratch.shape[0]
    inflow, slot_of, free, top = _inflow_pool(downstream_indices, n_steps * n_substeps)
    out = np.empty((block, n_steps), dtype=np.float32)
    for r0 in range(0, n_rivers, block):
        r1 = min(r0 + block, n_rivers)
        for r in range(r0, r1):
            aggregate_river(
                runoff_by_cell,
                indptr,
                cell,
                weight,
                scale,
                zero,
                cumulative,
                force_positive,
                replace_nan,
                r,
                scratch[r - r0],
            )
        for r in range(r0, r1):
            up = _upstream_row(r, inflow, slot_of)
            down = _downstream_row(downstream_indices[r], inflow, slot_of, free, top)
            q_t[r] = _route_static_river(
                q_t[r], c1[r], c2[r], c3[r], c4_dt[r], scratch[r - r0], up, down, out[r - r0], n_steps, n_substeps
            )
            _release_row(r, slot_of, free, top)
        _write_block(out, r0, r1, discharge_array)
    return


@numba.njit(cache=True, nogil=True)
def dynamic_vlateral(
    *,
    q_t,  # Array (n_rivers,) of discharge state, updated in-place to the final state
    discharge_array,  # Array (n_steps, n_rivers) to write discharge time series into
    downstream_indices,  # Array (n_rivers,) of downstream river indices, -1 for none; upstream always below downstream
    alpha,  # Array (n_rivers,) of alpha for river at index i, K_i = alpha_i * q^beta_i
    beta,  # Array (n_rivers,) of beta for river at index i
    x,  # Array (n_rivers,) of x for river at index i
    dt_routing,  # float32 routing timestep in seconds
    dt_runoff,  # float32 runoff timestep in seconds, scales lateral volumes to a rate
    n_substeps,  # integer number of routing substeps per runoff value
    vlateral,  # Array of lateral volumes, C-order (river, >= n_steps) if by_river else (>= n_steps, river)
    by_river,  # bool, vlateral is (river, time) and each river's row is read in place
    block,  # integer rivers per block of the transposes
):
    """Nonlinear Muskingum over the whole network river by river; see the module header."""
    n_steps, n_rivers = discharge_array.shape
    inflow, slot_of, free, top = _inflow_pool(downstream_indices, n_steps * n_substeps)
    inv_dt_runoff = np.float32(1.0 / dt_runoff)
    scratch = np.empty((0 if by_river else block, n_steps), dtype=np.float32)
    out = np.empty((block, n_steps), dtype=np.float32)
    for r0 in range(0, n_rivers, block):
        r1 = min(r0 + block, n_rivers)
        if not by_river:
            _read_block(vlateral, r0, r1, n_steps, scratch)
        for r in range(r0, r1):
            lateral = vlateral[r] if by_river else scratch[r - r0]
            up = _upstream_row(r, inflow, slot_of)
            down = _downstream_row(downstream_indices[r], inflow, slot_of, free, top)
            q_t[r] = _route_dynamic_river(
                q_t[r],
                alpha[r],
                beta[r],
                x[r],
                dt_routing,
                inv_dt_runoff,
                lateral,
                up,
                down,
                out[r - r0],
                n_steps,
                n_substeps,
            )
            _release_row(r, slot_of, free, top)
        _write_block(out, r0, r1, discharge_array)
    return
