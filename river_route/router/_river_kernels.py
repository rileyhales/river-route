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
# the step starts, and _recurrence advances it eight steps per serial multiply-add. The trade is that every step of a
# river's forcing must be available when that river is routed.
#
# Passes. The network is routed in passes over blocks of rivers, the same partition the time order kernels use. A
# region is a contiguous upstream-closed block; its outlet's whole unclamped series is written to its row of `boundary`
# instead of an inflow row, so regions route concurrently. One pass may hold many regions as separate blocks, which
# keeps per-pass overhead to once per thread. The main stem pass runs last and adds each boundary row in `cut_target`
# into the inflow of the river it drains into before routing that river. One pass over the whole network with no
# outlet and no cuts is the single threaded case.
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
# Inputs are NaN free: RunoffGaussianGrid and RunoffVlateral replace NaN with zero when they prepare the forcing.


@numba.njit(cache=True, nogil=True)
def _cuts_by_target(cut_target):
    """The cuts that drain into a river, as positions into cut_target sorted by target, basin outlets (-1) dropped."""
    order = np.argsort(cut_target, kind='mergesort')
    first = 0
    while first < order.shape[0] and cut_target[order[first]] < 0:
        first += 1
    return order[first:]


@numba.njit(cache=True, nogil=True)
def plan_inflow_rows(downstream_indices, block_starts, block_stops, block_outlet, cut_target):
    """
    Peak number of inflow rows live at once when a pass routes its blocks in order. A river's row opens when its first
    upstream is routed, or when a region's buffered outlet series is injected into it, and closes once the river
    itself is routed. A block's outlet drains outside the pass, so it never opens a row.
    """
    first, span = _pass_span(block_starts, block_stops)
    is_open = np.zeros(span, dtype=np.bool_)
    cuts = _cuts_by_target(cut_target)
    next_cut = 0
    live = 0
    peak = 0
    for b in range(block_starts.shape[0]):
        for r in range(block_starts[b], block_stops[b]):
            while next_cut < cuts.shape[0] and cut_target[cuts[next_cut]] == r:
                if not is_open[r - first]:
                    is_open[r - first] = True
                    live += 1
                next_cut += 1
            d = downstream_indices[r]
            if d >= 0 and r != block_outlet[b] and not is_open[d - first]:
                is_open[d - first] = True
                live += 1
            if live > peak:
                peak = live
            if is_open[r - first]:
                live -= 1
    return peak


@numba.njit(cache=True, nogil=True)
def _pass_span(block_starts, block_stops):
    """The first river a pass routes and how many indices its blocks span, which sizes its per-river bookkeeping."""
    first = block_starts.min() if block_starts.shape[0] else 0
    last = block_stops.max() if block_stops.shape[0] else 0
    return first, max(last - first, 0)


@numba.njit(cache=True, nogil=True)
def _inflow_pool(span, n_rows, n_routing_steps):
    """The inflow rows, the row each river holds (-1 for none), and the stack of free rows with its top index."""
    inflow = np.empty((n_rows + 2, n_routing_steps + 1), dtype=np.float32)
    inflow[n_rows, :] = 0.0  # headwaters read this row as their inflow
    slot_of = np.full(span, -1, dtype=np.int64)
    free = np.arange(n_rows)
    top = np.full(1, n_rows, dtype=np.int64)
    return inflow, slot_of, free, top


@numba.njit(cache=True, nogil=True)
def _upstream_row(r, inflow, slot_of):
    """The inflow row river r reads: the one its upstreams filled, or the zero row for a headwater. Every index into
    slot_of is relative to the pass's first river."""
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


@numba.njit(cache=True, nogil=True, fastmath={'contract'})
def _recurrence(q, c3, work, n):
    """
    Replace work[g] with q_g = c3 q_(g-1) + work[g] in place and return the last q.

    Steps are taken eight at a time: within a group the partial sums p_k = sum_j c3^(k-j) work[j] do not depend on q,
    so they overlap with the previous group, and only q_(g+7) = c3^8 q_(g-1) + p_7 waits on the one before it. The
    serial chain is one multiply-add per eight steps instead of one per step.
    """
    c3_2 = c3 * c3
    c3_3 = c3_2 * c3
    c3_4 = c3_2 * c3_2
    c3_5 = c3_4 * c3
    c3_6 = c3_4 * c3_2
    c3_7 = c3_4 * c3_3
    c3_8 = c3_4 * c3_4
    g = 0
    while g + 8 <= n:
        a0 = work[g]
        p1 = c3 * a0 + work[g + 1]
        p2 = c3 * p1 + work[g + 2]
        p3 = c3 * p2 + work[g + 3]
        p4 = c3 * p3 + work[g + 4]
        p5 = c3 * p4 + work[g + 5]
        p6 = c3 * p5 + work[g + 6]
        p7 = c3 * p6 + work[g + 7]
        work[g] = c3 * q + a0
        work[g + 1] = c3_2 * q + p1
        work[g + 2] = c3_3 * q + p2
        work[g + 3] = c3_4 * q + p3
        work[g + 4] = c3_5 * q + p4
        work[g + 5] = c3_6 * q + p5
        work[g + 6] = c3_7 * q + p6
        q = c3_8 * q + p7
        work[g + 7] = q
        g += 8
    while g < n:
        q = c3 * q + work[g]
        work[g] = q
        g += 1
    return q


@numba.njit(cache=True, nogil=True, fastmath={'contract'})
def _route_static_river(q, c1, c2, c3, c4_dt, lateral, up, down, out, work, n_steps, n_substeps):
    """
    Route one river's whole series. Its unclamped series is added into ``down`` for the downstream river, its
    clamped per-step mean is written into ``out``, and its final state is returned. ``lateral`` is None for channel
    routing. ``work`` holds n_steps * n_substeps values of scratch.

    Only the recurrence q = c3 q + forcing is serial, so it runs alone in its own loop; the passes before and after
    it are independent per step and vectorize. Fusing multiply-adds is the only fastmath flag, which keeps NaN
    semantics and makes each serial step one instruction.
    """
    zero = np.float32(0.0)
    n_routing = n_steps * n_substeps

    # forcing of every step, all known before the recurrence starts
    for g in range(n_routing):
        work[g] = c1 * up[g + 1] + c2 * up[g]
    if lateral is not None:
        if n_substeps == 1:
            for t in range(n_steps):
                work[t] += c4_dt * np.float32(lateral[t])
        else:
            for t in range(n_steps):
                external = c4_dt * np.float32(lateral[t])
                for s in range(t * n_substeps, (t + 1) * n_substeps):
                    work[s] += external

    down[0] += q
    q = _recurrence(q, c3, work, n_routing)

    for g in range(n_routing):
        down[g + 1] += work[g]
    # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
    if n_substeps == 1:
        for t in range(n_steps):
            val = work[t]
            out[t] = val if val > zero else zero
    else:
        inv_substeps = np.float32(1.0 / n_substeps)
        for t in range(n_steps):
            interval_sum = zero
            for s in range(t * n_substeps, (t + 1) * n_substeps):
                interval_sum += work[s]
            val = interval_sum * inv_substeps
            out[t] = val if val > zero else zero
    return q
    for g in range(n_routing):
        down[g + 1] += work[g]
    inv_substeps = np.float32(1.0 / n_substeps)
    for t in range(n_steps):
        interval_sum = zero
        for s in range(t * n_substeps, (t + 1) * n_substeps):
            interval_sum += work[s]
        val = interval_sum * inv_substeps
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
def _route_pass(
    q_t,
    discharge_array,
    downstream_indices,
    c1,
    c2,
    c3,
    c4_dt,
    dynamic,
    alpha,
    beta,
    x,
    dt_routing,
    inv_dt_runoff,
    n_substeps,
    lateral_source,
    vlateral,
    runoff_by_cell,
    indptr,
    cell,
    weight,
    scale,
    cumulative,
    force_positive,
    block_starts,
    block_stops,
    block_outlet,
    block_region,
    cut_target,
    boundary,
    block,
    discharge_by_river,
):
    """
    Route every river in the given blocks, river by river; see the module header for the pass contract. A block's
    outlet (-1 for none) hands its unclamped series to ``boundary[block_region]`` instead of an inflow row, and each
    boundary row in ``cut_target`` is injected into its target river before that river is routed.

    ``lateral_source`` is 0 for none, 1 for (time, river) vlateral, 2 for (river, time) vlateral read in place, and 3
    for gridded runoff aggregated per block. ``discharge_by_river`` means discharge_array is C-order (river, time) and
    each river's series is written into its own row in place rather than transposed into (time, river) by block.
    """
    n_rivers, n_steps = discharge_array.shape if discharge_by_river else discharge_array.shape[::-1]
    n_routing = n_steps * n_substeps
    n_rows = plan_inflow_rows(downstream_indices, block_starts, block_stops, block_outlet, cut_target)
    first, span = _pass_span(block_starts, block_stops)
    inflow, slot_of, free, top = _inflow_pool(span, n_rows, n_routing)
    cuts = _cuts_by_target(cut_target)
    next_cut = 0
    scratch = np.zeros((block if lateral_source in (0, 1, 3) else 0, n_steps), dtype=np.float32)
    out = np.empty((0 if discharge_by_river else block, n_steps), dtype=np.float32)
    work = np.empty(n_routing, dtype=np.float32)
    zero = np.float32(0.0)

    for b in range(block_starts.shape[0]):
        outlet = block_outlet[b]
        for r0 in range(block_starts[b], block_stops[b], block):
            r1 = min(r0 + block, block_stops[b])
            if lateral_source == 1:
                _read_block(vlateral, r0, r1, n_steps, scratch)
            elif lateral_source == 3:
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
                        r,
                        scratch[r - r0],
                    )
            for r in range(r0, r1):
                while next_cut < cuts.shape[0] and cut_target[cuts[next_cut]] == r:
                    injected = _downstream_row(r - first, inflow, slot_of, free, top)
                    injected += boundary[cuts[next_cut]]
                    next_cut += 1
                up = _upstream_row(r - first, inflow, slot_of)
                if r == outlet:
                    down = boundary[block_region[b]]
                    down[:] = 0.0
                else:
                    d = downstream_indices[r]
                    down = _downstream_row(d - first if d >= 0 else -1, inflow, slot_of, free, top)
                series = discharge_array[r] if discharge_by_river else out[r - r0]
                lateral = vlateral[r] if lateral_source == 2 else scratch[r - r0]
                if dynamic:
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
                        series,
                        n_steps,
                        n_substeps,
                    )
                else:
                    q_t[r] = _route_static_river(
                        q_t[r], c1[r], c2[r], c3[r], c4_dt[r], lateral, up, down, series, work, n_steps, n_substeps
                    )
                _release_row(r - first, slot_of, free, top)
            if not discharge_by_river:
                _write_block(out, r0, r1, discharge_array)
    return


_NO_FLOATS = np.zeros(0, dtype=np.float32)
_NO_GRID = np.zeros((0, 0), dtype=np.float32)
_NO_INTS = np.zeros(0, dtype=np.int32)
_UNUSED = dict(
    c4_dt=_NO_FLOATS,
    dynamic=False,
    alpha=_NO_FLOATS,
    beta=_NO_FLOATS,
    x=_NO_FLOATS,
    dt_routing=np.float32(0),
    inv_dt_runoff=np.float32(0),
    vlateral=_NO_GRID,
    runoff_by_cell=_NO_GRID,
    indptr=_NO_INTS,
    cell=_NO_INTS,
    weight=_NO_FLOATS,
    scale=_NO_FLOATS,
    cumulative=False,
    force_positive=False,
)


def _route(q_t, schedule: dict, **arguments) -> None:
    """Call _route_pass with the unused arguments of the kernel filled in and the pass defaulted to the whole
    network with no regions."""
    n_rivers = q_t.shape[0]
    pass_arguments = dict(
        block_starts=schedule.get('block_starts', np.array([0], dtype=np.int32)),
        block_stops=schedule.get('block_stops', np.array([n_rivers], dtype=np.int32)),
        block_outlet=schedule.get('block_outlet', np.array([-1], dtype=np.int32)),
        block_region=schedule.get('block_region', np.array([0], dtype=np.int32)),
        cut_target=schedule.get('cut_target', _NO_INTS),
        boundary=schedule.get('boundary', np.zeros((1, 0), dtype=np.float32)),
    )
    discharge, by_river = _discharge_layout(arguments.pop('discharge_array'))
    _route_pass(
        q_t=q_t, discharge_array=discharge, discharge_by_river=by_river, **{**_UNUSED, **arguments}, **pass_arguments
    )


def _discharge_layout(discharge_array):
    """
    The array the kernel writes and whether it writes it by river. A (time, river) discharge array that is the
    transpose of a C-order (river, time) buffer takes each river's series in place; a C-order (time, river) array is
    written by block. Output is never copied, so any other layout is refused.
    """
    if discharge_array.flags.c_contiguous:
        return discharge_array, False
    if discharge_array.T.flags.c_contiguous:
        return discharge_array.T, True
    raise ValueError('discharge_array must be C-order (time, river) or the transpose of a C-order (river, time) array')


def static_channel(*, q_t, discharge_array, downstream_indices, c1, c2, c3, n_substeps, block, **schedule):
    """Route the initial state through the rivers of one pass with no lateral inflow."""
    _route(
        q_t,
        schedule,
        discharge_array=discharge_array,
        downstream_indices=downstream_indices,
        c1=c1,
        c2=c2,
        c3=c3,
        c4_dt=np.zeros_like(c3),
        n_substeps=n_substeps,
        lateral_source=0,
        block=block,
    )


def static_vlateral(
    *, q_t, discharge_array, downstream_indices, c1, c2, c3, c4_dt, n_substeps, vlateral, by_river, block, **schedule
):
    """Route the rivers of one pass with lateral inflow from a vlateral array, (river, time) if by_river."""
    _route(
        q_t,
        schedule,
        discharge_array=discharge_array,
        downstream_indices=downstream_indices,
        c1=c1,
        c2=c2,
        c3=c3,
        c4_dt=c4_dt,
        n_substeps=n_substeps,
        lateral_source=2 if by_river else 1,
        vlateral=vlateral,
        block=block,
    )


def static_grid(
    *,
    q_t,
    discharge_array,
    downstream_indices,
    c1,
    c2,
    c3,
    c4_dt,
    n_substeps,
    runoff_by_cell,
    indptr,
    cell,
    weight,
    scale,
    cumulative,
    force_positive,
    block,
    **schedule,
):
    """
    Aggregate gridded runoff onto the rivers of one pass and route it in the same sweep. No vlateral array is built:
    each block of rivers is aggregated into scratch and routed while it is still in cache. The runoff must be NaN
    free, as RunoffGaussianGrid prepares it.
    """
    _route(
        q_t,
        schedule,
        discharge_array=discharge_array,
        downstream_indices=downstream_indices,
        c1=c1,
        c2=c2,
        c3=c3,
        c4_dt=c4_dt,
        n_substeps=n_substeps,
        lateral_source=3,
        runoff_by_cell=runoff_by_cell,
        indptr=indptr,
        cell=cell,
        weight=weight,
        scale=scale,
        cumulative=cumulative,
        force_positive=force_positive,
        block=block,
    )


def dynamic_vlateral(
    *,
    q_t,
    discharge_array,
    downstream_indices,
    alpha,
    beta,
    x,
    dt_routing,
    dt_runoff,
    n_substeps,
    vlateral,
    by_river,
    block,
    **schedule,
):
    """Nonlinear Muskingum over the rivers of one pass with lateral inflow from a vlateral array."""
    _route(
        q_t,
        schedule,
        discharge_array=discharge_array,
        downstream_indices=downstream_indices,
        c1=_NO_FLOATS,
        c2=_NO_FLOATS,
        c3=_NO_FLOATS,
        dynamic=True,
        alpha=alpha,
        beta=beta,
        x=x,
        dt_routing=np.float32(dt_routing),
        inv_dt_runoff=np.float32(1.0 / dt_runoff),
        n_substeps=n_substeps,
        lateral_source=2 if by_river else 1,
        vlateral=vlateral,
        block=block,
    )
