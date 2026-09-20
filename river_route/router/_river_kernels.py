import numba
import numpy as np
from llvmlite import ir
from numba import types
from numba.extending import intrinsic, overload

from ..runoff._numba_kernels import aggregate_river

__all__ = [
    'BLOCK',
    'TRANSPOSE_TILE',
    'decode_discharge',
    'to_time_major',
    'plan_inflow_rows',
    'static_channel',
    'static_vlateral',
    'static_grid',
    'dynamic_vlateral',
]

BLOCK = 64  # rivers per block of the (time, river) transposes, the block size aggregate_to_rivers was tuned with

# Discharge storage dtype. The routing math is always float32; only the stored copy is narrowed, so rounding never
# feeds back into the recurrence. numba has no float16 type on the CPU, so a narrowed buffer is carried as uint16 bit
# patterns the kernels fill through an LLVM fptrunc, and viewed as float16 by whoever reads it. float16 keeps 11
# significand bits but only covers 6.1e-5 to 65,504, so it suits networks whose flows stay inside that range.


@intrinsic
def _f32_to_f16_bits(typingctx, x):
    """float32 -> the 16 bits of its float16, as one native conversion instruction."""

    def codegen(context, builder, signature, args):
        return builder.bitcast(builder.fptrunc(args[0], ir.HalfType()), ir.IntType(16))

    return types.uint16(types.float32), codegen


def store_discharge(out, t, val):
    """Store one float32 discharge value into out[t], narrowed to float16 bits when out is uint16."""
    out[t] = val


@overload(store_discharge, inline='always')
def _store_discharge(out, t, val):
    """
    Compile the store for the dtype of the output array: a plain store for float, a narrowing one for uint16.

    This dispatch is what keeps a narrowed buffer correct. numba would otherwise compile ``out[t] = val`` on a uint16
    array as a NUMERIC cast, silently storing int(discharge).
    """
    if isinstance(out.dtype, types.Float):

        def impl(out, t, val):
            out[t] = val

        return impl
    if isinstance(out.dtype, types.Integer) and out.dtype.bitwidth == 16:

        def impl(out, t, val):
            out[t] = _f32_to_f16_bits(np.float32(val))

        return impl
    raise TypeError('discharge_array must hold float32, or uint16 for a narrowed dtype')


def decode_discharge(stored):
    """The stored discharge as the float dtype it holds: float32 passes through, uint16 is float16 bit patterns."""
    return stored if stored.dtype == np.float32 else stored.view(np.float16)


TRANSPOSE_TILE = 32  # rows and columns per tile of the discharge transpose, sized so a tile stays in cache


@numba.njit(cache=True, nogil=True)
def _transpose_tiles(by_river, out, tile):
    """Copy a (river, time) array into a C-order (time, river) one, a square tile at a time."""
    n_rivers, n_steps = by_river.shape
    for r0 in range(0, n_rivers, tile):
        r1 = min(r0 + tile, n_rivers)
        for t0 in range(0, n_steps, tile):
            t1 = min(t0 + tile, n_steps)
            for t in range(t0, t1):
                row = out[t]
                for r in range(r0, r1):
                    row[r] = by_river[r, t]


def to_time_major(discharge_array):
    """
    A C-order (time, river) copy of a routed (river, time) discharge array, for a writer or consumer whose format
    needs each time step's rivers contiguous.

    The kernels always write (river, time), because that is the layout they solve in. Transposing is a scatter
    however it is done, so this does it a tile at a time to keep both sides in cache rather than striding the whole
    array per column. An array that is already the transpose of a C-order (time, river) buffer is returned as that
    buffer's view, with nothing copied.
    """
    if discharge_array.T.flags.c_contiguous:
        return discharge_array.T
    by_river = np.ascontiguousarray(discharge_array)
    n_rivers, n_steps = by_river.shape
    out = np.empty((n_steps, n_rivers), dtype=by_river.dtype)
    _transpose_tiles(by_river, out, TRANSPOSE_TILE)
    return out


# Routing kernels: each river's whole time series is routed before the next river. The recurrence they solve, how
# passes and regions are scheduled, how inflow rows are pooled, and which layouts the lateral inflow may arrive in
# are described in docs/references/kernels.md, under "How routing works".
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
def _route_static_river(
    q_t, p0, n_pieces, m, c1, c2, c3, c4_dt, lateral, up, down, out, work, chain, n_steps, n_substeps
):
    """
    Route one river's whole series through its ``n_pieces`` equal sub-reaches in series, whose states are
    ``q_t[p0:p0 + n_pieces]``: one piece on a standard network, and as many as stability needs on a stabilized one.
    The last piece's unclamped series is added into ``down`` for the downstream river and its clamped per-step mean is
    written into ``out``; the states are updated in place. ``lateral`` is None for channel routing, and every piece
    takes an equal share of it through ``c4_dt``.

    ``m`` sub-cycles the river: each routing step is taken as m equal steps of its own, which is how a river too short
    for the routing step is kept stable. The upstream series arrives once per routing step and is interpolated
    linearly between those levels, and the downstream river is handed the level at the end of each routing step,
    since Muskingum inflow is the discharge at the time levels. With m = 1 the river is routed at the routing step.

    ``work`` holds n_steps * n_substeps * m values of scratch and ``chain`` one more, which carries a piece's series
    to the next piece in the same layout as an inflow row.

    Only the recurrence q = c3 q + forcing is serial, so it runs alone in its own loop; the passes before and after
    it are independent per step and vectorize. Fusing multiply-adds is the only fastmath flag, which keeps NaN
    semantics and makes each serial step one instruction.
    """
    zero = np.float32(0.0)
    n_routing = n_steps * n_substeps
    n_per_step = n_substeps * m  # the river's own steps per runoff step
    n_fine = n_steps * n_per_step
    source = up
    for j in range(n_pieces):
        # forcing of every step, all known before the recurrence starts
        if j == 0 and m > 1:
            inv_m = np.float32(1.0 / m)
            for g in range(n_routing):
                start = up[g]
                rise = (up[g + 1] - start) * inv_m
                before = start
                for s in range(m):
                    after = start + rise * np.float32(s + 1)
                    work[g * m + s] = c1 * after + c2 * before
                    before = after
        else:
            for h in range(n_fine):
                work[h] = c1 * source[h + 1] + c2 * source[h]
        if lateral is not None:
            if n_per_step == 1:
                for t in range(n_steps):
                    work[t] += c4_dt * np.float32(lateral[t])
            else:
                for t in range(n_steps):
                    external = c4_dt * np.float32(lateral[t])
                    for h in range(t * n_per_step, (t + 1) * n_per_step):
                        work[h] += external

        # every step of the source has been read, so the chain row can take this piece's series in its place
        q = q_t[p0 + j]
        if j == n_pieces - 1:
            down[0] += q
            q_t[p0 + j] = _recurrence(q, c3, work, n_fine)
            if m == 1:
                for g in range(n_routing):
                    down[g + 1] += work[g]
            else:
                for g in range(n_routing):
                    down[g + 1] += work[(g + 1) * m - 1]
        else:
            chain[0] = q
            q_t[p0 + j] = _recurrence(q, c3, work, n_fine)
            for h in range(n_fine):
                chain[h + 1] = work[h]
            source = chain

    # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
    if n_per_step == 1:
        for t in range(n_steps):
            val = work[t]
            store_discharge(out, t, val if val > zero else zero)
    else:
        inv_per_step = np.float32(1.0 / n_per_step)
        for t in range(n_steps):
            interval_sum = zero
            for h in range(t * n_per_step, (t + 1) * n_per_step):
                interval_sum += work[h]
            val = interval_sum * inv_per_step
            store_discharge(out, t, val if val > zero else zero)
    return


@numba.njit(cache=True, nogil=True)
def _route_dynamic_river(q, alpha, beta, x, dt_routing, inv_dt_runoff, lateral, up, down, out, n_steps, n_substeps):
    """
    Nonlinear Muskingum for one river's whole series. The coefficients are rebuilt from the river's own discharge
    every substep, and the upstream series is weighted by them.
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
        store_discharge(out, t, val if val > zero else zero)
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
    reach_indptr,
    substeps,
):
    """
    Route every river in the given blocks, river by river; see the module header for the pass contract. A block's
    outlet (-1 for none) hands its unclamped series to ``boundary[block_region]`` instead of an inflow row, and each
    boundary row in ``cut_target`` is injected into its target river before that river is routed.

    ``lateral_source`` is 0 for none, 1 for (time, river) vlateral, 2 for (river, time) vlateral read in place, and 3
    for gridded runoff aggregated per block. ``discharge_array`` is C-order (river, time): each river's series is
    written into its own row in place.

    ``reach_indptr`` is empty on a standard network, where ``q_t`` holds one state per river. On a stabilized network
    river r is the static sub-reaches ``reach_indptr[r]:reach_indptr[r + 1]`` in series, and ``q_t`` holds one state
    per sub-reach. ``substeps`` is empty when no river is sub-cycled, and otherwise gives the steps each river takes
    per routing step.
    """
    n_rivers, n_steps = discharge_array.shape
    n_routing = n_steps * n_substeps
    n_rows = plan_inflow_rows(downstream_indices, block_starts, block_stops, block_outlet, cut_target)
    first, span = _pass_span(block_starts, block_stops)
    inflow, slot_of, free, top = _inflow_pool(span, n_rows, n_routing)
    cuts = _cuts_by_target(cut_target)
    next_cut = 0
    scratch = np.zeros((block if lateral_source in (0, 1, 3) else 0, n_steps), dtype=np.float32)
    expanded = reach_indptr.shape[0] > 0
    cycled = substeps.shape[0] > 0
    most = substeps.max() if cycled else 1
    work = np.empty(n_routing * most, dtype=np.float32)
    chain = np.empty(n_routing * most + 1 if expanded else 0, dtype=np.float32)
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
                series = discharge_array[r]
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
                    p0 = reach_indptr[r] if expanded else r
                    n_pieces = reach_indptr[r + 1] - p0 if expanded else 1
                    _route_static_river(
                        q_t,
                        p0,
                        n_pieces,
                        substeps[r] if cycled else 1,
                        c1[r],
                        c2[r],
                        c3[r],
                        c4_dt[r],
                        lateral,
                        up,
                        down,
                        series,
                        work,
                        chain,
                        n_steps,
                        n_substeps,
                    )
                _release_row(r - first, slot_of, free, top)
    return


_NO_FLOATS = np.zeros(0, dtype=np.float32)
_NO_GRID = np.zeros((0, 0), dtype=np.float32)
_NO_INTS = np.zeros(0, dtype=np.int32)
_NO_REACHES = np.zeros(0, dtype=np.int64)  # a standard network: one reach per river
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
    reach_indptr=_NO_REACHES,
    substeps=_NO_REACHES,
)


def _route(q_t, schedule: dict, **arguments) -> None:
    """Call _route_pass with the unused arguments of the kernel filled in and the pass defaulted to the whole
    network with no regions."""
    n_rivers = arguments['downstream_indices'].shape[0]
    pass_arguments = dict(
        block_starts=schedule.get('block_starts', np.array([0], dtype=np.int32)),
        block_stops=schedule.get('block_stops', np.array([n_rivers], dtype=np.int32)),
        block_outlet=schedule.get('block_outlet', np.array([-1], dtype=np.int32)),
        block_region=schedule.get('block_region', np.array([0], dtype=np.int32)),
        cut_target=schedule.get('cut_target', _NO_INTS),
        boundary=schedule.get('boundary', np.zeros((1, 0), dtype=np.float32)),
    )
    discharge = arguments.pop('discharge_array')
    # each river's series is written into its own contiguous row, which is the layout the kernels solve in. Output
    # is never copied, so any other layout is refused rather than transposed.
    if not discharge.flags.c_contiguous:
        raise ValueError('discharge_array must be a C-order (river, time) array')
    _route_pass(q_t=q_t, discharge_array=discharge, **{**_UNUSED, **arguments}, **pass_arguments)


def static_channel(
    *,
    q_t,
    discharge_array,
    downstream_indices,
    c1,
    c2,
    c3,
    n_substeps,
    block,
    reach_indptr=_NO_REACHES,
    substeps=_NO_REACHES,
    **schedule,
):
    """Route the initial state through the rivers of one pass with no lateral inflow. A nonempty reach_indptr routes
    each river as its sub-reaches, with one state per sub-reach in q_t, and a nonempty substeps
    sub-cycles each river in that many steps per routing step."""
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
        reach_indptr=reach_indptr,
        substeps=substeps,
    )


def static_vlateral(
    *,
    q_t,
    discharge_array,
    downstream_indices,
    c1,
    c2,
    c3,
    c4_dt,
    n_substeps,
    vlateral,
    by_river,
    block,
    reach_indptr=_NO_REACHES,
    substeps=_NO_REACHES,
    **schedule,
):
    """Route the rivers of one pass with lateral inflow from a vlateral array, (river, time) if by_river. A nonempty
    reach_indptr routes each river as its sub-reaches, with one state per sub-reach in q_t, and a nonempty substeps
    sub-cycles each river in that many steps per routing step."""
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
        reach_indptr=reach_indptr,
        substeps=substeps,
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
    reach_indptr=_NO_REACHES,
    substeps=_NO_REACHES,
    **schedule,
):
    """
    Aggregate gridded runoff onto the rivers of one pass and route it in the same sweep. No vlateral array is built:
    each block of rivers is aggregated into scratch and routed while it is still in cache. The runoff must be NaN
    free, as RunoffGaussianGrid prepares it. A nonempty reach_indptr routes each river as its sub-reaches, with one
    state per sub-reach in q_t.
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
        reach_indptr=reach_indptr,
        substeps=substeps,
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
