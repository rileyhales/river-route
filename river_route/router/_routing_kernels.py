"""
Routing kernels: each river's whole time series is routed before the next river. The recurrence they solve, how
passes and regions are scheduled, how inflow rows are pooled, and which layouts the lateral inflow may arrive in are
described in docs/references/kernels.md, under "How routing works". Inputs are NaN free: GaussianGridRunoff and
CatchmentRunoff replace NaN with zero when they prepare the forcing.

The arguments of a routing pass are the NamedTuples in _routing_kernel_arguments.py. An input a pass does not use is
None, and numba compiles a version of the pass for each combination with the branches that read the missing inputs
removed.

The runoff a pass reads is one argument whose type chooses how it is read, through the overloads of
_count_runoff_scratch_rows, _count_rivers_per_runoff_group, _prepare_runoff_for_river_group, and
_get_river_runoff_series: None for channel routing, a CatchmentByRiver array of catchment runoff read in place, or the
gaussian grid runoff as the CellRunoff its reader builds. Adding a kind of runoff is a type and its overloads. The
overloads do not use inline='always': inlining into route_scheduled_rivers at the numba IR level (numba 0.67) made the
whole-array in-place add of a region's outlet series compile as a new array rather than into the inflow row, which
dropped that series from the main stem of threaded routing.
"""

import numba
import numpy as np
from numba import types
from numba.extending import overload

from ..runoff._numba_kernels import aggregate_river
from ..runoff.GaussianGridRunoff import CellRunoff
from ._inflow_rows import (
    allocate_inflow_row_pool,
    count_peak_live_inflow_rows,
    first_river_and_span_of_pass,
    get_or_open_downstream_inflow_row,
    get_upstream_inflow_row,
    release_inflow_row_to_pool,
    sort_region_cuts_by_target_river,
)
from ._routing_kernel_arguments import CatchmentByRiver

__all__ = ['route_scheduled_rivers']

_GRID_RIVERS_AGGREGATED_TOGETHER = 64


@numba.njit(cache=True, nogil=True)
def _aggregate_grid_runoff_for_rivers(grid, first, stop, scratch):
    """Aggregate the gridded runoff of rivers first:stop into scratch rows 0:stop - first."""
    for r in range(first, stop):
        aggregate_river(
            grid.runoff_by_cell,
            grid.indptr,
            grid.cell,
            grid.weight,
            grid.scale,
            np.float32(0.0),
            grid.cumulative,
            grid.force_positive,
            r,
            scratch[r - first],
        )


@numba.njit(cache=True, nogil=True, fastmath={'contract'})
def _solve_muskingum_recurrence_in_place(q, c3, work, n):
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


def _is_runoff_type(source, kind):
    """Whether the numba type of a runoff argument is the NamedTuple class ``kind``."""
    return isinstance(source, types.BaseNamedTuple) and source.instance_class is kind


def _count_runoff_scratch_rows(runoff):
    """Rows of scratch the runoff is prepared in: one per river aggregated together for gridded runoff, else none."""


@overload(_count_runoff_scratch_rows)
def _count_runoff_scratch_rows_overload(runoff):
    if _is_runoff_type(runoff, CellRunoff):
        return lambda runoff: _GRID_RIVERS_AGGREGATED_TOGETHER
    if isinstance(runoff, types.NoneType) or _is_runoff_type(runoff, CatchmentByRiver):
        return lambda runoff: 0
    raise TypeError(f'no routing kernel reads runoff of type {runoff}')


def _count_rivers_per_runoff_group(runoff, n_rivers):
    """Rivers a pass prepares the runoff of together, out of the n_rivers of a schedule block: the rivers aggregated
    together for gridded runoff, else the whole block, which then needs no preparing."""


@overload(_count_rivers_per_runoff_group)
def _count_rivers_per_runoff_group_overload(runoff, n_rivers):
    if _is_runoff_type(runoff, CellRunoff):
        return lambda runoff, n_rivers: _GRID_RIVERS_AGGREGATED_TOGETHER
    return lambda runoff, n_rivers: max(n_rivers, 1)


def _prepare_runoff_for_river_group(runoff, group_start, group_stop, scratch):
    """Prepare the runoff of rivers group_start:group_stop in scratch before they are routed, if it needs it."""


@overload(_prepare_runoff_for_river_group)
def _prepare_runoff_for_river_group_overload(runoff, group_start, group_stop, scratch):
    if _is_runoff_type(runoff, CellRunoff):
        return lambda runoff, group_start, group_stop, scratch: _aggregate_grid_runoff_for_rivers(
            runoff, group_start, group_stop, scratch
        )
    return lambda runoff, group_start, group_stop, scratch: None


def _get_river_runoff_series(runoff, r, group_start, scratch):
    """River r's runoff series: its row of a (river, time) array, its scratch row of the gridded runoff its group
    aggregated, or None for channel routing."""


@overload(_get_river_runoff_series)
def _get_river_runoff_series_overload(runoff, r, group_start, scratch):
    if isinstance(runoff, types.NoneType):
        return lambda runoff, r, group_start, scratch: None
    if _is_runoff_type(runoff, CatchmentByRiver):
        return lambda runoff, r, group_start, scratch: runoff.runoff[r]
    return lambda runoff, r, group_start, scratch: scratch[r - group_start]


@numba.njit(cache=True, nogil=True)
def _interpolate_upstream_inflow_for_substeps(up, m, c1, c2, work, n_routing):
    """
    The upstream forcing of a river sub-cycled in m steps per routing step, with the upstream series interpolated
    linearly between the routing step levels: work[g * m + s] = c1 * after + c2 * before for each of its own steps.

    Compiled without fastmath contract, unlike the rest of the river routine. With contract, whether these multiply-
    adds were fused differed between freshly compiled code and the same code loaded from the numba cache, so a
    sub-cycled river's discharge changed in the last bit depending on whether the kernels had just been compiled.
    """
    inv_m = np.float32(1.0 / m)
    for g in range(n_routing):
        start = up[g]
        rise = (up[g + 1] - start) * inv_m
        before = start
        for s in range(m):
            after = start + rise * np.float32(s + 1)
            work[g * m + s] = c1 * after + c2 * before
            before = after


@numba.njit(cache=True, nogil=True, fastmath={'contract'})
def _route_river_with_static_coefficients(
    q_t, r, coefficients, layout, lateral, up, down, out, work, chain, n_steps, n_substeps
):
    """
    Route river r's whole series through its ``n_pieces`` equal sub-reaches in series, whose states are
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
    expanded = layout.reach_indptr.shape[0] > 0
    p0 = layout.reach_indptr[r] if expanded else r
    n_pieces = layout.reach_indptr[r + 1] - p0 if expanded else 1
    m = layout.substeps[r] if layout.substeps.shape[0] > 0 else 1
    c1 = coefficients.c1[r]
    c2 = coefficients.c2[r]
    c3 = coefficients.c3[r]
    c4_dt = coefficients.c4_dt[r]
    zero = np.float32(0.0)
    n_routing = n_steps * n_substeps
    n_per_step = n_substeps * m  # the river's own steps per runoff step
    n_fine = n_steps * n_per_step
    source = up
    for j in range(n_pieces):
        # forcing of every step, all known before the recurrence starts
        if j == 0 and m > 1:
            _interpolate_upstream_inflow_for_substeps(up, m, c1, c2, work, n_routing)
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
            q_t[p0 + j] = _solve_muskingum_recurrence_in_place(q, c3, work, n_fine)
            if m == 1:
                for g in range(n_routing):
                    down[g + 1] += work[g]
            else:
                for g in range(n_routing):
                    down[g + 1] += work[(g + 1) * m - 1]
        else:
            chain[0] = q
            q_t[p0 + j] = _solve_muskingum_recurrence_in_place(q, c3, work, n_fine)
            for h in range(n_fine):
                chain[h + 1] = work[h]
            source = chain

    # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
    if n_per_step == 1:
        for t in range(n_steps):
            val = work[t]
            out[t] = val if val > zero else zero
    else:
        inv_per_step = np.float32(1.0 / n_per_step)
        for t in range(n_steps):
            interval_sum = zero
            for h in range(t * n_per_step, (t + 1) * n_per_step):
                interval_sum += work[h]
            val = interval_sum * inv_per_step
            out[t] = val if val > zero else zero
    return


@numba.njit(cache=True, nogil=True)
def _route_river_with_dynamic_coefficients(q, r, dynamic, lateral, up, down, out, n_steps, n_substeps):
    """
    Nonlinear Muskingum for river r's whole series. The coefficients are rebuilt from the river's own discharge
    every substep, and the upstream series is weighted by them.
    """
    alpha = dynamic.alpha[r]
    beta = dynamic.beta[r]
    x = dynamic.x[r]
    dt_routing = dynamic.dt_routing
    inv_dt_runoff = dynamic.inv_dt_runoff
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
def route_scheduled_rivers(
    q_t, discharge_array, downstream_indices, n_substeps, coefficients, dynamic, layout, runoff, schedule
):
    """
    Route every river in the blocks of ``schedule``, river by river.

    One of ``coefficients`` and ``dynamic`` is given and the other is None. ``layout`` describes how each river is
    routed, and ``runoff`` is the runoff routed into the rivers: None, a CatchmentByRiver, or a CellRunoff.
    ``discharge_array`` is C-order (river, time): each river's series is written into its own row in place.
    """
    n_rivers, n_steps = discharge_array.shape
    n_routing = n_steps * n_substeps
    block_starts, block_stops, block_outlet, block_region, cut_target, boundary, out_row = schedule
    n_rows = count_peak_live_inflow_rows(downstream_indices, block_starts, block_stops, block_outlet, cut_target)
    first, span = first_river_and_span_of_pass(block_starts, block_stops)
    inflow, slot_of, free, top = allocate_inflow_row_pool(span, n_rows, n_routing)
    cuts = sort_region_cuts_by_target_river(cut_target)
    next_cut = 0
    scratch = np.zeros((_count_runoff_scratch_rows(runoff), n_steps), dtype=np.float32)
    expanded = layout.reach_indptr.shape[0] > 0
    most = layout.substeps.max() if layout.substeps.shape[0] > 0 else 1
    work = np.empty(n_routing * most, dtype=np.float32)
    chain = np.empty(n_routing * most + 1 if expanded else 0, dtype=np.float32)
    renumbered = out_row.shape[0] > 0
    discarded = np.empty(n_steps if renumbered else 0, dtype=discharge_array.dtype)

    for b in range(block_starts.shape[0]):
        outlet = block_outlet[b]
        start, stop = block_starts[b], block_stops[b]
        group_size = _count_rivers_per_runoff_group(runoff, stop - start)
        for group_start in range(start, stop, group_size):
            group_stop = min(group_start + group_size, stop)
            _prepare_runoff_for_river_group(runoff, group_start, group_stop, scratch)
            for r in range(group_start, group_stop):
                while next_cut < cuts.shape[0] and cut_target[cuts[next_cut]] == r:
                    injected = get_or_open_downstream_inflow_row(r - first, inflow, slot_of, free, top)
                    source = boundary[cuts[next_cut]]
                    for g in range(injected.shape[0]):  # a loop, not injected += source; see the module docstring
                        injected[g] += source[g]
                    next_cut += 1
                up = get_upstream_inflow_row(r - first, inflow, slot_of)
                if r == outlet:
                    down = boundary[block_region[b]]
                    down[:] = 0.0
                else:
                    d = downstream_indices[r]
                    down = get_or_open_downstream_inflow_row(d - first if d >= 0 else -1, inflow, slot_of, free, top)
                row = out_row[r] if renumbered else r
                series = discharge_array[row] if row >= 0 else discarded
                lateral = _get_river_runoff_series(runoff, r, group_start, scratch)
                if dynamic is None:
                    _route_river_with_static_coefficients(
                        q_t, r, coefficients, layout, lateral, up, down, series, work, chain, n_steps, n_substeps
                    )
                else:
                    q_t[r] = _route_river_with_dynamic_coefficients(
                        q_t[r], r, dynamic, lateral, up, down, series, n_steps, n_substeps
                    )
                release_inflow_row_to_pool(r - first, slot_of, free, top)
    return
