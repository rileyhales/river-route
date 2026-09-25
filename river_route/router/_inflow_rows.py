"""
The pool of inflow rows a routing pass accumulates upstream series into. A river's row is taken from the pool when its
first upstream is routed, or when a region's buffered outlet series is injected into it, and returned once the river
itself is routed. Every river index into ``slot_of`` is relative to the pass's first river.
"""

import numba
import numpy as np

__all__ = [
    'sort_region_cuts_by_target_river',
    'first_river_and_span_of_pass',
    'count_peak_live_inflow_rows',
    'allocate_inflow_row_pool',
    'get_upstream_inflow_row',
    'get_or_open_downstream_inflow_row',
    'release_inflow_row_to_pool',
]


@numba.njit(cache=True, nogil=True)
def sort_region_cuts_by_target_river(cut_target):
    """The cuts that drain into a river, as positions into cut_target sorted by target, basin outlets (-1) dropped."""
    order = np.argsort(cut_target, kind='mergesort')
    first = 0
    while first < order.shape[0] and cut_target[order[first]] < 0:
        first += 1
    return order[first:]


@numba.njit(cache=True, nogil=True)
def first_river_and_span_of_pass(block_starts, block_stops):
    """The first river a pass routes and how many indices its blocks span, which sizes its per-river bookkeeping."""
    first = block_starts.min() if block_starts.shape[0] else 0
    last = block_stops.max() if block_stops.shape[0] else 0
    return first, max(last - first, 0)


@numba.njit(cache=True, nogil=True)
def count_peak_live_inflow_rows(downstream_indices, block_starts, block_stops, block_outlet, cut_target):
    """
    Peak number of inflow rows live at once when a pass routes its blocks in order. A block's outlet drains outside
    the pass, so it never opens a row.
    """
    first, span = first_river_and_span_of_pass(block_starts, block_stops)
    is_open = np.zeros(span, dtype=np.bool_)
    cuts = sort_region_cuts_by_target_river(cut_target)
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
def allocate_inflow_row_pool(span, n_rows, n_routing_steps):
    """The inflow rows, the row each river holds (-1 for none), and the stack of free rows with its top index."""
    inflow = np.empty((n_rows + 2, n_routing_steps + 1), dtype=np.float32)
    inflow[n_rows, :] = 0.0  # headwaters read this row as their inflow
    slot_of = np.full(span, -1, dtype=np.int64)
    free = np.arange(n_rows)
    top = np.full(1, n_rows, dtype=np.int64)
    return inflow, slot_of, free, top


@numba.njit(cache=True, nogil=True)
def get_upstream_inflow_row(r, inflow, slot_of):
    """The inflow row river r reads: the one its upstreams filled, or the zero row for a headwater."""
    s = slot_of[r]
    return inflow[s] if s >= 0 else inflow[inflow.shape[0] - 2]


@numba.njit(cache=True, nogil=True)
def get_or_open_downstream_inflow_row(d, inflow, slot_of, free, top):
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
def release_inflow_row_to_pool(r, slot_of, free, top):
    """River r has been routed, so the inflow row it read goes back to the pool."""
    s = slot_of[r]
    if s >= 0:
        free[top[0]] = s
        top[0] += 1
        slot_of[r] = -1
