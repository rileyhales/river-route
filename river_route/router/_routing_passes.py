"""
How a network is routed. A pass routes the rivers of a schedule's blocks one at a time, each river's whole series before
the next, and takes every river through three stages: finding its catchment runoff, transforming that runoff, and
routing it with the routing method. The stages are the functions without a body below. Each is implemented by a numba
overload registered in the module that defines the type of argument it reads, so numba compiles one version of
route_scheduled_rivers for each combination of argument types it is given:

    runoff     None for channel routing (here), CatchmentRunoffVolumes (runoff/CatchmentRunoff.py), or
               GridCellRunoff (runoff/GaussianGridRunoff.py)
    transform  None for the uniform transform, which changes nothing, so the stage is removed at compile time
    method     StaticMuskingum (router/static_muskingum.py) or DynamicMuskingum (router/dynamic_muskingum.py)

An overload returns None for types it does not implement, so numba tries the next one. The overloads do not use
inline='always': inlining into route_scheduled_rivers at the numba IR level (numba 0.67) made the whole-array in-place
add of a region's outlet series compile as a new array rather than into the inflow row, which dropped that series from
the main stem of threaded routing. numba removes a branch on ``argument is None`` only when the argument is None, so a
stage chooses between types by overload, and a branch on None only skips a stage, as for the uniform transform: an
identity overload for it measured 2% slower on the Amazon, since each call passes and returns arrays.

route_network runs the passes over a Router's schedule. The recurrence, the regions, and the inflow row pool are
described in docs/references/kernels.md under "How routing works". Inputs are NaN free: GaussianGridRunoff and
CatchmentRunoff replace NaN with zero when they read the runoff.
"""

import heapq
from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING, NamedTuple

import numba
import numpy as np
from numba import types
from numba.extending import overload

from ..types import FloatArray

if TYPE_CHECKING:
    from ..runoff import CatchmentRunoffVolumes, GridCellRunoff
    from .Router import Router

__all__ = [
    # the stages a routing pass takes each river through, implemented next to the types they read
    'count_rivers_prepared_together',
    'prepare_runoff_of_rivers',
    'get_river_catchment_runoff',
    'transform_catchment_runoff',
    'route_river',
    'is_argument_type',
    # what a pass reads, the pass, and running the passes over a Router's schedule
    'Layout',
    'STANDARD_LAYOUT',
    'Schedule',
    'route_scheduled_rivers',
    'route_network',
]


################################################
# The stages of routing one river, implemented next to the types they read
################################################


def count_rivers_prepared_together(runoff):
    """How many rivers' runoff is prepared at once in scratch rows before they are routed, 0 when read in place."""


def prepare_runoff_of_rivers(runoff, first_river, stop_river, scratch):
    """Prepare the runoff of rivers first_river:stop_river in the scratch rows before they are routed."""


def get_river_catchment_runoff(runoff, r, first_river, scratch):
    """River r's catchment runoff volume in each runoff step, or None for channel routing. first_river is the first
    river whose runoff was prepared in scratch."""


def transform_catchment_runoff(transform, catchment_runoff, r, out):
    """River r's catchment runoff volume in each runoff step after the runoff transform, which may write it into out."""


def route_river(method, layout, q_t, r, catchment_runoff, inflow, downstream_inflow, discharge, work, chain):
    """
    Route river r's whole series with the parameters of a routing method. ``inflow`` is the summed discharge of its
    upstreams before the first routing step and at every step. The river adds its own series, in the same form, into
    ``downstream_inflow``, writes its discharge averaged over each runoff step into ``discharge``, and updates its
    state in ``q_t``. ``work`` and ``chain`` are scratch sized for the ``layout``.
    """


def is_argument_type(numba_type, kind: type) -> bool:
    """Whether a stage argument's numba type is the NamedTuple class ``kind``."""
    return getattr(numba_type, 'instance_class', None) is kind


@overload(count_rivers_prepared_together)
def _channel_routing_prepares_no_runoff(runoff):
    if isinstance(runoff, types.NoneType):
        return lambda runoff: 0


@overload(prepare_runoff_of_rivers)
def _channel_routing_has_no_runoff_to_prepare(runoff, first_river, stop_river, scratch):
    if isinstance(runoff, types.NoneType):
        return lambda runoff, first_river, stop_river, scratch: None


@overload(get_river_catchment_runoff)
def _channel_routing_has_no_catchment_runoff(runoff, r, first_river, scratch):
    if isinstance(runoff, types.NoneType):
        return lambda runoff, r, first_river, scratch: None


################################################
# What a pass reads
################################################


class Layout(NamedTuple):
    """
    How each river is routed. Both are empty on a standard network, where ``q_t`` holds one state per river. On a
    stabilized network river r is the sub-reaches ``reach_indptr[r]:reach_indptr[r + 1]`` in series, with one state
    per sub-reach in ``q_t``, and ``substeps`` gives the steps each river takes per routing step, or is empty when no
    river is sub-cycled.
    """

    reach_indptr: np.ndarray  # (n + 1,) int64, or empty
    substeps: np.ndarray  # (n,) int64, or empty


STANDARD_LAYOUT = Layout(reach_indptr=np.zeros(0, dtype=np.int64), substeps=np.zeros(0, dtype=np.int64))


class Schedule(NamedTuple):
    """
    The blocks of rivers one pass routes, in order. A block's outlet (-1 for none) hands its unclamped series to
    ``boundary[block_region]`` instead of an inflow row, and each boundary row in ``cut_target`` is injected into its
    target river before that river is routed. ``out_row`` is empty when every river has a discharge row, and otherwise
    gives each river's row, -1 for a synthetic river whose series is written to a scratch row and discarded.
    """

    block_starts: np.ndarray  # (n_blocks,) int32
    block_stops: np.ndarray  # (n_blocks,) int32
    block_outlet: np.ndarray  # (n_blocks,) int32
    block_region: np.ndarray  # (n_blocks,) int32
    cut_target: np.ndarray  # (n_regions,) int32, or empty
    boundary: np.ndarray  # (n_regions, n_routing + 1) float32
    out_row: np.ndarray  # (n,) int32, or empty


################################################
# The pool of inflow rows a pass accumulates upstream series into. A river's row is taken from the pool when its first
# upstream is routed, or when a region's buffered outlet series is injected into it, and returned once the river itself
# is routed. Every river index into slot_of is relative to the pass's first river.
################################################


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
            peak = max(peak, live)
            if is_open[r - first]:
                live -= 1
    return peak


@numba.njit(cache=True, nogil=True)
def allocate_inflow_row_pool(span, n_rows, n_routing_steps):
    """The inflow rows, the row each river holds (-1 for none), and the stack of free rows with its top index."""
    inflow_rows = np.empty((n_rows + 2, n_routing_steps + 1), dtype=np.float32)
    inflow_rows[n_rows, :] = 0.0  # headwaters read this row as their inflow
    slot_of = np.full(span, -1, dtype=np.int64)
    free = np.arange(n_rows)
    top = np.full(1, n_rows, dtype=np.int64)
    return inflow_rows, slot_of, free, top


@numba.njit(cache=True, nogil=True)
def get_upstream_inflow_row(r, inflow_rows, slot_of):
    """The inflow row river r reads: the one its upstreams filled, or the zero row for a headwater."""
    s = slot_of[r]
    return inflow_rows[s] if s >= 0 else inflow_rows[inflow_rows.shape[0] - 2]


@numba.njit(cache=True, nogil=True)
def get_or_open_downstream_inflow_row(d, inflow_rows, slot_of, free, top):
    """The inflow row river d accumulates into, taken from the pool and zeroed on first use; the sink if d < 0."""
    if d < 0:
        return inflow_rows[inflow_rows.shape[0] - 1]
    s = slot_of[d]
    if s < 0:
        top[0] -= 1
        s = free[top[0]]
        slot_of[d] = s
        inflow_rows[s, :] = 0.0
    return inflow_rows[s]


@numba.njit(cache=True, nogil=True)
def release_inflow_row_to_pool(r, slot_of, free, top):
    """River r has been routed, so the inflow row it read goes back to the pool."""
    s = slot_of[r]
    if s >= 0:
        free[top[0]] = s
        top[0] += 1
        slot_of[r] = -1


################################################
# The pass
################################################


@numba.njit(cache=True, nogil=True)
def route_scheduled_rivers(
    q_t, discharge_array, downstream_indices, n_substeps, method, layout, runoff, transform, schedule
):
    """
    Route every river in the blocks of ``schedule``, river by river, through the stages: its catchment runoff from
    ``runoff``, transformed by ``transform``, and routed with the routing ``method``'s parameters. ``discharge_array``
    is C-order (river, time): each river's series is written into its own row in place.

    Rivers are routed in groups whose runoff is prepared together first. Neighboring catchments share grid cells, so
    aggregating gridded runoff for a group back to back reads each shared cell's series while it is still in cache.
    """
    n_steps = discharge_array.shape[1]
    n_routing = n_steps * n_substeps
    block_starts, block_stops, block_outlet, block_region, cut_target, boundary, out_row = schedule
    n_rows = count_peak_live_inflow_rows(downstream_indices, block_starts, block_stops, block_outlet, cut_target)
    first, span = first_river_and_span_of_pass(block_starts, block_stops)
    inflow_rows, slot_of, free, top = allocate_inflow_row_pool(span, n_rows, n_routing)
    cuts = sort_region_cuts_by_target_river(cut_target)
    next_cut = 0
    rivers_prepared_together = count_rivers_prepared_together(runoff)
    runoff_scratch = np.zeros((rivers_prepared_together, n_steps), dtype=np.float32)
    transform_scratch = np.empty(0 if transform is None else n_steps, dtype=np.float32)
    expanded = layout.reach_indptr.shape[0] > 0
    most = layout.substeps.max() if layout.substeps.shape[0] > 0 else 1
    work = np.empty(n_routing * most, dtype=np.float32)
    chain = np.empty(n_routing * most + 1 if expanded else 0, dtype=np.float32)
    renumbered = out_row.shape[0] > 0
    discarded = np.empty(n_steps if renumbered else 0, dtype=discharge_array.dtype)

    for b in range(block_starts.shape[0]):
        outlet = block_outlet[b]
        start, stop = block_starts[b], block_stops[b]
        group_size = rivers_prepared_together or max(stop - start, 1)
        for group_start in range(start, stop, group_size):
            group_stop = min(group_start + group_size, stop)
            prepare_runoff_of_rivers(runoff, group_start, group_stop, runoff_scratch)
            for r in range(group_start, group_stop):
                while next_cut < cuts.shape[0] and cut_target[cuts[next_cut]] == r:
                    injected = get_or_open_downstream_inflow_row(r - first, inflow_rows, slot_of, free, top)
                    region_outlet_series = boundary[cuts[next_cut]]
                    for g in range(injected.shape[0]):  # a loop, not +=; see the module docstring
                        injected[g] += region_outlet_series[g]
                    next_cut += 1
                inflow = get_upstream_inflow_row(r - first, inflow_rows, slot_of)
                if r == outlet:
                    downstream_inflow = boundary[block_region[b]]
                    downstream_inflow[:] = 0.0
                else:
                    d = downstream_indices[r]
                    downstream_inflow = get_or_open_downstream_inflow_row(
                        d - first if d >= 0 else -1, inflow_rows, slot_of, free, top
                    )
                row = out_row[r] if renumbered else r
                discharge = discharge_array[row] if row >= 0 else discarded
                catchment_runoff = get_river_catchment_runoff(runoff, r, group_start, runoff_scratch)
                if transform is not None:
                    catchment_runoff = transform_catchment_runoff(transform, catchment_runoff, r, transform_scratch)
                route_river(method, layout, q_t, r, catchment_runoff, inflow, downstream_inflow, discharge, work, chain)
                release_inflow_row_to_pool(r - first, slot_of, free, top)
    return


################################################
# Running the passes over a Router's schedule
################################################

_NO_CUTS = np.zeros(0, dtype=np.int32)  # a region pass injects nothing; only the main stem does


def route_network(
    router: Router,
    q_t: FloatArray,
    discharge_array: FloatArray,
    runoff: CatchmentRunoffVolumes | GridCellRunoff | None,
    thread_pool: ThreadPoolExecutor | None,
) -> None:
    """
    Route ``runoff`` through the whole network with the Router's routing method and layout: every region, concurrently
    on ``thread_pool`` when given, then the main stem, which consumes the boundary buffer the regions filled. ``runoff``
    is None for channel routing, or what the Router's Runoff generator yielded.

    A boundary row holds a region outlet's whole series plus its initial state. One pass takes many regions as blocks,
    each with its own outlet, so the regions are packed longest first into one pass per thread. That keeps the
    per-pass cost in python, and the buffers each pass allocates, to once per thread rather than once per region.
    """
    _check_everything_a_pass_reads(router, q_t, discharge_array, runoff)
    downstream_indices = router.network.downstream_indices
    n_substeps = router.num_routing_steps_per_runoff
    n_regions = len(router.routing_jobs) - 1
    boundary = np.zeros((max(n_regions, 1), router.num_runoff_steps * n_substeps + 1), dtype=np.float32)
    synthetic = router.network.synthetic  # each river's discharge row, -1 for a synthetic river that has none
    out_row = _NO_CUTS if synthetic is None else np.where(synthetic, -1, np.cumsum(~synthetic) - 1).astype(np.int32)

    method, layout, transform = router.routing_parameters, router.layout, None  # None is the uniform transform

    def route_pass(block_starts, block_stops, block_outlet, block_region, cut_target=_NO_CUTS) -> None:
        schedule = Schedule(block_starts, block_stops, block_outlet, block_region, cut_target, boundary, out_row)
        route_scheduled_rivers(
            q_t, discharge_array, downstream_indices, n_substeps, method, layout, runoff, transform, schedule
        )

    n_passes = router.threads if thread_pool is not None else 1
    region_passes = _pack_regions_into_passes(router.routing_jobs[:-1], n_passes)
    if thread_pool is None:
        for blocks in region_passes:
            route_pass(*blocks)
    else:
        list(thread_pool.map(lambda blocks: route_pass(*blocks), region_passes))  # list() so a worker exception raises

    # the one barrier of the simulation: every region has finished before the main stem consumes its buffer
    stem_starts, stem_stops, _, _ = router.routing_jobs[-1]
    no_outlet = np.full(stem_starts.shape[0], -1, dtype=np.int32)
    route_pass(stem_starts, stem_stops, no_outlet, np.zeros(stem_starts.shape[0], dtype=np.int32), router.cut_target)
    return


def _pack_regions_into_passes(regions: tuple, n_passes: int) -> list[tuple[np.ndarray, ...]]:
    """
    Pack region jobs, longest first, onto the least loaded of ``n_passes`` passes. Each pass lists its regions as
    blocks in index order: (block_starts, block_stops, block_outlet, block_region).
    """
    loads = [(0, p) for p in range(max(1, n_passes))]
    members: list[list] = [[] for _ in loads]
    for starts, stops, outlet, region in regions:
        load, p = heapq.heappop(loads)
        members[p].append((int(starts[0]), int(stops[0]), outlet, region))
        heapq.heappush(loads, (load + int(stops[0] - starts[0]), p))
    passes = []
    for jobs in sorted(members, key=lambda jobs: -sum(stop - start for start, stop, _, _ in jobs)):
        if jobs:
            jobs.sort()
            passes.append(tuple(np.array(v, dtype=np.int32) for v in zip(*jobs, strict=True)))
    return passes


def _check_everything_a_pass_reads(
    router: Router, q_t: FloatArray, discharge_array: FloatArray, runoff: CatchmentRunoffVolumes | GridCellRunoff | None
) -> None:
    """
    Check every array a pass reads, because this is the last point where a wrong one can be caught. The passes are
    compiled with ``numba.njit`` and therefore do no bounds checking: an array shorter than the network is read and
    written past its end rather than raising IndexError. A malformed schedule is worse: a cut_target inside a region,
    or blocks that overlap, make one thread write into rivers another thread is routing, which corrupts results
    without raising and without reproducing reliably.

    Raises:
        ValueError: if any array does not match the network, the number of routing steps, or the schedule
    """
    # each river's series is written into its own contiguous row, which is the layout the passes solve in. Output is
    # never copied, so any other layout is refused rather than transposed.
    if not discharge_array.flags.c_contiguous:
        raise ValueError('discharge_array must be a C-order (river, time) array')
    n_rivers = router.network.river_ids.shape[0]
    n_steps = router.num_runoff_steps
    reach_indptr, substeps = router.layout
    if reach_indptr.shape[0] and (
        reach_indptr.shape != (n_rivers + 1,) or reach_indptr[0] != 0 or np.any(np.diff(reach_indptr) < 1)
    ):
        raise ValueError(f'reach_indptr must be ({n_rivers + 1},) offsets from 0 with at least one reach per river')
    if substeps.shape[0] and (substeps.shape != (n_rivers,) or substeps.min() < 1):
        raise ValueError(f'substeps must be ({n_rivers},) counts of at least 1, or empty')
    n_states = int(reach_indptr[-1]) if reach_indptr.shape[0] else n_rivers  # one state per reach
    synthetic = router.network.synthetic
    n_out = n_rivers if synthetic is None else int(np.count_nonzero(~synthetic))  # synthetic rivers have no row
    expected_shapes = {'q_t': ((n_states,), q_t.shape), 'discharge_array': ((n_out, n_steps), discharge_array.shape)}
    for name, (want, got) in expected_shapes.items():
        if got != want:
            raise ValueError(f'{name} has shape {got}, expected {want} for {n_rivers} rivers and {n_steps} steps')
    per_river = {'downstream_indices': router.network.downstream_indices, **router.routing_parameters._asdict()}
    for name, array in per_river.items():
        if isinstance(array, np.ndarray) and array.shape != (n_rivers,):
            raise ValueError(f'{name} has shape {array.shape}, expected ({n_rivers},)')
    if runoff is not None:
        runoff.check(n_rivers, n_steps)

    covered = 0
    for starts, stops, _outlet, _region in router.routing_jobs:
        if starts.shape != stops.shape:
            raise ValueError(f'block starts/stops shapes differ: {starts.shape} vs {stops.shape}')
        if starts.size and (starts.min() < 0 or stops.max() > n_rivers or np.any(stops <= starts)):
            raise ValueError(f'blocks must be non-empty ranges within [0, {n_rivers})')
        covered += int((stops - starts).sum())
    if covered != n_rivers:
        raise ValueError(f'region schedule covers {covered} rivers, expected every one of {n_rivers} exactly once')
    targets = router.cut_target[router.cut_target >= 0]
    if targets.size and targets.max() >= n_rivers:
        raise ValueError(f'every cut_target must fall inside [0, {n_rivers}); found {int(targets.max())}')
    return
