"""
Customizable routing kernels in numba.

Terminology:

- The full watershed or region worth of work to be routed is given by the config file. If threading is used, then:
- The region can be split into jobs, one for each thread. To balance work between jobs/threads, then:
- Each job is assigned sub-watersheds, or "blocks" of rivers, so the amount of work for each thread is balanced.

1. The user calls rr.Router.route() once per simulation which calls these relevant functions.

    Router.route()
    ├── Router._execute_routing_channel()             forcing 'channel', called once
    │   └── route_region(router, q_t, discharge_array, None, thread_pool)
    └── Router._execute_routing_forced()              once per file the Runoff's generator yields
        └── route_region(router, q_t, discharge_array, runoff, thread_pool)

2. route_region maps all work to be done to sub pieces, called jobs, that can run concurrently.

    route_region()                                    python
    ├── Network.routing_blocks()                      cached per thread count
    │   └── streams.assign_blocks()                   the blocks of each job, read from the upstreamCount column
    │       └── pack_blocks_into_jobs()               network/_numba_kernels.py, one job per thread
    ├── _check_everything_a_job_reads()
    ├── route_job(job)                                each sub-watershed job, concurrently on thread_pool when given
    └── route_job(main stem job)                      after every other job has finished

3. Route each job. When threading is used, 1 block from every job is routing concurrently.

    route_job()                                       numba.njit
    ├── count_peak_live_inflow_rows()                 sizes the inflow row pool
    ├── first_river_and_span_of_job()
    ├── allocate_inflow_row_pool()
    ├── sort_block_cuts_by_target_river()
    └── for each block, for each of its rivers r in index order:
        ├── get_or_open_downstream_inflow_row()       main stem job: add the block outlet series that drain into r
        ├── get_upstream_inflow_row()                 the row r's upstream rivers added their series into
        ├── get_or_open_downstream_inflow_row()       the row r adds its series into, or boundary at a block outlet
        ├── get_river_forcing(runoff, r)              stage
        ├── transform_runoff(transform, ...)          stage, skipped for the uniform transform
        ├── route_river(method, ...)                  stage
        │   └── add_river_forcing(forcing, ...)       stage, called by the routing method
        └── release_inflow_row_to_pool()

4. There are 4 stages of routing 1 river which can be overloaded for different types of inputs or methods.

    stage              argument type           what the overload does                      registered in
    get_river_forcing  None (channel)          returns None, for channel routing           this module
                       CatchmentRunoffVolumes  returns river r's row of catchment runoff   runoff/bases.py
                       GridCellRunoff          returns river r's RiverGridCells            runoff/bases.py
    transform_runoff   None                    no overload: the uniform transform is skipped
    route_river        StaticMuskingum         routes river r with static coefficients     router/static_muskingum.py
                       DynamicMuskingum        routes river r with dynamic coefficients    router/dynamic_muskingum.py
    add_river_forcing  1D array                adds a catchment runoff series into work    runoff/bases.py
                       RiverGridCells          adds each grid cell's runoff into work      runoff/bases.py

The stages are the functions without a body below. numba picks the overload that matches the type of the argument and
compiles one version of route_job for each combination of argument types it is given. transform_runoff has no
overload yet: the only transform is the uniform one, which is None, so route_job skips that call.

To add a kind of runoff, add overloads of get_river_forcing and add_river_forcing next to its type; to add a routing
method, add an overload of route_river. Nothing else in this module changes. Overloads should not use inline='always'.
"""

from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING, NamedTuple

import numba
import numpy as np
from numba import types
from numba.extending import overload

from ..types import FloatArray, Int32Array, JobBlocks

if TYPE_CHECKING:
    from ..runoff import CatchmentRunoffVolumes, GridCellRunoff
    from .Router import Router

__all__ = [
    # the stages a job takes each river through, implemented next to the types they read
    'get_river_forcing',
    'transform_runoff',
    'add_river_forcing',
    'route_river',
    'is_argument_type',
    # what a job reads, routing a job, and running the jobs of a Router's blocks
    'Layout',
    'STANDARD_LAYOUT',
    'Job',
    'route_job',
    'route_region',
]


class Layout(NamedTuple):
    """
    How each river is routed. Both are empty on a standard network, where ``q_t`` holds one state per river. On a
    stabilized network river r is routed in substeps, the sub-reaches ``reach_indptr[r]:reach_indptr[r + 1]`` in
    series, with one state per sub-reach in ``q_t``, and ``subcycles`` gives the steps of its own each river takes per
    routing step, or is empty when no river is subcycled.
    """

    reach_indptr: np.ndarray  # (n + 1,) int64, or empty
    subcycles: np.ndarray  # (n,) int64, or empty


class Job(NamedTuple):
    """
    The blocks of rivers one thread routes, in order. A block's outlet (-1 for none) hands its unclamped series to
    ``boundary[block_number]`` instead of an inflow row, and each boundary row in ``cut_target`` is injected into its
    target river before that river is routed. ``out_row`` is empty when every river has a discharge row, and otherwise
    gives each river's row, -1 for a synthetic river whose series is written to a scratch row and discarded.
    """

    block_starts: np.ndarray  # (n_job_blocks,) int32
    block_stops: np.ndarray  # (n_job_blocks,) int32
    block_outlet: np.ndarray  # (n_job_blocks,) int32
    block_number: np.ndarray  # (n_job_blocks,) int32
    cut_target: np.ndarray  # (n_blocks,) int32, or empty
    boundary: np.ndarray  # (n_blocks, n_routing + 1) float32
    out_row: np.ndarray  # (n,) int32, or empty


STANDARD_LAYOUT = Layout(reach_indptr=np.zeros(0, dtype=np.int64), subcycles=np.zeros(0, dtype=np.int64))
_NO_CUTS = np.zeros(0, dtype=np.int32)  # a job of sub-watershed blocks injects nothing; only the main stem job does


def get_river_forcing(runoff, r):
    """River r's catchment runoff, in the form add_river_forcing reads, or None for channel routing."""


def transform_runoff(transform, catchment_runoff, r, out):
    """River r's catchment runoff volume in each runoff step after the runoff transform, which may write it into out."""


def add_river_forcing(catchment_runoff, multiplier, n_steps, n_per_step, work):
    """
    Add ``multiplier`` times a river's catchment runoff volume in each of the first ``n_steps`` runoff steps into each
    of the ``n_per_step`` steps of ``work`` that the runoff step spans.
    """


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


@overload(get_river_forcing)
def _channel_routing_has_no_catchment_runoff(runoff, r):
    if isinstance(runoff, types.NoneType):
        return lambda runoff, r: None


@numba.njit(cache=True, nogil=True)
def sort_block_cuts_by_target_river(cut_target):
    """The cuts that drain into a river, as positions into cut_target sorted by target, basin outlets (-1) dropped."""
    order = np.argsort(cut_target, kind='mergesort')
    first = 0
    while first < order.shape[0] and cut_target[order[first]] < 0:
        first += 1
    return order[first:]


@numba.njit(cache=True, nogil=True)
def first_river_and_span_of_job(block_starts, block_stops):
    """The first river a job routes and how many indices its blocks span, which sizes its per-river bookkeeping."""
    first = block_starts.min() if block_starts.shape[0] else 0
    last = block_stops.max() if block_stops.shape[0] else 0
    return first, max(last - first, 0)


@numba.njit(cache=True, nogil=True)
def count_peak_live_inflow_rows(downstream_indices, block_starts, block_stops, block_outlet, cut_target):
    """
    Peak number of inflow rows live at once when a job routes its blocks in order. A block's outlet drains outside
    the job, so it never opens a row.
    """
    first, span = first_river_and_span_of_job(block_starts, block_stops)
    is_open = np.zeros(span, dtype=np.bool_)
    cuts = sort_block_cuts_by_target_river(cut_target)
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


@numba.njit(cache=True, nogil=True)
def route_job(
    q_t, discharge_array, downstream_indices, routing_steps_per_runoff_step, method, layout, runoff, transform, job
):
    """
    Route every river in the blocks of ``job``, river by river, through the stages: its catchment runoff from
    ``runoff``, transformed by ``transform``, and routed with the routing ``method``'s parameters. ``discharge_array``
    is C-order (river, time): each river's series is written into its own row in place.
    """
    n_steps = discharge_array.shape[1]
    n_routing = n_steps * routing_steps_per_runoff_step
    block_starts, block_stops, block_outlet, block_number, cut_target, boundary, out_row = job
    n_rows = count_peak_live_inflow_rows(downstream_indices, block_starts, block_stops, block_outlet, cut_target)
    first, span = first_river_and_span_of_job(block_starts, block_stops)
    inflow_rows, slot_of, free, top = allocate_inflow_row_pool(span, n_rows, n_routing)
    cuts = sort_block_cuts_by_target_river(cut_target)
    next_cut = 0
    transform_scratch = np.empty(0 if transform is None else n_steps, dtype=np.float32)
    expanded = layout.reach_indptr.shape[0] > 0
    most_subcycles = layout.subcycles.max() if layout.subcycles.shape[0] > 0 else 1
    work = np.empty(n_routing * most_subcycles, dtype=np.float32)
    chain = np.empty(n_routing * most_subcycles + 1 if expanded else 0, dtype=np.float32)
    renumbered = out_row.shape[0] > 0
    discarded = np.empty(n_steps if renumbered else 0, dtype=discharge_array.dtype)

    for b in range(block_starts.shape[0]):
        outlet = block_outlet[b]
        for r in range(block_starts[b], block_stops[b]):
            while next_cut < cuts.shape[0] and cut_target[cuts[next_cut]] == r:
                injected = get_or_open_downstream_inflow_row(r - first, inflow_rows, slot_of, free, top)
                block_outlet_series = boundary[cuts[next_cut]]
                for g in range(injected.shape[0]):
                    injected[g] += block_outlet_series[g]
                next_cut += 1
            inflow = get_upstream_inflow_row(r - first, inflow_rows, slot_of)
            if r == outlet:
                downstream_inflow = boundary[block_number[b]]
                downstream_inflow[:] = 0.0
            else:
                d = downstream_indices[r]
                downstream_inflow = get_or_open_downstream_inflow_row(
                    d - first if d >= 0 else -1, inflow_rows, slot_of, free, top
                )
            row = out_row[r] if renumbered else r
            discharge = discharge_array[row] if row >= 0 else discarded
            catchment_runoff = get_river_forcing(runoff, r)
            if transform is not None:
                catchment_runoff = transform_runoff(transform, catchment_runoff, r, transform_scratch)
            route_river(method, layout, q_t, r, catchment_runoff, inflow, downstream_inflow, discharge, work, chain)
            release_inflow_row_to_pool(r - first, slot_of, free, top)
    return


def route_region(
    router: Router,
    q_t: FloatArray,
    discharge_array: FloatArray,
    runoff: CatchmentRunoffVolumes | GridCellRunoff | None,
    thread_pool: ThreadPoolExecutor | None,
) -> None:
    """
    Route ``runoff`` through the whole region with the Router's routing method and layout: every job of sub-watershed
    blocks, concurrently on ``thread_pool`` when given, then the main stem, which consumes the boundary buffer the
    blocks filled. ``runoff`` is None for channel routing, or what the Router's Runoff generator yielded.

    A boundary row holds a block outlet's whole series plus its initial state. The Network packs the blocks into one
    job per thread, which keeps the per-job cost in python, and the buffers each job allocates, to once per thread
    rather than once per block.
    """
    # the Network derives the blocks once per thread count and caches them, so this is a lookup after the first file
    job_blocks, cut_target = router.network.routing_blocks(router.threads if thread_pool is not None else 1)
    _check_everything_a_job_reads(router, q_t, discharge_array, runoff, job_blocks)
    downstream_indices = router.network.downstream_indices
    routing_steps_per_runoff_step = router.num_routing_steps_per_runoff
    n_routing_steps = router.num_runoff_steps * routing_steps_per_runoff_step
    boundary = np.zeros((max(cut_target.shape[0], 1), n_routing_steps + 1), dtype=np.float32)
    synthetic = router.network.synthetic  # each river's discharge row, -1 for a synthetic river that has none
    out_row = _NO_CUTS if synthetic is None else np.where(synthetic, -1, np.cumsum(~synthetic) - 1).astype(np.int32)

    method, layout, transform = router.routing_parameters, router.layout, None  # None is the uniform transform

    def route_one_job(blocks: JobBlocks, cut_target_of_job: Int32Array = _NO_CUTS) -> None:
        route_job(
            q_t,
            discharge_array,
            downstream_indices,
            routing_steps_per_runoff_step,
            method,
            layout,
            runoff,
            transform,
            Job(*blocks, cut_target_of_job, boundary, out_row),
        )

    *sub_watershed_jobs, main_stem = job_blocks
    if thread_pool is None:
        for blocks in sub_watershed_jobs:
            route_one_job(blocks)
    else:
        list(thread_pool.map(route_one_job, sub_watershed_jobs))  # list() so a worker exception raises

    # the one barrier of the simulation: every block has finished before the main stem consumes its buffer
    route_one_job(main_stem, cut_target)
    return


def _check_everything_a_job_reads(
    router: Router,
    q_t: FloatArray,
    discharge_array: FloatArray,
    runoff: CatchmentRunoffVolumes | GridCellRunoff | None,
    job_blocks: tuple[JobBlocks, ...],
) -> None:
    """Last chance to validate before routing. The numba kernels will not check or raise useful errors."""
    # each river's series is written into its own contiguous row, which is the layout the jobs solve in.
    # Output is never copied, so any other layout is refused rather than transposed.
    if not discharge_array.flags.c_contiguous:
        raise ValueError('discharge_array must be a C-order (river, time) array')
    n_rivers = router.network.river_ids.shape[0]
    n_steps = router.num_runoff_steps
    reach_indptr, subcycles = router.layout
    if reach_indptr.shape[0] and (
        reach_indptr.shape != (n_rivers + 1,) or reach_indptr[0] != 0 or np.any(np.diff(reach_indptr) < 1)
    ):
        raise ValueError(f'reach_indptr must be ({n_rivers + 1},) offsets from 0 with at least one reach per river')
    if subcycles.shape[0] and (subcycles.shape != (n_rivers,) or subcycles.min() < 1):
        raise ValueError(f'subcycles must be ({n_rivers},) counts of at least 1, or empty')
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
        runoff.check(router.network.river_ids, n_steps)

    # sorted by start, the blocks tile [0, n_rivers) exactly when each ends where the next begins: no gap, no overlap
    starts = np.concatenate([blocks[0] for blocks in job_blocks])
    stops = np.concatenate([blocks[1] for blocks in job_blocks])
    order = np.argsort(starts)
    starts, stops = starts[order], stops[order]
    if np.any(stops <= starts) or starts[0] != 0 or stops[-1] != n_rivers or np.any(stops[:-1] != starts[1:]):
        raise ValueError(f'the blocks must cover every one of the {n_rivers} rivers exactly once')
    return
