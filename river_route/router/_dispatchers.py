from __future__ import annotations

import heapq
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING

import numpy as np

from ..runoff import CellRunoff
from ..types import FloatArray
from ._routing_kernel_arguments import CatchmentByRiver, DynamicCoefficients, Layout, Schedule, StaticCoefficients
from ._routing_kernels import route_scheduled_rivers

if TYPE_CHECKING:
    from ..configs import Configs
    from .Router import Router

__all__ = ['DISPATCHERS', 'resolve_dispatcher', 'dispatch_channel', 'dispatch_catchment', 'dispatch_gaussian_grid']

# One dispatcher per kernel, registered under every combination of configs it can route. A key names the kernel's
# abilities: (coefficients, forcing, transform, runoff_type, network_type). Channel routing reads no runoff, so its
# transform and runoff_type are None. A stabilized network is described by the reach_indptr and substeps arrays, so a
# kernel that routes one registers for both network types rather than having a second kernel.
type DispatchKey = tuple[str, str, str | None, str | None, str]
DISPATCHERS: dict[DispatchKey, Callable[..., None]] = {}


_NO_CUTS: FloatArray = np.zeros(0, dtype=np.int32)  # a region pass injects nothing; only the main stem does


def _describe(key: DispatchKey) -> str:
    coefficients, forcing, transform, runoff_type, network_type = key
    return (
        f'coefficients={coefficients},'
        f'forcing={forcing},'
        f'transform={transform},'
        f'runoff_type={runoff_type},'
        f'network_type={network_type}'
    )


def register(coefficients: str, forcing: str, transform: str | None, runoff_type: str | None, network_type: str):
    """Register the decorated dispatcher for one combination of configs."""
    key = (coefficients, forcing, transform, runoff_type, network_type)

    def decorate(dispatcher: Callable[..., None]) -> Callable[..., None]:
        if key in DISPATCHERS:
            raise ValueError(f'Duplicate dispatcher registration for {_describe(key)}')
        DISPATCHERS[key] = dispatcher
        return dispatcher

    return decorate


def resolve_dispatcher(configs: Configs) -> Callable[..., None]:
    """The dispatcher for the configs, or NotImplementedError naming the combination and listing those that exist."""
    channel = configs.forcing == 'channel'
    key = (
        configs.coefficients,
        configs.forcing,
        None if channel else configs.transform,
        None if channel else configs.runoff_type,
        configs.network_type,
    )
    dispatcher = DISPATCHERS.get(key)
    if dispatcher is None:
        available = '\n  '.join(sorted(_describe(k) for k in DISPATCHERS))
        raise NotImplementedError(f'No routing kernel yet for {_describe(key)}.\nImplemented:\n  {available}')
    return dispatcher


def _check_kernel_arrays(
    router: Router, q_t: FloatArray, discharge_array: FloatArray, catchment_runoff: FloatArray | None
) -> None:
    """
    Check the shapes of every array handed to a kernel.

    The kernels are compiled with ``numba.njit`` and therefore do no bounds checking: an array that is shorter
    than ``n_rivers`` is read and written past its end rather than raising IndexError. This is the last point
    where that can be caught, so it is checked here for every routing pass.

    Raises:
        ValueError: if any array does not match the network size or the number of routing steps
    """
    n_rivers = router.network.river_ids.shape[0]
    n_steps = router.num_runoff_steps
    reach_indptr = router.reach_indptr
    if reach_indptr.shape[0] and (
        reach_indptr.shape != (n_rivers + 1,) or reach_indptr[0] != 0 or np.any(np.diff(reach_indptr) < 1)
    ):
        raise ValueError(f'reach_indptr must be ({n_rivers + 1},) offsets from 0 with at least one reach per river')
    n_states = int(reach_indptr[-1]) if reach_indptr.shape[0] else n_rivers  # one state per reach
    substeps = router.substeps
    if substeps.shape[0] and (substeps.shape != (n_rivers,) or substeps.min() < 1):
        raise ValueError(f'substeps must be ({n_rivers},) counts of at least 1, or empty')
    synthetic = router.network.synthetic
    n_out = n_rivers if synthetic is None else int(np.count_nonzero(~synthetic))  # synthetic rivers have no row
    expected = {'q_t': ((n_states,), q_t.shape), 'discharge_array': ((n_out, n_steps), discharge_array.shape)}
    if catchment_runoff is not None:
        expected['catchment_runoff'] = ((n_rivers, n_steps), catchment_runoff.shape)
    for name, (want, got) in expected.items():
        if got != want:
            raise ValueError(f'{name} has shape {got}, expected {want} for {n_rivers} rivers and {n_steps} steps')

    per_river = {
        'downstream_indices': router.network.downstream_indices,
        'c1': getattr(router, 'c1', None),
        'c2': getattr(router, 'c2', None),
        'c3': getattr(router, 'c3', None),
        'c4_dt': getattr(router, 'c4_dt', None),
        'dynamicAlpha': router.network.dynamicAlpha,
        'dynamicBeta': router.network.dynamicBeta,
        'x': router.network.x,
    }
    for name, array in per_river.items():
        if array is not None and array.shape != (n_rivers,):
            raise ValueError(f'{name} has shape {array.shape}, expected ({n_rivers},)')

    _check_partition_arrays(router, n_rivers)
    return


def _check_partition_arrays(router: Router, n_rivers: int) -> None:
    """
    Check the region schedule before it reaches a kernel, for the same reason as _check_kernel_arrays.

    A malformed schedule is worse than a wrong-sized array: a cut_target inside a region, or blocks that overlap,
    makes one thread write into a range another thread is sweeping. That corrupts results without raising and
    without reproducing reliably. None of it is checkable once inside the njit kernel.
    """
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
    if targets.size and (targets.min() < 0 or targets.max() >= n_rivers):
        raise ValueError(f'every cut_target must fall inside [0, {n_rivers}); found {int(targets.max())}')
    return


@register('static', 'channel', None, None, 'standard')
@register('static', 'channel', None, None, 'stabilized')
def dispatch_channel(
    router: Router, q_t: FloatArray, discharge_array: FloatArray, thread_pool: ThreadPoolExecutor | None = None
) -> None:
    """
    Route the initial channel state through the whole network with no runoff, concurrently across regions on
    ``thread_pool`` when given. The schedule is built once by the Router and only consumed here, so nothing about the
    partition is recomputed per call.
    """
    _check_kernel_arrays(router, q_t=q_t, discharge_array=discharge_array, catchment_runoff=None)
    _route(router, q_t, discharge_array, None, thread_pool)
    return


@register('static', 'runoff', 'uniform', 'catchment', 'standard')
@register('static', 'runoff', 'uniform', 'catchment', 'stabilized')
@register('dynamic', 'runoff', 'uniform', 'catchment', 'standard')
def dispatch_catchment(
    router: Router,
    q_t: FloatArray,
    discharge_array: FloatArray,
    catchment_runoff: FloatArray,
    thread_pool: ThreadPoolExecutor | None = None,
) -> None:
    """
    Route catchment runoff volumes, a (river, time) array, through the whole network with the uniform transform,
    concurrently across regions on ``thread_pool`` when given. Each river's row is read in place, so the rows must be
    contiguous; a row view of a longer file, as the Router makes to route part of one, is. Static or dynamic
    coefficients follow the ``coefficients`` config.
    """
    catchment_runoff = catchment_runoff.astype(np.float32, copy=False)
    if catchment_runoff.ndim != 2 or catchment_runoff.strides[1] != catchment_runoff.itemsize:
        raise ValueError('catchment_runoff must be a (river, time) array with each river row contiguous')
    _check_kernel_arrays(router, q_t=q_t, discharge_array=discharge_array, catchment_runoff=catchment_runoff)
    _route(router, q_t, discharge_array, CatchmentByRiver(catchment_runoff), thread_pool)
    return


def _route(
    router: Router,
    q_t: FloatArray,
    discharge_array: FloatArray,
    runoff: CatchmentByRiver | CellRunoff | None,
    thread_pool: ThreadPoolExecutor | None,
) -> None:
    """
    Run ``route_scheduled_rivers`` over the Router's schedule with the Router's coefficients, network layout, and
    the given runoff, None for channel routing. The type of the runoff chooses how the pass reads it.
    """
    # each river's series is written into its own contiguous row, which is the layout the kernels solve in. Output
    # is never copied, so any other layout is refused rather than transposed.
    if not discharge_array.flags.c_contiguous:
        raise ValueError('discharge_array must be a C-order (river, time) array')
    coefficients, dynamic = None, None
    if router.configs.coefficients == 'dynamic':
        dynamic = DynamicCoefficients(
            alpha=router.network.dynamicAlpha,
            beta=router.network.dynamicBeta,
            x=router.network.x,
            dt_routing=np.float32(router.dt_routing),
            inv_dt_runoff=np.float32(1.0 / router.dt_runoff),
        )
    else:
        coefficients = StaticCoefficients(c1=router.c1, c2=router.c2, c3=router.c3, c4_dt=router.c4_dt)
    layout = Layout(reach_indptr=router.reach_indptr, substeps=router.substeps)
    downstream_indices = router.network.downstream_indices
    n_substeps = router.num_routing_steps_per_runoff

    def run(schedule: Schedule) -> None:
        route_scheduled_rivers(
            q_t, discharge_array, downstream_indices, n_substeps, coefficients, dynamic, layout, runoff, schedule
        )

    _run_schedule(router, run, thread_pool)
    return


def _run_schedule(router: Router, run: Callable[[Schedule], None], thread_pool: ThreadPoolExecutor | None) -> None:
    """
    Run the Router's schedule: every region, concurrently on ``thread_pool`` when given, then the main stem, which
    consumes the boundary buffer the regions filled.

    A boundary row holds a region outlet's whole series plus its initial state. One pass takes many regions as
    blocks, each with its own outlet, so the regions are packed longest first into one pass per thread. That keeps
    the per-pass cost in python, and the buffers each pass allocates, to once per thread rather than once per region.
    """
    n_regions = len(router.routing_jobs) - 1
    n_routing_steps = router.num_runoff_steps * router.num_routing_steps_per_runoff
    boundary = np.zeros((max(n_regions, 1), n_routing_steps + 1), dtype=np.float32)
    regions = router.routing_jobs[:-1]  # already ordered longest first by the Router
    stem_starts, stem_stops, _, _ = router.routing_jobs[-1]
    synthetic = router.network.synthetic  # each river's discharge row, -1 for a synthetic river that has none
    out_row = _NO_CUTS if synthetic is None else np.where(synthetic, -1, np.cumsum(~synthetic) - 1).astype(np.int32)
    passes = [
        Schedule(*blocks, cut_target=_NO_CUTS, boundary=boundary, out_row=out_row)
        for blocks in _pack_regions(regions, router.threads if thread_pool is not None else 1)
    ]
    if passes:
        if thread_pool is None:
            for schedule in passes:
                run(schedule)
        else:
            list(thread_pool.map(run, passes))  # list() so a worker exception propagates

    # the one barrier of the simulation: every region has finished before the main stem consumes its buffer
    run(
        Schedule(
            block_starts=stem_starts,
            block_stops=stem_stops,
            block_outlet=np.full(stem_starts.shape[0], -1, dtype=np.int32),
            block_region=np.zeros(stem_starts.shape[0], dtype=np.int32),
            cut_target=router.cut_target,
            boundary=boundary,
            out_row=out_row,
        )
    )
    return


def _pack_regions(regions: tuple, n_bins: int) -> list[tuple[np.ndarray, ...]]:
    """
    Pack region jobs, longest first, onto the least loaded of ``n_bins`` passes. Each pass lists its regions as
    blocks in index order: (block_starts, block_stops, block_outlet, block_region).
    """
    loads = [(0, b) for b in range(max(1, n_bins))]
    members: list[list] = [[] for _ in loads]
    for starts, stops, outlet, region in regions:
        load, b = heapq.heappop(loads)
        members[b].append((int(starts[0]), int(stops[0]), outlet, region))
        heapq.heappush(loads, (load + int(stops[0] - starts[0]), b))
    passes = []
    for jobs in sorted(members, key=lambda jobs: -sum(stop - start for start, stop, _, _ in jobs)):
        if not jobs:
            continue
        jobs.sort()
        passes.append(tuple(np.array(v, dtype=np.int32) for v in zip(*jobs, strict=True)))
    return passes


@register('static', 'runoff', 'uniform', 'gaussian_grid', 'standard')
@register('static', 'runoff', 'uniform', 'gaussian_grid', 'stabilized')
def dispatch_gaussian_grid(
    router: Router,
    q_t: FloatArray,
    discharge_array: FloatArray,
    runoff: CellRunoff | FloatArray,
    thread_pool: ThreadPoolExecutor | None = None,
) -> None:
    """
    Aggregate gaussian grid runoff onto the rivers and route it in the same sweep, which never builds a catchment
    runoff array, concurrently across regions on ``thread_pool`` when given. Static coefficients with the uniform
    transform. A file whose timesteps were resampled arrives already aggregated, as a catchment runoff array, and is
    routed by ``dispatch_catchment``.
    """
    if isinstance(runoff, np.ndarray):
        dispatch_catchment(router, q_t, discharge_array, runoff, thread_pool)
        return
    _check_kernel_arrays(router, q_t=q_t, discharge_array=discharge_array, catchment_runoff=None)
    n_rivers, n_steps = router.network.river_ids.shape[0], discharge_array.shape[1]
    if runoff.runoff_by_cell.ndim != 2 or runoff.runoff_by_cell.shape[1] < n_steps:
        raise ValueError(f'runoff_by_cell has shape {runoff.runoff_by_cell.shape}, expected (n_cells, >= {n_steps})')
    if not runoff.runoff_by_cell.flags.c_contiguous:
        raise ValueError('runoff_by_cell must be C-contiguous, each cell series contiguous')
    if runoff.indptr.shape[0] != n_rivers + 1 or runoff.weight.shape != runoff.cell.shape:
        raise ValueError(f'the weight table does not describe the {n_rivers} rivers being routed')
    if runoff.cell.size and runoff.cell.max() >= runoff.runoff_by_cell.shape[0]:
        raise ValueError(f'the weight table references cells beyond the {runoff.runoff_by_cell.shape[0]} read')
    if runoff.scale.shape[0] not in (0, n_rivers):
        raise ValueError(f'scale has shape {runoff.scale.shape}, expected ({n_rivers},) or (0,)')
    _route(router, q_t, discharge_array, runoff, thread_pool)
    return
