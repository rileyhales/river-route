from __future__ import annotations

import heapq
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING

import numpy as np

from ..types import FloatArray
from . import _river_kernels as river_kernels

if TYPE_CHECKING:
    from ..runoff import CellRunoff
    from .Router import Router

__all__ = ['KERNEL_REGISTRY', 'resolve_kernel', 'dispatch', 'dispatch_grid']

KERNEL_REGISTRY: dict[tuple[str, str, str, str], Callable[..., None]] = {}


def _describe_kernel(coeff: str, forcing: str, transform: str, network_conditioning: str) -> str:
    return f'coeff={coeff}, forcing={forcing}, transform={transform}, network_conditioning={network_conditioning}'


def register(coeff: str, forcing: str, transform: str, network_conditioning: str):
    """Register the decorated function as the kernel-runner for one (coeff, forcing, transform,
    network_conditioning) combination."""
    key = (coeff, forcing, transform, network_conditioning)

    def deco(run: Callable[..., None]):
        if key in KERNEL_REGISTRY:
            raise ValueError(f'Duplicate kernel registration for {_describe_kernel(*key)}')
        KERNEL_REGISTRY[key] = run
        return run

    return deco


def _lookup_transform(forcing: str, transform: str) -> str:
    """Channel-only routing applies no runoff transformation, so its kernels are keyed under 'uniform'
    regardless of the configured transform. Lateral forcing uses the configured transform verbatim."""
    return 'uniform' if forcing == 'channel' else transform


def resolve_kernel(coeff: str, forcing: str, transform: str, network_conditioning: str) -> Callable[..., None]:
    """Look up the kernel-runner for a combination, or raise ``NotImplementedError`` listing what is implemented."""
    key = (coeff, forcing, _lookup_transform(forcing, transform), network_conditioning)
    kernel = KERNEL_REGISTRY.get(key)
    if kernel is None:
        available = '\n  '.join(sorted(_describe_kernel(*k) for k in KERNEL_REGISTRY))
        raise NotImplementedError(
            f'No routing kernel for {_describe_kernel(*key)}.\nImplemented combinations:\n  {available}'
        )
    return kernel


def _check_kernel_arrays(
    router: Router, q_t: FloatArray, discharge_array: FloatArray, vlateral: FloatArray | None
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
    expected = {'q_t': ((n_states,), q_t.shape), 'discharge_array': ((n_rivers, n_steps), discharge_array.shape)}
    if vlateral is not None:
        expected['vlateral'] = ((n_steps, n_rivers), vlateral.shape)
    for name, (want, got) in expected.items():
        if got != want:
            raise ValueError(f'{name} has shape {got}, expected {want} for {n_rivers} rivers and {n_steps} steps')
    _check_discharge_dtype(discharge_array)

    per_river = {
        'downstream_indices': router.network.downstream_indices,
        'c1': getattr(router, 'c1', None),
        'c2': getattr(router, 'c2', None),
        'c3': getattr(router, 'c3', None),
        'c4_dt': getattr(router, 'c4_dt', None),
        'alpha': router.network.alpha,
        'beta': router.network.beta,
        'x': router.network.x,
    }
    for name, array in per_river.items():
        if array is not None and array.shape != (n_rivers,):
            raise ValueError(f'{name} has shape {array.shape}, expected ({n_rivers},)')

    _check_partition_arrays(router, n_rivers)
    return


def _check_discharge_dtype(discharge_array: FloatArray) -> None:
    """
    Check that the discharge array is one the kernels can store into: float32, or the uint16 bit patterns a narrowed
    ``Configs.discharge_dtype`` fills, since numba has no float16 type. Router._discharge_buffer is what chooses
    between them; this only catches an array that came from somewhere else, for which numba would compile a float
    store as a numeric cast and silently write int(discharge).
    """
    if discharge_array.dtype not in (np.float32, np.uint16):
        raise ValueError(f'discharge_array has dtype {discharge_array.dtype}, expected float32 or uint16')
    return


_NO_CUTS: FloatArray = np.zeros(0, dtype=np.int32)  # a region pass injects nothing; only the main stem does


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


def dispatch(
    router: Router,
    q_t: FloatArray,
    discharge_array: FloatArray,
    vlateral: FloatArray | None = None,
    thread_pool: ThreadPoolExecutor | None = None,
) -> None:
    """
    Route the whole network, concurrently across regions on ``thread_pool`` when the schedule and config allow it.

    Threading is deliberately not part of the registry key. That key says which MATH to run; how the sweep is
    scheduled is orthogonal to it, and folding concurrency in would double every combination that has to be
    registered. The schedule itself is built once by the Router and only consumed here, so nothing about the
    partition is recomputed per input file.
    """
    # todo is there a way to do this without if branching? probably not without failing static checks
    # todo expand to match based on network also?
    kernel = resolve_kernel(
        coeff=router.configs.coeff,
        forcing=router.configs.forcing,
        transform=router.configs.transform,
        network_conditioning=router.configs.network_conditioning,
    )
    if router.configs.forcing == 'vlateral' and vlateral is None:
        raise ValueError('vlateral array must be provided for vlateral forcing')
    _check_kernel_arrays(router, q_t=q_t, discharge_array=discharge_array, vlateral=vlateral)

    forcing = {'vlateral': vlateral} if router.configs.forcing == 'vlateral' else {}
    _run_schedule(
        router, lambda **job: kernel(router, q_t=q_t, discharge_array=discharge_array, **forcing, **job), thread_pool
    )
    return


def _run_schedule(router: Router, run_pass: Callable[..., None], thread_pool: ThreadPoolExecutor | None) -> None:
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
    passes = _pack_regions(regions, router.threads if thread_pool is not None else 1)
    stem = dict(
        block_outlet=np.full(stem_starts.shape[0], -1, dtype=np.int32),
        block_region=np.zeros(stem_starts.shape[0], dtype=np.int32),
    )

    def run(job: dict) -> None:
        run_pass(cut_target=_NO_CUTS, boundary=boundary, **job)

    if passes:
        if thread_pool is None:
            for job in passes:
                run(job)
        else:
            list(thread_pool.map(run, passes))  # list() so a worker exception propagates

    # the one barrier of the simulation: every region has finished before the main stem consumes its buffer
    run_pass(block_starts=stem_starts, block_stops=stem_stops, cut_target=router.cut_target, boundary=boundary, **stem)
    return


def _pack_regions(regions: tuple, n_bins: int) -> list[dict]:
    """
    Pack region jobs, longest first, onto the least loaded of ``n_bins`` passes. Each pass lists its regions as
    blocks in index order with the outlet and boundary row of each.
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
        block_starts, block_stops, block_outlet, block_region = (
            np.array(v, dtype=np.int32) for v in zip(*jobs, strict=True)
        )
        passes.append(
            dict(
                block_starts=block_starts, block_stops=block_stops, block_outlet=block_outlet, block_region=block_region
            )
        )
    return passes


def dispatch_grid(
    router: Router,
    q_t: FloatArray,
    discharge_array: FloatArray,
    runoff: CellRunoff,
    thread_pool: ThreadPoolExecutor | None = None,
) -> None:
    """
    Aggregate gridded runoff and route it with the fused kernel, which never builds a vlateral array, concurrently
    across regions on ``thread_pool`` when given. Supports static coefficients with uniform lateral forcing on a
    standard or stabilized network; the Router chooses this path only for that combination.
    """
    _check_kernel_arrays(router, q_t=q_t, discharge_array=discharge_array, vlateral=None)
    n_rivers, n_steps = discharge_array.shape
    grid = router.runoff
    if runoff.runoff_by_cell.ndim != 2 or runoff.runoff_by_cell.shape[1] < n_steps:
        raise ValueError(f'runoff_by_cell has shape {runoff.runoff_by_cell.shape}, expected (n_cells, >= {n_steps})')
    if not runoff.runoff_by_cell.flags.c_contiguous:
        raise ValueError('runoff_by_cell must be C-contiguous, each cell series contiguous')
    if grid.indptr.shape[0] != n_rivers + 1 or runoff.weight.shape != grid.cell.shape:
        raise ValueError(f'the weight table does not describe the {n_rivers} rivers being routed')
    if grid.cell.size and grid.cell.max() >= runoff.runoff_by_cell.shape[0]:
        raise ValueError(f'the weight table references cells beyond the {runoff.runoff_by_cell.shape[0]} read')
    if runoff.scale.shape[0] not in (0, n_rivers):
        raise ValueError(f'scale has shape {runoff.scale.shape}, expected ({n_rivers},) or (0,)')

    def run_pass(**job):
        river_kernels.static_grid(
            q_t=q_t,
            discharge_array=discharge_array,
            downstream_indices=router.network.downstream_indices,
            c1=router.c1,
            c2=router.c2,
            c3=router.c3,
            c4_dt=router.c4_dt,
            n_substeps=router.num_routing_steps_per_runoff,
            runoff_by_cell=runoff.runoff_by_cell,
            indptr=grid.indptr,
            cell=grid.cell,
            weight=runoff.weight,
            scale=runoff.scale,
            cumulative=grid.cumulative,
            force_positive=grid.force_positive_runoff,
            block=river_kernels.BLOCK,
            reach_indptr=router.reach_indptr,
            substeps=router.substeps,
            **job,
        )

    _run_schedule(router, run_pass, thread_pool)
    return


# ──────────────────────────────────────────────────────────────────────────────
# Kernel cells. Each runner reads persistent vectors off the router and takes the per-pass arrays as keywords. A
# region's boundary row holds its outlet's whole series, so one buffer serves static and dynamic coefficients alike.
# ──────────────────────────────────────────────────────────────────────────────


def _river_layout(vlateral: FloatArray) -> tuple[FloatArray, bool]:
    """
    The array a river order kernel reads and whether it reads it by river. A (time, river) vlateral whose transpose
    is C-order, which is what a reader yields when it builds (river, time) arrays and hands over their ``.T``, is read
    one contiguous row per river; anything else is read as C-order (time, river).
    """
    if vlateral.T.flags.c_contiguous and not vlateral.flags.c_contiguous:
        return vlateral.T, True
    return np.ascontiguousarray(vlateral), False


@register('static', 'channel', 'uniform', 'stabilized')
@register('static', 'channel', 'uniform', 'standard')
def static_channel(r: Router, *, q_t: FloatArray, discharge_array: FloatArray, **job) -> None:
    river_kernels.static_channel(
        q_t=q_t,
        discharge_array=discharge_array,
        downstream_indices=r.network.downstream_indices,
        c1=r.c1,
        c2=r.c2,
        c3=r.c3,
        n_substeps=r.num_routing_steps_per_runoff,
        block=river_kernels.BLOCK,
        reach_indptr=r.reach_indptr,
        substeps=r.substeps,
        **job,
    )


@register('static', 'vlateral', 'uniform', 'stabilized')
@register('static', 'vlateral', 'uniform', 'standard')
def static_vlateral(r: Router, *, q_t: FloatArray, discharge_array: FloatArray, vlateral: FloatArray, **job) -> None:
    vlateral, by_river = _river_layout(vlateral)
    river_kernels.static_vlateral(
        q_t=q_t,
        discharge_array=discharge_array,
        downstream_indices=r.network.downstream_indices,
        c1=r.c1,
        c2=r.c2,
        c3=r.c3,
        c4_dt=r.c4_dt,
        n_substeps=r.num_routing_steps_per_runoff,
        vlateral=vlateral,
        by_river=by_river,
        block=river_kernels.BLOCK,
        reach_indptr=r.reach_indptr,
        substeps=r.substeps,
        **job,
    )


@register('dynamic', 'vlateral', 'uniform', 'standard')
def dynamic_vlateral(r: Router, *, q_t: FloatArray, discharge_array: FloatArray, vlateral: FloatArray, **job) -> None:
    vlateral, by_river = _river_layout(vlateral)
    river_kernels.dynamic_vlateral(
        q_t=q_t,
        discharge_array=discharge_array,
        downstream_indices=r.network.downstream_indices,
        alpha=r.network.alpha,
        beta=r.network.beta,
        x=r.network.x,
        dt_routing=r.dt_routing,
        dt_runoff=r.dt_runoff,
        n_substeps=r.num_routing_steps_per_runoff,
        vlateral=vlateral,
        by_river=by_river,
        block=river_kernels.BLOCK,
        **job,
    )
