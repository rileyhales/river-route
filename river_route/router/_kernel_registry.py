from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING

import numpy as np

from ..types import FloatArray
from . import _numba_kernels as kernels

if TYPE_CHECKING:
    from .Router import Router

__all__ = ['KERNEL_REGISTRY', 'resolve_kernel', 'dispatch']

KERNEL_REGISTRY: dict[tuple[str, str, str, str], Callable[..., None]] = {}


def _describe_kernel(coeff: str, forcing: str, transform: str, network: str) -> str:
    return f'coeff={coeff}, forcing={forcing}, transform={transform}, network={network}'


def register(coeff: str, forcing: str, transform: str, network: str, boundary_buffers: int = 1):
    """
    Register the decorated function as the kernel-runner for one (coeff, forcing, transform, network) combination.

    ``boundary_buffers`` is how many (n_regions, n_routing_steps) arrays the kernel needs to hand each region's
    outlet contribution to the sequential remainder. Static coefficients need one (the finished contribution);
    dynamic coefficients need two, because the weights belong to the river being pushed into and are rebuilt
    every substep, so the two discharges are handed over instead of the product.
    """
    key = (coeff, forcing, transform, network)

    def deco(run: Callable[..., None]):
        if key in KERNEL_REGISTRY:
            raise ValueError(f'Duplicate kernel registration for {_describe_kernel(*key)}')
        run.boundary_buffers = boundary_buffers
        KERNEL_REGISTRY[key] = run
        return run

    return deco


def _lookup_transform(forcing: str, transform: str) -> str:
    """Channel-only routing applies no runoff transformation, so its kernels are keyed under 'uniform'
    regardless of the configured transform. Lateral forcing uses the configured transform verbatim."""
    return 'uniform' if forcing == 'channel' else transform


def resolve_kernel(coeff: str, forcing: str, transform: str, network: str) -> Callable[..., None]:
    """Look up the kernel-runner for a combination, or raise ``NotImplementedError`` listing what is implemented."""
    key = (coeff, forcing, _lookup_transform(forcing, transform), network)
    kernel = KERNEL_REGISTRY.get(key)
    if kernel is None:
        available = '\n  '.join(sorted(_describe_kernel(*k) for k in KERNEL_REGISTRY))
        raise NotImplementedError(
            f'No routing kernel for {_describe_kernel(coeff, forcing, transform, network)}.\n'
            f'Implemented combinations:\n  {available}'
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
    expected = {'q_t': ((n_rivers,), q_t.shape), 'discharge_array': ((n_steps, n_rivers), discharge_array.shape)}
    if vlateral is not None:
        expected['vlateral'] = ((n_steps, n_rivers), vlateral.shape)
    for name, (want, got) in expected.items():
        if got != want:
            raise ValueError(f'{name} has shape {got}, expected {want} for {n_rivers} rivers and {n_steps} steps')

    per_river = {
        'downstream_indices': router.network.downstream_indices,
        'c3': getattr(router, 'c3', None),
        'c4_dt': getattr(router, 'c4_dt', None),
        'downstream_c1': getattr(router, 'downstream_c1', None),
        'downstream_c2': getattr(router, 'downstream_c2', None),
        'alpha': router.network.alpha,
        'beta': router.network.beta,
        'x': router.network.x,
    }
    for name, array in per_river.items():
        if array is not None and array.shape != (n_rivers,):
            raise ValueError(f'{name} has shape {array.shape}, expected ({n_rivers},)')

    _check_partition_arrays(router, n_rivers)
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
        network=router.configs.network,
    )
    if router.configs.forcing == 'vlateral' and vlateral is None:
        raise ValueError('vlateral array must be provided for vlateral forcing')
    _check_kernel_arrays(router, q_t=q_t, discharge_array=discharge_array, vlateral=vlateral)

    n_regions = len(router.routing_jobs) - 1
    n_routing_steps = router.num_runoff_steps * router.num_routing_steps_per_runoff
    boundary = tuple(
        np.zeros((max(n_regions, 1), n_routing_steps), dtype=np.float32) for _ in range(kernel.boundary_buffers)
    )
    forcing = {'vlateral': vlateral} if router.configs.forcing == 'vlateral' else {}

    def run(job, cut_target: FloatArray) -> None:
        starts, stops, outlet, region = job
        kernel(
            router,
            q_t=q_t,
            discharge_array=discharge_array,
            block_starts=starts,
            block_stops=stops,
            outlet=outlet,
            region=region,
            cut_target=cut_target,
            boundary=boundary,
            **forcing,
        )

    if n_regions:
        regions = router.routing_jobs[:-1]  # already ordered longest first by the Router
        if thread_pool is None:
            for job in regions:
                run(job, _NO_CUTS)
        else:
            list(thread_pool.map(lambda job: run(job, _NO_CUTS), regions))  # list() so a worker exception propagates

    # the one barrier of the simulation: every region has finished before the main stem consumes its buffer
    run(router.routing_jobs[-1], router.cut_target)
    return


# ──────────────────────────────────────────────────────────────────────────────
# Kernel cells. Each runner reads persistent vectors off the router and takes the per-pass arrays as keywords.
# ──────────────────────────────────────────────────────────────────────────────


@register('static', 'channel', 'uniform', 'standard')
def static_channel(
    r: Router,
    *,
    q_t: FloatArray,
    discharge_array: FloatArray,
    block_starts: FloatArray,
    block_stops: FloatArray,
    outlet: int,
    region: int,
    cut_target: FloatArray,
    boundary: tuple[FloatArray, ...],
) -> None:
    kernels.static_channel(
        block_starts=block_starts,
        block_stops=block_stops,
        outlet=outlet,
        region=region,
        cut_target=cut_target,
        boundary=boundary[0],
        q_t=q_t,
        discharge_array=discharge_array,
        downstream_indices=r.network.downstream_indices,
        downstream_c1=r.downstream_c1,
        downstream_c2=r.downstream_c2,
        c3=r.c3,
        n_rivers=r.network.river_ids.shape[0],
        n_steps=r.num_runoff_steps,
        n_substeps=r.num_routing_steps_per_runoff,
    )


@register('static', 'vlateral', 'uniform', 'standard')
def static_vlateral(
    r: Router,
    *,
    q_t: FloatArray,
    discharge_array: FloatArray,
    vlateral: FloatArray,
    block_starts: FloatArray,
    block_stops: FloatArray,
    outlet: int,
    region: int,
    cut_target: FloatArray,
    boundary: tuple[FloatArray, ...],
) -> None:
    kernels.static_vlateral(
        block_starts=block_starts,
        block_stops=block_stops,
        outlet=outlet,
        region=region,
        cut_target=cut_target,
        boundary=boundary[0],
        q_t=q_t,
        discharge_array=discharge_array,
        downstream_indices=r.network.downstream_indices,
        downstream_c1=r.downstream_c1,
        downstream_c2=r.downstream_c2,
        c3=r.c3,
        n_rivers=r.network.river_ids.shape[0],
        n_steps=r.num_runoff_steps,
        n_substeps=r.num_routing_steps_per_runoff,
        vlateral=vlateral,
        c4_dt=r.c4_dt,
    )


@register('dynamic', 'vlateral', 'uniform', 'standard', boundary_buffers=2)
def dynamic_vlateral(
    r: Router,
    *,
    q_t: FloatArray,
    discharge_array: FloatArray,
    vlateral: FloatArray,
    block_starts: FloatArray,
    block_stops: FloatArray,
    outlet: int,
    region: int,
    cut_target: FloatArray,
    boundary: tuple[FloatArray, ...],
) -> None:
    kernels.dynamic_vlateral(
        block_starts=block_starts,
        block_stops=block_stops,
        outlet=outlet,
        region=region,
        cut_target=cut_target,
        boundary_old=boundary[0],
        boundary_new=boundary[1],
        q_t=q_t,
        discharge_array=discharge_array,
        downstream_indices=r.network.downstream_indices,
        alpha=r.network.alpha,
        beta=r.network.beta,
        x=r.network.x,
        dt_routing=np.float32(r.dt_routing),
        dt_runoff=np.float32(r.dt_runoff),
        n_rivers=r.network.river_ids.shape[0],
        n_steps=r.num_runoff_steps,
        n_substeps=r.num_routing_steps_per_runoff,
        vlateral=vlateral,
    )
