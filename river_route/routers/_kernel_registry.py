from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np

from ..types import FloatArray
from . import _numba_kernels as kernels

if TYPE_CHECKING:
    from .Router import Router

__all__ = ['KERNEL_REGISTRY', 'resolve_kernel', 'dispatch', ]

KERNEL_REGISTRY: dict[tuple[str, str, str], Callable[..., None]] = {}


def _describe_kernel(coeff: str, forcing: str, network: str) -> str:
    return f'coeff={coeff}, forcing={forcing}, network={network}'


def register(coeff: str, forcing: str, network: str):
    """Register the decorated function as the kernel-runner for one (coeff, forcing, network) combination."""
    key = (coeff, forcing, network)

    def deco(run: Callable[..., None]):
        if key in KERNEL_REGISTRY:
            raise ValueError(f'Duplicate kernel registration for {_describe_kernel(*key)}')
        KERNEL_REGISTRY[key] = run
        return run

    return deco


def resolve_kernel(coeff: str, forcing: str, network: str) -> Callable[..., None]:
    """Look up the kernel-runner for a combination, or raise ``NotImplementedError`` listing what is implemented."""
    kernel = KERNEL_REGISTRY.get((coeff, forcing, network))
    if kernel is None:
        available = '\n  '.join(sorted(_describe_kernel(*k) for k in KERNEL_REGISTRY))
        raise NotImplementedError(
            f'No routing kernel for {_describe_kernel(coeff, forcing, network)}.\n'
            f'Implemented combinations:\n  {available}'
        )
    return kernel


def dispatch(router: Router, q_t: FloatArray, discharge_array: FloatArray,
             vlateral: FloatArray | None = None) -> None:
    # todo is there a way to do this without if branching? probably not without failing static checks
    # todo expand to match based on network also?
    kernel = resolve_kernel(coeff=router.cfg.coeff, forcing=router.cfg.forcing, network=router.cfg.network)
    if router.cfg.forcing == 'lateral':
        if vlateral is None:
            raise ValueError('vlateral array must be provided for lateral forcing')
        kernel(router, q_t=q_t, discharge_array=discharge_array, vlateral=vlateral)
    else:
        kernel(router, q_t=q_t, discharge_array=discharge_array)


# ──────────────────────────────────────────────────────────────────────────────
# Kernel cells. Each runner reads persistent vectors off the router and takes the per-pass arrays as keywords.
# ──────────────────────────────────────────────────────────────────────────────

@register('static', 'channel', 'standard')
def static_channel(r: Router, *, q_t: FloatArray, discharge_array: FloatArray) -> None:
    kernels.static_channel(
        q_t=q_t, discharge_array=discharge_array,
        downstream_indices=r.downstream_indices,
        downstream_c1=r.downstream_c1, downstream_c2=r.downstream_c2, c3=r.c3,
        n_rivers=r.river_ids.shape[0], n_steps=r.num_runoff_steps, n_substeps=r.num_routing_steps_per_runoff,
    )


@register('static', 'lateral', 'standard')
def static_qlateral(r: Router, *, q_t: FloatArray, discharge_array: FloatArray, vlateral: FloatArray) -> None:
    kernels.static_qlateral(
        q_t=q_t, discharge_array=discharge_array,
        downstream_indices=r.downstream_indices,
        downstream_c1=r.downstream_c1, downstream_c2=r.downstream_c2, c3=r.c3,
        n_rivers=r.river_ids.shape[0], n_steps=r.num_runoff_steps, n_substeps=r.num_routing_steps_per_runoff,
        vlateral=vlateral, c4_dt=r.c4_dt,
    )


@register('dynamic', 'lateral', 'standard')
def dynamic_qlateral(r: Router, *, q_t: FloatArray, discharge_array: FloatArray, vlateral: FloatArray) -> None:
    kernels.dynamic_qlateral(
        q_t=q_t, discharge_array=discharge_array,
        downstream_indices=r.downstream_indices,
        alpha=r.alpha, beta=r.beta, x=r.x,
        dt_routing=np.float32(r.dt_routing), dt_runoff=np.float32(r.dt_runoff),
        n_rivers=r.river_ids.shape[0], n_steps=r.num_runoff_steps, n_substeps=r.num_routing_steps_per_runoff,
        vlateral=vlateral,
    )
