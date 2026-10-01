"""
Muskingum routing with dynamic coefficients, the nonlinear Muskingum method: each river's travel time is rebuilt from
its own discharge every routing step as K = dynamicAlpha * Q ** dynamicBeta, so its coefficients change as it routes.
Standard networks only: stabilizing a river whose K moves needs a rule for choosing its substeps and subcycles.
"""

from typing import TYPE_CHECKING, Literal, NamedTuple

import numpy as np
from numba.extending import overload

from ._numba_kernels import STANDARD_LAYOUT, Layout, add_river_forcing, is_argument_type, route_river

if TYPE_CHECKING:
    from ..network.Network import Network

__all__ = ['DynamicMuskingum', 'prepare_routing']


class DynamicMuskingum(NamedTuple):
    """Nonlinear Muskingum parameters per river; the coefficients are rebuilt from K = alpha * Q ** beta each step."""

    alpha: np.ndarray  # (n,)
    beta: np.ndarray  # (n,)
    x: np.ndarray  # (n,)
    dt_routing: np.float32
    inv_dt_runoff: np.float32


def prepare_routing(
    network: Network,
    network_type: Literal['standard'],
    dt_routing: int,
    dt_runoff: int,
    unstable_coefficients: Literal['warn', 'raise', 'ignore'],
) -> tuple[Layout, DynamicMuskingum]:
    """
    The standard layout and each river's nonlinear parameters. The coefficients change as each river routes, so there
    is nothing for ``unstable_coefficients`` to check before routing starts.

    Args:
        network: the network whose dynamicAlpha, dynamicBeta, and x the parameters are read from
        network_type: always ``standard``; Configs refuses dynamic coefficients on a stabilized network
        dt_routing: the routing step in seconds
        dt_runoff: the runoff step in seconds, which turns each catchment runoff volume into an inflow rate
        unstable_coefficients: unused, since there are no fixed coefficients to check

    Returns:
        tuple: (the standard Layout, the DynamicMuskingum parameters of every river)

    Raises:
        ValueError: if the network file has no dynamicAlpha or dynamicBeta column
    """
    alpha, beta = network.dynamic_alpha, network.dynamic_beta
    if alpha is None or beta is None:
        raise ValueError('dynamic coefficients need dynamicAlpha and dynamicBeta columns in the network file')
    return STANDARD_LAYOUT, DynamicMuskingum(
        alpha, beta, network.x, np.float32(dt_routing), np.float32(1.0 / dt_runoff)
    )


def _route_river_with_dynamic_coefficients(
    method, layout, q_t, r, catchment_runoff, inflow, downstream_inflow, discharge, work, chain
):
    """
    Nonlinear Muskingum for river r's whole series. The coefficients are rebuilt from the river's own discharge every
    routing step, and the upstream series is weighted by them. The layout is always standard, so ``chain`` goes unused,
    and ``work`` holds the river's catchment runoff series, which each step reads as it is routed.
    """
    alpha = method.alpha[r]
    beta = method.beta[r]
    x = method.x[r]
    dt_routing = method.dt_routing
    inv_dt_runoff = method.inv_dt_runoff
    zero = np.float32(0.0)
    qmin = np.float32(1e-6)
    two_x = np.float32(2.0) * x
    two_one_minus_x = np.float32(2.0) * (np.float32(1.0) - x)
    n_steps = discharge.shape[0]
    # an inflow row holds the level before the first step and at every step
    routing_steps_per_runoff_step = (inflow.shape[0] - 1) // n_steps
    inv_routing_steps_per_runoff_step = np.float32(1.0 / routing_steps_per_runoff_step)
    if catchment_runoff is not None:
        for t in range(n_steps):
            work[t] = zero
        add_river_forcing(catchment_runoff, np.float32(1.0), n_steps, 1, work)
    q = q_t[r]
    u_prev = inflow[0]
    downstream_inflow[0] += q
    g = 1
    for t in range(n_steps):
        external = zero
        if catchment_runoff is not None:
            external = work[t] * inv_dt_runoff
        interval_sum = zero
        for _ in range(routing_steps_per_runoff_step):
            k = alpha * max(qmin, q) ** beta
            dt_div_k = dt_routing / k
            denominator = dt_div_k + two_one_minus_x
            c1 = (dt_div_k - two_x) / denominator
            c2 = (dt_div_k + two_x) / denominator
            c3 = (two_one_minus_x - dt_div_k) / denominator
            u = inflow[g]
            q = c3 * q + ((c1 + c2) * external + c1 * u + c2 * u_prev)
            u_prev = u
            downstream_inflow[g] += q
            interval_sum += q
            g += 1
        # coefficients that are not Muskingum-stable give negative discharge, which is clamped to zero here and does not
        # conserve mass; the unclamped value is what flows downstream and carries to the next step
        discharge[t] = max(zero, interval_sum * inv_routing_steps_per_runoff_step)
    q_t[r] = q
    return


@overload(route_river, jit_options={'nogil': True})
def _route_river_with_dynamic_muskingum(
    method, layout, q_t, r, catchment_runoff, inflow, downstream_inflow, discharge, work, chain
):
    if is_argument_type(method, DynamicMuskingum):
        return _route_river_with_dynamic_coefficients
