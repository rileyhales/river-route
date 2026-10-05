"""
Muskingum routing with dynamic coefficients, the variable parameter Muskingum method: each river's k is recomputed from
its own discharge every routing step by the parameterization the configs choose, so its coefficients change as it
routes. Each parameterization reads its own columns of the network file, PARAMETERIZATION_COLUMNS, into its own
parameters, and is an overload of the muskingum_k stage:

- power_law: k = powerLawAlpha * Q ** powerLawBeta

Validity cannot be forced yet: routing a river whose k moves in substeps and subcycles needs a rule for choosing them.
"""

from typing import TYPE_CHECKING, NamedTuple

import numpy as np
from numba.extending import overload

from ._numba_kernels import STANDARD_LAYOUT, Layout, add_river_forcing, is_argument_type, muskingum_k, route_river

if TYPE_CHECKING:
    from ..configs.Configs import Configs
    from ..network.Network import Network

__all__ = ['DynamicMuskingum', 'PowerLaw', 'PARAMETERS_FOR_PARAMETERIZATION', 'prepare_routing']


class PowerLaw(NamedTuple):
    """The power law parameterization, k = alpha * Q ** beta, from the powerLawAlpha and powerLawBeta columns."""

    alpha: np.ndarray  # (n,)
    beta: np.ndarray  # (n,)


class DynamicMuskingum(NamedTuple):
    """Muskingum parameters per river whose k the parameterization recomputes from its discharge every step."""

    parameterization: PowerLaw
    x: np.ndarray  # (n,)
    dt_routing: np.float32
    inv_dt_runoff: np.float32


# the parameters of each parameterization, built from its PARAMETERIZATION_COLUMNS in the order of its fields
PARAMETERS_FOR_PARAMETERIZATION: dict[str, type[PowerLaw]] = {'power_law': PowerLaw}


def prepare_routing(
    network: Network, configs: Configs, dt_routing: int, dt_runoff: int
) -> tuple[Layout, DynamicMuskingum]:
    """
    The standard layout and the parameters of each river. The coefficients change as each river routes, so there is
    nothing for ``invalid_coefficients`` to check before routing starts.

    Args:
        network: the network whose x and parameterization columns the parameters are read from
        configs: ``parameterization`` chooses how k is recomputed; Configs refuses every parameterization but
            ``static`` with ``force_validity``
        dt_routing: the routing step in seconds
        dt_runoff: the runoff step in seconds, which turns each catchment runoff volume into an inflow rate

    Returns:
        tuple: (the standard Layout, the DynamicMuskingum parameters of every river)

    Raises:
        ValueError: if the network file lacks a column of the parameterization
    """
    columns = network.parameterization_columns(configs.parameterization)
    parameters = PARAMETERS_FOR_PARAMETERIZATION[configs.parameterization](*columns)
    return STANDARD_LAYOUT, DynamicMuskingum(parameters, network.x, np.float32(dt_routing), np.float32(1.0 / dt_runoff))


def _route_river_with_dynamic_coefficients(
    method, layout, q_t, r, catchment_runoff, inflow, downstream_inflow, discharge, work, chain
):
    """
    Variable parameter Muskingum for river r's whole series. The parameterization recomputes k from the river's own
    discharge every routing step, the coefficients are rebuilt from it, and the upstream series is weighted by them. The
    layout is always standard, so ``chain`` goes unused, and ``work`` holds the river's catchment runoff series, which
    each step reads as it is routed. The discharge of a runoff step is the trapezoidal mean of the river's discharges
    over it, as for static coefficients.
    """
    parameterization = method.parameterization
    x = method.x[r]
    dt_routing = method.dt_routing
    inv_dt_runoff = method.inv_dt_runoff
    zero = np.float32(0.0)
    half = np.float32(0.5)
    two_x = np.float32(2.0) * x
    two_one_minus_x = np.float32(2.0) * (np.float32(1.0) - x)
    n_steps = discharge.shape[0]
    # an inflow row holds the discharge before the first step and at every step
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
        before = q  # the discharge at the start of the runoff step
        interval_sum = zero
        for _ in range(routing_steps_per_runoff_step):
            dt_div_k = dt_routing / muskingum_k(parameterization, r, q)
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
        discharge[t] = (interval_sum + half * (before - q)) * inv_routing_steps_per_runoff_step
    q_t[r] = q
    return


@overload(route_river, jit_options={'nogil': True})
def _route_river_with_dynamic_muskingum(
    method, layout, q_t, r, catchment_runoff, inflow, downstream_inflow, discharge, work, chain
):
    if is_argument_type(method, DynamicMuskingum):
        return _route_river_with_dynamic_coefficients


def _muskingum_k_by_power_law(parameterization, r, q):
    """k = alpha * Q ** beta, with Q held at 1e-6 m³/s or more so that a negative beta never divides by zero."""
    return parameterization.alpha[r] * max(np.float32(1e-6), q) ** parameterization.beta[r]


@overload(muskingum_k, jit_options={'nogil': True})
def _muskingum_k_of_power_law(parameterization, r, q):
    if is_argument_type(parameterization, PowerLaw):
        return _muskingum_k_by_power_law
