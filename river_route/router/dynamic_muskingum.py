"""
Muskingum routing with dynamic coefficients, the variable parameter Muskingum method: every routing step, a dynamic
method gives each river's k and x at its own discharge, and its coefficients are rebuilt from them. Each method reads
its columns of the network file, COLUMNS_FOR_PARAMETERS, into its own NamedTuple of parameters, and is an overload of
the k_and_x_at_discharge stage:

- power_law: k = powerLawAlpha * Q ** powerLawBeta, with the fixed x of muskingumX

Validity cannot be forced yet: routing a river whose k and x move in sub-reaches and subcycles needs a rule to choose
them.
"""

from typing import TYPE_CHECKING, NamedTuple

import numpy as np
from numba.extending import overload

from ..types import DynamicParameters, FloatArray, RiverForcing
from ._numba_kernels import ComputePlan, add_forcing_terms, is_argument_type, k_and_x_at_discharge, route_river

if TYPE_CHECKING:
    from ..configs.Configs import Configs
    from ..network.Network import Network

__all__ = ['DynamicCoefficients', 'PowerLaw', 'PARAMETERS_CLASS_FOR_PARAMETERS', 'prepare_routing']


class PowerLaw(NamedTuple):
    """The power law parameters, k = alpha * Q ** beta with a fixed x, from powerLawAlpha, powerLawBeta, muskingumX."""

    alpha: FloatArray  # (n,)
    beta: FloatArray  # (n,)
    x: FloatArray  # (n,)


class DynamicCoefficients(NamedTuple):
    """What each river's coefficients are rebuilt from at every routing step: the parameters of a dynamic method."""

    parameters: DynamicParameters
    dt_routing: np.float32
    inverse_dt_runoff: np.float32


# the NamedTuple of each dynamic method's parameters, built from its COLUMNS_FOR_PARAMETERS in the order of its fields
PARAMETERS_CLASS_FOR_PARAMETERS: dict[str, type[DynamicParameters]] = {'power_law': PowerLaw}


def prepare_routing(
    network: Network, configs: Configs, dt_routing: int, dt_runoff: int
) -> tuple[ComputePlan, DynamicCoefficients]:
    """
    A compute plan that routes every river whole, and what each river's coefficients are rebuilt from. The coefficients
    change as each river routes, so there is nothing for ``invalid_coefficients`` to check before routing starts.

    Raises:
        ValueError: if the network file lacks a column of the parameters
    """
    columns = network.columns_for_parameters(configs.parameters)
    parameters = PARAMETERS_CLASS_FOR_PARAMETERS[configs.parameters](*columns)
    whole_rivers = np.ones(network.size, dtype=np.int64)
    compute_plan = ComputePlan.from_counts(whole_rivers, whole_rivers)
    return compute_plan, DynamicCoefficients(parameters, np.float32(dt_routing), np.float32(1.0 / dt_runoff))


def _route_river_with_dynamic_coefficients(
    river_index: int,
    compute_plan: ComputePlan,
    coefficients: DynamicCoefficients,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    mean_discharge: FloatArray,
    known_terms: FloatArray,
    sub_reach_inflow: FloatArray,
) -> None:
    """
    Route the river at river_index whole with c1, c2, and c3 rebuilt from k_and_x_at_discharge at every routing step.
    They depend on Q, so each step is solved before the terms of the next are known.
    """
    parameters = coefficients.parameters
    dt_routing = coefficients.dt_routing
    zero = np.float32(0.0)
    half = np.float32(0.5)
    runoff_steps = mean_discharge.shape[0]
    routing_steps_per_runoff_step = (upstream_inflow.shape[0] - 1) // runoff_steps
    inverse_routing_steps_per_runoff_step = np.float32(1.0 / routing_steps_per_runoff_step)
    runoff_volumes = known_terms  # the scratch holds the river's runoff volume (m³) of each runoff step
    if river_forcing is not None:
        for runoff_step in range(runoff_steps):
            runoff_volumes[runoff_step] = zero
        add_forcing_terms(river_forcing, np.float32(1.0), runoff_steps, 1, runoff_volumes)
    outflow = channel_state[river_index]
    previous_inflow = upstream_inflow[0]
    downstream_inflow[0] += outflow
    routing_step = 1
    for runoff_step in range(runoff_steps):
        lateral_inflow = zero  # m³/s over the runoff step
        if river_forcing is not None:
            lateral_inflow = runoff_volumes[runoff_step] * coefficients.inverse_dt_runoff
        outflow_at_runoff_step_start = outflow
        step_sum = zero
        for _ in range(routing_steps_per_runoff_step):
            k, x = k_and_x_at_discharge(parameters, river_index, outflow)
            dt_div_k = dt_routing / k
            two_x = np.float32(2.0) * x
            two_one_minus_x = np.float32(2.0) * (np.float32(1.0) - x)
            denominator = dt_div_k + two_one_minus_x
            c1 = (dt_div_k - two_x) / denominator
            c2 = (dt_div_k + two_x) / denominator
            c3 = (two_one_minus_x - dt_div_k) / denominator
            inflow = upstream_inflow[routing_step]
            outflow = c3 * outflow + ((c1 + c2) * lateral_inflow + c1 * inflow + c2 * previous_inflow)
            previous_inflow = inflow
            downstream_inflow[routing_step] += outflow
            step_sum += outflow
            routing_step += 1
        mean_discharge[runoff_step] = (
            step_sum + half * (outflow_at_runoff_step_start - outflow)
        ) * inverse_routing_steps_per_runoff_step
    channel_state[river_index] = outflow
    return


@overload(route_river, jit_options={'nogil': True})
def _overload_route_river_for_dynamic_coefficients(
    river_index: int,
    compute_plan: ComputePlan,
    coefficients: DynamicCoefficients,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    mean_discharge: FloatArray,
    known_terms: FloatArray,
    sub_reach_inflow: FloatArray,
):
    if is_argument_type(coefficients, DynamicCoefficients):
        return _route_river_with_dynamic_coefficients


def _k_and_x_by_power_law(parameters: PowerLaw, river_index: int, q: np.float32) -> tuple[np.float32, np.float32]:
    """k = alpha * Q ** beta, with Q held at 1e-6 m³/s or more so that a negative beta never divides by zero, and x."""
    k = parameters.alpha[river_index] * max(np.float32(1e-6), q) ** parameters.beta[river_index]
    return k, parameters.x[river_index]


@overload(k_and_x_at_discharge, jit_options={'nogil': True})
def _overload_k_and_x_at_discharge_for_power_law(parameters: PowerLaw, river_index: int, q: np.float32):
    if is_argument_type(parameters, PowerLaw):
        return _k_and_x_by_power_law
