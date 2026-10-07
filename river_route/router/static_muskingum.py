"""
Muskingum routing with the following characteristics:
- Static coefficients per river. k, x, and dt are fixed for the whole run for a given river
- Forcing data may optionally be provided
- Handles reach-splitting and subcycling to force valid conditions
"""

from typing import TYPE_CHECKING, NamedTuple

import numba
import numpy as np
from numba.extending import overload

from ..types import FloatArray, RiverForcing
from ._numba_kernels import ComputePlan, add_forcing_terms, is_argument_type, route_river

if TYPE_CHECKING:
    from ..configs.Configs import Configs
    from ..network.Network import Network

__all__ = ['StaticCoefficients', 'prepare_routing']


class StaticCoefficients(NamedTuple):
    """Each river's Muskingum coefficients, fixed for the whole run and shared by all of its sub-reaches."""

    c1: FloatArray  # (n,)
    c2: FloatArray  # (n,)
    c3: FloatArray  # (n,)
    c4: FloatArray  # (n,) (c1 + c2) / (dt_runoff * sub_reaches), the coefficient of a runoff volume (m³)


def prepare_routing(
    network: Network, configs: Configs, dt_routing: int, dt_runoff: int
) -> tuple[ComputePlan, StaticCoefficients]:
    """
    Each river's compute plan and coefficients at dt_routing. With force_validity a river too long for dt_routing is
    split into sub-reaches of k / sub_reaches with the x of Network.sub_reach_x, and a river too short for it takes
    subcycles of dt_routing / subcycles.

    Raises:
        ValueError: if the coefficients do not sum to 1, or invalid_coefficients is raise and a coefficient is negative
    """
    sub_reaches = np.ones(network.size, dtype=np.int64)
    subcycles = np.ones(network.size, dtype=np.int64)
    if configs.force_validity:
        sub_reaches, subcycles, _ = network.conditioning(dt_routing)
    k = (network.k / sub_reaches).astype(np.float32)
    x = network.sub_reach_x(sub_reaches).astype(np.float32)
    dt = (dt_routing / subcycles).astype(np.float32)  # each river's own routing step
    dt_div_k = dt / k
    two_x = 2 * x
    denominator = dt_div_k + (2 * (1 - x))
    c1 = (dt_div_k - two_x) / denominator
    c2 = (dt_div_k + two_x) / denominator
    c3 = ((2 * (1 - x)) - dt_div_k) / denominator
    c4 = ((c1 + c2) / dt_runoff / sub_reaches).astype(np.float32)  # each sub-reach takes an equal share of runoff
    network.check_validity(
        dt_routing, action=configs.invalid_coefficients, sub_reaches=sub_reaches, subcycles=subcycles
    )
    if not np.allclose(c1 + c2 + c3, 1):
        raise ValueError('Muskingum coefficients do not sum to 1, check k and x of the network file and dt_routing')
    return ComputePlan.from_counts(sub_reaches, subcycles), StaticCoefficients(c1, c2, c3, c4)


@numba.njit(cache=True, nogil=True, fastmath={'contract'})
def _write_upstream_terms(c1: np.float32, c2: np.float32, inflow: FloatArray, known_terms: FloatArray) -> None:
    """known_terms[t] = c1 I_{t+1} + c2 I_t at each routing step t, where inflow[0] is I before the first step."""
    for routing_step in range(inflow.shape[0] - 1):
        known_terms[routing_step] = c1 * inflow[routing_step + 1] + c2 * inflow[routing_step]


@numba.njit(cache=True, nogil=True)
def _write_interpolated_upstream_terms(
    c1: np.float32, c2: np.float32, inflow: FloatArray, subcycles: int, known_terms: FloatArray
) -> None:
    """known_terms[s] = c1 I_{s+1} + c2 I_s at each subcycle s, with I interpolated linearly between routing steps."""
    # Compiled without fastmath contract: with it, whether these multiply-adds were fused differed between freshly
    # compiled code and the same code loaded from the numba cache, which changed the last bit of the discharge.
    inverse_subcycles = np.float32(1.0 / subcycles)
    for routing_step in range(inflow.shape[0] - 1):
        start = inflow[routing_step]
        rise = (inflow[routing_step + 1] - start) * inverse_subcycles
        before = start
        for subcycle in range(subcycles):
            after = start + rise * np.float32(subcycle + 1)
            known_terms[routing_step * subcycles + subcycle] = c1 * after + c2 * before
            before = after


@numba.njit(cache=True, nogil=True, fastmath={'contract'})
def _solve_outflow_recurrence(outflow: np.float32, c3: np.float32, known_terms: FloatArray, steps: int) -> np.float32:
    """Replace known_terms[t] with Q_{t+1} = c3 Q_t + known_terms[t] for the first steps, from Q_0 = outflow, and
    return the last Q."""
    # Eight steps at a time: within a group the partial sums of c3 powers times the known terms do not depend on Q, so
    # they overlap with the group before, and only Q_(t+8) = c3^8 Q_t + partial_7 waits on it. The serial chain is one
    # multiply-add per eight steps instead of one per step.
    c3_2 = c3 * c3  # c3_n is c3 to the power n
    c3_3 = c3_2 * c3
    c3_4 = c3_2 * c3_2
    c3_5 = c3_4 * c3
    c3_6 = c3_4 * c3_2
    c3_7 = c3_4 * c3_3
    c3_8 = c3_4 * c3_4
    step = 0
    while step + 8 <= steps:
        term_0 = known_terms[step]
        partial_1 = c3 * term_0 + known_terms[step + 1]
        partial_2 = c3 * partial_1 + known_terms[step + 2]
        partial_3 = c3 * partial_2 + known_terms[step + 3]
        partial_4 = c3 * partial_3 + known_terms[step + 4]
        partial_5 = c3 * partial_4 + known_terms[step + 5]
        partial_6 = c3 * partial_5 + known_terms[step + 6]
        partial_7 = c3 * partial_6 + known_terms[step + 7]
        known_terms[step] = c3 * outflow + term_0
        known_terms[step + 1] = c3_2 * outflow + partial_1
        known_terms[step + 2] = c3_3 * outflow + partial_2
        known_terms[step + 3] = c3_4 * outflow + partial_3
        known_terms[step + 4] = c3_5 * outflow + partial_4
        known_terms[step + 5] = c3_6 * outflow + partial_5
        known_terms[step + 6] = c3_7 * outflow + partial_6
        outflow = c3_8 * outflow + partial_7
        known_terms[step + 7] = outflow
        step += 8
    while step < steps:
        outflow = c3 * outflow + known_terms[step]
        known_terms[step] = outflow
        step += 1
    return outflow


@numba.njit(cache=True, nogil=True)
def _add_outflow_to_downstream_inflow(
    outflow_before_first_step: np.float32, outflow: FloatArray, downstream_inflow: FloatArray
) -> None:
    """Add the river's Q before the first routing step and at each routing step into the downstream river's inflow."""
    downstream_inflow[0] += outflow_before_first_step
    for routing_step in range(downstream_inflow.shape[0] - 1):
        downstream_inflow[routing_step + 1] += outflow[routing_step]


@numba.njit(cache=True, nogil=True, fastmath={'contract'})
def _write_mean_discharge(
    outflow_before_first_step: np.float32,
    outflow: FloatArray,
    river_steps_per_runoff_step: int,
    mean_discharge: FloatArray,
) -> None:
    """
    The mean Q over each runoff step of m of the river's steps, (Q_0/2 + Q_1 + ... + Q_(m-1) + Q_m/2) / m, since each Q
    averaged with the Q before it is the mean over the step between them.
    """
    half = np.float32(0.5)
    before = outflow_before_first_step
    if river_steps_per_runoff_step == 1:
        for runoff_step in range(mean_discharge.shape[0]):
            mean_discharge[runoff_step] = half * (before + outflow[runoff_step])
            before = outflow[runoff_step]
        return
    inverse_steps = np.float32(1.0 / river_steps_per_runoff_step)
    for runoff_step in range(mean_discharge.shape[0]):
        first_step = runoff_step * river_steps_per_runoff_step
        last_step = first_step + river_steps_per_runoff_step - 1
        step_sum = np.float32(0.0)
        for step in range(first_step, last_step + 1):
            step_sum += outflow[step]
        mean_discharge[runoff_step] = (step_sum + half * (before - outflow[last_step])) * inverse_steps
        before = outflow[last_step]


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, inline='always')
def _route_river_whole(
    state_index: int,
    c1: np.float32,
    c2: np.float32,
    c3: np.float32,
    c4: np.float32,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    mean_discharge: FloatArray,
    known_terms: FloatArray,
) -> None:
    """Q_{t+1} = c1 I_{t+1} + c2 I_t + c3 Q_t + c4 V_t at each routing step t, the river routed as one reach."""
    routing_steps = upstream_inflow.shape[0] - 1  # upstream_inflow[0] is I before the first routing step
    runoff_steps = mean_discharge.shape[0]
    routing_steps_per_runoff_step = routing_steps // runoff_steps
    _write_upstream_terms(c1, c2, upstream_inflow, known_terms)
    if river_forcing is not None:  # numba compiles channel routing without this
        add_forcing_terms(river_forcing, c4, runoff_steps, routing_steps_per_runoff_step, known_terms)
    outflow_before_first_step = channel_state[state_index]
    channel_state[state_index] = _solve_outflow_recurrence(outflow_before_first_step, c3, known_terms, routing_steps)
    outflow = known_terms  # the recurrence replaced each step's known terms with its Q
    _add_outflow_to_downstream_inflow(outflow_before_first_step, outflow, downstream_inflow)
    _write_mean_discharge(outflow_before_first_step, outflow, routing_steps_per_runoff_step, mean_discharge)


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, inline='always')
def _route_river_in_sub_reaches(
    first_sub_reach: int,
    sub_reaches: int,
    c1: np.float32,
    c2: np.float32,
    c3: np.float32,
    c4: np.float32,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    mean_discharge: FloatArray,
    known_terms: FloatArray,
    sub_reach_inflow: FloatArray,
) -> None:
    """Route a river too long for the routing step (c1 < 0) as equal sub-reaches in series, each routed whole."""
    routing_steps = upstream_inflow.shape[0] - 1
    runoff_steps = mean_discharge.shape[0]
    routing_steps_per_runoff_step = routing_steps // runoff_steps
    last_sub_reach = first_sub_reach + sub_reaches - 1
    outflow_before_first_step = channel_state[last_sub_reach]  # the river's outflow is its last sub-reach's
    inflow = upstream_inflow
    for sub_reach in range(first_sub_reach, last_sub_reach + 1):
        _write_upstream_terms(c1, c2, inflow, known_terms)
        if river_forcing is not None:  # c4 gives each sub-reach an equal share of the river's runoff
            add_forcing_terms(river_forcing, c4, runoff_steps, routing_steps_per_runoff_step, known_terms)
        sub_reach_outflow_before_first_step = channel_state[sub_reach]
        channel_state[sub_reach] = _solve_outflow_recurrence(
            sub_reach_outflow_before_first_step, c3, known_terms, routing_steps
        )
        if sub_reach < last_sub_reach:  # this sub-reach's outflow is the next one's inflow
            sub_reach_inflow[0] = sub_reach_outflow_before_first_step
            for routing_step in range(routing_steps):
                sub_reach_inflow[routing_step + 1] = known_terms[routing_step]
            inflow = sub_reach_inflow
    outflow = known_terms  # the last sub-reach's Q at each routing step
    _add_outflow_to_downstream_inflow(outflow_before_first_step, outflow, downstream_inflow)
    _write_mean_discharge(outflow_before_first_step, outflow, routing_steps_per_runoff_step, mean_discharge)


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, inline='always')
def _route_river_in_subcycles(
    state_index: int,
    subcycles: int,
    c1: np.float32,
    c2: np.float32,
    c3: np.float32,
    c4: np.float32,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    mean_discharge: FloatArray,
    known_terms: FloatArray,
) -> None:
    """Route a river too short for the routing step (c3 < 0) in equal subcycles of its own per routing step."""
    routing_steps = upstream_inflow.shape[0] - 1
    runoff_steps = mean_discharge.shape[0]
    river_steps = routing_steps * subcycles
    river_steps_per_runoff_step = (routing_steps // runoff_steps) * subcycles
    _write_interpolated_upstream_terms(c1, c2, upstream_inflow, subcycles, known_terms)
    if river_forcing is not None:
        add_forcing_terms(river_forcing, c4, runoff_steps, river_steps_per_runoff_step, known_terms)
    outflow_before_first_step = channel_state[state_index]
    channel_state[state_index] = _solve_outflow_recurrence(outflow_before_first_step, c3, known_terms, river_steps)
    outflow = known_terms  # the recurrence replaced each subcycle's known terms with its Q
    # the downstream river reads Q at the routing steps, the last subcycle of each
    downstream_inflow[0] += outflow_before_first_step
    for routing_step in range(routing_steps):
        downstream_inflow[routing_step + 1] += outflow[(routing_step + 1) * subcycles - 1]
    _write_mean_discharge(outflow_before_first_step, outflow, river_steps_per_runoff_step, mean_discharge)


def _route_river_with_static_coefficients(
    river_index: int,
    compute_plan: ComputePlan,
    coefficients: StaticCoefficients,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    mean_discharge: FloatArray,
    known_terms: FloatArray,
    sub_reach_inflow: FloatArray,
) -> None:
    """Route the river at river_index whole, split into sub-reaches where c1 < 0, or in subcycles where c3 < 0."""
    first_sub_reach = compute_plan.sub_reach_starts[river_index]
    sub_reaches = compute_plan.sub_reach_starts[river_index + 1] - first_sub_reach
    subcycles = compute_plan.subcycles[river_index]
    c1 = coefficients.c1[river_index]
    c2 = coefficients.c2[river_index]
    c3 = coefficients.c3[river_index]
    c4 = coefficients.c4[river_index]
    if sub_reaches > 1:
        _route_river_in_sub_reaches(
            first_sub_reach,
            sub_reaches,
            c1,
            c2,
            c3,
            c4,
            river_forcing,
            channel_state,
            upstream_inflow,
            downstream_inflow,
            mean_discharge,
            known_terms,
            sub_reach_inflow,
        )
    elif subcycles > 1:
        _route_river_in_subcycles(
            first_sub_reach,
            subcycles,
            c1,
            c2,
            c3,
            c4,
            river_forcing,
            channel_state,
            upstream_inflow,
            downstream_inflow,
            mean_discharge,
            known_terms,
        )
    else:
        _route_river_whole(
            first_sub_reach,
            c1,
            c2,
            c3,
            c4,
            river_forcing,
            channel_state,
            upstream_inflow,
            downstream_inflow,
            mean_discharge,
            known_terms,
        )


@overload(route_river, jit_options={'nogil': True, 'fastmath': {'contract'}})
def _overload_route_river_for_static_coefficients(
    river_index: int,
    compute_plan: ComputePlan,
    coefficients: StaticCoefficients,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    mean_discharge: FloatArray,
    known_terms: FloatArray,
    sub_reach_inflow: FloatArray,
):
    """
    The static routine is itself the implementation numba compiles into route_job, and numba inlines its three
    routines into it. A call to a separately compiled function increments and decrements the reference count of every
    array it is passed, two atomic operations per array per river: a lambda calling the routine measured 1% slower, and
    calls to the three routines, 16 atomic operations per river, 2.6% slower on the Columbia.
    """
    if is_argument_type(coefficients, StaticCoefficients):
        return _route_river_with_static_coefficients
