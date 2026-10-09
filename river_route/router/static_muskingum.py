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
from ._numba_kernels import ComputePlan, add_forcing_of_four_steps, add_forcing_terms, is_argument_type, route_river
from ._simd import (
    FourLanes,
    add,
    broadcast,
    four_lanes,
    last_lane_in_every_lane,
    load_up_to_four,
    multiply,
    multiply_add,
    shift_in_one_lane,
    shift_in_two_lanes,
    store_up_to_four,
)

if TYPE_CHECKING:
    from ..configs.Configs import Configs
    from ..network.Network import Network

__all__ = ['StaticCoefficients', 'prepare_routing']


class StaticCoefficients(NamedTuple):
    """The Muskingum coefficients of each sub-reach, one per river that is not split, fixed for the whole run."""

    c1: FloatArray  # (sub-reaches,)
    c2: FloatArray  # (sub-reaches,)
    c3: FloatArray  # (sub-reaches,)
    c4: FloatArray  # (sub-reaches,) (c1 + c2) * share of k / dt_runoff, the coefficient of a runoff volume (m³)


def prepare_routing(
    network: Network, configs: Configs, dt_routing: int, dt_runoff: int
) -> tuple[ComputePlan, StaticCoefficients]:
    """
    Each river's compute plan and the coefficients of each sub-reach at dt_routing. With force_validity a river too
    long for dt_routing is split into the sub-reaches of Network.sub_reach_k_and_x, each taking its share of the
    river's runoff, and a river too short for it takes subcycles of dt_routing / subcycles.

    Raises:
        ValueError: if the coefficients do not sum to 1, or invalid_coefficients is raise and a coefficient is negative
    """
    sub_reaches = np.ones(network.size, dtype=np.int64)
    subcycles = np.ones(network.size, dtype=np.int64)
    if configs.force_validity:
        sub_reaches, subcycles, _ = network.conditioning(dt_routing)
    k_reach, x_reach = network.sub_reach_k_and_x(sub_reaches)
    share = k_reach / np.repeat(network.k.astype(np.float64), sub_reaches)  # of the river's k, and of its runoff
    k = k_reach.astype(np.float32)
    x = x_reach.astype(np.float32)
    dt = np.repeat(dt_routing / subcycles, sub_reaches).astype(np.float32)  # each sub-reach's own routing step
    dt_div_k = dt / k
    two_x = 2 * x
    denominator = dt_div_k + (2 * (1 - x))
    c1 = (dt_div_k - two_x) / denominator
    c2 = (dt_div_k + two_x) / denominator
    c3 = ((2 * (1 - x)) - dt_div_k) / denominator
    c4 = ((c1 + c2) * share / dt_runoff).astype(np.float32)
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


@numba.njit(cache=True, nogil=True)
def _write_discharge_at_runoff_steps(
    outflow: FloatArray, river_steps_per_runoff_step: int, discharge_at_runoff_steps: FloatArray
) -> None:
    """The river's Q at the end of each runoff step, the last of the m steps it takes in each."""
    for runoff_step in range(discharge_at_runoff_steps.shape[0]):
        discharge_at_runoff_steps[runoff_step] = outflow[(runoff_step + 1) * river_steps_per_runoff_step - 1]


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, inline='always')
def _scan_four_steps(terms: FourLanes, c3: FourLanes, c3_squared: FourLanes, zero: FourLanes) -> FourLanes:
    """
    The partial sums of Q_{t+1} = c3 Q_t + terms_t over four steps from Q = 0, terms_0, terms_1 + c3 terms_0, and so
    on: each lane plus c3 times the lane before it, then plus c3^2 times the sum two lanes before it.
    """
    sums_of_two = multiply_add(c3, shift_in_one_lane(zero, terms), terms)
    return multiply_add(c3_squared, shift_in_two_lanes(zero, sums_of_two), sums_of_two)


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, inline='always')
def _route_eight_steps(
    c1: np.float32,
    c2: np.float32,
    c3: np.float32,
    c4: np.float32,
    river_forcing: RiverForcing | None,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    outflow_at_steps: FloatArray,
    first_step: int,
    first_lanes: int,
    second_lanes: int,
    before: FourLanes,
) -> FourLanes:
    """
    The routing steps first_step to first_step + 7 of _route_river_eight_steps_at_a_time, from Q = before in every lane,
    in two vectors of first_lanes and second_lanes steps, four unless the series ends among them. Returns the Q at the
    end of the eight steps in every lane. The powers of c3 are the same at every call, so the compiler computes them
    once per river.
    """
    second_step = first_step + 4
    c3_squared = c3 * c3
    c3_fourth = c3_squared * c3_squared
    powers_1_to_4 = four_lanes(c3, c3_squared, c3_squared * c3, c3_fourth)
    powers_5_to_8 = multiply(powers_1_to_4, broadcast(c3_fourth))
    c1_lanes, c2_lanes, c3_lanes, c3_squared_lanes = broadcast(c1), broadcast(c2), broadcast(c3), broadcast(c3_squared)
    zero = broadcast(np.float32(0.0))
    first_terms = multiply_add(
        c1_lanes,
        load_up_to_four(upstream_inflow, first_step + 1, first_lanes),
        multiply(c2_lanes, load_up_to_four(upstream_inflow, first_step, first_lanes)),
    )
    second_terms = multiply_add(
        c1_lanes,
        load_up_to_four(upstream_inflow, second_step + 1, second_lanes),
        multiply(c2_lanes, load_up_to_four(upstream_inflow, second_step, second_lanes)),
    )
    if river_forcing is not None:  # numba compiles channel routing without this
        first_terms = add_forcing_of_four_steps(river_forcing, c4, first_step, first_lanes, first_terms)
        second_terms = add_forcing_of_four_steps(river_forcing, c4, second_step, second_lanes, second_terms)
    first_sums = _scan_four_steps(first_terms, c3_lanes, c3_squared_lanes, zero)
    second_sums = multiply_add(
        powers_1_to_4,
        last_lane_in_every_lane(first_sums),
        _scan_four_steps(second_terms, c3_lanes, c3_squared_lanes, zero),
    )
    first_outflow = multiply_add(powers_1_to_4, before, first_sums)
    second_outflow = multiply_add(powers_5_to_8, before, second_sums)
    first_downstream = add(load_up_to_four(downstream_inflow, first_step + 1, first_lanes), first_outflow)
    second_downstream = add(load_up_to_four(downstream_inflow, second_step + 1, second_lanes), second_outflow)
    store_up_to_four(downstream_inflow, first_step + 1, first_downstream, first_lanes)
    store_up_to_four(downstream_inflow, second_step + 1, second_downstream, second_lanes)
    store_up_to_four(outflow_at_steps, first_step, first_outflow, first_lanes)
    store_up_to_four(outflow_at_steps, second_step, second_outflow, second_lanes)
    return last_lane_in_every_lane(second_outflow)


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, inline='always')
def _route_river_eight_steps_at_a_time(
    c1: np.float32,
    c2: np.float32,
    c3: np.float32,
    c4: np.float32,
    river_forcing: RiverForcing | None,
    outflow: np.float32,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    outflow_at_steps: FloatArray,
) -> np.float32:
    """
    Q_{t+1} = c1 I_{t+1} + c2 I_t + c3 Q_t + c4 V_t at each routing step t from Q_0 = outflow, with a routing step per
    runoff step, in one loop over eight steps at a time held in two vectors of four. The partial sums of the recurrence
    over each vector from Q = 0 do not depend on Q, so only Q_(t+i) = c3^i Q_t + partial_i waits on the steps before:
    one vector multiply-add per eight steps. Each Q_{t+1} is added into downstream_inflow[t + 1] and written into
    outflow_at_steps[t], and the last is returned.

    jsrr's C kernel routes this way, written with WebAssembly SIMD. On the Amazon the same loop here routed 1.55 times
    faster than taking each river through the upstream terms, runoff, recurrence, downstream inflow, and discharge in
    loops of their own. Every eight steps but the last few take four lanes a vector, a count the compiler knows, so
    their loads and stores are single vector instructions: passing the count at every step cost the Amazon 0.55 s.
    """
    steps = outflow_at_steps.shape[0]
    steps_in_whole_blocks = steps - steps % 8
    before = broadcast(outflow)  # Q at the end of the eight steps before, in every lane
    for first_step in range(0, steps_in_whole_blocks, 8):
        before = _route_eight_steps(
            c1,
            c2,
            c3,
            c4,
            river_forcing,
            upstream_inflow,
            downstream_inflow,
            outflow_at_steps,
            first_step,
            4,
            4,
            before,
        )
    last_steps = steps - steps_in_whole_blocks
    if last_steps > 0:
        _route_eight_steps(
            c1,
            c2,
            c3,
            c4,
            river_forcing,
            upstream_inflow,
            downstream_inflow,
            outflow_at_steps,
            steps_in_whole_blocks,
            min(last_steps, 4),
            max(last_steps - 4, 0),
            before,
        )
    return outflow_at_steps[steps - 1] if steps > 0 else outflow


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
    discharge_at_runoff_steps: FloatArray,
    runoff_terms: FloatArray,
    outflow_at_steps: FloatArray,
) -> None:
    """
    Q_{t+1} = c1 I_{t+1} + c2 I_t + c3 Q_t + c4 V_t at each routing step t, the river routed as one reach. With a
    routing step per runoff step its runoff is read as it is routed and each Q is written straight into
    discharge_at_runoff_steps. With m routing steps per runoff step, c4 V_t is first spread over runoff_terms, which are
    then routed, and the Q at the end of each runoff step is written.
    """
    routing_steps = upstream_inflow.shape[0] - 1  # upstream_inflow[0] is I before the first routing step
    runoff_steps = discharge_at_runoff_steps.shape[0]
    routing_steps_per_runoff_step = routing_steps // runoff_steps
    outflow = channel_state[state_index]
    downstream_inflow[0] += outflow
    if routing_steps_per_runoff_step == 1:
        channel_state[state_index] = _route_river_eight_steps_at_a_time(
            c1, c2, c3, c4, river_forcing, outflow, upstream_inflow, downstream_inflow, discharge_at_runoff_steps
        )
        return
    runoff_terms[:] = 0.0
    if river_forcing is not None:  # numba compiles channel routing without this
        add_forcing_terms(river_forcing, c4, runoff_steps, routing_steps_per_runoff_step, runoff_terms)
    channel_state[state_index] = _route_river_eight_steps_at_a_time(
        c1, c2, c3, np.float32(1.0), runoff_terms, outflow, upstream_inflow, downstream_inflow, outflow_at_steps
    )
    _write_discharge_at_runoff_steps(outflow_at_steps, routing_steps_per_runoff_step, discharge_at_runoff_steps)


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, inline='always')
def _route_river_in_sub_reaches(
    first_sub_reach: int,
    sub_reaches: int,
    coefficients: StaticCoefficients,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    discharge_at_runoff_steps: FloatArray,
    known_terms: FloatArray,
    sub_reach_inflow: FloatArray,
) -> None:
    """Route a river too long for the routing step (c1 < 0) as sub-reaches in series, each routed whole."""
    routing_steps = upstream_inflow.shape[0] - 1
    runoff_steps = discharge_at_runoff_steps.shape[0]
    routing_steps_per_runoff_step = routing_steps // runoff_steps
    last_sub_reach = first_sub_reach + sub_reaches - 1
    outflow_before_first_step = channel_state[last_sub_reach]  # the river's outflow is its last sub-reach's
    inflow = upstream_inflow
    for sub_reach in range(first_sub_reach, last_sub_reach + 1):
        _write_upstream_terms(coefficients.c1[sub_reach], coefficients.c2[sub_reach], inflow, known_terms)
        if river_forcing is not None:  # c4 gives each sub-reach its share of the river's runoff
            add_forcing_terms(
                river_forcing, coefficients.c4[sub_reach], runoff_steps, routing_steps_per_runoff_step, known_terms
            )
        sub_reach_outflow_before_first_step = channel_state[sub_reach]
        channel_state[sub_reach] = _solve_outflow_recurrence(
            sub_reach_outflow_before_first_step, coefficients.c3[sub_reach], known_terms, routing_steps
        )
        if sub_reach < last_sub_reach:  # this sub-reach's outflow is the next one's inflow
            sub_reach_inflow[0] = sub_reach_outflow_before_first_step
            for routing_step in range(routing_steps):
                sub_reach_inflow[routing_step + 1] = known_terms[routing_step]
            inflow = sub_reach_inflow
    outflow = known_terms  # the last sub-reach's Q at each routing step
    _add_outflow_to_downstream_inflow(outflow_before_first_step, outflow, downstream_inflow)
    _write_discharge_at_runoff_steps(outflow, routing_steps_per_runoff_step, discharge_at_runoff_steps)


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
    discharge_at_runoff_steps: FloatArray,
    known_terms: FloatArray,
) -> None:
    """Route a river too short for the routing step (c3 < 0) in equal subcycles of its own per routing step."""
    routing_steps = upstream_inflow.shape[0] - 1
    runoff_steps = discharge_at_runoff_steps.shape[0]
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
    _write_discharge_at_runoff_steps(outflow, river_steps_per_runoff_step, discharge_at_runoff_steps)


def _route_river_with_static_coefficients(
    river_index: int,
    compute_plan: ComputePlan,
    coefficients: StaticCoefficients,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    discharge_at_runoff_steps: FloatArray,
    scratch: FloatArray,
) -> None:
    """Route the river at river_index whole, split into sub-reaches where c1 < 0, or in subcycles where c3 < 0."""
    known_terms = scratch[0]
    sub_reach_inflow = scratch[1, : upstream_inflow.shape[0]]  # a river split in sub-reaches has no subcycles
    routing_steps = upstream_inflow.shape[0] - 1
    first_sub_reach = compute_plan.sub_reach_starts[river_index]
    sub_reaches = compute_plan.sub_reach_starts[river_index + 1] - first_sub_reach
    subcycles = compute_plan.subcycles[river_index]
    c1 = coefficients.c1[first_sub_reach]  # a river routed whole or in subcycles is its one sub-reach
    c2 = coefficients.c2[first_sub_reach]
    c3 = coefficients.c3[first_sub_reach]
    c4 = coefficients.c4[first_sub_reach]
    if sub_reaches > 1:
        _route_river_in_sub_reaches(
            first_sub_reach,
            sub_reaches,
            coefficients,
            river_forcing,
            channel_state,
            upstream_inflow,
            downstream_inflow,
            discharge_at_runoff_steps,
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
            discharge_at_runoff_steps,
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
            discharge_at_runoff_steps,
            scratch[0, :routing_steps],
            scratch[1, :routing_steps],
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
    discharge_at_runoff_steps: FloatArray,
    scratch: FloatArray,
):
    """
    The static routine is itself the implementation numba compiles into route_job, and numba inlines its three
    routines into it. A call to a separately compiled function increments and decrements the reference count of every
    array it is passed, two atomic operations per array per river: a lambda calling the routine measured 1% slower, and
    calls to the three routines, 16 atomic operations per river, 2.6% slower on the Columbia.
    """
    if is_argument_type(coefficients, StaticCoefficients):
        return _route_river_with_static_coefficients
