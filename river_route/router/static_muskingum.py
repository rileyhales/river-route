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

from ..types import Float64Array, FloatArray, IntArray, RiverForcing
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
    """
    The Muskingum coefficients of each sub-reach, one per river that is not split, fixed for the whole run. Those of a
    river of one reach are of a whole runoff step, its subcycles combined into one step.
    """

    c1: FloatArray  # (sub-reaches,)
    c2: FloatArray  # (sub-reaches,)
    c3: FloatArray  # (sub-reaches,)
    c4: FloatArray  # (sub-reaches,) the coefficient of the river's runoff volume (m³) of a runoff step


def prepare_routing(
    network: Network, configs: Configs, dt_runoff: int, dt_routing: int | None
) -> tuple[ComputePlan, StaticCoefficients]:
    """
    Each river's compute plan and the coefficients of each of its sub-reaches. Every river hands its discharge
    downstream at the end of each runoff step and takes subcycles of its own within it. With dt_routing auto (None), a
    river too long for the runoff step is split into the sub-reaches of Network.sub_reach_k_and_x, each taking its
    share of the river's runoff, and a river too short for it takes the fewest subcycles that keep its coefficients
    non-negative. With dt_routing given, every river is routed whole in subcycles of exactly dt_routing, and a river
    with a negative coefficient at it is reported.

    The subcycles of a river of one reach are combined into one step by _combine_subcycles, which gives the discharge
    of taking them one at a time, so only a river split into sub-reaches takes its subcycles one by one.

    Raises:
        ValueError: if the coefficients do not sum to 1, or invalid_coefficients is raise and a coefficient is negative
    """
    if dt_routing:
        sub_reaches = np.ones(network.size, dtype=np.int64)
        subcycles = np.full(network.size, dt_runoff // dt_routing, dtype=np.int64)
    else:
        sub_reaches, subcycles, _ = network.conditioning(dt_runoff)
    network.check_validity(dt_runoff, action=configs.invalid_coefficients, sub_reaches=sub_reaches, subcycles=subcycles)
    k_reach, x_reach = network.sub_reach_k_and_x(sub_reaches)
    share = k_reach / np.repeat(network.k.astype(np.float64), sub_reaches)  # of the river's k, and of its runoff
    dt_div_k = np.repeat(dt_runoff / subcycles, sub_reaches) / k_reach  # each sub-reach's own step over its k
    denominator = dt_div_k + 2 * (1 - x_reach)
    c1 = (dt_div_k - 2 * x_reach) / denominator
    c2 = (dt_div_k + 2 * x_reach) / denominator
    c3 = (2 * (1 - x_reach) - dt_div_k) / denominator
    c4 = (c1 + c2) * share / dt_runoff
    one_reach = sub_reaches == 1
    c1, c2, c3, c4 = _combine_subcycles(c1, c2, c3, c4, np.repeat(np.where(one_reach, subcycles, 1), sub_reaches))
    if not np.allclose(c1 + c2 + c3, 1):
        raise ValueError('Muskingum coefficients do not sum to 1, check k and x of the network file and dt_routing')
    coefficients = StaticCoefficients(
        c1.astype(np.float32), c2.astype(np.float32), c3.astype(np.float32), c4.astype(np.float32)
    )
    return ComputePlan.from_counts(sub_reaches, np.where(one_reach, 1, subcycles)), coefficients


def _combine_subcycles(
    c1: Float64Array, c2: Float64Array, c3: Float64Array, c4: Float64Array, subcycles: IntArray
) -> tuple[Float64Array, Float64Array, Float64Array, Float64Array]:
    """
    The coefficients of one step that is s subcycles of a reach, its inflow rising linearly from I_t to I_{t+1} across
    them and its runoff spread evenly over them. Subcycle j of s ends at the fraction j/s of the step, with
    Q_j = c3 Q_(j-1) + c1 I(j/s) + c2 I((j-1)/s) + c4 V, so Q_s = c3^s Q_0 plus the terms of each subcycle j times
    c3^(s-j), whose inflows split between I_t and I_{t+1}: Q_(t+1) = c3' Q_t + c1' I_{t+1} + c2' I_t + c4' V with
    c3' = c3^s. Like each subcycle's, c1' + c2' + c3' is 1, and none is negative where none of the subcycle's is.
    """
    combined_c1 = np.zeros_like(c1)
    combined_c2 = np.zeros_like(c2)
    combined_c4 = np.zeros_like(c4)
    weight = np.ones_like(c3)  # c3^(s - j) of subcycle j, from the last subcycle back to the first
    for subcycles_after in range(int(subcycles.max())):
        subcycle = subcycles - subcycles_after  # j, from 1
        taken = subcycle >= 1
        end = subcycle / subcycles  # the fraction of the step at the end of subcycle j
        start = (subcycle - 1) / subcycles
        combined_c1 += np.where(taken, weight * (c1 * end + c2 * start), 0.0)
        combined_c2 += np.where(taken, weight * (c1 * (1 - end) + c2 * (1 - start)), 0.0)
        combined_c4 += np.where(taken, weight * c4, 0.0)
        weight = np.where(taken, weight * c3, weight)
    return combined_c1, combined_c2, weight, combined_c4


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
    downstream_inflow: FloatArray | None,
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
    if downstream_inflow is not None:  # numba compiles the sub-reaches that hand their outflow on without this
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
    downstream_inflow: FloatArray | None,
    outflow_at_steps: FloatArray,
) -> np.float32:
    """
    Q_{t+1} = c1 I_{t+1} + c2 I_t + c3 Q_t + c4 V_t at each step t from Q_0 = outflow, in one loop over eight steps at a
    time held in two vectors of four. The partial sums of the recurrence over each vector from Q = 0 do not depend on Q,
    so only Q_(t+i) = c3^i Q_t + partial_i waits on the steps before: one vector multiply-add per eight steps. Each
    Q_{t+1} is added into downstream_inflow[t + 1], unless it is None, and written into outflow_at_steps[t], and the
    last is returned.

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


@numba.njit(cache=True, nogil=True)
def _write_inflow_at_subcycles(inflow: FloatArray, subcycles: int, inflow_at_subcycles: FloatArray) -> None:
    """
    The inflow at the end of each subcycle, from the inflow before the first runoff step and at the end of each,
    interpolated linearly across the subcycles of each runoff step: a weighted mean of the inflow at the start and end
    of the step, which never leaves the range between them however it rounds.
    """
    # Compiled without fastmath contract: with it, whether these multiply-adds were fused differed between freshly
    # compiled code and the same code loaded from the numba cache, which changed the last bit of the discharge.
    one = np.float32(1.0)
    inverse_subcycles = one / np.float32(subcycles)
    inflow_at_subcycles[0] = inflow[0]
    for runoff_step in range(inflow.shape[0] - 1):
        for subcycle in range(1, subcycles + 1):
            after = np.float32(subcycle) * inverse_subcycles
            inflow_at_subcycles[runoff_step * subcycles + subcycle] = (one - after) * inflow[runoff_step] + after * (
                inflow[runoff_step + 1]
            )


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
) -> None:
    """
    Q_{t+1} = c1 I_{t+1} + c2 I_t + c3 Q_t + c4 V_t at each runoff step t, a reach in one step of each runoff step, its
    runoff read as it is routed and each Q written straight into discharge_at_runoff_steps.
    """
    outflow = channel_state[state_index]
    downstream_inflow[0] += outflow
    channel_state[state_index] = _route_river_eight_steps_at_a_time(
        c1, c2, c3, c4, river_forcing, outflow, upstream_inflow, downstream_inflow, discharge_at_runoff_steps
    )


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, inline='always')
def _route_sub_reaches_in_series(
    first_sub_reach: int,
    stop_sub_reach: int,
    coefficients: StaticCoefficients,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    inflow: FloatArray,
    scratch: FloatArray,
    steps: int,
) -> FloatArray:
    """
    Route the sub-reaches from first_sub_reach up to stop_sub_reach in series over steps steps, each taking the outflow
    of the one before as its inflow, and return the outflow row of the last: its Q before the first step and at the end
    of each. Two scratch rows take turns, since a sub-reach reads its inflow row while it writes its outflow row.
    """
    for sub_reach in range(first_sub_reach, stop_sub_reach):
        outflow = scratch[2 + (sub_reach - first_sub_reach) % 2, : steps + 1]
        outflow[0] = channel_state[sub_reach]
        channel_state[sub_reach] = _route_river_eight_steps_at_a_time(
            coefficients.c1[sub_reach],
            coefficients.c2[sub_reach],
            coefficients.c3[sub_reach],
            coefficients.c4[sub_reach],
            river_forcing,
            outflow[0],
            inflow,
            None,
            outflow[1:],
        )
        inflow = outflow
    return inflow


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, inline='always')
def _route_river_in_sub_reaches(
    first_sub_reach: int,
    sub_reaches: int,
    subcycles: int,
    coefficients: StaticCoefficients,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    discharge_at_runoff_steps: FloatArray,
    scratch: FloatArray,
) -> None:
    """
    Route a river too long for its step as sub-reaches in series, each in the subcycles of the river. With one a runoff
    step, each sub-reach reads the river's runoff as it is routed and the last is routed like a reach of one step. With
    more, the river's inflow is first interpolated across its subcycles and its runoff spread over them, and its Q at
    the end of each runoff step is taken from those of its subcycles.
    """
    runoff_steps = discharge_at_runoff_steps.shape[0]
    last_sub_reach = first_sub_reach + sub_reaches - 1
    if subcycles == 1:
        inflow = _route_sub_reaches_in_series(
            first_sub_reach,
            last_sub_reach,
            coefficients,
            river_forcing,
            channel_state,
            upstream_inflow,
            scratch,
            runoff_steps,
        )
        _route_river_whole(
            last_sub_reach,
            coefficients.c1[last_sub_reach],
            coefficients.c2[last_sub_reach],
            coefficients.c3[last_sub_reach],
            coefficients.c4[last_sub_reach],
            river_forcing,
            channel_state,
            inflow,
            downstream_inflow,
            discharge_at_runoff_steps,
        )
        return
    steps = runoff_steps * subcycles
    inflow_at_subcycles = scratch[0, : steps + 1]
    _write_inflow_at_subcycles(upstream_inflow, subcycles, inflow_at_subcycles)
    runoff_at_subcycles = scratch[1, :steps]  # each subcycle's runoff volume, that of its runoff step
    runoff_at_subcycles[:] = 0.0
    if river_forcing is not None:  # numba compiles channel routing without this
        add_forcing_terms(river_forcing, np.float32(1.0), runoff_steps, subcycles, runoff_at_subcycles)
    outflow = _route_sub_reaches_in_series(
        first_sub_reach,
        last_sub_reach + 1,
        coefficients,
        runoff_at_subcycles,
        channel_state,
        inflow_at_subcycles,
        scratch,
        steps,
    )
    for runoff_step in range(runoff_steps + 1):
        downstream_inflow[runoff_step] += outflow[runoff_step * subcycles]
    _write_discharge_at_runoff_steps(outflow[1:], subcycles, discharge_at_runoff_steps)


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
    """Route the river at river_index as one reach, or as sub-reaches in series where it is too long for its step."""
    first_sub_reach = compute_plan.sub_reach_starts[river_index]
    sub_reaches = compute_plan.sub_reach_starts[river_index + 1] - first_sub_reach
    if sub_reaches > 1:
        _route_river_in_sub_reaches(
            first_sub_reach,
            sub_reaches,
            compute_plan.subcycles[river_index],
            coefficients,
            river_forcing,
            channel_state,
            upstream_inflow,
            downstream_inflow,
            discharge_at_runoff_steps,
            scratch,
        )
    else:
        _route_river_whole(
            first_sub_reach,
            coefficients.c1[first_sub_reach],
            coefficients.c2[first_sub_reach],
            coefficients.c3[first_sub_reach],
            coefficients.c4[first_sub_reach],
            river_forcing,
            channel_state,
            upstream_inflow,
            downstream_inflow,
            discharge_at_runoff_steps,
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
    The static routine is itself the implementation numba compiles into route_job, and numba inlines its two
    routines into it. A call to a separately compiled function increments and decrements the reference count of every
    array it is passed, two atomic operations per array per river: a lambda calling the routine measured 1% slower, and
    calls to the three routines, 16 atomic operations per river, 2.6% slower on the Columbia.
    """
    if is_argument_type(coefficients, StaticCoefficients):
        return _route_river_with_static_coefficients
