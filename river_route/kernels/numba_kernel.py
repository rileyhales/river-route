"""
The routing kernel of routing.c in numba, routine for routine and under the same names, returning the same
problems. Runoff and discharge are flat arrays indexed with their strides, since a view of a row made inside a kernel
costs numba a reference count update.

route_river is inlined: a call between separately compiled numba functions updates the reference count of every array
it passes, atomically, and the threads routing jobs at once contend for those counts. route_sub_reach is compiled on
its own: inlined, LLVM no longer makes copies of its loops for sub-reaches with and without inflow rivers or weights.
"""

import numba
import numpy as np

from ._simd import (
    add,
    broadcast,
    four_lanes,
    last_lane,
    last_lane_in_every_lane,
    load_four,
    multiply,
    multiply_add,
    shift_in_one_lane,
    shift_in_two_lanes,
    store_four,
)

__all__ = ['PROBLEMS', 'RR_OK', 'route_rivers']

RR_OK = 0
RR_PROBLEM_MISSING_ARRAY = -1
RR_PROBLEM_SIZES = -2
RR_PROBLEM_ORDER = -3
RR_PROBLEM_SUB_REACHES = -4
RR_PROBLEM_WEIGHTS = -5
RR_PROBLEM_BLOCKS = -6
PROBLEMS = {  # checks.c's rr_describe_problem
    RR_OK: 'no problem',
    RR_PROBLEM_MISSING_ARRAY: 'an array was not given',
    RR_PROBLEM_SIZES: 'a count of rivers, sub-reaches, cells, steps, or blocks, or a stride, is wrong',
    RR_PROBLEM_ORDER: 'a river drains into a river before it in routing order',
    RR_PROBLEM_SUB_REACHES: 'a river has no sub-reaches, too many, or sub-reaches outside the network',
    RR_PROBLEM_WEIGHTS: 'a river has weights out of order or a weight reads a cell outside the runoff',
    RR_PROBLEM_BLOCKS: 'the blocks of rivers are not in routing order or leave the network',
}
MAX_RIVERS = 16777216  # checks.c's limits
MAX_SUB_REACHES = 67108864
MAX_CELLS = 67108864
MAX_STEPS = 16777216
MAX_SUB_REACHES_PER_RIVER = 256
MAX_INFLOW_RIVERS_PER_RIVER = 65536
MAX_WEIGHTS_PER_RIVER = 65536
ZERO = np.float32(0.0)


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, boundscheck=False)
def check_job(
    network,
    n_cells,
    depth_size,
    cell_stride,
    n_steps,
    block_starts,
    block_stops,
    channel_state,
    discharge_size,
    river_stride,
):
    """The counts and arrays of a call are in range, and every series of runoff and discharge has room for its steps."""
    sub_reach_starts, inflow_river_starts, inflow_rivers, c1, c2, c3, c4, weight_starts, weight_cells, weights = network
    n_rivers = sub_reach_starts.shape[0] - 1
    n_sub_reaches = channel_state.shape[0]
    rivers_given = inflow_river_starts.shape[0] == n_rivers + 1 and weight_starts.shape[0] == n_rivers + 1
    if (
        n_rivers < 1
        or n_rivers > MAX_RIVERS
        or n_sub_reaches < n_rivers
        or n_sub_reaches > MAX_SUB_REACHES
        or sub_reach_starts[n_rivers] != n_sub_reaches
        or not rivers_given
        or weight_cells.shape[0] != weights.shape[0]
        or n_cells < 0
        or n_cells > MAX_CELLS
        or n_steps < 1
        or n_steps > MAX_STEPS
        or block_starts.shape[0] != block_stops.shape[0]
        or block_starts.shape[0] > MAX_RIVERS
        or river_stride <= n_steps
        or (n_cells > 0 and cell_stride < n_steps)
    ):
        return RR_PROBLEM_SIZES
    for coefficients in (c1, c2, c3, c4):
        if coefficients.shape[0] != n_sub_reaches:
            return RR_PROBLEM_SIZES
    if (n_rivers - 1) * river_stride + n_steps >= discharge_size:
        return RR_PROBLEM_SIZES
    if n_cells > 0 and (n_cells - 1) * cell_stride + n_steps - 1 >= depth_size:
        return RR_PROBLEM_SIZES
    routed_to = 0  # the blocks come in routing order, none overlapping the one before
    for block in range(block_starts.shape[0]):
        if block_starts[block] < routed_to or block_stops[block] < block_starts[block] or block_stops[block] > n_rivers:
            return RR_PROBLEM_BLOCKS
        routed_to = block_stops[block]
    return RR_OK


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, boundscheck=False)
def check_river(network, n_cells, river):
    """A river's sub-reaches, inflow rivers, and weights can be read, and it drains into no river before it."""
    sub_reach_starts, inflow_river_starts, inflow_rivers, c1, c2, c3, c4, weight_starts, weight_cells, weights = network
    n_all_sub_reaches = sub_reach_starts[sub_reach_starts.shape[0] - 1]
    first_sub_reach = sub_reach_starts[river]
    n_sub_reaches = sub_reach_starts[river + 1] - first_sub_reach
    if (
        first_sub_reach < 0
        or n_sub_reaches < 1
        or n_sub_reaches > MAX_SUB_REACHES_PER_RIVER
        or first_sub_reach > n_all_sub_reaches - n_sub_reaches
    ):
        return RR_PROBLEM_SUB_REACHES
    first_entry = inflow_river_starts[river]
    n_inflow_rivers = inflow_river_starts[river + 1] - first_entry
    if (
        first_entry < 0
        or n_inflow_rivers < 0
        or n_inflow_rivers > MAX_INFLOW_RIVERS_PER_RIVER
        or first_entry > inflow_rivers.shape[0] - n_inflow_rivers
    ):
        return RR_PROBLEM_ORDER
    for entry in range(first_entry, first_entry + n_inflow_rivers):
        if inflow_rivers[entry] < 0 or inflow_rivers[entry] >= river:
            return RR_PROBLEM_ORDER
    first_weight = weight_starts[river]
    n_weights = weight_starts[river + 1] - first_weight
    if (
        first_weight < 0
        or n_weights < 0
        or n_weights > MAX_WEIGHTS_PER_RIVER
        or first_weight > weight_cells.shape[0] - n_weights
    ):
        return RR_PROBLEM_WEIGHTS
    for weight in range(first_weight, first_weight + n_weights):
        if weight_cells[weight] < 0 or weight_cells[weight] >= n_cells:
            return RR_PROBLEM_WEIGHTS
    return RR_OK


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, boundscheck=False)
def start_job(
    network,
    n_cells,
    depth_size,
    cell_stride,
    n_steps,
    block_starts,
    block_stops,
    channel_state,
    discharge,
    river_stride,
):
    """Check every river of a call and start it: each river's discharge row at the discharge of its last sub-reach."""
    problem = check_job(
        network,
        n_cells,
        depth_size,
        cell_stride,
        n_steps,
        block_starts,
        block_stops,
        channel_state,
        discharge.shape[0],
        river_stride,
    )
    if problem != RR_OK:
        return problem
    sub_reach_starts = network[0]
    for block in range(block_starts.shape[0]):
        for river in range(block_starts[block], block_stops[block]):
            problem = check_river(network, n_cells, river)
            if problem != RR_OK:
                return problem
            discharge[river * river_stride] = channel_state[sub_reach_starts[river + 1] - 1]
    return RR_OK


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, boundscheck=False, inline='always')
def inflow_at_step(discharge, river_stride, inflow_rivers, first_entry, stop_entry, step):
    """The summed discharge at time step of the rivers inflow_rivers[first_entry:stop_entry]."""
    inflow = ZERO
    for entry in range(first_entry, stop_entry):
        inflow += discharge[inflow_rivers[entry] * river_stride + step]
    return inflow


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, boundscheck=False, inline='always')
def runoff_volume_at_step(depth, cell_stride, weight_cells, weights, first_weight, stop_weight, step):
    """The river's runoff volume V in step."""
    volume = ZERO
    for weight in range(first_weight, stop_weight):
        volume += weights[weight] * depth[weight_cells[weight] * cell_stride + step]
    return volume


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, boundscheck=False, inline='always')
def powers_of_c3(c3):
    """c3 and c3^2 in every lane, then c3 to c3^4 and c3^5 to c3^8, a power a lane: routing.c's struct powers."""
    squared = c3 * c3
    fourth = squared * squared
    powers_1_to_4 = four_lanes(c3, squared, squared * c3, fourth)
    return broadcast(c3), broadcast(squared), powers_1_to_4, multiply(powers_1_to_4, broadcast(fourth))


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, boundscheck=False, inline='always')
def scan_four_steps(terms, powers):
    """
    The partial sums of Q[i+1] = c3 Q[i] + terms[i] over four steps from Q = 0: terms0, terms1 + c3 terms0, and so on:
    each lane plus c3 times the lane before it, then plus c3^2 times the sum two lanes before it.
    """
    c3, c3_squared = powers[0], powers[1]
    sums_of_two = multiply_add(c3, shift_in_one_lane(broadcast(ZERO), terms), terms)
    return multiply_add(c3_squared, shift_in_two_lanes(broadcast(ZERO), sums_of_two), sums_of_two)


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, boundscheck=False)
def route_sub_reach(
    c1,
    c2,
    c3,
    c4,
    inflow_rivers,
    first_entry,
    stop_entry,
    weight_cells,
    weights,
    first_weight,
    stop_weight,
    depth,
    cell_stride,
    discharge,
    river_stride,
    row_start,
    n_steps,
    discharge_at_start,
    inflow_at_start,
):
    """
    Route a sub-reach over every step, as routing.c's routine of the same name: the steps of whole eights eight at a
    time, then any left one at a time. Its inflow rivers are inflow_rivers[first_entry:stop_entry]. Returns its
    discharge at the end.
    """
    powers = powers_of_c3(c3)
    _, _, powers_1_to_4, powers_5_to_8 = powers
    c1_lanes, c2_lanes, c4_lanes = broadcast(c1), broadcast(c2), broadcast(c4)
    discharge_before = broadcast(discharge_at_start)  # Q at the end of the eight steps before
    inflow_before = broadcast(inflow_at_start)  # I at the end of the eight steps before, in its last lane
    step = 0
    while step + 8 <= n_steps:
        first_inflow_after = second_inflow_after = broadcast(ZERO)
        for entry in range(first_entry, stop_entry):
            inflow = inflow_rivers[entry] * river_stride + step
            first_inflow_after = add(first_inflow_after, load_four(discharge, inflow + 1))
            second_inflow_after = add(second_inflow_after, load_four(discharge, inflow + 5))
        first_volume = second_volume = broadcast(ZERO)
        for weight in range(first_weight, stop_weight):
            depth_at = weight_cells[weight] * cell_stride + step
            cell_weight = broadcast(weights[weight])
            first_volume = multiply_add(cell_weight, load_four(depth, depth_at), first_volume)
            second_volume = multiply_add(cell_weight, load_four(depth, depth_at + 4), second_volume)
        first_inflow_before = shift_in_one_lane(inflow_before, first_inflow_after)
        second_inflow_before = shift_in_one_lane(first_inflow_after, second_inflow_after)
        inflow_before = second_inflow_after
        first_terms = multiply_add(
            c1_lanes, first_inflow_after, multiply_add(c2_lanes, first_inflow_before, multiply(c4_lanes, first_volume))
        )
        second_terms = multiply_add(
            c1_lanes,
            second_inflow_after,
            multiply_add(c2_lanes, second_inflow_before, multiply(c4_lanes, second_volume)),
        )
        first_sums = scan_four_steps(first_terms, powers)
        second_sums = multiply_add(
            powers_1_to_4, last_lane_in_every_lane(first_sums), scan_four_steps(second_terms, powers)
        )
        first_discharge = multiply_add(powers_1_to_4, discharge_before, first_sums)
        second_discharge = multiply_add(powers_5_to_8, discharge_before, second_sums)
        store_four(discharge, row_start + step + 1, first_discharge)
        store_four(discharge, row_start + step + 5, second_discharge)
        discharge_before = last_lane_in_every_lane(second_discharge)
        step += 8
    discharge_now, inflow_now = last_lane(discharge_before), last_lane(inflow_before)
    for step in range(n_steps - n_steps % 8, n_steps):  # the steps left, one at a time
        inflow_after = inflow_at_step(discharge, river_stride, inflow_rivers, first_entry, stop_entry, step + 1)
        volume = runoff_volume_at_step(depth, cell_stride, weight_cells, weights, first_weight, stop_weight, step)
        discharge_now = c1 * inflow_after + c2 * inflow_now + c4 * volume + c3 * discharge_now
        discharge[row_start + step + 1] = discharge_now
        inflow_now = inflow_after
    return discharge_now


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, boundscheck=False, inline='always')
def route_river(network, river_itself, depth, cell_stride, discharge, river_stride, channel_state, river, n_steps):
    """Route a river over every step, its sub-reaches in turn, in place in its row: routing.c's."""
    sub_reach_starts, inflow_river_starts, inflow_rivers, c1, c2, c3, c4, weight_starts, weight_cells, weights = network
    row_start = river * river_stride
    reach_inflow_rivers = inflow_rivers
    first_entry, stop_entry = inflow_river_starts[river], inflow_river_starts[river + 1]
    first_weight, stop_weight = weight_starts[river], weight_starts[river + 1]
    inflow_at_start = inflow_at_step(discharge, river_stride, reach_inflow_rivers, first_entry, stop_entry, 0)
    for sub_reach in range(sub_reach_starts[river], sub_reach_starts[river + 1]):
        discharge_at_start = channel_state[sub_reach]
        channel_state[sub_reach] = route_sub_reach(
            c1[sub_reach],
            c2[sub_reach],
            c3[sub_reach],
            c4[sub_reach],
            reach_inflow_rivers,
            first_entry,
            stop_entry,
            weight_cells,
            weights,
            first_weight,
            stop_weight,
            depth,
            cell_stride,
            discharge,
            river_stride,
            row_start,
            n_steps,
            discharge_at_start,
            inflow_at_start,
        )
        inflow_at_start = discharge_at_start
        river_itself[0] = river  # every sub-reach after the first takes its inflow from the river's own row
        reach_inflow_rivers, first_entry, stop_entry = river_itself, 0, 1


@numba.njit(cache=True, nogil=True, fastmath={'contract'}, boundscheck=False)
def route_blocks(
    network, depth, cell_stride, n_steps, block_starts, block_stops, channel_state, discharge, river_stride
):
    """Route each river of the blocks over every step before the next river: routing.c's rr_route_rivers."""
    river_itself = np.zeros(1, dtype=np.int32)  # routing.c's &river: the inflow river of a river's later sub-reaches
    for block in range(block_starts.shape[0]):
        for river in range(block_starts[block], block_stops[block]):
            route_river(
                network, river_itself, depth, cell_stride, discharge, river_stride, channel_state, river, n_steps
            )


def flat_view(array: np.ndarray) -> tuple[np.ndarray, int]:
    """A 2D array whose rows hold their values one after another as a flat view of its memory, and its row stride in
    floats."""
    rows = array.strides[0] // array.itemsize
    extent = (array.shape[0] - 1) * rows + array.shape[1] if array.size else 0
    return np.lib.stride_tricks.as_strided(array, shape=(extent,), strides=(array.itemsize,)), rows


def route_rivers(network, depth, n_steps, block_starts, block_stops, channel_state, discharge) -> int:
    """routing.c's rr_route_rivers: runoff and discharge must have each series' steps one after another."""
    depth_flat, cell_stride = flat_view(depth)
    discharge_flat, river_stride = flat_view(discharge)
    problem = start_job(
        network,
        depth.shape[0],
        depth_flat.shape[0],
        cell_stride,
        n_steps,
        block_starts,
        block_stops,
        channel_state,
        discharge_flat,
        river_stride,
    )
    if problem == RR_OK:
        route_blocks(
            network,
            depth_flat,
            cell_stride,
            n_steps,
            block_starts,
            block_stops,
            channel_state,
            discharge_flat,
            river_stride,
        )
    return problem
