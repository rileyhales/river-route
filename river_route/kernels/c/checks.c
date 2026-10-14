/*
 * The checks of a call of rr_route_rivers, made before any river is routed: every array is given, every count is
 * within the limits below, every index routing.c reads is inside its array, and every river drains into a river after
 * it in routing order. Every loop of the kernel is bounded by a count checked here. rr_describe_problem names the
 * problem a check finds.
 */
#include "checks.h"

enum {
    MAX_RIVERS = 16777216,      /* 2^24 */
    MAX_SUB_REACHES = 67108864, /* 2^26 */
    MAX_CELLS = 67108864,
    MAX_STEPS = 16777216,
    MAX_SUB_REACHES_PER_RIVER = 256,
    MAX_INFLOW_RIVERS_PER_RIVER = 65536,
    MAX_WEIGHTS_PER_RIVER = 65536,
};

const char *rr_describe_problem(int problem) {
    switch (problem) {
    case RR_OK: return "no problem";
    case RR_PROBLEM_MISSING_ARRAY: return "an array was not given";
    case RR_PROBLEM_SIZES: return "a count of rivers, sub-reaches, cells, steps, or blocks, or a stride, is wrong";
    case RR_PROBLEM_ORDER: return "a river drains into a river before it in routing order";
    case RR_PROBLEM_SUB_REACHES: return "a river has no sub-reaches, too many, or sub-reaches outside the network";
    case RR_PROBLEM_WEIGHTS: return "a river has weights out of order or a weight reads a cell outside the runoff";
    case RR_PROBLEM_BLOCKS: return "the blocks of rivers are not in routing order or leave the network";
    default: return "an unknown problem";
    }
}

/* The counts and arrays of a call are in range, and every series of runoff and discharge has room for its steps. */
static int check_job(struct rr_network network, struct rr_runoff runoff, struct rr_job job) {
    if (network.sub_reach_starts == 0 || network.inflow_river_starts == 0 || network.inflow_rivers == 0
        || network.c1 == 0 || network.c2 == 0 || network.c3 == 0 || network.c4 == 0 || network.weight_starts == 0
        || network.weight_cells == 0 || network.weights == 0 || job.block_starts == 0 || job.block_stops == 0
        || job.channel_state == 0 || job.discharge == 0 || (runoff.n_cells > 0 && runoff.depth == 0)) {
        return RR_PROBLEM_MISSING_ARRAY;
    }
    if (network.n_rivers < 1 || network.n_rivers > MAX_RIVERS || network.n_sub_reaches < network.n_rivers
        || network.n_sub_reaches > MAX_SUB_REACHES || runoff.n_cells < 0 || runoff.n_cells > MAX_CELLS
        || job.n_steps < 1 || job.n_steps > MAX_STEPS || job.n_blocks < 0 || job.n_blocks > MAX_RIVERS
        || job.river_stride <= job.n_steps || (runoff.n_cells > 0 && runoff.cell_stride < job.n_steps)
        || network.n_inflow_entries < 0 || network.n_weights < 0) {
        return RR_PROBLEM_SIZES;
    }
    if (network.sub_reach_starts[network.n_rivers] != network.n_sub_reaches) {
        return RR_PROBLEM_SIZES;
    }
    /*
     * The last river's row ends inside discharge, and the last cell's inside depth. They are compared by division,
     * since a count times a stride can overflow the 32 bit ptrdiff_t of WebAssembly, and the strides are positive.
     */
    if (job.discharge_size <= job.n_steps
        || (job.discharge_size - 1 - job.n_steps) / job.river_stride < network.n_rivers - 1) {
        return RR_PROBLEM_SIZES;
    }
    if (runoff.n_cells > 0
        && (runoff.depth_size < job.n_steps
            || (runoff.depth_size - job.n_steps) / runoff.cell_stride < runoff.n_cells - 1)) {
        return RR_PROBLEM_SIZES;
    }
    int32_t routed_to = 0; /* the blocks come in routing order, none overlapping the one before */
    for (int32_t block = 0; block < job.n_blocks; block++) {
        if (job.block_starts[block] < routed_to || job.block_stops[block] < job.block_starts[block]
            || job.block_stops[block] > network.n_rivers) {
            return RR_PROBLEM_BLOCKS;
        }
        routed_to = job.block_stops[block];
    }
    return RR_OK;
}

/* A river's sub-reaches, inflow rivers, and weights can be read, and it drains into no river before it. */
static int check_river(struct rr_network network, int32_t n_cells, int32_t river) {
    const int32_t first_sub_reach = network.sub_reach_starts[river];
    const int32_t n_sub_reaches = network.sub_reach_starts[river + 1] - first_sub_reach;
    if (first_sub_reach < 0 || n_sub_reaches < 1 || n_sub_reaches > MAX_SUB_REACHES_PER_RIVER
        || first_sub_reach > network.n_sub_reaches - n_sub_reaches) {
        return RR_PROBLEM_SUB_REACHES;
    }
    const int32_t first_entry = network.inflow_river_starts[river];
    const int32_t n_inflow_rivers = network.inflow_river_starts[river + 1] - first_entry;
    if (first_entry < 0 || n_inflow_rivers < 0 || n_inflow_rivers > MAX_INFLOW_RIVERS_PER_RIVER
        || first_entry > network.n_inflow_entries - n_inflow_rivers) {
        return RR_PROBLEM_ORDER;
    }
    for (int32_t entry = first_entry; entry < first_entry + n_inflow_rivers; entry++) {
        if (network.inflow_rivers[entry] < 0 || network.inflow_rivers[entry] >= river) {
            return RR_PROBLEM_ORDER;
        }
    }
    const int32_t first_weight = network.weight_starts[river];
    const int32_t n_weights = network.weight_starts[river + 1] - first_weight;
    if (first_weight < 0 || n_weights < 0 || n_weights > MAX_WEIGHTS_PER_RIVER
        || first_weight > network.n_weights - n_weights) {
        return RR_PROBLEM_WEIGHTS;
    }
    for (int32_t weight = first_weight; weight < first_weight + n_weights; weight++) {
        if (network.weight_cells[weight] < 0 || network.weight_cells[weight] >= n_cells) {
            return RR_PROBLEM_WEIGHTS;
        }
    }
    return RR_OK;
}

int start_job(const struct rr_network *network_pointer, const struct rr_runoff *runoff_pointer,
              const struct rr_job *job_pointer) {
    if (network_pointer == 0 || runoff_pointer == 0 || job_pointer == 0) {
        return RR_PROBLEM_MISSING_ARRAY;
    }
    const struct rr_network network = *network_pointer;
    const struct rr_job job = *job_pointer;
    int problem = check_job(network, *runoff_pointer, job);
    for (int32_t block = 0; block < job.n_blocks && problem == RR_OK; block++) {
        for (int32_t river = job.block_starts[block]; river < job.block_stops[block] && problem == RR_OK; river++) {
            problem = check_river(network, runoff_pointer->n_cells, river);
            if (problem == RR_OK) {
                job.discharge[(ptrdiff_t)river * job.river_stride] =
                    job.channel_state[network.sub_reach_starts[river + 1] - 1];
            }
        }
    }
    return problem;
}
