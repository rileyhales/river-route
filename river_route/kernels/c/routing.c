/*
 * The routing kernel of river-route in C17, the interface of routing.h, routed in the vectors of simd.h. It allocates
 * nothing and calls no library, so it compiles to a native library for Python and to WebAssembly for the browser.
 * Before routing a call, rr_route_rivers has start_job of checks.c check every index it reads, so every loop here is
 * bounded by a count checked there.
 *
 * route_river routes a river over every step of a call: each of its sub-reaches in turn over all of them, eight steps
 * at a time and any last steps one at a time, in the river's own discharge row, where each sub-reach reads the
 * discharge the sub-reach before it left and writes its own over it. rr_route_rivers routes each river of a job so
 * before the next river.
 */
#include "routing.h"
#include "checks.h"
#include "simd.h"

/* What routing one sub-reach reads: its coefficients, its river's runoff weights, and the rivers of its inflow. */
struct sub_reach {
    float c1;
    float c2;
    float c3;
    float c4;
    int32_t n_weights;
    const int32_t *weight_cells;
    const float *weights;
    int32_t n_inflow_rivers;
    const int32_t *inflow_rivers; /* the river's inflow rivers for its first sub-reach, and the river itself after */
};

/* A sub-reach's discharge and inflow at the start of the steps routed. */
struct discharge_and_inflow {
    float discharge;
    float inflow;
};

/* The powers of a sub-reach's c3 that solve eight steps at once. */
struct powers {
    four_lanes c3;            /* in every lane */
    four_lanes c3_squared;    /* in every lane */
    four_lanes powers_1_to_4; /* c3, c3^2, c3^3, c3^4 */
    four_lanes powers_5_to_8; /* c3^5 to c3^8 */
};

/* A river's discharge, its steps one after another. */
static inline float *discharge_row(struct rr_job job, int32_t river) {
    return job.discharge + ((ptrdiff_t)river * job.river_stride);
}

/* A cell's runoff depths, its steps one after another. */
static inline const float *cell_depths(struct rr_runoff runoff, int32_t cell) {
    return runoff.depth + ((ptrdiff_t)cell * runoff.cell_stride);
}

/* The summed discharge, at time step, of the rows a sub-reach takes its inflow from. */
static inline float inflow_at_step(struct rr_job job, struct sub_reach reach, int32_t step) {
    float inflow = 0.0F;
    for (int32_t entry = 0; entry < reach.n_inflow_rivers; entry++) {
        inflow += discharge_row(job, reach.inflow_rivers[entry])[step];
    }
    return inflow;
}

/* The river's runoff volume V in step. */
static inline float runoff_volume_at_step(struct rr_runoff runoff, struct sub_reach reach, int32_t step) {
    float volume = 0.0F;
    for (int32_t weight = 0; weight < reach.n_weights; weight++) {
        volume += reach.weights[weight] * cell_depths(runoff, reach.weight_cells[weight])[step];
    }
    return volume;
}

static inline struct powers powers_of_c3(float c3) {
    const float squared = c3 * c3;
    const float fourth = squared * squared;
    const four_lanes powers_1_to_4 = {c3, squared, squared * c3, fourth};
    const struct powers powers = {broadcast(c3), broadcast(squared), powers_1_to_4,
                                  multiply(powers_1_to_4, broadcast(fourth))};
    return powers;
}

/*
 * The partial sums of Q[i+1] = c3 Q[i] + terms[i] over four steps from Q = 0: terms0, terms1 + c3 terms0, and so on:
 * each lane plus c3 times the lane before it, then plus c3^2 times the sum two lanes before it.
 */
static inline four_lanes scan_four_steps(four_lanes terms, struct powers powers) {
    const four_lanes sums_of_two = multiply_add(powers.c3, shift_in_one_lane(broadcast(0.0F), terms), terms);
    return multiply_add(powers.c3_squared, shift_in_two_lanes(broadcast(0.0F), sums_of_two), sums_of_two);
}

/*
 * Route a sub-reach over every step of the job into row, from its discharge and inflow at the start, and return its
 * discharge at the end: the steps of whole eights eight at a time, then any steps left one at a time. Over eight steps
 * the recurrence Q[i+1] = c3 Q[i] + terms[i] is Q[i+1] = c3^(i+1) Q[0] + the partial sum of terms to i, which does not
 * depend on Q, so only one multiply-add per eight steps waits on the steps before. The steps run in two vectors of
 * four, the partial sums of the second taking the last of the first times the powers of c3. Each eight steps read the
 * inflow after each of their steps, and the inflow before each is the same moved up a lane behind the last inflow of
 * the eight steps before, so the sub-reach may write over its inflow row; the steps left go on from that last inflow
 * rather than read it again, since the row may hold the sub-reach's own discharge by then.
 */
static float route_sub_reach(struct rr_runoff runoff, struct rr_job job, struct sub_reach reach, float *row,
                             struct discharge_and_inflow start) {
    const struct powers powers = powers_of_c3(reach.c3);
    const four_lanes c1 = broadcast(reach.c1);
    const four_lanes c2 = broadcast(reach.c2);
    const four_lanes c4 = broadcast(reach.c4);
    four_lanes discharge_before = broadcast(start.discharge); /* Q at the end of the eight steps before */
    four_lanes inflow_before = broadcast(start.inflow);       /* I at the end of the eight steps before, last lane */
    for (int32_t step = 0; step + 8 <= job.n_steps; step += 8) {
        four_lanes first_inflow_after = broadcast(0.0F);
        four_lanes second_inflow_after = broadcast(0.0F);
        for (int32_t entry = 0; entry < reach.n_inflow_rivers; entry++) {
            const float *inflow = discharge_row(job, reach.inflow_rivers[entry]) + step;
            first_inflow_after = add(first_inflow_after, load_four(inflow + 1));
            second_inflow_after = add(second_inflow_after, load_four(inflow + 5));
        }
        four_lanes first_volume = broadcast(0.0F);
        four_lanes second_volume = broadcast(0.0F);
        for (int32_t weight = 0; weight < reach.n_weights; weight++) {
            const float *depth = cell_depths(runoff, reach.weight_cells[weight]) + step;
            const four_lanes cell_weight = broadcast(reach.weights[weight]);
            first_volume = multiply_add(cell_weight, load_four(depth), first_volume);
            second_volume = multiply_add(cell_weight, load_four(depth + 4), second_volume);
        }
        const four_lanes first_inflow_before = shift_in_one_lane(inflow_before, first_inflow_after);
        const four_lanes second_inflow_before = shift_in_one_lane(first_inflow_after, second_inflow_after);
        inflow_before = second_inflow_after;
        const four_lanes first_terms =
            multiply_add(c1, first_inflow_after, multiply_add(c2, first_inflow_before, multiply(c4, first_volume)));
        const four_lanes second_terms =
            multiply_add(c1, second_inflow_after, multiply_add(c2, second_inflow_before, multiply(c4, second_volume)));
        const four_lanes first_sums = scan_four_steps(first_terms, powers);
        const four_lanes second_sums = multiply_add(powers.powers_1_to_4, last_lane_in_every_lane(first_sums),
                                                    scan_four_steps(second_terms, powers));
        const four_lanes first_discharge = multiply_add(powers.powers_1_to_4, discharge_before, first_sums);
        const four_lanes second_discharge = multiply_add(powers.powers_5_to_8, discharge_before, second_sums);
        store_four(row + step + 1, first_discharge);
        store_four(row + step + 5, second_discharge);
        discharge_before = last_lane_in_every_lane(second_discharge);
    }
    float discharge_now = last_lane(discharge_before);
    float inflow_now = last_lane(inflow_before);
    for (int32_t step = job.n_steps - (job.n_steps % 8); step < job.n_steps; step++) { /* any left, one at a time */
        const float inflow_after = inflow_at_step(job, reach, step + 1);
        discharge_now = (reach.c1 * inflow_after) + (reach.c2 * inflow_now)
                        + (reach.c4 * runoff_volume_at_step(runoff, reach, step)) + (reach.c3 * discharge_now);
        row[step + 1] = discharge_now;
        inflow_now = inflow_after;
    }
    return discharge_now;
}

/*
 * Route a river over every step of the job: each sub-reach over them all in turn, in the river's row, its first
 * sub-reach from the rows of its inflow rivers and every other from the discharge the sub-reach before it left in the
 * row.
 */
static void route_river(struct rr_network network, struct rr_runoff runoff, struct rr_job job, int32_t river) {
    const int32_t first_entry = network.inflow_river_starts[river];
    const int32_t first_weight = network.weight_starts[river];
    struct sub_reach reach = {0.0F,
                              0.0F,
                              0.0F,
                              0.0F,
                              network.weight_starts[river + 1] - first_weight,
                              network.weight_cells + first_weight,
                              network.weights + first_weight,
                              network.inflow_river_starts[river + 1] - first_entry,
                              network.inflow_rivers + first_entry};
    float *row = discharge_row(job, river);
    float inflow_at_start = inflow_at_step(job, reach, 0);
    for (int32_t sub_reach = network.sub_reach_starts[river]; sub_reach < network.sub_reach_starts[river + 1];
         sub_reach++) {
        const struct discharge_and_inflow start = {job.channel_state[sub_reach], inflow_at_start};
        reach.c1 = network.c1[sub_reach];
        reach.c2 = network.c2[sub_reach];
        reach.c3 = network.c3[sub_reach];
        reach.c4 = network.c4[sub_reach];
        job.channel_state[sub_reach] = route_sub_reach(runoff, job, reach, row, start);
        inflow_at_start = start.discharge;
        reach.n_inflow_rivers = 1; /* every sub-reach after the first takes its inflow from the river's own row */
        reach.inflow_rivers = &river;
    }
}

int rr_route_rivers(const struct rr_network *network, const struct rr_runoff *runoff, const struct rr_job *job) {
    const int problem = start_job(network, runoff, job);
    if (problem != RR_OK) {
        return problem;
    }
    for (int32_t block = 0; block < job->n_blocks; block++) {
        for (int32_t river = job->block_starts[block]; river < job->block_stops[block]; river++) {
            route_river(*network, *runoff, *job, river);
        }
    }
    return RR_OK;
}
