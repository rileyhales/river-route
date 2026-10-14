/*
 * The interface of river-route's routing kernel, routing.c with the checks of checks.c, which Python loads as a native
 * library and the browser as WebAssembly. numba_kernel.py routes the same arrays the same way, routine for routine.
 *
 * A network is rivers in routing order, every river after the rivers draining into it. Each river is one or more
 * sub-reaches, upstream to downstream, and every sub-reach takes one Muskingum step per runoff step:
 *
 *     Q[t+1] = c1 I[t+1] + c2 I[t] + c3 Q[t] + c4 V[t]
 *
 * where I is the sub-reach's inflow, the discharge of the sub-reach before it or, for a river's first sub-reach, the
 * summed discharge of the river's inflow rivers, and V[t] is the river's runoff volume in step t: the sum over its
 * weights of a weight times the runoff depth of the weight's cell. Runoff already aggregated to catchments is a grid
 * whose cells are the rivers, and channel routing is a network without weights. The coefficients are the same in
 * every step.
 *
 * Every river's discharge is a row: its value at the start of the steps routed, then at the end of each step. A river
 * reads its inflow from the rows of its inflow rivers, so a job writes nothing another job reads, and a job needs only
 * the jobs upstream of it to have finished.
 *
 * Every function returns RR_OK or a negative problem that rr_describe_problem names.
 */
#ifndef RIVER_ROUTE_ROUTING_H
#define RIVER_ROUTE_ROUTING_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum rr_problem {
    RR_OK = 0,
    RR_PROBLEM_MISSING_ARRAY = -1,
    RR_PROBLEM_SIZES = -2,
    RR_PROBLEM_ORDER = -3,
    RR_PROBLEM_SUB_REACHES = -4,
    RR_PROBLEM_WEIGHTS = -5,
    RR_PROBLEM_BLOCKS = -6,
};

/*
 * The rivers of a network in routing order and the sub-reaches they are routed as. The counts are the lengths of the
 * arrays, which the kernel reads nothing beyond.
 */
struct rr_network {
    int32_t n_rivers;
    int32_t n_sub_reaches;
    int32_t n_inflow_entries; /* the length of inflow_rivers */
    int32_t n_weights;        /* the length of weight_cells and of weights */
    const int32_t *sub_reach_starts;    /* n_rivers + 1: each river's first sub-reach, then n_sub_reaches */
    const int32_t *inflow_river_starts; /* n_rivers + 1: each river's first entry in inflow_rivers, then their count */
    const int32_t *inflow_rivers;       /* n_inflow_entries: the rivers draining into each river, each before it */
    const float *c1;                    /* n_sub_reaches: the Muskingum coefficients of each sub-reach */
    const float *c2;
    const float *c3;
    const float *c4;              /* n_sub_reaches: the coefficient of the river's runoff volume V */
    const int32_t *weight_starts; /* n_rivers + 1: each river's first weight, then the count of weights */
    const int32_t *weight_cells;  /* n_weights: the runoff cell each weight reads */
    const float *weights;         /* n_weights: the runoff volume (m3) a meter of depth in the cell gives the river */
};

/*
 * Runoff depth of every cell in each step routed, each cell's steps one after another: cell c in step t is
 * depth[c * cell_stride + t].
 */
struct rr_runoff {
    int32_t n_cells;
    ptrdiff_t cell_stride;
    ptrdiff_t depth_size; /* the length of depth */
    const float *depth;
};

/*
 * What one call routes: the rivers of a job, as blocks of rivers in routing order, over n_steps steps. A caller that
 * routes a long run a window of steps at a time calls once per window, with the runoff and discharge of that window.
 */
struct rr_job {
    int32_t n_steps;
    int32_t n_blocks;
    const int32_t *block_starts; /* n_blocks: the first river of each block */
    const int32_t *block_stops;  /* n_blocks: one past the last river of each block */
    float *channel_state;        /* n_sub_reaches: each sub-reach's discharge, at the start and then at the end */
    float *discharge; /* river r at time t is discharge[r * river_stride + t], t from 0 to n_steps */
    ptrdiff_t river_stride;
    ptrdiff_t discharge_size; /* the length of discharge */
};

/* Route each river of the job over every step before the next river, eight steps at a time. */
int rr_route_rivers(const struct rr_network *network, const struct rr_runoff *runoff, const struct rr_job *job);

const char *rr_describe_problem(int problem);

#ifdef __cplusplus
}
#endif

#endif
