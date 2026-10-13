/*
 * Four values of consecutive steps, one per lane of a vector register, and the operations routing.c needs on them: a
 * vector type of GCC and clang, which compiles to vector instructions natively and in WebAssembly alike, so no
 * processor's instructions are named here. The functions are those of river_route's _simd.py, under the same names.
 */
#ifndef RIVER_ROUTE_SIMD_H
#define RIVER_ROUTE_SIMD_H

typedef float four_lanes __attribute__((vector_size(16)));
typedef float unaligned_four_lanes __attribute__((vector_size(16), aligned(4), may_alias)); /* four floats of a row */

/* value in every lane. */
static inline four_lanes broadcast(float value) {
    const four_lanes lanes = {value, value, value, value};
    return lanes;
}

/* values[0] to values[3], from anywhere in a row: rows are read and written four steps from any step. */
static inline four_lanes load_four(const float *values) {
    return *(const unaligned_four_lanes *)values;
}

static inline void store_four(float *values, four_lanes lanes) {
    *(unaligned_four_lanes *)values = lanes;
}

static inline four_lanes add(four_lanes left, four_lanes right) {
    return left + right;
}

static inline four_lanes multiply(four_lanes left, four_lanes right) {
    return left * right;
}

/* left * right + addend in each lane, fused where the processor can. */
static inline four_lanes multiply_add(four_lanes left, four_lanes right, four_lanes addend) {
    return (left * right) + addend;
}

/*
 * Each lane of lanes moved up one lane, and the last lane of before in the first: (b3, l0, l1, l2). A shuffle picks
 * lanes by number, 0 to 3 from its first vector and 4 to 7 from its second.
 */
static inline four_lanes shift_in_one_lane(four_lanes before, four_lanes lanes) {
    return __builtin_shufflevector(before, lanes, 3, 4, 5, 6);
}

/* Each lane of lanes moved up two lanes, and the last two lanes of before in the first two: (b2, b3, l0, l1). */
static inline four_lanes shift_in_two_lanes(four_lanes before, four_lanes lanes) {
    return __builtin_shufflevector(before, lanes, 2, 3, 4, 5);
}

static inline four_lanes last_lane_in_every_lane(four_lanes lanes) {
    return __builtin_shufflevector(lanes, lanes, 3, 3, 3, 3);
}

static inline float last_lane(four_lanes lanes) {
    return lanes[3];
}

#endif
