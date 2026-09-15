import numba

__all__ = ['aggregate_to_rivers']

# The aggregation kernel covers a continuous range of rivers, which could be the whole network, so that code does not
# change between the single-threaded and multithreaded cases. Disjoint ranges can run concurrently because each
# writes only its own output columns, provided each is given its own scratch block.


# no fastmath: it lets the compiler assume values are never NaN, which would remove the NaN replacement
@numba.njit(cache=True, nogil=True)
def aggregate_to_rivers(
    *,
    runoff_by_cell,  # Array (n_cells, n_steps) of runoff depths, each cell's time series contiguous
    indptr,  # Array (n_rivers + 1,) of sparse row pointers, the weights of river r are indptr[r]:indptr[r + 1]
    cell,  # Array (n_weights,) of the row of runoff_by_cell each weight applies to
    weight,  # Array (n_weights,) of area proportions, already multiplied by the depth unit conversion factor
    scale,  # Array (n_rivers,) of per-river multipliers (catchment area for volumes); EMPTY to skip scaling
    zero,  # scalar zero in the dtype of scratch
    cumulative,  # bool, de-accumulate each river's series from cumulative to incremental
    force_positive,  # bool, clip negative values to zero
    replace_nan,  # bool, replace NaN values with zero
    r_start,  # first river index this call aggregates
    r_stop,  # one past the last river index this call aggregates
    scratch,  # Array (rivers_per_block, n_steps) of working space owned by this call
    out,  # Array (n_steps, n_rivers) C-order vlateral to write columns r_start:r_stop into
):
    """
    Area weighted sum of each river's cells, with every per-value conversion applied in the same pass:
    de-accumulation, clipping, NaN replacement, and the per-river scale.

    Each river's whole time series is summed from contiguous (cell, time) series, which the compiler vectorizes,
    into a scratch block of rivers small enough to stay in cache. Each finished block is then written into the
    C-order (time, river) output one short contiguous run per row.
    """
    n_steps = runoff_by_cell.shape[1]
    block = scratch.shape[0]
    for r0 in range(r_start, r_stop, block):
        r1 = min(r0 + block, r_stop)
        for r in range(r0, r1):
            row = scratch[r - r0]
            for t in range(n_steps):
                row[t] = zero
            for k in range(indptr[r], indptr[r + 1]):
                w = weight[k]
                series = runoff_by_cell[cell[k]]
                for t in range(n_steps):
                    row[t] += w * series[t]
            if cumulative:
                # descending, so each step subtracts the previous step's cumulative total before it is replaced
                for t in range(n_steps - 1, 0, -1):
                    row[t] -= row[t - 1]
            if force_positive:
                for t in range(n_steps):
                    if row[t] < zero:
                        row[t] = zero
            if replace_nan:
                for t in range(n_steps):
                    if row[t] != row[t]:
                        row[t] = zero
            if scale.shape[0]:
                s = scale[r]
                for t in range(n_steps):
                    row[t] *= s
        for t in range(n_steps):
            destination = out[t]
            for b in range(r1 - r0):
                destination[r0 + b] = scratch[b, t]
    return
