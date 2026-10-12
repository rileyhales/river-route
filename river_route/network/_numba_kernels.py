import numba
import numpy as np

__all__ = ['depth_first_search_order', 'pack_blocks_into_jobs']


@numba.njit(cache=True, nogil=True)
def depth_first_search_order(
    downstream_indices,  # Array (n_rivers,) of the row of each river's downstream river, -1 at a basin outlet
):
    """
    The rows of the rivers in DFS order, each river after every river upstream of it, and the number of rivers upstream
    of each river in that order. Rivers draining into the same river, and basin outlets, keep their relative order.
    Rivers that drain into a cycle are downstream of no basin outlet and are left out, so fewer rows are returned.
    """
    n_rivers = downstream_indices.shape[0]
    root = n_rivers  # the river every basin outlet drains into
    drains_into = np.where(downstream_indices < 0, root, downstream_indices)
    inflow_rivers = np.argsort(drains_into, kind='mergesort')  # every river, grouped by the river it drains into
    inflow_river_starts = np.zeros(n_rivers + 2, dtype=np.int64)
    inflow_river_starts[1:] = np.cumsum(np.bincount(drains_into, minlength=n_rivers + 1))
    next_inflow_entry = inflow_river_starts[:-1].copy()
    first_row_upstream = np.zeros(n_rivers + 1, dtype=np.int64)  # the row of the first river of each watershed
    stack = np.empty(n_rivers + 1, dtype=np.int64)
    stack[0] = root
    stack_size = 1
    order = np.empty(n_rivers, dtype=np.int64)
    upstream_counts = np.empty(n_rivers, dtype=np.int64)
    n_ordered = 0
    # each river is pushed onto the stack once and popped once, and so is the root
    for _ in range(2 * n_rivers + 2):
        if stack_size == 0:
            break
        river = stack[stack_size - 1]
        if next_inflow_entry[river] < inflow_river_starts[river + 1]:
            inflow_river = inflow_rivers[next_inflow_entry[river]]
            next_inflow_entry[river] += 1
            first_row_upstream[inflow_river] = n_ordered
            stack[stack_size] = inflow_river
            stack_size += 1
            continue
        stack_size -= 1
        if river != root:
            order[n_ordered] = river
            upstream_counts[n_ordered] = n_ordered - first_row_upstream[river]
            n_ordered += 1
    return order[:n_ordered], upstream_counts[:n_ordered]


@numba.njit(cache=True, nogil=True)
def pack_blocks_into_jobs(
    block_sizes,  # Array (n_blocks,) of the number of rivers in each block
    n_jobs,  # number of jobs to pack the blocks into, one per thread
):
    """
    Pack blocks into jobs largest first: each block, from the largest down, goes to the job with the fewest rivers so
    far. Returns the job of each block and the number of rivers in each job.
    """
    job_of_block = np.empty(block_sizes.shape[0], dtype=np.int32)
    job_sizes = np.zeros(n_jobs, dtype=np.int64)
    for block in np.argsort(-block_sizes, kind='mergesort'):
        job = np.argmin(job_sizes)
        job_of_block[block] = job
        job_sizes[job] += block_sizes[block]
    return job_of_block, job_sizes
