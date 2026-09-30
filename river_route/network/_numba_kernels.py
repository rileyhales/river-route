import numba
import numpy as np

__all__ = ['pack_blocks_into_jobs']


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
