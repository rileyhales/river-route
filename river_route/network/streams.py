"""
Functions of a network table in DFS order: the divisors of a time step, which Network.largest_stable_dt searches, the
sub-watershed blocks threads route and how evenly they divide a network, and the subset of a network file, and of its
grid weights, to one river and every river upstream of it. The blocks and the subset only read the contiguous run of
rows that a river's riverIndex and upstreamCount give.
"""

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import xarray as xr

from ..types import Int32Array, IntArray, JobBlocks, PathInput
from ._numba_kernels import pack_blocks_into_jobs

if TYPE_CHECKING:
    from .Network import Network

__all__ = ['divisors_of', 'assign_blocks', 'analyze_partitioning', 'subset_network_to_river']

# the caps on the rivers in a block assign_blocks tries, as fractions of an even share of the rivers per job
_BLOCK_CAP_FRACTIONS: tuple[float, ...] = (0.125, 0.25, 0.375, 0.5, 0.65, 0.8, 1.0)


def divisors_of(n: int) -> np.ndarray:
    """All positive integer divisors of n, sorted ascending."""
    small = [d for d in range(1, int(n**0.5) + 1) if n % d == 0]
    return np.array(sorted(set(small + [n // d for d in small])), dtype=np.int64)


def assign_blocks(
    upstream_counts: IntArray, downstream_indices: Int32Array, threads: int
) -> tuple[tuple[JobBlocks, ...], Int32Array]:
    """
    Divide a region in DFS order into sub-watershed blocks, pack the blocks into ``threads`` jobs that route
    concurrently, and leave every other river to the main stem, the last job, which routes after them.

    In DFS order the watershed of river r is the rows ``r - upstream_counts[r]`` through r, so a block is known by its
    outlet alone. Under a cap on the rivers in a block, the block outlets are the rivers whose watershed fits the cap
    and whose downstream river's watershed does not, or that are a basin outlet. A river's watershed always holds more
    rivers than the watershed of any river upstream of it, so these blocks never overlap, and each drains into a river
    over the cap, which is in no block: the main stem.

    A loose cap leaves a short main stem but packs unevenly, a tight one packs evenly but leaves more rivers to the main
    stem, and the trade turns over at a different cap per network. So each cap of _BLOCK_CAP_FRACTIONS is tried, and the
    one whose busiest job and main stem together route the fewest rivers is kept. Blocks are kept only when they beat
    routing every river in one job, which they never do on one thread.

    Args:
        upstream_counts: (n,) number of rivers upstream of each river, the upstreamCount column
        downstream_indices: (n,) index of each river's downstream river, -1 at a basin outlet
        threads: number of jobs to pack the blocks into

    Returns:
        tuple: (job_blocks, cut_target). Each entry of job_blocks is (block_starts, block_stops, block_outlet,
            block_number) for one job, with its blocks in index order. The last entry is the main stem, whose blocks
            are the runs of rows between the sub-watershed blocks and have no outlet (-1). The outlet of block b drains
            into cut_target[b], -1 at a basin outlet.
    """
    n_rivers = upstream_counts.shape[0]
    watershed_sizes = upstream_counts.astype(np.int64) + 1
    downstream_watershed_sizes = np.where(downstream_indices >= 0, watershed_sizes[downstream_indices], n_rivers + 1)
    block_outlets = np.zeros(0, dtype=np.int32)
    job_of_block = np.zeros(0, dtype=np.int32)
    fewest_rivers_in_busiest_job_and_main_stem = n_rivers
    for fraction in _BLOCK_CAP_FRACTIONS:
        cap = max(1, int(fraction * n_rivers / threads))
        outlets = np.flatnonzero((watershed_sizes <= cap) & (downstream_watershed_sizes > cap))
        jobs, job_sizes = pack_blocks_into_jobs(watershed_sizes[outlets], threads)
        rivers_in_busiest_job_and_main_stem = int(job_sizes.max()) + n_rivers - int(watershed_sizes[outlets].sum())
        if rivers_in_busiest_job_and_main_stem < fewest_rivers_in_busiest_job_and_main_stem:
            block_outlets, job_of_block = outlets.astype(np.int32), jobs
            fewest_rivers_in_busiest_job_and_main_stem = rivers_in_busiest_job_and_main_stem

    block_starts = (block_outlets - upstream_counts[block_outlets]).astype(np.int32)
    block_stops = (block_outlets + 1).astype(np.int32)
    job_blocks: list[JobBlocks] = []
    for job in range(threads):
        block_numbers = np.flatnonzero(job_of_block == job).astype(np.int32)
        if block_numbers.shape[0]:
            job_blocks.append(
                (block_starts[block_numbers], block_stops[block_numbers], block_outlets[block_numbers], block_numbers)
            )
    stem_starts = np.concatenate(([0], block_stops)).astype(np.int32)
    stem_stops = np.concatenate((block_starts, [n_rivers])).astype(np.int32)
    in_main_stem = stem_stops > stem_starts
    no_outlets = np.full(np.count_nonzero(in_main_stem), -1, dtype=np.int32)
    job_blocks.append((stem_starts[in_main_stem], stem_stops[in_main_stem], no_outlets, np.zeros_like(no_outlets)))
    return tuple(job_blocks), downstream_indices[block_outlets].astype(np.int32)


def analyze_partitioning(network: Network, threads: int = 4) -> dict:
    """
    Static analysis of dividing a network into blocks for multithreaded routing (see assign_blocks).

    Reports the block count, how evenly they pack into ``threads`` jobs, how much of the network falls into
    the sequential main stem, and the speedup that implies. The speedup is an upper bound from work alone: it
    assumes perfect scaling and ignores memory bandwidth, which this kernel is largely bound by. Also prints a
    human-readable report.

    Args:
        network: the network to divide
        threads: thread count the blocks are sized for

    Returns:
        dict: the summary statistics
    """
    *sub_watershed_jobs, main_stem_job = network.routing_blocks(threads)[0]
    block_sizes = [int(size) for starts, stops, _, _ in sub_watershed_jobs for size in stops - starts]
    span = max((int(np.sum(stops - starts)) for starts, stops, _, _ in sub_watershed_jobs), default=0)
    main_stem = int(np.sum(main_stem_job[1] - main_stem_job[0]))
    n_rivers = network.size

    summary = {
        'threads': threads,
        'n_rivers': n_rivers,
        'n_blocks': len(block_sizes),
        'largest_block': max(block_sizes, default=0),
        'smallest_block': min(block_sizes, default=0),
        'parallel_span': span,
        'main_stem': main_stem,
        'main_stem_fraction': main_stem / n_rivers if n_rivers else 0.0,
        'predicted_speedup': n_rivers / (span + main_stem) if n_rivers else 1.0,
    }

    print(f'Partitioning analysis for threads={threads}')
    print(f'  rivers in:              {summary["n_rivers"]:,}')
    print(
        f'  sub-watershed blocks:   {summary["n_blocks"]:,} '
        f'(largest {summary["largest_block"]:,}, smallest {summary["smallest_block"]:,})'
    )
    print(f'  parallel span:          {summary["parallel_span"]:,} rivers on the busiest thread')
    print(f'  sequential main stem:   {summary["main_stem"]:,} ({summary["main_stem_fraction"]:.1%})')
    print(f'  predicted speedup:      {summary["predicted_speedup"]:.2f}x (work only; ignores memory bandwidth)')
    return summary


def subset_network_to_river(
    target_river: int,
    network_file: PathInput,
    out_network_file: PathInput,
    weights_file: PathInput | None = None,
    out_weights_file: PathInput | None = None,
) -> None:
    """
    Subset a network file, and optionally its grid weights, to the target river and every river upstream of it.

    In DFS order a river's upstream watershed is the run of rows ending at it, so the subset is the rows whose
    riverIndex lies from the target's riverIndex minus its upstreamCount to its riverIndex. The target river becomes
    the outlet (nextRiverId set to -1) of the subset. Every river keeps its riverIndex and upstreamCount, which
    describe the subset as well as the whole network. With weight table paths the weight table is also filtered to
    the rivers of the subset.

    Args:
        target_river: riverId of the river to subset to (becomes the outlet)
        network_file: path to the full network file parquet
        out_network_file: path to write the subset network file parquet
        weights_file: path to the full grid weights netCDF file (optional)
        out_weights_file: path to write the subset grid weights netCDF file, given exactly when ``weights_file`` is

    Raises:
        ValueError: if the target river is not in the network file, or only one of the weights files is given
    """
    if (weights_file is None) != (out_weights_file is None):
        raise ValueError('weights_file and out_weights_file must be given together')
    table = pd.read_parquet(network_file)
    target = table['riverId'].to_numpy() == target_river
    if not target.any():
        raise ValueError(f'riverId {target_river} is not in the network file: {network_file}')
    last = int(table['riverIndex'].to_numpy()[target][0])
    first = last - int(table['upstreamCount'].to_numpy()[target][0])
    subset = table[table['riverIndex'].between(first, last)].copy()
    subset.loc[subset['riverId'] == target_river, 'nextRiverId'] = -1
    subset.to_parquet(out_network_file, index=False)

    if weights_file is not None:
        with xr.open_dataset(weights_file) as ds:
            ds.isel(index=np.isin(ds['riverId'].values, subset['riverId'].to_numpy())).to_netcdf(out_weights_file)
    return
