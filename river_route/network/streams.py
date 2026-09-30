"""
Static analysis of a network table (riverId, nextRiverId, muskingumK, muskingumX).

The analysis answers three questions for a given routing time step dt:
1. Is the network connectivity valid (unique ids, downstreams exist, topologically sorted, no cycles)?
2. For Muskingum routing, are the k/x parameters numerically stable for dt?
3. If not, how many substeps must each river too long for dt be divided into, or how many subcycles must each river
   too short for it be routed in, and therefore how much work does the stabilized network with the fewest possible
   stability errors take?

Stability comes from requiring non-negative Muskingum coefficients (see Network.stability_window):
    c1 >= 0  <=>  dt >= 2*k*x          (else "too long":  river travel time too long for dt)
    c3 >= 0  <=>  dt <= 2*k*(1-x)       (else "too short": river travel time too short for dt)
so a single reach is stable when  2*k*x <= dt <= 2*k*(1-x).

A river too long for dt is routed in N substeps: N equal sub-reaches in series, each with k' = k/N and the same x.
Substituting k -> k/N into the inequalities gives the integer window of valid substeps:
    N >= 2*k*x / dt          (lower bound, from c1 >= 0)
    N <= 2*k*(1-x) / dt      (upper bound, from c3 >= 0)
The fewest-reaches choice is the smallest valid N, i.e. N_lo = ceil(2*k*x/dt), provided N_lo <= floor(2*k*(1-x)/dt).
If that window is empty the reach cannot be made stable by substeps (it is too short for dt, or x is so close to 0.5
that no integer fits the window). A river too short for dt is instead routed in S subcycles: S equal steps of its own
within each dt, overriding dt with the shorter routing step dt/S.

The module also divides a network into the sub-watershed blocks threads route, and subsets a network file (and its
grid weights) to one river and every river upstream of it. Both only read the contiguous run of rows that a river's
riverIndex and upstreamCount give.
"""

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import xarray as xr

from ..types import Int32Array, IntArray, JobBlocks, PathInput
from ._numba_kernels import pack_blocks_into_jobs

if TYPE_CHECKING:
    from .Network import Network

__all__ = [
    'connectivity_is_valid',
    'divisors_of',
    'required_substep_reaches',
    'assign_stable_dt',
    'analyze_dt_assignment',
    'optimize_network_compute',
    'stable_static_network',
    'expand_network',
    'analyze_min_compute',
    'analyze_stability',
    'assign_blocks',
    'analyze_partitioning',
    'subset_network_to_river',
]

# the caps on the rivers in a block assign_blocks tries, as fractions of an even share of the rivers per job
_BLOCK_CAP_FRACTIONS: tuple[float, ...] = (0.125, 0.25, 0.375, 0.5, 0.65, 0.8, 1.0)


def connectivity_is_valid(df: pd.DataFrame) -> bool:
    """
    Check that river ids are unique, all downstreams exist as rivers except -1, and that the table is
    topologically sorted from upstream to downstream (which also rules out cycles). Prints the first problem found.

    Args:
        df: network table with columns riverId and nextRiverId

    Returns:
        bool: True when the connectivity is valid
    """
    total_rows = df.shape[0]
    unique_rivers = df['riverId'].nunique()
    if unique_rivers != total_rows:
        print(f'riverId column must be unique: {total_rows} rows, {unique_rivers} unique river ids')
        return False

    downstreams = set(df['nextRiverId'])
    river_ids = set(df['riverId'])
    downstreams_not_in_rivers = downstreams - river_ids
    if downstreams_not_in_rivers != {-1}:  # only -1 doesn't need to be in river_ids
        print(
            f'all nextRiverId values must be in riverId column or -1. This might have been intentional or '
            f'may indicate a hole in the topology. Downstream ids not in river ids: {downstreams_not_in_rivers}'
        )
        return False

    # check that the rivers are topologically sorted from upstream to downstream
    river_id_to_index = {river_id: idx for idx, river_id in enumerate(df['riverId'])}
    river_id_to_downstream_id = dict(zip(df['riverId'], df['nextRiverId'], strict=True))
    for river_id in df['riverId']:
        downstream_id = river_id_to_downstream_id[river_id]
        if downstream_id == -1:
            continue
        if river_id_to_index[downstream_id] <= river_id_to_index[river_id]:
            print(
                f'Rivers must be topologically sorted in the input file (upstream before downstream). '
                f'River {river_id} has downstream {downstream_id} which appears before it in the file.'
            )
            return False

    return True


def required_substep_reaches(k: np.ndarray, x: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """
    For each river, compute its substeps: the smallest number of equal sub-reaches that makes it Muskingum-stable for
    dt, using the closed-form valid window  ceil(2kx/dt) <= N <= floor(2k(1-x)/dt).

    Args:
        k: array of Muskingum k values (travel time, same units as dt)
        x: array of Muskingum x weighting factors (0 <= x <= 0.5)
        dt: routing time step

    Returns:
        substeps: integer array of sub-reaches per river. For resolvable rivers this is the smallest valid N
            (1 when the river is already stable). For unresolvable rivers it is 1 (the reach is kept as-is and
            remains an error).
        resolvable: boolean array, True where substeps yield stable sub-reaches.
    """
    k = np.asarray(k, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)

    # smallest N satisfying the c1 >= 0 (too-long) bound; at least 1 reach must always exist
    fewest_substeps = np.maximum(1, np.ceil(2 * k * x / dt)).astype(np.int64)
    # largest N still satisfying the c3 >= 0 (too-short) bound; smaller k/N reduces this, so more substeps can break it
    most_substeps = np.floor(2 * k * (1 - x) / dt).astype(np.int64)

    resolvable = fewest_substeps <= most_substeps
    substeps = np.where(resolvable, fewest_substeps, 1)
    return substeps, resolvable


def divisors_of(n: int) -> np.ndarray:
    """All positive integer divisors of n, sorted ascending."""
    small = [d for d in range(1, int(n**0.5) + 1) if n % d == 0]
    return np.array(sorted(set(small + [n // d for d in small])), dtype=np.int64)


def assign_stable_dt(k: np.ndarray, x: np.ndarray, period: int = 3600) -> tuple[np.ndarray, np.ndarray]:
    """
    For each river pick a routing time step dt that (a) evenly divides ``period`` and (b) keeps the river
    Muskingum-stable, i.e. 2*k*x <= dt <= 2*k*(1-x). The largest valid divisor is chosen so the river takes the
    fewest subcycles (period/dt) per outer step. Unlike substeps, which divide a river into sub-reaches, lowering dt
    can also resolve "too short" reaches, so the only unresolvable rivers are those whose stability window contains
    no divisor of ``period``.

    Args:
        k: array of Muskingum k values (travel time, same units as period)
        x: array of Muskingum x weighting factors (0 <= x <= 0.5)
        period: outer time step that each per-river dt must divide evenly (default 3600 s = 1 hour)

    Returns:
        dt: integer array of chosen time steps. For resolvable rivers this is the largest valid divisor of
            ``period``; for unresolvable rivers it is 0.
        resolvable: boolean array, True where a valid divisor exists in the stability window.
    """
    k = np.asarray(k, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)

    divisors = divisors_of(period)
    window_lo = 2 * k * x
    window_hi = 2 * k * (1 - x)

    # largest divisor <= window_hi (index of last divisor not exceeding the upper bound)
    idx = np.searchsorted(divisors, window_hi, side='right') - 1
    candidate = np.where(idx >= 0, divisors[np.clip(idx, 0, len(divisors) - 1)], 0)

    resolvable = (idx >= 0) & (candidate >= window_lo)
    dt = np.where(resolvable, candidate, 0).astype(np.int64)
    return dt, resolvable


def analyze_dt_assignment(df: pd.DataFrame, period: int = 3600) -> dict:
    """
    Static analysis of assigning each river its own stable, period-dividing time step (see assign_stable_dt).

    Reports how many rivers are resolvable, the distribution of chosen dt, the total subcycles the network would
    take per outer ``period``, and how many rivers have no admissible dt. Also prints a human-readable report.

    Args:
        df: network table with columns muskingumK and muskingumX
        period: outer time step in seconds each per-river dt must divide evenly

    Returns:
        dict: the summary statistics
    """
    k = df['muskingumK'].to_numpy()
    x = df['muskingumX'].to_numpy()
    dt, resolvable = assign_stable_dt(k, x, period)

    n_rivers = df.shape[0]
    subcycles = np.where(resolvable, period // np.where(dt > 0, dt, 1), 0)
    dt_counts = (
        pd.Series(dt[resolvable]).value_counts().sort_index(ascending=False)
        if resolvable.any()
        else pd.Series(dtype=np.int64)
    )

    summary = {
        'period': period,
        'n_rivers': n_rivers,
        'n_resolvable': int(resolvable.sum()),
        'n_unresolvable': int((~resolvable).sum()),
        'total_subcycles': int(subcycles.sum()),
        'min_dt': int(dt[resolvable].min()) if resolvable.any() else 0,
        'max_dt': int(dt[resolvable].max()) if resolvable.any() else 0,
        'dt_counts': dt_counts.to_dict(),
    }

    print(f'Per-river dt assignment for period={period}')
    print(f'  rivers in:           {summary["n_rivers"]:,}')
    print(f'  resolvable:          {summary["n_resolvable"]:,} (dt range {summary["min_dt"]}..{summary["max_dt"]} s)')
    print(f'  unresolvable:        {summary["n_unresolvable"]:,} (no divisor of {period} fits the window)')
    print(f'  total subcycles/{period}s: {summary["total_subcycles"]:,}')
    return summary


def optimize_network_compute(k: np.ndarray, x: np.ndarray, period: int = 3600, cap: int | None = 10) -> dict:
    """
    For each river choose the combination of substeps N (equal sub-reaches in series) and subcycles S (equal steps of
    its own per period) that keeps every sub-reach Muskingum-stable while minimizing compute cost N*S (each sub-reach
    routed once per subcycle is one compute cycle). Both levers are searched jointly over the divisors of ``period``:

    With N substeps (k' = k/N) and S subcycles (dt' = period/S), stability of each sub-reach requires
        ceil(2*k*x*S / period)  <=  N  <=  floor(2*k*(1-x)*S / period).
    For each candidate S (a divisor of ``period``) the cheapest feasible N is the lower bound, giving cost N*S.
    The minimum over all S is the per-river optimum. This collapses to substeps alone (S=1) for "too long" reaches
    and subcycles alone (N=1) for "too short" reaches; combining only resolves narrow (x near 0.5) windows.

    Args:
        k: per-river Muskingum k in seconds
        x: per-river Muskingum x
        period: outer time step the per-river dt must divide (default 3600 s)
        cap: maximum allowed substeps N and subcycles S (default 10). Rivers needing more of either to be stable are
            reported as unresolvable. Pass None for no cap.

    Returns a dict of arrays: substeps, subcycles, dt (= period/subcycles), cost (= substeps*subcycles, 0 if
    unresolvable), and resolvable.
    """
    k = np.asarray(k, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    divisors = divisors_of(period)
    if cap is not None:
        divisors = divisors[divisors <= cap]

    inf = np.iinfo(np.int64).max
    best_cost = np.full(k.shape, inf, dtype=np.int64)
    best_substeps = np.ones(k.shape, dtype=np.int64)
    best_subcycles = np.ones(k.shape, dtype=np.int64)
    resolvable = np.zeros(k.shape, dtype=bool)

    for subcycles in divisors:
        dt = period / subcycles
        fewest_substeps = np.maximum(1, np.ceil(2 * k * x / dt)).astype(np.int64)
        most_substeps = np.floor(2 * k * (1 - x) / dt).astype(np.int64)
        within_cap = fewest_substeps <= cap if cap is not None else True
        feasible = (fewest_substeps <= most_substeps) & within_cap
        cost = fewest_substeps * subcycles
        improve = feasible & (cost < best_cost)
        best_cost = np.where(improve, cost, best_cost)
        best_substeps = np.where(improve, fewest_substeps, best_substeps)
        best_subcycles = np.where(improve, subcycles, best_subcycles)
        resolvable |= feasible

    cost = np.where(resolvable, best_cost, 0)
    return {
        'substeps': best_substeps,
        'subcycles': best_subcycles,
        'dt': (period // best_subcycles).astype(np.int64),
        'cost': cost,
        'resolvable': resolvable,
    }


def stable_static_network(df: pd.DataFrame, period: int = 3600, cap: int | None = 10) -> pd.DataFrame:
    """
    Annotate a network table with the two stability levers a stability-aware kernel needs, without
    changing the network topology or adding any rows. Each original river keeps its row, k, x, and connectivity;
    two independent-axis columns are added describing how to route that river to stay Muskingum-stable for ``period``.

    The two levers are orthogonal and both come straight from optimize_network_compute:

        substeps (N): route the river as N equal sub-reaches in SERIES, each with k/N and catchment runoff/N and the
            same x (for "too long" reaches). The reported flow is the instantaneous outflow of the final sub-reach.
            N == 1 means one reach.
        subcycles (S): route the river S times at its own dt = period/S (for "too short" reaches) and report the
            average of the S subcycle outflows. S == 1 means a single route.

    The optimizer's least-cost choice is always a single lever (N>1 xor S>1), but keeping them as separate columns
    expresses each axis explicitly, lets a kernel treat "report instantaneous outlet" as the S==1 case of "average
    over S", and leaves room for the rare reach that needs both. There are no synthetic river ids in the file (the
    substeps are materialized only in memory by expand_network); the file is exactly the input plus two columns.
    Rivers with no stable (N, S) within ``cap`` are left at N=S=1 and remain an error; see the ``stable`` column.

    Args:
        df: network table with columns riverId, nextRiverId, muskingumK, muskingumX (extra columns are preserved)
        period: outer time step each per-river dt must divide (default 3600 s)
        cap: maximum substeps or subcycles per river (default 10); see optimize_network_compute

    Returns:
        A copy of df with added columns:
            substeps  -- int, equal sub-reaches in series; k and catchment runoff are divided by it
            subcycles -- int, steps of its own to route and average over
            stable    -- bool, False where no stable routing exists within cap (kept at N=S=1, an error)
    """
    opt = optimize_network_compute(df['muskingumK'].to_numpy(), df['muskingumX'].to_numpy(), period, cap=cap)
    out = df.copy()
    # int32 (not int16): with cap=None a very long reach can need > 32767 substeps and would wrap negative
    out['substeps'] = opt['substeps'].astype(np.int32)
    out['subcycles'] = opt['subcycles'].astype(np.int32)
    out['stable'] = opt['resolvable']
    return out


def expand_network(df: pd.DataFrame, period: int = 3600, cap: int | None = 10) -> dict:
    """
    Materialize the stability-expanded routing network in memory as flat, CSR-indexed arrays ready for a routing
    kernel. Every river is replaced by a contiguous block of ``substeps`` reaches in series (the expanded
    reaches follow the inlet); a river routed only in subcycles keeps a single reach. k and the catchment runoff
    scale are pre-divided by substeps so the kernel never needs the substep count.

    The arrays are laid out in original (topological) order: river i's reaches occupy reaches[indptr[i]:indptr[i+1]],
    the first is its inlet and the last is its outlet (the only reach whose discharge is reported for river i).
    Connectivity stays a strictly-lower-triangular DAG: within a block each reach feeds the next, and the outlet
    feeds the inlet (head) of its original downstream river's block.

    Args:
        df: network table with riverId, nextRiverId, muskingumK, muskingumX. If substeps/subcycles columns are not
            present they are computed via stable_static_network(df, period, cap).
        period: outer time step each per-river dt must divide (default 3600 s), passed to stable_static_network
        cap: maximum substeps or subcycles per river (default 10), passed to stable_static_network

    Returns a dict of arrays (n = expanded reach count, m = original river count):
        n_reaches      -- int, total expanded reaches
        n_rivers       -- int, original river count
        k              -- float32 (n,), per-reach k (k_river / substeps)
        x              -- float32 (n,), per-reach x (unchanged within a river)
        catchment_runoff_scale -- float32 (n,), per-reach catchment runoff multiplier (1 / substeps)
        downstream_index -- int32 (n,), expanded downstream reach index, -1 at the network outlet
        parent_index   -- int32 (n,), original river index a reach belongs to (for runoff lookup / output grouping)
        reach_river_id -- int32 (n,), the original river id R each reach belongs to (first identity column)
        subreach_number -- int32 (n,), 0 at the outlet, 1..substeps-1 upstream (second identity column); the
                           deterministic (reach_river_id, subreach_number) pair identifies a reach for state I/O
        reach_indptr   -- int32 (m+1,), CSR offsets: river i's reaches are [reach_indptr[i], reach_indptr[i+1])
        outlet_index   -- int32 (m,), reach index whose discharge is reported for each original river
        subcycles_per_reach -- int32 (n,), subcycles to route+average for each reach (the sub-reaches of a river share
                               its subcycles; that is 1 unless the river needs both levers)
        subcycles      -- int32 (m,), subcycles of each original river
        river_id       -- int32 (m,), original river ids in order (for labeling output)
        stable         -- bool (m,), per original river, False if no stable routing within cap
    """
    if 'substeps' not in df.columns or 'subcycles' not in df.columns:
        df = stable_static_network(df, period=period, cap=cap)

    substeps = df['substeps'].to_numpy(dtype=np.int64)
    subcycles = df['subcycles'].to_numpy(dtype=np.int64)
    orig_id = df['riverId'].to_numpy(dtype=np.int32)
    orig_down = df['nextRiverId'].to_numpy(dtype=np.int32)
    k = df['muskingumK'].to_numpy(dtype=np.float64)
    x = df['muskingumX'].to_numpy(dtype=np.float64)
    stable = df['stable'].to_numpy() if 'stable' in df.columns else np.ones(orig_id.shape[0], dtype=bool)
    n_orig = orig_id.shape[0]

    # the kernel cannot detect a malformed network, so enforce its preconditions here at the build boundary
    if np.any(substeps < 1):
        raise ValueError('substeps must be >= 1 for every river')
    if np.any(subcycles < 1):
        raise ValueError('subcycles must be >= 1 for every river')
    total = int(substeps.sum())

    # CSR offsets: each river's block of substeps reaches, laid out in topological order
    reach_indptr = np.empty(n_orig + 1, dtype=np.int64)
    reach_indptr[0] = 0
    np.cumsum(substeps, out=reach_indptr[1:])
    first_reach = reach_indptr[:-1]
    outlet_index = reach_indptr[1:] - 1

    # per-reach arrays via run-length expansion of the per-river values
    parent_index = np.repeat(np.arange(n_orig, dtype=np.int64), substeps)
    reach_river_id = np.repeat(orig_id, substeps)
    k_exp = np.repeat(k / substeps, substeps)
    x_exp = np.repeat(x, substeps)
    catchment_runoff_scale = np.repeat(1.0 / substeps, substeps)
    block_size = np.repeat(substeps, substeps)
    pos_in_block = np.arange(total) - np.repeat(first_reach, substeps)
    is_outlet = pos_in_block == (block_size - 1)
    # deterministic reach identity is the two-column pair (reach_river_id, subreach_number): the outlet is
    # subreach_number 0 and the n-1 upstream subreaches are numbered 1..n-1 by distance upstream. Re-expanding the
    # same (river, substeps) always yields identical pairs, so saved per-reach state can be matched back.
    subreach_number = (block_size - 1 - pos_in_block).astype(np.int32)

    # connectivity: non-outlet reaches feed the next reach; outlets feed the head of the downstream river's block
    downstream_index = np.empty(total, dtype=np.int64)
    non_outlet_idx = np.nonzero(~is_outlet)[0]
    downstream_index[non_outlet_idx] = non_outlet_idx + 1
    id_to_index = pd.Series(np.arange(n_orig, dtype=np.int64), index=orig_id)
    down_orig_index = np.full(n_orig, -1, dtype=np.int64)
    has_down = orig_down != -1
    mapped = id_to_index.reindex(orig_down[has_down]).to_numpy()
    if np.isnan(mapped).any():  # a nextRiverId absent from riverId is a topology hole, not a valid -1 outlet
        raise ValueError(
            'nextRiverId values reference ids not present in riverId (topology hole); '
            'validate with connectivity_is_valid before expanding'
        )
    down_orig_index[has_down] = mapped.astype(np.int64)
    outlet_downstream = np.where(down_orig_index >= 0, first_reach[np.clip(down_orig_index, 0, n_orig - 1)], -1)
    downstream_index[is_outlet] = outlet_downstream  # outlets are emitted in original order

    # the kernel routes reaches in array order and requires a topologically sorted DAG (each reach feeds a later one
    # or -1); an unsorted input would silently push flow into an already-finalized reach and lose mass
    if not np.all((downstream_index < 0) | (downstream_index > np.arange(total))):
        raise ValueError(
            'input rivers must be topologically sorted upstream-before-downstream '
            '(nextRiverId must appear after riverId); run connectivity_is_valid'
        )

    return {
        'n_reaches': total,
        'n_rivers': n_orig,
        'k': k_exp.astype(np.float32),
        'x': x_exp.astype(np.float32),
        'catchment_runoff_scale': catchment_runoff_scale.astype(np.float32),
        'downstream_index': downstream_index.astype(np.int32),
        'parent_index': parent_index.astype(np.int32),
        'reach_river_id': reach_river_id,
        'subreach_number': subreach_number,
        'reach_indptr': reach_indptr.astype(np.int32),
        'outlet_index': outlet_index.astype(np.int32),
        'subcycles_per_reach': np.repeat(subcycles, substeps).astype(np.int32),
        'subcycles': subcycles.astype(np.int32),
        'river_id': orig_id,
        'stable': stable,
    }


def analyze_min_compute(df: pd.DataFrame, period: int = 3600, cap: int | None = 10) -> dict:
    """
    Static analysis choosing, per river, the combination of substeps and subcycles that minimizes compute cost (see
    optimize_network_compute). Reports the total compute cycles per ``period`` for the optimized network versus
    the base case of one cycle per river, and how each river is resolved. Also prints a human-readable report.

    Args:
        df: network table with columns muskingumK and muskingumX
        period: outer time step in seconds each per-river dt must divide evenly
        cap: maximum substeps or subcycles per river, or None for no cap

    Returns:
        dict: the summary statistics
    """
    res = optimize_network_compute(df['muskingumK'].to_numpy(), df['muskingumX'].to_numpy(), period, cap=cap)
    substeps, subcycles, cost, resolvable = res['substeps'], res['subcycles'], res['cost'], res['resolvable']

    n_rivers = df.shape[0]
    stable = resolvable & (substeps == 1) & (subcycles == 1)
    via_subcycles = resolvable & (substeps == 1) & (subcycles > 1)
    via_substeps = resolvable & (substeps > 1) & (subcycles == 1)
    via_both = resolvable & (substeps > 1) & (subcycles > 1)

    summary = {
        'period': period,
        'n_rivers': n_rivers,
        'n_stable': int(stable.sum()),
        'n_via_subcycles': int(via_subcycles.sum()),
        'n_via_substeps': int(via_substeps.sum()),
        'n_via_both': int(via_both.sum()),
        'n_unresolvable': int((~resolvable).sum()),
        'base_cost': n_rivers,
        'total_cost': int(cost.sum()),
        'extra_cost': int(cost.sum()) - int(resolvable.sum()),
        'new_subreaches': int((substeps[resolvable] - 1).sum()),
        'max_cost': int(cost.max()) if n_rivers else 0,
    }

    print(f'Minimum-compute analysis for period={period} (cost = substeps x subcycles per period)')
    print(f'  rivers in:              {summary["n_rivers"]:,}')
    print(f'  stable as-is:           {summary["n_stable"]:,}')
    print(f'  fixed by subcycles:     {summary["n_via_subcycles"]:,}')
    print(f'  fixed by substeps:      {summary["n_via_substeps"]:,}')
    print(f'  fixed by both:          {summary["n_via_both"]:,}')
    print(f'  unresolvable:           {summary["n_unresolvable"]:,}')
    print(f'  base cost (1/river):    {summary["base_cost"]:,}')
    print(f'  optimized total cost:   {summary["total_cost"]:,} cycles/period')
    print(f'  extra cost over base:   {summary["extra_cost"]:,}')
    print(f'  new sub-reaches to make:{summary["new_subreaches"]:,}')
    return summary


def analyze_stability(df: pd.DataFrame, dt: float) -> dict:
    """
    Static stability analysis of a network table for a given dt.

    Reports how many rivers are already stable, how many can be fixed by substeps (and into how many sub-reaches),
    how many cannot be fixed by substeps, and the total number of reaches the fewest-errors fixed network would
    contain.

    Args:
        df: network table with columns muskingumK and muskingumX
        dt: routing time step in seconds

    Returns:
        dict: the summary statistics. Also prints a human-readable report.
    """
    k = df['muskingumK'].to_numpy()
    x = df['muskingumX'].to_numpy()
    substeps, resolvable = required_substep_reaches(k, x, dt)

    n_rivers = df.shape[0]
    already_stable = resolvable & (substeps == 1)
    resolved_by_substeps = resolvable & (substeps > 1)
    unresolvable = ~resolvable

    total_reaches = int(substeps.sum())
    summary = {
        'dt': dt,
        'n_rivers': n_rivers,
        'n_already_stable': int(already_stable.sum()),
        'n_resolved_by_substeps': int(resolved_by_substeps.sum()),
        'n_unresolvable': int(unresolvable.sum()),
        'total_reaches': total_reaches,
        'reaches_to_create': total_reaches - n_rivers,
        'max_substeps': int(substeps.max()) if n_rivers else 0,
    }

    print(f'Static stability analysis for dt={dt}')
    print(f'  rivers in:              {summary["n_rivers"]:,}')
    print(f'  already stable:         {summary["n_already_stable"]:,}')
    print(
        f'  fixable by substeps:    {summary["n_resolved_by_substeps"]:,} '
        f'(up to {summary["max_substeps"]} sub-reaches each)'
    )
    print(
        f'  unresolvable errors:    {summary["n_unresolvable"]:,} (too short for dt or window empty; kept as 1 reach)'
    )
    print(f'  total reaches in fixed network: {summary["total_reaches"]:,}')
    print(f'  new reaches to create:          {summary["reaches_to_create"]:,}')
    return summary


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
            ds.isel(index=np.isin(ds['river_id'].values, subset['riverId'].to_numpy())).to_netcdf(out_weights_file)
    return
