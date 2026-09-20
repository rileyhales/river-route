"""
Static analysis of a routing parameter table (river_id, next_river_id, k, x).

The analysis answers three questions for a given routing time step dt:
1. Is the network connectivity valid (unique ids, downstreams exist, topologically sorted, no cycles)?
2. For Muskingum routing, are the k/x parameters numerically stable for dt?
3. If not, how many sub-reaches must each unstable river be split into, and therefore how many total
   streams are required to build the "fixed" network that has the fewest possible stability errors?

Stability comes from requiring non-negative Muskingum coefficients (see Muskingum._set_muskingum_coefficients):
    c1 >= 0  <=>  dt >= 2*k*x          (else "too long":  river travel time too long for dt)
    c3 >= 0  <=>  dt <= 2*k*(1-x)       (else "too short": river travel time too short for dt)
so a single reach is stable when  2*k*x <= dt <= 2*k*(1-x).

A reach can be modeled as N equal sub-reaches in series, each with k' = k/N and the same x. Substituting
k -> k/N into the inequalities gives the integer window of valid split counts:
    N >= 2*k*x / dt          (lower bound, from c1 >= 0)
    N <= 2*k*(1-x) / dt      (upper bound, from c3 >= 0)
The fewest-streams choice is the smallest valid N, i.e. N_lo = ceil(2*k*x/dt), provided N_lo <= floor(2*k*(1-x)/dt).
If that window is empty the reach cannot be made stable by uniform splitting (it is "too short" / Case 2, or x is
so close to 0.5 that no integer fits the window) and it is kept as a single reach that remains an error.

The module also holds the graph utilities that operate on the same table: the NetworkX and sparse adjacency
representations of the connectivity and subsetting a parameter table (and its grid weights) to one river.
"""

import heapq
import logging

import networkx as nx
import numba
import numpy as np
import pandas as pd
import scipy
import xarray as xr

from ..types import PathInput

logger = logging.getLogger(__name__)

__all__ = [
    'connectivity_is_valid',
    'required_subreaches',
    'assign_stable_dt',
    'analyze_dt_assignment',
    'optimize_network_compute',
    'stable_static_network',
    'expand_network',
    'broadcast_state_to_reaches',
    'analyze_min_compute',
    'analyze_stability',
    'is_dfs_ordered',
    'shreve_order',
    'assign_regions',
    'partition_network',
    'analyze_partitioning',
    'regions_to_layout',
    'subset_configs_to_river',
    'connectivity_to_digraph',
    'adjacency_matrix',
]


def connectivity_is_valid(df: pd.DataFrame) -> bool:
    """
    checks that river ids are unique, all downstreams exist as rivers except -1, and that the table is
    topologically sorted from upstream to downstream (which also rules out cycles).
    """
    total_rows = df.shape[0]
    unique_rivers = df['river_id'].nunique()
    if unique_rivers != total_rows:
        print(f'river_id column must be unique: {total_rows} rows, {unique_rivers} unique river ids')
        return False

    downstreams = set(df['next_river_id'])
    river_ids = set(df['river_id'])
    downstreams_not_in_rivers = downstreams - river_ids
    if downstreams_not_in_rivers != {-1}:  # only -1 doesn't need to be in river_ids
        print(
            f'all next_river_id values must be in river_id column or -1. This might have been intentional or '
            f'may indicate a hole in the topology. Downstream ids not in river ids: {downstreams_not_in_rivers}'
        )
        return False

    # check that the rivers are topologically sorted from upstream to downstream
    river_id_to_index = {river_id: idx for idx, river_id in enumerate(df['river_id'])}
    river_id_to_downstream_id = dict(zip(df['river_id'], df['next_river_id'], strict=True))
    for river_id in df['river_id']:
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


def required_subreaches(k: np.ndarray, x: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """
    For each river, compute the smallest number of equal sub-reaches that makes it Muskingum-stable for dt,
    using the closed-form valid window  ceil(2kx/dt) <= N <= floor(2k(1-x)/dt).

    Args:
        k: array of Muskingum k values (travel time, same units as dt)
        x: array of Muskingum x weighting factors (0 <= x <= 0.5)
        dt: routing time step

    Returns:
        n_subreaches: integer array of sub-reaches per river. For resolvable rivers this is the smallest valid N
            (1 when the river is already stable). For unresolvable rivers it is 1 (the reach is kept as-is and
            remains an error).
        resolvable: boolean array, True where uniform splitting yields a stable reach.
    """
    k = np.asarray(k, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)

    # smallest N satisfying the c1 >= 0 (too-long) bound; at least 1 reach must always exist
    n_lo = np.maximum(1, np.ceil(2 * k * x / dt)).astype(np.int64)
    # largest N still satisfying the c3 >= 0 (too-short) bound; smaller k/N reduces this, so splitting can break it
    n_hi = np.floor(2 * k * (1 - x) / dt).astype(np.int64)

    resolvable = n_lo <= n_hi
    n_subreaches = np.where(resolvable, n_lo, 1)
    return n_subreaches, resolvable


def _divisors(n: int) -> np.ndarray:
    """All positive integer divisors of n, sorted ascending."""
    small = [d for d in range(1, int(n**0.5) + 1) if n % d == 0]
    return np.array(sorted(set(small + [n // d for d in small])), dtype=np.int64)


def assign_stable_dt(k: np.ndarray, x: np.ndarray, period: int = 3600) -> tuple[np.ndarray, np.ndarray]:
    """
    For each river pick a routing time step dt that (a) evenly divides ``period`` and (b) keeps the river
    Muskingum-stable, i.e. 2*k*x <= dt <= 2*k*(1-x). The largest valid divisor is chosen so the river takes the
    fewest substeps (period/dt) per outer step. Unlike sub-reach splitting, lowering dt can also resolve
    "too short" (Case 2) reaches, so the only unresolvable rivers are those whose stability window contains no
    divisor of ``period``.

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

    divisors = _divisors(period)
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

    Reports how many rivers are resolvable, the distribution of chosen dt, the total substeps the network would
    take per outer ``period``, and how many rivers have no admissible dt. Also prints a human-readable report.
    """
    k = df['k'].to_numpy()
    x = df['x'].to_numpy()
    dt, resolvable = assign_stable_dt(k, x, period)

    n_rivers = df.shape[0]
    substeps = np.where(resolvable, period // np.where(dt > 0, dt, 1), 0)
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
        'total_substeps': int(substeps.sum()),
        'min_dt': int(dt[resolvable].min()) if resolvable.any() else 0,
        'max_dt': int(dt[resolvable].max()) if resolvable.any() else 0,
        'dt_counts': dt_counts.to_dict(),
    }

    print(f'Per-river dt assignment for period={period}')
    print(f'  rivers in:           {summary["n_rivers"]:,}')
    print(f'  resolvable:          {summary["n_resolvable"]:,} (dt range {summary["min_dt"]}..{summary["max_dt"]} s)')
    print(f'  unresolvable:        {summary["n_unresolvable"]:,} (no divisor of {period} fits the window)')
    print(f'  total substeps/{period}s: {summary["total_substeps"]:,}')
    return summary


def optimize_network_compute(k: np.ndarray, x: np.ndarray, period: int = 3600, cap: int = 10) -> dict:
    """
    For each river choose the combination of equal sub-reach count N and substeps-per-period S that keeps every
    sub-reach Muskingum-stable while minimizing compute cost N*S (each sub-reach routed once per substep is one
    compute cycle). Both levers are searched jointly over the divisors of ``period``:

    With N sub-reaches (k' = k/N) and substeps S (dt' = period/S), stability of each sub-reach requires
        ceil(2*k*x*S / period)  <=  N  <=  floor(2*k*(1-x)*S / period).
    For each candidate S (a divisor of ``period``) the cheapest feasible N is the lower bound, giving cost N*S.
    The minimum over all S is the per-river optimum. This collapses to pure splitting (S=1) for "too long" reaches
    and pure substepping (N=1) for "too short" reaches; combining only resolves narrow (x near 0.5) windows.

    Args:
        k, x: per-river Muskingum parameters
        period: outer time step the per-river dt must divide (default 3600 s)
        cap: maximum allowed splits N and substeps S (default 10). Rivers needing more of either to be stable are
            reported as unresolvable. Pass None for no cap.

    Returns a dict of arrays: n_reaches, substeps, dt (= period/substeps), cost (= n_reaches*substeps, 0 if
    unresolvable), and resolvable.
    """
    k = np.asarray(k, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    divisors = _divisors(period)
    if cap is not None:
        divisors = divisors[divisors <= cap]

    inf = np.iinfo(np.int64).max
    best_cost = np.full(k.shape, inf, dtype=np.int64)
    best_n = np.ones(k.shape, dtype=np.int64)
    best_s = np.ones(k.shape, dtype=np.int64)
    resolvable = np.zeros(k.shape, dtype=bool)

    for s in divisors:
        dt = period / s
        n_lo = np.maximum(1, np.ceil(2 * k * x / dt)).astype(np.int64)
        n_hi = np.floor(2 * k * (1 - x) / dt).astype(np.int64)
        feasible = (n_lo <= n_hi) & (n_lo <= (cap if cap is not None else n_lo))
        cost = n_lo * s
        improve = feasible & (cost < best_cost)
        best_cost = np.where(improve, cost, best_cost)
        best_n = np.where(improve, n_lo, best_n)
        best_s = np.where(improve, s, best_s)
        resolvable |= feasible

    cost = np.where(resolvable, best_cost, 0)
    return {
        'n_reaches': best_n,
        'substeps': best_s,
        'dt': (period // best_s).astype(np.int64),
        'cost': cost,
        'resolvable': resolvable,
    }


def stable_static_network(df: pd.DataFrame, period: int = 3600, cap: int = 10) -> pd.DataFrame:
    """
    Annotate a routing parameter table with the two stability levers a stability-aware kernel needs, without
    changing the network topology or adding any rows. Each original river keeps its row, k, x, and connectivity;
    two independent-axis columns are added describing how to route that river to stay Muskingum-stable for ``period``.

    The two levers are orthogonal and both come straight from optimize_network_compute:

        subdivisions (N): route the river as N equal sub-reaches in SERIES, each with k/N and vlateral/N and the
            same x (a spatial split for "too long" reaches). The reported flow is the instantaneous outflow of the
            final sub-reach. N == 1 means no split.
        substeps (S): sub-cycle the river S times at dt = period/S (a temporal refinement for "too short" reaches)
            and report the average of the S substep outflows. S == 1 means a single route.

    The optimizer's least-cost choice is always a single lever (N>1 xor S>1), but keeping them as separate columns
    expresses each axis explicitly, lets a kernel treat "report instantaneous outlet" as the S==1 case of "average
    over S", and leaves room for the rare reach that needs both. There are no synthetic river ids in the file (the
    split is materialized only in memory by expand_network); the file is exactly the input plus two columns.
    Rivers with no stable (N, S) within ``cap`` are left at N=S=1 and remain an error; see the ``stable`` column.

    Args:
        df: parameter table with columns river_id, next_river_id, k, x (extra columns are preserved)
        period: outer time step each per-river dt must divide (default 3600 s)
        cap: maximum subdivisions or substeps per river (default 10); see optimize_network_compute

    Returns:
        A copy of df with added columns:
            subdivisions -- int, equal sub-reaches in series (spatial split); k and vlateral are divided by it
            substeps     -- int, temporal substeps to route and average over
            stable        -- bool, False where no stable routing exists within cap (kept at N=S=1, an error)
    """
    opt = optimize_network_compute(df['k'].to_numpy(), df['x'].to_numpy(), period, cap=cap)
    out = df.copy()
    # int32 (not int16): with cap=None a very long reach can need > 32767 subdivisions and would wrap negative
    out['subdivisions'] = opt['n_reaches'].astype(np.int32)
    out['substeps'] = opt['substeps'].astype(np.int32)
    out['stable'] = opt['resolvable']
    return out


def expand_network(df: pd.DataFrame, period: int = 3600, cap: int = 10) -> dict:
    """
    Materialize the stability-expanded routing network in memory as flat, CSR-indexed arrays ready for a routing
    kernel. Every river is replaced by a contiguous block of ``subdivisions`` reaches in series (the expanded
    reaches follow the inlet); a river that is only substepped keeps a single reach. k and the lateral-inflow
    scale are pre-divided by subdivisions so the kernel never needs the subdivision count.

    The arrays are laid out in original (topological) order: river i's reaches occupy reaches[indptr[i]:indptr[i+1]],
    the first is its inlet and the last is its outlet (the only reach whose discharge is reported for river i).
    Connectivity stays a strictly-lower-triangular DAG: within a block each reach feeds the next, and the outlet
    feeds the inlet (head) of its original downstream river's block.

    Args:
        df: parameter table with river_id, next_river_id, k, x. If subdivisions/substeps columns are not
            present they are computed via stable_static_network(df, period, cap).

    Returns a dict of arrays (n = expanded reach count, m = original river count):
        n_reaches      -- int, total expanded reaches
        n_rivers       -- int, original river count
        k              -- float32 (n,), per-reach k (k_river / subdivisions)
        x              -- float32 (n,), per-reach x (unchanged within a river)
        lateral_scale  -- float32 (n,), per-reach lateral multiplier (1 / subdivisions)
        downstream_index -- int32 (n,), expanded downstream reach index, -1 at the network outlet
        parent_index   -- int32 (n,), original river index a reach belongs to (for vlateral lookup / output grouping)
        reach_river_id -- int64 (n,), the original river id R each reach belongs to (first identity column)
        subreach_number -- int32 (n,), 0 at the outlet, 1..subdivisions-1 upstream (second identity column); the
                           deterministic (reach_river_id, subreach_number) pair identifies a reach for state I/O
        reach_indptr   -- int32 (m+1,), CSR offsets: river i's reaches are [reach_indptr[i], reach_indptr[i+1])
        outlet_index   -- int32 (m,), reach index whose discharge is reported for each original river
        substeps_per_reach -- int32 (n,), temporal substeps to route+average for each reach (a subdivided river's
                              reaches share its substeps; that is 1 unless the river needs both levers)
        substeps       -- int32 (m,), temporal substeps for each original river
        river_id       -- int64 (m,), original river ids in order (for labeling output)
        stable         -- bool (m,), per original river, False if no stable routing within cap
    """
    if 'subdivisions' not in df.columns or 'substeps' not in df.columns:
        df = stable_static_network(df, period=period, cap=cap)

    n_subdiv = df['subdivisions'].to_numpy(dtype=np.int64)
    substeps = df['substeps'].to_numpy(dtype=np.int64)
    orig_id = df['river_id'].to_numpy(dtype=np.int64)
    orig_down = df['next_river_id'].to_numpy(dtype=np.int64)
    k = df['k'].to_numpy(dtype=np.float64)
    x = df['x'].to_numpy(dtype=np.float64)
    stable = df['stable'].to_numpy() if 'stable' in df.columns else np.ones(orig_id.shape[0], dtype=bool)
    n_orig = orig_id.shape[0]

    # the kernel cannot detect a malformed network, so enforce its preconditions here at the build boundary
    if np.any(n_subdiv < 1):
        raise ValueError('subdivisions must be >= 1 for every river')
    if np.any(substeps < 1):
        raise ValueError('substeps must be >= 1 for every river')
    total = int(n_subdiv.sum())

    # CSR offsets: each river's block of subdivisions reaches, laid out in topological order
    reach_indptr = np.empty(n_orig + 1, dtype=np.int64)
    reach_indptr[0] = 0
    np.cumsum(n_subdiv, out=reach_indptr[1:])
    group_start = reach_indptr[:-1]
    outlet_index = reach_indptr[1:] - 1

    # per-reach arrays via run-length expansion of the per-river values
    parent_index = np.repeat(np.arange(n_orig, dtype=np.int64), n_subdiv)
    reach_river_id = np.repeat(orig_id, n_subdiv)
    k_exp = np.repeat(k / n_subdiv, n_subdiv)
    x_exp = np.repeat(x, n_subdiv)
    lateral_scale = np.repeat(1.0 / n_subdiv, n_subdiv)
    block_size = np.repeat(n_subdiv, n_subdiv)
    pos_in_block = np.arange(total) - np.repeat(group_start, n_subdiv)
    is_outlet = pos_in_block == (block_size - 1)
    # deterministic reach identity is the two-column pair (reach_river_id, subreach_number): the outlet is
    # subreach_number 0 and the n-1 upstream subreaches are numbered 1..n-1 by distance upstream. Re-expanding the
    # same (river, subdivisions) always yields identical pairs, so saved per-reach state can be matched back.
    subreach_number = (block_size - 1 - pos_in_block).astype(np.int32)

    # connectivity: non-outlet reaches feed the next reach; outlets feed the head of the downstream river's block
    downstream_index = np.empty(total, dtype=np.int64)
    non_outlet_idx = np.nonzero(~is_outlet)[0]
    downstream_index[non_outlet_idx] = non_outlet_idx + 1
    id_to_index = pd.Series(np.arange(n_orig, dtype=np.int64), index=orig_id)
    down_orig_index = np.full(n_orig, -1, dtype=np.int64)
    has_down = orig_down != -1
    mapped = id_to_index.reindex(orig_down[has_down]).to_numpy()
    if np.isnan(mapped).any():  # a next_river_id absent from river_id is a topology hole, not a valid -1 outlet
        raise ValueError(
            'next_river_id values reference ids not present in river_id (topology hole); '
            'validate with connectivity_is_valid before expanding'
        )
    down_orig_index[has_down] = mapped.astype(np.int64)
    outlet_downstream = np.where(down_orig_index >= 0, group_start[np.clip(down_orig_index, 0, n_orig - 1)], -1)
    downstream_index[is_outlet] = outlet_downstream  # outlets are emitted in original order

    # the kernel sweeps in array order and requires a topologically sorted DAG (each reach feeds a later one or -1);
    # an unsorted input would silently push flow into an already-finalized reach and lose mass
    if not np.all((downstream_index < 0) | (downstream_index > np.arange(total))):
        raise ValueError(
            'input rivers must be topologically sorted upstream-before-downstream '
            '(next_river_id must appear after river_id); run connectivity_is_valid'
        )

    return {
        'n_reaches': total,
        'n_rivers': n_orig,
        'k': k_exp.astype(np.float32),
        'x': x_exp.astype(np.float32),
        'lateral_scale': lateral_scale.astype(np.float32),
        'downstream_index': downstream_index.astype(np.int32),
        'parent_index': parent_index.astype(np.int32),
        'reach_river_id': reach_river_id,
        'subreach_number': subreach_number,
        'reach_indptr': reach_indptr.astype(np.int32),
        'outlet_index': outlet_index.astype(np.int32),
        'substeps_per_reach': np.repeat(substeps, n_subdiv).astype(np.int32),
        'substeps': substeps.astype(np.int32),
        'river_id': orig_id,
        'stable': stable,
    }


def broadcast_state_to_reaches(expanded: dict, channel_state: np.ndarray) -> np.ndarray:
    """
    Seed an expanded per-reach state array from a per-river channel state by broadcasting each river's single value
    across all of its subreaches (the order expand_network lays them out). Returns a float32 array of length
    expanded['n_reaches'] suitable as the kernel's ``q``.
    """
    channel_state = np.asarray(channel_state, dtype=np.float32)
    return channel_state[expanded['parent_index']]


def analyze_min_compute(df: pd.DataFrame, period: int = 3600, cap: int = 10) -> dict:
    """
    Static analysis choosing, per river, the split/substep combination that minimizes compute cost (see
    optimize_network_compute). Reports the total compute cycles per ``period`` for the optimized network versus
    the base case of one cycle per river, and how each river is resolved. Also prints a human-readable report.
    """
    res = optimize_network_compute(df['k'].to_numpy(), df['x'].to_numpy(), period, cap=cap)
    n_reaches, substeps, cost, resolvable = res['n_reaches'], res['substeps'], res['cost'], res['resolvable']

    n_rivers = df.shape[0]
    stable = resolvable & (n_reaches == 1) & (substeps == 1)
    via_substeps = resolvable & (n_reaches == 1) & (substeps > 1)
    via_split = resolvable & (n_reaches > 1) & (substeps == 1)
    via_combined = resolvable & (n_reaches > 1) & (substeps > 1)

    summary = {
        'period': period,
        'n_rivers': n_rivers,
        'n_stable': int(stable.sum()),
        'n_via_substeps': int(via_substeps.sum()),
        'n_via_split': int(via_split.sum()),
        'n_via_combined': int(via_combined.sum()),
        'n_unresolvable': int((~resolvable).sum()),
        'base_cost': n_rivers,
        'total_cost': int(cost.sum()),
        'extra_cost': int(cost.sum()) - int(resolvable.sum()),
        'new_subreaches': int((n_reaches[resolvable] - 1).sum()),
        'max_cost': int(cost.max()) if n_rivers else 0,
    }

    print(f'Minimum-compute analysis for period={period} (cost = reaches x substeps per period)')
    print(f'  rivers in:              {summary["n_rivers"]:,}')
    print(f'  stable as-is:           {summary["n_stable"]:,}')
    print(f'  fixed by substeps:      {summary["n_via_substeps"]:,}')
    print(f'  fixed by splitting:     {summary["n_via_split"]:,}')
    print(f'  fixed by both:          {summary["n_via_combined"]:,}')
    print(f'  unresolvable:           {summary["n_unresolvable"]:,}')
    print(f'  base cost (1/river):    {summary["base_cost"]:,}')
    print(f'  optimized total cost:   {summary["total_cost"]:,} cycles/period')
    print(f'  extra cost over base:   {summary["extra_cost"]:,}')
    print(f'  new sub-reaches to make:{summary["new_subreaches"]:,}')
    return summary


def analyze_stability(df: pd.DataFrame, dt: float) -> dict:
    """
    Static stability analysis of a routing parameter table for a given dt.

    Reports how many rivers are already stable, how many can be fixed by splitting (and into how many pieces),
    how many cannot be fixed by uniform splitting, and the total number of streams the fewest-errors fixed
    network would contain.

    Returns a dict of the summary statistics. Also prints a human-readable report.
    """
    k = df['k'].to_numpy()
    x = df['x'].to_numpy()
    n_subreaches, resolvable = required_subreaches(k, x, dt)

    n_rivers = df.shape[0]
    already_stable = resolvable & (n_subreaches == 1)
    resolved_by_split = resolvable & (n_subreaches > 1)
    unresolvable = ~resolvable

    total_streams = int(n_subreaches.sum())
    summary = {
        'dt': dt,
        'n_rivers': n_rivers,
        'n_already_stable': int(already_stable.sum()),
        'n_resolved_by_split': int(resolved_by_split.sum()),
        'n_unresolvable': int(unresolvable.sum()),
        'total_streams': total_streams,
        'streams_to_create': total_streams - n_rivers,
        'max_subreaches': int(n_subreaches.max()) if n_rivers else 0,
    }

    print(f'Static stability analysis for dt={dt}')
    print(f'  rivers in:              {summary["n_rivers"]:,}')
    print(f'  already stable:         {summary["n_already_stable"]:,}')
    print(
        f'  fixable by splitting:   {summary["n_resolved_by_split"]:,} '
        f'(up to {summary["max_subreaches"]} sub-reaches each)'
    )
    print(
        f'  unresolvable errors:    {summary["n_unresolvable"]:,} (too short for dt or window empty; kept as 1 reach)'
    )
    print(f'  total streams in fixed network: {summary["total_streams"]:,}')
    print(f'  new streams to create:          {summary["streams_to_create"]:,}')
    return summary


# ──────────────────────────────────────────────────────────────────────────────
# Network partitioning for multi-threaded routing
#
# Everything here assumes the parameter table is in DFS computation order: a depth-first walk from each basin
# outlet upstream, emitted so that every river follows all of its own upstream rivers. That ordering makes every
# subtree a CONTIGUOUS block ending at its own outlet, which is what lets a thread be handed a plain index range
# instead of a list of rivers. Nothing here ever reorders a parameter table: the order it is given in is the
# order it is routed and written in, so the forcing, the state and the routed discharge always line up with the
# file the user provided. DFS order is a property of the input, checked with is_dfs_ordered and required by
# assign_regions; a table that lacks it is reported, never rewritten.
# ──────────────────────────────────────────────────────────────────────────────


def _downstream_indices(df: pd.DataFrame) -> np.ndarray:
    """
    Map the next_river_id column to positional indices, -1 where a river has no downstream.

    Raises:
        ValueError: if a next_river_id is absent from river_id (a topology hole), or if the table is not
            sorted upstream-before-downstream (which every routine here relies on).
    """
    n = df.shape[0]
    river_id = df['river_id'].to_numpy(dtype=np.int64)
    next_river_id = df['next_river_id'].to_numpy(dtype=np.int64)
    downstream_index = np.full(n, -1, dtype=np.int64)
    has_downstream = next_river_id != -1
    mapped = pd.Series(np.arange(n, dtype=np.int64), index=river_id).reindex(next_river_id[has_downstream]).to_numpy()
    if np.isnan(mapped).any():
        raise ValueError(
            'next_river_id values reference ids not present in river_id (topology hole); '
            'validate with connectivity_is_valid before partitioning'
        )
    downstream_index[has_downstream] = mapped.astype(np.int64)
    if not np.all((downstream_index < 0) | (downstream_index > np.arange(n))):
        raise ValueError(
            'rivers must be topologically sorted upstream-before-downstream '
            '(next_river_id must appear after river_id); run connectivity_is_valid'
        )
    return downstream_index


@numba.njit(cache=True)
def _subtree_extent(downstream_index: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    One forward pass giving each river's upstream count and the lowest index in its subtree.

    Correct only because the table is topologically sorted: every contributor of river i sits at a lower index,
    so both values are final by the time the loop reaches i and can be folded into its downstream.
    """
    n = downstream_index.shape[0]
    size = np.ones(n, dtype=np.int64)
    lowest = np.arange(n)
    for i in range(n):
        d = downstream_index[i]
        if d >= 0:
            size[d] += size[i]
            if lowest[i] < lowest[d]:
                lowest[d] = lowest[i]
    return size, lowest


@numba.njit(cache=True)
def shreve_order(downstream_index: np.ndarray) -> np.ndarray:
    """
    Shreve magnitude of every river: 1 for a headwater, and the sum of its upstream rivers' magnitudes otherwise, which
    is the number of headwaters upstream of and including it. One forward pass, because the table is topologically
    sorted and every upstream river is final before its downstream is reached.
    """
    n = downstream_index.shape[0]
    magnitude = np.zeros(n, dtype=np.int64)
    for i in range(n):
        if magnitude[i] == 0:
            magnitude[i] = 1
        d = downstream_index[i]
        if d >= 0:
            magnitude[d] += magnitude[i]
    return magnitude


def is_dfs_ordered(downstream_index: np.ndarray) -> np.ndarray:
    """
    Per-river check that the table is in DFS computation order.

    A river's subtree is contiguous exactly when the lowest index it contains is ``i - upstream_count``. Any
    river failing this has some of its upstream water sitting outside its own block, so the block cannot be
    routed as an independent range. Returns the per-river boolean so callers can report how far off a table is.
    """
    size, lowest = _subtree_extent(downstream_index)
    return lowest == (np.arange(downstream_index.shape[0]) - size + 1)


@numba.njit(cache=True)
def _claim_regions(size: np.ndarray, lowest: np.ndarray, order: np.ndarray, cap: int):
    """
    Claim maximal subtrees of at most ``cap`` rivers as concurrent regions, largest first.

    In DFS computation order a subtree is just the range [lowest[v], v], so claiming one is marking a slice --
    no traversal. Visiting candidates in descending size is what makes the result safe to route concurrently.
    Any ancestor of a river is at least as large as it, so ancestors are considered first: once a river is
    unclaimed at its turn, no ancestor of it has been claimed either. That yields the two properties the kernels
    depend on:
        1. regions are disjoint contiguous ranges, each holding every river upstream of its own outlet, so a
           region needs nothing from outside its range;
        2. a region's outlet always drains into unclaimed water, never into another region, so the only
           cross-region write in a routing sweep is the single push at each region's outlet.

    Returns the per-river region id (-1 for rivers left to the sequential main stem) and the region count.
    """
    n = size.shape[0]
    region = np.full(n, -1, dtype=np.int64)
    n_regions = 0
    for oi in range(order.shape[0]):
        v = order[oi]
        if region[v] >= 0 or size[v] > cap:
            continue
        for i in range(lowest[v], v + 1):
            region[i] = n_regions
        n_regions += 1
    return region, n_regions


def _parallel_cost(region: np.ndarray, n_regions: int, threads: int) -> tuple[int, int, float]:
    """
    Predict the cost of a partition as (parallel span, sequential main stem, speedup over one thread).

    Regions are packed onto ``threads`` workers longest-first, which is how a work-stealing pool schedules them
    when the largest regions are submitted first. Cost is counted in rivers swept, since the kernel does a fixed
    amount of work per river per sweep. The main stem is added whole because it runs single-threaded afterwards.
    """
    counts = np.bincount(region[region >= 0], minlength=max(n_regions, 1))
    main_stem = int(np.count_nonzero(region < 0))
    bins = [(0, w) for w in range(threads)]
    heapq.heapify(bins)
    for count in sorted(counts.tolist(), reverse=True):
        load, w = heapq.heappop(bins)
        heapq.heappush(bins, (load + int(count), w))
    span = max(load for load, _ in bins) if threads else 0
    total = span + main_stem
    return span, main_stem, (region.shape[0] / total if total else 1.0)


# fractions of the ideal per-thread share to try as the largest allowed region; the best is chosen by
# _parallel_cost. A loose cap leaves a short main stem but packs unevenly, a tight one packs evenly but
# pushes more water into the sequential main stem, and the trade turns over at a different point per network.
_CAP_MULTIPLIERS: tuple[float, ...] = (0.25, 0.5, 0.75, 1.0, 1.3, 1.6, 2.0)


def assign_regions(
    downstream_index: np.ndarray,
    threads: int = 4,
    granularity: int = 2,
    cap_multiplier: float | None = None,
    measure: str = 'rivers',
) -> tuple[np.ndarray, int]:
    """
    Partition a DFS-ordered network into contiguous regions that can be routed concurrently (see _claim_regions).

    Args:
        downstream_index: positional index of each river's downstream, -1 where there is none
        threads: worker count the partition is being sized for
        granularity: regions to aim for per thread. More regions than threads let a work-stealing pool even out
            the uneven region sizes a river network produces; 2 is the measured sweet spot.
        cap_multiplier: largest allowed region as a multiple of the ideal per-thread share. None searches
            _CAP_MULTIPLIERS and keeps whichever partition _parallel_cost rates fastest.
        measure: what sizes a subtree when regions are claimed. 'rivers' counts the rivers in it; 'shreve' uses its
            Shreve magnitude, the number of headwaters it drains. Either never shrinks downstream, which is what makes
            claiming the largest subtrees first safe. The partitions are always rated by rivers, the work routed.

    Returns:
        region: per-river region id, -1 for rivers in the sequential main stem
        n_regions: number of concurrent regions

    Raises:
        ValueError: if the network is not in DFS computation order, since regions would not be contiguous
    """
    if threads < 1:
        raise ValueError(f'threads must be >= 1, got {threads}')
    if granularity < 1:
        raise ValueError(f'granularity must be >= 1, got {granularity}')
    n = downstream_index.shape[0]
    size, lowest = _subtree_extent(downstream_index)
    contiguous = lowest == (np.arange(n) - size + 1)
    if not contiguous.all():
        raise ValueError(
            f'{int((~contiguous).sum()):,} of {n:,} rivers do not have a contiguous subtree, so the parameter '
            f'table is not in DFS computation order and cannot be split into concurrent ranges. Provide a '
            f'parameter table sorted in depth-first computation order, or route with threads=1, which places '
            f'no ordering requirement beyond upstream-before-downstream.'
        )
    if measure == 'rivers':
        claim_size = size
    elif measure == 'shreve':
        claim_size = shreve_order(downstream_index)
    else:
        raise ValueError(f"measure must be 'rivers' or 'shreve', got {measure!r}")
    # largest first, and among equal sizes the downstream river first: a chain of rivers shares one Shreve magnitude,
    # and its most downstream river must be claimed before any subtree inside it
    order = np.lexsort((-np.arange(n), -claim_size))
    share = claim_size[downstream_index < 0].sum() / (threads * granularity)  # the whole network's size

    multipliers = _CAP_MULTIPLIERS if cap_multiplier is None else (cap_multiplier,)
    best: tuple[float, np.ndarray, int] | None = None
    for multiplier in multipliers:
        region, n_regions = _claim_regions(claim_size, lowest, order, max(1, int(share * multiplier)))
        _, _, speedup = _parallel_cost(region, n_regions, threads)
        if best is None or speedup > best[0]:
            best = (speedup, region, n_regions)
    _, region, n_regions = best
    return region, n_regions


def partition_network(
    df: pd.DataFrame, threads: int = 4, granularity: int = 2, cap_multiplier: float | None = None
) -> pd.DataFrame:
    """
    Annotate a routing parameter table with the region each river belongs to for multi-threaded routing.

    Topology and row order are unchanged; one column is added. The table must already be in DFS computation
    order, which is what makes each region a contiguous range of rows, and is a property of the input this
    function checks rather than imposes.

    Routing a partitioned network is a two-stage pass with a single barrier between them, not a barrier per time
    step: coupling between rivers only ever runs downstream, so a region can be routed for the whole simulation
    in isolation. Each region's outlet contribution is buffered per routing step and injected when the main stem
    is swept.

    The partition depends only on connectivity and ``threads``, never on the forcing, dt, or coefficients, so it
    is computed once and stored with the parameters rather than rebuilt per simulation.

    Args:
        df: parameter table with columns river_id and next_river_id (other columns are preserved)
        threads: worker count the partition is sized for
        granularity: regions to aim for per thread; see assign_regions
        cap_multiplier: largest allowed region as a multiple of the ideal per-thread share; None auto-selects

    Returns:
        A copy of df with one added column:
            region -- int32, the concurrent region a river belongs to, or -1 if it is routed in the sequential
                      main stem after the barrier
    """
    region, _ = assign_regions(_downstream_indices(df), threads, granularity, cap_multiplier)
    out = df.copy()
    out['region'] = region.astype(np.int32)
    return out


def analyze_partitioning(
    df: pd.DataFrame, threads: int = 4, granularity: int = 2, cap_multiplier: float | None = None
) -> dict:
    """
    Static analysis of partitioning a network for multi-threaded routing (see partition_network).

    Reports the region count, how evenly they pack onto ``threads`` workers, how much of the network falls into
    the sequential main stem, and the speedup that implies. The speedup is an upper bound from work alone: it
    assumes perfect scaling and ignores memory bandwidth, which this kernel is largely bound by. Also prints a
    human-readable report.
    """
    region, n_regions = assign_regions(_downstream_indices(df), threads, granularity, cap_multiplier)
    span, main_stem, speedup = _parallel_cost(region, n_regions, threads)
    counts = np.bincount(region[region >= 0], minlength=max(n_regions, 1))
    n_rivers = df.shape[0]

    summary = {
        'threads': threads,
        'n_rivers': n_rivers,
        'n_regions': n_regions,
        'largest_region': int(counts.max()) if n_regions else 0,
        'smallest_region': int(counts.min()) if n_regions else 0,
        'parallel_span': span,
        'main_stem': main_stem,
        'main_stem_fraction': main_stem / n_rivers if n_rivers else 0.0,
        'predicted_speedup': speedup,
    }

    print(f'Partitioning analysis for threads={threads}')
    print(f'  rivers in:              {summary["n_rivers"]:,}')
    print(
        f'  concurrent regions:     {summary["n_regions"]:,} '
        f'(largest {summary["largest_region"]:,}, smallest {summary["smallest_region"]:,})'
    )
    print(f'  parallel span:          {summary["parallel_span"]:,} rivers on the busiest thread')
    print(f'  sequential main stem:   {summary["main_stem"]:,} ({summary["main_stem_fraction"]:.1%})')
    print(f'  predicted speedup:      {summary["predicted_speedup"]:.2f}x (work only; ignores memory bandwidth)')
    return summary


def regions_to_layout(region: np.ndarray, downstream_index: np.ndarray) -> dict:
    """
    Turn a per-river region assignment into the index ranges a region-parallel routing kernel consumes.

    Nothing is reordered. In DFS computation order each region is already a contiguous run of rows, so this only
    locates the run boundaries. The main stem is what is left between the regions, which is generally SEVERAL
    ranges rather than one: a main stem river sits between the tributary subtrees that feed it. The kernel
    therefore takes a list of blocks and sweeps them in increasing index order, which is still a valid
    topological sweep because the whole table is topologically sorted.

    Args:
        region: per-river region id from assign_regions, -1 for the sequential main stem
        downstream_index: positional index of each river's downstream, -1 where there is none

    Returns a dict of arrays:
        region_starts / region_stops -- int32 (n_regions,), region r is [start, stop)
        region_outlet -- int32 (n_regions,), the river whose push is buffered; always stop - 1 in DFS order
        cut_target    -- int32 (n_regions,), index each region's outlet drains into, -1 at a basin outlet
        stem_starts / stem_stops     -- int32 (m,), the main stem blocks, in increasing index order
        n_regions     -- int, number of concurrent regions

    Raises:
        ValueError: if a region is not one contiguous run, or if a region's outlet drains into another region
            rather than the main stem, either of which would make the regions unsafe to route concurrently
    """
    n = region.shape[0]
    if downstream_index.shape[0] != n:
        raise ValueError(f'region has {n} values but downstream_index has {downstream_index.shape[0]}')
    concurrent = region[region >= 0]
    n_regions = int(concurrent.max()) + 1 if concurrent.size else 0
    if concurrent.size and not np.array_equal(np.unique(concurrent), np.arange(n_regions)):
        raise ValueError(f'region ids must be a contiguous range 0..{n_regions - 1} (plus -1 for the main stem)')

    # boundaries of every run of equal region id, then split them into concurrent regions and main stem blocks
    edges = np.nonzero(np.diff(region))[0] + 1
    starts = np.concatenate(([0], edges))
    stops = np.concatenate((edges, [n]))
    labels = region[starts]

    region_starts = np.full(n_regions, -1, dtype=np.int64)
    region_stops = np.full(n_regions, -1, dtype=np.int64)
    for start, stop, label in zip(starts[labels >= 0], stops[labels >= 0], labels[labels >= 0], strict=True):
        if region_starts[label] >= 0:
            raise ValueError(f'region {int(label)} is split across more than one range of rows; it must be one block')
        region_starts[label] = start
        region_stops[label] = stop

    # a region's outlet is the last river of its block, and its downstream must land in the main stem
    region_outlet = region_stops - 1
    cut_target = downstream_index[region_outlet] if n_regions else np.zeros(0, dtype=np.int64)
    if n_regions and np.any(region[cut_target[cut_target >= 0]] >= 0):
        raise ValueError(
            'a region outlet drains into another region rather than the sequential main stem; '
            'the partition is not safe to route concurrently (build it with assign_regions)'
        )

    return {
        'region_starts': region_starts.astype(np.int32),
        'region_stops': region_stops.astype(np.int32),
        'region_outlet': region_outlet.astype(np.int32),
        'cut_target': np.asarray(cut_target).astype(np.int32),
        'stem_starts': starts[labels < 0].astype(np.int32),
        'stem_stops': stops[labels < 0].astype(np.int32),
        'n_regions': n_regions,
    }


def subset_configs_to_river(
    target_river: int,
    params: PathInput,
    out_params: PathInput,
    weights: PathInput | None = None,
    out_weights: PathInput | None = None,
) -> None:
    """
    Subset routing parameters and weight tables to the target river and all rivers upstream of it.

    The target river becomes the outlet (next_river_id set to -1) in the subset. If weight
    table paths are provided, the weight table is also filtered to only include matching rivers.

    Args:
        target_river: river_id of the river to subset to (becomes the outlet)
        params: path to the full routing parameters parquet file
        out_params: path to write the subsetted parameters parquet file
        weights: path to the full grid weights netCDF file (optional)
        out_weights: path to write the subsetted grid weights netCDF file (optional, required if weights is given)
    """
    pdf = pd.read_parquet(params)
    if target_river not in set(pdf['river_id'].values.tolist()):
        raise ValueError(f'river_id {target_river} is not in the parameter table: {params}')

    graph = connectivity_to_digraph(pdf['river_id'].values, pdf['next_river_id'].values)
    upstreams = list(nx.ancestors(graph, target_river))
    upstreams.append(target_river)

    subset = pdf[pdf['river_id'].isin(upstreams)].copy()
    subset.loc[subset['river_id'] == target_river, 'next_river_id'] = -1
    subset.to_parquet(out_params)

    if weights is not None and out_weights is not None:
        with xr.open_dataset(weights) as ds:
            mask = np.isin(ds['river_id'].values, list(upstreams))
            ds.isel(index=mask).to_netcdf(out_weights)
    return


def connectivity_to_digraph(river_ids: np.ndarray, downstream_ids: np.ndarray) -> nx.DiGraph:
    """
    Build a NetworkX DiGraph from river connectivity arrays.

    Each edge goes from a river to its downstream river (including the -1 sentinel for outlets).

    Args:
        river_ids: 1D array of river ID integers
        downstream_ids: 1D array of downstream river ID integers (-1 for outlets)

    Returns:
        Directed graph with edges from each river_id to its next_river_id
    """
    graph = nx.DiGraph()
    graph.add_edges_from(zip(river_ids, downstream_ids, strict=True))
    return graph


def adjacency_matrix(river_ids: np.ndarray, downstream_ids: np.ndarray) -> scipy.sparse.csc_matrix:
    """
    Build a sparse adjacency matrix for the river network.

    Entry A[downstream_idx, upstream_idx] = 1 for each river that flows into a downstream river.
    Outlet rivers (downstream_id < 0) have no outgoing edges. The input arrays must be topologically
    sorted (upstream before downstream) — a ValueError is raised otherwise.

    Args:
        river_ids: 1D array of river ID integers, topologically sorted upstream to downstream
        downstream_ids: 1D array of downstream river ID integers (-1 for outlets)

    Returns:
        Sparse CSC matrix of shape (n, n) where n = len(river_ids)

    Raises:
        ValueError: if a downstream_id is not found in river_ids, or if the arrays are not
            topologically sorted
    """
    river_index = {int(river_id): idx for idx, river_id in enumerate(river_ids.tolist())}
    row_indices: list[int] = []
    col_indices: list[int] = []
    for upstream_idx, next_river_id in enumerate(downstream_ids.tolist()):
        if next_river_id < 0:
            continue
        if next_river_id not in river_index:
            raise ValueError(f'Unknown next_river_id: {next_river_id}')
        downstream_idx = river_index[int(next_river_id)]
        if downstream_idx <= upstream_idx:
            raise ValueError('params_file must be topologically sorted upstream to downstream')
        row_indices.append(downstream_idx)
        col_indices.append(upstream_idx)

    data = np.ones(len(row_indices), dtype=np.float32)
    return scipy.sparse.csc_matrix((data, (row_indices, col_indices)), shape=(river_ids.shape[0], river_ids.shape[0]))
