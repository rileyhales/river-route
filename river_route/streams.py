"""
Static analysis of a routing parameter table (river_id, downstream_river_id, k, x).

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

eventually this module should merge with the present tools.py
"""
from typing import Tuple

import numpy as np
import pandas as pd


def river_connectivity_is_valid(df: pd.DataFrame) -> bool:
    """
    checks that river ids are unique, all downstreams exist as rivers except -1, and that the table is
    topologically sorted from upstream to downstream (which also rules out cycles).
    """
    total_rows = df.shape[0]
    unique_rivers = df['river_id'].nunique()
    if unique_rivers != total_rows:
        print(f'river_id column must be unique: {total_rows} rows, {unique_rivers} unique river ids')
        return False

    downstreams = set(df['downstream_river_id'])
    river_ids = set(df['river_id'])
    downstreams_not_in_rivers = downstreams - river_ids
    if downstreams_not_in_rivers != {-1}:  # only -1 doesn't need to be in river_ids
        print(
            f'all downstream_river_id values must be in river_id column or -1. This might have been intentional or '
            f'may indicate a hole in the topology. Downstream ids not in river ids: {downstreams_not_in_rivers}'
        )
        return False

    # check that the rivers are topologically sorted from upstream to downstream
    river_id_to_index = {river_id: idx for idx, river_id in enumerate(df['river_id'])}
    river_id_to_downstream_id = dict(zip(df['river_id'], df['downstream_river_id']))
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


def required_subreaches(k: np.ndarray, x: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
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
    small = [d for d in range(1, int(n ** 0.5) + 1) if n % d == 0]
    return np.array(sorted(set(small + [n // d for d in small])), dtype=np.int64)


def assign_stable_dt(k: np.ndarray, x: np.ndarray, period: int = 3600) -> Tuple[np.ndarray, np.ndarray]:
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
        if resolvable.any() else pd.Series(dtype=np.int64)
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
    print(f'  resolvable:          {summary["n_resolvable"]:,} '
          f'(dt range {summary["min_dt"]}..{summary["max_dt"]} s)')
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

    INF = np.iinfo(np.int64).max
    best_cost = np.full(k.shape, INF, dtype=np.int64)
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
        df: parameter table with columns river_id, downstream_river_id, k, x (extra columns are preserved)
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
        df: parameter table with river_id, downstream_river_id, k, x. If subdivisions/substeps columns are not
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
    orig_down = df['downstream_river_id'].to_numpy(dtype=np.int64)
    k = df['k'].to_numpy(dtype=np.float64)
    x = df['x'].to_numpy(dtype=np.float64)
    stable = (df['stable'].to_numpy() if 'stable' in df.columns
              else np.ones(orig_id.shape[0], dtype=bool))
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
    if np.isnan(mapped).any():  # a downstream_river_id absent from river_id is a topology hole, not a valid -1 outlet
        raise ValueError('downstream_river_id values reference ids not present in river_id (topology hole); '
                         'validate with river_connectivity_is_valid before expanding')
    down_orig_index[has_down] = mapped.astype(np.int64)
    outlet_downstream = np.where(down_orig_index >= 0,
                                 group_start[np.clip(down_orig_index, 0, n_orig - 1)], -1)
    downstream_index[is_outlet] = outlet_downstream  # outlets are emitted in original order

    # the kernel sweeps in array order and requires a topologically sorted DAG (each reach feeds a later one or -1);
    # an unsorted input would silently push flow into an already-finalized reach and lose mass
    if not np.all((downstream_index < 0) | (downstream_index > np.arange(total))):
        raise ValueError('input rivers must be topologically sorted upstream-before-downstream '
                         '(downstream_river_id must appear after river_id); run river_connectivity_is_valid')

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
    print(f'  fixable by splitting:   {summary["n_resolved_by_split"]:,} '
          f'(up to {summary["max_subreaches"]} sub-reaches each)')
    print(f'  unresolvable errors:    {summary["n_unresolvable"]:,} '
          f'(too short for dt or window empty; kept as 1 reach)')
    print(f'  total streams in fixed network: {summary["total_streams"]:,}')
    print(f'  new streams to create:          {summary["streams_to_create"]:,}')
    return summary