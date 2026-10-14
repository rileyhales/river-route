"""
The fewest sub-reaches a river splits into when equal sub-reaches cannot make every coefficient of it non-negative at a
step. The sub-reaches may differ in travel time, each takes subcycles of the step of its own, and every one keeps the
river's diffusion: a sub-reach of travel time k_i takes x_i = 1/2 - (1/2 - x) k / k_i, which keeps the river's travel
time variance k**2 (1 - 2x) across them whatever their travel times. The window of steps at which a sub-reach is valid,
2 k_i x_i <= step <= 2 k_i (1 - x_i), is then |step - k_i| <= (1 - 2x) k, as wide for every sub-reach of the river.

In m subcycles of a step dt, a sub-reach may take any travel time within (1 - 2x) k of dt/m and no shorter than
(1 - 2x) k, where its x is 0, so the travel times a sub-reach may take are a union of windows about dt, dt/2, dt/3, and
so on. A river splits into n sub-reaches when its k is a sum of n of them. The totals of 1, 2, 3, ... sub-reaches are
grown as sorted intervals, and n sub-reaches make k when a total of the first n - n//2 of them and a total of the last
n//2 do, so the totals are grown only to half the fewest n. The window of each sub-reach is then read back from them.
"""

import numba
import numpy as np

from ..types import Float64Array, IntArray

__all__ = ['split_rivers', 'fewest_sub_reaches']

# the share of a river's diffusion that its windows are held inside their edges by, so that rounding, as in a network
# file stored in float32, cannot flip the sign of a coefficient
MARGIN = 1e-3
# the shortest sub-reach, as a share of the shorter of the river's travel time and the step. Near x = 1/2 the windows
# narrow, and a split can need a sub-reach a tiny share of the river long, taking tens of thousands of subcycles only
# to make the travel times sum to k. A river whose every split needs one is kept whole.
SHORTEST_SHARE = 1 / 64
# the most totals one more sub-reach is added to at once. A river that needs more is kept whole.
MAX_PAIRS = 1 << 22


@numba.njit(cache=True, nogil=True)
def _windows(k: float, diffusion: float, dt: float) -> tuple[Float64Array, Float64Array]:
    """
    The travel times a sub-reach of a river may take in subcycles of dt, as windows from the longest to the shortest:
    one for each count of subcycles m from the first whose window reaches down to k, and last the window that every
    count joins into from the first whose window overlaps the next one's. None is shorter than the diffusion, where a
    sub-reach's x is 0, or than SHORTEST_SHARE of the river or the step.
    """
    half = diffusion * (1 - MARGIN)
    shortest = max(diffusion, SHORTEST_SHARE * min(k, dt))
    first = max(1, int(np.ceil(dt / (k + half))))
    # the last count whose window reaches up to the shortest, or past which every window lies within its window
    last = max(first, int(np.floor(dt / max(shortest - half, half))))
    # windows m and m + 1 overlap when dt / (m (m + 1)) <= 2 half, and so does every pair of windows after them
    joined = min(max(first, int(np.ceil((np.sqrt(1 + 2 * dt / half) - 1) / 2))), last)
    lower = np.empty(joined - first + 1)
    upper = np.empty(joined - first + 1)
    for window in range(joined - first):
        lower[window] = max(shortest, dt / (first + window) - half)
        upper[window] = dt / (first + window) + half
    lower[-1] = max(shortest, dt / last - half)  # the windows of every count from joined to last, joined into one
    upper[-1] = dt / joined + half
    return lower, upper


@numba.njit(cache=True, nogil=True)
def _merge(lower: Float64Array, upper: Float64Array) -> tuple[Float64Array, Float64Array]:
    """The union of the intervals from lower to upper, as sorted, disjoint intervals."""
    merged_lower = np.empty(lower.size)
    merged_upper = np.empty(lower.size)
    count = 0
    for interval in np.argsort(lower):
        if count > 0 and lower[interval] <= merged_upper[count - 1]:
            merged_upper[count - 1] = max(merged_upper[count - 1], upper[interval])
        else:
            merged_lower[count] = lower[interval]
            merged_upper[count] = upper[interval]
            count += 1
    return merged_lower[:count].copy(), merged_upper[:count].copy()


@numba.njit(cache=True, nogil=True)
def _holds(lower: Float64Array, upper: Float64Array, total: float) -> bool:
    """Whether the sorted, disjoint intervals from lower to upper hold total."""
    interval = np.searchsorted(lower, total, side='right') - 1
    return bool(interval >= 0 and upper[interval] >= total)


@numba.njit(cache=True, nogil=True)
def _add(
    total_lower: Float64Array, total_upper: Float64Array, lower: Float64Array, upper: Float64Array, k: float
) -> tuple[Float64Array, Float64Array]:
    """The totals, up to k, of one more sub-reach in any of the windows on any of the totals given."""
    sum_lower = np.empty(total_lower.size * lower.size)
    sum_upper = np.empty(total_lower.size * lower.size)
    count = 0
    for total in range(total_lower.size):
        for window in range(lower.size):
            if total_lower[total] + lower[window] <= k:
                sum_lower[count] = total_lower[total] + lower[window]
                sum_upper[count] = min(total_upper[total] + upper[window], k)
                count += 1
    return _merge(sum_lower[:count], sum_upper[:count])


@numba.njit(cache=True, nogil=True)
def _last_window(
    total_lower: Float64Array, total_upper: Float64Array, lower: Float64Array, upper: Float64Array, total: float
) -> tuple[int, float]:
    """
    The longest window whose sub-reach, added to one of the totals given, makes total, and the total it is added to,
    tested as _add builds the totals so that rounding cannot lose the one that made total. -1 if none does.
    """
    for window in range(lower.size):
        for interval in range(total_lower.size):
            if total_lower[interval] + lower[window] <= total <= total_upper[interval] + upper[window]:
                rest = min(total_upper[interval], total - lower[window])
                return window, max(total_lower[interval], rest)
    return -1, total


@numba.njit(cache=True, nogil=True)
def _halves(levels: list, count: int, k: float) -> tuple[bool, float, float]:
    """
    Whether count sub-reaches make k, as a total of the first count - count//2 of them and a total of the last
    count//2, and those two totals. levels[n - 1] are the totals of n sub-reaches, up to k.
    """
    first_lower, first_upper = levels[count - count // 2 - 1]
    if count == 1:
        return _holds(first_lower, first_upper, k), k, 0.0
    last_lower, last_upper = levels[count // 2 - 1]
    for interval in range(first_lower.size):
        # the last totals' interval reaching furthest down to k minus this one, and the one after it, against rounding
        match = np.searchsorted(last_lower, k - first_lower[interval], side='right') - 1
        for other in range(max(match, 0), min(match + 2, last_lower.size)):
            if first_lower[interval] + last_lower[other] <= k <= first_upper[interval] + last_upper[other]:
                first = max(first_lower[interval], min(first_upper[interval], k - last_lower[other]))
                return True, first, max(last_lower[other], min(last_upper[other], k - first))
    return False, k, 0.0


@numba.njit(cache=True, nogil=True)
def _read_back(levels: list, count: int, total: float, lower: Float64Array, upper: Float64Array, k: float) -> IntArray:
    """The windows of count sub-reaches whose totals of every count are levels, which make total: from the last
    sub-reach to the first, the longest window that leaves a total the others make."""
    chosen = np.empty(count, dtype=np.int64)
    if count == 0:
        return chosen
    for level in range(count - 1, 0, -1):
        total_lower, total_upper = levels[level - 1]
        chosen[level], total = _last_window(total_lower, total_upper, lower, upper, total)
    chosen[0] = -1
    for window in range(lower.size):
        if lower[window] <= total <= min(upper[window], k):
            chosen[0] = window
            break
    return chosen


@numba.njit(cache=True, nogil=True)
def _travel_times(chosen: IntArray, lower: Float64Array, upper: Float64Array, k: float) -> Float64Array:
    """
    The travel times of sub-reaches in the windows chosen, longest first, that sum to k: every one the same share of
    the way across its window, so that none is nearer an edge of its window than it must be.
    """
    shortest = lower[chosen]
    longest = upper[chosen]
    span = longest.sum() - shortest.sum()
    across = 0.5 if span <= 0 else min(1.0, max(0.0, (k - shortest.sum()) / span))
    times = shortest + across * (longest - shortest)
    times *= k / times.sum()
    return np.sort(times)[::-1].copy()


@numba.njit(cache=True, nogil=True)
def fewest_sub_reaches(k: float, x: float, dt: float, most: int) -> tuple[Float64Array, IntArray]:
    """
    The travel times, upstream to downstream, of the fewest sub-reaches of a river, at most ``most``, whose
    coefficients are all non-negative at dt, and the subcycles of dt each takes: the fewest whose step is in its
    window, as Network.subcycles gives a river. Empty arrays where no split of at most ``most`` sub-reaches is valid,
    as at x = 1/2, whose windows close to single steps.
    """
    none = (np.empty(0), np.empty(0, dtype=np.int64))
    diffusion = (1 - 2 * x) * k
    if not (k > 0 and dt > 0 and diffusion > 0) or np.ceil(k / (dt + diffusion)) > most:
        return none
    lower, upper = _windows(k, diffusion, dt)
    reachable = lower <= k
    levels = [_merge(lower[reachable], np.minimum(upper[reachable], k))]  # the totals of 1, 2, ... sub-reaches
    for count in range(1, most + 1):
        if count - count // 2 > len(levels):
            total_lower, total_upper = levels[-1]
            if total_lower.size * lower.size > MAX_PAIRS:
                return none
            levels.append(_add(total_lower, total_upper, lower, upper, k))
        found, first, last = _halves(levels, count, k)
        if not found:
            continue
        chosen = np.concatenate(
            (
                _read_back(levels, count - count // 2, first, lower, upper, k),
                _read_back(levels, count // 2, last, lower, upper, k),
            )
        )
        if np.any(chosen < 0):
            return none
        times = _travel_times(chosen, lower, upper, k)
        return times, np.ceil(dt / (times + diffusion)).astype(np.int64)
    return none


@numba.njit(cache=True, nogil=True)
def _grown(values: np.ndarray, size: int) -> np.ndarray:
    """values in an array of at least size, at least twice as long."""
    grown = np.empty(max(size, 2 * values.size), dtype=values.dtype)
    grown[: values.size] = values
    return grown


@numba.njit(cache=True, nogil=True)
def split_rivers(
    k: Float64Array, x: Float64Array, dt: Float64Array, most: int
) -> tuple[IntArray, Float64Array, IntArray]:
    """
    The fewest_sub_reaches of each river.

    Returns:
        tuple: (sub_reaches, travel_times, subcycles). ``sub_reaches`` is the count of each river, 0 where none is
            valid, and ``travel_times`` and ``subcycles`` have a value for each sub-reach, river after river.
    """
    sub_reaches = np.zeros(k.size, dtype=np.int64)
    travel_times = np.empty(4 * k.size)
    subcycles = np.empty(4 * k.size, dtype=np.int64)
    filled = 0
    for river in range(k.size):
        river_times, river_subcycles = fewest_sub_reaches(k[river], x[river], dt[river], most)
        count = river_times.size
        if filled + count > travel_times.size:
            travel_times = _grown(travel_times, filled + count)
            subcycles = _grown(subcycles, filled + count)
        travel_times[filled : filled + count] = river_times
        subcycles[filled : filled + count] = river_subcycles
        sub_reaches[river] = count
        filled += count
    return sub_reaches, travel_times[:filled].copy(), subcycles[:filled].copy()
