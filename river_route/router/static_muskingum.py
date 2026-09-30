"""
Muskingum routing with static coefficients: each river's c1, c2, and c3 are fixed for the whole run by its k and x and
the routing step. On a stabilized network a river too long for the routing step is routed in substeps, equal
sub-reaches in series, and a river too short for it is routed in subcycles, shorter steps of its own, each as
``Network.conditioning`` finds it needs to be stable.
"""

from typing import TYPE_CHECKING, Literal, NamedTuple

import numba
import numpy as np
from numba.extending import overload

from ._numba_kernels import STANDARD_LAYOUT, Layout, add_river_forcing, is_argument_type, route_river

if TYPE_CHECKING:
    from ..network.Network import Network
    from ..types import FloatArray, IntArray

__all__ = ['StaticMuskingum', 'prepare_routing']


class StaticMuskingum(NamedTuple):
    """Muskingum coefficients per river, fixed for the whole run."""

    c1: np.ndarray  # (n,) float32
    c2: np.ndarray  # (n,) float32
    c3: np.ndarray  # (n,) float32
    c4_dt: np.ndarray  # (n,) float32, c4 / dt_runoff per sub-reach: turns a catchment runoff volume into an inflow rate


def prepare_routing(
    network: Network,
    network_type: Literal['standard', 'stabilized'],
    dt_routing: int,
    dt_runoff: int,
    unstable_coefficients: Literal['warn', 'raise', 'ignore'],
) -> tuple[Layout, StaticMuskingum]:
    """
    Each river's layout and coefficients at ``dt_routing``. On a stabilized network a river too long for dt_routing is
    divided into N substeps, sub-reaches that are each k/N long and take 1/N of its runoff, so one set of coefficients
    serves all of them, and a river too short for it is routed in m subcycles of dt_routing/m, so its coefficients are
    built for that step. Rivers that are still unstable are reported as ``unstable_coefficients`` says.

    Args:
        network: the network whose k and x the coefficients are built from
        network_type: ``standard`` routes each river whole, ``stabilized`` in substeps and subcycles
        dt_routing: the routing step in seconds
        dt_runoff: the runoff step in seconds, which turns each catchment runoff volume into an inflow rate
        unstable_coefficients: what to do about rivers that are not Muskingum-stable: warn, raise, or ignore

    Returns:
        tuple: (the Layout each river is routed with, the StaticMuskingum coefficients of every river)

    Raises:
        ValueError: if the coefficients do not sum to 1, or unstable_coefficients is raise and a river is unstable
    """
    substeps: IntArray | None = None
    subcycles: IntArray | None = None
    dt_river: FloatArray | int = dt_routing
    layout, k = STANDARD_LAYOUT, network.k
    if network_type == 'stabilized':
        substeps, subcycles, _ = network.conditioning(dt_routing)
        reach_indptr = np.zeros(network.size + 1, dtype=np.int64)
        np.cumsum(substeps, out=reach_indptr[1:])
        layout = Layout(reach_indptr, subcycles if np.any(subcycles > 1) else STANDARD_LAYOUT.subcycles)
        k = (network.k / substeps).astype(np.float32)
        dt_river = (dt_routing / subcycles).astype(np.float32)  # float32 like the standard network's math
    dt_div_k = dt_river / k
    denominator = dt_div_k + (2 * (1 - network.x))
    two_x = 2 * network.x
    # contiguous arrays iterate faster in kernels due to cpu and ram access patterns
    c1 = np.ascontiguousarray((dt_div_k - two_x) / denominator, dtype=np.float32)
    c2 = np.ascontiguousarray((dt_div_k + two_x) / denominator, dtype=np.float32)
    c3 = np.ascontiguousarray(((2 * (1 - network.x)) - dt_div_k) / denominator, dtype=np.float32)
    c4_dt = np.ascontiguousarray((c1 + c2) / dt_runoff, dtype=np.float32)
    if substeps is not None:
        c4_dt = np.ascontiguousarray(c4_dt / substeps, dtype=np.float32)
    network.check_stability(dt_routing, action=unstable_coefficients, substeps=substeps, subcycles=subcycles)
    if not np.allclose(c1 + c2 + c3, 1):
        raise ValueError('Muskingum coefficients do not sum to 1, check k and x of the network file and dt_routing')
    return layout, StaticMuskingum(c1, c2, c3, c4_dt)


@numba.njit(cache=True, nogil=True, fastmath={'contract'})
def _solve_muskingum_recurrence_in_place(q, c3, work, n):
    """
    Replace work[g] with q_g = c3 q_(g-1) + work[g] in place and return the last q.

    Steps are taken eight at a time: within a group the partial sums p_k = sum_j c3^(k-j) work[j] do not depend on q,
    so they overlap with the previous group, and only q_(g+7) = c3^8 q_(g-1) + p_7 waits on the one before it. The
    serial chain is one multiply-add per eight steps instead of one per step.
    """
    c3_2 = c3 * c3
    c3_3 = c3_2 * c3
    c3_4 = c3_2 * c3_2
    c3_5 = c3_4 * c3
    c3_6 = c3_4 * c3_2
    c3_7 = c3_4 * c3_3
    c3_8 = c3_4 * c3_4
    g = 0
    while g + 8 <= n:
        a0 = work[g]
        p1 = c3 * a0 + work[g + 1]
        p2 = c3 * p1 + work[g + 2]
        p3 = c3 * p2 + work[g + 3]
        p4 = c3 * p3 + work[g + 4]
        p5 = c3 * p4 + work[g + 5]
        p6 = c3 * p5 + work[g + 6]
        p7 = c3 * p6 + work[g + 7]
        work[g] = c3 * q + a0
        work[g + 1] = c3_2 * q + p1
        work[g + 2] = c3_3 * q + p2
        work[g + 3] = c3_4 * q + p3
        work[g + 4] = c3_5 * q + p4
        work[g + 5] = c3_6 * q + p5
        work[g + 6] = c3_7 * q + p6
        q = c3_8 * q + p7
        work[g + 7] = q
        g += 8
    while g < n:
        q = c3 * q + work[g]
        work[g] = q
        g += 1
    return q


@numba.njit(cache=True, nogil=True)
def _interpolate_upstream_inflow_for_subcycles(inflow, subcycles, c1, c2, work, n_routing):
    """
    The upstream forcing of a river routed in subcycles, steps of its own, per routing step, with the upstream series
    interpolated linearly between the routing step levels: work[g * subcycles + s] = c1 * after + c2 * before for each
    of its own steps.

    Compiled without fastmath contract, unlike the rest of the river routine. With contract, whether these multiply-
    adds were fused differed between freshly compiled code and the same code loaded from the numba cache, so a
    subcycled river's discharge changed in the last bit depending on whether the kernels had just been compiled.
    """
    inverse_subcycles = np.float32(1.0 / subcycles)
    for g in range(n_routing):
        start = inflow[g]
        rise = (inflow[g + 1] - start) * inverse_subcycles
        before = start
        for s in range(subcycles):
            after = start + rise * np.float32(s + 1)
            work[g * subcycles + s] = c1 * after + c2 * before
            before = after


def _route_river_with_static_coefficients(
    method, layout, q_t, r, catchment_runoff, inflow, downstream_inflow, discharge, work, chain
):
    """
    Route river r's whole series through its ``n_pieces`` equal sub-reaches in series, whose states are
    ``q_t[p0:p0 + n_pieces]``: one piece on a standard network, and as many as stability needs on a stabilized one.
    The last piece's unclamped series is added into ``downstream_inflow`` and its clamped per-step mean is written into
    ``discharge``; the states are updated in place. ``catchment_runoff`` is None for channel routing, and every piece
    takes an equal share of it through ``c4_dt``, added into the piece's forcing by add_river_forcing.

    ``subcycles`` overrides the routing step of the river: each routing step is taken as that many equal steps of its
    own, which is how a river too short for the routing step is kept stable. The upstream series arrives once per
    routing step and is interpolated linearly between those levels, and the downstream river is handed the level at the
    end of each routing step, since Muskingum inflow is the discharge at the time levels. With one subcycle the river
    is routed at the routing step.

    ``work`` holds every step the river takes and ``chain`` one more, which carries a piece's series to the next piece
    in the same layout as an inflow row.

    Only the recurrence q = c3 q + forcing is serial, so it runs alone in its own loop; the loops over the steps before
    and after it are independent per step and vectorize. Fusing multiply-adds is the only fastmath flag, which keeps NaN
    semantics and makes each serial step one instruction.
    """
    expanded = layout.reach_indptr.shape[0] > 0
    p0 = layout.reach_indptr[r] if expanded else r
    n_pieces = layout.reach_indptr[r + 1] - p0 if expanded else 1
    subcycles = layout.subcycles[r] if layout.subcycles.shape[0] > 0 else 1
    c1 = method.c1[r]
    c2 = method.c2[r]
    c3 = method.c3[r]
    c4_dt = method.c4_dt[r]
    zero = np.float32(0.0)
    n_steps = discharge.shape[0]
    n_routing = inflow.shape[0] - 1  # an inflow row holds the level before the first routing step and at every step
    n_per_step = (n_routing // n_steps) * subcycles  # the river's own steps per runoff step
    n_fine = n_steps * n_per_step
    piece_inflow = inflow
    for j in range(n_pieces):
        # forcing of every step, all known before the recurrence starts
        if j == 0 and subcycles > 1:
            _interpolate_upstream_inflow_for_subcycles(inflow, subcycles, c1, c2, work, n_routing)
        else:
            for h in range(n_fine):
                work[h] = c1 * piece_inflow[h + 1] + c2 * piece_inflow[h]
        if catchment_runoff is not None:
            add_river_forcing(catchment_runoff, c4_dt, n_steps, n_per_step, work)

        # every step of the piece's inflow has been read, so the chain row can take this piece's series in its place
        q = q_t[p0 + j]
        if j == n_pieces - 1:
            downstream_inflow[0] += q
            q_t[p0 + j] = _solve_muskingum_recurrence_in_place(q, c3, work, n_fine)
            if subcycles == 1:
                for g in range(n_routing):
                    downstream_inflow[g + 1] += work[g]
            else:
                for g in range(n_routing):
                    downstream_inflow[g + 1] += work[(g + 1) * subcycles - 1]
        else:
            chain[0] = q
            q_t[p0 + j] = _solve_muskingum_recurrence_in_place(q, c3, work, n_fine)
            for h in range(n_fine):
                chain[h + 1] = work[h]
            piece_inflow = chain

    # only coefficients that are not Muskingum-stable give negative discharge, which unstable_coefficients reports.
    # It is clamped to zero here, which does not conserve mass; the unclamped series is what flows downstream.
    if n_per_step == 1:
        for t in range(n_steps):
            discharge[t] = max(zero, work[t])
    else:
        inv_per_step = np.float32(1.0 / n_per_step)
        for t in range(n_steps):
            interval_sum = zero
            for h in range(t * n_per_step, (t + 1) * n_per_step):
                interval_sum += work[h]
            discharge[t] = max(zero, interval_sum * inv_per_step)
    return


@overload(route_river, jit_options={'nogil': True, 'fastmath': {'contract'}})
def _route_river_with_static_muskingum(
    method, layout, q_t, r, catchment_runoff, inflow, downstream_inflow, discharge, work, chain
):
    """
    Route rivers with StaticMuskingum parameters by the static coefficients routine. The routine is itself the
    implementation numba compiles into route_job, rather than a lambda calling a separately jitted function, which
    added a call per river that passes a dozen arrays and measured 1% slower on the Amazon.
    """
    if is_argument_type(method, StaticMuskingum):
        return _route_river_with_static_coefficients
