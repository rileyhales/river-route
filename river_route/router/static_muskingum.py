"""
Muskingum routing with static coefficients: k, x, and the step are fixed for the whole run for each river, and every
river is routed validly by splitting it into sub-reaches, taking its step in subcycles, or both. This module
settles each river's sub-reaches and the coefficients of each, in float64, and either kernel of river_route.kernels
routes them.
"""

from typing import TYPE_CHECKING, NamedTuple

import numpy as np

from ..types import Float64Array, FloatArray, Int32Array, IntArray

if TYPE_CHECKING:
    from ..configs.Configs import Configs
    from ..network.Network import Network

__all__ = ['Coefficients', 'prepare_routing']


class Coefficients(NamedTuple):
    """
    The sub-reaches each river is routed as, one for a river that is not split, and the Muskingum coefficients of each,
    fixed for the whole run. Each sub-reach's coefficients are of a whole runoff step, its subcycles combined into one
    step, so every sub-reach takes one step per runoff step.
    """

    sub_reach_starts: Int32Array  # (n_rivers + 1,) channel_state index of each river's first sub-reach, then the total
    c1: FloatArray  # (sub-reaches,)
    c2: FloatArray  # (sub-reaches,)
    c3: FloatArray  # (sub-reaches,)
    c4: FloatArray  # (sub-reaches,) the coefficient of the river's runoff volume (m³) of a runoff step


def prepare_routing(network: Network, configs: Configs, dt_runoff: int, dt_routing: int | None) -> Coefficients:
    """
    Each river's sub-reaches and the coefficients of each. Every river hands its discharge downstream at the end of
    each runoff step and computes within it as Network.conditioning lays it out at its step: dt_routing when given,
    and with dt_routing auto (None) the runoff step, at which a river too short for it takes the fewest subcycles that
    keep its coefficients non-negative, its own valid step. With use_subreaches a river too long for its step is split
    into sub-reaches, each taking its share of the river's runoff, with use_subcycles a river too short for it is
    routed in subcycles of it, and with both a river neither makes valid alone is split into sub-reaches that each take
    subcycles of their own. Without them a river is routed whole at its step, and a river with a negative coefficient
    at it is reported. Network.check_routable refuses sub-reaches the kernels cannot route before any coefficient is
    built.

    The subcycles of each sub-reach are combined into one step by _combine_subcycles, which gives the discharge of
    taking them one at a time, so every sub-reach takes one step per runoff step and hands its discharge to the next at
    the end of it, as a river hands its discharge downstream.

    Raises:
        ValueError: if the coefficients do not sum to 1, invalid_coefficients is raise and a coefficient is negative,
            or use_subreaches or use_subcycles is on and a river's x is above MAX_STABILIZED_X
    """
    step = dt_routing or dt_runoff
    if configs.use_subreaches or configs.use_subcycles:
        network.check_stabilizable()
    layout = network.conditioning(step, use_subreaches=configs.use_subreaches, use_subcycles=configs.use_subcycles)
    network.check_routable(layout.sub_reaches)
    network.check_validity(step, action=configs.invalid_coefficients, layout=layout)
    share = layout.k / np.repeat(network.k, layout.sub_reaches)  # of the river's k and runoff
    subcycles = dt_runoff // step * layout.subcycles  # the steps of its own each sub-reach takes in each runoff step
    dt_div_k = dt_runoff / subcycles / layout.k  # each sub-reach's own step over its k
    denominator = dt_div_k + 2 * (1 - layout.x)
    c1 = (dt_div_k - 2 * layout.x) / denominator
    c2 = (dt_div_k + 2 * layout.x) / denominator
    c3 = (2 * (1 - layout.x) - dt_div_k) / denominator
    c4 = (c1 + c2) * share / dt_runoff
    c1, c2, c3, c4 = _combine_subcycles(c1, c2, c3, c4, subcycles)
    if not np.allclose(c1 + c2 + c3, 1):
        raise ValueError('Muskingum coefficients do not sum to 1, check k and x of the network file and dt_routing')
    sub_reach_starts = np.zeros(network.size + 1, dtype=np.int32)
    sub_reach_starts[1:] = np.cumsum(layout.sub_reaches)
    return Coefficients(sub_reach_starts, *(c.astype(np.float32) for c in (c1, c2, c3, c4)))


def _combine_subcycles(
    c1: Float64Array, c2: Float64Array, c3: Float64Array, c4: Float64Array, subcycles: IntArray
) -> tuple[Float64Array, Float64Array, Float64Array, Float64Array]:
    """
    The coefficients of one step that is s subcycles of a reach, its inflow rising linearly from I_t to I_{t+1} across
    them and its runoff spread evenly over them. Subcycle j of s ends at the fraction j/s of the step, with
    Q_j = c3 Q_(j-1) + c1 I(j/s) + c2 I((j-1)/s) + c4 V, so Q_s = c3^s Q_0 plus the terms of each subcycle j times
    c3^(s-j), whose inflows split between I_t and I_{t+1}: Q_(t+1) = c3' Q_t + c1' I_{t+1} + c2' I_t + c4' V with
    c3' = c3^s. Like each subcycle's, c1' + c2' + c3' is 1, and none is negative where none of the subcycle's is.

    With G the sum of c3^i and A the sum of i c3^i over i from 0 to s - 1, c1' = c1 (G - A/s) + c2 (G - (A + G)/s),
    c2' = (c1 A + c2 (A + G))/s, and c4' = c4 G. G, A, and c3^s are built over the bits of s, doubling the count of
    terms at each bit and adding one more where the bit is set, so a reach of s subcycles costs log2(s) steps.
    """
    power = np.ones_like(c3)  # c3^n, of the n terms summed so far
    geometric = np.zeros_like(c3)  # G of n terms
    weighted = np.zeros_like(c3)  # A of n terms
    terms = np.zeros_like(subcycles)
    for bit in range(int(subcycles.max()).bit_length() - 1, -1, -1):
        # the second n terms are the first n times c3^n, each i of them n more
        weighted = weighted + power * (weighted + terms * geometric)
        geometric = geometric + power * geometric
        power = power * power
        terms = 2 * terms
        taken = (subcycles >> bit) & 1 == 1
        weighted = np.where(taken, weighted + terms * power, weighted)
        geometric = np.where(taken, geometric + power, geometric)
        power = np.where(taken, power * c3, power)
        terms = terms + taken
    combined_c1 = c1 * (geometric - weighted / subcycles) + c2 * (geometric - (weighted + geometric) / subcycles)
    combined_c2 = (c1 * weighted + c2 * (weighted + geometric)) / subcycles
    return combined_c1, combined_c2, power, c4 * geometric
