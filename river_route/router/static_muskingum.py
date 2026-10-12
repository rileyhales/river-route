"""
Muskingum routing with static coefficients: k, x, and the step are fixed for the whole run for each river, and
conditioning routes every river validly by splitting it into sub-reaches or taking its step in subcycles. This module
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
    each runoff step and takes subcycles of its own within it, settled in two jobs. The first sets each river's step:
    with dt_routing auto (None), its own valid step, the runoff step over the fewest subcycles that keep its
    coefficients non-negative, and with dt_routing given, exactly dt_routing. The second, conditioning, splits a river
    too long for its step into the sub-reaches of Network.sub_reach_k_and_x, each taking its share of the river's
    runoff, and routes a river too short for it in subcycles of its step. Without conditioning every river is routed
    whole at its step, and a river with a negative coefficient at it is reported.

    The subcycles of each sub-reach are combined into one step by _combine_subcycles, which gives the discharge of
    taking them one at a time, so every sub-reach takes one step per runoff step and hands its discharge to the next at
    the end of it, as a river hands its discharge downstream.

    Raises:
        ValueError: if the coefficients do not sum to 1, or invalid_coefficients is raise and a coefficient is negative
    """
    if dt_routing:
        subcycles = np.full(network.size, dt_runoff // dt_routing, dtype=np.int64)
    else:
        subcycles, _ = network.subcycles(dt_runoff)
    sub_reaches = np.ones(network.size, dtype=np.int64)
    if configs.conditioning:
        sub_reaches, conditioning_subcycles, _ = network.conditioning(dt_runoff / subcycles)
        subcycles = subcycles * conditioning_subcycles
    network.check_validity(dt_runoff, action=configs.invalid_coefficients, sub_reaches=sub_reaches, subcycles=subcycles)
    k_reach, x_reach = network.sub_reach_k_and_x(sub_reaches)
    share = k_reach / np.repeat(network.k.astype(np.float64), sub_reaches)  # of the river's k, and of its runoff
    dt_div_k = np.repeat(dt_runoff / subcycles, sub_reaches) / k_reach  # each sub-reach's own step over its k
    denominator = dt_div_k + 2 * (1 - x_reach)
    c1 = (dt_div_k - 2 * x_reach) / denominator
    c2 = (dt_div_k + 2 * x_reach) / denominator
    c3 = (2 * (1 - x_reach) - dt_div_k) / denominator
    c4 = (c1 + c2) * share / dt_runoff
    c1, c2, c3, c4 = _combine_subcycles(c1, c2, c3, c4, np.repeat(subcycles, sub_reaches))
    if not np.allclose(c1 + c2 + c3, 1):
        raise ValueError('Muskingum coefficients do not sum to 1, check k and x of the network file and dt_routing')
    sub_reach_starts = np.zeros(network.size + 1, dtype=np.int32)
    sub_reach_starts[1:] = np.cumsum(sub_reaches)
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
    """
    combined_c1 = np.zeros_like(c1)
    combined_c2 = np.zeros_like(c2)
    combined_c4 = np.zeros_like(c4)
    weight = np.ones_like(c3)  # c3^(s - j) of subcycle j, from the last subcycle back to the first
    for subcycles_after in range(int(subcycles.max())):
        subcycle = subcycles - subcycles_after  # j, from 1
        taken = subcycle >= 1
        end = subcycle / subcycles  # the fraction of the step at the end of subcycle j
        start = (subcycle - 1) / subcycles
        combined_c1 += np.where(taken, weight * (c1 * end + c2 * start), 0.0)
        combined_c2 += np.where(taken, weight * (c1 * (1 - end) + c2 * (1 - start)), 0.0)
        combined_c4 += np.where(taken, weight * c4, 0.0)
        weight = np.where(taken, weight * c3, weight)
    return combined_c1, combined_c2, weight, combined_c4
