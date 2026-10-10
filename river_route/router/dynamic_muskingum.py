"""
Muskingum routing with dynamic coefficients, the variable parameter Muskingum method: regression equations give each
river's k and x at its own discharge, and its coefficients are rebuilt from them at every runoff step. Regression
holds the dynamic columns of the network file, COLUMNS_FOR_PARAMETERS['dynamic'], and is an overload of the
k_and_x_at_discharges stage: k = dynamicKCoefficient * Q ** dynamicKExponent and
x = 1/2 - dynamicXCoefficient * Q ** dynamicXExponent, with Q held within dynamicMinQ and dynamicMaxQ and x held from 0
to 0.5, or to RESTRICTED_LARGEST_X with restrict_x_range.

Each step's k and x are evaluated at the discharge entering the river over the step, its inflow plus its runoff, which
is known for the river's whole series before it is routed, since every river upstream of it is routed first. So k and x
are evaluated for the whole series at once, which the compiler vectorizes, and the outflow is one multiply-add per
step. The Dynamic Parameters reference of the documentation describes the method, its accuracy, and its cost.

With dt_routing auto, plan_sub_reaches splits each river into the fewest sub-reaches, fixed for the run, that leave at
every discharge of its range a whole number of subcycles keeping every coefficient non-negative, and each runoff step of
each sub-reach takes the fewest subcycles its own k and x need. With dt_routing given, every river is routed whole and
each runoff step is taken in subcycles of exactly dt_routing.
"""

import warnings
from typing import TYPE_CHECKING, NamedTuple

import numba
import numpy as np
from numba.extending import overload

from ..network.Network import RESTRICTED_LARGEST_X, sub_reach_shares_and_x
from ..types import DynamicParameters, Float64Array, FloatArray, Int32Array, IntArray, RiverForcing
from ._numba_kernels import ComputePlan, add_forcing_terms, is_argument_type, k_and_x_at_discharges, route_river
from .static_muskingum import _add_outflow_to_downstream_inflow, _write_discharge_at_runoff_steps

if TYPE_CHECKING:
    from ..configs.Configs import Configs
    from ..network.Network import Network

__all__ = [
    'DynamicCoefficients',
    'Regression',
    'PARAMETERS_CLASS_FOR_PARAMETERS',
    'plan_sub_reaches',
    'prepare_routing',
]

# the discharges across each river's range, spread evenly in log, at which plan_sub_reaches checks its coefficients
PLAN_DISCHARGES = 64
# the most sub-reaches plan_sub_reaches tries before reporting a river it cannot make valid
MOST_SUB_REACHES = 256

# Polynomials, highest power first, fit to log2(m) of a float32 mantissa m in [1, 2) as a polynomial in m - 1.5, within
# 9e-7, and to 2^f of a fraction f in [0, 1) as a polynomial in f - 0.5, within 3e-7 of its value.
LOG2_OF_MANTISSA = np.array(
    [0.0145985, -0.0248256, 0.0375328, -0.0706862, 0.1425188, -0.3206248, 0.9617962, 0.5849627], dtype=np.float32
)
EXP2_OF_FRACTION = np.array([0.0018944, 0.0136765, 0.0784937, 0.3397254, 0.9802582, 1.4142137], dtype=np.float32)


class Regression(NamedTuple):
    """
    The regression equations of each river at the step dt its coefficients are built for, dt_routing when it is given
    and otherwise the runoff step, which the kernel divides into subcycles, as base 2 logs: dt/k = 2^(log2_dt_div_k -
    k_exponent * log2(Q)) and 1 - 2x = 2^(log2_one_minus_two_x + x_exponent * log2(Q)), held from
    smallest_one_minus_two_x to 1, with Q held within min_q and max_q. They are k = dynamicKCoefficient * Q **
    dynamicKExponent and x = 1/2 - dynamicXCoefficient * Q ** dynamicXExponent, held from 0 to the largest x routed.
    """

    log2_dt_div_k: FloatArray  # (n,) log2 of dt/k at a discharge of 1 m³/s
    k_exponent: FloatArray  # (n,)
    log2_one_minus_two_x: FloatArray  # (n,) log2 of 1 - 2x at a discharge of 1 m³/s
    x_exponent: FloatArray  # (n,)
    min_q: FloatArray  # (n,) m³/s
    max_q: FloatArray  # (n,) m³/s
    smallest_one_minus_two_x: np.float32  # 1 - 2x at the largest x routed

    @classmethod
    def at_step(cls, columns: tuple[FloatArray, ...], dt: int, largest_x: float) -> Regression:
        """
        The regression equations of the dynamic columns, in COLUMNS_FOR_PARAMETERS order, at the step dt, routing x no
        larger than largest_x.
        """
        k_coefficient, k_exponent, x_coefficient, x_exponent, min_q, max_q = columns
        # an x coefficient of 0, an x of 1/2 at every discharge, has a log2 of -inf, held at the smallest exponent
        # of a float32 so that 2 to its power is still a number
        log2_one_minus_two_x = np.maximum(np.log2(2 * x_coefficient.astype(np.float64)), -126)
        return cls(
            np.log2(dt / k_coefficient.astype(np.float64)).astype(np.float32),
            k_exponent.astype(np.float32),
            log2_one_minus_two_x.astype(np.float32),
            x_exponent.astype(np.float32),
            min_q.astype(np.float32),
            max_q.astype(np.float32),
            np.float32(1 - 2 * largest_x),
        )

    def dt_div_k_and_one_minus_two_x_over_range(self, discharges: int) -> tuple[Float64Array, Float64Array]:
        """(n, discharges) dt/k and 1 - 2x of each river at discharges spread evenly in log across its range."""
        log2_q = np.linspace(
            np.log2(self.min_q, dtype=np.float64), np.log2(self.max_q, dtype=np.float64), discharges, axis=1
        )
        dt_div_k = 2.0 ** (self.log2_dt_div_k[:, None] - self.k_exponent[:, None] * log2_q)
        one_minus_two_x = 2.0 ** (self.log2_one_minus_two_x[:, None] + self.x_exponent[:, None] * log2_q)
        return dt_div_k, np.clip(one_minus_two_x, self.smallest_one_minus_two_x, 1.0)


class DynamicCoefficients(NamedTuple):
    """
    What the coefficients of each river are rebuilt from at every runoff step: the parameters of a dynamic method, and
    the share of the river's k and runoff each of its sub-reaches takes, with the scale of its 1 - 2x that keeps the
    river's variance.
    """

    parameters: DynamicParameters
    share: FloatArray  # (sub-reaches,) of the river's k and of its runoff
    one_minus_two_x_scale: FloatArray  # (sub-reaches,) 1 / sum(w**2) of the shares w of the sub-reach's river
    inverse_dt_runoff: np.float32
    choose_subcycles: bool  # dt_routing auto: each step takes the fewest subcycles that keep it valid
    fixed_subcycles: np.float32  # the subcycles of every step when they are not chosen: dt_runoff / dt_routing, or 1


# the NamedTuple of each dynamic method's parameters, built from its COLUMNS_FOR_PARAMETERS in the order of its fields
PARAMETERS_CLASS_FOR_PARAMETERS: dict[str, type[DynamicParameters]] = {'dynamic': Regression}


def plan_sub_reaches(
    dt_div_k: Float64Array, one_minus_two_x: Float64Array
) -> tuple[IntArray, Float64Array, Float64Array, IntArray, np.ndarray, np.ndarray]:
    """
    The fewest sub-reaches each river must be split into, fixed for the run, so that at every one of its discharges
    some whole number of subcycles of the step dt keeps every coefficient of every sub-reach non-negative.

    A sub-reach of share w of k, at a discharge of dt/k = u and 1 - 2x = v, is routed at u/w and at the 1 - 2x of its
    layout, and its coefficients are non-negative at s subcycles when 1 - v <= u/s <= 1 + v: the fewest s that keep c3
    non-negative, ceil(u / (1 + v)), must also keep c1 non-negative. The sub-reaches are laid out one of two ways:

    - The shares of Network.sub_reach_k_and_x at the river's smallest x, with 1 - 2x scaled by 1/sum(w**2), keep its
      mean travel time and its variance at every discharge with an x of 0 or more.
    - Where no count of those is valid at every discharge, as for a river whose x is small at some discharges and
      large at others, equal shares with 1 - 2x scaled by their count and held at 1 or less keep the variance at the
      discharges where the x of the sub-reaches need not be held at 0.

    A river neither layout makes valid with up to MOST_SUB_REACHES sub-reaches keeps one reach and is reported as
    unresolvable, and the kernel keeps it valid by lowering its x at the steps that need it.

    Args:
        dt_div_k: (n, discharges) dt/k of each river at discharges across its range
        one_minus_two_x: (n, discharges) 1 - 2x of each river at the same discharges

    Returns:
        tuple: (sub_reaches, share, one_minus_two_x_scale, most_subcycles, resolvable, keeps_variance). share and
            one_minus_two_x_scale have a value per sub-reach, each river's upstream to downstream, most_subcycles is
            the most subcycles any discharge of each river needs, and keeps_variance is False for the rivers of equal
            shares.
    """
    n_rivers = dt_div_k.shape[0]
    sub_reaches = np.ones(n_rivers, dtype=np.int64)
    most_subcycles = np.ones(n_rivers, dtype=np.int64)
    resolvable = np.zeros(n_rivers, dtype=bool)
    equal_shares = np.zeros(n_rivers, dtype=bool)
    smallest_x = (1 - one_minus_two_x.max(axis=1)) / 2
    for keep_variance in (True, False):
        for count in range(1 if keep_variance else 2, MOST_SUB_REACHES + 1):
            rivers = np.flatnonzero(~resolvable)
            if rivers.size == 0:
                break
            if keep_variance:
                larger_share, last_share, _ = sub_reach_shares_and_x(np.full(rivers.size, count), smallest_x[rivers])
                scale = 1 / ((count - 1) * larger_share**2 + last_share**2)
            else:
                larger_share = last_share = np.full(rivers.size, 1 / count)
                scale = np.full(rivers.size, float(count))
            valid, subcycles = _valid_at_every_discharge(
                dt_div_k[rivers], one_minus_two_x[rivers], (larger_share, last_share), scale
            )
            if (
                count == 1
            ):  # a river no count makes valid keeps one reach and the subcycles that keep its c3 non-negative
                most_subcycles[rivers] = subcycles
            planned = rivers[valid]
            sub_reaches[planned] = count
            most_subcycles[planned] = subcycles[valid]
            resolvable[planned] = True
            equal_shares[planned] = not keep_variance
    larger_share, last_share, _ = sub_reach_shares_and_x(sub_reaches, smallest_x)
    larger_share[equal_shares] = last_share[equal_shares] = 1 / sub_reaches[equal_shares]
    scale = np.where(equal_shares, sub_reaches, 1 / ((sub_reaches - 1) * larger_share**2 + last_share**2))
    share = np.repeat(larger_share, sub_reaches)
    share[np.cumsum(sub_reaches) - 1] = last_share
    return sub_reaches, share, np.repeat(scale, sub_reaches), most_subcycles, resolvable, resolvable & ~equal_shares


def _valid_at_every_discharge(
    dt_div_k: Float64Array,
    one_minus_two_x: Float64Array,
    shares: tuple[Float64Array, Float64Array],
    scale: Float64Array,
) -> tuple[np.ndarray, Float64Array]:
    """
    Whether the fewest subcycles that keep c3 of the longest and the shortest sub-reach non-negative keep their c1
    non-negative at every discharge, with a margin for the rounding of float32, and the most subcycles they take. A
    share of 0, which keeping the variance of a river whose x is 0 at some discharge can ask for, is never valid.
    """
    valid = np.ones(dt_div_k.shape[0], dtype=bool)
    subcycles = np.ones(dt_div_k.shape[0])
    v = np.minimum(one_minus_two_x * scale[:, np.newaxis], 1.0)
    for share in shares:
        positive = share > 1e-9
        u = dt_div_k / np.where(positive, share, 1.0)[:, np.newaxis]
        fewest = np.maximum(1, np.ceil(u / (1 + v)))
        valid &= positive & np.all(u / fewest >= (1 - v) * (1 + 1e-6), axis=1)
        subcycles = np.maximum(subcycles, fewest.max(axis=1))
    return valid, subcycles


def prepare_routing(
    network: Network, configs: Configs, dt_runoff: int, dt_routing: int | None
) -> tuple[ComputePlan, DynamicCoefficients]:
    """
    Each river's compute plan and what its coefficients are rebuilt from at every runoff step. With dt_routing auto
    (None), each river is split into the sub-reaches of plan_sub_reaches and each step takes the fewest subcycles its
    own k and x need. With dt_routing given, every river is routed whole and every step is taken in subcycles of it.
    Rivers with a negative coefficient at some discharge of their range are reported as invalid_coefficients says.

    Raises:
        ValueError: if the network file lacks a column of the parameters, or invalid_coefficients is raise and a river
            has a negative coefficient
    """
    columns = network.columns_for_parameters(configs.parameters)
    largest_x = RESTRICTED_LARGEST_X if configs.restrict_x_range else 0.5
    step = dt_routing or dt_runoff  # the step the coefficients are built for, which chosen subcycles divide
    choose_subcycles = not dt_routing
    fixed_subcycles = dt_runoff // dt_routing if dt_routing else 1
    parameters = PARAMETERS_CLASS_FOR_PARAMETERS[configs.parameters].at_step(columns, step, largest_x)
    dt_div_k, one_minus_two_x = parameters.dt_div_k_and_one_minus_two_x_over_range(PLAN_DISCHARGES)
    if choose_subcycles:
        sub_reaches, share, scale, most_subcycles, resolvable, keeps_variance = plan_sub_reaches(
            dt_div_k, one_minus_two_x
        )
    else:  # a river routed whole at dt_routing is valid where dt/k is within 1 - 2x of 1 at every discharge
        sub_reaches = np.ones(network.size, dtype=np.int64)
        most_subcycles = np.full(network.size, fixed_subcycles, dtype=np.int64)
        share = scale = np.ones(network.size)
        resolvable = np.all((1 - one_minus_two_x <= dt_div_k) & (dt_div_k <= 1 + one_minus_two_x), axis=1)
        keeps_variance = np.ones(network.size, dtype=bool)
    _report_invalid_rivers(network, resolvable, keeps_variance, step, configs, choose_subcycles)
    coefficients = DynamicCoefficients(
        parameters,
        share.astype(np.float32),
        scale.astype(np.float32),
        np.float32(1.0 / dt_runoff),
        choose_subcycles,
        np.float32(fixed_subcycles),
    )
    return ComputePlan.from_counts(sub_reaches, most_subcycles), coefficients


def _report_invalid_rivers(
    network: Network,
    resolvable: np.ndarray,
    keeps_variance: np.ndarray,
    step: int,
    configs: Configs,
    choose_subcycles: bool,
) -> None:
    """
    Warn of, raise on, or ignore, as invalid_coefficients says, the rivers with a negative coefficient at some
    discharge, and with dt_routing auto the rivers whose sub-reaches do not keep their variance at every discharge.
    """
    invalid = int(np.count_nonzero(~resolvable))
    equal_shares = int(np.count_nonzero(resolvable & ~keeps_variance))
    if invalid + equal_shares == 0 or configs.invalid_coefficients == 'ignore':
        return
    if choose_subcycles:
        message = (
            f'At a step of {step} s, {equal_shares + invalid} of {network.size} rivers have no sub-reaches that keep '
            f'their variance and every coefficient non-negative at every discharge. {equal_shares} are split into '
            f'equal sub-reaches whose x is held at 0 where keeping the variance would take it lower. {invalid} have '
            f'no valid split into up to {MOST_SUB_REACHES} sub-reaches either; they are routed whole, and at the '
            f'steps where no subcycles keep every coefficient non-negative their x is lowered until they do.'
        )
    else:
        message = (
            f'{invalid} of {network.size} rivers have a negative Muskingum coefficient at some discharge at '
            f'dt_routing={step} s, at which every river is routed whole; set dt_routing to auto to split them into '
            f'sub-reaches with subcycles that keep every coefficient non-negative.'
        )
    message += " Set invalid_coefficients to 'ignore' to silence this."
    if configs.invalid_coefficients == 'raise':
        raise ValueError(message)
    warnings.warn(message, stacklevel=3)


@numba.njit(cache=True, nogil=True, fastmath=True)
def _log2(values: FloatArray, out: FloatArray, bits: Int32Array) -> None:
    """
    log2 of each positive float32 of values into out: its exponent plus LOG2_OF_MANTISSA of its mantissa. Unlike
    np.log2, which numba calls one value at a time, the compiler vectorizes it.
    """
    value_bits = values.view(np.int32)
    for i in range(values.shape[0]):
        bits[i] = (value_bits[i] & 0x007FFFFF) | 0x3F800000  # the mantissa with the exponent of 1
    mantissas = bits.view(np.float32)
    for i in range(values.shape[0]):
        exponent = np.float32(((value_bits[i] >> 23) & 0xFF) - 127)
        f = mantissas[i] - np.float32(1.5)
        polynomial = LOG2_OF_MANTISSA[0]
        for power in range(1, LOG2_OF_MANTISSA.shape[0]):
            polynomial = polynomial * f + LOG2_OF_MANTISSA[power]
        out[i] = exponent + polynomial


@numba.njit(cache=True, nogil=True, fastmath=True)
def _exp2(powers: FloatArray, out: FloatArray, bits: Int32Array) -> None:
    """
    2 to each of powers into out: EXP2_OF_FRACTION of its fraction times 2 to its whole part, built in the exponent
    bits of a float32. Powers are held within the exponents of a normal float32, -126 to 127.
    """
    for i in range(powers.shape[0]):
        power = min(max(powers[i], np.float32(-126.0)), np.float32(127.0))
        whole = np.floor(power)
        f = power - whole - np.float32(0.5)
        polynomial = EXP2_OF_FRACTION[0]
        for power in range(1, EXP2_OF_FRACTION.shape[0]):
            polynomial = polynomial * f + EXP2_OF_FRACTION[power]
        out[i] = polynomial
        bits[i] = (np.int32(whole) + 127) << 23
    scales = bits.view(np.float32)
    for i in range(powers.shape[0]):
        out[i] *= scales[i]


@numba.njit(cache=True, nogil=True, fastmath={'contract'})
def _write_step_coefficients(
    share: np.float32,
    one_minus_two_x_scale: np.float32,
    choose_subcycles: bool,
    fixed_subcycles: np.float32,
    dt_div_k: FloatArray,
    one_minus_two_x: FloatArray,
    c1: FloatArray,
    c3: FloatArray,
    subcycles: FloatArray,
) -> None:
    """
    c1, c3, and the subcycles of each runoff step of a sub-reach. The river's dt/k and 1 - 2x at the step, written u
    and v, are the sub-reach's u / share and v * one_minus_two_x_scale, held at 1 or less, which is an x of 0 or more.
    Each subcycle routes with c1 = (u - 1 + v) / D, c2 = (u + 1 - v) / D, and c3 = (1 + v - u) / D, with D = u + 1 + v
    and u divided by the subcycles; c2 is 1 - c1 - c3, held at 0 or more against rounding.

    Without choose_subcycles a step takes fixed_subcycles, of the step its u is of. With choose_subcycles, a step of
    the runoff step takes the fewest subcycles s that keep c3 non-negative, ceil(u / (1 + v)). Where
    those leave c1 negative, at a discharge whose window of valid steps holds no whole number of subcycles, the step's
    1 - 2x is raised to 1 - u/s, the least that makes c1 0, which lowers its x. plan_sub_reaches chooses sub-reaches
    that need this at none of the discharges it checks, so it is rare and small. The numerators are held at 0 or more
    against rounding, and D is their sum, so the coefficients still sum to 1. No step depends on another, so the
    compiler evaluates several at a time.
    """
    zero = np.float32(0.0)
    one = np.float32(1.0)
    inverse_share = one / share
    for step in range(dt_div_k.shape[0]):
        u = dt_div_k[step] * inverse_share
        v = min(one_minus_two_x[step] * one_minus_two_x_scale, one)
        steps = fixed_subcycles
        if choose_subcycles:
            steps = max(one, np.ceil(u / (one + v)))
            u = u / steps
            v = max(v, one - u)
            c1_numerator = max(u - one + v, zero)
            c3_numerator = max(one + v - u, zero)
        else:
            c1_numerator = u - one + v
            c3_numerator = one + v - u
        inverse = one / (c1_numerator + (u + one - v) + c3_numerator)
        c1[step] = c1_numerator * inverse
        c3[step] = c3_numerator * inverse
        subcycles[step] = steps


@numba.njit(cache=True, nogil=True, fastmath={'contract'})
def _route_sub_reach(
    outflow: np.float32,
    share: np.float32,
    c1: FloatArray,
    c3: FloatArray,
    subcycles: FloatArray,
    inflow: FloatArray,
    lateral_inflow: FloatArray,
    outflow_at_steps: FloatArray,
) -> np.float32:
    """
    Route a sub-reach from Q = outflow with the coefficients and subcycles of each runoff step, its inflow
    interpolated linearly across the subcycles, writing its Q at the end of each step into outflow_at_steps, and return
    its last Q.
    """
    zero = np.float32(0.0)
    one = np.float32(1.0)
    for step in range(c1.shape[0]):
        step_c1 = c1[step]
        step_c3 = c3[step]
        step_c2 = max(one - step_c1 - step_c3, zero)  # held at 0 or more against rounding
        runoff_term = (step_c1 + step_c2) * share * lateral_inflow[step]
        inflow_before = inflow[step]
        inflow_at_end = inflow[step + 1]
        steps = int(subcycles[step])
        if steps == 1:
            outflow = step_c3 * outflow + (step_c1 * inflow_at_end + step_c2 * inflow_before + runoff_term)
        else:
            inverse_steps = one / np.float32(steps)
            for subcycle in range(steps):
                # a weighted mean of the inflow at the start and end of the step, which never leaves the range between
                # them however it rounds
                after = np.float32(subcycle + 1) * inverse_steps
                inflow_after = (one - after) * inflow[step] + after * inflow_at_end
                outflow = step_c3 * outflow + (step_c1 * inflow_after + step_c2 * inflow_before + runoff_term)
                inflow_before = inflow_after
        outflow_at_steps[step] = outflow
    return outflow


def _route_river_with_dynamic_coefficients(
    river_index: int,
    compute_plan: ComputePlan,
    coefficients: DynamicCoefficients,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    discharge_at_runoff_steps: FloatArray,
    scratch: FloatArray,
) -> None:
    """
    Route the river at river_index in its sub-reaches in series, each sub-reach's whole series before the next, with
    the coefficients of each step at the discharge entering the sub-reach over the step: the mean of its inflow at the
    start and end of the step plus its share of the river's runoff.
    """
    parameters = coefficients.parameters
    runoff_steps = discharge_at_runoff_steps.shape[0]  # upstream_inflow also holds I before the first runoff step
    lateral_inflow = scratch[0, :runoff_steps]  # m³/s of the river's runoff over each step
    discharge = scratch[1, :runoff_steps]  # the discharge each step's k and x are evaluated at
    dt_div_k = scratch[2, :runoff_steps]
    one_minus_two_x = scratch[3, :runoff_steps]
    work = scratch[6:8]  # the dynamic method's own rows, whole: numba types a 2D slice of their columns as any layout
    half = np.float32(0.5)
    lateral_inflow[:] = 0.0
    if river_forcing is not None:
        add_forcing_terms(river_forcing, coefficients.inverse_dt_runoff, runoff_steps, 1, lateral_inflow)
    first_sub_reach = compute_plan.sub_reach_starts[river_index]
    last_sub_reach = compute_plan.sub_reach_starts[river_index + 1] - 1
    inflow = upstream_inflow
    for sub_reach in range(first_sub_reach, last_sub_reach + 1):
        share = coefficients.share[sub_reach]
        for step in range(runoff_steps):
            discharge[step] = half * (inflow[step] + inflow[step + 1]) + share * lateral_inflow[step]
        k_and_x_at_discharges(parameters, river_index, discharge, dt_div_k, one_minus_two_x, work)
        # the sub-reach's Q before the first step and at the end of each, the inflow of the next sub-reach. Two rows
        # take turns, since a sub-reach reads its inflow row while it writes its outflow row.
        c1 = discharge  # k and x have been evaluated at it
        c3 = scratch[6, :runoff_steps]
        subcycles = scratch[7, :runoff_steps]
        _write_step_coefficients(
            share,
            coefficients.one_minus_two_x_scale[sub_reach],
            coefficients.choose_subcycles,
            coefficients.fixed_subcycles,
            dt_div_k,
            one_minus_two_x,
            c1,
            c3,
            subcycles,
        )
        outflow = scratch[4 + (sub_reach - first_sub_reach) % 2, : runoff_steps + 1]
        outflow[0] = channel_state[sub_reach]
        channel_state[sub_reach] = _route_sub_reach(
            outflow[0], share, c1, c3, subcycles, inflow, lateral_inflow, outflow[1:]
        )
        inflow = outflow
    _add_outflow_to_downstream_inflow(inflow[0], inflow[1:], downstream_inflow)
    _write_discharge_at_runoff_steps(inflow[1:], 1, discharge_at_runoff_steps)


@overload(route_river, jit_options={'nogil': True, 'fastmath': {'contract'}})
def _overload_route_river_for_dynamic_coefficients(
    river_index: int,
    compute_plan: ComputePlan,
    coefficients: DynamicCoefficients,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    discharge_at_runoff_steps: FloatArray,
    scratch: FloatArray,
):
    if is_argument_type(coefficients, DynamicCoefficients):
        return _route_river_with_dynamic_coefficients


def _k_and_x_by_regression(
    parameters: Regression,
    river_index: int,
    discharge: FloatArray,
    dt_div_k: FloatArray,
    one_minus_two_x: FloatArray,
    work: FloatArray,
) -> None:
    """
    The dt/k and 1 - 2x of the regression equations at each discharge Q, held within the discharges the regression was
    fit over, so neither equation is extrapolated and the log of Q is that of a positive number. 1 - 2x is held at 1 or
    less, since Muskingum x is never negative, and at the largest x routed or more. The discharges are overwritten.

    Each of _log2 and _exp2 writes into a row other than the one it reads: the compiler vectorizes a loop over two
    arrays only after checking at run time that they do not overlap, and one array passed as both fails the check.
    """
    steps = discharge.shape[0]
    log2_discharge = work[0][:steps]
    bits = work[1][:steps].view(np.int32)
    min_q = parameters.min_q[river_index]
    max_q = parameters.max_q[river_index]
    for step in range(discharge.shape[0]):
        discharge[step] = min(max(discharge[step], min_q), max_q)
    _log2(discharge, log2_discharge, bits)
    log2_dt_div_k = parameters.log2_dt_div_k[river_index]
    k_exponent = parameters.k_exponent[river_index]
    log2_one_minus_two_x = parameters.log2_one_minus_two_x[river_index]
    x_exponent = parameters.x_exponent[river_index]
    power_of_dt_div_k = discharge
    power_of_one_minus_two_x = log2_discharge
    for step in range(discharge.shape[0]):
        log2_q = log2_discharge[step]
        power_of_dt_div_k[step] = log2_dt_div_k - k_exponent * log2_q
        power_of_one_minus_two_x[step] = log2_one_minus_two_x + x_exponent * log2_q
    _exp2(power_of_dt_div_k, dt_div_k, bits)
    _exp2(power_of_one_minus_two_x, one_minus_two_x, bits)
    smallest_one_minus_two_x = parameters.smallest_one_minus_two_x
    for step in range(discharge.shape[0]):
        one_minus_two_x[step] = min(max(one_minus_two_x[step], smallest_one_minus_two_x), np.float32(1.0))


@overload(k_and_x_at_discharges, jit_options={'nogil': True})
def _overload_k_and_x_at_discharges_for_regression(
    parameters: Regression,
    river_index: int,
    discharge: FloatArray,
    dt_div_k: FloatArray,
    one_minus_two_x: FloatArray,
    work: FloatArray,
):
    if is_argument_type(parameters, Regression):
        return _k_and_x_by_regression
