"""
The arguments of route_scheduled_rivers, grouped by what they describe. numba compiles NamedTuple fields as plain
arguments, so grouping them costs nothing at run time.
"""

from typing import NamedTuple

import numpy as np

__all__ = ['StaticCoefficients', 'DynamicCoefficients', 'Layout', 'Schedule', 'CatchmentByRiver']


class StaticCoefficients(NamedTuple):
    """Muskingum coefficients per river, fixed for the whole run."""

    c1: np.ndarray  # (n,) float32
    c2: np.ndarray  # (n,) float32
    c3: np.ndarray  # (n,) float32
    c4_dt: np.ndarray  # (n,) float32, c4 / dt_runoff: turns a catchment runoff volume into a uniform inflow rate


class DynamicCoefficients(NamedTuple):
    """Nonlinear Muskingum parameters per river; the coefficients are rebuilt from K = alpha * Q ** beta each step."""

    alpha: np.ndarray  # (n,)
    beta: np.ndarray  # (n,)
    x: np.ndarray  # (n,)
    dt_routing: np.float32
    inv_dt_runoff: np.float32


class Layout(NamedTuple):
    """
    How each river is routed. Both are empty on a standard network, where ``q_t`` holds one state per river. On a
    stabilized network river r is the sub-reaches ``reach_indptr[r]:reach_indptr[r + 1]`` in series, with one state
    per sub-reach in ``q_t``, and ``substeps`` gives the steps each river takes per routing step, or is empty when no
    river is sub-cycled.
    """

    reach_indptr: np.ndarray  # (n + 1,) int64, or empty
    substeps: np.ndarray  # (n,) int64, or empty


class Schedule(NamedTuple):
    """
    The blocks of rivers one pass routes, in order. A block's outlet (-1 for none) hands its unclamped series to
    ``boundary[block_region]`` instead of an inflow row, and each boundary row in ``cut_target`` is injected into its
    target river before that river is routed. ``out_row`` is empty when every river has a discharge row, and otherwise
    gives each river's row, -1 for a synthetic river whose series is written to a scratch row and discarded.
    """

    block_starts: np.ndarray  # (n_blocks,) int32
    block_stops: np.ndarray  # (n_blocks,) int32
    block_outlet: np.ndarray  # (n_blocks,) int32
    block_region: np.ndarray  # (n_blocks,) int32
    cut_target: np.ndarray  # (n_regions,) int32, or empty
    boundary: np.ndarray  # (n_regions, n_routing + 1) float32
    out_row: np.ndarray  # (n,) int32, or empty


class CatchmentByRiver(NamedTuple):
    """Catchment runoff volumes as a C-order (river, time) array, each river's row read in place with no copy."""

    runoff: np.ndarray  # (n, n_steps) float32
