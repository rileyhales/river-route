"""
Routes a region river by river. Rivers hand their discharge downstream at the end of every runoff step, and each river
computes in a plan of its own, the ComputePlan of its sub-reaches and of the subcycles it takes in each runoff step. Two
configs decide how a river is routed, and each is the type of an argument numba compiles route_job for: ``forcing``
(None for channel routing) and ``parameters`` (the coefficients).
"""

from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING, NamedTuple

import numba
import numpy as np
from numba import types
from numba.extending import overload

from ..types import (
    Coefficients,
    DynamicParameters,
    FloatArray,
    Forcing,
    Int32Array,
    IntArray,
    JobBlocks,
    RiverForcing,
    Transform,
)
from ._simd import FourLanes

if TYPE_CHECKING:
    from ..network.Network import Network

__all__ = [
    # the stages each river is routed through, implemented next to the types they read
    'forcing_of_river',
    'transform_runoff',
    'add_forcing_terms',
    'add_forcing_of_four_steps',
    'route_river',
    'k_and_x_at_discharges',
    'is_argument_type',
    # how each river is computed, what a job routes, and routing the jobs of a region
    'ComputePlan',
    'Job',
    'route_job',
    'route_region',
]

# the rows of scratch route_job gives every river, each as long as the river's steps: two for static coefficients, eight
# for dynamic ones
SCRATCH_ROWS = 8


class ComputePlan(NamedTuple):
    """How each river is computed: split into sub-reaches in series, and in subcycles of each runoff step."""

    sub_reach_starts: Int32Array  # (n + 1,) channel_state index of each river's first sub-reach, then the total
    subcycles: Int32Array  # (n,) steps each river takes per runoff step, the most it takes with dynamic coefficients

    @classmethod
    def from_counts(cls, sub_reaches: IntArray, subcycles: IntArray) -> ComputePlan:
        """The plan of rivers split into ``sub_reaches`` sub-reaches that take ``subcycles`` steps per runoff step."""
        sub_reach_starts = np.zeros(sub_reaches.shape[0] + 1, dtype=np.int32)
        sub_reach_starts[1:] = np.cumsum(sub_reaches)
        return cls(sub_reach_starts, subcycles.astype(np.int32))


class Job(NamedTuple):
    """
    The blocks of rivers one thread routes, in order. A block's outlet (-1 for none) adds its series into its boundary
    row, ``boundary[block_number]``, instead of an inflow row, and the main stem job adds each boundary row into the
    inflow of its ``boundary_target`` river before routing that river.
    """

    block_starts: Int32Array  # (n_job_blocks,)
    block_stops: Int32Array  # (n_job_blocks,)
    block_outlet: Int32Array  # (n_job_blocks,)
    block_number: Int32Array  # (n_job_blocks,)
    boundary_target: Int32Array  # (n_blocks,) river each boundary row drains into, -1 at a basin outlet, or empty
    boundary: FloatArray  # (n_blocks, runoff_steps + 1)
    discharge_row: Int32Array  # (n,) row of discharge_array of each river, -1 for a synthetic river, or empty for none


_EMPTY_INT32 = np.zeros(0, dtype=np.int32)  # no boundary targets, or no synthetic rivers


def forcing_of_river(forcing: Forcing, river_index: int) -> RiverForcing | None:
    """The forcing of the river at river_index, in the form add_forcing_terms reads, or None for channel routing."""
    raise NotImplementedError('stages run only inside kernels, where numba compiles them from their overloads')


def transform_runoff(
    transform: Transform, river_forcing: RiverForcing | None, river_index: int, out: FloatArray
) -> RiverForcing | None:
    """
    The forcing of the river at river_index after a runoff transform, which may write it into out. No transform but
    the uniform one exists yet.
    """
    raise NotImplementedError('stages run only inside kernels, where numba compiles them from their overloads')


def add_forcing_terms(
    river_forcing: RiverForcing,
    c4: np.float32,
    runoff_steps: int,
    river_steps_per_runoff_step: int,
    known_terms: FloatArray,
) -> None:
    """Add c4 V_t, c4 times the river's runoff volume (m³) of runoff step t, to known_terms at each step t spans."""
    raise NotImplementedError('stages run only inside kernels, where numba compiles them from their overloads')


def add_forcing_of_four_steps(
    river_forcing: RiverForcing, c4: np.float32, first_step: int, steps: int, terms: FourLanes
) -> FourLanes:
    """
    terms plus c4 V_t in each lane, for the up to four runoff steps t from first_step, a lane each. Lanes past the last
    step of the series are left as they are.
    """
    raise NotImplementedError('stages run only inside kernels, where numba compiles them from their overloads')


def route_river(
    river_index: int,
    compute_plan: ComputePlan,
    coefficients: Coefficients,
    river_forcing: RiverForcing | None,
    channel_state: FloatArray,
    upstream_inflow: FloatArray,
    downstream_inflow: FloatArray,
    discharge_at_runoff_steps: FloatArray,
    scratch: FloatArray,
) -> None:
    """
    Route the river at river_index over every runoff step: add its Q at the end of each into downstream_inflow, write
    it into discharge_at_runoff_steps, and leave its last Q in channel_state. scratch is (SCRATCH_ROWS, runoff steps
    times the most subcycles + 1) rows the routing method works in.
    """
    raise NotImplementedError('stages run only inside kernels, where numba compiles them from their overloads')


def k_and_x_at_discharges(
    parameters: DynamicParameters,
    river_index: int,
    discharge: FloatArray,
    dt_div_k: FloatArray,
    one_minus_two_x: FloatArray,
    work: FloatArray,
) -> None:
    """
    dt/k and 1 - 2x of the river at river_index at each of the discharges, by the dynamic method's parameters, the
    forms its coefficients are built from. The discharges and the rows of work are scratch the method may overwrite.
    """
    raise NotImplementedError('stages run only inside kernels, where numba compiles them from their overloads')


def is_argument_type(numba_type: object, kind: type | None) -> bool:
    """
    Whether numba typed a stage argument as ``kind``: None, np.ndarray for a 1D array, or a NamedTuple class. An
    overload is handed numba types, though it is annotated like the implementation it returns, as numba requires.
    """
    if kind is None:
        return isinstance(numba_type, types.NoneType)
    if kind is np.ndarray:
        return isinstance(numba_type, types.Array) and numba_type.ndim == 1
    return getattr(numba_type, 'instance_class', None) is kind


def _no_forcing_for_channel_routing(forcing: None, river_index: int) -> None:
    """Channel routing routes the inflow from upstream alone."""
    return None


@overload(forcing_of_river)
def _overload_forcing_of_river_for_channel_routing(forcing: None, river_index: int):
    if is_argument_type(forcing, None):
        return _no_forcing_for_channel_routing


@numba.njit(cache=True, nogil=True)
def sort_boundaries_by_target_river(boundary_target: Int32Array) -> IntArray:
    """Positions into boundary_target in order of the river each boundary drains into, basin outlets (-1) dropped."""
    order = np.argsort(boundary_target, kind='mergesort')
    first = 0
    while first < order.shape[0] and boundary_target[order[first]] < 0:
        first += 1
    return order[first:]


@numba.njit(cache=True, nogil=True)
def first_river_and_span_of_job(block_starts: Int32Array, block_stops: Int32Array) -> tuple[int, int]:
    """The first river a job routes and how many river indices its blocks span, which sizes its per-river arrays."""
    first_river = int(block_starts.min()) if block_starts.shape[0] else 0
    last_river = int(block_stops.max()) if block_stops.shape[0] else 0
    return first_river, max(last_river - first_river, 0)


@numba.njit(cache=True, nogil=True)
def count_peak_live_inflow_rows(
    downstream_indices: Int32Array,
    block_starts: Int32Array,
    block_stops: Int32Array,
    block_outlet: Int32Array,
    boundary_target: Int32Array,
) -> int:
    """The most inflow rows open at once while a job routes its blocks in order, which sizes its pool of rows."""
    first_river, river_span = first_river_and_span_of_job(block_starts, block_stops)
    has_inflow_row = np.zeros(river_span, dtype=np.bool_)
    boundaries_in_target_order = sort_boundaries_by_target_river(boundary_target)
    next_boundary = 0
    open_rows = 0
    most_open_rows = 0
    for block in range(block_starts.shape[0]):
        for river_index in range(block_starts[block], block_stops[block]):
            # the main stem job opens a row for each river a block outlet drains into
            while (
                next_boundary < boundaries_in_target_order.shape[0]
                and boundary_target[boundaries_in_target_order[next_boundary]] == river_index
            ):
                if not has_inflow_row[river_index - first_river]:
                    has_inflow_row[river_index - first_river] = True
                    open_rows += 1
                next_boundary += 1
            # a block outlet adds its series into a boundary row, so it never opens an inflow row
            downstream = downstream_indices[river_index]
            if downstream >= 0 and river_index != block_outlet[block] and not has_inflow_row[downstream - first_river]:
                has_inflow_row[downstream - first_river] = True
                open_rows += 1
            most_open_rows = max(most_open_rows, open_rows)
            if has_inflow_row[river_index - first_river]:  # the river has been routed, so its row closes
                open_rows -= 1
    return most_open_rows


@numba.njit(cache=True, nogil=True)
def allocate_inflow_row_pool(
    river_span: int, n_rows: int, runoff_steps: int
) -> tuple[FloatArray, IntArray, IntArray, IntArray]:
    """
    The inflow rows, the row each river holds (-1 for none), and the stack of free rows with its count. Two rows follow
    the pool: the zeros a headwater reads as its inflow, and the row basin outlets add their series into.
    """
    inflow_rows = np.empty((n_rows + 2, runoff_steps + 1), dtype=np.float32)
    inflow_rows[n_rows, :] = 0.0
    inflow_row_of_river = np.full(river_span, -1, dtype=np.int64)
    free_inflow_rows = np.arange(n_rows)
    free_inflow_row_count = np.full(1, n_rows, dtype=np.int64)  # an array, so taking and releasing a row changes it
    return inflow_rows, inflow_row_of_river, free_inflow_rows, free_inflow_row_count


@numba.njit(cache=True, nogil=True)
def get_upstream_inflow_row(river: int, inflow_rows: FloatArray, inflow_row_of_river: IntArray) -> FloatArray:
    """The inflow row a river reads: the one its upstream rivers added into, or the zeros for a headwater."""
    row = inflow_row_of_river[river]
    return inflow_rows[row] if row >= 0 else inflow_rows[inflow_rows.shape[0] - 2]


@numba.njit(cache=True, nogil=True)
def get_or_open_downstream_inflow_row(
    river: int,
    inflow_rows: FloatArray,
    inflow_row_of_river: IntArray,
    free_inflow_rows: IntArray,
    free_inflow_row_count: IntArray,
) -> FloatArray:
    """The inflow row a river collects its upstream series in, taken from the pool and zeroed on first use, or the row
    basin outlets add into if river < 0."""
    if river < 0:
        return inflow_rows[inflow_rows.shape[0] - 1]
    row = inflow_row_of_river[river]
    if row < 0:
        free_inflow_row_count[0] -= 1
        row = free_inflow_rows[free_inflow_row_count[0]]
        inflow_row_of_river[river] = row
        inflow_rows[row, :] = 0.0
    return inflow_rows[row]


@numba.njit(cache=True, nogil=True)
def release_inflow_row_to_pool(
    river: int, inflow_row_of_river: IntArray, free_inflow_rows: IntArray, free_inflow_row_count: IntArray
) -> None:
    """A river has been routed, so the inflow row it read goes back to the pool."""
    row = inflow_row_of_river[river]
    if row >= 0:
        free_inflow_rows[free_inflow_row_count[0]] = row
        free_inflow_row_count[0] += 1
        inflow_row_of_river[river] = -1


@numba.njit(cache=True, nogil=True)
def route_job(
    channel_state: FloatArray,
    discharge_array: FloatArray,
    downstream_indices: Int32Array,
    coefficients: Coefficients,
    compute_plan: ComputePlan,
    forcing: Forcing,
    transform: Transform,
    job: Job,
) -> None:
    """
    Route the rivers of each block of job in index order, each river's whole series before the next river. Each river's
    discharge at the end of each runoff step is written into its own row of the C-order (river, time) discharge_array.
    """
    runoff_steps = discharge_array.shape[1]
    block_starts, block_stops, block_outlet, block_number, boundary_target, boundary, discharge_row = job
    n_rows = count_peak_live_inflow_rows(downstream_indices, block_starts, block_stops, block_outlet, boundary_target)
    first_river, river_span = first_river_and_span_of_job(block_starts, block_stops)
    inflow_rows, inflow_row_of_river, free_inflow_rows, free_inflow_row_count = allocate_inflow_row_pool(
        river_span, n_rows, runoff_steps
    )
    boundaries_in_target_order = sort_boundaries_by_target_river(boundary_target)
    next_boundary = 0
    transform_scratch = np.empty(0 if transform is None else runoff_steps, dtype=np.float32)
    scratch = np.empty((SCRATCH_ROWS, runoff_steps * compute_plan.subcycles.max() + 1), dtype=np.float32)
    has_synthetic_rivers = discharge_row.shape[0] > 0
    synthetic_river_discharge = np.empty(runoff_steps if has_synthetic_rivers else 0, dtype=discharge_array.dtype)

    for block in range(block_starts.shape[0]):
        for river_index in range(block_starts[block], block_stops[block]):
            # the main stem job adds the series of each block outlet that drains into this river before routing it
            while (
                next_boundary < boundaries_in_target_order.shape[0]
                and boundary_target[boundaries_in_target_order[next_boundary]] == river_index
            ):
                inflow_of_river = get_or_open_downstream_inflow_row(
                    river_index - first_river, inflow_rows, inflow_row_of_river, free_inflow_rows, free_inflow_row_count
                )
                block_outlet_series = boundary[boundaries_in_target_order[next_boundary]]
                for routing_step in range(inflow_of_river.shape[0]):
                    inflow_of_river[routing_step] += block_outlet_series[routing_step]
                next_boundary += 1
            upstream_inflow = get_upstream_inflow_row(river_index - first_river, inflow_rows, inflow_row_of_river)
            if river_index == block_outlet[block]:
                downstream_inflow = boundary[block_number[block]]
                downstream_inflow[:] = 0.0
            else:
                downstream = downstream_indices[river_index]
                downstream_inflow = get_or_open_downstream_inflow_row(
                    downstream - first_river if downstream >= 0 else -1,
                    inflow_rows,
                    inflow_row_of_river,
                    free_inflow_rows,
                    free_inflow_row_count,
                )
            row = discharge_row[river_index] if has_synthetic_rivers else river_index
            discharge_at_runoff_steps = discharge_array[row] if row >= 0 else synthetic_river_discharge
            river_forcing = forcing_of_river(forcing, river_index)
            if transform is not None:  # numba compiles the uniform transform, None, without this
                river_forcing = transform_runoff(transform, river_forcing, river_index, transform_scratch)
            route_river(
                river_index,
                compute_plan,
                coefficients,
                river_forcing,
                channel_state,
                upstream_inflow,
                downstream_inflow,
                discharge_at_runoff_steps,
                scratch,
            )
            release_inflow_row_to_pool(
                river_index - first_river, inflow_row_of_river, free_inflow_rows, free_inflow_row_count
            )
    return


def route_region(
    network: Network,
    coefficients: Coefficients,
    compute_plan: ComputePlan,
    channel_state: FloatArray,
    discharge_array: FloatArray,
    forcing: Forcing,
    thread_pool: ThreadPoolExecutor | None,
    threads: int,
) -> None:
    """
    Route every job of sub-watershed blocks, concurrently on thread_pool when given, then the main stem job, which adds
    the series the block outlets left in their boundary rows into the rivers they drain into. discharge_array has a
    column per runoff step, and forcing is None for channel routing or what a Runoff generator yielded.
    """
    # the Network derives the blocks once per thread count and caches them, so this is a lookup after the first file
    job_blocks, boundary_target = network.routing_blocks(threads if thread_pool is not None else 1)
    _check_everything_a_job_reads(
        network, coefficients, compute_plan, channel_state, discharge_array, forcing, job_blocks
    )
    downstream_indices = network.downstream_indices
    # a boundary row holds a block outlet's Q before the first runoff step and at the end of each, like an inflow row
    boundary = np.zeros((max(boundary_target.shape[0], 1), discharge_array.shape[1] + 1), dtype=np.float32)
    synthetic = network.synthetic
    discharge_row = (
        _EMPTY_INT32 if synthetic is None else np.where(synthetic, -1, np.cumsum(~synthetic) - 1).astype(np.int32)
    )
    transform = None  # the uniform transform

    def route_one_job(blocks: JobBlocks, boundary_target_of_job: Int32Array = _EMPTY_INT32) -> None:
        route_job(
            channel_state,
            discharge_array,
            downstream_indices,
            coefficients,
            compute_plan,
            forcing,
            transform,
            Job(*blocks, boundary_target_of_job, boundary, discharge_row),
        )

    *sub_watershed_jobs, main_stem = job_blocks
    if thread_pool is None:
        for blocks in sub_watershed_jobs:
            route_one_job(blocks)
    else:
        list(thread_pool.map(route_one_job, sub_watershed_jobs))  # list() so a worker exception raises

    # the one barrier of the simulation: every block has finished before the main stem reads their boundary rows
    route_one_job(main_stem, boundary_target)
    return


def _check_everything_a_job_reads(
    network: Network,
    coefficients: Coefficients,
    compute_plan: ComputePlan,
    channel_state: FloatArray,
    discharge_array: FloatArray,
    forcing: Forcing,
    job_blocks: tuple[JobBlocks, ...],
) -> None:
    """Last chance to validate before routing. The numba kernels will not check or raise useful errors."""
    # each river's series is written into its own contiguous row, which is the layout the jobs solve in.
    # Output is never copied, so any other layout is refused rather than transposed.
    if not discharge_array.flags.c_contiguous:
        raise ValueError('discharge_array must be a C-order (river, time) array')
    n_rivers = network.river_ids.shape[0]
    n_steps = discharge_array.shape[1]
    sub_reach_starts, subcycles = compute_plan
    sub_reaches = np.diff(sub_reach_starts)
    if sub_reach_starts.shape != (n_rivers + 1,) or sub_reach_starts[0] != 0 or np.any(sub_reaches < 1):
        raise ValueError(
            f'sub_reach_starts must be ({n_rivers + 1},) offsets from 0 with a sub-reach per river or more'
        )
    if subcycles.shape != (n_rivers,) or np.any(subcycles < 1):
        raise ValueError(f'subcycles must be ({n_rivers},) counts of at least 1')
    n_states = int(sub_reach_starts[-1])  # one state per sub-reach
    n_out = network.original_river_ids.shape[0]  # synthetic rivers have no discharge row
    expected_shapes = {
        'channel_state': ((n_states,), channel_state.shape),
        'discharge_array': ((n_out, n_steps), discharge_array.shape),
    }
    for name, (want, got) in expected_shapes.items():
        if got != want:
            raise ValueError(f'{name} has shape {got}, expected {want} for {n_rivers} rivers and {n_steps} steps')
    # static coefficients have a value per sub-reach, and the parameters of a dynamic method one per river
    per_river = {'downstream_indices': network.downstream_indices}
    per_sub_reach = {}
    for name, value in coefficients._asdict().items():
        if hasattr(value, '_asdict'):
            per_river |= value._asdict()
        elif isinstance(value, np.ndarray):
            per_sub_reach[name] = value
    for expected, arrays in ((n_rivers, per_river), (n_states, per_sub_reach)):
        for name, array in arrays.items():
            if isinstance(array, np.ndarray) and array.shape != (expected,):
                raise ValueError(f'{name} has shape {array.shape}, expected ({expected},)')
    if forcing is not None:
        forcing.check(network.river_ids, n_steps)

    # sorted by start, the blocks tile [0, n_rivers) exactly when each ends where the next begins: no gap, no overlap
    starts = np.concatenate([blocks[0] for blocks in job_blocks])
    stops = np.concatenate([blocks[1] for blocks in job_blocks])
    order = np.argsort(starts)
    starts, stops = starts[order], stops[order]
    if np.any(stops <= starts) or starts[0] != 0 or stops[-1] != n_rivers or np.any(stops[:-1] != starts[1:]):
        raise ValueError(f'the blocks must cover every one of the {n_rivers} rivers exactly once')
    return
