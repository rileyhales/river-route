import ctypes
from concurrent.futures import Executor
from dataclasses import dataclass, fields
from functools import cache
from typing import TYPE_CHECKING

import numpy as np

from ..types import FloatArray, Forcing, Int32Array
from . import numba_kernel
from .build import LIBRARY

if TYPE_CHECKING:
    from ..network.Network import Network
    from ..router.static_muskingum import Coefficients

__all__ = ['KERNELS', 'KernelNetwork', 'choose_kernel', 'kernel_arrays', 'route', 'route_region']

KERNELS = ('auto', 'numba', 'c')


@dataclass(frozen=True)
class KernelNetwork:
    """The arrays of routing.h's rr_network, in the order and dtypes both kernels read them."""

    sub_reach_starts: Int32Array  # (n_rivers + 1,) each river's first sub-reach, then the count of sub-reaches
    inflow_river_starts: Int32Array  # (n_rivers + 1,) each river's first entry in inflow_rivers, then their count
    inflow_rivers: Int32Array  # the rivers draining directly into each river, every one before it
    c1: FloatArray  # (sub-reaches,) the Muskingum coefficients of each sub-reach
    c2: FloatArray
    c3: FloatArray
    c4: FloatArray  # the coefficient of the river's runoff volume V
    weight_starts: Int32Array  # (n_rivers + 1,) each river's first weight, then the count of weights
    weight_cells: Int32Array  # the runoff cell each weight reads
    weights: FloatArray  # the runoff volume (m³) a meter of runoff depth in the cell gives the river

    def arrays(self) -> tuple[np.ndarray, ...]:
        """The arrays themselves in field order; dataclasses.astuple would copy them."""
        return tuple(getattr(self, field.name) for field in fields(self))


class _Network(ctypes.Structure):
    _fields_ = [
        ('n_rivers', ctypes.c_int32),
        ('n_sub_reaches', ctypes.c_int32),
        *((field.name, ctypes.c_void_p) for field in fields(KernelNetwork)),
    ]


class _Runoff(ctypes.Structure):
    _fields_ = [('n_cells', ctypes.c_int32), ('cell_stride', ctypes.c_ssize_t), ('depth', ctypes.c_void_p)]


class _Job(ctypes.Structure):
    _fields_ = [
        ('n_steps', ctypes.c_int32),
        ('n_blocks', ctypes.c_int32),
        ('block_starts', ctypes.c_void_p),
        ('block_stops', ctypes.c_void_p),
        ('channel_state', ctypes.c_void_p),
        ('discharge', ctypes.c_void_p),
        ('river_stride', ctypes.c_ssize_t),
    ]


def choose_kernel(kernel: str) -> str:
    """
    The kernel a route with ``kernel`` routes with: 'auto' is the C once ``rr compile`` has compiled it on this machine,
    and numba until then.

    Raises:
        ValueError: if the kernel is not one of KERNELS
        FileNotFoundError: if the C kernel is asked for and has not been compiled
    """
    if kernel not in KERNELS:
        raise ValueError(f'kernel must be one of {KERNELS}, got {kernel!r}')
    if kernel == 'auto':
        return 'c' if LIBRARY.exists() else 'numba'
    if kernel == 'c' and not LIBRARY.exists():
        raise FileNotFoundError(f'the C kernel has not been compiled at {LIBRARY}: compile it with rr compile')
    return kernel


@cache
def _library() -> ctypes.CDLL:
    library = ctypes.CDLL(str(LIBRARY))
    library.rr_describe_problem.argtypes = [ctypes.c_int]
    library.rr_describe_problem.restype = ctypes.c_char_p
    library.rr_route_rivers.argtypes, library.rr_route_rivers.restype = [ctypes.c_void_p] * 3, ctypes.c_int
    return library


def route(
    kernel: str,
    network: KernelNetwork,
    depth: FloatArray,
    n_steps: int,
    block_starts: Int32Array,
    block_stops: Int32Array,
    channel_state: FloatArray,
    discharge: FloatArray,
) -> None:
    """
    Route the rivers of the blocks, each the rivers from a start up to a stop in routing order, over the first n_steps
    runoff steps of depth, a (cells, steps) array, from the channel_state of each sub-reach, which is left at its end,
    into discharge, a (rivers, n_steps + 1) array whose first column is each river's discharge at the start. Each river
    is routed over every step before the next river.

    Raises:
        ValueError: if depth or discharge do not hold each series' steps one after another, or the kernel refuses the
            arrays, naming the problem it found
    """
    if any(array.dtype != np.float32 for array in (depth, channel_state, discharge)):
        raise ValueError('depth, channel_state, and discharge must be float32')
    # an array with no values, or one step a row, has its steps one after another whatever numpy gives as its stride
    if any(array.size and array.shape[1] > 1 and array.strides[1] != array.itemsize for array in (depth, discharge)):
        raise ValueError('depth and discharge must hold each series of steps one after another')
    if kernel == 'numba':
        problem = numba_kernel.route_rivers(
            network.arrays(), depth, n_steps, block_starts, block_stops, channel_state, discharge
        )
        description = numba_kernel.PROBLEMS.get(problem, 'an unknown problem')
    else:
        problem, description = _route_with_c(
            network, depth, n_steps, block_starts, block_stops, channel_state, discharge
        )
    if problem != numba_kernel.RR_OK:
        raise ValueError(f'the {kernel} kernel refused to route: {description}')


def _route_with_c(
    network: KernelNetwork,
    depth: FloatArray,
    n_steps: int,
    block_starts: Int32Array,
    block_stops: Int32Array,
    channel_state: FloatArray,
    discharge: FloatArray,
) -> tuple[int, str]:
    """The problem rr_route_rivers returned, and its description."""
    block_starts = np.ascontiguousarray(block_starts, dtype=np.int32)
    block_stops = np.ascontiguousarray(block_stops, dtype=np.int32)
    n_rivers = network.sub_reach_starts.shape[0] - 1
    n_sub_reaches = int(network.sub_reach_starts[n_rivers])
    network_struct = _Network(n_rivers, n_sub_reaches, *(array.ctypes.data for array in network.arrays()))
    runoff_struct = _Runoff(depth.shape[0], depth.strides[0] // depth.itemsize, depth.ctypes.data)
    job_struct = _Job(
        n_steps,
        block_starts.shape[0],
        block_starts.ctypes.data,
        block_stops.ctypes.data,
        channel_state.ctypes.data,
        discharge.ctypes.data,
        discharge.strides[0] // discharge.itemsize,
    )
    library = _library()
    structs = (ctypes.byref(struct) for struct in (network_struct, runoff_struct, job_struct))
    problem = library.rr_route_rivers(*structs)
    return problem, library.rr_describe_problem(problem).decode()


def kernel_arrays(
    network: Network, coefficients: Coefficients, forcing: Forcing, n_steps: int
) -> tuple[KernelNetwork, FloatArray]:
    """
    The arrays both kernels route one forcing file, or the channel alone, with: the network's sub-reaches with their
    coefficients and the weights of the runoff's cells, and the (cells, steps) runoff depths, which are none for
    channel routing.

    Raises:
        ValueError: if the runoff does not describe the network
    """
    if forcing is None:
        depth = np.zeros((0, n_steps), dtype=np.float32)
        weight_starts = np.zeros(network.size + 1, dtype=np.int32)
        weight_cells, weights = np.zeros(0, dtype=np.int32), np.zeros(0, dtype=np.float32)
    else:
        forcing.check(network.river_ids, n_steps)
        depth, weight_starts, weight_cells, weights = forcing.cells_and_weights()
    kernel_network = KernelNetwork(
        coefficients.sub_reach_starts,
        network.inflow_river_starts,
        network.inflow_rivers,
        coefficients.c1,
        coefficients.c2,
        coefficients.c3,
        coefficients.c4,
        weight_starts,
        weight_cells,
        weights,
    )
    return kernel_network, depth


def route_region(
    kernel: str,
    network: Network,
    coefficients: Coefficients,
    channel_state: FloatArray,
    discharge: FloatArray,
    forcing: Forcing,
    thread_pool: Executor | None,
    threads: int,
) -> None:
    """
    Route one forcing file, or the channel alone, over every river of the network: every job of sub-watershed blocks,
    concurrently on thread_pool when given, then the main stem job, which reads the discharge rows of the blocks that
    drain into it. discharge is a (river, time) array with a column for the discharge at the start and one for the end
    of each runoff step, and forcing is None for channel routing or what a Runoff generator yielded.
    """
    n_steps = discharge.shape[1] - 1
    kernel_network, depth = kernel_arrays(network, coefficients, forcing, n_steps)
    job_blocks = network.routing_blocks(threads if thread_pool is not None else 1)

    def route_job(blocks: tuple[Int32Array, ...]) -> None:
        block_starts, block_stops = blocks[0], blocks[1]
        route(kernel, kernel_network, depth, n_steps, block_starts, block_stops, channel_state, discharge)

    *sub_watershed_jobs, main_stem = job_blocks
    if thread_pool is None:
        for blocks in sub_watershed_jobs:
            route_job(blocks)
    else:
        list(thread_pool.map(route_job, sub_watershed_jobs))  # list() so a worker exception raises
    # the one barrier of the simulation: every block has finished before the main stem reads their rows
    route_job(main_stem)
