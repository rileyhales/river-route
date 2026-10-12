from collections.abc import Generator
from concurrent.futures import Executor
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from .forcings import CatchmentRunoffVolumes, GridCellRunoff
    from .router.Router import Router

__all__ = [
    'PathInput',
    'PathList',
    'FloatArray',
    'Float64Array',
    'IntArray',
    'Int32Array',
    'DatetimeArray',
    'JobBlocks',
    'RunoffGenerator',
    'Forcing',
    'WriteDischargesFn',
]

PathInput = str | Path  # used at runtime for validation so it can't be a lazy type alias
type PathList = list[PathInput]
type FloatArray = NDArray[np.float32]
type Float64Array = NDArray[np.float64]
type IntArray = NDArray[np.int64]
type Int32Array = NDArray[np.int32]
type DatetimeArray = NDArray[np.datetime64]
# the blocks one job routes, as (block_starts, block_stops), each the rivers from a start up to a stop; see
# Network.routing_blocks
type JobBlocks = tuple[Int32Array, Int32Array]
# what Runoff.generator yields for each forcing file: its dates, its runoff in the form routing reads, and the file
type RunoffGenerator = Generator[tuple[DatetimeArray, CatchmentRunoffVolumes | GridCellRunoff, PathInput]]
# the runoff the kernels route, chosen by the forcing config: None is channel routing
type Forcing = CatchmentRunoffVolumes | GridCellRunoff | None


class WriteDischargesFn(Protocol):
    """A discharge writer. It is handed a (river, time) array with each river's row contiguous, a view of the memory
    the kernels route in, and the thread pool and thread count Router.route was given; see river_route.router.writers.
    The memory is routed into again for the next file, so a writer that keeps the array past the call must copy it."""

    def __call__(
        self,
        router: Router,
        dates: DatetimeArray,
        discharge_array: FloatArray,
        discharge_file: PathInput,
        forcing_file: PathInput = '',
        *,
        thread_pool: Executor | None = None,
        threads: int = 1,
    ) -> None: ...
