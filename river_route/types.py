from collections.abc import Generator
from concurrent.futures import Executor
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from .forcings import CatchmentRunoffVolumes, GridCellRunoff
    from .forcings.bases import RiverGridCells
    from .router.dynamic_muskingum import DynamicCoefficients, PowerLaw
    from .router.Router import Router
    from .router.static_muskingum import StaticCoefficients

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
    'RiverForcing',
    'DynamicParameters',
    'Coefficients',
    'Transform',
    'WriteDischargesFn',
]

PathInput = str | Path  # used at runtime for validation so it can't be a lazy type alias
type PathList = list[PathInput]
type FloatArray = NDArray[np.float32]
type Float64Array = NDArray[np.float64]
type IntArray = NDArray[np.int64]
type Int32Array = NDArray[np.int32]
type DatetimeArray = NDArray[np.datetime64]
# the blocks one job routes, as (block_starts, block_stops, block_outlet, block_number); see Network.routing_blocks
type JobBlocks = tuple[Int32Array, Int32Array, Int32Array, Int32Array]
# what Runoff.generator yields for each input: its dates, its runoff in the form routing reads, and the file
type RunoffGenerator = Generator[tuple[DatetimeArray, CatchmentRunoffVolumes | GridCellRunoff, PathInput]]
# the forcing route_job is compiled for, chosen by the forcing config: None is channel routing
type Forcing = CatchmentRunoffVolumes | GridCellRunoff | None
# one river's forcing, as forcing_of_river reads it: its row of catchment runoff volumes, or its grid cells
type RiverForcing = FloatArray | RiverGridCells
# the parameters of a dynamic method, which give each river's k and x at its discharge in every routing step
type DynamicParameters = PowerLaw
# the coefficients route_job is compiled for, chosen by the parameters config
type Coefficients = StaticCoefficients | DynamicCoefficients
# the runoff transform route_job is compiled for, chosen by the transform config: None is the uniform transform, and
# any other will be a NamedTuple of its parameters
type Transform = tuple | None


class WriteDischargesFn(Protocol):
    """A discharge writer. It is handed a C-order (river, time) array, the layout the kernels route in, and the thread
    pool and thread count Router.route was given; see river_route.router.writers."""

    def __call__(
        self,
        router: Router,
        dates: DatetimeArray,
        discharge_array: FloatArray,
        discharge_file: PathInput,
        runoff_file: PathInput = '',
        *,
        thread_pool: Executor | None = None,
        threads: int = 1,
    ) -> None: ...
