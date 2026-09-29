from collections.abc import Generator
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from .router.Router import Router
    from .runoff import CatchmentRunoffVolumes, GridCellRunoff

PathInput = str | Path  # used at runtime for validation so it can't be a lazy type alias
type PathList = list[PathInput]
type FloatArray = NDArray[np.float32]
type Float64Array = NDArray[np.float64]
type IntArray = NDArray[np.int64]
type Int32Array = NDArray[np.int32]
type DatetimeArray = NDArray[np.datetime64]
# what Runoff.generator yields for each input: its dates, its runoff in the form routing reads, and the file
type RunoffGenerator = Generator[tuple[DatetimeArray, CatchmentRunoffVolumes | GridCellRunoff, PathInput]]


class WriteDischargesFn(Protocol):
    """A discharge writer. It is handed a C-order (river, time) array, the layout the kernels route in; see
    river_route.router.writers."""

    def __call__(
        self,
        router: Router,
        dates: DatetimeArray,
        discharge_array: FloatArray,
        discharge_file: PathInput,
        runoff_file: PathInput = '',
    ) -> None: ...


__all__ = [
    'PathInput',
    'PathList',
    'FloatArray',
    'Float64Array',
    'IntArray',
    'Int32Array',
    'DatetimeArray',
    'RunoffGenerator',
    'WriteDischargesFn',
]
