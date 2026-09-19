from collections.abc import Generator
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from .router import Router

PathInput = str | Path  # used at runtime for validation so it can't be a lazy type alias
type PathList = list[PathInput]
type FloatArray = NDArray[np.float32]
type IntArray = NDArray[np.int64]
type DatetimeArray = NDArray[np.datetime64]
type VlateralGenerator = Generator[tuple[DatetimeArray, FloatArray, PathInput]]


class WriteDischargesFn(Protocol):
    """A discharge writer. It may also set ``discharge_layout = 'river'``; see river_route.router.writers."""

    def __call__(
        self,
        router: Router,
        dates: DatetimeArray,
        discharge_array: FloatArray,
        discharge_file: PathInput,
        runoff_file: PathInput = '',
    ) -> None: ...


__all__ = ['PathInput', 'PathList', 'FloatArray', 'IntArray', 'DatetimeArray', 'VlateralGenerator', 'WriteDischargesFn']
