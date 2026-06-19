from collections.abc import Generator
from pathlib import Path
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

PathInput = str | Path  # used at runtime for validation so it can't be a lazy type alias
type PathList = list[PathInput]
type FloatArray = NDArray[np.float32]
type IntArray = NDArray[np.int64]
type DatetimeArray = NDArray[np.datetime64]
type VlateralGeneratorSignature = Generator[tuple[DatetimeArray, FloatArray, PathInput, PathInput], None, None]


class WriteDischargesFn(Protocol):
    def __call__(
            self,
            dates: DatetimeArray,
            q_array: FloatArray,
            q_file: PathInput,
            routed_file: PathInput = '',
    ) -> None: ...


__all__ = [
    'PathInput',
    'PathList',
    'FloatArray',
    'IntArray',
    'DatetimeArray',
    'VlateralGeneratorSignature',
    'WriteDischargesFn',
]
