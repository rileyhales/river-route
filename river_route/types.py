from collections.abc import Generator, Callable
from pathlib import Path
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

PathInput = str | Path
PathTypes = (str, Path)  # for isinstance checks rather than type hinting
PathList = list[PathInput]
FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
DatetimeArray = NDArray[np.datetime64]
QlateralGeneratorSignature = Generator[tuple[DatetimeArray, FloatArray, FloatArray, PathInput, PathInput], None, None]


class WriteDischargesFn(Protocol):
    def __call__(
            self,
            dates: DatetimeArray,
            q_array: FloatArray,
            q_file: PathInput,
            routed_file: PathInput = '',
    ) -> None: ...


FactorizedSolveFn = Callable[[FloatArray], FloatArray]

__all__ = [
    'PathInput',
    'PathTypes',
    'PathList',
    'FloatArray',
    'IntArray',
    'DatetimeArray',
    'QlateralGeneratorSignature',
    'WriteDischargesFn',
    'FactorizedSolveFn',
]
