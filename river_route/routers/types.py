from typing import Callable, Generator, Tuple

from ..types import ConfigDict
from ..types import DatetimeArray
from ..types import FloatArray
from ..types import IntArray
from ..types import PathInput

RunoffGeneratorSignature = Generator[Tuple[DatetimeArray, FloatArray, PathInput, PathInput], None, None]
WriteDischargesFn = Callable[[DatetimeArray, FloatArray, PathInput, PathInput], None]
FactorizedSolveFn = Callable[[FloatArray], FloatArray]

__all__ = [
    'ConfigDict',
    'PathInput',
    'FloatArray',
    'IntArray',
    'DatetimeArray',
    'RunoffGeneratorSignature',
    'WriteDischargesFn',
    'FactorizedSolveFn',
]
