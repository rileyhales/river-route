from typing import Callable

from ..types import ConfigDict
from ..types import DatetimeArray
from ..types import FloatArray
from ..types import IntArray
from ..types import PathInput

WriteDischargesFn = Callable[[DatetimeArray, FloatArray, PathInput, PathInput], None]
FactorizedSolveFn = Callable[[FloatArray], FloatArray]

__all__ = [
    'ConfigDict',
    'PathInput',
    'FloatArray',
    'IntArray',
    'DatetimeArray',
    'WriteDischargesFn',
    'FactorizedSolveFn',
]
