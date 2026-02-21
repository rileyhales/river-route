from typing import Callable

from ..typing import ConfigDict
from ..typing import DatetimeArray
from ..typing import FloatArray
from ..typing import IntArray
from ..typing import PathInput

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
