from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

ConfigDict = dict[str, Any]
PathInput = str | Path
FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
DatetimeArray = NDArray[np.datetime64]

__all__ = [
    'ConfigDict',
    'PathInput',
    'FloatArray',
    'IntArray',
    'DatetimeArray',
]
