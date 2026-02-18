from pathlib import Path
from typing import Any
from typing import Callable

import numpy as np
from numpy.typing import NDArray

# type hints for checkers
ConfigDict = dict[str, Any]
PathInput = str | Path
FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
DatetimeArray = NDArray[np.datetime64]
WriteDischargesFn = Callable[[DatetimeArray, FloatArray, PathInput, PathInput], None]
FactorizedSolveFn = Callable[[FloatArray], FloatArray]
