import warnings
from abc import ABC, abstractmethod
from typing import Self

import numpy as np
import pandas as pd

from ..typing import FloatArray, PathInput

__all__ = ['AbstractBaseTransformer', ]


class AbstractBaseTransformer(ABC):
    """
    Abstract base class for stateful runoff transformers used by UnitMuskingum.

    Implements the shared lifecycle for any transformer to integrate with UnitMuskingum, including:
    - Initialization with time step (dt)
    - Precomputation of the kernel matrix of shape (n_time_steps, n_basins) 1 row per time step, 1 column per basin
    - State management of the unit hydrographs for each basin (the "kernel")
    - a transform method that applies the kernel to a runoff vector and advances the state by one time step.

    Subclasses only need to implement kernel generation in `_build_kernel`. that method should do the following:
    - create a 2D array of shape (n_time_steps, n_basins)
    - each column should constitute the unit hydrograph for the basin in that column (R=1)
    - the row value should be the average flow of the unit hydrograph discretized over the specified time step (dt)

    Explanation of unit hydrograph discretization:
    the unit hydrograph is a function that describes how runoff from a single time step of rainfall (R=1) translates
    into flow at the outlet over time. When we discretize this function into time steps of size dt, we need to compute
    the average flow that would occur in each time step due to that initial runoff. This means integrating the unit
    hydrograph function over each time step interval and dividing by dt to get the average flow. For example, if the
    unit hydrograph is defined as a continuous function u(t), then the value in row i of column j should equal the
    integral of u(t) from t=i*dt to t=(i+1)*dt, divided by dt to get the average flow over that time step.

    Notes:
    1. The kernel may be very sparse if there is a big difference between the Tc of the longest and shortest basins.
       The kernel will have the shape of n_time_steps of the longest basin which is close to Tc / dt
    2. Times should all always be given in seconds.
    """

    dt: float
    kernel: FloatArray
    state: FloatArray

    def __init__(self, dt: float) -> None:
        if dt <= 0:
            raise ValueError('dt must be > 0')
        self.dt = float(dt)
        self.precompute()

    def precompute(self) -> None:
        self.kernel = self._build_kernel()
        if self.kernel.ndim != 2:
            raise ValueError('Transformer kernel must be a 2D array')
        self.reset_state()

    def reset_state(self) -> None:
        self.state = np.zeros_like(self.kernel, dtype=np.float64)

    @abstractmethod
    def _build_kernel(self) -> FloatArray:
        ...

    def validate_runoff_vector(self, runoff_vector: FloatArray) -> FloatArray:
        runoff_vector = np.asarray(runoff_vector, dtype=np.float64)
        if runoff_vector.ndim != 1:
            raise ValueError('runoff_vector must be a 1D array')
        if runoff_vector.shape[0] != self.kernel.shape[1]:
            raise ValueError(
                f'runoff_vector length {runoff_vector.shape[0]} does not match number of basins '
                f'{self.kernel.shape[1]}'
            )
        if not np.all(np.isfinite(runoff_vector)):
            raise ValueError('runoff_vector contains non-finite values')
        return runoff_vector

    def save_kernel(self, filename: PathInput | None = None) -> Self:
        """
        Save the kernel in tall orientation so that it reads and writes as parquet better (n_basins, n_time_steps)
        The file name is {class_name}_kernel.parquet in the current working directory.
        """
        warnings.warn('it is strongly recommended that you label the time step of the kernel matrix')
        pd.DataFrame(self.kernel.T).to_parquet(filename)
        return self

    @classmethod
    def from_kernel(cls, dt: float, kernel: PathInput) -> Self:
        """
        Load a pre-computed kernel from a parquet file instead of running a subclass __init__.

        Parameters
        ----------
        dt     : time step in seconds.
        kernel : path to a parquet file of shape (n_basins, n_time_steps) — tall format.

        The kernel is transposed from tall (n_basins, n_time_steps) to the internal
        (n_time_steps, n_basins) layout.  Call set_state() after to warm-start the state.
        """
        instance = cls.__new__(cls)
        if float(dt) <= 0:
            raise ValueError('dt must be > 0')
        instance.dt = float(dt)
        instance.kernel = pd.read_parquet(kernel).T.to_numpy(dtype=np.float64)
        if instance.kernel.ndim != 2:
            raise ValueError('kernel must be a 2D array')
        instance.reset_state()
        return instance

    def set_state(self, state: PathInput) -> Self:
        """
        Ingest a state from a parquet file, replacing the current state.

        Parameters
        ----------
        state : path to a parquet file of shape (n_basins, n_time_steps) — tall format.

        The state is transposed from tall (n_basins, n_time_steps) to the internal
        (n_time_steps, n_basins) layout and must match the kernel shape.
        """
        _state = pd.read_parquet(state).T.to_numpy(dtype=np.float64)
        if _state.shape != self.kernel.shape:
            raise ValueError(
                f'state shape {_state.shape} does not match kernel shape {self.kernel.shape}'
            )
        self.state = _state
        return self

    def transform(self, runoff_vector: FloatArray) -> FloatArray:
        runoff_vector = self.validate_runoff_vector(runoff_vector)
        self.state += self.kernel * runoff_vector
        out = self.state[0, :].copy()  # the vector returned is the 1st time step (row)
        self.state[:-1, :] = self.state[1:, :]  # then slide up all the values 1 time step (up the rows)
        self.state[-1, :] = 0.0  # and clear the last time step (row) which is now empty after the slide
        return out
