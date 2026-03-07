from typing import Self

import numpy as np
import pandas as pd
import scipy.sparse
from scipy.signal import fftconvolve

from ..types import FloatArray, PathInput

__all__ = ['UnitHydrograph', ]


class UnitHydrograph:
    """
    Stateful runoff transformer driven by a precomputed kernel loaded from file.

    Two convolution modes are provided:

    - ``convolve`` — full timeseries FFT convolution for a 2D array of runoff timeseries
    - ``convolve_incrementally`` — incremental convolution for one timestep of runoff from a 1D vector

    Both modes maintain internal state so that carryover between successive calls is handled automatically.

    Parameters
    ----------
    kernel_file : path to a scipy sparse npz file

    Notes
    -----
    All times should be given in seconds; areas in m².
    See the Math Derivations page in the documentation for kernel structure and convolution details.
    """

    kernel: FloatArray  # (n_kernel_steps, n_basins)
    state: FloatArray  # (n_kernel_steps, n_basins)

    def __init__(self, kernel_file: PathInput) -> None:
        self.kernel = scipy.sparse.load_npz(kernel_file).toarray().astype(np.float64, copy=False)
        if self.kernel.ndim != 2:
            raise ValueError('kernel must be a 2D array')
        self.reset_state()

    def reset_state(self) -> None:
        self.state = np.zeros_like(self.kernel, dtype=np.float64)
        return

    def set_state(self, path: PathInput) -> Self:
        """
        Read convolution carryover state from a parquet file.

        The file should have shape (n_basins, n_kernel_steps) — basins as rows.
        """
        _state = pd.read_parquet(path).T.to_numpy(dtype=np.float64, copy=True)
        if _state.shape != self.kernel.shape:
            raise ValueError(f'state shape {_state.shape} does not match kernel shape {self.kernel.shape}')
        self.state = np.ascontiguousarray(_state)
        return self

    def write_state(self, path: PathInput) -> None:
        """Write convolution carryover state to a parquet file with shape (n_basins, n_kernel_steps)."""
        pd.DataFrame(self.state.T).to_parquet(path)
        return

    def convolve_incrementally(self, runoff_vector: FloatArray) -> FloatArray:
        """
        Incremental single-step convolution.

        Convolves one timestep of runoff with the kernel, advances the internal
        state, and returns the convolved output for this step.
        """
        self.state += self.kernel * runoff_vector
        out = self.state[0, :].copy()
        self.state[:-1, :] = self.state[1:, :]
        self.state[-1, :] = 0.0
        return out

    def convolve(self, lateral: FloatArray) -> FloatArray:
        """
        Full timeseries convolution in one pass.

        Convolves a (t, n_basins) lateral inflow array with the kernel,
        incorporates carryover state from previous calls, and updates the
        internal state with the new tail carryover.

        Parameters
        ----------
        lateral : array of shape (t, n_basins), runoff depth per timestep per basin

        Returns
        -------
        convolved : array of shape (t, n_basins), convolved lateral inflow (m³/s)
        """
        t = lateral.shape[0]
        n_ks = self.kernel.shape[0]

        # fast fourier transform to convolve along the time axis for all basins at once
        buf = fftconvolve(lateral, self.kernel, axes=0, mode='full')  # (t + n_ks - 1, n_basins)

        # Inject carryover state from previous call
        buf[:n_ks] += self.state

        # Save tail as carryover state for next call
        self.state[:] = 0
        if n_ks > 1:
            self.state[:n_ks - 1] = buf[t:]

        return buf[:t]
