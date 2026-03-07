from abc import ABC

import numpy as np
import scipy.sparse

from ..types import FloatArray, PathInput

__all__ = ['_SCSBase', ]


class _SCSBase(ABC):
    """
    Shared base for SCS unit hydrograph kernels.

    Subclasses define ``_scalar_t`` and ``_scalar_q`` class attributes representing the dimensionless UH shape.
    The base class integrates the table, scales to physical units, and builds the kernel.

    See the Math Derivations page in the documentation for the SCS parameterization equations.
    """

    _scalar_t: FloatArray
    _scalar_q: FloatArray

    kernel: FloatArray  # (n_kernel_steps, n_basins)

    tr: float  # duration of runoff, equal to dt when computing on discrete time steps
    tc: FloatArray
    area: FloatArray

    # derived parameters
    tl: FloatArray  # lag time
    tp: FloatArray  # time to peak, also called time of rise
    tb: FloatArray  # time to base (duration of the UH)
    qp: FloatArray  # peak flow (m²/s)

    def __init__(self, *, tc: FloatArray, area: FloatArray, tr: float, ) -> None:
        if float(tr) <= 0:
            raise ValueError('tr must be > 0')
        tc = np.asarray(tc, dtype=np.float64)
        area = np.asarray(area, dtype=np.float64)
        if tc.ndim != 1:
            raise ValueError('tc must be a 1D float array')
        if area.ndim != 1:
            raise ValueError('area must be a 1D float array')
        if tc.shape != area.shape:
            raise ValueError('tc and area must have the same length')
        self.tr = float(tr)
        self.tc = tc
        self.area = area

        dimensionless_integral = np.concatenate([
            [0.0],
            np.cumsum((self._scalar_q[:-1] + self._scalar_q[1:]) / 2.0 * np.diff(self._scalar_t)),
        ])
        dimensionless_uh_area = dimensionless_integral[-1]

        self.tl = 0.6 * self.tc
        self.tp = self.tl + self.tr / 2.0
        self.tb = self._scalar_t[-1] * self.tp
        self.qp = self.area / (dimensionless_uh_area * self.tp)

        n_steps = np.ceil(self.tb / self.tr).astype(int)
        max_steps = int(n_steps.max())

        i_row = np.arange(max_steps + 1, dtype=np.float64)[:, np.newaxis]
        t_edges = np.minimum(i_row * self.tr, self.tb)  # (max_steps+1, n_basins)

        t_over_tp = t_edges / self.tp
        i_cum_dim = np.interp(
            t_over_tp.ravel(), self._scalar_t, dimensionless_integral, right=dimensionless_uh_area
        ).reshape(t_over_tp.shape)

        u_cum = i_cum_dim * (self.qp * self.tp)

        # Average flow over each interval = Δu_cum / tr
        self.kernel = np.diff(u_cum, axis=0) / self.tr  # (max_steps, n_basins)

    def save(self, path: PathInput) -> None:
        """Save the kernel as a scipy sparse npz file."""
        scipy.sparse.save_npz(path, scipy.sparse.csr_matrix(self.kernel))
        return
