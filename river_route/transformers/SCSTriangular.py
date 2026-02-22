import numpy as np

from .AbstractTransformer import AbstractBaseTransformer
from ..types import FloatArray

__all__ = ['SCSUnitHydrograph', ]


class SCSUnitHydrograph(AbstractBaseTransformer):
    """
    SCS triangular dimensionless unit hydrograph transformer.

    This method needs to be modified to implement the following logic

    tp = time to peak
    tl = lag time, time from center of mass of runoff to peak.
    tr = time of runoff (same as dt. runoff is uniform over dt in this implementation)
    tc = time of concentration (provided)
    tb = time to base, the duration of the unit hydrograph
    qp = peak flow occurring at tp
    area = area of the basin in *meters squared*
    r = vector of runoff depths in meters

    This class builds its kernel and unit hydrograph base in the following way:
    1. accept a vector of tc values, vector of area in meters squared, and one time step (dt)
    2. calculate tl = 0.6 * tc
    3. calculate tp = tl + dt/2
    4. calculate tb = 2.67 * tp
    5. calculate qp = 2 * area / tb
    5. determine the number of kernel time steps needed as ceiling of (tb / dt) for each basin
    6. calculate the unit hydrograph kernel discretized to the average over each dt.
    7. store the kernel and a sparse empty state array of the same size.

    A notes on averaging the UH over dt. we're using the triangular unit hydrograph so we only need tb and qp.

    Initialization
    --------------
    Compute mode  __init__(tc, area, dt):
        Builds the kernel from the SCS triangular unit hydrograph equations above.

    Load mode  AbstractBaseTransformer.from_cached(dt, kernel, state=None):
        Inherited class method.  Accepts a pre-computed kernel (and optional state) in tall format
        (n_basins, n_time_steps) as arrays or parquet paths, skipping the SCS computation.
    """

    tc: FloatArray
    area: FloatArray

    def __init__(self, dt: float, tc: FloatArray, area: FloatArray) -> None:
        # Step 1: accept a vector of tc values, vector of area in meters squared, and one time step (dt)
        self.tc = np.asarray(tc, dtype=np.float64)
        self.area = np.asarray(area, dtype=np.float64)
        if self.tc.ndim != 1:
            raise ValueError('tc must be a 1D float array')
        if self.area.ndim != 1:
            raise ValueError('area must be a 1D float array')
        if self.tc.shape != self.area.shape:
            raise ValueError('tc and area must have the same length')
        super().__init__(dt=dt)  # validates dt > 0, calls precompute() -> _build_kernel() then reset_state()

    def _build_kernel(self) -> FloatArray:
        tl = 0.6 * self.tc
        tp = tl + self.dt / 2.0
        tb = 2.67 * tp
        qp = 2.0 * self.area / tb

        # Step 5b: number of kernel time steps per basin = ceiling(tb / dt)
        n_steps = np.ceil(tb / self.dt).astype(int)
        max_steps = int(n_steps.max())

        # Step 6: build a (max_steps+1, n_basins) matrix of time-edge values.
        # Row i holds the left edge of interval i; the extra row is the right edge of the last interval.
        # Clamping to tb[j] makes t_edges[i,j] == tb[j] for all i >= n_steps[j],
        # so np.diff produces zero for those padding rows automatically.
        i_row = np.arange(max_steps + 1, dtype=np.float64)[:, np.newaxis]  # (max_steps+1, 1)
        t_edges = np.minimum(i_row * self.dt, tb)  # (max_steps+1, n_basins)

        # Cumulative integral u_cum(t) = ∫₀ᵗ u(s) ds of the triangular UH, vectorized over all basins:
        #   rising limb  (0 ≤ t ≤ tp): u(t) = qp * t / tp
        #     → u_cum(t) = qp * t² / (2 * tp)
        #   falling limb (tp < t ≤ tb): u(t) = qp * (tb - t) / (tb - tp)
        #     → u_cum(t) = qp*tp/2 + qp/(tb-tp) * (tb*(t-tp) - (t²-tp²)/2)
        # tp, tb, qp are shape (n_basins,) and broadcast across rows automatically
        rising = qp * t_edges ** 2 / (2.0 * tp)
        falling = qp * tp / 2.0 + qp / (tb - tp) * (tb * (t_edges - tp) - (t_edges ** 2 - tp ** 2) / 2.0)
        u_cum = np.where(t_edges <= tp, rising, falling)  # (max_steps+1, n_basins)

        # Average flow over each interval = (u_cum[i+1] - u_cum[i]) / dt
        # This satisfies the AbstractBaseTransformer discretization contract (see base class docstring).
        # Step 7: zero padding from the clamp means shorter-basin columns are already zero-padded.
        return np.diff(u_cum, axis=0) / self.dt  # (max_steps, n_basins)
