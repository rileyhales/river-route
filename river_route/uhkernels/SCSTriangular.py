import numpy as np

from ._SCSBase import _SCSBase

__all__ = ['SCSTriangular', ]


class SCSTriangular(_SCSBase):
    """
    SCS triangular dimensionless unit hydrograph transformer.

    Builds a kernel from basin time-of-concentration and area values using the
    SCS triangular UH equations, then delegates all state management to Transformer.

    Parameters
    ----------
    tr   : duration of runoff generation in seconds (tr > 0)
    tc   : 1D array of time-of-concentration values in seconds, one per basin
    area : 1D array of basin areas in m², one per basin

    Notes
    -----
    The triangular UH is parameterized as follows:

        tl = 0.6 * tc               (lag time)
        tp = tl + tr / 2            (time to peak)
        tb = 2.67 * tp              (base time)
        qp = 2 * area / tb          (peak flow, m²/s)

    The kernel column for each basin is the average flow (m²/s) of the unit
    hydrograph over each tr interval, obtained by integrating the piecewise-linear
    UH analytically. Volume is conserved exactly: sum(kernel[:, j] * tr) == area[j].
    """
    _scalar_t = np.array([0.0, 1.0, 2.67])
    _scalar_q = np.array([0.0, 1.0, 0.0])
