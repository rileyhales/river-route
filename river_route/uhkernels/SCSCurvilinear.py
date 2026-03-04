import numpy as np

from ._SCSBase import _SCSBase

__all__ = ['SCSCurvilinear', ]


class SCSCurvilinear(_SCSBase):
    """
    SCS curvilinear dimensionless unit hydrograph transformer.

    Standard NRCS (SCS) dimensionless unit hydrograph table.
    Source: NRCS National Engineering Handbook (NEH) Part 630, Chapter 16, Table 16-1.
    Abscissa: t/tp (dimensionless time)
    Ordinate: q/qp (dimensionless discharge)

    Parameters
    ----------
    tr   : duration of runoff generation in seconds (tr > 0)
    tc   : 1D array of time-of-concentration values in seconds, one per basin
    area : 1D array of basin areas in m², one per basin

    Notes
    -----
    Parameterization (same relationships as SCSTriangular)::

        tl = 0.6 * tc               (lag time)
        tp = tl + tr / 2            (time to peak)
        tb = 5 * tp                 (base time: t/tp table ends at 5.0)
        qp = area / (I_dim * tp)    (peak flow, m²/s; I_dim ≈ 1.333 from table integral)

    Volume conservation is exact: sum(kernel[:, j] * tr) == area[j].
    The kernel is built by evaluating the precomputed cumulative integral of the
    piecewise-linear dimensionless UH at each time-edge point and dividing by tr.
    """
    _scalar_t = np.array([
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
        1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
        2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
        4.0, 4.5, 5.0,
    ], dtype=np.float64)

    _scalar_q = np.array([
        0.000, 0.015, 0.075, 0.160, 0.280, 0.430, 0.600, 0.770, 0.890, 0.970,
        1.000, 0.980, 0.920, 0.840, 0.750, 0.660, 0.560, 0.460, 0.390, 0.330,
        0.280, 0.207, 0.147, 0.107, 0.077, 0.055, 0.040, 0.029, 0.021, 0.015,
        0.011, 0.005, 0.000,
    ], dtype=np.float64)
