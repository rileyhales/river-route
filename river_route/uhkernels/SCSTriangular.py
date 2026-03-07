import numpy as np

from ._SCSBase import _SCSBase

__all__ = ['SCSTriangular', ]


class SCSTriangular(_SCSBase):
    """
    SCS triangular dimensionless unit hydrograph transformer.

    Source: NRCS National Engineering Handbook (NEH) Part 630, Chapter 16, Table 16-1.

    Parameters
    ----------
    tc   : 1D array of time-of-concentration values in seconds, one per basin
    area : 1D array of basin areas in m², one per basin
    tr   : duration of runoff generation in seconds (tr > 0)

    Notes
    -----
    See the Math Derivations page in the documentation for the parameterization equations.
    """
    _scalar_t = np.array([0.0, 1.0, 2.67])
    _scalar_q = np.array([0.0, 1.0, 0.0])
