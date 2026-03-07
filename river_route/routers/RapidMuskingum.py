import numpy as np

from .TransformMuskingum import TransformMuskingum
from ._numba_kernels import rapid_route
from ..types import FloatArray

__all__ = ['RapidMuskingum', ]


class RapidMuskingum(TransformMuskingum):
    """
    Muskingum channel routing with direct lateral inflow. Lateral flow is the runoff volume divided by the runoff
    timestep — all runoff enters the channel in the interval it is generated, ignoring overland flow delay.

    See the Math Derivations page in the documentation for the full equations.
    """
    _as_volumes = True

    def _router(self, qlateral: FloatArray) -> tuple[FloatArray, FloatArray]:
        """Execute the core routing math for one runoff file and return the discharge array."""
        self.logger.debug('Getting initial state arrays')
        n = self.A.shape[0]
        discharge_array = np.zeros((self.num_runoff_steps, n), dtype=np.float64)
        q_t = self.channel_state.astype(np.float64, copy=True)
        c4_dt = self.c4 / self.dt_runoff

        rapid_route(
            self._csc_indptr, self._csc_indices, self._lhs_off_data,
            self.c2, self.c3, c4_dt, q_t, qlateral,
            discharge_array,
            self.num_routing_steps_per_runoff,
        )
        return q_t, discharge_array
