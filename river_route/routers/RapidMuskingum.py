import numpy as np

from .TransformMuskingum import TransformMuskingum
from ._numba_kernels import rapid_substeps
from ..types import FloatArray

__all__ = ['RapidMuskingum', ]


class RapidMuskingum(TransformMuskingum):
    """
    Muskingum channel routing with a direct lateral inflow component, Ql. Lateral flow is the runoff volume divided
    by the runoff timestep. All runoff enters the channel in the interval it is generated which ignores the overland
    flow time, or in other words the runoff transformation.

    Q_t+1 = (c1 * I_t+1) + (c2 * I_t) + (c3 * Q_t) + (c4 * Ql_t)
    c1 = (dt/k - 2x) / (dt/k + 2(1-x))
    c2 = (dt/k + 2x) / (dt/k + 2(1-x))
    c3 = (2(1-x) - dt/k) / (dt/k + 2(1-x))
    c4 = c1 + c2

    In matrix form, the router needs to solve the equation:
    (I - c1 * A) @ Q_t+1 = c2 * (A @ Q_t) + c3 * Q_t + c4 * Ql_t
    """
    _as_volumes = True

    def transform_runoff(self, r_t: FloatArray) -> FloatArray:
        """Convert a runoff volume (m³) to the lateral term for the routing equation."""
        return r_t / self.dt_runoff * self.c4

    def _router(self, qlateral: FloatArray) -> tuple[FloatArray, FloatArray]:
        """Execute the core routing math for one runoff file and return the discharge array."""
        self.logger.debug('Getting initial state arrays')
        q_init = self.channel_state

        n = self.A.shape[0]
        discharge_array = np.zeros((self.num_runoff_steps, n), dtype=np.float64)
        q_t = q_init.astype(np.float64, copy=True)
        rhs = np.zeros(n, dtype=np.float64)
        interval_sum = np.zeros(n, dtype=np.float64)

        for runoff_time_step in range(self.num_runoff_steps):
            ql_t = self.transform_runoff(qlateral[runoff_time_step, :])
            interval_sum.fill(0.0)
            rapid_substeps(
                self._csc_indptr, self._csc_indices, self._lhs_off_data,
                self.c2, self.c3, q_t, ql_t,
                rhs, interval_sum,
                self.num_routing_steps_per_runoff,
            )
            discharge_array[runoff_time_step, :] = interval_sum / self.num_routing_steps_per_runoff
        return q_t, discharge_array
