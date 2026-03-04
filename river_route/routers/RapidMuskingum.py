import numpy as np
from tqdm import tqdm

from .TransformMuskingum import TransformMuskingum
from ..types import FloatArray

__all__ = ['RapidMuskingum', ]


class RapidMuskingum(TransformMuskingum):
    """
    Solves the Muskingum Cunge routing equation for channel routing with a lateral inflow component, Ql.

    Q_t+1 = (c1 * I_t+1) + (c2 * I_t) + (c3 * Q_t) + (c4 * Ql_t)
    c1 = (dt/k - 2x) / (dt/k + 2(1-x))
    c2 = (dt/k + 2x) / (dt/k + 2(1-x))
    c3 = (2(1-x) - dt/k) / (dt/k + 2(1-x))
    c4 = dt/k / (dt/k + 2(1-x)) = c1 + c2
    Ql_t is the lateral inflow volume generated at time t divided by dt, the time step of the runoff file.

    (I - c1 @ A) Q_t+1 = c2 * (A @ q_t) + c3 * q_t + c4 * Ql_t

    Because of the definition of Ql_t used here, all runoff is assumed to make it into the channel and is routed
    during the interval it is generated. This assumption is best for long time steps, such as daily averages, where the
    runoff distribution matters less compared to the time step and travel times. It can have the affect of speeding up
    and amplifying the peak discharges compared to using an overland flow runoff transformation.
    """

    @property
    def _catchment_runoff_as_volume(self) -> bool:
        return True  # True -> means as volume

    def transform_runoff(self, r_t: FloatArray) -> FloatArray:
        """Convert a runoff volume (m³) to the lateral term for the routing equation."""
        return r_t / self.dt_runoff * self.c4

    def _router(self, lateral: FloatArray) -> tuple[FloatArray, FloatArray]:
        """Execute the core routing math for one runoff file and return the discharge array."""
        self.logger.debug('Getting initial state arrays')
        q_init = self.channel_state

        n = self.A.shape[0]
        discharge_array = np.zeros((self.num_runoff_steps, n), dtype=np.float64)
        q_t = q_init.astype(np.float64, copy=True)
        rhs = np.zeros(n, dtype=np.float64)
        buffer = np.zeros(n, dtype=np.float64)
        interval_sum = np.zeros(n, dtype=np.float64)

        runoff_iter = range(self.num_runoff_steps)  # default to no progress bar
        self.logger.info('Performing routing computation iterations')
        if self.cfg.progress_bar:
            runoff_iter = tqdm(runoff_iter, desc='Runoff Routed')

        for runoff_time_step in runoff_iter:
            ql_t = self.transform_runoff(lateral[runoff_time_step, :])
            interval_sum.fill(0.0)
            for _ in range(self.num_routing_steps_per_runoff):
                # rhs = c2*(A @ q_t) + c3*q_t + c4*ql_t
                buffer[:] = self.A @ q_t
                np.multiply(self.c2, buffer, out=rhs)
                np.multiply(self.c3, q_t, out=buffer)
                np.add(rhs, buffer, out=rhs)
                np.add(rhs, ql_t, out=rhs)
                q_t[:] = self.lhs_factorized(rhs)
                interval_sum += q_t
            discharge_array[runoff_time_step, :] = interval_sum / self.num_routing_steps_per_runoff
        return q_t, discharge_array
