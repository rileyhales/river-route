import numpy as np
from tqdm import tqdm

from .TransformMuskingum import TransformMuskingum
from ..types import FloatArray
from ..uhkernels import UnitHydrograph

__all__ = ['UnitMuskingum', ]


class UnitMuskingum(TransformMuskingum):
    """
    Muskingum channel routing with unit hydrograph lateral inflow. Lateral inflow at each segment is determined by
    convolving runoff depths (m) with a precomputed unit hydrograph kernel. Discharge is the superposition of overland
    flow by unit hydrograph and channel routing which are merged together at each time step.

    c1 = (dt/k - 2x) / (dt/k + 2(1-x))
    c2 = (dt/k + 2x) / (dt/k + 2(1-x))
    c3 = (2(1-x) - dt/k) / (dt/k + 2(1-x))
    Ql_t = unit hydrograph convolution output at time t (m³/s)

    In matrix form, the router needs to solve the equations:
    (I - c1 * A) @ Q_channel_t+1 = c1 * (A @ Ql_t+1) + c2 * (A @ Q_full_t) + c3 * Q_channel_t
    Q_full_t+1 = Q_channel_t+1 + Ql_t+1

    The channel and full discharge are tracked separately so that headwater segments report exactly
    the UH output without Muskingum amplification, while downstream segments receive correctly
    routed contributions from all upstream sources.
    """
    _ROUTER_REQUIRED_CONFIGS = ('transformer_kernel_file',)
    _uh: UnitHydrograph | None = None
    _as_volumes = False

    def _hook_before_route(self) -> None:
        if self._uh is not None:
            return
        self.logger.info('Loading UH kernel')
        self._uh = UnitHydrograph(self.cfg.transformer_kernel_file)
        if self.cfg.transformer_state_init_file:
            self.logger.debug('Reading convolution state from parquet')
            self._uh.set_state(self.cfg.transformer_state_init_file)

    def _router(self, qlateral: FloatArray) -> tuple[FloatArray, FloatArray]:
        """Route with UH lateral superimposed on Muskingum channel routing.

        Two state vectors are maintained:
        - q_channel: flow from Muskingum channel routing only (no lateral accumulated via c3)
        - q_full: q_channel + ql — the reported discharge and what flows downstream

        The routing equation is:
            rhs = c1*(A @ ql_t) + c2*(A @ q_full) + c3*q_channel
            q_channel = (I - c1*A)^{-1} @ rhs
            q_full = q_channel + ql_t

        This ensures headwater discharge equals the UH output exactly (no Muskingum amplification)
        while correctly routing all flow downstream through the channel network.
        """
        self.logger.info('Precomputing UH convolution for full timeseries')
        convolved_lateral = self._uh.convolve(qlateral)

        n = self.A.shape[0]
        discharge_array = np.zeros((self.num_runoff_steps, n), dtype=np.float64)
        q_full = self.channel_state.astype(np.float64, copy=True)
        q_channel = q_full.copy()
        rhs = np.zeros(n, dtype=np.float64)
        buffer = np.zeros(n, dtype=np.float64)
        interval_sum = np.zeros(n, dtype=np.float64)

        runoff_iter = range(self.num_runoff_steps)
        self.logger.info('Performing routing computation iterations')
        if self.cfg.progress_bar:
            runoff_iter = tqdm(runoff_iter, desc='Runoff Routed')

        for runoff_time_step in runoff_iter:
            ql_t = convolved_lateral[runoff_time_step]
            interval_sum.fill(0.0)
            for _ in range(self.num_routing_steps_per_runoff):
                buffer[:] = self.A @ ql_t
                np.multiply(self.c1, buffer, out=rhs)
                buffer[:] = self.A @ q_full
                np.multiply(self.c2, buffer, out=buffer)
                np.add(rhs, buffer, out=rhs)
                np.multiply(self.c3, q_channel, out=buffer)
                np.add(rhs, buffer, out=rhs)
                q_channel[:] = self.lhs_factorized(rhs)
                np.add(q_channel, ql_t, out=q_full)
                interval_sum += q_full
            discharge_array[runoff_time_step, :] = interval_sum / self.num_routing_steps_per_runoff

        return q_full, discharge_array

    def _write_final_state(self) -> None:
        super()._write_final_state()
        if self.cfg.transformer_state_final_file and self._uh is not None:
            self.logger.debug('Writing final convolution state to parquet')
            self._uh.write_state(self.cfg.transformer_state_final_file)
