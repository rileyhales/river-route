import numpy as np

from ._numba_kernels import unit_substeps
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

    Notes on notation used in the code and solvers:
    q_full = qchannel + qlateral
    qchannel is the routed flow from muskingum routing
    qlateral is the unit hydrograph convolved runoff transformation
    _hw marks the headwater segments where only UH convolution matters
    _inner means all the non-headwater segments needing convolution and routing

    Headwater streams (no upstream connections) are excluded from the matrix solve.
    Their discharge is entirely the UH convolution output. Their outflow enters the
    reduced system as a known RHS contribution, cutting the sparse linear solve size
    roughly in half.
    """
    _ROUTER_REQUIRED_CONFIGS = ('transformer_kernel_file',)
    _uh: UnitHydrograph | None = None
    _as_volumes = False

    def _hook_before_route(self) -> None:
        if self._uh is None:
            self.logger.debug('Loading UH kernel')
            self._uh = UnitHydrograph(self.cfg.transformer_kernel_file)
            if self.cfg.transformer_state_init_file:
                self.logger.debug('Reading convolution state from parquet')
                self._uh.set_state(self.cfg.transformer_state_init_file)

        if not hasattr(self, 'hw_idx'):
            self._setup_headwater_split()

    def _setup_headwater_split(self) -> None:
        """Identify headwater streams and build reduced matrices for the inner (non-headwater) system."""
        incoming = np.asarray(self.A.sum(axis=1)).flatten()
        hw_mask = incoming == 0
        self.hw_idx = np.where(hw_mask)[0]
        self.inner_idx = np.where(~hw_mask)[0]
        self.A_inner = self.A[np.ix_(self.inner_idx, self.inner_idx)].tocsc()
        self.A_hw_to_inner = self.A[np.ix_(self.inner_idx, self.hw_idx)].tocsc()
        self.logger.info(
            f'Headwater split: {len(self.hw_idx)} headwater, {len(self.inner_idx)} inner '
            f'({len(self.hw_idx) / self.A.shape[0] * 100:.0f}% excluded from solve)'
        )

    def _set_muskingum_coefficients(self, dt_routing: float) -> None:
        super()._set_muskingum_coefficients(dt_routing)
        # Cache coefficient slices
        self._c1_inner = self.c1[self.inner_idx]
        self._c2_inner = self.c2[self.inner_idx]
        self._c3_inner = self.c3[self.inner_idx]

        self._inner_csc_indptr = self.A_inner.indptr
        self._inner_csc_indices = self.A_inner.indices
        self._inner_lhs_off_data = np.ascontiguousarray(-self._c1_inner[self.A_inner.indices])

    def _router(self, qlateral: FloatArray) -> tuple[FloatArray, FloatArray]:
        """Route with UH lateral superimposed on Muskingum channel routing."""
        self.logger.debug('Precomputing UH convolution for full timeseries')
        convolved_lateral = self._uh.convolve(qlateral)

        ql_hw = np.zeros(len(self.hw_idx), dtype=np.float64)
        ql_inner = np.zeros(len(self.inner_idx), dtype=np.float64)

        discharge_array = np.zeros((self.num_runoff_steps, self.river_ids.shape[0]), dtype=np.float64)

        # Inner channel state from initial conditions
        q_ch_inner = self.channel_state[self.inner_idx].astype(np.float64, copy=True)
        q_full_inner = q_ch_inner.copy()
        rhs = np.zeros_like(self.inner_idx, dtype=np.float64)
        interval_sum_inner = np.zeros_like(self.inner_idx, dtype=np.float64)

        runoff_iter = range(self.num_runoff_steps)
        self.logger.debug('Performing routing computation iterations')

        for runoff_time_step in runoff_iter:
            ql_t = convolved_lateral[runoff_time_step]
            ql_hw[:] = ql_t[self.hw_idx]
            ql_inner[:] = ql_t[self.inner_idx]

            discharge_array[runoff_time_step, self.hw_idx] = ql_hw  # headwater segments are entirely UH convolution
            c1_A_ql = self._c1_inner * (self.A_inner @ ql_inner + self.A_hw_to_inner @ ql_hw)
            hw_contrib = np.asarray(self.A_hw_to_inner @ ql_hw).ravel()
            interval_sum_inner.fill(0.0)
            unit_substeps(
                self._inner_csc_indptr, self._inner_csc_indices, self._inner_lhs_off_data,
                self._c2_inner, self._c3_inner,
                c1_A_ql, hw_contrib, ql_inner,
                q_ch_inner, q_full_inner,
                rhs, interval_sum_inner,
                self.num_routing_steps_per_runoff,
            )
            discharge_array[runoff_time_step, self.inner_idx] = interval_sum_inner / self.num_routing_steps_per_runoff

        # Recombine final state
        q_final = np.empty(self.river_ids.shape[0], dtype=np.float64)
        q_final[self.hw_idx] = ql_hw
        q_final[self.inner_idx] = q_full_inner
        return q_final, discharge_array

    def _write_final_state(self) -> None:
        super()._write_final_state()
        if self.cfg.transformer_state_final_file and self._uh is not None:
            self.logger.debug('Writing final convolution state to parquet')
            self._uh.write_state(self.cfg.transformer_state_final_file)
