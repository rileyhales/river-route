import numpy as np
from scipy.sparse import diags, eye
from scipy.sparse.linalg import factorized
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
            self.logger.info('Loading UH kernel')
            self._uh = UnitHydrograph(self.cfg.transformer_kernel_file)
            if self.cfg.transformer_state_init_file:
                self.logger.debug('Reading convolution state from parquet')
                self._uh.set_state(self.cfg.transformer_state_init_file)

        if not hasattr(self, '_hw_idx'):
            self._setup_headwater_split()

    def _setup_headwater_split(self) -> None:
        """Identify headwater streams and build reduced matrices for the inner (non-headwater) system."""
        incoming = np.asarray(self.A.sum(axis=1)).flatten()
        hw_mask = incoming == 0
        self._hw_idx = np.where(hw_mask)[0]
        self._inner_idx = np.where(~hw_mask)[0]
        self._A_inner = self.A[np.ix_(self._inner_idx, self._inner_idx)].tocsc()
        self._A_hw_to_inner = self.A[np.ix_(self._inner_idx, self._hw_idx)].tocsc()
        self.logger.info(
            f'Headwater split: {len(self._hw_idx)} headwater, {len(self._inner_idx)} inner '
            f'({len(self._hw_idx) / self.A.shape[0] * 100:.0f}% excluded from solve)'
        )

    def _set_muskingum_coefficients(self, dt_routing: float) -> None:
        super()._set_muskingum_coefficients(dt_routing)
        # Cache coefficient slices
        self._c1_inner = self.c1[self._inner_idx]
        self._c2_inner = self.c2[self._inner_idx]
        self._c3_inner = self.c3[self._inner_idx]
        # Re-factorize with the reduced inner matrix
        n_inner = len(self._inner_idx)
        lhs = eye(n_inner) - (diags(self._c1_inner) @ self._A_inner)
        self.logger.info(f'Factorizing reduced LHS matrix ({n_inner} x {n_inner})')
        self.lhs_factorized = factorized(lhs.tocsc())

    def _router(self, qlateral: FloatArray) -> tuple[FloatArray, FloatArray]:
        """Route with UH lateral superimposed on Muskingum channel routing.

        Headwater discharge is entirely the UH convolution output. Their outflow
        is injected as a known RHS contribution to the reduced inner system.
        """
        self.logger.info('Precomputing UH convolution for full timeseries')
        convolved_lateral = self._uh.convolve(qlateral)

        n = self.A.shape[0]
        hw_idx = self._hw_idx
        inner_idx = self._inner_idx
        n_inner = len(inner_idx)

        discharge_array = np.zeros((self.num_runoff_steps, n), dtype=np.float64)

        # Inner channel state from initial conditions
        q_ch_inner = self.channel_state[inner_idx].astype(np.float64, copy=True)
        q_full_inner = q_ch_inner.copy()

        rhs = np.zeros(n_inner, dtype=np.float64)
        buffer = np.zeros(n_inner, dtype=np.float64)
        interval_sum_inner = np.zeros(n_inner, dtype=np.float64)

        runoff_iter = range(self.num_runoff_steps)
        self.logger.info('Performing routing computation iterations')
        if self.cfg.progress_bar:
            runoff_iter = tqdm(runoff_iter, desc='Runoff Routed')

        for runoff_time_step in runoff_iter:
            ql_t = convolved_lateral[runoff_time_step]
            ql_hw = ql_t[hw_idx]
            ql_inner = ql_t[inner_idx]

            # Headwater discharge is purely the UH output
            discharge_array[runoff_time_step, hw_idx] = ql_hw

            # Precompute lateral upstream contribution (constant within routing substeps)
            c1_A_ql = self._c1_inner * (self._A_inner @ ql_inner + self._A_hw_to_inner @ ql_hw)

            interval_sum_inner.fill(0.0)

            for _ in range(self.num_routing_steps_per_runoff):
                # Inner RHS: c1*(A @ ql) + c2*(A @ q_full) + c3*q_channel
                rhs[:] = c1_A_ql
                buffer[:] = self._A_inner @ q_full_inner
                buffer += self._A_hw_to_inner @ ql_hw
                np.multiply(self._c2_inner, buffer, out=buffer)
                np.add(rhs, buffer, out=rhs)
                np.multiply(self._c3_inner, q_ch_inner, out=buffer)
                np.add(rhs, buffer, out=rhs)

                q_ch_inner[:] = self.lhs_factorized(rhs)
                np.add(q_ch_inner, ql_inner, out=q_full_inner)

                interval_sum_inner += q_full_inner

            discharge_array[runoff_time_step, inner_idx] = interval_sum_inner / self.num_routing_steps_per_runoff

        # Recombine final state
        q_final = np.empty(n, dtype=np.float64)
        q_final[hw_idx] = ql_hw
        q_final[inner_idx] = q_full_inner
        return q_final, discharge_array

    def _write_final_state(self) -> None:
        super()._write_final_state()
        if self.cfg.transformer_state_final_file and self._uh is not None:
            self.logger.debug('Writing final convolution state to parquet')
            self._uh.write_state(self.cfg.transformer_state_final_file)
