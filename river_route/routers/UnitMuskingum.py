import numpy as np

from .TransformMuskingum import TransformMuskingum
from ._numba_kernels import unit_route
from ..types import FloatArray
from ..uhkernels import UnitHydrograph

__all__ = ['UnitMuskingum', ]


class UnitMuskingum(TransformMuskingum):
    """
    Muskingum channel routing with unit hydrograph lateral inflow. Lateral inflow at each segment is determined by
    convolving runoff depths (m) with a precomputed unit hydrograph kernel. Discharge is the superposition of overland
    flow by unit hydrograph and channel routing merged at each time step.

    Notation used in the code and solvers:

    - q_full = q_channel + q_lateral
    - q_channel is the routed flow from Muskingum routing
    - q_lateral is the unit hydrograph convolved runoff transformation
    - _hw marks headwater segments (UH convolution only, excluded from matrix solve)
    - _inner marks non-headwater segments needing both convolution and routing

    See the Math Derivations page in the documentation for the full equations.
    """
    _ROUTER_REQUIRED_CONFIGS = ('uh_kernel_file',)
    _uh: UnitHydrograph | None = None
    _as_volumes = False

    def _hook_before_route(self) -> None:
        if self._uh is None:
            self.logger.debug('Loading UH kernel')
            self._uh = UnitHydrograph(self.cfg.uh_kernel_file)
            if self.cfg.uh_state_init_file:
                self.logger.debug('Reading convolution state from parquet')
                self._uh.set_state(self.cfg.uh_state_init_file)

        # prepare to split arrays for headwater vs inner segments
        if not hasattr(self, 'hw_idx'):
            incoming = np.asarray(self.A.sum(axis=1)).flatten()
            hw_mask = incoming == 0
            self.hw_idx = np.where(hw_mask)[0]
            self.inner_idx = np.where(~hw_mask)[0]
            self.A_inner = self.A[np.ix_(self.inner_idx, self.inner_idx)].tocsc()
            self.A_hw_to_inner = self.A[np.ix_(self.inner_idx, self.hw_idx)].tocsc()

            # Store CSC arrays for numba kernel
            self._a_inner_indptr = self.A_inner.indptr
            self._a_inner_indices = self.A_inner.indices
            self._a_inner_data = np.ascontiguousarray(self.A_inner.data.astype(np.float64, copy=False))
            self._a_hw_indptr = self.A_hw_to_inner.indptr
            self._a_hw_indices = self.A_hw_to_inner.indices
            self._a_hw_data = np.ascontiguousarray(self.A_hw_to_inner.data.astype(np.float64, copy=False))

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
        convolved_lateral = np.ascontiguousarray(self._uh.convolve(qlateral), dtype=np.float64)

        discharge_array = np.zeros((self.num_runoff_steps, self.river_ids.shape[0]), dtype=np.float64)
        q_ch_inner = self.channel_state[self.inner_idx].astype(np.float64, copy=True)
        q_full_inner = q_ch_inner.copy()

        self.logger.debug('Performing routing computation iterations')
        unit_route(
            self._inner_csc_indptr, self._inner_csc_indices, self._inner_lhs_off_data,
            self._a_inner_indptr, self._a_inner_indices, self._a_inner_data,
            self._a_hw_indptr, self._a_hw_indices, self._a_hw_data,
            self._c1_inner, self._c2_inner, self._c3_inner,
            self.hw_idx, self.inner_idx,
            q_ch_inner, q_full_inner,
            convolved_lateral,
            discharge_array,
            self.num_routing_steps_per_runoff,
        )

        # Recombine final state
        q_final = np.empty(self.river_ids.shape[0], dtype=np.float64)
        q_final[self.hw_idx] = convolved_lateral[-1][self.hw_idx]
        q_final[self.inner_idx] = q_full_inner
        return q_final, discharge_array

    def _write_final_state(self) -> None:
        super()._write_final_state()
        if self.cfg.uh_state_final_file and self._uh is not None:
            self.logger.debug('Writing final convolution state to parquet')
            self._uh.write_state(self.cfg.uh_state_final_file)
