import traceback

import numpy as np
import pandas as pd
from tqdm import tqdm

from .MuskingumVlateral import MuskingumVlateral
from ._numba_kernels import dynamic_muskingum_vlateral
from ..types import FloatArray, DatetimeArray

__all__ = ['DynamicMuskingumVlateral', ]


class DynamicMuskingumVlateral(MuskingumVlateral):
    """
    Muskingum with Dynamically determined k and x with lateral inflow. K is recomputed at every routing
    substep from the most recent outflow as:

        K_i = alpha_i * max(Q_i, qmin) ** beta_i

    Muskingum coefficients (c1, c2, c3) are rebuilt from K each substep before
    the forward-substitution sweep. Required per-segment columns in params_file:
    river_id, downstream_river_id, x, alpha, beta. The 'k' column read by the
    base class is not used by this router.
    """
    alpha: FloatArray
    beta: FloatArray

    def _set_network_dependent_vectors(self) -> None:
        super()._set_network_dependent_vectors()
        try:
            df = pd.read_parquet(self.cfg.params_file, columns=['alpha', 'beta'])
        except Exception as e:
            self.logger.error(f'Error reading alpha and beta columns from params_file: {e}')
            self.logger.debug(traceback.format_exc())
            raise
        self.alpha = np.ascontiguousarray(df['alpha'].to_numpy(copy=False), dtype=np.float32)
        self.beta = np.ascontiguousarray(df['beta'].to_numpy(copy=False), dtype=np.float32)
        if np.any(self.alpha <= 0):
            raise ValueError('alpha column in params_file must be strictly positive')

    def _set_network_and_time_dependent_vectors(self, dates: DatetimeArray) -> None:
        self.logger.debug('Setting and validating time parameters')
        self.dt_runoff = self.cfg.dt_runoff or (dates[1] - dates[0]).astype('timedelta64[s]').astype(int)
        self.dt_discharge = self.cfg.dt_discharge or self.dt_runoff
        self.dt_total = self.cfg.dt_total or self.dt_runoff * dates.shape[0]
        if not self.cfg.dt_routing:
            self.logger.warning('dt_routing was not provided or is Null/False, defaulting to dt_runoff')
        self.dt_routing = self.cfg.dt_routing or self.dt_runoff

        signature = (self.dt_total, self.dt_runoff, self.dt_discharge, self.dt_routing,)
        if self._network_time_signature == signature:
            return

        if self.dt_total < self.dt_runoff:
            raise ValueError('dt_total must be >= dt_runoff')
        if self.dt_total < self.dt_discharge:
            raise ValueError('dt_total must be >= dt_discharge')
        if self.dt_discharge < self.dt_runoff:
            raise ValueError('dt_discharge must be >= dt_runoff')
        if self.dt_runoff < self.dt_routing:
            raise ValueError('dt_runoff must be >= dt_routing')
        if self.dt_total % self.dt_runoff != 0:
            raise ValueError('dt_total must be an integer multiple of dt_runoff')
        if self.dt_total % self.dt_discharge != 0:
            raise ValueError('dt_total must be an integer multiple of dt_discharge')
        if self.dt_discharge % self.dt_runoff != 0:
            raise ValueError('dt_discharge must be an integer multiple of dt_runoff')
        if self.dt_runoff % self.dt_routing != 0:
            raise ValueError('dt_runoff must be an integer multiple of dt_routing')

        self.num_runoff_steps = int(self.dt_total / self.dt_runoff)
        self.num_runoff_steps_per_discharge = int(self.dt_discharge / self.dt_runoff)
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)

        self._network_time_signature = signature
        return

    def _execute_routing(self) -> None:
        self._ensemble_member_states = []

        total_files = len(self.cfg.qlateral_files or self.cfg.grid_runoff_files)
        file_iter = self._qlateral_generator()
        if self.cfg.progress_bar:
            file_iter = tqdm(file_iter, total=total_files, desc='Files Routed')

        for dates, qlateral, runoff_file, discharge_file in file_iter:
            self.logger.info(f'Routing qlateral: {runoff_file}')
            self._set_network_and_time_dependent_vectors(dates)
            self.logger.debug('Starting routing computation')
            q_t, q_array = self._router(qlateral)
            if self.cfg.runoff_processing_mode == 'sequential':
                self.channel_state = q_t
            elif self.cfg.runoff_processing_mode == 'ensemble':
                self._ensemble_member_states.append(q_t.copy())

            if self.dt_discharge > self.dt_runoff:
                q_array = (
                    q_array
                    .reshape((
                        int(self.dt_total / self.dt_discharge),
                        int(self.dt_discharge / self.dt_runoff),
                        self.river_ids.shape[0],
                    ))
                    .mean(axis=1)
                )
                dates = dates[::self.num_runoff_steps_per_discharge]

            q_array = q_array.astype(np.float32, copy=False)
            self._write_discharges(dates, q_array, discharge_file, runoff_file)

        if self.cfg.runoff_processing_mode == 'ensemble':
            self.channel_state = np.array(self._ensemble_member_states).mean(axis=0)
        self.logger.info('-' * 60)
        return

    def _router(self, qlateral: FloatArray) -> tuple[FloatArray, FloatArray]:
        n = self.river_ids.shape[0]
        discharge_array = np.zeros((self.num_runoff_steps, n), dtype=np.float32)
        q_t = self.channel_state.astype(np.float32, copy=True)
        qlateral = np.ascontiguousarray(qlateral, dtype=np.float32)

        dynamic_muskingum_vlateral(
            q_t=q_t,
            discharge_array=discharge_array,
            downstream_indices=self.downstream_indices,
            alpha=self.alpha,
            beta=self.beta,
            x=self.x,
            dt_routing=np.float32(self.dt_routing),
            dt_runoff=np.float32(self.dt_runoff),
            n_rivers=n,
            n_steps=self.num_runoff_steps,
            n_substeps=self.num_routing_steps_per_runoff,
            vlateral=qlateral,
        )
        return q_t, discharge_array
