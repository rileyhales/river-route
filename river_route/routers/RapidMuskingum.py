import datetime

import numpy as np
import tqdm

from .TransformRouter import TransformRouter
from ..types import FloatArray, DatetimeArray

__all__ = ['RapidMuskingum', ]


class RapidMuskingum(TransformRouter):
    """
    Muskingum router that accepts pre-aggregated catchment volumes as lateral inflow.

    Runoff volumes (m³) are divided by the runoff time step to convert to flow rates (m³/s)
    before being routed through the channel network using the Muskingum method. This matches
    the RAPID model convention where overland flow is assumed to enter the channel uniformly
    over the runoff time step.

    Required configs:
    - routing_params_file: path to routing parameters parquet file.

    Lateral input — one of:
    - catchment_runoff_files: list of paths to netCDF files containing per-catchment volumes (m³).
    - runoff_grid_files + grid_weights_file: gridded runoff depth inputs remapped to catchments.

    - discharge_files: list of output paths, one per input file.
    - dt_routing: routing sub-step in seconds (required).
    """

    @property
    def _catchment_runoff_as_volume(self) -> bool:
        return True  # True -> means as volume

    def _validate_router_configs(self) -> None:
        self._validate_lateral_runoff_configs()

    def _prepare_qlateral(self, array: FloatArray) -> FloatArray:
        """Convert catchment volumes (m³) to flow rates (m³/s) using the runoff time step."""
        array = array.astype(np.float64, copy=False)
        np.divide(array, self.dt_runoff, out=array)
        return array

    def _router(self, dates: DatetimeArray, lateral: FloatArray) -> FloatArray:
        self.logger.debug('Getting initial state arrays')
        self._read_initial_state()
        q_init = self.channel_state

        n = self.A.shape[0]
        discharge_array = np.zeros((lateral.shape[0], n), dtype=np.float64)
        q_t = q_init.astype(np.float64, copy=True)
        rhs = np.zeros(n, dtype=np.float64)
        buffer = np.zeros(n, dtype=np.float64)
        c2_r_t = np.zeros(n, dtype=np.float64)
        interval_sum = np.zeros(n, dtype=np.float64)

        t1 = datetime.datetime.now()
        runoff_iter = tqdm.tqdm(dates, desc='Runoff Volumes Routed') if self.cfg.progress_bar else dates
        if not self.cfg.progress_bar:
            self.logger.info('Performing routing computation iterations')

        for runoff_time_step, _ in enumerate(runoff_iter):
            r_t = lateral[runoff_time_step, :]
            np.multiply(self.c2, r_t, out=c2_r_t)
            interval_sum.fill(0.0)
            for _ in range(self.num_routing_steps_per_runoff):
                # rhs = c1*(A @ q_t) + c2*r_t + c3*q_t
                buffer[:] = self.A @ q_t
                np.multiply(self.c1, buffer, out=rhs)
                np.multiply(self.c3, q_t, out=buffer)
                np.add(rhs, buffer, out=rhs)
                np.add(rhs, c2_r_t, out=rhs)
                q_t[:] = self.lhs_factorized(rhs)
                interval_sum += q_t
            discharge_array[runoff_time_step, :] = interval_sum / self.num_routing_steps_per_runoff

        discharge_array[discharge_array < 0] = 0

        if self.cfg.runoff_processing_mode == 'sequential':
            self.logger.debug('Updating Channel State for Next Sequential Computation')
            self.channel_state = q_t
        elif self.cfg.runoff_processing_mode == 'ensemble':
            self.logger.debug('Recording Member State for Final State Aggregation')
            self._ensemble_member_states.append(q_t.copy())

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return discharge_array
