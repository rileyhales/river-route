import datetime
from typing import List

import numpy as np
import tqdm
import xarray as xr

from .AbstractTransformRouter import AbstractTransformRouter
from .types import FloatArray, DatetimeArray, PathInput

__all__ = ['RapidMuskingum', ]


class RapidMuskingum(AbstractTransformRouter):
    """
    A class for creating a RAPID style Muskingum model.

    In this model, runoff depths/volumes are assumed to uniformly discharge to the river segment during the time step
    of the runoff data. That is, for runoff data in hourly increments, the full volume is added to the inlet at
    a rate of (depth * area) / 3600 seconds over that same hour. It is added to the inflow from the upstream
    catchment(s). Overland flow is entirely ignored. The discharge is routed between catchments using the Muskingum
    method.

    Required configs:
    - routing_params_file: path to a parquet file containing routing parameters for each river segment. See docs.

    # option 1
    - lateral_volume_files: list of 1 or more paths to netcdf files containing catchment volume time series for each
        river segment at each time step. See docs.
    # option 2
    - runoff_depth_grids: list of 1 or more paths to netcdf files containing runoff depth time series for each river
        segment at each time step. Must also have a corresponding grid_weights_file. See docs.
    - grid_weights_file: path to a netCDF file containing the weights for converting runoff depths to volumes for each
        river segment. Must accompany runoff_depth_grids. See docs.

    - discharge_files: list of 1 or more paths to netCDF files where the routed discharge time series will be written
        for each river segment. Must match the number of input files. See docs.
    """
    # State variables
    _ensemble_member_states: List[FloatArray]  # for ensemble routing, stores member states for computing final state

    # Time options
    num_routing_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_discharge: int

    def _validate_router_configs(self) -> None:
        self._validate_lateral_runoff_configs()

    def _read_lateral_file(self, file_path: PathInput):
        self.logger.info(f'Reading catchment volume file: {file_path}')
        with xr.open_dataset(file_path) as ds:
            return (ds['time'].values.astype('datetime64[s]'),
                    ds[self.conf['var_catchment_volume']].values)

    def _prepare_lateral_array(self, array: FloatArray) -> FloatArray:
        """Convert catchment volumes to flow rates (m³ -> m³/s) using the routing timestep."""
        array = array.astype(np.float64, copy=False)
        np.divide(array, self.dt_runoff, out=array)
        return array

    def _router(self, dates: DatetimeArray, volumes: FloatArray) -> FloatArray:
        self.logger.debug('Getting initial state arrays')
        self._read_initial_state()
        q_init = self.channel_state
        self._ensemble_member_states = []

        # declare arrays for routing computations once to avoid repeated and duplicate allocations in the loop
        discharge_array = np.zeros((volumes.shape[0], self.A.shape[0]))
        q_t = np.empty_like(q_init, dtype=np.float64)
        r_prev = np.zeros_like(q_init, dtype=np.float64)
        rhs = np.zeros(self.A.shape[0], dtype=np.float64)
        buffer = np.zeros(self.A.shape[0], dtype=np.float64)
        c2_r_t = np.zeros(self.A.shape[0], dtype=np.float64)
        interval_sum = np.zeros(self.A.shape[0], dtype=np.float64)
        q_t[:] = q_init
        # r_prev[:] = r_init  todo: correct math to no longer use r_prev

        t1 = datetime.datetime.now()
        if self.conf['progress_bar']:
            dates = tqdm.tqdm(dates, desc='Runoff Volumes Routed')
        else:
            self.logger.info('Performing routing computation iterations')
        for runoff_time_step, runoff_end_date in enumerate(dates):
            r_t = volumes[runoff_time_step, :]
            np.multiply(self.c2, r_t, out=c2_r_t)
            interval_sum.fill(0.0)  # add then divide to avoid accumulating large array and averaging
            for routing_sub_iteration_num in range(self.num_routing_steps_per_runoff):
                # rhs = (self.c1 * ((self.A @ q_t) + r_prev)) + (self.c2 * r_t) + (self.c3 * q_t)
                buffer[:] = self.A @ q_t
                np.add(buffer, r_prev, out=buffer)
                np.multiply(self.c1, buffer, out=rhs)
                np.multiply(self.c3, q_t, out=buffer)
                np.add(rhs, buffer, out=rhs)
                np.add(rhs, c2_r_t, out=rhs)
                # solve Ax=b, LHS q_t = rhs
                q_t[:] = self.lhs_factorized(rhs)
                # add to interval accumulator for averaging at the end of the runoff time step
                interval_sum += q_t
            discharge_array[runoff_time_step, :] = interval_sum / self.num_routing_steps_per_runoff
            r_prev[:] = r_t

        # Enforce positive flows in case of negative numerical results
        discharge_array[discharge_array < 0] = 0

        # if simulation type is ensemble, then do not overwrite the initial state
        if self.conf['input_type'] == 'sequential':
            self.logger.debug('Updating Channel State for Next Sequential Computation')
            self.channel_state = q_t
        if self.conf['input_type'] == 'ensemble':
            self.logger.debug('Recording Member State for Final State Aggregation')
            self._ensemble_member_states.append(q_t.copy())

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return discharge_array
