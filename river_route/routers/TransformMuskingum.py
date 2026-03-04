from abc import ABC, abstractmethod

import numpy as np
import xarray as xr

from .Muskingum import Muskingum
from ..types import DatetimeArray, FloatArray, RunoffGeneratorSignature

__all__ = ['TransformMuskingum', ]

from ..runoff import grid_to_catchment


class TransformMuskingum(Muskingum, ABC):
    """
    Intermediate abstract router class adding routing methods that require pre-processing of the lateral inflow
    """
    _ROUTER_REQUIRED_CONFIGS = ()

    # State variables
    _ensemble_member_states: list[FloatArray]  # for ensemble routing

    # Time options
    num_runoff_steps: int
    num_routing_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_discharge: int

    @property
    @abstractmethod
    def _catchment_runoff_as_volume(self) -> bool:
        """Whether the catchment runoff generator yields volumes (m³) or depths (m)"""

    def _catchment_runoff_generator(self) -> RunoffGeneratorSignature:
        if self.cfg.catchment_runoff_files:
            for lateral_file, discharge_file in zip(self.cfg.catchment_runoff_files, self.cfg.discharge_files):
                self.logger.info('-' * 60)
                with xr.open_dataset(lateral_file) as ds:
                    dates = ds['time'].values.astype('datetime64[s]')
                    array = ds[self.cfg.var_catchment_runoff_variable].values.astype(np.float64, copy=False)
                    yield dates, array, lateral_file, discharge_file
        elif self.cfg.runoff_grid_files and self.cfg.grid_weights_file:
            for runoff_file, discharge_file in zip(self.cfg.runoff_grid_files, self.cfg.discharge_files):
                self.logger.info('-' * 60)
                self.logger.info(f'Calculating catchment volumes from runoff depth grid: {runoff_file}')
                ds = grid_to_catchment(
                    runoff_file,
                    weight_table=self.cfg.grid_weights_file,
                    runoff_var=self.cfg.var_runoff_depth,
                    x_var=self.cfg.var_x,
                    y_var=self.cfg.var_y,
                    time_var=self.cfg.var_t,
                    river_id_var=self.cfg.var_river_id,
                    cumulative=self.cfg.runoff_accumulation_type == 'cumulative',
                    as_volumes=self._catchment_runoff_as_volume
                )
                yield (
                    ds['time'].values.astype('datetime64[s]'),
                    ds[self.cfg.var_catchment_runoff_variable].values.astype(np.float64, copy=False),
                    runoff_file, discharge_file
                )

    def _validate_router_configs(self) -> None:
        lateral = self.cfg.catchment_runoff_files
        grids = self.cfg.runoff_grid_files and self.cfg.grid_weights_file

        if lateral and grids:
            raise ValueError('Provide catchment_runoff_files or runoff_grid_files with grid_weights_file, not both')
        if not lateral and not grids:
            raise ValueError('Provide catchment_runoff_files or runoff_grid_files with grid_weights_file')
        n_inputs = len(lateral) + len(self.cfg.runoff_grid_files or [])
        if len(self.cfg.discharge_files) != n_inputs:
            raise ValueError('Number of resolved discharge output files must match number of input files')
        return

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

        # check that time options have the correct sizes
        if self.dt_total < self.dt_runoff:
            raise ValueError('dt_total must be >= dt_runoff')
        if self.dt_total < self.dt_discharge:
            raise ValueError('dt_total must be >= dt_discharge')
        if self.dt_discharge < self.dt_runoff:
            raise ValueError('dt_discharge must be >= dt_runoff')
        if self.dt_runoff < self.dt_routing:
            raise ValueError('dt_runoff must be >= dt_routing')
        # check that time options are evenly divisible
        if self.dt_total % self.dt_runoff != 0:
            raise ValueError('dt_total must be an integer multiple of dt_runoff')
        if self.dt_total % self.dt_discharge != 0:
            raise ValueError('dt_total must be an integer multiple of dt_discharge')
        if self.dt_discharge % self.dt_runoff != 0:
            raise ValueError('dt_discharge must be an integer multiple of dt_runoff')
        if self.dt_runoff % self.dt_routing != 0:
            raise ValueError('dt_runoff must be an integer multiple of dt_routing')

        # set derived datetime parameters for computation cycles later
        self.num_runoff_steps = int(self.dt_total / self.dt_runoff)
        self.num_runoff_steps_per_discharge = int(self.dt_discharge / self.dt_runoff)
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)
        self.num_routing_steps = int(self.dt_total / self.dt_routing)

        self._set_muskingum_coefficients(self.dt_routing)
        self._network_time_signature = signature
        return

    def _execute_routing(self) -> None:
        self._ensemble_member_states = []

        for dates, qlateral, runoff_file, discharge_file in self._catchment_runoff_generator():
            self.logger.info(f'Routing lateral inflow: {runoff_file}')
            self._set_network_and_time_dependent_vectors(dates)
            self.logger.info('Starting routing computation')
            q_t, q_array = self._router(qlateral)
            q_array[q_array < 0] = 0
            if self.cfg.runoff_processing_mode == 'sequential':
                self.logger.debug('Updating Channel State for Next Sequential Computation')
                self.channel_state = q_t
            elif self.cfg.runoff_processing_mode == 'ensemble':
                self.logger.debug('Recording Member State for Final State Aggregation')
                self._ensemble_member_states.append(q_t.copy())

            if self.dt_discharge > self.dt_runoff:
                self.logger.info('Resampling dates and discharges to specified timestep')
                q_array = (
                    q_array
                    .reshape((
                        int(self.dt_total / self.dt_discharge),
                        int(self.dt_discharge / self.dt_runoff),
                        self.A.shape[0],
                    ))
                    .mean(axis=1)
                )
                dates = dates[::self.num_runoff_steps_per_discharge]

            self.logger.info('Writing Discharge Array to File')
            np.round(q_array, decimals=2, out=q_array)
            q_array = q_array.astype(np.float32, copy=False)
            self._write_discharges(dates, q_array, discharge_file, runoff_file)

        if self.cfg.runoff_processing_mode == 'ensemble':
            self.channel_state = np.array(self._ensemble_member_states).mean(axis=0)
        self.logger.info('-' * 60)
        return

    @abstractmethod
    def _router(self, lateral: FloatArray) -> tuple[FloatArray, FloatArray]:
        ...
