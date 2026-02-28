import datetime
from abc import ABC, abstractmethod
from typing import Self, List

import numpy as np
import xarray as xr

from .Router import Router
from ..types import DatetimeArray, FloatArray, RunoffGeneratorSignature

__all__ = ['TransformRouter', ]

from ..runoff import grid_to_catchment


class TransformRouter(Router, ABC):
    """
    Intermediate abstract router class adding methods routing methods that require pre-processing of the lateral inflow
    """
    # State variables
    _ensemble_member_states: List[FloatArray]  # for ensemble routing

    # Time options
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
                with xr.open_dataset(lateral_file) as ds:
                    dates = ds['time'].values.astype('datetime64[s]')
                    array = ds[self.cfg.var_catchment_runoff_variable].values.astype(np.float64, copy=False)
                    yield dates, array, lateral_file, discharge_file
        elif self.cfg.runoff_grid_files and self.cfg.grid_weights_file:
            for runoff_file, discharge_file in zip(self.cfg.runoff_grid_files, self.cfg.discharge_files):
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

    def _prepare_qlateral(self, array: FloatArray) -> FloatArray:
        """Pre-process the lateral array, precomputed or generated, before passing to _router"""
        return array

    def _validate_router_configs(self) -> None:
        lateral = self.cfg.catchment_runoff_files
        grids = self.cfg.runoff_grid_files and self.cfg.grid_weights_file

        if lateral and grids:
            raise ValueError('Provide catchment_runoff_files or runoff_grid_files with grid_weights_file, not both')
        if not lateral and not grids:
            raise ValueError('Provide catchment_runoff_files or runoff_grid_files with grid_weights_file')
        if len(self.cfg.discharge_files) != len(lateral) + len(grids):
            raise ValueError('Number of discharge_files must match number of input files')
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

        try:
            # check that time options have the correct sizes
            assert self.dt_total >= self.dt_runoff, 'dt_total must be >= dt_runoff'
            assert self.dt_total >= self.dt_discharge, 'dt_total must be >= dt_discharge'
            assert self.dt_discharge >= self.dt_runoff, 'dt_discharge must be >= dt_runoff'
            assert self.dt_runoff >= self.dt_routing, 'dt_runoff must be >= dt_routing'
            # check that time options are evenly divisible
            assert self.dt_total % self.dt_runoff == 0, 'dt_total must be an integer multiple of dt_runoff'
            assert self.dt_total % self.dt_discharge == 0, 'dt_total must be an integer multiple of dt_discharge'
            assert self.dt_discharge % self.dt_runoff == 0, 'dt_discharge must be an integer multiple of dt_runoff'
            assert self.dt_runoff % self.dt_routing == 0, 'dt_runoff must be an integer multiple of dt_routing'
        except AssertionError as e:
            self.logger.error(e)
            raise AssertionError('Time options are not valid')

        # set derived datetime parameters for computation cycles later
        self.num_runoff_steps_per_discharge = int(self.dt_discharge / self.dt_runoff)
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)
        self.num_routing_steps = int(self.dt_total / self.dt_routing)

        self._set_muskingum_coefficients(self.dt_routing)
        self._network_time_signature = signature
        return

    def route(self) -> Self:
        self.logger.info('Beginning routing')
        t1 = datetime.datetime.now()

        # validate configuration options
        self._validate_configs()
        self.logger.debug(self)

        # set arrays for routing
        self._set_network_dependent_vectors()
        self._read_initial_state()

        # hooks
        self._hook_before_route()
        self._ensemble_member_states = []

        for dates, qlateral, runoff_file, discharge_file in self._catchment_runoff_generator():
            self.logger.info('-' * 80)
            self._set_network_and_time_dependent_vectors(dates)
            qlateral = self._prepare_qlateral(qlateral)
            discharge_array = self._router(dates, qlateral)

            if self.dt_discharge > self.dt_runoff:
                self.logger.info('Resampling dates and discharges to specified timestep')
                discharge_array = (
                    discharge_array
                    .reshape((
                        int(self.dt_total / self.dt_discharge),
                        int(self.dt_discharge / self.dt_runoff),
                        self.A.shape[0],
                    ))
                    .mean(axis=1)
                )
                dates = dates[::self.num_runoff_steps_per_discharge]

            self.logger.info('Writing Discharge Array to File')
            np.round(discharge_array, decimals=2, out=discharge_array)
            discharge_array = discharge_array.astype(np.float32, copy=False)
            self._write_discharges(dates, discharge_array, discharge_file, runoff_file)

        if self.cfg.runoff_processing_mode == 'ensemble':
            self.channel_state = np.array(self._ensemble_member_states).mean(axis=0)
        self._write_final_state()
        self._hook_after_route()

        t2 = datetime.datetime.now()
        self.logger.info('All runoff files routed')
        self.logger.info(f'Total job time: {(t2 - t1).total_seconds()}')
        return self

    @abstractmethod
    def _router(self, dates: DatetimeArray, lateral: FloatArray) -> FloatArray:
        """Execute the core routing math for one runoff file and return the discharge array."""
