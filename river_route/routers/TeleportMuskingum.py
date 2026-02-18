import datetime
from typing import Any, Generator

import numpy as np
import pandas as pd
import tqdm
import xarray as xr
from scipy.sparse import diags
from scipy.sparse import eye
from scipy.sparse.linalg import factorized

from ..runoff import depth_to_volume
from .typing import FloatArray, DatetimeArray, PathInput, WriteDischargesFn
from .Muskingum import Muskingum

__all__ = ['TeleportMuskingum', ]

GeneratorSignature = Generator[tuple[DatetimeArray, FloatArray, PathInput, PathInput], None, None]


class TeleportMuskingum(Muskingum):
    """
    A class for creating a "Teleport" Muskingum model.

    In this model, runoff depths/volumes are assumed to uniformly discharge to the river segment during the time step
    of the runoff data. That is, for runoff data in hourly increments, the full volume is "teleported" to the inlet at
    a rate of (depth * area) / 3600 seconds over that same hour.  It is added to the inflow from the upstream
    catchment(s). Overland flow is entirely ignored. The discharge is routed between catchments using the Muskingum
    method.

    Required configs:
    - routing_params_file: path to a parquet file containing routing parameters for each river segment. See docs.

    # option 1
    - catchment_volumes_files: list of 1 or more paths to netcdf files containing catchment volume time series for each
        river segment at each time step. See docs.
    # option 2
    - runoff_depths_files: list of 1 or more paths to netcdf files containing runoff depth time series for each river
        segment at each time step. Must also have a corresponding weight_table_file. See docs.
    - weight_table_file: path to a netCDF file containing the weights for converting runoff depths to volumes for each
        river segment. Must accompany runoff_depths_files. See docs.

    - discharge_files: list of 1 or more paths to netCDF files where the routed discharge time series will be written
        for each river segment. Must match the number of input files (catchment_volumes_files or runoff_depths_files).
        See docs.
    """
    _network_time_signature: tuple[Any, ...] | None = None  # track if sequential runoffs use the same routing matrices

    # State variables
    _ensemble_member_states: list[FloatArray]  # for ensemble routing, stores member states for computing final state

    # Time options
    dt_runoff: float
    num_routing_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_discharge: int

    # Methods
    _write_discharges: WriteDischargesFn

    def _write_final_state(self) -> None:
        final_state_file = self.conf.get('final_state_file', '')
        if final_state_file == '':
            return
        if self.conf['input_type'] == 'sequential':
            array = np.array(self.initial_state).T
        elif self.conf['input_type'] == 'ensemble':
            array = np.array(self._ensemble_member_states).mean(axis=0)
        else:
            raise ValueError('Unexpected input type when writing final states')
        self.logger.debug('Writing Final State to Parquet')
        pd.DataFrame(array, columns=['Q', 'R']).to_parquet(self.conf['final_state_file'])
        return

    def _set_network_and_time_dependent_vectors(self, dates: DatetimeArray) -> None:
        self.logger.debug('Setting and validating time parameters')
        self.dt_runoff = self.conf.get('dt_runoff', (dates[1] - dates[0]).astype('timedelta64[s]').astype(int))
        self.dt_discharge = self.conf.get('dt_discharge', self.dt_runoff)
        self.dt_total = self.conf.get('dt_total', self.dt_runoff * dates.shape[0])
        if not self.conf.get('dt_routing', 0):
            self.logger.warning('dt_routing was not provided or is Null/False, defaulting to dt_runoff')
        self.dt_routing = self.conf.get('dt_routing', self.dt_runoff)

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

        self.logger.debug('Calculating network and time dependent vectors')
        dt_div_k = self.dt_routing / self.k
        denominator = dt_div_k + (2 * (1 - self.x))
        _2x = 2 * self.x
        self.c1 = (dt_div_k + _2x) / denominator
        self.c2 = (dt_div_k - _2x) / denominator
        self.c3 = ((2 * (1 - self.x)) - dt_div_k) / denominator

        # sum of coefficients should be 1 (or arbitrarily close) for all segments
        if not np.allclose(self.c1 + self.c2 + self.c3, 1):
            self.logger.warning('Muskingum coefficients do not sum to 1')
            self.logger.debug(f'c1: {self.c1}')
            self.logger.debug(f'c2: {self.c2}')
            self.logger.debug(f'c3: {self.c3}')
            raise ValueError('Muskingum coefficients do not sum to 1, check routing parameters and time step')

        self.lhs = eye(self.A.shape[0]) - (diags(self.c2) @ self.A)
        self.lhs = self.lhs.tocsc()
        self.logger.info('Calculating factorized LHS matrix')
        self.lhs_factorized = factorized(self.lhs)
        self._network_time_signature = signature
        return

    def _volumes_generator(self) -> GeneratorSignature:
        if self.conf.get('catchment_volumes_files', False):
            for volume_file, discharge_file in zip(self.conf['catchment_volumes_files'], self.conf['discharge_files']):
                self.logger.info(f'Reading catchment volumes file: {volume_file}')
                with xr.open_dataset(volume_file) as runoff_ds:
                    self.logger.debug('Reading time array')
                    dates = runoff_ds['time'].values.astype('datetime64[s]')
                    self.logger.debug('Reading volume array')
                    volumes_array = runoff_ds[self.conf['var_catchment_volume']].values
                yield dates, volumes_array, volume_file, discharge_file
        elif self.conf.get('runoff_depths_files', False):
            for runoff_file, discharge_file in zip(self.conf['runoff_depths_files'], self.conf['discharge_files']):
                self.logger.info(f'Calculating catchment volumes from runoff depths: {runoff_file}')
                volumes_ds = depth_to_volume(
                    runoff_file,
                    weight_table=self.conf['weight_table_file'],
                    runoff_var=self.conf['var_runoff_depth'],
                    x_var=self.conf['var_x'],
                    y_var=self.conf['var_y'],
                    time_var=self.conf['var_t'],
                    river_id_var=self.conf['var_river_id'],
                    cumulative=self.conf['runoff_type'] == 'cumulative'
                )
                yield volumes_ds['time'].values, volumes_ds[
                    self.conf['var_catchment_volume']].values, runoff_file, discharge_file
        else:
            raise ValueError('No runoff data found in configs. Provide catchment volumes or runoff depths.')

    def route(self) -> 'TeleportMuskingum':
        """
        Performs time-iterative runoff routing through the river network

        Returns:
            river_route.LumpedMuskingum
        """
        self.logger.info(f'Beginning routing')
        t1 = datetime.datetime.now()

        self._validate_configs()
        self._set_network_dependent_vectors()
        self.logger.debug(self)

        for dates, volumes_array, runoff_file, discharge_file in self._volumes_generator():
            self.logger.info('-' * 80)
            self._set_network_and_time_dependent_vectors(dates)
            volumes_array = volumes_array.astype(np.float64, copy=False)
            np.divide(volumes_array, self.dt_runoff, out=volumes_array)  # convert volume -> volume/time
            discharge_array = self._router(dates, volumes_array)

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

            self.logger.info('Writing Discharge Array to File')
            np.round(discharge_array, decimals=2, out=discharge_array)
            discharge_array = discharge_array.astype(np.float32, copy=False)
            self._write_discharges(dates, discharge_array, discharge_file, runoff_file)

        # write the final state to disc
        self._write_final_state()

        t2 = datetime.datetime.now()
        self.logger.info('All runoff files routed')
        self.logger.info(f'Total job time: {(t2 - t1).total_seconds()}')
        return self

    def _router(self, dates: DatetimeArray, volumes: FloatArray) -> FloatArray:
        self.logger.debug('Getting initial state arrays')
        self._read_initial_state()
        q_init, r_init = self.initial_state  # n x 1 arrays of initial discharges and runoff in each segment
        self._ensemble_member_states = []

        # declare arrays for routing computations once to avoid repeated and duplicate allocations in the loop
        discharge_array = np.zeros((volumes.shape[0], self.A.shape[0]))
        q_t = np.empty_like(q_init, dtype=np.float64)
        r_prev = np.empty_like(r_init, dtype=np.float64)
        rhs = np.zeros(self.A.shape[0], dtype=np.float64)
        buffer = np.zeros(self.A.shape[0], dtype=np.float64)
        c2_r_t = np.zeros(self.A.shape[0], dtype=np.float64)
        interval_sum = np.zeros(self.A.shape[0], dtype=np.float64)
        q_t[:] = q_init
        r_prev[:] = r_init

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
            self.logger.debug('Updating Initial State for Next Sequential Computation')
            self.initial_state = (q_t, r_prev)
        if self.conf['input_type'] == 'ensemble':
            self.logger.debug('Recording Member State for Final State Aggregation')
            self._ensemble_member_states.append(np.array([q_t, r_prev]).T)

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return discharge_array
