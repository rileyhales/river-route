import datetime
import os
from typing import Any, Self

import numpy as np
import tqdm

from .AbstractRouter import AbstractRouter
from .types import FloatArray, PathInput

__all__ = ['Muskingum', ]


class Muskingum(AbstractRouter):
    """
    Muskingum channel routing of discharge without lateral inflow.

    Required:
        routing_params_file: parquet file with river_id, downstream_river_id, k, x columns
        initial_state_file: parquet file with Q and R columns for initial conditions
        dt_routing: routing computation timestep in seconds
        dt_discharge: output timestep in seconds (defaults to dt_routing)
        dt_total: total simulation duration in seconds
        discharge_file: output netCDF file path
        final_state_file: parquet file path to write final state

    Optional config keys:
        start_datetime: simulation start datetime string (defaults to '2000-01-01')
    """

    def __init__(self, config_file: PathInput | None = None, **kwargs: Any):
        super().__init__(config_file, **kwargs)

    def _validate_router_configs(self) -> None:
        if not self.conf.get('discharge_file'):
            raise ValueError('discharge_file is required for Muskingum routing')
        directory = os.path.dirname(os.path.abspath(self.conf['discharge_file']))
        if not os.path.exists(directory):
            raise NotADirectoryError(f'Output file directory not found at: {directory}')
        if not self.conf.get('dt_routing'):
            raise ValueError('dt_routing is required for Muskingum routing')
        if not self.conf.get('dt_total'):
            raise ValueError('dt_total is required for Muskingum routing')

    def route(self) -> Self:
        """
        Execute the simulation described by the provided configs and routing parameters. All configs, file paths,
        parameters, and options must be set when the object is initialized so that validation is performed before the
        simulation.

        Returns:
            'Muskingum': the class instance with updated initial_state and output files written to disk
        """
        self.logger.info('Beginning routing')
        t1 = datetime.datetime.now()

        self._validate_configs()
        self._set_network_dependent_vectors()
        self.logger.debug(self)

        # Time parameters from config
        self.dt_routing = self.conf['dt_routing']
        self.dt_total = self.conf['dt_total']
        self.dt_discharge = self.conf.get('dt_discharge', self.dt_routing)

        try:
            assert self.dt_total >= self.dt_discharge >= self.dt_routing, 'Must have dt_total >= dt_discharge >= dt_routing'
            assert self.dt_total % self.dt_discharge == 0, 'dt_total must be an integer multiple of dt_discharge'
            assert self.dt_discharge % self.dt_routing == 0, 'dt_discharge must be an integer multiple of dt_routing'
        except AssertionError as e:
            self.logger.error(e)
            raise AssertionError('Time options are not valid')

        num_output_steps = int(self.dt_total / self.dt_discharge)
        num_routing_per_output = int(self.dt_discharge / self.dt_routing)

        self._set_muskingum_coefficients(self.dt_routing)
        self._read_initial_state()

        discharge_array = self._router(num_output_steps, num_routing_per_output)

        # Generate date array for output
        start = np.datetime64(self.conf.get('start_datetime', '2000-01-01'))
        dt_sec = int(self.dt_discharge)
        dates = start + (np.arange(num_output_steps) * np.timedelta64(dt_sec, 's'))

        self.logger.info('Writing Discharge Array to File')
        np.round(discharge_array, decimals=2, out=discharge_array)
        discharge_array = discharge_array.astype(np.float32, copy=False)
        self._write_discharges(dates, discharge_array, self.conf['discharge_file'])
        self._write_final_state()

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return self

    def _router(self, num_output_steps: int, num_routing_per_output: int) -> FloatArray:
        """
        Route discharge without lateral inflow.

        Muskingum channel routing equation:
            (I - c2*A) @ Q(t+1) = c1*(A @ Q(t)) + c3*Q(t)
        """
        self.logger.debug('Getting initial state arrays')
        q_init, _ = self.initial_state

        n = self.A.shape[0]
        discharge_array = np.zeros((num_output_steps, n), dtype=np.float64)
        q_t = q_init.astype(np.float64, copy=True)
        rhs = np.zeros(n, dtype=np.float64)
        buffer = np.zeros(n, dtype=np.float64)
        interval_sum = np.zeros(n, dtype=np.float64)

        t1 = datetime.datetime.now()
        output_iter = range(num_output_steps)
        if self.conf['progress_bar']:
            output_iter = tqdm.tqdm(output_iter, desc='Channel Routing')
        else:
            self.logger.info('Performing routing computation iterations')

        for output_step in output_iter:
            interval_sum.fill(0.0)
            for _ in range(num_routing_per_output):
                # rhs = c1*(A @ q_t) + c3*q_t
                buffer[:] = self.A @ q_t
                np.multiply(self.c1, buffer, out=rhs)
                np.multiply(self.c3, q_t, out=buffer)
                np.add(rhs, buffer, out=rhs)
                q_t[:] = self.lhs_factorized(rhs)
                interval_sum += q_t
            discharge_array[output_step, :] = interval_sum / num_routing_per_output

        discharge_array[discharge_array < 0] = 0

        self.logger.debug('Updating Initial State')
        self.initial_state = (q_t, np.zeros_like(q_t))

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return discharge_array
