import datetime
import json
import logging
import os
import random
import sys
import traceback
from pathlib import Path
from typing import Any, Self
from typing import Tuple

import netCDF4 as nc
import numpy as np
import pandas as pd
import tqdm
import yaml
from scipy.sparse import csc_matrix, diags, eye
from scipy.sparse.linalg import factorized

from .types import ConfigDict, IntArray, FactorizedSolveFn, WriteDischargesFn, DatetimeArray
from .types import FloatArray, PathInput
from ..__metadata__ import __version__
from ..tools import adjacency_matrix

__all__ = ['Router', ]


class Router:
    conf: ConfigDict
    logger: logging.Logger

    # Network dependent matrices and vectors from routing parameters file
    A: csc_matrix  # n x n - adjacency matrix
    river_ids: IntArray  # n x 1 - river ID for each segment
    downstream_river_ids: IntArray  # n x 1 - downstream river ID for each segment
    k: FloatArray  # n x 1 - K values for each segment
    x: FloatArray  # n x 1 - X values for each segment

    # Network and routing timestep dependent matrices and vectors
    c1: FloatArray  # n x 1 - C1 values for each segment => f(k, x, dt_routing)
    c2: FloatArray  # n x 1 - C2 values for each segment => f(k, x, dt_routing)
    c3: FloatArray  # n x 1 - C3 values for each segment => f(k, x, dt_routing)
    lhs: csc_matrix  # n x n - left hand side of matrix form routing equation => f(A, c1)
    lhs_factorized: FactorizedSolveFn  # callable factorized left hand side matrix for direct => f(lhs)

    # State variables
    channel_state: FloatArray
    _network_time_signature: Tuple[Any, ...] | None = None  # check if time params change between computes

    # Time options
    dt_routing: int  # compute time step, must be divisible into and <= min(dt_runoff, dt_discharge)
    dt_runoff: int  # time between lateral inflows, must be 1) constant, divisible into and <= dt_total
    dt_discharge: int  # time to average routed flows and save them, must be <=
    dt_total: int  # how long to simulate,

    # methods that are overridable via dependency injection
    _write_discharges: WriteDischargesFn

    def __init__(self, config_file: PathInput | None = None, **kwargs: Any) -> None:
        self.logger = logging.getLogger(f'river_route.{''.join([str(random.randint(0, 9)) for _ in range(8)])}')
        self._set_configs(config_file, **kwargs)
        return

    def __repr__(self):
        messages = ['Configs:', ] + [f'\t{k}: {v}' for k, v in self.conf.items()]
        return '\n'.join(messages)

    def _set_configs(self, config_file: PathInput | None, **kwargs: Any) -> None:
        if config_file is None or config_file == '':
            self.conf = {}
        elif str(config_file).endswith('.json'):
            with open(config_file, 'r') as f:
                self.conf = json.load(f)
        elif str(config_file).endswith('.yml') or str(config_file).endswith('.yaml'):
            with open(config_file, 'r') as f:
                self.conf = yaml.load(f, Loader=yaml.FullLoader)
        else:
            raise RuntimeError('Unrecognized simulation config file type. Must be .json or .yaml')

        self.conf.update(kwargs)
        self.conf['river-route-version'] = __version__
        self.conf['log'] = bool(self.conf.get('log', True))
        self.conf['progress_bar'] = self.conf.get('progress_bar', self.conf['log'])
        self.conf['log_level'] = self.conf.get('log_level', 'INFO')

        # compute and routing options (time is validated separately at compute step)
        self.conf['input_type'] = self.conf.get('input_type', 'sequential')
        self.conf['runoff_type'] = self.conf.get('runoff_type', 'incremental')

        # expected variable names in input/output files
        self.conf['var_river_id'] = self.conf.get('var_river_id', 'river_id')
        self.conf['var_discharge'] = self.conf.get('var_discharge', 'Q')
        self.conf['var_x'] = self.conf.get('var_x', 'x')
        self.conf['var_y'] = self.conf.get('var_y', 'y')
        self.conf['var_t'] = self.conf.get('var_t', 'time')
        self.conf['var_catchment_volume'] = self.conf.get('var_catchment_volume', 'volume')
        self.conf['var_runoff_depth'] = self.conf.get('var_runoff_depth', 'ro')

        # configure logging
        self.logger.disabled = not self.conf.get('log', True)
        self.logger.setLevel(self.conf.get('log_level', 'INFO'))
        log_destination = self.conf.get('log_stream', 'stdout')
        log_format = self.conf.get('log_format', '%(asctime)s - %(levelname)s - %(message)s')
        if log_destination == 'stdout':
            self.logger.addHandler(logging.StreamHandler(sys.stdout))
        elif isinstance(log_destination, str):
            self.logger.addHandler(logging.FileHandler(log_destination))
        self.logger.handlers[0].setFormatter(logging.Formatter(log_format))
        self.logger.debug('Logger initialized')
        return

    def _validate_configs(self) -> None:
        self.logger.debug('Validating configs file')
        self._validate_router_configs()
        if not self.conf.get('routing_params_file'):
            raise ValueError('routing_params_file is required')
        if not os.path.exists(self.conf['routing_params_file']):
            raise FileNotFoundError('Routing params file not found')

        # check for valid options
        if self.conf['input_type'] not in ['sequential', 'ensemble']:
            raise ValueError('Input type not recognized')
        if self.conf['runoff_type'] not in ['incremental', 'cumulative']:
            raise ValueError('Runoff type not recognized')

        # if the initial state was provided but is falsey then remove it so future checks behave properly
        if 'channel_state_file' in self.conf and not self.conf['channel_state_file']:
            del self.conf['channel_state_file']
        if self.conf.get('channel_state_file', ''):
            if not os.path.exists(self.conf['channel_state_file']):
                raise FileNotFoundError('channel_state_file not found')

        # convert all relative paths to absolute paths (after subclass has populated all file keys)
        for arg in [k for k in self.conf.keys() if 'file' in k]:
            if isinstance(self.conf[arg], list):
                self.conf[arg] = [os.path.abspath(path) for path in self.conf[arg]]
            elif isinstance(self.conf[arg], (str, Path)):
                self.conf[arg] = os.path.abspath(self.conf[arg])

        if not self.conf.get('dt_routing'):
            raise ValueError('dt_routing is required for Muskingum routing')

        self._hook_after_validate_configs()
        return

    def _validate_router_configs(self) -> None:
        """Called at the start of _validate_configs() for subclass-specific validation."""
        discharge_files = self.conf.get('discharge_files', [])
        if not discharge_files:
            raise ValueError('discharge_files is required for Muskingum routing')
        if len(discharge_files) != 1:
            raise ValueError('Muskingum routing requires exactly one entry in discharge_files')
        directory = os.path.dirname(os.path.abspath(discharge_files[0]))
        if not os.path.exists(directory):
            raise NotADirectoryError(f'Output file directory not found at: {directory}')
        if not self.conf.get('dt_total'):
            raise ValueError('dt_total is required for Muskingum routing')
        return

    ################################################
    # Computation state handling
    ################################################

    def _read_initial_state(self) -> None:
        if hasattr(self, 'channel_state'):
            return

        state_file = self.conf.get('channel_state_file', '')
        if state_file == '':
            self.logger.warning('channel_state_file not provided. Defaulting to zero initial conditions')
            self.channel_state = np.zeros(self.A.shape[0], dtype=np.float64)
            return
        self.logger.debug('Reading Initial State from Parquet')
        self.channel_state = pd.read_parquet(state_file).values.flatten().astype(np.float64, copy=False)
        return

    def _write_final_state(self) -> None:
        final_state_file = self.conf.get('final_channel_state_file', '')
        if final_state_file == '':
            return
        self.logger.debug('Writing Final State to Parquet')
        pd.DataFrame({'Q': self.channel_state}).to_parquet(self.conf['final_channel_state_file'])
        return

    ################################################
    # Prepare arrays for routing
    ################################################

    def _set_network_dependent_vectors(self) -> None:
        self.logger.debug('Calculating network dependent vectors')
        try:
            df = pd.read_parquet(
                self.conf['routing_params_file'],
                columns=[self.conf['var_river_id'], 'k', 'x', 'downstream_river_id']
            )
        except Exception as e:
            self.logger.error(f'Error reading required parameter columns from routing_params_file: {e}')
            self.logger.debug(traceback.format_exc())
            raise

        if df[self.conf['var_river_id']].duplicated().any():
            raise ValueError('routing_params_file contains duplicate river IDs.')

        self.river_ids = df[self.conf['var_river_id']].to_numpy(dtype=np.int64, copy=False)
        self.downstream_river_ids = df['downstream_river_id'].to_numpy(dtype=np.int64, copy=False)
        self.k = df['k'].to_numpy(dtype=np.float64, copy=False)
        self.x = df['x'].to_numpy(dtype=np.float64, copy=False)

        river_id_set = set(self.river_ids.tolist())
        downstream_ids = {d for d in self.downstream_river_ids.tolist() if d > 0}
        unknown_downstream_ids = sorted(downstream_ids - river_id_set)
        if unknown_downstream_ids:
            raise ValueError(
                f'routing_params_file has downstream IDs not in river_id column: {unknown_downstream_ids[:10]}')
        self.A = adjacency_matrix(self.river_ids, self.downstream_river_ids)
        self._hook_after_read_network()
        return

    def _set_muskingum_coefficients(self, dt_routing: float) -> None:
        self.logger.debug('Calculating Muskingum coefficients')
        dt_div_k = dt_routing / self.k
        denominator = dt_div_k + (2 * (1 - self.x))
        _2x = 2 * self.x
        self.c1 = (dt_div_k + _2x) / denominator
        self.c2 = (dt_div_k - _2x) / denominator
        self.c3 = ((2 * (1 - self.x)) - dt_div_k) / denominator

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
        return

    ################################################
    # Main method to execute routing simulation
    ################################################

    def route(self) -> Self:
        """
        Execute the simulation described by the provided configs and routing parameters. All configs, file paths,
        parameters, and options must be set when the object is initialized so that validation is performed before the
        simulation.

        Returns:
            Self: the class instance with updated channel_state and output files written to disk
        """
        self.logger.info('Beginning routing')
        t1 = datetime.datetime.now()

        # hooks
        self._hook_before_route()

        # validate configuration options
        self._validate_configs()
        self.logger.debug(self)

        # time parameters
        self.dt_routing = self.conf['dt_routing']
        self.dt_total = self.conf['dt_total']
        self.dt_discharge = self.conf.get('dt_discharge', self.dt_routing)
        assert self.dt_total >= self.dt_discharge >= self.dt_routing, 'Need dt_total >= dt_discharge >= dt_routing'
        assert self.dt_total % self.dt_discharge == 0, 'dt_total must be an integer multiple of dt_discharge'
        assert self.dt_discharge % self.dt_routing == 0, 'dt_discharge must be an integer multiple of dt_routing'
        num_output_steps = int(self.dt_total / self.dt_discharge)
        num_routing_per_output = int(self.dt_discharge / self.dt_routing)

        # set arrays for routing
        self._set_network_dependent_vectors()
        self._set_muskingum_coefficients(self.dt_routing)
        self._read_initial_state()

        discharge_array = self._router(num_output_steps, num_routing_per_output)

        # Generate date array for output
        start = np.datetime64(self.conf.get('start_datetime', '2000-01-01'))
        dt_sec = int(self.dt_discharge)
        dates = start + (np.arange(num_output_steps) * np.timedelta64(dt_sec, 's'))

        # write outputs and states
        self.logger.info('Writing Discharge Array to File')
        np.round(discharge_array, decimals=2, out=discharge_array)
        discharge_array = discharge_array.astype(np.float32, copy=False)
        self._write_discharges(dates, discharge_array, self.conf['discharge_files'][0])  # todo update signature
        self._write_final_state()

        # end hook
        self._hook_after_route()

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
        q_init = self.channel_state
        if not np.any(q_init):
            self.logger.warning(
                'Initial channel state is all zeros. Muskingum routing without lateral inflow requires a '
                'non-zero initial state to produce meaningful results. Provide channel_state_file.'
            )

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

        self.logger.debug('Updating Channel State')
        self.channel_state = q_t

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return discharge_array

    ################################################
    # Hooks so that subclasses can cleanly inject behavior with less overriding or duplicating
    ################################################

    def _more_config_validation(self) -> None:
        """Called at the end of _validate_configs() for subclass-specific validation."""
        pass

    def _hook_before_route(self) -> None:
        """Called at the start of route(), before validation. Default: no-op."""
        pass

    def _hook_after_route(self) -> None:
        """Called at the end of route(), after writing final state. Default: no-op."""
        pass

    def _hook_after_validate_configs(self) -> None:
        """Called after _validate_configs(). Default: no-op."""
        pass

    def _hook_after_read_network(self) -> None:
        """Called after _set_network_dependent_vectors(). Default: no-op."""
        pass

    ################################################
    # Class promotion methods for users to add behavior without calling the subclass
    ################################################

    # TODO: implement class-promotion method for rapid-style transformer
    def enable_transformer_rapid(self):
        raise NotImplementedError

    # TODO: implement class-promotion method for unit-hydrograph transformer
    def enable_transformer_unit_hydrograph(self):
        raise NotImplementedError

    ################################################
    # Dependency injection methods for users to overwrite default behaviors without subclassing
    ################################################

    def set_write_discharges(self, func: WriteDischargesFn) -> Self:
        """
        Overwrites the default write_discharges method to a custom function and returns the class instance so that you
        can chain the method with the constructor.

        Args:
            func (callable): function that takes dates, discharge_array, discharge_file, runoff_file and returns None
        """
        self._write_discharges = func
        return self

    def _write_discharges(self,
                          dates: DatetimeArray,
                          q_array: FloatArray,
                          q_file: PathInput,
                          routed_file: PathInput = '', ) -> None:
        """
        Writes routed discharge from a routing simulation to a netcdf file.
        You can overwrite this method with a custom handler using set_write_discharges.

        Args:
            dates: datetime array corresponding to the discharge rows
            q_array: routed discharge values with shape (time, river)
            q_file: path to write the discharge data to
            routed_file: path to the lateral inflow used to generate the discharge values, if applicable.

        Returns:
            None
        """
        with nc.Dataset(str(q_file), mode='w', format='NETCDF4') as ds:
            ds.createDimension('time', size=q_array.shape[0])
            ds.createDimension(self.conf['var_river_id'], size=q_array.shape[1])
            ds.runoff_file = str(routed_file)
            time_var = ds.createVariable('time', 'f8', ('time',))
            time_var.units = f'seconds since {pd.Timestamp(dates[0]).strftime("%Y-%m-%d %H:%M:%S")}'
            time_var[:] = (dates - dates[0]).astype('timedelta64[s]').astype(np.int64)
            id_var = ds.createVariable(self.conf['var_river_id'], 'i4', (self.conf['var_river_id']), )
            id_var[:] = self.river_ids
            flow_var = ds.createVariable(self.conf['var_discharge'], 'f4', ('time', self.conf['var_river_id']))
            flow_var[:] = q_array
            flow_var.long_name = 'Discharge at catchment outlet'
            flow_var.standard_name = 'discharge'
            flow_var.aggregation_method = 'mean'
            flow_var.units = 'm3 s-1'
        return
