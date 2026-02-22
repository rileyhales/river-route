import json
import logging
import os
import random
import sys
import traceback
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Self

import netCDF4 as nc
import numpy as np
import pandas as pd
import yaml
from scipy.sparse import csc_matrix, diags, eye
from scipy.sparse.linalg import factorized

from .types import ConfigDict, PathInput, FloatArray, IntArray, FactorizedSolveFn, WriteDischargesFn, DatetimeArray
from ..__metadata__ import __version__
from ..tools import adjacency_matrix

__all__ = ['AbstractRouter', ]


class AbstractRouter(ABC):
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
    lhs: csc_matrix  # n x n - left hand side of matrix form of routing equation
    lhs_factorized: FactorizedSolveFn  # callable factorized left hand side matrix for direct

    # State variables
    initial_state: tuple[FloatArray, FloatArray]  # Q init, R init

    # Time options
    dt_total: float
    dt_discharge: float
    dt_routing: float

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

    @abstractmethod
    def _validate_router_configs(self) -> None:
        """Router-specific config validation. Called at the start of _validate_configs before common checks."""

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
        if 'initial_state_file' in self.conf and not self.conf['initial_state_file']:
            del self.conf['initial_state_file']
        if self.conf.get('initial_state_file', ''):
            if not os.path.exists(self.conf['initial_state_file']):
                raise FileNotFoundError('Initial state file not found')

        # convert all relative paths to absolute paths (after subclass has populated all file keys)
        for arg in [k for k in self.conf.keys() if 'file' in k]:
            if isinstance(self.conf[arg], list):
                self.conf[arg] = [os.path.abspath(path) for path in self.conf[arg]]
            elif isinstance(self.conf[arg], (str, Path)):
                self.conf[arg] = os.path.abspath(self.conf[arg])
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

    def _read_initial_state(self) -> None:
        if hasattr(self, 'initial_state'):
            return

        state_file = self.conf.get('initial_state_file', '')
        if state_file == '':
            # todo raise error on muskingum but warning on others? but avoid duplicating code?
            # raise RuntimeError('initial_state_file not provided. Cannot route without water in the channels')
            self.logger.warning('initial_state_file not provided. Defaulting to zero initial conditions')
            initial_state = (np.zeros(self.A.shape[0], dtype=np.float64), np.zeros(self.A.shape[0], dtype=np.float64))
            self.initial_state = initial_state
            return
        self.logger.debug('Reading Initial State from Parquet')
        initial_state = pd.read_parquet(state_file).values
        initial_state = (initial_state[:, 0].flatten(), initial_state[:, 1].flatten())
        self.initial_state = initial_state
        return

    def _write_final_state(self) -> None:
        final_state_file = self.conf.get('final_state_file', '')
        if final_state_file == '':
            return
        self.logger.debug('Writing Final State to Parquet')
        array = np.array(self.initial_state).T
        pd.DataFrame(array, columns=['Q', 'R']).to_parquet(self.conf['final_state_file'])
        return

    @abstractmethod
    def route(self) -> Self:
        """Execute routing end-to-end and return self."""

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
        return

    def _write_discharges(self,
                          dates: DatetimeArray,
                          discharge_array: FloatArray,
                          discharge_file: PathInput,
                          routed_file: PathInput = '', ) -> None:
        """
        Writes routed discharge from a routing simulation to a netcdf file.
        You can overwrite this method with a custom handler using set_write_discharges.

        Args:
            dates: datetime array corresponding to the discharge rows
            discharge_array: routed discharge values with shape (time, river)
            discharge_file: path to write the discharge data to
            routed_file: path to the lateral inflow used to generate the discharge values, if applicable.

        Returns:
            None
        """
        with nc.Dataset(str(discharge_file), mode='w', format='NETCDF4') as ds:
            ds.createDimension('time', size=discharge_array.shape[0])
            ds.createDimension(self.conf['var_river_id'], size=discharge_array.shape[1])
            ds.runoff_file = str(routed_file)
            time_var = ds.createVariable('time', 'f8', ('time',))
            time_var.units = f'seconds since {pd.Timestamp(dates[0]).strftime("%Y-%m-%d %H:%M:%S")}'
            time_var[:] = (dates - dates[0]).astype('timedelta64[s]').astype(np.int64)
            id_var = ds.createVariable(self.conf['var_river_id'], 'i4', (self.conf['var_river_id']), )
            id_var[:] = self.river_ids
            flow_var = ds.createVariable(self.conf['var_discharge'], 'f4', ('time', self.conf['var_river_id']))
            flow_var[:] = discharge_array
            flow_var.long_name = 'Discharge at catchment outlet'
            flow_var.standard_name = 'discharge'
            flow_var.aggregation_method = 'mean'
            flow_var.units = 'm3 s-1'
        return

    def set_write_discharges(self, func: WriteDischargesFn) -> Self:
        """
        Overwrites the default write_discharges method to a custom function and returns the class instance so that you
        can chain the method with the constructor.

        Args:
            func (callable): function that takes dates, discharge_array, discharge_file, runoff_file and returns None
        """
        self._write_discharges = func
        return self
