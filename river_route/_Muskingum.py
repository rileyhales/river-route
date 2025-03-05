import datetime
import json
import logging
import os
import random
import sys
from functools import partial
from typing import Tuple

import netCDF4 as nc
import networkx as nx
import numpy as np
import pandas as pd
import tqdm
import xarray as xr
import yaml
from scipy.optimize import minimize_scalar
from scipy.sparse import csc_matrix
from scipy.sparse import diags
from scipy.sparse import eye
from scipy.sparse.linalg import cgs
from scipy.sparse.linalg import factorized

from .__metadata__ import __version__
from .metrics import kge2012
from .runoff import calc_catchment_volumes
from .tools import connectivity_to_digraph
from .tools import get_adjacency_matrix

__all__ = ['Muskingum', ]

lvl_calibrate = logging.INFO + 5
logging.addLevelName(lvl_calibrate, 'CALIBRATE')


class Muskingum:
    # Given configs
    conf: dict

    # Routing matrices and vectors
    A: csc_matrix  # n x n - adjacency matrix => f(n, connectivity)
    lhs: csc_matrix  # n x n - left hand side of matrix form of routing equation => f(n, dt_routing)
    lhs_factorized: callable  # n x n - factorized left hand side matrix for direct solutions
    k: np.array  # n x 1 - K values for each segment
    x: np.array  # n x 1 - X values for each segment
    c1: np.array  # n x 1 - C1 values for each segment => f(k, x, dt_routing)
    c2: np.array  # n x 1 - C2 values for each segment => f(k, x, dt_routing)
    c3: np.array  # n x 1 - C3 values for each segment => f(k, x, dt_routing)
    initial_state: (np.array, np.array)  # Q init, R init

    # Time options
    dt_total: float
    dt_runoff: float
    dt_outflow: float
    dt_routing: float
    num_routing_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_outflow: int

    # Solver options
    _solver_method: str = 'direct'  # changed by set_iterative_solver or set_direct_solver methods
    _solver_atol: float = 1e-5

    # Calibration variables
    _calibration_iteration_number: int

    # Methods
    write_outflows: callable

    def __init__(self, config_file: str = None, **kwargs, ):
        self.LOG = logging.getLogger(f'river_route.{''.join([str(random.randint(0, 9)) for _ in range(8)])}')
        self.set_configs(config_file, **kwargs)
        return

    def set_configs(self, config_file: str, **kwargs) -> None:
        """
        Validate simulation configs given by json, yaml, or kwargs

        You can also override the values in the config file or entirely ignore the config file by specifying config
        options using keyword arguments. See README or Docs for list of all recognized configuration options. If you do
        not want to use a config file, pass None to the first positional argument and then give kwargs.

        Args:
            config_file (str): path to a json or yaml file of configs. See README for all recognized options.

        Returns:
            None
        """
        # read the config file
        if config_file is None or config_file == '':
            self.conf = {}
        elif config_file.endswith('.json'):
            with open(config_file, 'r') as f:
                self.conf = json.load(f)
        elif config_file.endswith('.yml') or config_file.endswith('.yaml'):
            with open(config_file, 'r') as f:
                self.conf = yaml.load(f, Loader=yaml.FullLoader)
        else:
            raise RuntimeError('Unrecognized simulation config file type. Must be .json or .yaml')

        # overwrite config file values with kwargs
        self.conf.update(kwargs)

        # set default values for configs when possible
        self.conf['river-route-version'] = __version__
        self.conf['log'] = bool(self.conf.get('log', True))
        self.conf['progress_bar'] = self.conf.get('progress_bar', self.conf['log'])
        self.conf['log_level'] = self.conf.get('log_level', 'INFO')

        # compute, routing, solver options (time is validated separately at compute step)
        self.conf['input_type'] = self.conf.get('input_type', 'sequential')
        self.conf['runoff_type'] = self.conf.get('runoff_type', 'incremental')
        self._solver_atol = self.conf.get('solver_atol', self._solver_atol)

        # expected variable names in input/output files
        self.conf['var_river_id'] = self.conf.get('var_river_id', 'river_id')
        self.conf['var_discharge'] = self.conf.get('var_discharge', 'Q')
        self.conf['var_x'] = self.conf.get('var_x', 'x')
        self.conf['var_y'] = self.conf.get('var_y', 'y')
        self.conf['var_t'] = self.conf.get('var_t', 'time')
        self.conf['var_catchment_volume'] = self.conf.get('var_catchment_volume', 'volume')
        self.conf['var_runoff_depth'] = self.conf.get('var_runoff_depth', 'ro')

        self.LOG.disabled = not self.conf.get('log', True)
        self.LOG.setLevel(self.conf.get('log_level', 'INFO'))
        log_destination = self.conf.get('log_stream', 'stdout')
        log_format = self.conf.get('log_format', '%(asctime)s - %(levelname)s - %(message)s')
        if log_destination == 'stdout':
            self.LOG.addHandler(logging.StreamHandler(sys.stdout))
        elif isinstance(log_destination, str):
            self.LOG.addHandler(logging.FileHandler(log_destination))
        self.LOG.handlers[0].setFormatter(logging.Formatter(log_format))
        self.LOG.debug('Logger initialized')
        return

    def _validate_configs(self) -> None:
        self.LOG.debug('Validating configs file')
        # these 2 files must always exist
        assert os.path.exists(self.conf['routing_params_file']), FileNotFoundError('Routing params file not found')
        assert os.path.exists(self.conf['connectivity_file']), FileNotFoundError('Connectivity file not found')

        # check for valid options
        assert self.conf['input_type'] in ['sequential', 'ensemble'], 'Input type not recognized'
        assert self.conf['runoff_type'] in ['incremental', 'cumulative'], 'Runoff type not recognized'

        # format conversion for the inputs
        if isinstance(self.conf.get('catchment_volumes_file', []), str):
            self.conf['catchment_volumes_file'] = [self.conf['catchment_volumes_file'], ]
        if isinstance(self.conf.get('runoff_depths_file', []), str):
            self.conf['runoff_depths_file'] = [self.conf['runoff_depths_file'], ]
        if isinstance(self.conf['outflow_file'], str):
            self.conf['outflow_file'] = [self.conf['outflow_file'], ]

        # if the weight table file is provided, it must exist
        if self.conf.get('weight_table_file', ''):
            assert os.path.exists(self.conf['weight_table_file']), FileNotFoundError('Weight table file not found')

        # if the initial state file is provided, it must exist
        if self.conf.get('initial_state_file', ''):
            assert os.path.exists(self.conf['initial_state_file']), FileNotFoundError('Initial state file not found')
        # if the key exists and the path is '' then delete the key
        if 'initial_state_file' in self.conf and not self.conf['initial_state_file']:
            del self.conf['initial_state_file']

        # either the catchment volumes or runoff depths files must be provided but not both
        assert self.conf.get('catchment_volumes_file', False) or self.conf.get('runoff_depths_file', False), \
            'Either catchment volumes or runoff depths files must be provided'

        # if catchment volumes files are provided, they must exist
        if self.conf.get('catchment_volumes_file', False):
            for file in self.conf['catchment_volumes_file']:
                assert os.path.exists(file), FileNotFoundError(f'Runoff volumes file not found at: {file}')

        # if runoff depths files are provided, they must exist
        if self.conf.get('runoff_depths_file', False):
            for file in self.conf['runoff_depths_file']:
                assert os.path.exists(file), FileNotFoundError(f'Runoff depths file not found at: {file}')
            assert self.conf.get('weight_table_file', ''), 'weight_table_file should be provided'

        # the directory for each outflow file must exist
        for outflow_file in self.conf['outflow_file']:
            assert os.path.exists(os.path.dirname(outflow_file)), NotADirectoryError('Output directory not found')

        # convert all relative paths to absolute paths
        for arg in [k for k in self.conf.keys() if 'file' in k]:
            if isinstance(self.conf[arg], list):
                self.conf[arg] = [os.path.abspath(path) for path in self.conf[arg]]
            elif isinstance(self.conf[arg], str):
                self.conf[arg] = os.path.abspath(self.conf[arg])
        return

    def _log_configs(self) -> None:
        self.LOG.debug(f'river-route version: {__version__}')
        self.LOG.debug(f'Number of Rivers: {self._read_river_ids().shape[0]}')
        self.LOG.debug('Configs:')
        for k, v in self.conf.items():
            self.LOG.debug(f'\t{k}: {v}')
        return

    def _read_river_ids(self) -> np.array:
        return pd.read_parquet(self.conf['routing_params_file'], columns=[self.conf['var_river_id'], ]).values.flatten()

    def _set_linear_routing_params(self) -> None:
        if hasattr(self, 'k') and hasattr(self, 'x'):
            return
        self.k = pd.read_parquet(self.conf['routing_params_file'], columns=['k', ]).values.flatten()
        self.x = pd.read_parquet(self.conf['routing_params_file'], columns=['x', ]).values.flatten()
        return

    def _read_initial_state(self) -> None:
        if hasattr(self, 'initial_state'):
            return

        state_file = self.conf.get('initial_state_file', '')
        if state_file == '':
            self.LOG.debug('Setting initial state to zeros')
            self.initial_state = (np.zeros(self.A.shape[0]), np.zeros(self.A.shape[0]))
            return
        assert os.path.exists(state_file), FileNotFoundError(f'Initial state file not found at: {state_file}')
        self.LOG.debug('Reading Initial State from Parquet')
        initial_state = pd.read_parquet(state_file).values
        initial_state = (initial_state[:, 0].flatten(), initial_state[:, 1].flatten())
        self.initial_state = initial_state
        return

    def _write_final_state(self) -> None:
        final_state_file = self.conf.get('final_state_file', '')
        if final_state_file == '':
            return
        self.LOG.debug('Writing Final State to Parquet')
        pd.DataFrame(np.array(self.initial_state).T, columns=['Q', 'R']).to_parquet(self.conf['final_state_file'])
        return

    def _get_digraph(self) -> nx.DiGraph:
        return connectivity_to_digraph(self.conf['connectivity_file'])

    def _set_adjacency_matrix(self) -> None:
        if hasattr(self, 'A'):
            return
        self.LOG.debug('Calculating Network Adjacency Matrix (A)')
        self.A = get_adjacency_matrix(self.conf['routing_params_file'], self.conf['connectivity_file'])
        return

    def _set_lhs_matrix(self) -> None:
        if not hasattr(self, 'lhs'):
            self.lhs = eye(self.A.shape[0]) - (diags(self.c2) @ self.A)
            self.lhs = self.lhs.tocsc()
            if hasattr(self, 'lhs_factorized'):
                del self.lhs_factorized
        if not hasattr(self, 'lhs_factorized') and self._solver_method == 'direct':
            self.LOG.info('Calculating factorized LHS matrix')
            self.lhs_factorized = factorized(self.lhs)
        return

    def _set_time_params(self, dates: np.array) -> None:
        """
        Set time parameters for the simulation
        """
        self.LOG.debug('Setting and validating time parameters')
        self.dt_runoff = self.conf.get('dt_runoff', (dates[1] - dates[0]).astype('timedelta64[s]').astype(int))
        self.dt_outflow = self.conf.get('dt_outflow', self.dt_runoff)
        self.dt_total = self.conf.get('dt_total', self.dt_runoff * dates.shape[0])
        if not self.conf.get('dt_routing', 0):
            self.LOG.warning('dt_routing was not provided or is Null/False, defaulting to dt_runoff')
        self.dt_routing = self.conf.get('dt_routing', self.dt_runoff)

        try:
            # check that time options have the correct sizes
            assert self.dt_total >= self.dt_runoff, 'dt_total !>= dt_runoff'
            assert self.dt_total >= self.dt_outflow, 'dt_total !>= dt_outflow'
            assert self.dt_outflow >= self.dt_runoff, 'dt_outflow !>= dt_runoff'
            assert self.dt_runoff >= self.dt_routing, 'dt_runoff !>= dt_routing'
            # check that time options are evenly divisible
            assert self.dt_total % self.dt_runoff == 0, 'dt_total must be an integer multiple of dt_runoff'
            assert self.dt_total % self.dt_outflow == 0, 'dt_total must be an integer multiple of dt_outflow'
            assert self.dt_outflow % self.dt_runoff == 0, 'dt_outflow must be an integer multiple of dt_runoff'
            assert self.dt_runoff % self.dt_routing == 0, 'dt_runoff must be an integer multiple of dt_routing'
        except AssertionError as e:
            self.LOG.error(e)
            raise AssertionError('Time options are not valid')

        # set derived datetime parameters for computation cycles later
        self.num_runoff_steps_per_outflow = int(self.dt_outflow / self.dt_runoff)
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)
        self.num_routing_steps = int(self.dt_total / self.dt_routing)
        return

    def _set_muskingum_coefficients(self, k: np.ndarray = None, x: np.ndarray = None) -> None:
        """
        Calculate the 3 Muskingum routing coefficients for each segment using given k and x
        """
        self.LOG.debug('Calculating Muskingum coefficients')

        self._set_linear_routing_params()
        k = k if k is not None else self.k
        x = x if x is not None else self.x

        dt_div_k = self.dt_routing / k
        denom = dt_div_k + (2 * (1 - x))
        _2x = 2 * x
        self.c1 = (dt_div_k + _2x) / denom
        self.c2 = (dt_div_k - _2x) / denom
        self.c3 = ((2 * (1 - x)) - dt_div_k) / denom

        # sum of coefficients should be 1 (or arbitrarily close) for all segments
        if not np.allclose(self.c1 + self.c2 + self.c3, 1):
            self.LOG.warning('Muskingum coefficients do not sum to 1')
            self.LOG.debug(f'c1: {self.c1}')
            self.LOG.debug(f'c2: {self.c2}')
            self.LOG.debug(f'c3: {self.c3}')
            raise AssertionError('Muskingum coefficients do not sum to 1')
        return

    def _volumes_output_generator(self) -> Tuple[pd.DataFrame, str]:
        if self.conf.get('catchment_volumes_file', False):
            for volume_file, outflow_file in zip(self.conf['catchment_volumes_file'], self.conf['outflow_file']):
                self.LOG.info('-' * 80)
                self.LOG.info(f'Reading catchment volumes file: {volume_file}')
                with xr.open_dataset(volume_file) as runoff_ds:
                    self.LOG.debug('Reading time array')
                    dates = runoff_ds['time'].values.astype('datetime64[s]')
                    self.LOG.debug('Reading volume array')
                    volumes_array = runoff_ds[self.conf['var_catchment_volume']].values
                yield dates, volumes_array, volume_file, outflow_file
        elif self.conf.get('runoff_depths_file', False):
            for runoff_file, outflow_file in zip(self.conf['runoff_depths_file'], self.conf['outflow_file']):
                self.LOG.info('-' * 80)
                self.LOG.info(f'Calculating catchment volumes from runoff depths: {runoff_file}')
                volumes_df = calc_catchment_volumes(
                    runoff_file,
                    weight_table=self.conf['weight_table_file'],
                    river_id_var=self.conf['var_river_id'],
                    runoff_var=self.conf['var_runoff_depth'],
                    x_var=self.conf['var_x'],
                    y_var=self.conf['var_y'],
                    time_var=self.conf['var_t'],
                    cumulative=self.conf['runoff_type'] == 'cumulative',
                )
                yield volumes_df.index.values, volumes_df.values, runoff_file, outflow_file
        else:
            raise ValueError('No runoff data found in configs. Provide catchment volumes or runoff depths.')

    def route(self, **kwargs) -> 'Muskingum':
        """
        Performs time-iterative runoff routing through the river network

        Keyword Args:

        Returns:
            river_route.Muskingum
        """
        self.LOG.info(f'Beginning routing')
        t1 = datetime.datetime.now()

        if len(kwargs) > 0:
            self.LOG.info('Updating configs with kwargs')
            self.conf.update(kwargs)

        self._validate_configs()
        self._log_configs()
        self._set_adjacency_matrix()

        for dates, volumes_array, runoff_file, outflow_file in self._volumes_output_generator():
            self._set_time_params(dates)
            volumes_array = volumes_array / self.dt_runoff  # convert volume -> volume/time
            self._set_muskingum_coefficients()
            outflow_array = self._router(dates, volumes_array)

            if self.dt_outflow > self.dt_runoff:
                self.LOG.info('Resampling dates and outflows to specified timestep')
                outflow_array = (
                    outflow_array
                    .reshape((
                        int(self.dt_total / self.dt_outflow),
                        int(self.dt_outflow / self.dt_runoff),
                        self.A.shape[0],
                    ))
                    .mean(axis=1)
                )

            self.LOG.info('Writing Outflow Array to File')
            outflow_array = np.round(outflow_array, decimals=2)
            self._write_outflows(dates, outflow_array, outflow_file, runoff_file)

        # write the final state to disc
        self._write_final_state()

        t2 = datetime.datetime.now()
        self.LOG.info('All runoff files routed')
        self.LOG.info(f'Total job time: {(t2 - t1).total_seconds()}')
        return self

    def calibrate(self, observed: pd.DataFrame, overwrite_params_file: bool = False) -> 'Muskingum':
        """
        Calibrate K and X to given measured_values using optimization algorithms

        Args:
            observed: A pandas DataFrame with a datetime index, river id column names, and discharge values
            overwrite_params_file: If True, the calibration will overwrite the routing_params_file with the new values

        Returns:
            river_route.Muskingum
        """
        self.LOG.info(f'Beginning optimization')
        t1 = datetime.datetime.now()

        self._set_linear_routing_params()
        self._validate_configs()
        self._log_configs()
        self._set_adjacency_matrix()

        # find the indices of the rivers which have observed flows
        g = self._get_digraph()
        subgraph_rivers = set()
        for river_id in observed.columns:
            subgraph_rivers.update([river_id, ])
            subgraph_rivers.update(nx.ancestors(g, river_id))
        # sort the river IDs to match the order they appear in the routing parameters file
        river_ids = pd.Series(self._read_river_ids())
        subgraph_rivers = pd.Series(list(subgraph_rivers))
        sorted_indices = subgraph_rivers.map(pd.Series(river_ids.index, index=river_ids.values)).sort_values().index
        subgraph_rivers = subgraph_rivers[sorted_indices].reset_index(drop=True).values
        # make boolean selectors to subset the inputs and select rivers to compare with observed values
        riv_idxes = river_ids.isin(subgraph_rivers).values
        river_ids = river_ids[riv_idxes].reset_index(drop=True)
        obs_idxes = river_ids[river_ids.isin(observed.columns)].index.values

        self.k = self.k[riv_idxes]
        self.x = self.x[riv_idxes]
        self.A = self.A[riv_idxes, :][:, riv_idxes]

        self._calibration_iteration_number = 0
        self.LOG.setLevel('CALIBRATE')

        obj = partial(self._calibration_objective,
                      observed=observed, riv_idxes=riv_idxes, obs_idxes=obs_idxes, var='k')
        result_k = minimize_scalar(obj, bounds=(0.75, 1.2), method='bounded', options={'xatol': 1e-2})
        self.k = self.k * result_k.x
        self.LOG.log(lvl_calibrate, f'Optimal k scalar: {result_k.x}')
        self.LOG.log(lvl_calibrate, result_k)
        self.LOG.log(lvl_calibrate, '-' * 80)

        self._calibration_iteration_number = 0
        obj = partial(self._calibration_objective,
                      observed=observed, riv_idxes=riv_idxes, obs_idxes=obs_idxes, var='x')
        result_x = minimize_scalar(obj, bounds=(0.9, 1.1), method='bounded', options={'xatol': 1e-2})
        self.x = self.x * result_x.x
        self.LOG.log(lvl_calibrate, f'Optimal x scalar: {result_x.x}')
        self.LOG.log(lvl_calibrate, result_x)

        if overwrite_params_file:
            df = pd.read_parquet(self.conf['routing_params_file'])
            df['k'] = df['k'] * result_k.x
            df['x'] = df['x'] * result_x.x
            df.to_parquet(self.conf['routing_params_file'])

        self.LOG.log(lvl_calibrate, '-' * 80)
        self.LOG.log(lvl_calibrate, f'Final k scalar: {result_k.x}')
        self.LOG.log(lvl_calibrate, f'Final x scalar: {result_x.x}')
        self.LOG.log(lvl_calibrate, f'Total iterations: {result_k.nit + result_x.nit}')
        self.LOG.log(lvl_calibrate, '-' * 80)

        self.LOG.setLevel(self.conf.get('log_level', 'INFO'))
        t2 = datetime.datetime.now()
        self.LOG.info(f'Total job time: {(t2 - t1).total_seconds()}')

        # delete calibration arrays to force recalculation
        self.LOG.debug('Deleting routing related arrays affected by calibration process')
        del self.k
        del self.x
        del self.A
        del self.lhs
        del self.c1
        del self.c2
        del self.c3
        self._set_linear_routing_params()
        return self

    def _router(self, dates: np.array, volumes: np.array) -> np.array:
        self.LOG.debug('Getting initial state arrays')
        self._read_initial_state()
        q_init, r_init = self.initial_state

        # set initial values
        self._set_lhs_matrix()
        outflow_array = np.zeros((volumes.shape[0], self.A.shape[0]))
        q_t = np.zeros_like(q_init)
        q_t[:] = q_init
        r_prev = np.zeros_like(r_init)
        r_prev[:] = r_init

        t1 = datetime.datetime.now()
        if self.conf['progress_bar']:
            dates = tqdm.tqdm(dates, desc='Runoff Volumes Routed')
        else:
            self.LOG.info('Performing routing computation iterations')
        for runoff_time_step, runoff_end_date in enumerate(dates):
            r_t = volumes[runoff_time_step, :]
            interval_flows = np.zeros((self.num_routing_steps_per_runoff, self.A.shape[0]))
            for routing_sub_iteration_num in range(self.num_routing_steps_per_runoff):
                rhs = (self.c1 * ((self.A @ q_t) + r_prev)) + (self.c2 * r_t) + (self.c3 * q_t)
                q_t[:] = self._solve(rhs, q_t)
                interval_flows[routing_sub_iteration_num, :] = q_t
            outflow_array[runoff_time_step, :] = np.mean(interval_flows, axis=0)
            r_prev[:] = r_t

        # Enforce positive flows in case of negative solver results
        outflow_array[outflow_array < 0] = 0

        # if simulation type is ensemble, then do not overwrite the initial state
        if self.conf['input_type'] == 'sequential':
            self.LOG.debug('Updating Initial State for Next Sequential Computation')
            self.initial_state = (q_t, r_prev)

        t2 = datetime.datetime.now()
        self.LOG.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return outflow_array

    def _solve(self, rhs: np.array, q_t: np.array) -> np.array:
        return self._solve_direct(rhs, q_t)

    def _solve_cgs(self, rhs: np.array, q_t: np.array) -> np.array:
        return cgs(self.lhs, rhs, x0=q_t, atol=self._solver_atol)[0]

    def _solve_direct(self, rhs: np.array, q_t: np.array) -> np.array:
        return self.lhs_factorized(rhs)

    def set_direct_solver(self, ) -> 'Muskingum':
        """
        Set using a direct solver of the factorized LHS of the matrix muskingum equation. Returns self.

        Returns:
            river_route.Muskingum
        """
        self._solver_method = 'direct'
        self._solve = self._solve_direct
        return self

    def set_iterative_solver(self) -> 'Muskingum':
        """
        Set using an iterative solver of the LHS of the matrix muskingum equation. Returns self.

        Returns:
            river_route.Muskingum
        """
        self._solver_method = 'iterative'
        self._solve = self._solve_cgs
        return self

    def _calibration_objective(self,
                               iteration: np.array, *, observed: pd.DataFrame,
                               riv_idxes: np.array, obs_idxes: np.ndarray, var: str = None) -> np.float64:
        """
        Objective function for calibration

        Args:
            iteration: a scalar value to multiply the k values by
            riv_idxes: array of the indices of rivers with observations in the river id list from params
            obs_idxes: array of the indices of observed rivers in the river id list from params

        Returns:
            np.float64
        """
        self._calibration_iteration_number += 1
        self.LOG.log(lvl_calibrate, f'Iteration {self._calibration_iteration_number} - Testing scalar: {iteration}')
        iter_df = pd.DataFrame(columns=riv_idxes)

        # set iteration k and x depending on the routing type in the configs and the size of the scalar
        assert var in ['k', 'x'], 'Calibration variable must be k or x for linear calibration'
        k = self.k if var != 'k' else self.k * iteration
        x = self.x if var != 'x' else self.x * iteration

        for dates, volumes_array, runoff_file, outflow_file in self._volumes_output_generator():
            self._set_time_params(dates)
            self._set_muskingum_coefficients(k=k, x=x)
            volumes_array = volumes_array / self.dt_runoff  # convert volume -> volume/time
            outflow_array = self._router(dates, volumes_array)[:, obs_idxes]
            if iter_df.empty:
                iter_df = pd.DataFrame(outflow_array, index=dates, columns=observed.columns)
            else:
                iter_df = pd.concat([iter_df, pd.DataFrame(outflow_array, index=dates, columns=observed.columns)])
        y_true, y_pred = observed.merge(iter_df, left_index=True, right_index=True,
                                        suffixes=('_true', '_pred')).dropna().values.T
        objective_function_value = np.absolute(kge2012(y_true, y_pred) - 1)
        self.LOG.log(lvl_calibrate,
                     f'Iteration {self._calibration_iteration_number} - ABS(KGE - 1): {objective_function_value}')
        del self.initial_state
        # return mse
        return objective_function_value

    def _write_outflows(self, dates: np.array, outflow_array: np.array, outflow_file: str, runoff_file: str) -> None:
        dates = dates[::self.num_runoff_steps_per_outflow].astype('datetime64[s]')
        df = pd.DataFrame(outflow_array, index=dates, columns=self._read_river_ids())
        self.write_outflows(df=df, outflow_file=outflow_file, runoff_file=runoff_file)
        return

    def write_outflows(self, df: pd.DataFrame, outflow_file: str, runoff_file: str) -> None:
        """
        Writes the outflows from a routing simulation to a netcdf file. You should overwrite this method with a custom
        handler that writes it in a format that fits your needs.

        Args:
            df: a Pandas DataFrame with a datetime Index, river_id column names, and discharge values
            outflow_file: the file path to write the outflows to
            runoff_file: the file path to the runoff file used to generate the outflows

        Returns:
            None
        """
        with nc.Dataset(outflow_file, mode='w', format='NETCDF4') as ds:
            ds.createDimension('time', size=df.shape[0])
            ds.createDimension(self.conf['var_river_id'], size=df.shape[1])

            time_var = ds.createVariable('time', 'f8', ('time',))
            time_var.units = f'seconds since {df.index[0].strftime("%Y-%m-%d %H:%M:%S")}'
            time_var[:] = df.index.values - df.index.values[0]

            id_var = ds.createVariable(self.conf['var_river_id'], 'i4', (self.conf['var_river_id']), )
            id_var[:] = df.columns.values

            flow_var = ds.createVariable(self.conf['var_discharge'], 'f4', ('time', self.conf['var_river_id']))
            flow_var[:] = df.values
            flow_var.long_name = 'Discharge at catchment outlet'
            flow_var.standard_name = 'discharge'
            flow_var.aggregation_method = 'mean'
            flow_var.units = 'm3 s-1'
        return

    def set_write_outflows(self, func: callable) -> 'Muskingum':
        """
        Overwrites the default write_outflows method to a custom function and returns the class instance so that you
        can chain the method with the constructor.

        Args:
            func (callable): a function that takes 3 keyword arguments: df, outflow_file, runoff_file and returns None

        Returns:
            river_route.Muskingum
        """
        self.write_outflows = func
        return self

    def hydrograph(self, river_id: int) -> pd.DataFrame:
        """
        Get the hydrograph for a given river id as a pandas dataframe

        Args:
            river_id: the ID of a river reach in the output files

        Returns:
            pandas.DataFrame
        """
        with xr.open_mfdataset(self.conf['outflow_file']) as ds:
            df = (
                ds
                [self.conf['var_discharge']]
                .sel(**{self.conf['var_river_id']: river_id})
                .to_dataframe()
                [[self.conf['var_discharge'], ]]
            )
            df.columns = [river_id, ]
        return df

    def mass_balance(self, river_id: int, ancestors: list = None) -> pd.DataFrame:
        """
        Get the mass balance for a given river id as a pandas dataframe

        Args:
            river_id: the ID of a river reach in the output files
            ancestors: a list of the given river_id and all rivers upstream of that river

        Returns:
            pandas.DataFrame
        """
        if type(river_id) is not int:
            raise TypeError(f'river_id should be an integer ID of a river to mass balance')
        if ancestors is None:
            g = connectivity_to_digraph(self.conf['connectivity_file'])
            ancestors = set(list(nx.ancestors(g, river_id)) + [river_id, ])
        with xr.open_mfdataset(self.conf['catchment_volumes_file']) as ds:
            vdf = (
                ds
                .sel(**{self.conf['var_river_id']: ancestors})
                [self.conf['var_catchment_volume']]
                .to_dataframe()
                [[self.conf['var_catchment_volume'], ]]
                .reset_index()
                .pivot(index='time', columns=self.conf['var_river_id'], values=self.conf['var_catchment_volume'])
                .sum(axis=1)
                .cumsum()
                .rename('runoff_volume')
            )
        with xr.open_mfdataset(self.conf['outflow_file']) as ds:
            qdf = (
                ds
                .sel(**{self.conf['var_river_id']: river_id})
                .to_dataframe()
                [[self.conf['var_discharge'], ]]
                .cumsum()
            )
            # convert to discharged volume - multiply by the time delta in seconds
            qdf = qdf * (qdf.index[1] - qdf.index[0]).total_seconds()
            qdf.columns = ['discharge_volume', ]

        df = qdf.join(vdf)
        df['runoff-discharge'] = df['runoff_volume'] - df['discharge_volume']
        if not df['runoff-discharge'].gt(0).all():
            self.LOG.warning(f'More discharge than runoff volume for river {river_id}')

        return df
