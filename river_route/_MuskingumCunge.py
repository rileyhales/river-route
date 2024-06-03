import datetime
import json
import logging
import os
import sys
from functools import partial

import netCDF4 as nc
import networkx as nx
import numpy as np
import pandas as pd
import scipy
import tqdm
import xarray as xr
import yaml
from scipy.optimize import minimize_scalar

from .__metadata__ import __version__ as VERSION
from .tools import connectivity_to_adjacency_matrix
from .tools import connectivity_to_digraph

__all__ = ['MuskingumCunge', ]

LOG = logging.getLogger('river_route')
LOG.disabled = True


class MuskingumCunge:
    # Given configs
    conf: dict

    # Routing matrices
    A: scipy.sparse.csc_matrix
    lhs: scipy.sparse.csc_matrix
    k: np.array
    x: np.array
    c1: np.array
    c2: np.array
    c3: np.array
    initial_state: (np.array, np.array)  # Q init, R init

    # Time options
    dt_total: float
    dt_runoff: float
    dt_outflow: float
    dt_routing: float
    num_routing_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_outflow: int

    # Calibration variables
    _calibration_iteration_number: int

    def __init__(self, config_file: str = None, **kwargs, ):
        self.set_configs(config_file, **kwargs)
        return

    def set_configs(self, config_file, **kwargs) -> None:
        """
        Validate simulation configs given by json, yaml, or kwargs

        Args:
            config_file: path to a json or yaml file containing configs. See README for list of all recognized options

        Keyword Args:
            See README for list of all recognized configuration options

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
        self.conf['river-route-version'] = VERSION
        self.conf['job_name'] = self.conf.get('job_name', 'untitled_job')
        self.conf['log'] = bool(self.conf.get('log', False))
        self.conf['progress_bar'] = self.conf.get('progress_bar', self.conf['log'])
        self.conf['log_level'] = self.conf.get('log_level', 'INFO')
        self.conf['var_runoff_volume'] = self.conf.get('var_runoff_volume', 'ro_vol')
        self.conf['var_river_id'] = self.conf.get('var_river_id', 'river_id')
        self.conf['var_discharge'] = self.conf.get('var_discharge', 'Q')

        # routing and solver options - time is validated at route time
        self.conf['routing'] = self.conf.get('routing', 'linear')
        assert self.conf['routing'] in ['linear', 'nonlinear'], 'Routing method not recognized'
        self.conf['positive_flow'] = self.conf.get('positive_flow', True)

        # check that the routing params files are given depending on the routing method
        if self.conf['routing'] == 'nonlinear':
            assert 'nonlinear_routing_params_file' in self.conf, 'Nonlinear routing requires nonlinear routing params'
            assert 'nonlinear_thresholds_file' in self.conf, 'Nonlinear routing requires nonlinear thresholds'
        else:
            assert 'routing_params_file' in self.conf, 'Linear routing requires linear routing params'

        # type and path checking on file paths
        if isinstance(self.conf['runoff_volumes_file'], str):
            self.conf['runoff_volumes_file'] = [self.conf['runoff_volumes_file'], ]
        if isinstance(self.conf['outflow_file'], str):
            self.conf['outflow_file'] = [self.conf['outflow_file'], ]
        for arg in [k for k in self.conf.keys() if 'file' in k]:
            if isinstance(self.conf[arg], list):
                self.conf[arg] = [os.path.abspath(path) for path in self.conf[arg]]
            elif isinstance(self.conf[arg], str):
                self.conf[arg] = os.path.abspath(self.conf[arg])

        # Enable logging options
        if self.conf['log']:
            LOG.disabled = False
            LOG.setLevel(self.conf.get('log_level', 'INFO'))
            log_destination = self.conf.get('log_stream', 'stdout')
            log_format = self.conf.get('log_format', '%(asctime)s - %(levelname)s - %(message)s')
            if log_destination == 'stdout':
                LOG.addHandler(logging.StreamHandler(sys.stdout))
            elif log_destination == 'stderr':
                LOG.addHandler(logging.StreamHandler(sys.stderr))
            elif isinstance(log_destination, str):
                LOG.addHandler(logging.FileHandler(log_destination))
            LOG.handlers[0].setFormatter(logging.Formatter(log_format))
            LOG.debug('Logger initialized')
        return

    def _validate_configs(self) -> None:
        LOG.info('Validating configs file')
        required_file_paths = ['connectivity_file',
                               'runoff_volumes_file',
                               'outflow_file', ]
        paths_should_exist = ['connectivity_file', ]

        if self.conf['routing'] == 'linear':
            required_file_paths.append('routing_params_file')
        else:
            required_file_paths.append('nonlinear_routing_params_file')
            required_file_paths.append('nonlinear_thresholds_file')

        for arg in required_file_paths:
            if arg not in self.conf:
                raise ValueError(f'{arg} not found in configs')
        for arg in paths_should_exist:
            if not os.path.exists(self.conf[arg]):
                raise FileNotFoundError(f'{arg} not found at given path')
        for path in self.conf['runoff_volumes_file']:
            assert os.path.exists(path), FileNotFoundError(f'runoff file not found at given path: {path}')

        return

    def _log_configs(self) -> None:
        LOG.debug(f'river-route version: {VERSION}')
        LOG.debug(f'Number of Rivers: {self._read_river_ids().shape[0]}')
        LOG.debug('Configs:')
        for k, v in self.conf.items():
            LOG.debug(f'\t{k}: {v}')
        return

    def _read_river_ids(self) -> np.array:
        """
        Reads river ids vector from parquet given in config file
        """
        return pd.read_parquet(self.conf['routing_params_file'], columns=[self.conf['var_river_id'], ]).values.flatten()

    def _read_linear_k(self) -> np.array:
        """
        Reads K vector from parquet given in config file
        """
        if hasattr(self, 'k'):
            return self.k
        self.k = pd.read_parquet(self.conf['routing_params_file'], columns=['k', ]).values.flatten()
        return self.k

    def _read_linear_x(self) -> np.array:
        """
        Reads X vector from parquet given in config file
        """
        if hasattr(self, 'x'):
            return self.x
        self.x = pd.read_parquet(self.conf['routing_params_file'], columns=['x', ]).values.flatten()
        return self.x

    def _read_initial_state(self) -> (np.array, np.array):
        if hasattr(self, 'initial_state'):
            return self.initial_state

        state_file = self.conf.get('initial_state_file', '')
        if state_file == '':
            LOG.debug('Setting initial state to zeros')
            return np.zeros(self.A.shape[0]), np.zeros(self.A.shape[0])
        assert os.path.exists(state_file), FileNotFoundError(f'Initial state file not found at: {state_file}')
        LOG.debug('Reading Initial State from Parquet')
        initial_state = pd.read_parquet(state_file).values
        initial_state = (initial_state[:, 0].flatten(), initial_state[:, 1].flatten())
        self.initial_state = initial_state
        return

    def _write_final_state(self) -> None:
        final_state_file = self.conf.get('final_state_file', '')
        if final_state_file == '':
            return
        LOG.debug('Writing Final State to Parquet')
        pd.DataFrame(self.initial_state, columns=['Q', 'R']).to_parquet(self.conf['final_state_file'])
        return

    def _get_digraph(self) -> nx.DiGraph:
        """
        Returns a directed graph of the river network
        """
        return connectivity_to_digraph(self.conf['connectivity_file'])

    def _set_adjacency_matrix(self) -> None:
        """
        Calculate the adjacency array from the connectivity file
        """
        if hasattr(self, 'A'):
            return
        LOG.debug('Calculating Network Adjacency Matrix (A)')
        self.A = connectivity_to_adjacency_matrix(self.conf['connectivity_file'])
        return

    def _set_lhs_matrix(self, c2: np.array = None) -> None:
        if c2 is None:
            c2 = self.c2
        self.lhs = scipy.sparse.eye(self.A.shape[0]) - (scipy.sparse.diags(self.c2) @ self.A)
        self.lhs = self.lhs.tocsc()
        return

    def _set_time_params(self, dates: np.array) -> None:
        """
        Set time parameters for the simulation
        """
        LOG.debug('Setting and validating time parameters')
        self.dt_runoff = self.conf.get('dt_runoff', (dates[1] - dates[0]).astype('timedelta64[s]').astype(int))
        self.dt_outflow = self.conf.get('dt_outflow', self.dt_runoff)
        self.dt_total = self.conf.get('dt_total', self.dt_runoff * dates.shape[0])
        if self.conf.get('dt_routing', 0):
            LOG.warning('dt_routing was not provided or is Null/False, defaulting to dt_runoff')
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
            LOG.error(e)
            raise AssertionError('Time options are not valid')

        # set derived datetime parameters for computation cycles later
        self.num_runoff_steps_per_outflow = int(self.dt_outflow / self.dt_runoff)
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)
        self.num_routing_steps = int(self.dt_total / self.dt_routing)
        return

    def _set_muskingum_coefficients(self, k: np.ndarray = None, x: np.ndarray = None) -> None:
        """
        Calculate the 3 MuskingumCunge routing coefficients for each segment using given k and x
        """
        LOG.debug('Calculating MuskingumCunge coefficients')

        if k is None:
            k = self._read_linear_k()
        if x is None:
            x = self._read_linear_x()

        dt_div_k = self.dt_routing / k
        denom = dt_div_k + (2 * (1 - x))
        _2x = 2 * x
        self.c1 = (dt_div_k + _2x) / denom
        self.c2 = (dt_div_k - _2x) / denom
        self.c3 = ((2 * (1 - x)) - dt_div_k) / denom

        # sum of coefficients should be 1 (or arbitrarily close) for all segments
        assert np.allclose(self.c1 + self.c2 + self.c3, 1), 'MuskingumCunge coefficients do not sum to 1'
        return

    def _set_nonlinear_muskingum_coefficients(self, q_t: np.array) -> None:
        if not hasattr(self, 'nonlinear_thresholds'):
            self.nonlinear_thresholds = pd.read_parquet(self.conf['nonlinear_thresholds_file'])
        threshold_columns = sorted(self.nonlinear_thresholds.columns)
        if not np.any(q_t > self.nonlinear_thresholds[threshold_columns[0]]):
            return

        if not hasattr(self, 'nonlinear_k_table'):
            self.nonlinear_k_table = pd.read_parquet(self.conf['nonlinear_routing_params_file'])

        k = self._read_linear_k()

        for threshold in threshold_columns:  # small to large
            values_to_change = (q_t > self.nonlinear_thresholds[threshold]).values
            k[values_to_change] = self.nonlinear_k_table[f'k_{threshold}'][values_to_change]

        self._set_muskingum_coefficients(k=k)
        self._set_lhs_matrix()
        return

    def route(self, **kwargs) -> 'MuskingumCunge':
        """
        Performs time-iterative runoff routing through the river network

        Args:
            **kwargs: optional keyword arguments to override and update previously calculated or given configs

        Returns:
            river_route.MuskingumCunge
        """
        LOG.info(f'Beginning routing: {self.conf["job_name"]}')
        t1 = datetime.datetime.now()

        if len(kwargs) > 0:
            LOG.info('Updating configs with kwargs')
            self.conf.update(kwargs)

        self._validate_configs()
        self._log_configs()
        self._set_adjacency_matrix()

        for runoff_file, outflow_file in zip(self.conf['runoff_volumes_file'], self.conf['outflow_file']):
            LOG.info('-' * 80)
            LOG.info(f'Reading runoff volumes file: {runoff_file}')
            with xr.open_dataset(runoff_file) as runoff_ds:
                LOG.debug('Reading time array')
                dates = runoff_ds['time'].values.astype('datetime64[s]')
                LOG.debug('Reading runoff array')
                runoffs = runoff_ds[self.conf['var_runoff_volume']].values

            self._set_time_params(dates)
            self._set_muskingum_coefficients()
            runoffs = runoffs / self.dt_runoff  # convert volume -> volume/time
            LOG.debug('Getting initial value arrays')
            q_t, r_t = self._read_initial_state()
            outflow_array = self._router(dates, runoffs, q_t, r_t)

            if self.dt_outflow > self.dt_runoff:
                LOG.info('Resampling dates and outflows to specified timestep')
                outflow_array = (
                    outflow_array
                    .reshape((
                        int(self.dt_total / self.dt_outflow),
                        int(self.dt_outflow / self.dt_runoff),
                        self.A.shape[0],
                    ))
                    .mean(axis=1)
                )

            LOG.info('Writing Outflow Array to File')
            outflow_array = np.round(outflow_array, decimals=2)
            self._write_outflows(outflow_file, dates, outflow_array)

        # write the final state to disc
        self._write_final_state()

        t2 = datetime.datetime.now()
        LOG.info('All runoff files routed')
        LOG.info(f'Total job time: {(t2 - t1).total_seconds()}')
        return self

    def calibrate(self, observed_df: pd.DataFrame, run_calibrated: bool = True) -> 'MuskingumCunge':
        """
        Calibrate K and X to given measured_values using optimization algorithms

        Args:
            observed_df: A pandas DataFrame with a datetime index, river id column names, and discharge values
            run_calibrated: whether to run the model with the calibrated parameters after optimization

        Returns:
            river_route.MuskingumCunge
        """
        LOG.info(f'Beginning optimization: {self.conf["job_name"]}')
        t1 = datetime.datetime.now()

        # todo nonlinear calibration
        if self.conf['routing'] == 'nonlinear':
            raise NotImplementedError('Nonlinear calibration not yet implemented')

        self._calibration_iteration_number = 0
        self._read_linear_k()
        self._read_linear_x()
        self._validate_configs()
        self._log_configs()
        self._set_adjacency_matrix()

        # find the indices of the rivers which have observed flows
        # todo create the smallest subgraph containing only the observed rivers for faster routing
        river_indices = self._read_river_ids()
        river_indices = [np.where(river_indices == r)[0][0] for r in observed_df.columns]
        objective = partial(self._calibration_objective, observed_df=observed_df, river_indices=river_indices)
        result = minimize_scalar(objective, method='bounded', bounds=(0.75, 1.25), options={'disp': True}, tol=1e-2)

        LOG.info('Optimization Results')
        LOG.info(result)
        self.k = self.k * result.x
        if self.conf['calibrated_linear_params_file']:
            LOG.info('Writing optimized k values to file')
            df = pd.read_parquet(self.conf['routing_params_file'])
            df['k'] = self.k
            df.to_parquet(self.conf['calibrated_linear_params_file'])

        t2 = datetime.datetime.now()
        LOG.info(f'Total job time: {(t2 - t1).total_seconds()}')

        if run_calibrated:
            LOG.info('Rerunning model with calibrated parameters')
            self.route()
        return self

    def _router(self, dates: np.array, runoffs: np.array, q_init: np.array, r_init: np.array, ) -> np.array:
        # set initial values
        self._set_lhs_matrix()
        outflow_array = np.zeros((runoffs.shape[0], self.A.shape[0]))
        q_t = np.zeros_like(q_init)
        q_t[:] = q_init
        r_prev = np.zeros_like(r_init)
        r_prev[:] = r_init

        LOG.info('Performing routing computation iterations')
        t1 = datetime.datetime.now()
        if self.conf['progress_bar']:
            dates = tqdm.tqdm(dates, desc='Runoff Routed')
        for runoff_time_step, runoff_end_date in enumerate(dates):
            r_t = runoffs[runoff_time_step, :]
            interval_flows = np.zeros((self.num_routing_steps_per_runoff, self.A.shape[0]))
            for routing_sub_iteration_num in range(self.num_routing_steps_per_runoff):
                if self.conf['routing'] == 'nonlinear':
                    self._set_nonlinear_muskingum_coefficients(q_t)
                rhs = (self.c1 * ((self.A @ q_t) + r_prev)) + (self.c2 * r_t) + (self.c3 * q_t)
                q_t[:] = self._solver(rhs)
                interval_flows[routing_sub_iteration_num, :] = q_t
            outflow_array[runoff_time_step, :] = np.mean(interval_flows, axis=0)
            r_prev[:] = r_t

        if self.conf['positive_flow']:
            outflow_array[outflow_array < 0] = 0

        t2 = datetime.datetime.now()
        LOG.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        self.initial_state = (q_t, r_prev)
        return outflow_array

    def _solver(self, rhs: np.array) -> np.array:
        return scipy.sparse.linalg.cgs(self.lhs, rhs, x0=rhs, atol=1, rtol=1e-2)[0]

    def _calibration_objective(self,
                               iteration_array: np.array,
                               observed_df: pd.DataFrame,
                               river_indices: np.array) -> np.float64:
        """
        Objective function for calibration

        Args:
            iteration_array: a scalar value to multiply the k values by
            river_indices: array of the indices of rivers with observations in the river id list from params

        Returns:
            np.float64
        """
        self._calibration_iteration_number += 1
        LOG.info(f'Iteration {self._calibration_iteration_number} - Testing scalar: {iteration_array}')
        iter_df = pd.DataFrame(columns=river_indices)
        LOG.disabled = True

        # set iteration k and x depending on the routing type in the configs and the size of the scalar
        k = self.k * iteration_array
        x = self.x

        self.conf['progress_bar'] = False
        for runoff_file, outflow_file in zip(self.conf['runoff_volumes_file'], self.conf['outflow_file']):
            with xr.open_dataset(runoff_file) as runoff_ds:
                dates = runoff_ds['time'].values.astype('datetime64[s]')
                runoffs = runoff_ds[self.conf['var_runoff_volume']].values
            self._set_time_params(dates)
            self._set_muskingum_coefficients(k=k, x=x)
            runoffs = runoffs / self.dt_runoff  # convert volume -> volume/time
            q_t, r_t = self._read_initial_state()
            outflow_array = self._router(dates, runoffs, q_t, r_t)
            if iter_df.empty:
                iter_df = pd.DataFrame(outflow_array, index=dates)[river_indices]
            else:
                iter_df = pd.concat([iter_df, pd.DataFrame(outflow_array, index=dates)[river_indices]])
        LOG.disabled = False
        self.conf['progress_bar'] = True
        # todo clip by date ranges, merge dataframes, fill with nans, calculate error metric
        mse = np.sum(np.mean((iter_df.values - observed_df.values) ** 2, axis=0))
        LOG.info(f'Iteration {self._calibration_iteration_number} - MSE: {mse}')
        return mse

    def _write_outflows(self, outflow_file: str, dates: np.array, outflow_array: np.array) -> None:
        reference_date = datetime.datetime.fromtimestamp(dates[0].astype(int), tz=datetime.timezone.utc)
        dates = dates[::self.num_runoff_steps_per_outflow].astype('datetime64[s]')
        dates = dates - dates[0]

        with nc.Dataset(outflow_file, mode='w', format='NETCDF4') as ds:
            ds.createDimension('time', size=dates.shape[0])
            ds.createDimension(self.conf['var_river_id'], size=outflow_array.shape[1])

            time_var = ds.createVariable('time', 'f8', ('time',))
            time_var.units = f'seconds since {reference_date.strftime("%Y-%m-%d %H:%M:%S")}'
            time_var[:] = dates

            id_var = ds.createVariable(self.conf['var_river_id'], 'i4', (self.conf['var_river_id']), )
            id_var[:] = self._read_river_ids()

            flow_var = ds.createVariable(self.conf['var_discharge'], 'f4', ('time', self.conf['var_river_id']))
            flow_var[:] = outflow_array
            flow_var.long_name = 'Discharge at catchment outlet'
            flow_var.standard_name = 'discharge'
            flow_var.aggregation_method = 'mean'
            flow_var.units = 'm3 s-1'
        return

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
            G = connectivity_to_digraph(self.conf['connectivity_file'])
            ancestors = set(list(nx.ancestors(G, river_id)) + [river_id, ])
        with xr.open_mfdataset(self.conf['runoff_volumes_file']) as ds:
            vdf = (
                ds
                .sel(**{self.conf['var_river_id']: ancestors})
                [self.conf['var_runoff_volume']]
                .to_dataframe()
                [[self.conf['var_runoff_volume'], ]]
                .reset_index()
                .pivot(index='time', columns=self.conf['var_river_id'], values=self.conf['var_runoff_volume'])
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
            LOG.warning(f'More discharge than runoff volume for river {river_id}')

        return df

    def save_configs(self, path: str) -> None:
        """
        Save the current configs of the class to a json file
        Args:
            path: the file path where the json will be written

        Returns:
            None
        """
        with open(path, 'w') as f:
            json.dump(self.conf, f)
        return
