import datetime
import json
import logging
import os
import sys

import netCDF4 as nc
import networkx as nx
import numpy as np
import pandas as pd
import scipy
import tqdm
import xarray as xr
import yaml

from ._meta import __version__ as VERSION
from .tools import connectivity_to_adjacency_matrix
from .tools import connectivity_to_digraph

__all__ = ['MuskingumCunge', ]

LOG = logging.getLogger('river_route')
LOG.disabled = True


class MuskingumCunge:
    # Given configs
    conf: dict
    log: logging.Logger

    # Routing matrices
    A: scipy.sparse.csc_matrix
    lhs: scipy.sparse.csc_matrix
    lhsinv: scipy.sparse.csc_matrix
    c1: np.array
    c2: np.array
    c3: np.array
    qinit: np.array
    rinit: np.array

    # Time options
    dt_total: float
    dt_runoff: float
    dt_outflow: float
    dt_routing: float
    num_routing_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_outflow: int

    def __init__(self, config_file: str = None, **kwargs, ):
        """
        Implements Matrix MuskingumCunge routing
        """
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
        self.conf['runoff_volume_var'] = self.conf.get('runoff_volume_var', 'm3_riv')
        self.conf['dt_routing'] = self.conf.get('dt_routing', 60)
        self.conf['log_level'] = self.conf.get('log_level', 'INFO')
        self.conf['min_q'] = self.conf.get('min_q', False)
        self.conf['max_q'] = self.conf.get('max_q', False)

        # type and path checking on file paths
        if isinstance(self.conf['runoff_file'], str):
            self.conf['runoff_file'] = [self.conf['runoff_file'], ]
        if isinstance(self.conf['outflow_file'], str):
            self.conf['outflow_file'] = [self.conf['outflow_file'], ]
        for arg in [k for k in self.conf.keys() if 'file' in k]:
            if isinstance(self.conf[arg], list):
                self.conf[arg] = [os.path.abspath(path) for path in self.conf[arg]]
            elif isinstance(self.conf[arg], str):
                self.conf[arg] = os.path.abspath(self.conf[arg])

        if self.conf['log']:
            LOG.disabled = False
            LOG.setLevel(self.conf.get('log_level', 'DEBUG'))
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
                               'runoff_file',
                               'outflow_file', ]
        paths_should_exist = ['connectivity_file', ]
        required_time_opts = ['dt_routing', ]

        for arg in required_file_paths + required_time_opts:
            if arg not in self.conf:
                raise ValueError(f'{arg} not found in configs')
        for arg in paths_should_exist:
            if not os.path.exists(self.conf[arg]):
                raise FileNotFoundError(f'{arg} not found at given path')
        for path in self.conf['runoff_file']:
            assert os.path.exists(path), FileNotFoundError(f'runoff file not found at given path: {path}')

        return

    def _log_configs(self) -> None:
        LOG.debug(f'river-route version: {VERSION}')
        LOG.debug('Configs:')
        for k, v in self.conf.items():
            LOG.debug(f'\t{k}: {v}')
        return

    def _read_riverids(self) -> np.array:
        """
        Reads river ids vector from parquet given in config file
        """
        return pd.read_parquet(self.conf['routing_params_file'], columns=['rivid', ]).values.flatten()

    def _read_k(self) -> np.array:
        """
        Reads K vector from parquet given in config file
        """
        return pd.read_parquet(self.conf['routing_params_file'], columns=['k', ]).values.flatten()

    def _read_x(self) -> np.array:
        """
        Reads X vector from parquet given in config file
        """
        return pd.read_parquet(self.conf['routing_params_file'], columns=['x', ]).values.flatten()

    def _read_qinit(self) -> np.array:
        if hasattr(self, 'qinit'):
            return self.qinit

        qinit = self.conf.get('qinit_file', None)
        if qinit is None or qinit == '':
            return np.zeros(self.A.shape[0])
        return pd.read_parquet(self.conf['qinit_file']).values.flatten()

    def _read_rinit(self) -> np.array:
        if hasattr(self, 'rinit'):
            return self.rinit

        rinit = self.conf.get('rinit_file', None)
        if rinit is None or rinit == '':
            return np.zeros(self.A.shape[0])
        return pd.read_parquet(self.conf['rinit_file']).values.flatten()

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

        if os.path.exists(self.conf.get('adj_file', '')):
            LOG.info('Loading adjacency matrix from file')
            self.A = scipy.sparse.load_npz(self.conf['adj_file'])
            return

        LOG.info('Calculating Network Adjacency Matrix (A)')
        self.A = connectivity_to_adjacency_matrix(self.conf['connectivity_file'])
        if self.conf.get('adj_file', ''):
            LOG.info('Saving adjacency matrix to file')
            scipy.sparse.save_npz(self.conf['adj_file'], self.A)
        return

    def _set_lhs_matrix(self) -> None:
        """
        Calculate the LHS matrix for the routing problem
        """
        if hasattr(self, 'lhs'):
            return

        if os.path.exists(self.conf.get('lhs_file', '')):
            LOG.info('Loading LHS matrix from file')
            self.lhs = scipy.sparse.load_npz(self.conf['lhs_file'])
            return

        LOG.info('Calculating LHS Matrix')
        self.lhs = scipy.sparse.eye(self.A.shape[0]) - scipy.sparse.diags(self.c2) @ self.A
        self.lhs = self.lhs.tocsc()
        if self.conf.get('lhs_file', ''):
            LOG.info('Saving LHS matrix to file')
            scipy.sparse.save_npz(self.conf['lhs_file'], self.lhs)
        return

    def _set_lhs_inv_matrix(self) -> None:
        """
        Calculate the LHS matrix for the routing problem
        """
        if hasattr(self, 'lhsinv'):
            return

        if os.path.exists(self.conf.get('lhsinv_file', '')):
            LOG.info('Loading LHS Inverse matrix from file')
            self.lhsinv = scipy.sparse.load_npz(self.conf['lhsinv_file'])
            return

        self._set_lhs_matrix()
        LOG.info('Inverting LHS Matrix')
        self.lhsinv = scipy.sparse.csc_matrix(scipy.sparse.linalg.inv(self.lhs))
        if self.conf.get('lhsinv_file', ''):
            LOG.info('Saving LHS Inverse matrix to file')
            scipy.sparse.save_npz(self.conf['lhsinv_file'], self.lhsinv)
        return

    def _set_time_params(self, dates: np.array):
        """
        Set time parameters for the simulation
        """
        LOG.info('Setting and validating time parameters')
        self.dt_routing = self.conf['dt_routing']
        self.dt_runoff = self.conf.get('dt_runoff', (dates[1] - dates[0]).astype('timedelta64[s]').astype(int))
        self.dt_outflow = self.conf.get('dt_outflow', self.dt_runoff)
        self.dt_total = self.conf.get('dt_total', self.dt_runoff * dates.shape[0])

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

    def _calculate_muskingum_coefficients(self, k: np.ndarray = None, x: np.ndarray = None) -> None:
        """
        Calculate the 3 MuskingumCunge routing coefficients for each segment using given k and x
        """
        LOG.info('Calculating MuskingumCunge coefficients')

        if k is None:
            k = self._read_k()
        if x is None:
            x = self._read_x()

        dt_div_k = self.conf['dt_routing'] / k
        denom = dt_div_k + (2 * (1 - x))
        self.c1 = (dt_div_k + (2 * x)) / denom
        self.c2 = (dt_div_k - (2 * x)) / denom
        self.c3 = ((2 * (1 - x)) - dt_div_k) / denom

        # sum of muskingum coefficients should be 1 for all segments
        assert np.allclose(self.c1 + self.c2 + self.c3, 1), 'MuskingumCunge coefficients do not approximately sum to 1'
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
        self._calculate_muskingum_coefficients()

        for runoff_file, outflow_file in zip(self.conf['runoff_file'], self.conf['outflow_file']):
            LOG.info(f'Reading runoff volumes file: {runoff_file}')
            with xr.open_dataset(runoff_file) as runoff_ds:
                LOG.info('Reading time array')
                dates = runoff_ds['time'].values.astype('datetime64[s]')
                LOG.info('Reading runoff array')
                runoffs = runoff_ds[self.conf['runoff_volume_var']].values

            self._set_time_params(dates)
            runoffs = runoffs / self.dt_runoff  # convert volume -> volume/time

            LOG.debug('Getting initial value arrays')
            q_t = self._read_qinit()
            r_t = self._read_rinit()

            # begin the analytical solution
            if self.conf['routing'] == 'linear':
                outflow_array = self._solver_linear(dates, runoffs, q_t, r_t)
            elif self.conf['routing'] == 'nonlinear':
                outflow_array = self._solver_nonlinear(dates, runoffs, q_t, r_t)
            else:
                raise ValueError('Routing method not recognized')

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
        if self.conf.get('qfinal_file', False):
            LOG.info('Writing Qfinal parquet')
            pd.DataFrame(self.qinit, columns=['Q', ]).astype(float).to_parquet(self.conf['qfinal_file'])
        if self.conf.get('rfinal_file', False):
            LOG.info('Writing Rfinal parquet')
            pd.DataFrame(self.rinit, columns=['R', ]).astype(float).to_parquet(self.conf['rfinal_file'])

        t2 = datetime.datetime.now()
        LOG.info('All runoff files routed')
        LOG.info(f'Total job time: {(t2 - t1).total_seconds()}')
        return self

    def _solver_linear(self,
                       dates: np.array,
                       runoffs: np.array,
                       q_init: np.array,
                       r_init: np.array, ) -> np.array:
        self._set_lhs_inv_matrix()
        outflow_array = np.zeros((runoffs.shape[0], self.A.shape[0]))

        # set initial values
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
            for routing_substep_iteration in range(self.num_routing_steps_per_runoff):
                q_t = self.lhsinv @ ((self.c1 * ((self.A @ q_t) + r_prev)) + (self.c2 * r_t) + (self.c3 * q_t))
                interval_flows[routing_substep_iteration, :] = q_t
            r_prev[:] = r_t
            interval_flows = np.round(np.mean(interval_flows, axis=0), decimals=2)
            outflow_array[runoff_time_step, :] = interval_flows
        t2 = datetime.datetime.now()
        LOG.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')

        self.qinit = q_t
        self.rinit = r_prev
        return outflow_array

    def _solver_nonlinear(self,
                          dates: np.array,
                          runoffs: np.array,
                          q_init: np.array,
                          r_init: np.array, ) -> np.array:
        self._set_lhs_matrix()
        outflow_array = np.zeros((runoffs.shape[0], self.A.shape[0]))

        # set initial values
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
            for routing_substep_iteration in range(self.num_routing_steps_per_runoff):
                rhs = ((self.c1 * ((self.A @ q_t) + r_prev)) + (self.c2 * r_t) + (self.c3 * q_t))
                q_t[:] = scipy.sparse.linalg.cgs(self.lhs, rhs, x0=q_t)[0]
                interval_flows[routing_substep_iteration, :] = q_t
            outflow_array[runoff_time_step, :] = np.mean(interval_flows, axis=0)
            r_prev[:] = r_t

        t2 = datetime.datetime.now()
        LOG.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        self.qinit = q_t
        self.rinit = r_prev
        return outflow_array

    def _write_outflows(self, outflow_file: str, dates: np.array, outflow_array: np.array) -> None:
        reference_date = datetime.datetime.fromtimestamp(dates[0].astype(int), tz=datetime.timezone.utc)
        dates = dates[::self.num_runoff_steps_per_outflow].astype('datetime64[s]')
        dates = dates - dates[0]

        with nc.Dataset(outflow_file, mode='w', format='NETCDF4') as ds:
            ds.createDimension('time', size=dates.shape[0])
            ds.createDimension('rivid', size=self.A.shape[0])

            ds.createVariable('time', 'f8', ('time',))
            ds['time'].units = f'seconds since {reference_date.strftime("%Y-%m-%d %H:%M:%S")}'
            ds['time'][:] = dates

            ds.createVariable('rivid', 'i4', ('rivid',))
            ds['rivid'][:] = self._read_riverids()

            ds.createVariable('Qout', 'f4', ('time', 'rivid'))
            ds['Qout'][:] = outflow_array
            ds['Qout'].long_name = 'Discharge at the outlet of each river reach'
            ds['Qout'].standard_name = 'discharge'
            ds['Qout'].aggregation_method = 'mean'
            ds['Qout'].units = 'm3 s-1'
        return

    def hydrograph(self, rivid: int) -> pd.DataFrame:
        """
        Get the hydrograph for a given river id as a pandas dataframe

        Args:
            rivid: the ID of a river reach in the output files

        Returns:
            pandas.DataFrame
        """
        with xr.open_mfdataset(self.conf['outflow_file']) as ds:
            df = ds.Qout.sel(rivid=rivid).to_dataframe()[['Qout', ]]
            df.columns = [rivid, ]
        return df

    def mass_balance(self, rivid: int, ancestors: list = None) -> pd.DataFrame:
        """
        Get the mass balance for a given river id as a pandas dataframe

        Args:
            rivid: the ID of a river reach in the output files
            ancestors: a list of the given rivid and all rivers upstream of that river

        Returns:
            pandas.DataFrame
        """
        if type(rivid) is not int:
            raise TypeError(f'rivid should be an integer ID of a river to mass balance')
        if ancestors is None:
            G = connectivity_to_digraph(self.conf['connectivity_file'])
            ancestors = list(nx.ancestors(G, rivid))
        with xr.open_mfdataset(self.conf['runoff_file']) as ds:
            vdf = (
                ds
                .sel(rivid=ancestors)
                .m3_riv
                .to_dataframe()
                [['m3_riv', ]]
                .reset_index()
                .set_index('time')
                .pivot(columns='rivid', values='m3_riv')
                .sum(axis=1)
                .cumsum()
                .rename('runoff_volume')
            )
        with xr.open_mfdataset(self.conf['outflow_file']) as ds:
            qdf = ds.sel(rivid=rivid).to_dataframe()[['Qout', ]].cumsum()
            # convert to discharged volume - multiply by the time delta in seconds
            qdf = qdf * (qdf.index[1] - qdf.index[0]).total_seconds()
            qdf.columns = ['discharge_volume', ]

        df = qdf.join(vdf)
        df['runoff-discharge'] = df['runoff_volume'] - df['discharge_volume']
        if not df['runoff-discharge'].gt(0).all():
            LOG.warning(f'More discharge than runoff volume for river {rivid}')

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