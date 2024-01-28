import datetime
import json
import logging
import os
import sys

import matplotlib.pyplot as plt
import netCDF4 as nc
import networkx as nx
import numpy as np
import pandas as pd
import scipy
import tqdm
import xarray as xr
import yaml

from .tools import connectivity_to_adjacency_matrix
from .tools import connectivity_to_digraph

__all__ = ['Muskingum', ]


class Muskingum:
    # Given configs
    conf: dict

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
    num_timesteps_resample: int

    def __init__(self, config_file: str = None, **kwargs, ):
        """
        Implements Matrix Muskingum routing
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
        self.conf['job_name'] = self.conf.get('job_name', 'untitled_job')
        self.conf['progress_bar'] = self.conf.get('progress_bar', True)
        self.conf['runoff_volume_var'] = self.conf.get('runoff_volume_var', 'm3_riv')
        self.conf['dt_routing'] = self.conf.get('dt_routing', 300)
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

        # start a logger
        log_basic_configs = {
            'stream': sys.stdout,
            'level': self.conf['log_level'],
            'format': '%(asctime)s - %(levelname)s - %(message)s',
        }
        if self.conf.get('log_file', ''):
            log_basic_configs['filename'] = self.conf['log_file']
            log_basic_configs['filemode'] = 'w'
            log_basic_configs.pop('stream')
        logging.basicConfig(**log_basic_configs)
        return

    def _validate_configs(self) -> None:
        logging.info('Validating configs file')
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
        logging.debug('Configs:')
        for k, v in self.conf.items():
            logging.debug(f'\t{k}: {v}')
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
            logging.info('Loading adjacency matrix from file')
            self.A = scipy.sparse.load_npz(self.conf['adj_file'])
            return

        logging.info('Calculating Network Adjacency Matrix (A)')
        self.A = connectivity_to_adjacency_matrix(self.conf['connectivity_file'])
        if self.conf.get('adj_file', ''):
            logging.info('Saving adjacency matrix to file')
            scipy.sparse.save_npz(self.conf['adj_file'], self.A)
        return

    def _set_lhs_matrix(self) -> None:
        """
        Calculate the LHS matrix for the routing problem
        """
        if hasattr(self, 'lhs'):
            return

        if os.path.exists(self.conf.get('lhs_file', '')):
            logging.info('Loading LHS matrix from file')
            self.lhs = scipy.sparse.load_npz(self.conf['lhs_file'])
            return

        logging.info('Calculating LHS Matrix')
        self.lhs = scipy.sparse.eye(self.A.shape[0]) - scipy.sparse.diags(self.c2) @ self.A
        self.lhs = self.lhs.tocsc()
        if self.conf.get('lhs_file', ''):
            logging.info('Saving LHS matrix to file')
            scipy.sparse.save_npz(self.conf['lhs_file'], self.lhs)
        return

    def _set_lhs_inv_matrix(self) -> None:
        """
        Calculate the LHS matrix for the routing problem
        """
        if hasattr(self, 'lhsinv'):
            return

        if os.path.exists(self.conf.get('lhsinv_file', '')):
            logging.info('Loading LHS Inverse matrix from file')
            self.lhsinv = scipy.sparse.load_npz(self.conf['lhsinv_file'])
            return

        self._set_lhs_matrix()
        logging.info('Inverting LHS Matrix')
        self.lhsinv = scipy.sparse.csc_matrix(scipy.sparse.linalg.inv(self.lhs))
        if self.conf.get('lhsinv_file', ''):
            logging.info('Saving LHS Inverse matrix to file')
            scipy.sparse.save_npz(self.conf['lhsinv_file'], self.lhsinv)
        return

    def _set_time_params(self, dates: np.array):
        """
        Set time parameters for the simulation
        """
        logging.info('Setting and validating time parameters')
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
            logging.error(e)
            raise AssertionError('Time options are not valid')

        # set derived datetime parameters for computation cycles later
        self.num_runoff_steps_per_outflow = int(self.dt_outflow / self.dt_runoff)
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)
        self.num_routing_steps = int(self.dt_total / self.dt_routing)
        self.num_timesteps_resample = int(self.dt_outflow / self.dt_runoff)
        return

    def _calculate_muskingum_coefficients(self) -> None:
        """
        Calculate the 3 Muskingum Cunge routing coefficients for each segment using given k and x
        """
        logging.info('Calculating Muskingum coefficients')

        k = self._read_k()
        x = self._read_x()

        dt_route_half = self.conf['dt_routing'] / 2
        kx = k * x
        denom = k - kx + dt_route_half

        self.c1 = (dt_route_half - kx) / denom
        self.c2 = (dt_route_half + kx) / denom
        self.c3 = (k - kx - dt_route_half) / denom

        # sum of muskingum coefficients should be 1 for all segments
        assert np.allclose(self.c1 + self.c2 + self.c3, 1), 'Muskingum coefficients do not approximately sum to 1'
        return

    def route(self, **kwargs) -> 'Muskingum':
        """
        Performs time-iterative runoff routing through the river network

        Args:
            **kwargs: optional keyword arguments to override and update previously calculated or given configs

        Returns:
            river_route.Muskingum
        """
        logging.info(f'Beginning routing: {self.conf["job_name"]}')
        t1 = datetime.datetime.now()

        if len(kwargs) > 0:
            logging.info('Updating configs with kwargs')
            self.conf.update(kwargs)

        self._validate_configs()
        self._log_configs()
        self._set_adjacency_matrix()
        self._calculate_muskingum_coefficients()

        for runoff_file, outflow_file in zip(self.conf['runoff_file'], self.conf['outflow_file']):
            logging.info(f'Reading runoff volumes file: {runoff_file}')
            with xr.open_dataset(runoff_file) as runoff_ds:
                logging.info('Reading time array')
                dates = runoff_ds['time'].values.astype('datetime64[s]')
                self._set_time_params(dates)
                logging.info(f'Reading runoff array: {runoff_file}')
                runoffs = runoff_ds[self.conf['runoff_volume_var']].values
                runoffs[runoffs < 0] = np.nan
                runoffs = np.nan_to_num(runoffs, nan=0.0)
                runoffs = runoffs / self.dt_runoff  # volume to volume/time

            logging.debug('Getting initial value arrays')
            q_t = self._read_qinit()
            r_t = self._read_rinit()
            inflow_t = (self.A @ q_t) + r_t

            # begin the analytical solution
            self._set_lhs_inv_matrix()
            outflow_array = self._analytical_solution(dates, runoffs, q_t, r_t, inflow_t)

            if self.dt_outflow > self.dt_runoff:
                logging.info('Resampling dates and outflows to specified timestep')
                outflow_array = (
                    outflow_array
                    .reshape((
                        int(self.dt_total / self.dt_outflow),
                        int(self.dt_outflow / self.dt_runoff),
                        self.A.shape[0],
                    ))
                    .mean(axis=1)
                )

            logging.info('Writing Outflow Array to File')
            outflow_array = np.round(outflow_array, decimals=2)
            self._write_outflows(outflow_file, dates, outflow_array)

        # write the final state to disc
        if self.conf.get('qfinal_file', False):
            logging.info('Writing Qfinal parquet')
            pd.DataFrame(self.qinit, columns=['Q', ]).astype(float).to_parquet(self.conf['qfinal_file'])
        if self.conf.get('rfinal_file', False):
            logging.info('Writing Rfinal parquet')
            pd.DataFrame(self.rinit, columns=['R', ]).astype(float).to_parquet(self.conf['rfinal_file'])

        t2 = datetime.datetime.now()
        logging.info('All runoff files routed')
        logging.info(f'Total job time: {(t2 - t1).total_seconds()}')
        return self

    def _analytical_solution(self,
                             dates: np.array,
                             runoffs: np.array,
                             q_t: np.array,
                             r_t: np.array,
                             inflow_t: np.array, ) -> np.array:
        outflow_array = np.zeros((runoffs.shape[0], self.A.shape[0]))

        force_min_max = self.conf['min_q'] or self.conf['max_q']

        logging.info('Performing routing computation iterations')
        t1 = datetime.datetime.now()
        if self.conf['progress_bar']:
            dates = tqdm.tqdm(dates, desc='Runoff Routed')
        for inflow_time_step, inflow_end_date in enumerate(dates):
            r_t = runoffs[inflow_time_step, :]
            interval_flows = np.zeros((self.num_routing_steps_per_runoff, self.A.shape[0]))
            for routing_substep_iteration in range(self.num_routing_steps_per_runoff):
                inflow_tnext = (self.A @ q_t) + r_t
                q_t = self.lhsinv @ ((self.c1 * inflow_t) + (self.c2 * r_t) + (self.c3 * q_t))
                q_t = q_t if not force_min_max else np.clip(q_t, a_min=self.conf['min_q'], a_max=self.conf['max_q'])
                interval_flows[routing_substep_iteration, :] = q_t
                inflow_t = inflow_tnext
            interval_flows = np.mean(interval_flows, axis=0)
            interval_flows = np.round(interval_flows, decimals=2)
            outflow_array[inflow_time_step, :] = interval_flows
        t2 = datetime.datetime.now()
        logging.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')

        self.qinit = q_t
        self.rinit = r_t
        return outflow_array

    def _write_outflows(self, outflow_file: str, dates: np.array, outflow_array: np.array) -> None:
        reference_date = datetime.datetime.fromtimestamp(dates[0].astype(int), tz=datetime.timezone.utc)
        dates = dates[::self.num_timesteps_resample].astype('datetime64[s]')
        dates = dates - dates[0]

        with nc.Dataset(outflow_file, mode='w') as ds:
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

    def plot(self, rivid: int, show: bool = False) -> plt.Figure:
        """
        Create a plot of the hydrograph for a given river id with matplotlib and optionally display it

        Args:
            rivid: the ID of a river reach in the output files
            show: boolean flag to display the figure or not

        Returns:
            matplotlib.pyplot.Figure
        """
        with xr.open_mfdataset(self.conf['outflow_file']) as ds:
            figure = ds['Qout'].sel(rivid=rivid).to_dataframe()['Qout'].plot()
        if show:
            plt.show()
        return figure

    def hydrograph_to_csv(self, rivid: int, csv_path: str = None) -> None:
        """
        Save the hydrograph for a given river id to a csv file

        Args:
            rivid: the ID of a river reach in the output files
            csv_path: the file path where the csv will be written

        Returns:
            None
        """
        with xr.open_mfdataset(self.conf['outflow_file']) as ds:
            df = ds.sel(rivid=rivid).to_dataframe()[['Qout', ]]
            df.columns = [rivid, ]
        df.to_csv(csv_path)
        return

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
