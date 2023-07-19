import json
import logging
import os
import sys

import numpy as np
import pandas as pd
import scipy
import xarray as xr
import yaml
from numba import float64, int64, types
from scipy.sparse.linalg import bicgstab

numba_class_specs = [
    ('sim_config_file', types.unicode_type),
    ('config', types.DictType(types.unicode_type, types.unicode_type)),
    ('N', float64[:]),
    ('c1', float64[:]),
    ('c2', float64[:]),
    ('c3', float64[:]),
    ('num_outflow_steps', int64),
    ('num_routing_substeps_per_outflow', int64),
    ('flow_temporal_aggregation_function', types.FunctionType),
]


class RouteMuskingum:
    sim_config_file: str

    configs: dict

    N: np.array
    c1: np.array
    c2: np.array
    c3: np.array

    num_outflow_steps: int
    num_routing_substeps_per_outflow: int

    flow_temporal_aggregation_function: callable

    def __init__(self,
                 sim_config_file: str,
                 **kwargs, ) -> None:
        """
        Read config files and prepare to execute simulation

        Args:
            sim_config_file (str): path to the simulation config file
        """
        self.sim_config_file = sim_config_file
        self.read_and_validate_configs()
        self.calculate_connections_array()
        self.calculate_muskingum_coefficients()
        self.define_flow_temporal_aggregation_function()
        return

    def read_and_validate_configs(self) -> None:
        """
        Validate simulation configs
        """
        if self.sim_config_file.endswith('.json'):
            with open(self.sim_config_file, 'r') as f:
                self.configs = json.load(f)
        elif self.sim_config_file.endswith('.yml') or self.sim_config_file.endswith('.yaml'):
            with open(self.sim_config_file, 'r') as f:
                self.configs = yaml.load(f, Loader=yaml.FullLoader)
        else:
            raise RuntimeError('Unrecognized simulation config file type')

        # check that required file paths are given and exist
        required_file_paths = ['routing_params_file',
                               'connectivity_file',
                               'inflow_file',
                               'outflow_file', ]
        required_time_opts = ['dt_total',
                              'dt_inflows',
                              'dt_outflows',
                              'dt_routing', ]
        optional_file_paths = ['qinit_file',
                               'qfinal_file',
                               'log_file', ]
        paths_should_exist = ['routing_params_file',
                              'connectivity_file',
                              'inflow_file', ]
        for arg in required_file_paths + required_time_opts:
            if arg not in self.configs:
                raise ValueError(f'{arg} not found in config file')
        for arg in paths_should_exist:
            if not os.path.exists(self.configs[arg]):
                raise FileNotFoundError(f'{arg} not found at given path')

        # start a logger
        log_basic_configs = {
            'stream': sys.stdout,
            'level': logging.INFO,
            'format': '%(asctime)s %(levelname)s %(message)s',
        }
        if self.configs.get('log_file', ''):
            log_basic_configs['filename'] = self.configs['log_file']
            log_basic_configs['filemode'] = 'w'
            log_basic_configs.pop('stream')
        logging.basicConfig(**log_basic_configs)

        # check that time options have the correct sizes
        assert self.configs['dt_total'] >= self.configs['dt_inflows'], 'dt_total !>= dt_inflows'
        assert self.configs['dt_inflows'] >= self.configs['dt_outflows'], 'dt_inflows !>= dt_outflows'
        assert self.configs['dt_outflows'] >= self.configs['dt_routing'], 'dt_outflows !>= dt_routing'

        # check that time options are evenly divisible
        assert self.configs['dt_total'] % self.configs['dt_inflows'] == 0, \
            'dt_total should be a whole number multiple of dt_inflows'
        assert self.configs['dt_total'] % self.configs['dt_outflows'] == 0, \
            'dt_total should be a whole number multiple of dt_outflows'
        assert self.configs['dt_inflows'] % self.configs['dt_outflows'] == 0, \
            'dt_inflows should be a whole number multiple of dt_outflows'
        assert self.configs['dt_outflows'] % self.configs['dt_routing'] == 0, \
            'dt_outflows should be a whole number multiple of dt_routing'

        # set derived datetime parameters for computation cycles later
        self.num_outflow_steps = int(self.configs.get('dt_total') / self.configs.get('dt_outflows'))
        self.num_routing_substeps_per_outflow = int(self.configs.get('dt_outflows') / self.configs.get('dt_routing'))

        return

    def calculate_muskingum_coefficients(self) -> None:
        """
        Calculate the 3 Muskingum Cunge routing coefficients for each segment using given k and x
        """
        logging.info('Calculating Muskingum coefficients')

        k = self.get_k()
        x = self.get_x()

        dt_route_half = np.ones(k.shape[0]) * self.configs['dt_routing'] / 2
        denom = k * (1 - x) + dt_route_half

        self.c1 = (dt_route_half - k * x) / denom
        self.c2 = (dt_route_half + k * x) / denom
        self.c3 = k * (1 - x) - dt_route_half
        return

    def calculate_connections_array(self) -> None:
        """
        Calculate the connections array from the connectivity file
        """
        logging.info('Calculating Connectivity array')

        df = pd.read_parquet(self.configs['connectivity_file'])
        self.N = np.zeros((df.values.shape[0], df.values.shape[0]))
        for idx, row in df.iterrows():
            upstreams = (
                pd
                .Series(row.values[2:])
                .replace(0, np.nan)
                .dropna()
                .apply(lambda x: df.loc[df['id'] == x, 'topo_order'].values[0])
                .values
            )
            if len(upstreams):
                self.N[idx, upstreams] = 1
        self.N = scipy.sparse.csc_matrix(self.N.astype(np.float32))
        return

    def define_flow_temporal_aggregation_function(self) -> None:
        """
        Define the function to use for temporal aggregation of flows calculated during the routing substeps between
        the output time steps

        Expects a config parameter 'aggregation_method' to be set to one of the following:
            'mean' - average of all substep flows
            'max' - maximum of all substep flows
            'min' - minimum of all substep flows
            'instantaneous' - use the flow from the last substep
        """

        if self.configs['aggregation_method'] == 'mean':
            self.flow_temporal_aggregation_function = np.mean
        elif self.configs['aggregation_method'] == 'max':
            self.flow_temporal_aggregation_function = np.max
        elif self.configs['aggregation_method'] == 'min':
            self.flow_temporal_aggregation_function = np.min
        elif self.configs['aggregation_method'] == 'instantaneous':
            self.flow_temporal_aggregation_function = lambda x: x[-1]
        else:
            # todo check that the aggregation method is valid in read_validate_configs
            raise RuntimeError('Unrecognized aggregation method bypassed config file validation')

    def get_k(self) -> np.array:
        """
        Reads K vector from parquet given in config file
        """
        return pd.read_parquet(self.configs['routing_params_file'], columns=['k', ]).values.flatten()

    def get_x(self) -> np.array:
        """
        Reads X vector from parquet given in config file
        """
        return pd.read_parquet(self.configs['routing_params_file'], columns=['x', ]).values.flatten()
        # return np.ones(self.N.shape[0]) * 0.25

    def get_riverids(self) -> np.array:
        """
        Reads riverids vector from parquet given in config file
        """
        return pd.read_parquet(self.configs['routing_params_file'], columns=['id', ]).values.flatten()

    def qinit(self) -> np.array:
        qinit = self.configs.get('qinit_file', None)
        if qinit is None or qinit == '':
            return np.zeros(self.N.shape[0])
        return pd.read_parquet(self.configs['qinit_file']).values.flatten()

    def route(self):
        # get the list of dates for runoff data
        logging.info('Reading Inflow Data')
        with xr.open_dataset(self.configs['inflow_file']) as inflow_ds:
            dates = inflow_ds['time'].values
            inflows = inflow_ds['m3_riv'].values
            inflows.fill(0)

        logging.info('Scaffolding Outflow File')
        self.scaffold_outflow_file(dates)

        logging.info('Beginning Computation loops')
        q_t = self.qinit()
        lhs = scipy.sparse.csc_matrix(np.diag((1 - self.c1).flatten()))
        # preconditioner = scipy.sparse.linalg.spilu(lhs)

        outflow_array = np.zeros((self.num_outflow_steps, self.N.shape[0]))
        for inflow_time_step, inflow_end_date in enumerate(dates):
            q_ro = inflows[inflow_time_step, :] / self.num_routing_substeps_per_outflow
            interval_flows = np.zeros((self.num_routing_substeps_per_outflow, self.N.shape[0]))
            for routing_substep_iteration in range(self.num_routing_substeps_per_outflow):
                rhs = (self.c1 * q_ro) + \
                      (self.c2 * ((self.N @ (self.N @ q_t)) + q_ro)) + \
                      (self.c3 * q_t)
                q_t = bicgstab(lhs, rhs, x0=q_t)[0]
                # rhs_pc = preconditioner.solve(rhs)
                # q_t = bicgstab(lhs, rhs_pc, x0=q_t)[0]
                interval_flows[routing_substep_iteration, :] = q_t

            interval_flows[interval_flows < 0] = 0
            interval_flows = self.flow_temporal_aggregation_function(np.array(interval_flows), axis=0)
            outflow_array[inflow_time_step, :] = interval_flows

        logging.info('Writing Outflow Array to File')
        self.write_outflows(outflow_array)

        # write the final outflows to disc
        if self.configs.get('qfinal_file', False):
            logging.info('Writing Qfinal parquet')
            pd.DataFrame(interval_flows, columns=['Q', ]).astype(float).to_parquet(self.configs.get('qfinal_file'))

        logging.info('Routing computation Loops completed')
        return q_t

    def scaffold_outflow_file(self, dates) -> None:
        xr.Dataset(
            data_vars={
                'Qout': (['time', 'rivid'], np.zeros((self.num_outflow_steps, self.N.shape[0]))),
            },
            coords={
                'time': dates,
                'rivid': self.get_riverids(),
            },
            attrs={
                'long_name': 'Discharge at the outlet of each river reach',
                'units': 'm3 s-1',
                'standard_name': 'discharge',
                'aggregation_method': 'mean',
            },
        ).to_netcdf(
            self.configs['outflow_file'],
            mode='w',
        )
        with xr.open_dataset(self.configs['outflow_file'], mode='a') as ds:
            ds['time'].attrs['units'] = \
                f'days since {dates[0].astype("datetime64[s]").tolist().strftime("%Y-%m-%d %H:%M:%S")}'

        return

    def write_dates_to_outflow(self, dates: np.array) -> None:
        with xr.open_dataset(self.configs['outflow_file'], mode='a') as ds:
            ds['time'].attrs['units'] = \
                f'days since {dates[0].astype("datetime64[s]").tolist().strftime("%Y-%m-%d %H:%M:%S")}'
        return

    def write_outflows(self, outflow_array) -> None:
        with xr.open_dataset(self.configs['outflow_file'], mode='a') as ds:
            ds['Qout'][:] = outflow_array
        return

    def write_incremental_outflow(self, q_t, inflow_time_step) -> None:
        with xr.open_dataset(self.configs['outflow_file'], mode='a') as ds:
            ds['Qout'][inflow_time_step, :] = q_t
        return
