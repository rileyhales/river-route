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
import xarray as xr
import yaml
from petsc4py import PETSc


class RouteMuskingum:
    sim_config_file: str

    conf: dict

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
        return

    def read_and_validate_configs(self) -> None:
        """
        Validate simulation conf
        """
        if self.sim_config_file.endswith('.json'):
            with open(self.sim_config_file, 'r') as f:
                self.conf = json.load(f)
        elif self.sim_config_file.endswith('.yml') or self.sim_config_file.endswith('.yaml'):
            with open(self.sim_config_file, 'r') as f:
                self.conf = yaml.load(f, Loader=yaml.FullLoader)
        else:
            raise RuntimeError('Unrecognized simulation config file type')

        # check that required file paths are given and exist
        required_file_paths = ['routing_params_file',
                               'connectivity_file',
                               'inflow_file',
                               'outflow_file', ]
        paths_should_exist = ['routing_params_file',
                              'connectivity_file',
                              'inflow_file', ]

        # check for required time options
        required_time_opts = ['dt_total',
                              'dt_inflows',
                              'dt_outflows',
                              'dt_routing', ]

        # check for optional file paths
        optional_file_paths = ['qinit_file',
                               'qfinal_file',
                               'log_file', ]

        for arg in required_file_paths + required_time_opts:
            if arg not in self.conf:
                raise ValueError(f'{arg} not found in config file')
        for arg in paths_should_exist:
            if not os.path.exists(self.conf[arg]):
                raise FileNotFoundError(f'{arg} not found at given path')

        # start a logger
        log_basic_configs = {
            'stream': sys.stdout,
            'level': logging.INFO,
            'format': '%(asctime)s %(levelname)s %(message)s',
        }
        if self.conf.get('log_file', ''):
            log_basic_configs['filename'] = self.conf['log_file']
            log_basic_configs['filemode'] = 'w'
            log_basic_configs.pop('stream')
        logging.basicConfig(**log_basic_configs)

        # check that time options have the correct sizes
        logging.info('Validating time paramters')
        assert self.conf['dt_total'] >= self.conf['dt_inflows'], 'dt_total !>= dt_inflows'
        assert self.conf['dt_inflows'] >= self.conf['dt_outflows'], 'dt_inflows !>= dt_outflows'
        assert self.conf['dt_outflows'] >= self.conf['dt_routing'], 'dt_outflows !>= dt_routing'

        # check that time options are evenly divisible
        assert self.conf['dt_total'] % self.conf['dt_inflows'] == 0, \
            'dt_total should be a whole number multiple of dt_inflows'
        assert self.conf['dt_total'] % self.conf['dt_outflows'] == 0, \
            'dt_total should be a whole number multiple of dt_outflows'
        assert self.conf['dt_inflows'] % self.conf['dt_outflows'] == 0, \
            'dt_inflows should be a whole number multiple of dt_outflows'
        assert self.conf['dt_outflows'] % self.conf['dt_routing'] == 0, \
            'dt_outflows should be a whole number multiple of dt_routing'

        # set derived datetime parameters for computation cycles later
        self.num_outflow_steps = int(self.conf.get('dt_total') / self.conf.get('dt_outflows'))
        self.num_routing_substeps_per_outflow = int(self.conf.get('dt_outflows') / self.conf.get('dt_routing'))

        return

    def get_connectivity(self) -> pd.DataFrame:
        """
        Reads connectivity matrix from parquet given in config file
        """
        return pd.read_parquet(self.conf['connectivity_file'])

    def get_riverids(self) -> np.array:
        """
        Reads riverids vector from parquet given in config file
        """
        return pd.read_parquet(self.conf['routing_params_file'], columns=['id', ]).values.flatten()

    def get_k(self) -> np.array:
        """
        Reads K vector from parquet given in config file
        """
        return pd.read_parquet(self.conf['routing_params_file'], columns=['k', ]).values.flatten()

    def get_x(self) -> np.array:
        """
        Reads X vector from parquet given in config file
        """
        return pd.read_parquet(self.conf['routing_params_file'], columns=['x', ]).values.flatten()

    def get_qinit(self) -> np.array:
        qinit = self.conf.get('qinit_file', None)
        if qinit is None or qinit == '':
            return np.zeros(self.N.shape[0])
        return pd.read_parquet(self.conf['qinit_file']).values.flatten()

    def make_adjacency_array(self) -> None:
        """
        Calculate the connections array from the connectivity file
        """
        logging.info('Calculating Network Adjacency Array (N)')
        df = self.get_connectivity()
        G = nx.DiGraph()
        G.add_edges_from(df[df.columns[:2]].values)
        sorted_order = list(nx.topological_sort(G))
        if -1 in sorted_order:
            sorted_order.remove(-1)
        self.N = scipy.sparse.csc_matrix(nx.to_numpy_array(G, nodelist=sorted_order).T)
        return

    def calculate_muskingum_coefficients(self) -> None:
        """
        Calculate the 3 Muskingum Cunge routing coefficients for each segment using given k and x
        """
        logging.info('Calculating Muskingum coefficients')

        k = self.get_k()
        x = self.get_x()

        dt_route_half = self.conf['dt_routing'] / 2
        kx = k * x
        denom = k - kx + dt_route_half

        self.c1 = (dt_route_half - kx) / denom
        self.c2 = (dt_route_half + kx) / denom
        self.c3 = (k - kx - dt_route_half) / denom

        # sum of muskingum coefficiencts should be 1 for all segments
        a = self.c1 + self.c2 + self.c3
        assert np.allclose(a, 1), 'Muskingum coefficients are not approximately equal to 1'

        self.c1 = scipy.sparse.csc_matrix(np.diag(self.c1))
        self.c2 = scipy.sparse.csc_matrix(np.diag(self.c2))
        self.c3 = scipy.sparse.csc_matrix(np.diag(self.c3))
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

        if self.conf['aggregation_method'] == 'mean':
            self.flow_temporal_aggregation_function = np.mean
        elif self.conf['aggregation_method'] == 'max':
            self.flow_temporal_aggregation_function = np.max
        elif self.conf['aggregation_method'] == 'min':
            self.flow_temporal_aggregation_function = np.min
        elif self.conf['aggregation_method'] == 'instantaneous':
            self.flow_temporal_aggregation_function = lambda x: x[-1]
        else:
            raise RuntimeError('Unrecognized aggregation method bypassed config file validation')

    def route(self) -> None:
        self.make_adjacency_array()
        self.calculate_muskingum_coefficients()
        self.define_flow_temporal_aggregation_function()

        logging.info('Reading Inflow Data')
        with xr.open_dataset(self.conf['inflow_file']) as inflow_ds:
            # read dates from the netcdf
            dates = inflow_ds['time'].values
            # read inflows from the netcdf
            inflows = inflow_ds['m3_riv'].values
            # convert to m3/s in each routing time step
            inflows = inflows / self.num_routing_substeps_per_outflow / self.conf['dt_routing']

        logging.info('Scaffolding Outflow File')
        self.scaffold_outflow_file(dates)

        logging.info('Preparing Arrays')
        lhs = scipy.sparse.csc_matrix(np.eye(self.N.shape[0]) - (self.N @ self.c1))
        outflow_array = np.zeros((self.num_outflow_steps, self.N.shape[0]))
        interval_flows = np.zeros((self.num_routing_substeps_per_outflow, self.N.shape[0]))
        q_t = self.get_qinit()

        logging.info('Init PETSc Objects and Options')
        A = PETSc.Mat().createAIJ(size=lhs.shape, csr=(lhs.indptr, lhs.indices, lhs.data))
        x = PETSc.Vec().createSeq(size=lhs.shape[0])
        b = PETSc.Vec().createSeq(size=lhs.shape[0])
        ksp = PETSc.KSP().create()
        ksp.setType('bicg')
        ksp.setTolerances(rtol=1e-2)
        ksp.setOperators(A)

        logging.info('Performing routing computation iterations')
        time_start = datetime.datetime.now()
        for inflow_time_step, inflow_end_date in enumerate(dates):
            q_ro = inflows[inflow_time_step, :]
            interval_flows = np.zeros((self.num_routing_substeps_per_outflow, self.N.shape[0]))
            c1_matmul_q_ro = self.c1 @ q_ro

            for routing_substep_iteration in range(self.num_routing_substeps_per_outflow):
                # solve the right hand side of the equation
                rhs = c1_matmul_q_ro + \
                      (self.c2 @ (self.N @ q_t + q_ro)) + \
                      (self.c3 @ q_t)

                # set current iteration values in PETSc objects and solve
                b.setArray(rhs)
                x.setArray(q_t)
                ksp.solve(b, x)
                q_t = x.getArray()

                # remove negatives before other iterations
                q_t[q_t < 0] = 0
                interval_flows[routing_substep_iteration, :] = q_t

            interval_flows = self.flow_temporal_aggregation_function(np.array(interval_flows), axis=0)
            outflow_array[inflow_time_step, :] = interval_flows
        time_end = datetime.datetime.now()
        logging.info(f'Routing completed in {(time_end - time_start).total_seconds()} seconds')

        logging.info('Cleaning up PETSc objects')
        A.destroy()
        x.destroy()
        b.destroy()
        ksp.destroy()

        logging.info('Writing Outflow Array to File')
        outflow_array = np.round(outflow_array, decimals=2)
        self.write_outflows(outflow_array)

        # write the final outflows to disc
        if self.conf.get('qfinal_file', False):
            logging.info('Writing Qfinal parquet')
            pd.DataFrame(interval_flows, columns=['Q', ]).astype(float).to_parquet(self.conf.get('qfinal_file'))

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
            self.conf['outflow_file'],
            mode='w',
        )
        return

    def write_outflows(self, outflow_array) -> None:
        with nc.Dataset(self.conf['outflow_file'], mode='a') as ds:
            ds['Qout'][:] = outflow_array
            ds.sync()
        return

    def plot(self, rivid: int) -> None:
        rivid_index = np.where(self.get_riverids() == rivid)[0][0]
        with xr.open_dataset(self.conf['outflow_file']) as ds:
            ds['Qout'][:, rivid_index].to_dataframe()['Qout'].plot()
            plt.show()
        return

    def mass_balance(self) -> None:
        outlet_ids = self.get_connectivity().values
        outlet_ids = outlet_ids[outlet_ids[:, 1] == -1, 0]

        with xr.open_dataset(self.conf['outflow_file']) as ds:
            out_df = ds.sel(rivid=outlet_ids).to_dataframe()[['Qout', ]].groupby('time').sum().cumsum()
            out_df = out_df * self.conf['dt_routing'] * self.num_routing_substeps_per_outflow
        with xr.open_dataset(self.conf['inflow_file']) as ds:
            in_df = ds.sel(rivid=outlet_ids).to_dataframe()[['m3_riv', ]].groupby('time').sum().cumsum()

        df = out_df.merge(in_df, left_index=True, right_index=True)
        logging.info(f'\n{df.sum()}')
        df.plot()
        plt.show()
