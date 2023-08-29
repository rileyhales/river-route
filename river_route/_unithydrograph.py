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


class UnitHydrograph:
    # Given configs
    conf: dict

    # Routing Matrices
    A: np.array or scipy.sparse.csc_matrix

    def __init__(self,
                 config_file: str = None,
                 **kwargs, ) -> None:
        """
        Read config files to initialize routing class
        """
        self.read_configs(config_file, **kwargs)
        return

    def read_connectivity(self) -> pd.DataFrame:
        """
        Reads connectivity matrix from parquet given in config file
        """
        return pd.read_parquet(self.conf['connectivity_file'])

    def get_directed_graph(self) -> nx.DiGraph:
        """
        Returns a directed graph of the river network
        """
        df = self.read_connectivity()
        G = nx.DiGraph()
        G.add_edges_from(df.values)
        return G

    def set_adjacency_matrix(self) -> None:
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
        G = self.get_directed_graph()
        sorted_order = list(nx.topological_sort(G))
        if -1 in sorted_order:
            sorted_order.remove(-1)
        self.A = scipy.sparse.csc_matrix(nx.to_scipy_sparse_array(G, nodelist=sorted_order).T)
        if self.conf.get('adj_file', ''):
            scipy.sparse.save_npz(self.conf['adj_file'], self.A)
        return

    def read_configs(self, config_file, **kwargs) -> None:
        """
        Validate simulation conf
        """
        # read the config file
        if config_file.endswith('.json'):
            with open(config_file, 'r') as f:
                self.conf = json.load(f)
        elif config_file.endswith('.yml') or config_file.endswith('.yaml'):
            with open(config_file, 'r') as f:
                self.conf = yaml.load(f, Loader=yaml.FullLoader)
        elif config_file is None or config_file == '':
            self.conf = {}
        else:
            raise RuntimeError('Unrecognized simulation config file type')

        # overwrite configs with kwargs
        self.conf.update(kwargs)

        if isinstance(self.conf['runoff_file'], str):
            self.conf['runoff_file'] = [self.conf['runoff_file'], ]
        if isinstance(self.conf['outflow_file'], str):
            self.conf['outflow_file'] = [self.conf['outflow_file'], ]

        for arg in [k for k in self.conf.keys() if 'file' in k]:
            if isinstance(self.conf[arg], list):
                self.conf[arg] = [
                    os.path.abspath(os.path.join(os.path.dirname(config_file), f))
                    if f.startswith('.') else f for f in self.conf[arg]
                ]
            elif self.conf[arg].startswith('.'):
                self.conf[arg] = os.path.abspath(os.path.join(os.path.dirname(config_file), self.conf[arg]))

        # start a logger
        log_basic_configs = {
            'stream': sys.stdout,
            'level': logging.INFO,
            'format': '%(asctime)s %(message)s',
        }
        if self.conf.get('log_file', ''):
            log_basic_configs['filename'] = self.conf['log_file']
            log_basic_configs['filemode'] = 'w'
            log_basic_configs.pop('stream')
        logging.basicConfig(**log_basic_configs)
        return

    def read_riverids(self) -> np.array:
        """
        Reads river ids vector from parquet given in config file
        """
        return pd.read_parquet(self.conf['routing_params_file'], columns=['rivid', ]).values.flatten()

    def read_lag_times(self) -> np.array:
        """
        Reads lag times vector from parquet given in config file
        """
        return pd.read_parquet(self.conf['routing_params_file'], columns=['lag_time', ]).values.flatten()

    def route(self) -> None:
        """
        Performs time-iterative runoff routing through the river network
        """
        self.set_adjacency_matrix()
        logging.info('Beginning UH Transform')
        t1 = datetime.datetime.now()

        # todo handle the time step of the runoff vs the unit hydrograph
        dt_uh = 900

        for runoff_file, outflow_file in zip(self.conf['runoff_file'], self.conf['outflow_file']):
            logging.info(f'Reading Inflow Data: {runoff_file}')
            with xr.open_dataset(runoff_file) as runoff_ds:
                dates = runoff_ds['time'].values.astype('datetime64[s]')
                dt_runoff = (dates[1] - dates[0]).astype(int)
                runoffs = runoff_ds['m3_riv'].values
                runoffs[runoffs < 0] = np.nan
                runoffs = np.nan_to_num(runoffs, nan=0.0)

            num_uh_steps_per_runoff = int(dt_runoff / dt_uh)

            # todo read UH from file
            # todo check the uh time deltas within
            uhseries = [0, .075, .15, .225, .20, .125, .075, .05, 0.25, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
            uh = np.array(uhseries).reshape(-1, 1) / sum(uhseries) / dt_uh
            uh = np.broadcast_to(uh, (uh.shape[0], runoffs.shape[1]))

            Qro = np.concatenate([
                np
                .convolve(np.repeat(runoffs[:, i], num_uh_steps_per_runoff, axis=0), uh[:, i], mode='same')
                .reshape(-1, 1)
                for i in tqdm.tqdm(range(runoffs.shape[1]), desc='River Segment Convolutions')
            ], axis=1)

            Qro = (
                Qro
                .reshape((
                    runoffs.shape[0],
                    num_uh_steps_per_runoff,
                    runoffs.shape[1],
                ))
                .mean(axis=1)
                .round(2)
                [:dates.shape[0], :]  # clip to the date of the last runoff
            )

            lag_index_steps = (self.read_lag_times() / dt_uh).astype(int)

            # for each row of Qro, add the upstreams rows from Qro but offset by the lag time for the current row
            for column in tqdm.tqdm(range(Qro.shape[1]), desc='Lagging Upstream Flows'):
                # select the rows from Qro which have a value of 1 in the adjacency matrix
                Qro[:, column] += (
                    np.roll(
                        np.pad(
                            Qro
                            [:, (self.A[column, :] == 1).toarray().squeeze()]  # select columns of upstream segments
                            .sum(axis=1),  # sum the upstream flows
                            lag_index_steps[column],
                            mode='constant',
                            constant_values=0  # pad with zeros
                        ),
                        lag_index_steps[column]
                    )
                    [:Qro.shape[0]]  # clip to the length of Qro
                )

            self.write_outflows(outflow_file, dates, Qro)

        t2 = datetime.datetime.now()
        logging.info('All UH Convolutions Completed')
        logging.info(f'Total Transformation time: {(t2 - t1).total_seconds()}')
        return

    def write_outflows(self, outflow_file: str, dates: np.array, outflow_array: np.array) -> None:
        pydates = list(map(datetime.datetime.utcfromtimestamp, dates.astype(int)))
        with nc.Dataset(outflow_file, mode='w') as ds:
            ds.createDimension('time', size=dates.shape[0])
            ds.createDimension('rivid', size=outflow_array.shape[1])

            ds.createVariable('time', 'f8', ('time',))
            ds['time'].units = f'seconds since {pydates[0].strftime("%Y-%m-%d %H:%M:%S")}'
            ds['time'][:] = nc.date2num(pydates, units=ds['time'].units)

            ds.createVariable('rivid', 'i4', ('rivid',))
            ds['rivid'][:] = self.read_riverids()

            ds.createVariable('Qout', 'f4', ('time', 'rivid'))
            ds['Qout'].units = 'm3 s-1'
            ds['Qout'].long_name = 'Discharge at the outlet of each river reach'
            ds['Qout'].standard_name = 'discharge'
            ds['Qout'].aggregation_method = 'mean'
            ds['Qout'][:] = outflow_array
        return

    def plot(self, rivid: int) -> None:
        with xr.open_mfdataset(self.conf['outflow_file']) as ds:
            ds['Qout'].sel(rivid=rivid).to_dataframe()['Qout'].plot()
            plt.show()
        return
