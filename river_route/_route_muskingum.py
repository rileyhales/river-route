import json
import logging

import numpy as np
import pandas as pd
import yaml
import numba as nb
from numba import jit, njit, float64

logger = logging.getLogger(__name__)


numba_class_specs = [
    ()
]


class RouteMuskingum:
    sim_config_file: str
    inflow_file: str
    q_init_file: str
    q_final_file: str

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
        self.read_validate_configs()

        logger.info('Calculating connections array')
        self.calculate_connections_array()
        logger.info('Calculating muskingum coefficients')
        self.calculate_muskingum_coefficients()

        self.define_flow_temporal_aggregation_function()

        return

    def read_validate_configs(self) -> None:
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
                               'qfinal_file', ]
        paths_should_exist = ['routing_params_file',
                              'connectivity_file',
                              'inflow_file', ]
        # for arg in required_file_paths + required_time_opts:
        #     if arg not in self.configs:
        #         raise ValueError(f'{arg} not found in config file')
        # for arg in paths_should_exist:
        #     if not os.path.exists(self.configs[arg]):
        #         raise FileNotFoundError(f'{arg} not found at given path')

        # todo check that number of rivers in the inputs matches the number in the inflows

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
        # todo calculate dates and iterate over these with enumerate later

        return

    def calculate_muskingum_coefficients(self) -> None:
        """
        Calculate the 3 Muskingum Cunge routing coefficients for each segment using given k and x
        """
        k = self.k()
        x = self.x()

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
        logger.info('Calculating connectivity array')

        # df = pd.read_parquet(self.configs['connectivity_file']).reset_index(names=['topo_order'])
        # todo replace with real connectivity file
        df = pd.DataFrame({
            'id': [1, 2, 3, 4, 5],
            'us1': [-1, -1, 1, -1, 3],
            'us2': [-1, -1, 2, -1, 4]
        })
        odf = df[['id', ]].copy().reset_index().rename(columns={'index': 'topo_order'})
        self.N = np.zeros((df.values.shape[0], df.values.shape[0]))
        for idx, row in df.iterrows():
            upstreams = (
                pd
                .Series(row.values[1:])
                .replace(-1, np.nan)
                .dropna()
                .apply(lambda x: odf.loc[odf['id'] == x, 'topo_order'].values[0])
                .values
            )
            if len(upstreams):
                self.N[idx, upstreams] = 1
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

    def k(self) -> np.array:
        """
        Reads K vector from parquet given in config file
        """
        # return pd.read_parquet(self.configs['routing_params_file'], columns=['k', ]).values
        return np.ones(self.N.shape[0]) * 100

    def x(self) -> np.array:
        """
        Reads X vector from parquet given in config file
        """
        # return pd.read_parquet(self.configs['routing_params_file'], columns=['x', ]).values
        return np.ones(self.N.shape[0]) * 0.25

    def qinit(self) -> np.array:
        qinit = self.configs.get('qinit_file', None)
        if qinit is None or qinit == '':
            return np.zeros(self.N.shape[0])
        return pd.read_parquet(self.configs['qinit_file']).values

    def route(self):
        q_t = self.qinit()

        left_hand_side = np.diag(1 - self.c1) @ self.N
        logger.info('Beginning computation loops')

        # todo this counts the number of loops to make but need a way to pair with dates
        for outflow_interval in range(int(self.configs['dt_total'] / self.configs['dt_outflows'])):
            q_ro = np.random.rand(self.N.shape[0], ) * 3
            interval_flows = []
            # todo pair with actual datetimes
            for routing_substep in range(self.num_routing_substeps_per_outflow):
                # todo get the correct q_ro based on the iteration and computation substep number
                right_hand_side = (self.c1 @ q_ro) + \
                                  (self.c2 * np.matmul(self.N, np.matmul(self.N, q_t) + q_ro)) + \
                                  (self.c3 * q_t)

                # q_tplusdt, residuals, rank, s = np.linalg.lstsq(left_hand_side, right_hand_side)
                q_tplusdt, _, _, _ = np.linalg.lstsq(left_hand_side, right_hand_side, rcond=None)
                q_tplusdt = np.clip(q_tplusdt, 0, None)
                q_t = q_tplusdt
                interval_flows.append(q_tplusdt)

            q = self.flow_temporal_aggregation_function(np.array(interval_flows))
            # todo write outflows to file

        # write the final outflows to disc
        if self.configs.get('qfinal_path', False):
            logger.info('Writing Qfinal parquet')
            pd.DataFrame(q, columns=['q', ]).astype(float).to_parquet(self.configs.get('qfinal_path'))

        logger.info('Routing computation Loops completed')
        return q_t

    def scaffold_outflow_zarr(self):
        # todo
        return

    def write_incremental_outflow(self):
        # todo
        return
