import datetime
import json
import logging
import os
import random
import sys
from typing import Any, Callable, Generator

import netCDF4 as nc
import networkx as nx
import numpy as np
import pandas as pd
import tqdm
import xarray as xr
import yaml
from numpy.typing import NDArray
from scipy.sparse import csc_matrix
from scipy.sparse import diags
from scipy.sparse import eye
from scipy.sparse.linalg import factorized

from .__metadata__ import __version__
from .runoff import calc_catchment_volumes
from .tools import connectivity_to_digraph

__all__ = ['Muskingum', ]

# type hints for checkers
FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
DatetimeArray = NDArray[np.datetime64]
ConfigDict = dict[str, Any]
WriteDischargesFn = Callable[[DatetimeArray, FloatArray, str, str], None]
FactorizedSolveFn = Callable[[FloatArray], FloatArray]


class Muskingum:
    # Given configs
    conf: ConfigDict

    # Network dependent matrices and vectors
    A: csc_matrix  # n x n - adjacency matrix => from connectivity file
    k: FloatArray  # n x 1 - K values for each segment => from routing params file
    x: FloatArray  # n x 1 - X values for each segment => from routing params file
    river_ids: IntArray  # n x 1 - river ID for each segment => from routing params file
    # Network and routing timestep dependent matrices and vectors
    c1: FloatArray  # n x 1 - C1 values for each segment => f(k, x, dt_routing)
    c2: FloatArray  # n x 1 - C2 values for each segment => f(k, x, dt_routing)
    c3: FloatArray  # n x 1 - C3 values for each segment => f(k, x, dt_routing)
    lhs: csc_matrix  # n x n - left hand side of matrix form of routing equation => f(n, dt_routing)
    lhs_factorized: FactorizedSolveFn  # n x n - factorized left hand side matrix for direct solutions
    _network_time_signature: tuple[Any, ...] | None  # used to track if sequential runoffs use the same routing matrices

    # State variables
    initial_state: tuple[FloatArray, FloatArray]  # Q init, R init
    _ensemble_member_states: list[FloatArray]  # used for ensemble routing, stores member states for final state aggregation

    # Time options
    dt_total: float
    dt_runoff: float
    dt_discharge: float
    dt_routing: float
    num_routing_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_discharge: int

    # Methods
    _write_discharges: WriteDischargesFn

    # Logger
    logger: logging.Logger

    def __init__(self, config_file: str | None = None, **kwargs: Any) -> None:
        self.logger = logging.getLogger(f'river_route.{''.join([str(random.randint(0, 9)) for _ in range(8)])}')
        self._set_configs(config_file, **kwargs)
        return

    def _set_configs(self, config_file: str | None, **kwargs: Any) -> None:
        """
        Validate simulation configs given by JSON, YAML, or kwargs

        You can also override the values in the config file or entirely ignore the config file by specifying config
        options using keyword arguments. See README or Docs for list of all recognized configuration options. If you do
        not want to use a config file, pass None to the first positional argument and then give kwargs.

        Args:
            config_file (str): path to a JSON or YAML file of configs. See README for all recognized options.

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
        # these 2 files must always exist
        assert os.path.exists(self.conf['routing_params_file']), FileNotFoundError('Routing params file not found')
        assert os.path.exists(self.conf['connectivity_file']), FileNotFoundError('Connectivity file not found')

        # check for valid options
        assert self.conf['input_type'] in ['sequential', 'ensemble'], 'Input type not recognized'
        assert self.conf['runoff_type'] in ['incremental', 'cumulative'], 'Runoff type not recognized'

        # format conversion for the inputs
        if isinstance(self.conf.get('catchment_volumes_files', []), str):
            self.conf['catchment_volumes_files'] = [self.conf['catchment_volumes_files'], ]
        if isinstance(self.conf.get('runoff_depths_files', []), str):
            self.conf['runoff_depths_files'] = [self.conf['runoff_depths_files'], ]
        if isinstance(self.conf.get('discharge_files', []), str):
            self.conf['discharge_files'] = [self.conf['discharge_files'], ]

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
        assert self.conf.get('catchment_volumes_files', False) or self.conf.get('runoff_depths_files', False), \
            'Either catchment volumes or runoff depths files must be provided'

        # if catchment volumes files are provided, they must exist
        if self.conf.get('catchment_volumes_files', False):
            for file in self.conf['catchment_volumes_files']:
                assert os.path.exists(file), FileNotFoundError(f'Runoff volumes file not found at: {file}')

        # if runoff depths files are provided, they must exist
        if self.conf.get('runoff_depths_files', False):
            for file in self.conf['runoff_depths_files']:
                assert os.path.exists(file), FileNotFoundError(f'Runoff depths file not found at: {file}')
            assert self.conf.get('weight_table_file', ''), 'weight_table_file should be provided'

        # the directory for each discharge file must exist
        assert self.conf.get('discharge_files', False), 'discharge_files must be provided'
        for discharge_file in self.conf['discharge_files']:
            assert os.path.exists(os.path.dirname(discharge_file)), NotADirectoryError('Output directory not found')

        # convert all relative paths to absolute paths
        for arg in [k for k in self.conf.keys() if 'file' in k]:
            if isinstance(self.conf[arg], list):
                self.conf[arg] = [os.path.abspath(path) for path in self.conf[arg]]
            elif isinstance(self.conf[arg], str):
                self.conf[arg] = os.path.abspath(self.conf[arg])
        return

    def _log_configs(self) -> None:
        self.logger.debug(f'river-route version: {__version__}')
        self.logger.debug(f'Number of Rivers: {self.river_ids.shape[0]}')
        self.logger.debug('Configs:')
        for k, v in self.conf.items():
            self.logger.debug(f'\t{k}: {v}')
        return

    def _read_initial_state(self) -> None:
        if hasattr(self, 'initial_state'):
            return

        state_file = self.conf.get('initial_state_file', '')
        if state_file == '':
            self.logger.debug('Setting initial state to zeros')
            self.initial_state = (np.zeros(self.A.shape[0]), np.zeros(self.A.shape[0]))
            return
        assert os.path.exists(state_file), FileNotFoundError(f'Initial state file not found at: {state_file}')
        self.logger.debug('Reading Initial State from Parquet')
        initial_state = pd.read_parquet(state_file).values
        initial_state = (initial_state[:, 0].flatten(), initial_state[:, 1].flatten())
        self.initial_state = initial_state
        return

    def _write_final_state(self) -> None:
        final_state_file = self.conf.get('final_state_file', '')
        if final_state_file == '':
            return
        if self.conf['input_type'] == 'sequential':
            array = np.array(self.initial_state).T
        elif self.conf['input_type'] == 'ensemble':
            array = np.array(self._ensemble_member_states).mean(axis=0)
        else:
            raise ValueError('Unexpected input type when writing final states')
        self.logger.debug('Writing Final State to Parquet')
        pd.DataFrame(array, columns=['Q', 'R']).to_parquet(self.conf['final_state_file'])
        return

    def _set_network_dependent_vectors(self) -> None:
        self.logger.debug('Calculating network dependent vectors')
        df = pd.read_parquet(self.conf['routing_params_file'])
        self.river_ids = df[self.conf['var_river_id']].to_numpy(dtype=np.int64, copy=False)
        self.k = df['k'].to_numpy(dtype=np.float64, copy=False)
        self.x = df['x'].to_numpy(dtype=np.float64, copy=False)
        graph = connectivity_to_digraph(self.conf['connectivity_file'])
        self.A = csc_matrix(nx.convert_matrix.to_scipy_sparse_array(graph, nodelist=self.river_ids.tolist()).T)
        return

    def _set_network_and_time_dependent_vectors(self, dates: DatetimeArray) -> None:
        self.logger.debug('Setting and validating time parameters')
        self.dt_runoff = self.conf.get('dt_runoff', (dates[1] - dates[0]).astype('timedelta64[s]').astype(int))
        self.dt_discharge = self.conf.get('dt_discharge', self.dt_runoff)
        self.dt_total = self.conf.get('dt_total', self.dt_runoff * dates.shape[0])
        if not self.conf.get('dt_routing', 0):
            self.logger.warning('dt_routing was not provided or is Null/False, defaulting to dt_runoff')
        self.dt_routing = self.conf.get('dt_routing', self.dt_runoff)

        try:
            # check that time options have the correct sizes
            assert self.dt_total >= self.dt_runoff, 'dt_total !>= dt_runoff'
            assert self.dt_total >= self.dt_discharge, 'dt_total !>= dt_discharge'
            assert self.dt_discharge >= self.dt_runoff, 'dt_discharge !>= dt_runoff'
            assert self.dt_runoff >= self.dt_routing, 'dt_runoff !>= dt_routing'
            # check that time options are evenly divisible
            assert self.dt_total % self.dt_runoff == 0, 'dt_total must be an integer multiple of dt_runoff'
            assert self.dt_total % self.dt_discharge == 0, 'dt_total must be an integer multiple of dt_discharge'
            assert self.dt_discharge % self.dt_runoff == 0, \
                'dt_discharge must be an integer multiple of dt_runoff'
            assert self.dt_runoff % self.dt_routing == 0, 'dt_runoff must be an integer multiple of dt_routing'
        except AssertionError as e:
            self.logger.error(e)
            raise AssertionError('Time options are not valid')

        # set derived datetime parameters for computation cycles later
        self.num_runoff_steps_per_discharge = int(self.dt_discharge / self.dt_runoff)
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)
        self.num_routing_steps = int(self.dt_total / self.dt_routing)

        signature = (self.dt_total, self.dt_runoff, self.dt_discharge, self.dt_routing,)
        if self._network_time_signature == signature:
            return

        self.logger.debug('Calculating network and time dependent vectors')
        dt_div_k = self.dt_routing / self.k
        denominator = dt_div_k + (2 * (1 - self.x))
        _2x = 2 * self.x
        self.c1 = (dt_div_k + _2x) / denominator
        self.c2 = (dt_div_k - _2x) / denominator
        self.c3 = ((2 * (1 - self.x)) - dt_div_k) / denominator

        # sum of coefficients should be 1 (or arbitrarily close) for all segments
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
        self._network_time_signature = signature
        return

    def _volumes_output_generator(self) -> Generator[tuple[DatetimeArray, FloatArray, str, str], None, None]:
        if self.conf.get('catchment_volumes_files', False):
            for volume_file, discharge_file in zip(self.conf['catchment_volumes_files'], self.conf['discharge_files']):
                self.logger.info(f'Reading catchment volumes file: {volume_file}')
                with xr.open_dataset(volume_file) as runoff_ds:
                    self.logger.debug('Reading time array')
                    dates = runoff_ds['time'].values.astype('datetime64[s]')
                    self.logger.debug('Reading volume array')
                    volumes_array = runoff_ds[self.conf['var_catchment_volume']].values
                yield dates, volumes_array, volume_file, discharge_file
        elif self.conf.get('runoff_depths_files', False):
            for runoff_file, discharge_file in zip(self.conf['runoff_depths_files'], self.conf['discharge_files']):
                self.logger.info(f'Calculating catchment volumes from runoff depths: {runoff_file}')
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
                yield volumes_df.index.values, volumes_df.values, runoff_file, discharge_file
        else:
            raise ValueError('No runoff data found in configs. Provide catchment volumes or runoff depths.')

    def route(self) -> 'Muskingum':
        """
        Performs time-iterative runoff routing through the river network

        Returns:
            river_route.Muskingum
        """
        self.logger.info(f'Beginning routing')
        t1 = datetime.datetime.now()

        self._validate_configs()
        self._set_network_dependent_vectors()
        self._log_configs()

        for dates, volumes_array, runoff_file, discharge_file in self._volumes_output_generator():
            self.logger.info('-' * 80)
            self._set_network_and_time_dependent_vectors(dates)
            volumes_array = volumes_array.astype(np.float64, copy=False)
            np.divide(volumes_array, self.dt_runoff, out=volumes_array)  # convert volume -> volume/time
            discharge_array = self._router(dates, volumes_array)

            if self.dt_discharge > self.dt_runoff:
                self.logger.info('Resampling dates and discharges to specified timestep')
                discharge_array = (
                    discharge_array
                    .reshape((
                        int(self.dt_total / self.dt_discharge),
                        int(self.dt_discharge / self.dt_runoff),
                        self.A.shape[0],
                    ))
                    .mean(axis=1)
                )

            self.logger.info('Writing Discharge Array to File')
            np.round(discharge_array, decimals=2, out=discharge_array)
            discharge_array = discharge_array.astype(np.float32, copy=False)
            self._write_discharges(dates, discharge_array, discharge_file, runoff_file)

        # write the final state to disc
        self._write_final_state()

        t2 = datetime.datetime.now()
        self.logger.info('All runoff files routed')
        self.logger.info(f'Total job time: {(t2 - t1).total_seconds()}')
        return self

    def _router(self, dates: DatetimeArray, volumes: FloatArray) -> FloatArray:
        self.logger.debug('Getting initial state arrays')
        self._read_initial_state()
        q_init, r_init = self.initial_state  # n x 1 arrays of initial discharges and runoff in each segment
        self._ensemble_member_states = []

        # declare arrays needed to do routing computations once to avoid repeated and duplicate allocations in the loop
        discharge_array = np.zeros((volumes.shape[0], self.A.shape[0]))
        q_t = np.empty_like(q_init, dtype=np.float64)
        r_prev = np.empty_like(r_init, dtype=np.float64)
        rhs = np.zeros(self.A.shape[0], dtype=np.float64)
        buffer = np.zeros(self.A.shape[0], dtype=np.float64)
        c2_r_t = np.zeros(self.A.shape[0], dtype=np.float64)
        interval_sum = np.zeros(self.A.shape[0], dtype=np.float64)

        q_t[:] = q_init
        r_prev[:] = r_init

        t1 = datetime.datetime.now()
        if self.conf['progress_bar']:
            dates = tqdm.tqdm(dates, desc='Runoff Volumes Routed')
        else:
            self.logger.info('Performing routing computation iterations')
        for runoff_time_step, runoff_end_date in enumerate(dates):
            r_t = volumes[runoff_time_step, :]
            np.multiply(self.c2, r_t, out=c2_r_t)
            interval_sum.fill(0.0) # add then divide to avoid accumulating large array and averaging
            for routing_sub_iteration_num in range(self.num_routing_steps_per_runoff):
                # rhs = (self.c1 * ((self.A @ q_t) + r_prev)) + (self.c2 * r_t) + (self.c3 * q_t)
                buffer[:] = self.A @ q_t
                np.add(buffer, r_prev, out=buffer)
                np.multiply(self.c1, buffer, out=rhs)
                np.multiply(self.c3, q_t, out=buffer)
                np.add(rhs, buffer, out=rhs)
                np.add(rhs, c2_r_t, out=rhs)
                # solve Ax=b, LHS q_t = rhs
                q_t[:] = self.lhs_factorized(rhs)
                # add to interval accumulator for averaging at the end of the runoff time step
                interval_sum += q_t
            discharge_array[runoff_time_step, :] = interval_sum / self.num_routing_steps_per_runoff
            r_prev[:] = r_t

        # Enforce positive flows in case of negative numerical results
        discharge_array[discharge_array < 0] = 0

        # if simulation type is ensemble, then do not overwrite the initial state
        if self.conf['input_type'] == 'sequential':
            self.logger.debug('Updating Initial State for Next Sequential Computation')
            self.initial_state = (q_t, r_prev)
        if self.conf['input_type'] == 'ensemble':
            self.logger.debug('Recording Member State for Final State Aggregation')
            self._ensemble_member_states.append(np.array([q_t, r_prev]).T)

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return discharge_array

    def _write_discharges(self, dates: DatetimeArray, discharge_array: FloatArray, discharge_file: str, runoff_file: str, ) -> None:
        """
        Writes routed discharge from a routing simulation to a netcdf file.
        You should overwrite this method with a custom handler using set_write_discharges.

        Args:
            dates: datetime array corresponding to the discharge rows
            discharge_array: routed discharge values with shape (time, river)
            discharge_file: the file path to write the discharge data to
            runoff_file: the file path to the runoff file used to generate the discharge values

        Returns:
            None
        """
        dates = dates[::self.num_runoff_steps_per_discharge].astype('datetime64[s]')

        with nc.Dataset(discharge_file, mode='w', format='NETCDF4') as ds:
            ds.createDimension('time', size=discharge_array.shape[0])
            ds.createDimension(self.conf['var_river_id'], size=discharge_array.shape[1])
            ds.runoff_file = runoff_file

            time_var = ds.createVariable('time', 'f8', ('time',))
            t0 = pd.Timestamp(dates[0]).strftime("%Y-%m-%d %H:%M:%S")
            time_var.units = f'seconds since {t0}'
            time_var[:] = dates - dates[0]

            id_var = ds.createVariable(self.conf['var_river_id'], 'i4', (self.conf['var_river_id']), )
            id_var[:] = self.river_ids

            flow_var = ds.createVariable(self.conf['var_discharge'], 'f4', ('time', self.conf['var_river_id']))
            flow_var[:] = discharge_array
            flow_var.long_name = 'Discharge at catchment outlet'
            flow_var.standard_name = 'discharge'
            flow_var.aggregation_method = 'mean'
            flow_var.units = 'm3 s-1'
        return

    def set_write_discharges(self, func: WriteDischargesFn) -> 'Muskingum':
        """
        Overwrites the default write_discharges method to a custom function and returns the class instance so that you
        can chain the method with the constructor.

        Args:
            func (callable): function that takes dates, discharge_array, discharge_file, runoff_file and returns None

        Returns:
            river_route.Muskingum
        """
        self._write_discharges = func
        return self

    def hydrograph(self, river_id: int) -> pd.DataFrame:
        """
        Get the hydrograph for a given river id as a pandas dataframe

        Args:
            river_id: the ID of a river reach in the output files

        Returns:
            pandas.DataFrame
        """
        with xr.open_mfdataset(self.conf['discharge_files']) as ds:
            df = (
                ds
                [self.conf['var_discharge']]
                .sel(**{self.conf['var_river_id']: river_id})
                .to_dataframe()
                [[self.conf['var_discharge'], ]]
            )
            df.columns = [river_id, ]
        return df

    def mass_balance(self, river_id: int, ancestors: list[int] | set[int] | None = None) -> pd.DataFrame:
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
        with xr.open_mfdataset(self.conf['catchment_volumes_files']) as ds:
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
        with xr.open_mfdataset(self.conf['discharge_files']) as ds:
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
            self.logger.warning(f'More discharge than runoff volume for river {river_id}')

        return df
