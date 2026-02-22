import datetime
import os
from pathlib import Path
from typing import Any, Generator, Self, List, Tuple

import numpy as np
import pandas as pd
import tqdm
import xarray as xr

from .AbstractRouter import AbstractRouter
from .types import DatetimeArray, FloatArray, PathInput
from ..transformers import AbstractBaseTransformer, Transformer
from ..runoff import depth_to_volume

__all__ = ['UnitMuskingum']

GeneratorSignature = Generator[Tuple[DatetimeArray, FloatArray, PathInput, PathInput], None, None]


class UnitMuskingum(AbstractRouter):
    _network_time_signature: Tuple[Any, ...] | None = None
    _ensemble_member_states: List[FloatArray]

    _runoff_transformer: AbstractBaseTransformer | None = None

    # Time options
    dt_runoff: float
    num_routing_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_discharge: int

    def __init__(self, config_file: PathInput | None = None, **kwargs: Any) -> None:
        super().__init__(config_file, **kwargs)

    def set_transformer(self, transformer: AbstractBaseTransformer) -> Self:
        """
        Inject a user instantiated transformer object including custom subclasses of AbstractBaseTransformer allowing
        users to implement their own transformation logic in a custom transformer. The transformer must be a
        AbstractBaseTransformer subclass with a kernel and state set.

        You must call this before route() to have any affect.
        """
        if not isinstance(transformer, AbstractBaseTransformer):
            raise TypeError(f'transformer must be a AbstractBaseTransformer, got {type(transformer).__name__}')
        self._runoff_transformer = transformer
        return self

    def _validate_router_configs(self) -> None:
        self.conf['lateral_depth_files'] = self.conf.get('lateral_depth_files', [])
        self.conf['runoff_depth_grids'] = self.conf.get('runoff_depth_grids', [])
        self.conf['discharge_files'] = self.conf.get('discharge_files', [])
        if isinstance(self.conf['lateral_depth_files'], (str, Path)):
            self.conf['lateral_depth_files'] = [self.conf['lateral_depth_files'], ]
        if isinstance(self.conf['runoff_depth_grids'], (str, Path)):
            self.conf['runoff_depth_grids'] = [self.conf['runoff_depth_grids'], ]
        if isinstance(self.conf['discharge_files'], (str, Path)):
            self.conf['discharge_files'] = [self.conf['discharge_files'], ]

        self.conf['runoff_depth_grids'] = [os.path.abspath(p) for p in self.conf['runoff_depth_grids']]

        n_files = len(self.conf['lateral_depth_files']) + len(self.conf['runoff_depth_grids'])
        if self.conf['lateral_depth_files'] and self.conf['runoff_depth_grids']:
            raise ValueError('Provide either lateral_depth_files or runoff_depth_grids, not both')
        if not len(self.conf['discharge_files']) == n_files:
            raise ValueError('Number of discharge files must match the number of input files')
        for file in self.conf['lateral_depth_files']:
            if not os.path.exists(file):
                raise FileNotFoundError(f'Lateral depth file not found at: {file}')
        for file in self.conf['runoff_depth_grids']:
            if not os.path.exists(file):
                raise FileNotFoundError(f'Runoff depth grid file not found at: {file}')
        for directory in set(os.path.dirname(os.path.abspath(file)) for file in self.conf['discharge_files']):
            if not os.path.exists(directory):
                raise NotADirectoryError(f'Output file directory not found at: {directory}')
        if self.conf['runoff_depth_grids']:
            if not self.conf.get('grid_weights_file'):
                raise ValueError('grid_weights_file is required when using runoff_depth_grids')
            if not os.path.exists(self.conf['grid_weights_file']):
                raise FileNotFoundError('grid_weights_file not found')

        if not self.conf['lateral_depth_files']:
            del self.conf['lateral_depth_files']
        if not self.conf['runoff_depth_grids']:
            del self.conf['runoff_depth_grids']

        if self._runoff_transformer is None and not self.conf.get('transformer_kernel_file'):
            raise RuntimeError(
                'No transformer configured. '
                'Provide transformer_kernel_file in config or call set_transformer() before route().'
            )

    def _initialize_runoff_transformer(self) -> None:
        if self._runoff_transformer is not None:
            return
        transformer_kernel_file = self.conf.get('transformer_kernel_file')
        if transformer_kernel_file:
            self._runoff_transformer = Transformer.from_kernel(dt=self.dt_runoff, kernel=transformer_kernel_file)
            if self.conf.get('transformer_state_file'):
                self._runoff_transformer.set_state(self.conf['transformer_state_file'])

    def _load_runoff(self) -> GeneratorSignature:
        if self.conf.get('lateral_depth_files', False):
            for depth_file, discharge_file in zip(self.conf['lateral_depth_files'], self.conf['discharge_files']):
                self.logger.info(f'Reading lateral depth file: {depth_file}')
                with xr.open_dataset(depth_file) as ds:
                    self.logger.debug('Reading time array')
                    dates = ds['time'].values.astype('datetime64[s]')
                    self.logger.debug('Reading depth array')
                    depths_array = ds[self.conf['var_runoff_depth']].values.astype(np.float64, copy=False)
                yield dates, depths_array, depth_file, discharge_file
        elif self.conf.get('runoff_depth_grids', False):
            for runoff_file, discharge_file in zip(self.conf['runoff_depth_grids'], self.conf['discharge_files']):
                self.logger.info(f'Calculating catchment volumes from runoff depth grid: {runoff_file}')
                volumes_ds = depth_to_volume(
                    runoff_file,
                    weight_table=self.conf['grid_weights_file'],
                    runoff_var=self.conf['var_runoff_depth'],
                    x_var=self.conf['var_x'],
                    y_var=self.conf['var_y'],
                    time_var=self.conf['var_t'],
                    river_id_var=self.conf['var_river_id'],
                    cumulative=self.conf['runoff_type'] == 'cumulative'
                )
                yield (volumes_ds['time'].values,
                       volumes_ds[self.conf['var_catchment_volume']].values.astype(np.float64, copy=False),
                       runoff_file, discharge_file)
        else:
            raise ValueError('No runoff data found in configs. Provide lateral_depth_files or runoff_depth_grids.')

    def route(self) -> Self:
        self.logger.info('Beginning routing')
        t1 = datetime.datetime.now()

        self._validate_configs()
        self._set_network_dependent_vectors()
        self.logger.debug(self)

        for dates, lateral_array, runoff_file, discharge_file in self._load_runoff():
            self.logger.info('-' * 80)
            self._set_network_and_time_dependent_vectors(dates)
            discharge_array = self._router(dates, lateral_array)

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
                dates = dates[::self.num_runoff_steps_per_discharge]

            self.logger.info('Writing Discharge Array to File')
            np.round(discharge_array, decimals=2, out=discharge_array)
            discharge_array = discharge_array.astype(np.float32, copy=False)
            self._write_discharges(dates, discharge_array, discharge_file, runoff_file)

        if self.conf['input_type'] == 'ensemble':
            self.channel_state = np.array(self._ensemble_member_states).mean(axis=0)
        self._write_final_state()
        if self.conf.get('final_transformer_state_file') and self._runoff_transformer is not None:
            self.logger.debug('Writing final transformer state to parquet')
            pd.DataFrame(self._runoff_transformer.state.T).to_parquet(self.conf['final_transformer_state_file'])

        t2 = datetime.datetime.now()
        self.logger.info('All runoff files routed')
        self.logger.info(f'Total job time: {(t2 - t1).total_seconds()}')
        return self

    def _router(self, dates: DatetimeArray, volumes: FloatArray) -> FloatArray:
        # todo why is dates not used?
        self.logger.debug('Getting initial state arrays')
        self._read_initial_state()
        q_init = self.channel_state
        self._ensemble_member_states = []

        T, n = volumes.shape
        discharge_array = np.zeros((T, n), dtype=np.float64)
        q_t = q_init.astype(np.float64, copy=True)
        rhs = np.zeros(n, dtype=np.float64)
        work = np.zeros(n, dtype=np.float64)
        interval_sum = np.zeros(n, dtype=np.float64)

        if self._runoff_transformer is None:
            raise RuntimeError('Runoff transformer has not been initialized')
        transformer = self._runoff_transformer

        t1 = datetime.datetime.now()
        runoff_iter = range(T)
        if self.conf['progress_bar']:
            runoff_iter = tqdm.tqdm(runoff_iter, desc='Unit-Muskingum Routing', total=T)
        else:
            self.logger.info('Performing Unit-Muskingum routing iterations')

        for t in runoff_iter:
            ql_t = transformer.transform(volumes[t, :])
            interval_sum.fill(0.0)
            for _ in range(self.num_routing_steps_per_runoff):
                work[:] = self.A @ q_t
                np.multiply(self.c1, work, out=rhs)
                np.multiply(self.c3, q_t, out=work)
                np.add(rhs, work, out=rhs)
                np.add(rhs, ql_t, out=rhs)
                q_t[:] = self.lhs_factorized(rhs)
                interval_sum += q_t
            discharge_array[t, :] = interval_sum / self.num_routing_steps_per_runoff

        discharge_array[discharge_array < 0] = 0

        if self.conf['input_type'] == 'sequential':
            self.logger.debug('Updating Channel State for Next Sequential Computation')
            self.channel_state = q_t
        if self.conf['input_type'] == 'ensemble':
            self.logger.debug('Recording Member State for Final State Aggregation')
            self._ensemble_member_states.append(q_t.copy())

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return discharge_array

    def _set_network_and_time_dependent_vectors(self, dates: DatetimeArray) -> None:
        self.logger.debug('Setting and validating time parameters')
        self.dt_runoff = self.conf.get('dt_runoff', (dates[1] - dates[0]).astype('timedelta64[s]').astype(int))
        self.dt_discharge = self.conf.get('dt_discharge', self.dt_runoff)
        self.dt_total = self.conf.get('dt_total', self.dt_runoff * dates.shape[0])
        if not self.conf.get('dt_routing', 0):
            self.logger.warning('dt_routing was not provided or is Null/False, defaulting to dt_runoff')
        self.dt_routing = self.conf.get('dt_routing', self.dt_runoff)

        signature = (self.dt_total, self.dt_runoff, self.dt_discharge, self.dt_routing,)
        if self._network_time_signature == signature:
            return

        try:
            assert self.dt_total >= self.dt_runoff, 'dt_total must be >= dt_runoff'
            assert self.dt_total >= self.dt_discharge, 'dt_total must be >= dt_discharge'
            assert self.dt_discharge >= self.dt_runoff, 'dt_discharge must be >= dt_runoff'
            assert self.dt_runoff >= self.dt_routing, 'dt_runoff must be >= dt_routing'
            assert self.dt_total % self.dt_runoff == 0, 'dt_total must be an integer multiple of dt_runoff'
            assert self.dt_total % self.dt_discharge == 0, 'dt_total must be an integer multiple of dt_discharge'
            assert self.dt_discharge % self.dt_runoff == 0, 'dt_discharge must be an integer multiple of dt_runoff'
            assert self.dt_runoff % self.dt_routing == 0, 'dt_runoff must be an integer multiple of dt_routing'
        except AssertionError as e:
            self.logger.error(e)
            raise AssertionError('Time options are not valid')

        self.num_runoff_steps_per_discharge = int(self.dt_discharge / self.dt_runoff)
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)
        self.num_routing_steps = int(self.dt_total / self.dt_routing)

        self._set_muskingum_coefficients(self.dt_routing)
        self._network_time_signature = signature
        self._initialize_runoff_transformer()
