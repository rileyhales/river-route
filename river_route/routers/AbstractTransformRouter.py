import datetime
import os
from abc import abstractmethod
from pathlib import Path
from typing import Generator, Self, Tuple, List

import numpy as np

from .AbstractRouter import AbstractRouter
from .types import DatetimeArray, FloatArray, PathInput
from ..runoff import depth_to_volume
from ..transformers import AbstractBaseTransformer

__all__ = ['AbstractTransformRouter', ]

GeneratorSignature = Generator[Tuple[DatetimeArray, FloatArray, PathInput, PathInput], None, None]


class AbstractTransformRouter(AbstractRouter):
    """
    Abstract intermediate class for routers that receive lateral water (runoff) as input alongside the channel network.

    Extends AbstractRouter and provides shared config validation, runoff data loading, and the route() method.
    Remains abstract because _validate_router_configs (from AbstractRouter), _router, and _read_lateral_file
    are left for concrete subclasses to implement.

    Concrete subclasses must implement:
      - _validate_router_configs()  — from AbstractRouter; call _validate_lateral_runoff_configs() inside
      - _router(dates, lateral_array) -> FloatArray  — the core routing math
      - _read_lateral_file(file_path) -> Tuple[DatetimeArray, FloatArray]  — file format handling

    Optional overrides (hooks with default no-ops or defaults):
      - _lateral_water_files_key  — config key for per-catchment lateral files (default: 'lateral_volume_files')
      - _prepare_lateral_array(array) -> FloatArray  — pre-processing before _router (default: pass-through)
      - _on_runoff_file_start()  — called before _router on each file iteration (default: no-op)
    """
    # State variables
    _ensemble_member_states: List[FloatArray]  # for ensemble routing, stores member states for computing final state

    # Time options
    num_routing_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_discharge: int

    # for unit hydrograph style convolutions
    _Transformer: AbstractBaseTransformer | None = None

    @property
    def _lateral_water_files_key(self) -> str:
        """Config key for pre-processed per-catchment lateral input files."""
        return 'lateral_volume_files'

    @abstractmethod
    def _read_lateral_file(self, file_path: PathInput) -> Tuple[DatetimeArray, FloatArray]:
        """Open a pre-processed lateral input file and return (dates, array)."""

    def _prepare_lateral_array(self, array: FloatArray) -> FloatArray:
        """Pre-process the lateral array before passing to _router. Default: pass-through."""
        return array

    def _on_runoff_file_start(self) -> None:
        """Called at the start of each per-file routing iteration, after time setup. Default: no-op."""
        pass

    @abstractmethod
    def _router(self, dates: DatetimeArray, lateral: FloatArray) -> FloatArray:
        """Execute the core routing math for one runoff file and return the discharge array."""

    def _validate_lateral_runoff_configs(self) -> None:
        key = self._lateral_water_files_key
        for k in (key, 'runoff_depth_grids', 'discharge_files'):
            self.conf[k] = self.conf.get(k, [])
            if isinstance(self.conf[k], (str, Path)):
                self.conf[k] = [self.conf[k]]

        self.conf['runoff_depth_grids'] = [os.path.abspath(p) for p in self.conf['runoff_depth_grids']]

        lateral = self.conf[key]
        grids = self.conf['runoff_depth_grids']
        n_files = len(lateral) + len(grids)

        if lateral and grids:
            raise ValueError(f'Provide either {key} or runoff_depth_grids, not both')
        if len(self.conf['discharge_files']) != n_files:
            raise ValueError('Number of discharge_files must match number of input files')
        for f in lateral:
            if not os.path.exists(f):
                raise FileNotFoundError(f'Lateral file not found: {f}')
        for f in grids:
            if not os.path.exists(f):
                raise FileNotFoundError(f'Runoff depth grid not found: {f}')
        for directory in {os.path.dirname(os.path.abspath(f)) for f in self.conf['discharge_files']}:
            if not os.path.exists(directory):
                raise NotADirectoryError(f'Output directory not found: {directory}')
        if grids:
            if not self.conf.get('grid_weights_file'):
                raise ValueError('grid_weights_file is required when using runoff_depth_grids')
            if not os.path.exists(self.conf['grid_weights_file']):
                raise FileNotFoundError('grid_weights_file not found')

        if not lateral:
            del self.conf[key]
        if not grids:
            del self.conf['runoff_depth_grids']

    def _lateral_runoff_generator(self) -> GeneratorSignature:
        key = self._lateral_water_files_key
        if self.conf.get(key):
            for lateral_file, discharge_file in zip(self.conf[key], self.conf['discharge_files']):
                dates, array = self._read_lateral_file(lateral_file)
                yield dates, array, lateral_file, discharge_file
        elif self.conf.get('runoff_depth_grids'):
            for runoff_file, discharge_file in zip(self.conf['runoff_depth_grids'], self.conf['discharge_files']):
                self.logger.info(f'Calculating catchment volumes from runoff depth grid: {runoff_file}')
                ds = depth_to_volume(
                    runoff_file,
                    weight_table=self.conf['grid_weights_file'],
                    runoff_var=self.conf['var_runoff_depth'],
                    x_var=self.conf['var_x'],
                    y_var=self.conf['var_y'],
                    time_var=self.conf['var_t'],
                    river_id_var=self.conf['var_river_id'],
                    cumulative=self.conf['runoff_type'] == 'cumulative',
                )
                yield (ds['time'].values,
                       ds[self.conf['var_catchment_volume']].values.astype(np.float64, copy=False),
                       runoff_file, discharge_file)
        else:
            raise ValueError(f'No runoff data in configs. Provide {key} or runoff_depth_grids.')

    def route(self) -> Self:
        self.logger.info('Beginning routing')
        t1 = datetime.datetime.now()

        self._on_before_route()
        self._validate_configs()
        self._set_network_dependent_vectors()
        self.logger.debug(self)

        self._ensemble_member_states = []

        for dates, lateral_array, runoff_file, discharge_file in self._lateral_runoff_generator():
            self.logger.info('-' * 80)
            self._set_network_and_time_dependent_vectors(dates)
            self._on_runoff_file_start()
            lateral_array = self._prepare_lateral_array(lateral_array)
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
        self._on_after_route()

        t2 = datetime.datetime.now()
        self.logger.info('All runoff files routed')
        self.logger.info(f'Total job time: {(t2 - t1).total_seconds()}')
        return self
