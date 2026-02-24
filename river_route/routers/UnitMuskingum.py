import datetime
from typing import Any, Self

import numpy as np
import pandas as pd
import tqdm
import xarray as xr

from .AbstractTransformRouter import AbstractTransformRouter
from .types import DatetimeArray, FloatArray, PathInput
from ..transformers import AbstractBaseTransformer, Transformer

__all__ = ['UnitMuskingum']


class UnitMuskingum(AbstractTransformRouter):
    """
    A Muskingum router that applies a unit hydrograph transformation to lateral depth inputs before routing.

    Runoff depths are first convolved with a unit hydrograph kernel (via a Transformer) to distribute overland
    flow through time before it enters the channel network. The transformed lateral inflow is then routed between
    river segments using the Muskingum method. This contrasts with RapidMuskingum, which assumes runoff enters
    the channel uniformly over its input timestep with no overland travel time.

    The Transformer can be supplied in two ways:
      1. Set ``transformer_kernel_file`` in config — a Transformer is built automatically from the kernel file.
      2. Call ``set_transformer()`` before ``route()`` — inject a pre-built or custom ``AbstractBaseTransformer``.

    Required configs:
    - routing_params_file: path to a parquet file containing Muskingum routing parameters for each river segment.

    Lateral input — one of:
    - lateral_depth_files: list of 1 or more paths to netCDF files containing per-catchment runoff depth time
        series for each river segment.
    - runoff_depth_grids: list of 1 or more paths to netCDF files containing gridded runoff depth time series.
        Must be accompanied by grid_weights_file.
    - grid_weights_file: path to a netCDF file of weights for remapping gridded depths to catchment depths.
        Required when using runoff_depth_grids.

    - discharge_files: list of paths where routed discharge netCDF files will be written. Must match the number
        of input files.

    Transformer configs:
    - transformer_kernel_file: path to a file defining the unit hydrograph kernel used to build the Transformer.
        Required unless set_transformer() is called.
    - transformer_state_file: (optional) path to a parquet file containing the initial transformer convolution
        state. Used to warm-start the transformer, e.g. when chaining sequential simulations.
    - final_transformer_state_file: (optional) path where the final transformer convolution state will be written
        as a parquet file after routing completes.
    """
    def __init__(self, config_file: PathInput | None = None, **kwargs: Any) -> None:
        super().__init__(config_file, **kwargs)

    @property
    def _lateral_water_files_key(self) -> str:
        return 'lateral_depth_files'

    def set_transformer(self, transformer: AbstractBaseTransformer) -> Self:
        """
        Inject a user instantiated transformer object including custom subclasses of AbstractBaseTransformer allowing
        users to implement their own transformation logic in a custom transformer. The transformer must be a
        AbstractBaseTransformer subclass with a kernel and state set.

        You must call this before route() to have any affect.
        """
        if not isinstance(transformer, AbstractBaseTransformer):
            raise TypeError(f'transformer must be a AbstractBaseTransformer, got {type(transformer).__name__}')
        self._Transformer = transformer
        return self

    def _validate_router_configs(self) -> None:
        self._validate_lateral_runoff_configs()
        if self._Transformer is None and not self.conf.get('transformer_kernel_file'):
            raise RuntimeError(
                'No transformer configured. '
                'Provide transformer_kernel_file in config or call set_transformer() before route().'
            )

    def _read_lateral_file(self, file_path: PathInput):
        self.logger.info(f'Reading lateral depth file: {file_path}')
        with xr.open_dataset(file_path) as ds:
            return (ds['time'].values.astype('datetime64[s]'),
                    ds[self.conf['var_runoff_depth']].values.astype(np.float64, copy=False))

    def _on_runoff_file_start(self) -> None:
        self._initialize_runoff_transformer()

    def _on_after_route(self) -> None:
        if self.conf.get('final_transformer_state_file') is not None and self._Transformer is not None:
            self.logger.debug('Writing final transformer state to parquet')
            pd.DataFrame(self._Transformer.state.T).to_parquet(self.conf['final_transformer_state_file'])

    def _initialize_runoff_transformer(self) -> None:
        if self._Transformer is not None:
            return
        transformer_kernel_file = self.conf.get('transformer_kernel_file')
        if transformer_kernel_file:
            self._Transformer = Transformer.from_kernel(dt=self.dt_runoff, kernel=transformer_kernel_file)
            if self.conf.get('transformer_state_file'):
                self._Transformer.set_state(self.conf['transformer_state_file'])

    def _router(self, dates: DatetimeArray, volumes: FloatArray) -> FloatArray:
        self.logger.debug('Getting initial state arrays')
        self._read_initial_state()
        q_init = self.channel_state
        self._ensemble_member_states = []

        discharge_array = np.zeros((volumes.shape[0], self.A.shape[0]), dtype=np.float64)
        q_t = np.empty_like(q_init, dtype=np.float64)
        rhs = np.zeros(self.A.shape[0], dtype=np.float64)
        work = np.zeros(self.A.shape[0], dtype=np.float64)
        interval_sum = np.zeros(self.A.shape[0], dtype=np.float64)
        q_t[:] = q_init

        if self._Transformer is None:
            raise RuntimeError('Runoff transformer has not been initialized')
        transformer = self._Transformer

        t1 = datetime.datetime.now()
        if self.conf['progress_bar']:
            dates = tqdm.tqdm(dates, desc='Lateral Depths Routed')
        else:
            self.logger.info('Performing routing computation iterations')

        for runoff_time_step, runoff_end_date in enumerate(dates):
            ql_t = transformer.transform(volumes[runoff_time_step, :])
            interval_sum.fill(0.0)
            for _ in range(self.num_routing_steps_per_runoff):
                work[:] = self.A @ q_t
                np.multiply(self.c1, work, out=rhs)
                np.multiply(self.c3, q_t, out=work)
                np.add(rhs, work, out=rhs)
                np.add(rhs, ql_t, out=rhs)
                q_t[:] = self.lhs_factorized(rhs)
                interval_sum += q_t
            discharge_array[runoff_time_step, :] = interval_sum / self.num_routing_steps_per_runoff

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
