import datetime
from typing import Tuple

import numpy as np
import pandas as pd
import tqdm
import xarray as xr

from .TransformRouter import TransformRouter
from ..types import DatetimeArray, FloatArray, PathInput

__all__ = ['UnitMuskingum', ]


class UnitMuskingum(TransformRouter):
    """
    Muskingum router that applies a unit hydrograph convolution to lateral depth inputs.

    Runoff depths are convolved with a precomputed unit hydrograph kernel each timestep to
    distribute overland flow through time before it enters the channel network. The transformed
    lateral inflow is then routed between river segments using the Muskingum method.

    Required configs:
    - params_file: path to routing parameters parquet file.
    - transformer_kernel_file: path to a parquet kernel file (n_basins × n_time_steps).

    Lateral input — one of:
    - catchment_runoff_files: list of paths to netCDF files containing per-catchment runoff depths.
    - runoff_grid_files + grid_weights_file: gridded runoff depth inputs remapped to catchments.

    - discharge_files: list of output paths, one per input file.
    - dt_routing: routing sub-step in seconds (required).

    State configs:
    - transformer_state_init_file: (optional) path to warm-start the convolution state.
    - transformer_state_final_file: (optional) path to write the final convolution state.
    """

    _uh_kernel: FloatArray | None = None
    _uh_state: FloatArray | None = None

    @property
    def _catchment_runoff_as_volume(self) -> bool:
        return False  # False -> means as depth

    def _validate_router_configs(self) -> None:
        self._validate_lateral_runoff_configs()
        if not self.cfg.transformer_kernel_file:
            raise RuntimeError('transformer_kernel_file is required for UnitMuskingum routing')

    def _read_lateral_file(self, file_path: PathInput) -> Tuple[DatetimeArray, FloatArray]:
        self.logger.info(f'Reading lateral depth file: {file_path}')
        with xr.open_dataset(file_path) as ds:
            return (ds['time'].values.astype('datetime64[s]'),
                    ds[self.cfg.var_runoff_depth].values.astype(np.float64, copy=False))

    def _on_runoff_file_start(self) -> None:
        self._read_uh_kernel()

    def _hook_after_route(self) -> None:
        if self.cfg.transformer_state_final_file and self._uh_state is not None:
            self.logger.debug('Writing final convolution state to parquet')
            pd.DataFrame(self._uh_state.T).to_parquet(self.cfg.transformer_state_final_file)

    def _read_uh_kernel(self) -> None:
        if self._uh_kernel is not None:
            return
        self.logger.debug('Reading UH kernel from parquet')
        self._uh_kernel = pd.read_parquet(self.cfg.transformer_kernel_file).T.to_numpy(dtype=np.float64)
        state_file = self.cfg.transformer_state_init_file
        if state_file:
            self.logger.debug('Reading convolution state from parquet')
            _state = pd.read_parquet(state_file).T.to_numpy(dtype=np.float64, copy=True)
            if _state.shape != self._uh_kernel.shape:
                raise ValueError(f'State shape {_state.shape} does not match kernel shape {self._uh_kernel.shape}')
            self._uh_state = _state
        else:
            self._uh_state = np.zeros_like(self._uh_kernel, dtype=np.float64)

    def _router(self, dates: DatetimeArray, lateral: FloatArray) -> FloatArray:
        self.logger.debug('Getting initial state arrays')
        self._read_initial_state()
        q_init = self.channel_state

        if self._uh_kernel is None or self._uh_state is None:
            raise RuntimeError('UH kernel has not been initialized')

        n = self.A.shape[0]
        discharge_array = np.zeros((lateral.shape[0], n), dtype=np.float64)
        q_t = q_init.astype(np.float64, copy=True)
        rhs = np.zeros(n, dtype=np.float64)
        buffer = np.zeros(n, dtype=np.float64)
        interval_sum = np.zeros(n, dtype=np.float64)

        t1 = datetime.datetime.now()
        runoff_iter = tqdm.tqdm(dates, desc='Lateral Depths Routed') if self.cfg.progress_bar else dates
        if not self.cfg.progress_bar:
            self.logger.info('Performing routing computation iterations')

        for runoff_time_step, _ in enumerate(runoff_iter):
            # convolve: accumulate kernel response for this runoff depth, then advance state
            self._uh_state += self._uh_kernel * lateral[runoff_time_step, :]
            ql_t = self._uh_state[0, :].copy()
            self._uh_state[:-1, :] = self._uh_state[1:, :]
            self._uh_state[-1, :] = 0.0

            interval_sum.fill(0.0)
            for _ in range(self.num_routing_steps_per_runoff):
                # rhs = c1*(A @ q_t) + c3*q_t + ql_t
                buffer[:] = self.A @ q_t
                np.multiply(self.c1, buffer, out=rhs)
                np.multiply(self.c3, q_t, out=buffer)
                np.add(rhs, buffer, out=rhs)
                np.add(rhs, ql_t, out=rhs)
                q_t[:] = self.lhs_factorized(rhs)
                interval_sum += q_t
            discharge_array[runoff_time_step, :] = interval_sum / self.num_routing_steps_per_runoff

        discharge_array[discharge_array < 0] = 0

        if self.cfg.runoff_processing_mode == 'sequential':
            self.logger.debug('Updating Channel State for Next Sequential Computation')
            self.channel_state = q_t
        elif self.cfg.runoff_processing_mode == 'ensemble':
            self.logger.debug('Recording Member State for Final State Aggregation')
            self._ensemble_member_states.append(q_t.copy())

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return discharge_array
