import numpy as np
import pandas as pd

from .TransformMuskingum import TransformMuskingum
from ..types import FloatArray

__all__ = ['UnitMuskingum', ]


class UnitMuskingum(TransformMuskingum):
    """
    Muskingum router that applies a unit hydrograph convolution to lateral depth inputs.

    Runoff depths are convolved with a precomputed unit hydrograph kernel each timestep to
    distribute overland flow through time before it enters the channel network. The transformed
    lateral inflow is then routed between river segments using the Muskingum method.

    Required configs:
    - params_file: path to routing parameters parquet file.
    - discharge_dir: directory where output netCDF files are written (named after input files).
    - transformer_kernel_file: path to a parquet kernel file (n_basins × n_time_steps).
    Lateral input — one of:
    - catchment_runoff_files: list of paths to netCDF files containing per-catchment runoff depths.
    - runoff_grid_files + grid_weights_file: gridded runoff depth inputs remapped to catchments.
    """
    _ROUTER_REQUIRED_CONFIGS = ('transformer_kernel_file',)
    _uh_kernel: FloatArray | None = None
    _uh_state: FloatArray | None = None

    @property
    def _catchment_runoff_as_volume(self) -> bool:
        return False  # False -> means as depth

    def _hook_before_route(self) -> None:
        # Load the unit hydrograph kernel before routing begins. Kernel is static so it only needs to load once.
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

    def _write_final_state(self) -> None:
        super()._write_final_state()
        if self.cfg.transformer_state_final_file and self._uh_state is not None:
            self.logger.debug('Writing final convolution state to parquet')
            pd.DataFrame(self._uh_state.T).to_parquet(self.cfg.transformer_state_final_file)

    def transform_runoff(self, r_t: FloatArray) -> FloatArray:
        """Convolve runoff depth with the unit hydrograph kernel and advance convolution state."""
        if self._uh_kernel is None or self._uh_state is None:
            raise RuntimeError('UH kernel has not been initialized')
        self._uh_state += self._uh_kernel * r_t
        ql_t = self._uh_state[0, :].copy()
        self._uh_state[:-1, :] = self._uh_state[1:, :]
        self._uh_state[-1, :] = 0.0
        return ql_t
