import numpy as np
import xarray as xr
from tqdm import tqdm

from .Muskingum import Muskingum
from ._numba_kernels import static_muskingum_vlateral
from ..runoff import runoff_to_qlateral
from ..types import FloatArray, DatetimeArray, QlateralGeneratorSignature

__all__ = ['MuskingumVlateral', ]


class MuskingumVlateral(Muskingum):
    """
    Muskingum channel routing with direct lateral inflow. Lateral flow is the runoff volume divided by the runoff
    timestep — all runoff enters the channel in the interval it is generated, ignoring overland flow delay.

    See the Math Derivations page in the documentation for the full equations.
    """
    _ROUTER_REQUIRED_CONFIGS = ()

    # State variables
    _ensemble_member_states: list[FloatArray]  # for ensemble routing

    # Time options
    num_runoff_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_discharge: int

    _as_volumes: bool = True

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._validate_router_configs()

    def _qlateral_generator(self) -> QlateralGeneratorSignature:
        if self.cfg.qlateral_files:
            for lateral_file, discharge_file in zip(self.cfg.qlateral_files, self.cfg.discharge_files):
                self.logger.info('-' * 60)
                with xr.open_dataset(lateral_file) as ds:
                    dates = ds['time'].values.astype('datetime64[s]')
                    array = ds['qlateral'].values.astype(np.float32, copy=False)
                    yield dates, array, lateral_file, discharge_file
        elif self.cfg.grid_runoff_files and self.cfg.grid_weights_file:
            for runoff_file, discharge_file in zip(self.cfg.grid_runoff_files, self.cfg.discharge_files):
                self.logger.info('-' * 60)
                self.logger.debug(f'Calculating qlateral: {runoff_file}')
                ds = runoff_to_qlateral(runoff_file, grid_weights_file=self.cfg.grid_weights_file,
                                        var_runoff=self.cfg.var_grid_runoff, var_x=self.cfg.var_x, var_y=self.cfg.var_y,
                                        var_t=self.cfg.var_t, var_river_id=self.cfg.var_river_id,
                                        cumulative=self.cfg.grid_accumulation_type == 'cumulative',
                                        as_volumes=self._as_volumes)
                yield (
                    ds['time'].values.astype('datetime64[s]'),
                    ds['qlateral'].values.astype(np.float32, copy=False),
                    runoff_file, discharge_file
                )

    def _validate_router_configs(self) -> None:
        qlateral = self.cfg.qlateral_files
        grids = self.cfg.grid_runoff_files and self.cfg.grid_weights_file

        if qlateral and grids:
            raise ValueError('Provide qlateral_files or grid_runoff_files with grid_weights_file, not both')
        if not qlateral and not grids:
            raise ValueError('Provide qlateral_files or grid_runoff_files with grid_weights_file')
        n_inputs = len(qlateral) + len(self.cfg.grid_runoff_files or [])
        if len(self.cfg.discharge_files) != n_inputs:
            raise ValueError('Number of resolved discharge output files must match number of input files')
        return

    def _set_network_and_time_dependent_vectors(self, dates: DatetimeArray) -> None:
        self.logger.debug('Setting and validating time parameters')
        self.dt_runoff = self.cfg.dt_runoff or (dates[1] - dates[0]).astype('timedelta64[s]').astype(int)
        self.dt_discharge = self.cfg.dt_discharge or self.dt_runoff
        self.dt_total = self.cfg.dt_total or self.dt_runoff * dates.shape[0]
        if not self.cfg.dt_routing:
            self.logger.warning('dt_routing was not provided or is Null/False, defaulting to dt_runoff')
        self.dt_routing = self.cfg.dt_routing or self.dt_runoff

        signature = (self.dt_total, self.dt_runoff, self.dt_discharge, self.dt_routing,)
        if self._network_time_signature == signature:
            return

        # check that time options have the correct sizes
        if self.dt_total < self.dt_runoff:
            raise ValueError('dt_total must be >= dt_runoff')
        if self.dt_total < self.dt_discharge:
            raise ValueError('dt_total must be >= dt_discharge')
        if self.dt_discharge < self.dt_runoff:
            raise ValueError('dt_discharge must be >= dt_runoff')
        if self.dt_runoff < self.dt_routing:
            raise ValueError('dt_runoff must be >= dt_routing')
        # check that time options are evenly divisible
        if self.dt_total % self.dt_runoff != 0:
            raise ValueError('dt_total must be an integer multiple of dt_runoff')
        if self.dt_total % self.dt_discharge != 0:
            raise ValueError('dt_total must be an integer multiple of dt_discharge')
        if self.dt_discharge % self.dt_runoff != 0:
            raise ValueError('dt_discharge must be an integer multiple of dt_runoff')
        if self.dt_runoff % self.dt_routing != 0:
            raise ValueError('dt_runoff must be an integer multiple of dt_routing')

        # set derived datetime parameters for computation cycles later
        self.num_runoff_steps = int(self.dt_total / self.dt_runoff)
        self.num_runoff_steps_per_discharge = int(self.dt_discharge / self.dt_runoff)
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)

        self._set_muskingum_coefficients(self.dt_routing)
        self.c4_dt = np.ascontiguousarray(self.c4 / self.dt_runoff, dtype=np.float32)
        self._network_time_signature = signature
        return

    def _execute_routing(self) -> None:
        self._ensemble_member_states = []

        total_files = len(self.cfg.qlateral_files or self.cfg.grid_runoff_files)
        file_iter = self._qlateral_generator()
        if self.cfg.progress_bar:
            file_iter = tqdm(file_iter, total=total_files, desc='Files Routed')

        for dates, qlateral, runoff_file, discharge_file in file_iter:
            self.logger.info(f'Routing qlateral: {runoff_file}')
            self._set_network_and_time_dependent_vectors(dates)
            self.logger.debug('Starting routing computation')
            q_t, q_array = self._router(qlateral)
            if self.cfg.runoff_processing_mode == 'sequential':
                self.logger.debug('Updating Channel State for Next Sequential Computation')
                self.channel_state = q_t
            elif self.cfg.runoff_processing_mode == 'ensemble':
                self.logger.debug('Recording Member State for Final State Aggregation')
                self._ensemble_member_states.append(q_t.copy())

            if self.dt_discharge > self.dt_runoff:
                self.logger.debug('Resampling dates and discharges to specified timestep')
                q_array = (
                    q_array
                    .reshape((
                        int(self.dt_total / self.dt_discharge),
                        int(self.dt_discharge / self.dt_runoff),
                        self.river_ids.shape[0],
                    ))
                    .mean(axis=1)
                )
                dates = dates[::self.num_runoff_steps_per_discharge]

            self.logger.debug('Writing Discharge Array to File')
            q_array = q_array.astype(np.float32, copy=False)
            self._write_discharges(dates, q_array, discharge_file, runoff_file)

        if self.cfg.runoff_processing_mode == 'ensemble':
            self.channel_state = np.array(self._ensemble_member_states).mean(axis=0)
        self.logger.info('-' * 60)
        return

    def _router(self, qlateral: FloatArray) -> tuple[FloatArray, FloatArray]:
        """Execute the core routing math for one runoff file and return the discharge array"""
        self.logger.debug('Getting initial state arrays')
        n = self.river_ids.shape[0]
        discharge_array = np.zeros((self.num_runoff_steps, n), dtype=np.float32)
        q_t = self.channel_state.astype(np.float32, copy=True)
        qlateral = np.ascontiguousarray(qlateral, dtype=np.float32)

        static_muskingum_vlateral(
            q_t=q_t,
            discharge_array=discharge_array,
            downstream_indices=self.downstream_indices,
            downstream_c1=self.downstream_c1,
            downstream_c2=self.downstream_c2,
            c3=self.c3,
            n_rivers=n,
            n_steps=self.num_runoff_steps,
            n_substeps=self.num_routing_steps_per_runoff,
            qlateral=qlateral,
            c4_dt=self.c4_dt
        )
        return q_t, discharge_array
