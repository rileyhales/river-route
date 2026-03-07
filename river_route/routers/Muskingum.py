import datetime
import json
import logging
import sys
import traceback
from typing import Any, Self

import netCDF4 as nc
import numpy as np
import pandas as pd
import yaml

from .Config import Configs
from ._numba_kernels import muskingum_route
from ..logging import PROGRESS
from ..tools import adjacency_matrix
from ..types import IntArray, FloatArray, PathInput, WriteDischargesFn, DatetimeArray

__all__ = ['Muskingum', ]


class Muskingum:
    """
    Solves Muskingum routing using a matrix form for simultaneous solution on all rivers. Channel routing only.
    Inflow to each river segment is the sum of discharge from upstream segments.

    Q_t+1 = (c1 * I_t+1) + (c2 * I_t) + (c3 * Q_t)
    c1 = (dt/k - 2x) / (dt/k + 2(1-x))
    c2 = (dt/k + 2x) / (dt/k + 2(1-x))
    c3 = (2(1-x) - dt/k) / (dt/k + 2(1-x))

    In matrix form, the router needs to solve the equation:
    (I - c1 * A) @ Q_t+1 = c2 * (A @ Q_t) + c3 * Q_t

    The LHS matrix (I - c1 * A) is unit lower triangular (rows are topologically sorted).
    It is solved via numba JIT-compiled forward substitution operating on the CSC sparse structure.
    """
    cfg: Configs
    logger: logging.Logger

    # any subclass that requires non-null values for configs should have an iterable of the keys needed
    _ROUTER_REQUIRED_CONFIGS = ('channel_state_init_file', 'dt_routing', 'dt_total')

    # Network dependent matrices and vectors from routing parameters file
    A: object  # n x n - adjacency matrix (scipy csc_matrix)
    river_ids: IntArray  # n x 1 - river ID for each segment
    k: FloatArray  # n x 1 - K values for each segment
    x: FloatArray  # n x 1 - X values for each segment

    # Network and routing timestep dependent matrices and vectors
    c1: FloatArray  # n x 1 - C1 values for each segment => f(k, x, dt_routing)
    c2: FloatArray  # n x 1 - C2 values for each segment => f(k, x, dt_routing)
    c3: FloatArray  # n x 1 - C3 values for each segment => f(k, x, dt_routing)

    # State variables
    channel_state: FloatArray
    _network_time_signature: tuple[Any, ...] | None = None  # check if time params change between computes

    # Time options
    dt_routing: int  # compute time step, must be divisible into and <= min(dt_runoff, dt_discharge)
    dt_runoff: int  # time between lateral inflows, must be 1) constant, divisible into and <= dt_total
    dt_discharge: int  # time to average routed flows and save them, must be <=
    dt_total: int  # how long to simulate,

    # methods that are overridable via dependency injection
    _write_discharges: WriteDischargesFn

    # indices and data for numba solvers
    _csc_indptr: IntArray
    _csc_indices: IntArray
    _lhs_off_data: FloatArray

    def __init__(self, configs: PathInput | Configs | None = None, **kwargs: Any) -> None:
        # parse and createe
        raw: dict[str, Any] = {}
        if configs is not None and configs != '':
            if str(configs).endswith('.json'):
                with open(configs, 'r') as f:
                    raw = json.load(f)
            elif str(configs).endswith(('.yml', '.yaml')):
                with open(configs, 'r') as f:
                    raw = yaml.load(f, Loader=yaml.FullLoader)
            else:
                raise RuntimeError('Unrecognized simulation config file type. Must be .json or .yaml')
        raw.update(kwargs)
        raw.pop('_router', None)
        self.cfg = Configs(**raw)

        # configure logging - progress bar and info/debug logs are mutually exclusive
        self.logger = logging.getLogger(f'river_route.{id(self):x}')
        self.logger.disabled = not self.cfg.log
        self.logger.setLevel(self.cfg.log_level)
        if self.cfg.log_stream == 'stdout':
            self.logger.addHandler(logging.StreamHandler(sys.stdout))
        else:
            self.logger.addHandler(logging.FileHandler(self.cfg.log_stream))
        self.logger.handlers[0].setFormatter(logging.Formatter(self.cfg.log_format))
        self.logger.debug('Logger initialized')
        return

    def __repr__(self) -> str:
        return f'{type(self).__name__}(params_file={self.cfg.params_file!r})'

    def _validate_configs(self) -> None:
        self.logger.debug('Validating configs file')
        for key in self._ROUTER_REQUIRED_CONFIGS:
            if not getattr(self.cfg, key, None):
                raise ValueError(f'{key} is required for {type(self).__name__}')
        self._validate_router_configs()
        return

    def _validate_router_configs(self) -> None:
        """Subclass hook for relational validation beyond _ROUTER_REQUIRED_CONFIGS."""
        if len(self.cfg.discharge_files) != 1:
            raise ValueError('Muskingum requires exactly one entry in discharge_files')
        return

    ################################################
    # Computation state handling
    ################################################

    def _read_initial_state(self) -> None:
        if hasattr(self, 'channel_state'):
            return

        state_file = self.cfg.channel_state_init_file
        if not state_file:
            self.logger.warning('channel_state_init_file not provided. Defaulting to zero initial conditions')
            self.channel_state = np.zeros(self.A.shape[0], dtype=np.float64)
            return
        self.logger.debug('Reading Initial State from Parquet')
        self.channel_state = pd.read_parquet(state_file).values.flatten().astype(np.float64, copy=False)
        return

    def _write_final_state(self) -> None:
        final_state_file = self.cfg.channel_state_final_file
        if not final_state_file:
            return
        self.logger.debug('Writing Final State to Parquet')
        pd.DataFrame({'Q': self.channel_state}).to_parquet(self.cfg.channel_state_final_file)
        return

    ################################################
    # Prepare arrays for routing
    ################################################

    def _set_network_dependent_vectors(self) -> None:
        self.logger.debug('Calculating network dependent vectors')
        try:
            df = pd.read_parquet(
                self.cfg.params_file,
                columns=[self.cfg.var_river_id, 'k', 'x', 'downstream_river_id']
            )
        except Exception as e:
            self.logger.error(f'Error reading required parameter columns from params_file: {e}')
            self.logger.debug(traceback.format_exc())
            raise

        if df[self.cfg.var_river_id].duplicated().any():
            raise ValueError('params_file contains duplicate river IDs.')

        self.river_ids = df[self.cfg.var_river_id].to_numpy(dtype=np.int64, copy=False)
        downstream_river_ids = df['downstream_river_id'].to_numpy(dtype=np.int64, copy=False)
        self.k = df['k'].to_numpy(dtype=np.float64, copy=False)
        self.x = df['x'].to_numpy(dtype=np.float64, copy=False)

        river_id_set = set(self.river_ids.tolist())
        downstream_ids = {d for d in downstream_river_ids.tolist() if d > 0}
        unknown_downstream_ids = sorted(downstream_ids - river_id_set)
        if unknown_downstream_ids:
            raise ValueError(
                f'params_file has downstream IDs not in river_id column: {unknown_downstream_ids[:10]}')
        self.A = adjacency_matrix(self.river_ids, downstream_river_ids)
        self.logger.log(PROGRESS, f'Network: {self.A.shape[0]} river segments')
        return

    def _set_muskingum_coefficients(self, dt_routing: float) -> None:
        self.logger.debug('Calculating Muskingum coefficients')
        dt_div_k = dt_routing / self.k
        denominator = dt_div_k + (2 * (1 - self.x))
        _2x = 2 * self.x
        self.c1 = (dt_div_k - _2x) / denominator
        self.c2 = (dt_div_k + _2x) / denominator
        self.c3 = ((2 * (1 - self.x)) - dt_div_k) / denominator
        if not np.allclose(self.c1 + self.c2 + self.c3, 1):
            self.logger.warning('Muskingum coefficients do not sum to 1')
            self.logger.debug(f'c1: {self.c1}')
            self.logger.debug(f'c2: {self.c2}')
            self.logger.debug(f'c3: {self.c3}')
            raise ValueError('Muskingum coefficients do not sum to 1, check routing parameters and time step')

        # LHS = I - diags(c1) @ A is unit lower triangular.
        # Off-diagonal entries share A's sparsity: data = -c1[row] per edge.
        A_csc = self.A.tocsc()
        self._csc_indptr = A_csc.indptr
        self._csc_indices = A_csc.indices
        self._lhs_off_data = np.ascontiguousarray(-self.c1[A_csc.indices])
        return

    ################################################
    # Methods to execute routing simulation
    ################################################

    def route(self) -> Self:
        """
        Execute the simulation described by the provided configs and routing parameters. All configs, file paths,
        parameters, and options must be set when the object is initialized so that validation is performed before the
        simulation.

        Returns:
            Self: the class instance with updated channel_state and output files written to disk
        """
        # start timer
        self.logger.log(PROGRESS, 'Beginning routing')
        t1 = datetime.datetime.now()
        # validate configuration options
        self._validate_configs()
        self.logger.debug(self)
        # set arrays for routing
        self._set_network_dependent_vectors()
        self._read_initial_state()
        # init hook
        self._hook_before_route()
        # routing handled by subclass routing logic
        self._execute_routing()
        # final hook
        self._write_final_state()
        self._hook_after_route()
        # log total time
        t2 = datetime.datetime.now()
        self.logger.log(PROGRESS, f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return self

    def _execute_routing(self) -> None:
        self.logger.info('-' * 60)
        # time parameters
        self.dt_routing = self.cfg.dt_routing
        self.dt_total = self.cfg.dt_total
        self.dt_discharge = self.cfg.dt_discharge or self.dt_routing
        if not (self.dt_total >= self.dt_discharge >= self.dt_routing):
            raise ValueError('Need dt_total >= dt_discharge >= dt_routing')
        if self.dt_total % self.dt_discharge != 0:
            raise ValueError('dt_total must be an integer multiple of dt_discharge')
        if self.dt_discharge % self.dt_routing != 0:
            raise ValueError('dt_discharge must be an integer multiple of dt_routing')
        num_output_steps = int(self.dt_total / self.dt_discharge)
        num_routing_per_output = int(self.dt_discharge / self.dt_routing)
        self._set_muskingum_coefficients(self.dt_routing)

        self.logger.debug('Starting routing computation')
        discharge_array = self._router(num_output_steps, num_routing_per_output)

        # Generate date array for output
        dates = pd.date_range(
            start=self.cfg.start_datetime,
            periods=num_output_steps,
            freq=pd.to_timedelta(self.dt_discharge, unit='s')
        ).to_numpy()

        # write outputs
        self.logger.debug('Writing Discharge Array to File')
        np.round(discharge_array, decimals=2, out=discharge_array)
        discharge_array = discharge_array.astype(np.float32, copy=False)
        self._write_discharges(dates, discharge_array, self.cfg.discharge_files[0])
        self.logger.info('-' * 60)
        return

    def _router(self, num_output_steps: int, num_routing_per_output: int) -> FloatArray:
        """
        Route discharge without lateral inflow.

        Muskingum channel routing equation:
            (I - c1*A) @ Q(t+1) = c2*(A @ Q(t)) + c3*Q(t)
        """
        self.logger.debug('Getting initial state arrays')
        q_init = self.channel_state
        if not np.any(q_init):
            self.logger.warning(
                'Initial channel state is all zeros. Muskingum routing without lateral inflow requires a '
                'non-zero initial state to produce meaningful results. Provide channel_state_init_file.'
            )

        n = self.A.shape[0]
        discharge_array = np.zeros((num_output_steps, n), dtype=np.float64)
        q_t = q_init.astype(np.float64, copy=True)

        muskingum_route(
            self._csc_indptr, self._csc_indices, self._lhs_off_data,
            self.c2, self.c3, q_t,
            discharge_array,
            num_output_steps, num_routing_per_output,
        )

        discharge_array[discharge_array < 0] = 0

        self.logger.debug('Updating Channel State')
        self.channel_state = q_t
        return discharge_array

    ################################################
    # Hooks so that subclasses can cleanly inject behavior with less overriding or duplicating
    ################################################

    def _hook_before_route(self) -> None:
        """Called after validation and network setup, before routing begins. Default: no-op."""
        return

    def _hook_after_route(self) -> None:
        """Called at the end of route(), after writing final state. Default: no-op."""
        return

    ################################################
    # Dependency injection methods for users to overwrite default behaviors without subclassing
    ################################################

    def set_write_discharges(self, func: WriteDischargesFn) -> Self:
        """
        Overwrites the default write_discharges method to a custom function and returns the class instance so that you
        can chain the method with the constructor.

        Args:
            func (callable): function that takes dates, discharge_array, discharge_file, runoff_file and returns None
        """
        self._write_discharges = func
        return self

    def _write_discharges(self,
                          dates: DatetimeArray,
                          q_array: FloatArray,
                          q_file: PathInput,
                          routed_file: PathInput = '', ) -> None:
        """
        Writes routed discharge from a routing simulation to a netcdf file.
        You can overwrite this method with a custom handler using set_write_discharges.

        Args:
            dates: datetime array corresponding to the discharge rows
            q_array: routed discharge values with shape (time, river)
            q_file: path to write the discharge data to
            routed_file: path to the lateral inflow used to generate the discharge values, if applicable.

        Returns:
            None
        """
        with nc.Dataset(str(q_file), mode='w', format='NETCDF4') as ds:
            ds.createDimension('time', size=q_array.shape[0])
            ds.createDimension(self.cfg.var_river_id, size=q_array.shape[1])
            ds.runoff_file = str(routed_file)
            time_var = ds.createVariable('time', 'f8', ('time',))
            time_var.units = f'seconds since {pd.Timestamp(dates[0]).strftime("%Y-%m-%d %H:%M:%S")}'
            time_var[:] = (dates - dates[0]).astype('timedelta64[s]').astype(np.int64)
            id_var = ds.createVariable(self.cfg.var_river_id, 'i4', self.cfg.var_river_id, )
            id_var[:] = self.river_ids
            flow_var = ds.createVariable(self.cfg.var_discharge, 'f4', ('time', self.cfg.var_river_id))
            flow_var[:] = q_array
            flow_var.long_name = 'Discharge at catchment outlet'
            flow_var.standard_name = 'discharge'
            flow_var.aggregation_method = 'mean'
            flow_var.units = 'm3 s-1'
        return
