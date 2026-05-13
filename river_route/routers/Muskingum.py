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
from ._numba_kernels import static_muskingum
from ..logging import PROGRESS
from ..types import IntArray, FloatArray, PathInput, WriteDischargesFn, DatetimeArray

__all__ = ['Muskingum', ]


class Muskingum:
    """
    Muskingum channel routing using a matrix formulation solved by forward substitution.
    Inflow to each river segment is the sum of discharge from upstream segments.
    No lateral inflow — channel routing only.

    See the Math Derivations and Forward Substitution pages in the documentation for the
    full equations and algorithm details.
    """
    cfg: Configs
    logger: logging.Logger

    # any subclass that requires non-null values for configs should have an iterable of the keys needed
    _ROUTER_REQUIRED_CONFIGS = ('channel_state_init_file', 'dt_routing', 'dt_total')

    # Network dependent matrices and vectors from routing parameters file
    river_ids: IntArray  # n x 1 - river ID for each segment
    downstream_river_ids: IntArray  # n x 1 - downstream river ID for each segment, -1 if no downstream segment
    compute_groups: IntArray  # n x 1 - segments with same ID can be computed concurrently (e.g. topological levels)
    k: FloatArray  # n x 1 - K values for each segment  # todo possibly can be deleted later
    x: FloatArray  # n x 1 - X values for each segment  # todo possibly can be deleted later
    c1: FloatArray  # n x 1 - C1 values for each segment => f(k, x, dt_routing)
    c2: FloatArray  # n x 1 - C2 values for each segment => f(k, x, dt_routing)
    c3: FloatArray  # n x 1 - C3 values for each segment => f(k, x, dt_routing)
    c4: FloatArray  # n x 1 - C4 values for each segment => f(k, x, dt_routing) - used if lateral inflow provided

    # derived indexing
    downstream_indices: IntArray
    downstream_c1: FloatArray
    downstream_c2: FloatArray
    upstream_indptr: IntArray
    upstream_indices: IntArray

    # State variables
    channel_state: FloatArray  # routing depends only on a channel state vector
    _network_time_signature: tuple[Any, ...] | None = None  # check if time params change between computes

    # Time options
    dt_routing: int  # compute time step, must be divisible into and <= min(dt_runoff, dt_discharge)
    dt_runoff: int  # time between lateral inflows, must be 1) constant, divisible into and <= dt_total
    dt_discharge: int  # time to average routed flows and save them, must be <=
    dt_total: int  # how long to simulate,

    # methods that are overridable via dependency injection
    _write_discharges: WriteDischargesFn

    def __init__(self, configs: PathInput | Configs | None = None, **kwargs: Any) -> None:
        # parse and create configs
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
            self.channel_state = np.zeros(self.river_ids.shape[0], dtype=np.float32)
            return
        self.logger.debug('Reading Initial State from Parquet')
        self.channel_state = pd.read_parquet(state_file).values.flatten().astype(np.float32, copy=False)
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
        """
        Creates several arrays derived from the routing parameters file used for computations

        Uses the columns:
            - river_ids: n x 1 - river ID for each segment (directly from params)
            - downstream_river_ids: n x 1 - downstream river ID for each segment (directly from params)
            - k: n x 1 - K values for each segment (directly from params)
            - x: n x 1 - X values for each segment (directly from params)
        To calculate the columns:
            - downstream_indices: n x 1 - index of downstream river, -1 if no downstream river
            - upstream_indptr: n + 1 x 1 - index pointer for start of upstream segments in upstream_indices
            - upstream_indices: m x 1 - indices of upstream rivers, m is the total upstream connections in all rivers
            - compute_groups: n x 1 - segments with same ID can be computed concurrently (e.g. topological levels)
        """
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

        self.river_ids = np.ascontiguousarray(df[self.cfg.var_river_id].to_numpy(copy=False), dtype=np.int64)
        downstream_river_ids = np.ascontiguousarray(df['downstream_river_id'].to_numpy(copy=False), dtype=np.int64)
        self.k = np.ascontiguousarray(df['k'].to_numpy(copy=False), dtype=np.float32)
        self.x = np.ascontiguousarray(df['x'].to_numpy(copy=False), dtype=np.float32)

        n = self.river_ids.shape[0]
        river_id_set = set(self.river_ids.tolist())
        downstream_ids = {d for d in downstream_river_ids.tolist() if d > 0}
        unknown_downstream_ids = sorted(downstream_ids - river_id_set)
        if unknown_downstream_ids:
            raise ValueError(f'params_file has downstream IDs not in river_id column: {unknown_downstream_ids}')

        # to avoid constant lookups and to allow carefully sorting arrays for better cpu-memory access patterns, make
        # a 1D array giving the index of the downstream river in parameter arrays, -1 if none downstream
        river_index = {int(river_id): idx for idx, river_id in enumerate(self.river_ids.tolist())}
        self.downstream_indices = np.full(n, -1, dtype=np.int32)
        counts = np.zeros(n, dtype=np.int32)
        for upstream_idx, downstream_river_id in enumerate(downstream_river_ids.tolist()):
            if downstream_river_id < 0:
                continue
            downstream_idx = river_index[int(downstream_river_id)]
            if downstream_idx <= upstream_idx:
                raise ValueError('params_file must be topologically sorted upstream to downstream')
            self.downstream_indices[upstream_idx] = downstream_idx
            counts[downstream_idx] += 1

        self.upstream_indptr = np.empty(n + 1, dtype=np.int32)
        self.upstream_indptr[0] = 0
        np.cumsum(counts, out=self.upstream_indptr[1:])

        self.upstream_indices = np.empty(int(self.upstream_indptr[-1]), dtype=np.int32)
        write_pos = self.upstream_indptr[:-1].copy()
        for upstream_idx, downstream_river_id in enumerate(downstream_river_ids.tolist()):
            if downstream_river_id < 0:
                continue
            downstream_idx = river_index[int(downstream_river_id)]
            pos = write_pos[downstream_idx]
            self.upstream_indices[pos] = upstream_idx
            write_pos[downstream_idx] += 1

        self.logger.log(PROGRESS, f'Network: {self.river_ids.shape[0]} river segments')
        return

    def _set_muskingum_coefficients(self, dt_routing: float) -> None:
        self.logger.debug('Calculating Muskingum coefficients')
        dt_div_k = dt_routing / self.k
        denominator = dt_div_k + (2 * (1 - self.x))
        _2x = 2 * self.x
        # when arrays are contiguous, they iterate much faster in numba kernels which sequentially iterate
        self.c1 = np.ascontiguousarray((dt_div_k - _2x) / denominator, dtype=np.float32)
        self.c2 = np.ascontiguousarray((dt_div_k + _2x) / denominator, dtype=np.float32)
        self.c3 = np.ascontiguousarray(((2 * (1 - self.x)) - dt_div_k) / denominator, dtype=np.float32)
        self.c4 = np.ascontiguousarray(self.c1 + self.c2, dtype=np.float32)
        if not np.allclose(self.c1 + self.c2 + self.c3, 1):
            self.logger.warning('Muskingum coefficients do not sum to 1')
            self.logger.debug(f'c1: {self.c1}')
            self.logger.debug(f'c2: {self.c2}')
            self.logger.debug(f'c3: {self.c3}')
            raise ValueError('Muskingum coefficients do not sum to 1, check routing parameters and time step')

        # todo: figure out how to check for invalid dt and parameters without needing to clamp later.
        # # Courant-like check: dt_routing >= 2*K*X guarantees C1 >= 0 which means negative discharge is impossible
        # n_violations = int(np.sum(self.c1 < 0))
        # if n_violations:
        #     worst = np.min(self.c1)
        #     self.logger.warning(
        #         f'Courant check: C1 < 0 in {n_violations} river segments, worst violation is {worst:.2f}.'
        #     )
        # shuffling arrays to list coefficient of the downstream increases performance of kernel which can
        # read sequentially when solving each river, rather than essentially randomly throughout the array
        self.downstream_c1 = np.zeros(self.downstream_indices.shape[0], dtype=np.float32)
        self.downstream_c2 = np.zeros(self.downstream_indices.shape[0], dtype=np.float32)
        valid = self.downstream_indices >= 0
        self.downstream_c1[valid] = self.c1[self.downstream_indices[valid]]
        self.downstream_c2[valid] = self.c2[self.downstream_indices[valid]]
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
        discharge_array = discharge_array.astype(np.float32, copy=False)
        self._write_discharges(dates, discharge_array, self.cfg.discharge_files[0])
        self.logger.info('-' * 60)
        return

    def _router(self, num_output_steps: int, num_routing_per_output: int) -> FloatArray:
        """Route discharge without lateral inflow"""
        self.logger.debug('Getting initial state arrays')
        q_init = self.channel_state
        if not np.any(q_init):
            self.logger.warning(
                'Initial channel state is all zeros. Muskingum routing without lateral inflow requires a '
                'non-zero initial state to produce meaningful results. Provide channel_state_init_file.'
            )

        n = self.river_ids.shape[0]
        discharge_array = np.zeros((num_output_steps, n), dtype=np.float32)
        q_t = q_init.astype(np.float32, copy=True)

        static_muskingum(
            q_t=q_t,
            discharge_array=discharge_array,
            downstream_indices=self.downstream_indices,
            downstream_c1=self.downstream_c1,
            downstream_c2=self.downstream_c2,
            c3=self.c3,
            n_rivers=n,
            n_steps=num_output_steps,
            n_substeps=num_routing_per_output
        )

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
