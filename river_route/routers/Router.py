import datetime
import itertools
import json
import logging
import sys
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Self

import numpy as np
import pandas as pd
import xarray as xr
import yaml
from tqdm import tqdm

from .. import streams, writers
from ..runoff import Runoff
from ..types import DatetimeArray, FloatArray, IntArray, PathInput, VlateralGeneratorSignature, WriteDischargesFn
from ._kernel_registry import dispatch
from .Configs import Configs

__all__ = ['Router']

PROGRESS = 25
logging.addLevelName(PROGRESS, 'PROGRESS')

_ROUTER_COUNT = itertools.count()  # id() is reused after garbage collection, so it cannot name a logger uniquely


class Router:
    """
    Muskingum style river routing allowing
    - static or dynamic coefficients (e.g. muskingum vs muskingum-cunge style)
    - channel-only or volumetric lateral water inputs (forcing)
    - uniform or unit-hydrograph runoff transformation (transform)
    - standard or stability-expanded networks (using reach subdivisions and substeps)
    """

    cfg: Configs
    logger: logging.Logger

    # Network Parameters
    river_ids: IntArray  # n x 1 - river ID for each segment
    next_river_ids: IntArray  # n x 1 - downstream river ID for each segment, -1 if no downstream
    # routing parameters
    k: FloatArray  # n x 1 - K values for each segment
    x: FloatArray  # n x 1 - X values for each segment
    alpha: FloatArray  # n x 1 - alpha values for each segment   (dynamic coeff only)
    beta: FloatArray  # n x 1 - beta values for each segment     (dynamic coeff only)
    # calculated muskingum coefficients - consumed by static kernels or modified by dynamic kernels
    c1: FloatArray  # n x 1 - C1 values for each segment => f(k, x, dt_routing)
    c2: FloatArray  # n x 1 - C2 values for each segment => f(k, x, dt_routing)
    c3: FloatArray  # n x 1 - C3 values for each segment => f(k, x, dt_routing)
    c4: FloatArray  # n x 1 - C4 values for each segment => f(k, x, dt_routing) - used if lateral inflow provided
    c4_dt: FloatArray  # n x 1 - c4 / dt_runoff, scales lateral inflow volumes to a rate

    # derived indexing
    downstream_indices: IntArray
    downstream_c1: FloatArray
    downstream_c2: FloatArray

    # Region schedule for threaded routing (see streams.assign_regions). Built once per Router and then only
    # read, never rebuilt per input file. These are index ranges into the parameter table as given: no river
    # vector is ever reordered and no forcing or discharge array is ever gathered.
    routing_jobs: tuple[tuple[IntArray, IntArray, int, int], ...]  # (block_starts, block_stops, outlet, region)
    cut_target: IntArray  # (n_regions,) river each region's outlet drains into, -1 at a basin outlet
    _stored_region: IntArray | None  # the params file 'region' column, when it has one

    # State variables
    channel_state: FloatArray  # routing depends only on a channel state vector
    _ensemble_member_states: list[FloatArray]  # for ensemble routing

    # Time options
    dt_routing: int  # compute time step, must be divisible into and <= min(dt_runoff, dt_discharge)
    dt_runoff: int  # time between lateral inflows, must be 1) constant, divisible into and <= dt_total
    dt_discharge: int  # time to average routed flows and save them, must be <=
    dt_total: int  # how long to simulate,
    num_runoff_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_discharge: int

    # methods overridable via dependency injection
    _write_discharges: WriteDischargesFn

    def __init__(self, configs: PathInput | None = None, **kwargs: Any) -> None:
        # parse and create configs
        raw: dict[str, Any] = {}
        if configs is not None and configs != '':
            configs = str(configs)
            if configs.endswith('.json'):
                with open(configs) as f:
                    raw = json.load(f)
            elif configs.endswith(('.yml', '.yaml')):
                with open(configs) as f:
                    raw = yaml.load(f, Loader=yaml.FullLoader)
            else:
                raise RuntimeError('Unrecognized simulation config file type. Must be .json or .yaml')
        raw.update(kwargs)
        self.cfg = Configs.from_mapping(raw)

        # configure logging - progress bar and info/debug logs are mutually exclusive
        self.logger = logging.getLogger(f'river_route.router{next(_ROUTER_COUNT)}')
        self.logger.propagate = False  # this logger owns its handler; propagating would double every message
        self.logger.disabled = not self.cfg.log
        self.logger.setLevel(self.cfg.log_level)
        if self.cfg.log_stream == 'stdout':
            self.logger.addHandler(logging.StreamHandler(sys.stdout))
        else:
            self.logger.addHandler(logging.FileHandler(self.cfg.log_stream))
        self.logger.handlers[0].setFormatter(logging.Formatter(self.cfg.log_format))
        self.logger.debug('Logger initialized')

        # default discharge writer; overridable via set_write_discharges
        self._write_discharges = writers.netcdf_writer
        return

    def __repr__(self) -> str:
        return f'{type(self).__name__}(params_file={self.cfg.params_file!r})'

    ################################################
    # Initial and final state handling
    ################################################

    def _read_initial_state(self) -> None:
        """Read the initial channel state from the config. Called on every route() so that repeated calls on
        the same object always start from the configured state instead of the previous run's final state."""
        n_rivers = self.river_ids.shape[0]
        state_file = self.cfg.channel_state_init_file
        if not state_file:
            self.logger.warning('channel_state_init_file not provided. Defaulting to zero initial conditions')
            self.channel_state = np.zeros(n_rivers, dtype=np.float32)
            return
        self.logger.debug('Reading Initial State from Parquet')
        state = pd.read_parquet(state_file).values.flatten().astype(np.float32, copy=False)
        if state.shape[0] != n_rivers:
            raise ValueError(
                f'channel_state_init_file has {state.shape[0]} values but {self.cfg.params_file} has '
                f'{n_rivers} rivers. The state file must have one row per river in the same order.'
            )
        self.channel_state = state
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

    def _set_vectors_from_params(self):
        self.logger.debug('Calculating network dependent vectors')
        # todo use a cache for the params file to speed up possible repeated reads
        df = pd.read_parquet(self.cfg.params_file)

        required = [self.cfg.var_river_id, 'next_river_id', 'k', 'x']
        if self.cfg.coeff == 'dynamic':
            required += ['alpha', 'beta']
        missing = [column for column in required if column not in df.columns]
        if missing:
            raise ValueError(f'params_file is missing required column(s): {", ".join(missing)}')

        if df[self.cfg.var_river_id].duplicated().any():
            raise ValueError('params_file contains duplicate river IDs.')

        # a stored partition (streams.partition_network) is reused as-is; absent, one is derived at setup
        self._stored_region = (
            np.ascontiguousarray(df['region'].to_numpy(copy=False), dtype=np.int64) if 'region' in df.columns else None
        )

        self.river_ids = np.ascontiguousarray(df[self.cfg.var_river_id].to_numpy(copy=False), dtype=np.int64)
        self.next_river_ids = np.ascontiguousarray(df['next_river_id'].to_numpy(copy=False), dtype=np.int64)
        self.k = np.ascontiguousarray(df['k'].to_numpy(copy=False), dtype=np.float32)
        self.x = np.ascontiguousarray(df['x'].to_numpy(copy=False), dtype=np.float32)

        # dynamic coefficients additionally need the alpha and beta columns
        if self.cfg.coeff == 'dynamic':
            self.alpha = np.ascontiguousarray(df['alpha'].to_numpy(copy=False), dtype=np.float32)
            self.beta = np.ascontiguousarray(df['beta'].to_numpy(copy=False), dtype=np.float32)
            if np.any(self.alpha <= 0):
                raise ValueError('alpha column in params_file must be strictly positive')
        return

    def _set_connectivity_vectors(self):
        """Build the index vectors describing network connectivity from river_ids and next_river_ids"""
        self.logger.debug('Calculating network connectivity vectors')
        n = self.river_ids.shape[0]
        river_index = {int(river_id): idx for idx, river_id in enumerate(self.river_ids.tolist())}

        # 1D array giving the index of the downstream river in the parameter arrays, -1 if none downstream
        self.downstream_indices = np.full(n, -1, dtype=np.int32)
        for upstream_idx, next_river_id in enumerate(self.next_river_ids.tolist()):
            if next_river_id < 0:
                continue
            downstream_idx = river_index.get(int(next_river_id))
            if downstream_idx is None:
                raise ValueError(f'params_file next_river_id {next_river_id} is not in the river_id column')
            if downstream_idx <= upstream_idx:
                raise ValueError('params_file must be topologically sorted upstream to downstream')
            self.downstream_indices[upstream_idx] = downstream_idx

        self.logger.log(logging.INFO, f'Network: {n} river segments')
        return

    def _set_region_schedule(self, thread_pool: ThreadPoolExecutor | None = None) -> None:
        """
        Build the list of index ranges the kernels sweep. Computed once per Router and reused for every file.

        The parameter table is always used in the order it is given. Nothing in the routing path reorders a
        river, so the forcing, the state and the routed discharge stay in parameter file order from end to end.
        Threaded routing additionally requires that order to be DFS computation order -- a river following all
        of its own upstream rivers, which makes every subtree a contiguous block that a worker can be handed as
        a plain index range. That is checked against the input and reported if absent, never corrected here.

        The jobs are ordered longest first so the pool packs the very uneven region sizes a river network
        produces, with the main stem last. That orders the work queue only; the rivers inside each block keep
        their file positions. Single-threaded routing, without a thread_pool or with threads=1, gets one job spanning
        the whole network, which is the same sweep the kernels did before any of this existed.
        """
        n = self.river_ids.shape[0]
        whole = (np.array([0], dtype=np.int32), np.array([n], dtype=np.int32), -1, 0)
        if thread_pool is None or self.cfg.threads < 2:
            self.routing_jobs = (whole,)
            self.cut_target = np.zeros(0, dtype=np.int32)
            return

        downstream_index = self.downstream_indices.astype(np.int64)
        if self._stored_region is not None:
            self.logger.debug('Using the region column stored in the params file')
            region = self._stored_region
        else:
            self.logger.debug('Deriving a network partition for threaded routing')
            region, _ = streams.assign_regions(downstream_index, threads=self.cfg.threads)
        layout = streams.regions_to_layout(region, downstream_index)

        n_regions = layout['n_regions']
        if not n_regions:
            self.logger.warning('The network did not split into any concurrent regions; routing single-threaded')
            self.routing_jobs = (whole,)
            self.cut_target = np.zeros(0, dtype=np.int32)
            return

        starts, stops = layout['region_starts'], layout['region_stops']
        order = np.argsort(starts - stops)  # submission order for the pool only; river order is untouched
        jobs = [(starts[r : r + 1], stops[r : r + 1], int(layout['region_outlet'][r]), int(r)) for r in order.tolist()]
        jobs.append((layout['stem_starts'], layout['stem_stops'], -1, 0))
        self.routing_jobs = tuple(jobs)
        self.cut_target = layout['cut_target']

        main_stem = int((layout['stem_stops'] - layout['stem_starts']).sum())
        self.logger.log(
            logging.INFO,
            f'Partition: {n_regions} concurrent regions on {self.cfg.threads} threads, '
            f'{main_stem} rivers ({main_stem / n:.2%}) routed sequentially as the main stem',
        )
        return

    def _set_channel_time_options(self) -> None:
        self.dt_routing = self.cfg.dt_routing
        self.dt_total = self.cfg.dt_total
        self.dt_discharge = self.cfg.dt_discharge or self.dt_routing
        self.dt_runoff = self.dt_discharge  # Assign to pass time validation. It is never used in channel routing
        self._validate_time_options()
        return

    def _set_forced_time_options(self, dates: DatetimeArray) -> None:
        self.dt_runoff = self.cfg.dt_runoff or (dates[1] - dates[0]).astype('timedelta64[s]').astype(int)
        self.dt_discharge = self.cfg.dt_discharge or self.dt_runoff
        self.dt_total = self.cfg.dt_total or self.dt_runoff * dates.shape[0]
        if not self.cfg.dt_routing:
            self.logger.warning('dt_routing was not provided or is Null/False, defaulting to dt_runoff')
        self.dt_routing = self.cfg.dt_routing or self.dt_runoff
        self._validate_time_options()
        return

    def _validate_time_options(self):
        # check that time options have the correct relative sizes
        if self.dt_total < self.dt_runoff:
            raise ValueError('dt_total must be >= dt_runoff')
        if self.dt_total < self.dt_discharge:
            raise ValueError('dt_total must be >= dt_discharge')
        if self.dt_discharge < self.dt_runoff:
            raise ValueError('dt_discharge must be >= dt_runoff')
        if self.dt_runoff < self.dt_routing:
            raise ValueError('dt_runoff must be >= dt_routing')
        if not (self.dt_total >= self.dt_discharge >= self.dt_runoff >= self.dt_routing):
            raise ValueError('Need dt_total >= dt_discharge >= dt_runoff >= dt_routing')
        # check that time options are evenly divisible
        if self.dt_total % self.dt_runoff != 0:
            raise ValueError('dt_total must be an integer multiple of dt_runoff')
        if self.dt_total % self.dt_discharge != 0:
            raise ValueError('dt_total must be an integer multiple of dt_discharge')
        if self.dt_discharge % self.dt_runoff != 0:
            raise ValueError('dt_discharge must be an integer multiple of dt_runoff')
        if self.dt_runoff % self.dt_routing != 0:
            raise ValueError('dt_runoff must be an integer multiple of dt_routing')

        # Now that we know time parameters are valid, set time-derived parameters for computation cycles
        self.num_runoff_steps = int(self.dt_total / self.dt_runoff)
        self.num_runoff_steps_per_discharge = int(self.dt_discharge / self.dt_runoff)  # to resample/reshape results
        # todo for synthesized networks, this will have to become a vector
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)
        return

    def _set_static_muskingum_coefficients(self):
        """
        implied dependency on having set time options and the network parameter vectors
        """
        self.logger.debug('Calculating Muskingum coefficients')
        dt_div_k = self.dt_routing / self.k
        denominator = dt_div_k + (2 * (1 - self.x))
        _2x = 2 * self.x
        # contiguous arrays iterate faster in kernels due to cpu and ram access patterns
        self.c1 = np.ascontiguousarray((dt_div_k - _2x) / denominator, dtype=np.float32)
        self.c2 = np.ascontiguousarray((dt_div_k + _2x) / denominator, dtype=np.float32)
        self.c3 = np.ascontiguousarray(((2 * (1 - self.x)) - dt_div_k) / denominator, dtype=np.float32)
        self.c4 = np.ascontiguousarray(self.c1 + self.c2, dtype=np.float32)
        self.c4_dt = np.ascontiguousarray(self.c4 / self.dt_runoff, dtype=np.float32)
        self._check_coefficient_stability()
        if not np.allclose(self.c1 + self.c2 + self.c3, 1):
            self.logger.warning('Muskingum coefficients do not sum to 1')
            self.logger.debug(f'c1: {self.c1}')
            self.logger.debug(f'c2: {self.c2}')
            self.logger.debug(f'c3: {self.c3}')
            raise ValueError('Muskingum coefficients do not sum to 1, check routing parameters and time step')

        # shuffling arrays to list coefficient of the downstream increases performance of kernel which can
        # read sequentially when solving each river, rather than essentially randomly throughout the array
        self.downstream_c1 = np.zeros(self.downstream_indices.shape[0], dtype=np.float32)
        self.downstream_c2 = np.zeros(self.downstream_indices.shape[0], dtype=np.float32)
        valid = self.downstream_indices >= 0
        self.downstream_c1[valid] = self.c1[self.downstream_indices[valid]]
        self.downstream_c2[valid] = self.c2[self.downstream_indices[valid]]
        return

    def _check_coefficient_stability(self) -> None:
        """
        Report rivers whose Muskingum coefficients are not stable for the current dt_routing.

        Non-negative coefficients require ``2*k*x <= dt_routing <= 2*k*(1-x)``. Outside that window the
        solution oscillates and the kernels clamp the resulting negative discharges to zero, which does not
        conserve mass. The c1+c2+c3 == 1 identity holds for negative coefficients too, so it cannot detect this.

        The action taken is set by the ``unstable_coefficients`` config: ``warn`` (default), ``raise``, or
        ``ignore``. See ``river_route.streams`` for the analysis tools that quantify how to fix a network.
        """
        if self.cfg.unstable_coefficients == 'ignore':
            return
        dt = self.dt_routing
        too_long = 2 * self.k * self.x > dt  # c1 < 0
        too_short = 2 * self.k * (1 - self.x) < dt  # c3 < 0
        n_unstable = int(np.count_nonzero(too_long | too_short))
        if not n_unstable:
            return
        message = (
            f'{n_unstable} of {self.k.shape[0]} rivers are not Muskingum-stable for dt_routing={dt} s '
            f'({int(np.count_nonzero(too_long))} need a larger dt_routing, '
            f'{int(np.count_nonzero(too_short))} need a smaller one). '
            f'Stability requires 2*k*x <= dt_routing <= 2*k*(1-x) for every river. '
            f'Routed discharge for these rivers oscillates and negative values are clamped to zero, '
            f'which does not conserve mass. Use river_route.streams.analyze_stability to inspect the network, '
            f"or set unstable_coefficients to 'ignore' to silence this."
        )
        if self.cfg.unstable_coefficients == 'raise':
            raise ValueError(message)
        self.logger.warning(message)
        return

    #################################################
    # Generator for lateral flow routing
    #################################################

    def _check_vlateral_alignment(self, ds: xr.Dataset, source: PathInput) -> None:
        """
        Verify that a lateral inflow dataset lines up with the routing network before its values reach the
        kernels. The kernels are compiled without bounds checking and index vlateral by river position, so a
        mismatched river count reads past the end of the array and a mismatched order silently routes each
        river's water down the wrong reach.

        Args:
            ds: dataset holding the ``vlateral`` variable with dimensions (time, river)
            source: path of the file the dataset came from, used in error messages

        Raises:
            ValueError: if the variable is missing, not 2D, ordered (river, time), has a river count that does
                not match the params file, or carries river ids that differ from the params file
        """
        rid = self.cfg.var_river_id
        if 'vlateral' not in ds:
            raise ValueError(f'{source} does not contain a vlateral variable')
        array = ds['vlateral']
        dims = tuple(array.dims)
        if len(dims) != 2:
            raise ValueError(f'{source} vlateral must have 2 dimensions (time, {rid}), found {dims}')
        if rid in dims and dims[1] != rid:
            raise ValueError(f'{source} vlateral dimensions must be ordered (time, {rid}), found {dims}')
        file_ids = ds[rid].values if rid in ds.variables else None
        self._check_river_alignment(file_ids, array.shape[1], source)
        return

    def _check_river_alignment(self, file_ids: IntArray | None, n_columns: int, source: PathInput) -> None:
        """
        Verify that the river columns of a lateral inflow source match the params file in count and order.

        Args:
            file_ids: river id of each column in order, or None when the source does not record them
            n_columns: number of river columns the source provides
            source: path of the file the columns came from, used in error messages

        Raises:
            ValueError: if the river count does not match the params file, or the river ids differ from it
        """
        rid = self.cfg.var_river_id
        n_rivers = self.river_ids.shape[0]
        if n_columns != n_rivers:
            raise ValueError(
                f'{source} provides lateral inflow for {n_columns} rivers but {self.cfg.params_file} has '
                f'{n_rivers} rivers. The two files must describe the same network.'
            )
        if file_ids is None:
            self.logger.warning(
                f'{source} has no {rid} variable so its column order cannot be verified against the params '
                f'file. Values are assumed to be in the same order as {self.cfg.params_file}.'
            )
            return
        file_ids = np.asarray(file_ids).astype(np.int64, copy=False)
        if not np.array_equal(file_ids, self.river_ids):
            n_diff = int(np.count_nonzero(file_ids != self.river_ids))
            same_set = set(file_ids.tolist()) == set(self.river_ids.tolist())
            reason = 'are in a different order than' if same_set else 'are not the same rivers as'
            raise ValueError(
                f'{source} {rid} values {reason} {self.cfg.params_file} ({n_diff} positions differ). '
                f'Sort the lateral inflow columns to match the params file river order.'
            )
        return

    def _vlateral_generator(self, thread_pool: ThreadPoolExecutor | None = None) -> VlateralGeneratorSignature:
        if self.cfg.vlateral_files:
            for lateral_file, discharge_file in zip(self.cfg.vlateral_files, self.cfg.discharge_files, strict=True):
                self.logger.info('-' * 60)
                with xr.open_dataset(lateral_file) as ds:
                    self._check_vlateral_alignment(ds, lateral_file)
                    dates = ds['time'].values.astype('datetime64[s]')
                    array = ds['vlateral'].values.astype(np.float32, copy=False)
                    yield dates, array, lateral_file, discharge_file
        elif self.cfg.grid_runoff_files and self.cfg.grid_weights_file:
            # The weight table is identical for every runoff file, so it is read and checked against the params file
            # once. Every file is then aggregated into one reused C-order buffer that the kernels read directly: no
            # per-file allocation or copy. The yielded array is overwritten by the next file.
            runoff = Runoff(
                self.cfg.grid_weights_file,
                var_river_id=self.cfg.var_river_id,
                var_runoff=self.cfg.var_grid_runoff,
                var_x=self.cfg.var_x,
                var_y=self.cfg.var_y,
                var_t=self.cfg.var_t,
                cumulative=self.cfg.grid_accumulation_type == 'cumulative',
                as_volumes=True,
                thread_pool=thread_pool,
                threads=self.cfg.threads,
            )
            self._check_river_alignment(runoff.river_ids, runoff.river_ids.shape[0], self.cfg.grid_weights_file)
            for runoff_file, discharge_file in zip(self.cfg.grid_runoff_files, self.cfg.discharge_files, strict=True):
                self.logger.info('-' * 60)
                self.logger.debug(f'Calculating vlateral: {runoff_file}')
                vlateral, dates = runoff.vlateral(runoff_file)
                yield (
                    dates.astype('datetime64[s]'),
                    vlateral.astype(np.float32, copy=False),
                    runoff_file,
                    discharge_file,
                )

    ################################################
    # Methods to execute routing simulation
    ################################################

    def route(self, thread_pool: ThreadPoolExecutor | None = None) -> Self:
        """
        Execute the simulation described by the provided configs and routing parameters. All configs, file paths,
        parameters, and options must be set when the object is initialized so that validation is performed before the
        simulation.

        Args:
            thread_pool: optional pool to route concurrently on, split into ``threads`` regions. It is used as given
                and never shut down here, so it can be shared and closed by the caller's with block. Without one,
                routing is single-threaded regardless of ``threads``.

        Returns:
            Self: the class instance with updated channel_state and output files written to disk
        """
        # start timer
        self.logger.log(PROGRESS, 'Beginning routing')
        t1 = datetime.datetime.now()
        # validate configuration options
        self.logger.debug('Validating configs')
        self.cfg.validate()
        if self.cfg.deep_validation:
            self.logger.debug('Validating contents of input files')
            self.cfg.deep_validate()
        self.logger.debug(self)
        # build network vectors
        self._set_vectors_from_params()  # river_id, next_river_id, k, x, and alpha and beta (dynamic only)
        self._set_connectivity_vectors()  # downstream_indices
        self._set_region_schedule(thread_pool)  # index ranges the kernels sweep; no vector is reordered
        # read state, route, write state
        self._read_initial_state()
        self._execute_routing(thread_pool)
        self._write_final_state()
        # log total time
        t2 = datetime.datetime.now()
        self.logger.log(PROGRESS, f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return self

    def _execute_routing(self, thread_pool: ThreadPoolExecutor | None) -> None:
        # there are two types of loops, one for channel only, one if lateral forcing(s) are provided.
        if self.cfg.forcing == 'channel':
            self._execute_routing_channel(thread_pool)
        else:
            self._execute_routing_forced(thread_pool)
        return

    def _execute_routing_channel(self, thread_pool: ThreadPoolExecutor | None) -> None:
        self.logger.info('-' * 60)

        # channel routing time options and (static) coefficients
        self._set_channel_time_options()
        self._set_static_muskingum_coefficients()

        self.logger.debug('Starting routing computation')
        q_t = self.channel_state.astype(np.float32, copy=True)
        discharge_array = np.zeros((self.num_runoff_steps, self.river_ids.shape[0]), dtype=np.float32)
        dispatch(self, q_t=q_t, discharge_array=discharge_array, thread_pool=thread_pool)
        self.channel_state = q_t

        # generate dates since they cannot be copied from an external forcing file
        dates = pd.date_range(
            start=self.cfg.start_datetime,
            periods=self.num_runoff_steps,
            freq=pd.to_timedelta(self.dt_discharge, unit='s'),
        ).to_numpy()

        # write outputs
        self.logger.debug('Writing Discharge Array to File')
        discharge_array = discharge_array.astype(np.float32, copy=False)
        self._write_discharges(self, dates, discharge_array, self.cfg.discharge_files[0])
        self.logger.info('-' * 60)
        return

    def _execute_routing_forced(self, thread_pool: ThreadPoolExecutor | None) -> None:
        self._ensemble_member_states = []

        total_files = len(self.cfg.vlateral_files or self.cfg.grid_runoff_files or [])
        file_iter = self._vlateral_generator(thread_pool)
        if self.cfg.progress_bar:
            file_iter = tqdm(file_iter, total=total_files, desc='Files Routed')

        coeff_dt: tuple[int, int] | None = None  # (dt_routing, dt_runoff) the static coefficients were built for
        for dates, vlateral, runoff_file, discharge_file in file_iter:
            self.logger.info(f'Routing vlateral: {runoff_file}')
            self._set_forced_time_options(dates)
            if self.num_runoff_steps > dates.shape[0]:
                raise ValueError(
                    f'dt_total={self.dt_total} s needs {self.num_runoff_steps} steps of lateral inflow but '
                    f'{runoff_file} provides {dates.shape[0]}. Lower dt_total or provide more input.'
                )
            if self.num_runoff_steps < dates.shape[0]:
                self.logger.debug(f'Using the first {self.num_runoff_steps} of {dates.shape[0]} steps for dt_total')
                dates = dates[: self.num_runoff_steps]
                vlateral = vlateral[: self.num_runoff_steps]
            # static coefficients depend only on dt_routing/dt_runoff, so rebuild them only when those change
            # across files (dynamic coefficients are rebuilt inside the kernel each substep)
            if self.cfg.coeff != 'dynamic' and (self.dt_routing, self.dt_runoff) != coeff_dt:
                self._set_static_muskingum_coefficients()
                coeff_dt = (self.dt_routing, self.dt_runoff)
            self.logger.debug('Starting routing computation')
            q_t = self.channel_state.astype(np.float32, copy=True)
            q_array = np.zeros((self.num_runoff_steps, self.river_ids.shape[0]), dtype=np.float32)
            dispatch(
                self,
                q_t=q_t,
                discharge_array=q_array,
                vlateral=np.ascontiguousarray(vlateral, dtype=np.float32),
                thread_pool=thread_pool,
            )
            if self.cfg.runoff_processing_mode == 'sequential':
                self.logger.debug('Updating Channel State for Next Sequential Computation')
                self.channel_state = q_t
            elif self.cfg.runoff_processing_mode == 'ensemble':
                self.logger.debug('Recording Member State for Final State Aggregation')
                self._ensemble_member_states.append(q_t.copy())

            if self.dt_discharge > self.dt_runoff:
                self.logger.debug('Resampling dates and discharges to specified timestep')
                q_array = q_array.reshape(
                    (
                        int(self.dt_total / self.dt_discharge),
                        int(self.dt_discharge / self.dt_runoff),
                        self.river_ids.shape[0],
                    )
                ).mean(axis=1)
                dates = dates[:: self.num_runoff_steps_per_discharge]

            self.logger.debug('Writing Discharge Array to File')
            q_array = q_array.astype(np.float32, copy=False)
            self._write_discharges(self, dates, q_array, discharge_file, runoff_file)

        if self.cfg.runoff_processing_mode == 'ensemble':
            self.channel_state = np.array(self._ensemble_member_states).mean(axis=0)
        self.logger.info('-' * 60)
        return

    ################################################
    # Dependency injection methods for users to overwrite default behaviors without subclassing
    ################################################

    def set_write_discharges(self, func: WriteDischargesFn) -> Self:
        """
        Replace the discharge writer, which defaults to ``river_route.writers.netcdf_writer``. Premade writers are in
        ``river_route.writers``, e.g. ``router.set_write_discharges(rr.writers.zarr_writer)``.

        Args:
            func (callable): function that takes router, dates, discharge_array, discharge_file, runoff_file and
                returns None. It is called once per routed input file with this Router as the first argument.
        """
        self._write_discharges = func
        return self
