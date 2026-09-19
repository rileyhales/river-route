import datetime
import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Self

import numpy as np
import pandas as pd
from tqdm import tqdm

from .._logging import PROGRESS, build_logger
from ..configs import Configs
from ..network import Network
from ..runoff import CellRunoff, Runoff, RunoffGaussianGrid, RunoffVlateral
from ..types import DatetimeArray, FloatArray, IntArray, WriteDischargesFn
from ._kernel_registry import dispatch, dispatch_grid
from .writers import netcdf_writer

__all__ = ['Router']


class Router:
    """
    Muskingum style river routing allowing
    - static or dynamic coefficients (e.g. muskingum vs muskingum-cunge style)
    - channel-only or volumetric lateral water inputs (forcing)
    - uniform or unit-hydrograph runoff transformation (transform)
    - standard or stability-expanded networks (using reach subdivisions and substeps)

    A Router owns one simulation: the coefficients derived from the network parameters and the time options, the
    channel state, and the routing loop. The network itself -- the ids, the topology, the k and x vectors, the
    concurrent partition, and the stability analysis over them -- belongs to ``Network``, which the Router holds
    and reads through. That split is what makes the network work reusable: a Network parses and partitions a
    parameter table once, and every simulation over it reuses the result.
    """

    configs: Configs
    logger: logging.Logger
    network: Network  # ids, topology, k, x, the partition, and the stability analysis
    runoff: Runoff | None  # the weight table gridded runoff from the configs is aggregated with

    # calculated muskingum coefficients - consumed by static kernels or modified by dynamic kernels
    c1: FloatArray  # n x 1 - C1 values for each segment => f(k, x, dt_routing)
    c2: FloatArray  # n x 1 - C2 values for each segment => f(k, x, dt_routing)
    c3: FloatArray  # n x 1 - C3 values for each segment => f(k, x, dt_routing)
    c4: FloatArray  # n x 1 - C4 values for each segment => f(k, x, dt_routing) - used if lateral inflow provided
    c4_dt: FloatArray  # n x 1 - c4 / dt_runoff, scales lateral inflow volumes to a rate
    downstream_c1: FloatArray  # n x 1 - c1 of the downstream river, -1 positions left at zero
    downstream_c2: FloatArray  # n x 1 - c2 of the downstream river, -1 positions left at zero

    # The schedule the kernels sweep, bound from the Network per route(). These are index ranges into the
    # parameter table as given: no river vector is ever reordered and no forcing or discharge array is
    # ever gathered.
    routing_jobs: tuple[tuple[IntArray, IntArray, int, int], ...]  # (block_starts, block_stops, outlet, region)
    cut_target: IntArray  # (n_regions,) river each region's outlet drains into, -1 at a basin outlet
    threads: int = 1  # the threads passed to the most recent route(), read by writers such as zarr_writer

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
    _discharge_writer: WriteDischargesFn

    def __init__(self, configs: Configs, network: Network | None = None, runoff: Runoff | None = None) -> None:
        """
        Args:
            configs: options describing the simulation, from ``Configs(...)`` or ``Configs.from_file(path)``. They
                are validated for routing when ``route`` is called. A Router takes its options only from this
                object, so change them with ``Configs.replace`` and build a Router from the result.
            network: the network to route over. One is built from the params file the first time it is needed
                when none is given, so pass one to reuse a parsed and partitioned network across Routers.
            runoff: the RunoffGaussianGrid gridded runoff is aggregated with. One is built from the weight table the
                first time it is needed when none is given, so pass one to reuse a read weight table, or to route
                runoff a custom RunoffGaussianGrid prepares.
        """
        if not isinstance(configs, Configs):
            raise TypeError('Router must be given an rr.Configs object')
        if not isinstance(network, Network) and network is not None:
            raise TypeError(f'network must be a Network or None, got {type(network).__name__}')
        if not isinstance(runoff, RunoffGaussianGrid) and runoff is not None:
            raise TypeError(f'runoff must be a RunoffGaussianGrid or None, got {type(runoff).__name__}')

        self.configs = configs
        self.network = network if network is not None else Network.from_configs(configs)
        self.runoff = runoff

        # configure logging - progress bar and info/debug logs are mutually exclusive
        self.logger = build_logger(self.configs, 'router')
        self.logger.debug('Logger initialized')

        # default discharge writer; overridable via set_discharge_writer
        self._discharge_writer = netcdf_writer
        return

    def __repr__(self) -> str:
        return f'{type(self).__name__}(params_file={self.configs.params_file!r})'

    ################################################
    # Initial and final state handling
    ################################################

    def _read_initial_state(self) -> None:
        """Read the initial channel state from the config. Called on every route() so that repeated calls on
        the same object always start from the configured state instead of the previous run's final state."""
        n_rivers = self.network.river_ids.shape[0]
        state_file = self.configs.channel_state_init_file
        if not state_file:
            self.logger.warning('channel_state_init_file not provided. Defaulting to zero initial conditions')
            self.channel_state = np.zeros(n_rivers, dtype=np.float32)
            return
        self.logger.debug('Reading Initial State from Parquet')
        state = pd.read_parquet(state_file).values.flatten().astype(np.float32, copy=False)
        if state.shape[0] != n_rivers:
            raise ValueError(
                f'channel_state_init_file has {state.shape[0]} values but {self.configs.params_file} has '
                f'{n_rivers} rivers. The state file must have one row per river in the same order.'
            )
        self.channel_state = state
        return

    def _write_final_state(self) -> None:
        final_state_file = self.configs.channel_state_final_file
        if not final_state_file:
            return
        self.logger.debug('Writing Final State to Parquet')
        pd.DataFrame({'Q': self.channel_state}).to_parquet(self.configs.channel_state_final_file)
        return

    ################################################
    # Prepare arrays for routing
    ################################################

    def _set_routing_schedule(self, thread_pool: ThreadPoolExecutor | None = None, threads: int = 1) -> None:
        """Bind the index ranges the kernels sweep, derived and cached by the Network so that repeated
        simulations over one network never re-partition it."""
        self.routing_jobs, self.cut_target = self.network.routing_schedule(
            threads=threads, concurrent=thread_pool is not None
        )
        return

    def _set_channel_time_options(self) -> None:
        self.dt_routing = self.configs.dt_routing
        self.dt_total = self.configs.dt_total
        self.dt_discharge = self.configs.dt_discharge or self.dt_routing
        self.dt_runoff = self.dt_discharge  # Assign to pass time validation. It is never used in channel routing
        self._validate_time_options()
        return

    def _set_forced_time_options(self, dates: DatetimeArray) -> None:
        self.dt_runoff = self.configs.dt_runoff or (dates[1] - dates[0]).astype('timedelta64[s]').astype(int)
        self.dt_discharge = self.configs.dt_discharge or self.dt_runoff
        self.dt_total = self.configs.dt_total or self.dt_runoff * dates.shape[0]
        if not self.configs.dt_routing:
            self.logger.warning('dt_routing was not provided or is Null/False, defaulting to dt_runoff')
        self.dt_routing = self.configs.dt_routing or self.dt_runoff
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
        dt_div_k = self.dt_routing / self.network.k
        denominator = dt_div_k + (2 * (1 - self.network.x))
        _2x = 2 * self.network.x
        # contiguous arrays iterate faster in kernels due to cpu and ram access patterns
        self.c1 = np.ascontiguousarray((dt_div_k - _2x) / denominator, dtype=np.float32)
        self.c2 = np.ascontiguousarray((dt_div_k + _2x) / denominator, dtype=np.float32)
        self.c3 = np.ascontiguousarray(((2 * (1 - self.network.x)) - dt_div_k) / denominator, dtype=np.float32)
        self.c4 = np.ascontiguousarray(self.c1 + self.c2, dtype=np.float32)
        self.c4_dt = np.ascontiguousarray(self.c4 / self.dt_runoff, dtype=np.float32)
        self.network.check_stability(self.dt_routing, action=self.configs.unstable_coefficients)
        if not np.allclose(self.c1 + self.c2 + self.c3, 1):
            self.logger.warning('Muskingum coefficients do not sum to 1')
            self.logger.debug(f'c1: {self.c1}')
            self.logger.debug(f'c2: {self.c2}')
            self.logger.debug(f'c3: {self.c3}')
            raise ValueError('Muskingum coefficients do not sum to 1, check routing parameters and time step')

        # shuffling arrays to list coefficient of the downstream increases performance of kernel which can
        # read sequentially when solving each river, rather than essentially randomly throughout the array
        self.downstream_c1 = np.zeros(self.network.downstream_indices.shape[0], dtype=np.float32)
        self.downstream_c2 = np.zeros(self.network.downstream_indices.shape[0], dtype=np.float32)
        valid = self.network.downstream_indices >= 0
        self.downstream_c1[valid] = self.c1[self.network.downstream_indices[valid]]
        self.downstream_c2[valid] = self.c2[self.network.downstream_indices[valid]]
        return

    ################################################
    # Methods to execute routing simulation
    ################################################

    def route(self, thread_pool: ThreadPoolExecutor | None = None, threads: int = 1) -> Self:
        """
        Execute the simulation described by the configs and routing parameters. The configs are validated for routing
        first. The thread pool and thread count are runtime resources, not configs, so they are given here.

        Args:
            thread_pool: optional pool to route and aggregate gridded runoff concurrently on. It is used as given and
                never shut down here, so it can be shared and closed by the caller's with block. Without one,
                routing is single-threaded regardless of ``threads``.
            threads: number of regions the network is split into, and river ranges gridded runoff is aggregated in,
                when ``thread_pool`` is given. Also the concurrency limit of writers that follow it, like zarr_writer.

        Returns:
            Self: the class instance with updated channel_state and output files written to disk
        """
        if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
            raise ValueError(f'threads must be an integer >= 1, got {threads!r}')
        self.threads = threads
        # start timer
        self.logger.log(PROGRESS, 'Beginning routing')
        t1 = datetime.datetime.now()
        self.logger.debug('Validating configs')
        self.configs.validate_routing()
        self.logger.debug(self)
        # the Network parses and partitions the parameter table; both are cached there and reused across runs
        self._set_routing_schedule(thread_pool, threads)  # index ranges the kernels sweep; no vector is reordered
        # read state, route, write state
        self._read_initial_state()
        self._execute_routing(thread_pool, threads)
        self._write_final_state()
        # log total time
        t2 = datetime.datetime.now()
        self.logger.log(PROGRESS, f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return self

    def _execute_routing(self, thread_pool: ThreadPoolExecutor | None, threads: int) -> None:
        # there are two types of loops, one for channel only, one if lateral forcing(s) are provided.
        if self.configs.forcing == 'channel':
            self._execute_routing_channel(thread_pool)
        else:
            self._execute_routing_forced(thread_pool, threads)
        return

    def _execute_routing_channel(self, thread_pool: ThreadPoolExecutor | None) -> None:
        self.logger.info('-' * 60)

        # channel routing time options and (static) coefficients
        self._set_channel_time_options()
        self._set_static_muskingum_coefficients()

        self.logger.debug('Starting routing computation')
        q_t = self.channel_state.astype(np.float32, copy=True)
        discharge_array = self._discharge_buffer()
        dispatch(self, q_t=q_t, discharge_array=discharge_array, thread_pool=thread_pool)
        self.channel_state = q_t

        # generate dates since they cannot be copied from an external forcing file
        dates = pd.date_range(
            start=self.configs.start_datetime,
            periods=self.num_runoff_steps,
            freq=pd.to_timedelta(self.dt_discharge, unit='s'),
        ).to_numpy()

        # write outputs
        self.logger.debug('Writing Discharge Array to File')
        self._discharge_writer(self, dates, discharge_array, self.configs.discharge_files[0])
        self.logger.info('-' * 60)
        return

    def _execute_routing_forced(self, thread_pool: ThreadPoolExecutor | None, threads: int) -> None:
        self._ensemble_member_states = []

        total_files = len(self.configs.vlateral_files or self.configs.grid_runoff_files or [])
        if self.configs.vlateral_files:
            runoff_iter = RunoffVlateral().reader(self.configs.vlateral_files)
        else:
            if self.runoff is None:
                self.runoff = RunoffGaussianGrid.from_configs(self.configs)
            if self._fuses_grid_runoff():
                self.logger.debug('Aggregating and routing gridded runoff in one pass with the fused kernel')
                runoff_iter = self.runoff.cell_reader(self.configs.grid_runoff_files)
            else:
                runoff_iter = self.runoff.reader(self.configs.grid_runoff_files, thread_pool, threads)
        file_iter = (
            (dates, forcing, runoff_file, discharge_file)
            for (dates, forcing, runoff_file), discharge_file in zip(
                runoff_iter, self.configs.discharge_files, strict=True
            )
        )
        if self.configs.progress_bar:
            file_iter = tqdm(file_iter, total=total_files, desc='Files Routed')

        coeff_dt: tuple[int, int] | None = None  # (dt_routing, dt_runoff) the static coefficients were built for
        for dates, forcing, runoff_file, discharge_file in file_iter:
            self.logger.info('-' * 60)
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
                if not isinstance(forcing, CellRunoff):  # the fused kernel reads only the steps it routes
                    forcing = forcing[: self.num_runoff_steps]
            # static coefficients depend only on dt_routing/dt_runoff, so rebuild them only when those change
            # across files (dynamic coefficients are rebuilt inside the kernel each substep)
            if self.configs.coeff != 'dynamic' and (self.dt_routing, self.dt_runoff) != coeff_dt:
                self._set_static_muskingum_coefficients()
                coeff_dt = (self.dt_routing, self.dt_runoff)
            self.logger.debug('Starting routing computation')
            q_t = self.channel_state.astype(np.float32, copy=True)
            q_array = self._discharge_buffer()
            if isinstance(forcing, CellRunoff):
                dispatch_grid(self, q_t=q_t, discharge_array=q_array, runoff=forcing, thread_pool=thread_pool)
            else:
                dispatch(
                    self,
                    q_t=q_t,
                    discharge_array=q_array,
                    vlateral=forcing.astype(np.float32, copy=False),  # the kernel runner picks the layout
                    thread_pool=thread_pool,
                )
            if self.configs.runoff_processing_mode == 'sequential':
                self.logger.debug('Updating Channel State for Next Sequential Computation')
                self.channel_state = q_t
            elif self.configs.runoff_processing_mode == 'ensemble':
                self.logger.debug('Recording Member State for Final State Aggregation')
                self._ensemble_member_states.append(q_t.copy())

            if self.dt_discharge > self.dt_runoff:
                self.logger.debug('Resampling dates and discharges to specified timestep')
                q_array = q_array.reshape(
                    (
                        int(self.dt_total / self.dt_discharge),
                        int(self.dt_discharge / self.dt_runoff),
                        self.network.river_ids.shape[0],
                    )
                ).mean(axis=1)
                dates = dates[:: self.num_runoff_steps_per_discharge]

            self.logger.debug('Writing Discharge Array to File')
            self._discharge_writer(self, dates, q_array, discharge_file, runoff_file)

        if self.configs.runoff_processing_mode == 'ensemble':
            self.channel_state = np.array(self._ensemble_member_states).mean(axis=0)
        self.logger.info('-' * 60)
        return

    def _discharge_buffer(self) -> FloatArray:
        """
        A zeroed (time, river) array for the routed discharge. River order kernels write each river's series in place
        when the array is the transpose of a C-order (river, time) buffer, which skips transposing into (time, river),
        so that layout is used whenever the discharge writer declares ``discharge_layout = 'river'`` (it reads that
        layout at least as fast). Otherwise, and always for time order, the array is C-order (time, river).
        """
        shape = (self.num_runoff_steps, self.network.river_ids.shape[0])
        river_layout = getattr(self._discharge_writer, 'discharge_layout', 'time') == 'river'
        if self.configs.routing_order == 'river' and river_layout:
            return np.zeros(shape[::-1], dtype=np.float32).T
        return np.zeros(shape, dtype=np.float32)

    def _fuses_grid_runoff(self) -> bool:
        """
        Whether gridded runoff is aggregated and routed in one pass by the fused river order kernel instead of being
        aggregated to a vlateral array first. The fused kernel implements river order routing with static
        coefficients on a standard network with uniform lateral forcing, so every other combination aggregates first.
        """
        return (
            self.configs.routing_order == 'river'
            and self.configs.coeff == 'static'
            and self.configs.transform == 'uniform'
            and self.configs.network == 'standard'
        )

    ################################################
    # Dependency injection methods for users to overwrite default behaviors without subclassing
    ################################################

    def set_discharge_writer(self, func: WriteDischargesFn) -> Self:
        """
        Replace the function which writes computed flows to disc. Defaults in rr.router.writers.

        Args:
            func (callable): function that takes router, dates, discharge_array, discharge_file, runoff_file and
                returns None. It is called once per routed input file with this Router as the first argument.
        """
        self._discharge_writer = func
        return self
