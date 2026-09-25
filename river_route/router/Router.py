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
from ..runoff import RUNOFF_CLASS_FOR_RUNOFF_TYPE, CellRunoff, Runoff
from ..types import DatetimeArray, FloatArray, IntArray, WriteDischargesFn
from ._dispatchers import resolve_dispatcher
from .writers import zarr_writer

__all__ = ['Router']


class Router:
    """
    Muskingum style river routing allowing
    - channel-only or runoff forcing
    - static or dynamic coefficients (e.g. muskingum vs muskingum-cunge style)
    - uniform or unit-hydrograph runoff transformation (e.g. runoff transform method)
    - standard or stabilized networks forcing k and x within Muskingum valid ranges
    """

    configs: Configs
    network: Network  # ids, topology, k, x, the partition, and the stability analysis
    runoff: Runoff  # the weight table gridded runoff from the configs is aggregated with

    # calculated muskingum coefficients - consumed by static kernels or modified by dynamic kernels
    c1: FloatArray  # n x 1 - C1 values for each segment => f(k, x, dt_routing)
    c2: FloatArray  # n x 1 - C2 values for each segment => f(k, x, dt_routing)
    c3: FloatArray  # n x 1 - C3 values for each segment => f(k, x, dt_routing)
    c4: FloatArray  # n x 1 - C4 values for each segment => f(k, x, dt_routing) - used if lateral inflow provided
    c4_dt: FloatArray  # n x 1 - c4 / dt_runoff, turns catchment runoff volumes into a uniform lateral inflow rate
    subdivisions: IntArray  # n x 1 - equal sub-reaches each river is routed as, all 1 on a standard network
    reach_indptr: IntArray  # (n + 1,) offsets of each river's sub-reach states when stabilized, else empty
    substeps: IntArray  # n x 1 - steps each river takes per routing step when stabilized, else empty

    # The parallelizable groups of rivers the kernels will solve
    routing_jobs: tuple[tuple[IntArray, IntArray, int, int], ...]  # (block_starts, block_stops, outlet, region)
    cut_target: IntArray  # (n_regions,) river each region's outlet drains into, -1 at a basin outlet
    threads: int = 1  # the threads passed to the most recent route(), read by writers such as zarr_writer

    # State variables
    channel_state: FloatArray  # one value per river, or per sub-reach on a stabilized network
    _ensemble_member_states: list[FloatArray]  # for ensemble routing

    # Time options
    dt_routing: int  # compute time step, must be divisible into and <= min(dt_runoff, dt_discharge)
    dt_runoff: int  # time between catchment runoff steps, must be 1) constant, divisible into and <= dt_total
    dt_discharge: int  # time to average routed flows and save them, must be <=
    dt_total: int  # how long to simulate,
    num_runoff_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_discharge: int
    logger: logging.Logger

    # methods overridable via dependency injection
    _discharge_writer: WriteDischargesFn

    def __init__(self, configs: Configs, network: Network | None = None, runoff: Runoff | None = None) -> None:
        """
        Args:
            configs: options describing the simulation, from ``Configs(...)`` or ``Configs.from_json(path)``. They
                are validated for routing when ``route`` is called. A Router takes its options only from this
                object, so build a new Configs and a Router from it to change them.
            network: the network to route over. One is built from the params file the first time it is needed
                when none is given, so pass one to reuse a parsed and partitioned network across Routers.
            runoff: the Runoff that aggregates ``runoff_files`` to catchments. It must be the class for the
                ``runoff_type`` config. One is built from the configs the first time it is needed when none is given,
                so pass one to reuse a read weight table across Routers.
        """
        if not isinstance(configs, Configs):
            raise TypeError('Router must be given an rr.Configs object')
        if not isinstance(network, Network) and network is not None:
            raise TypeError(f'network must be a Network or None, got {type(network).__name__}')
        if not isinstance(runoff, Runoff) and runoff is not None:
            raise TypeError(f'runoff must be a Runoff or None, got {type(runoff).__name__}')
        if runoff is not None and configs.forcing == 'runoff':
            expected = RUNOFF_CLASS_FOR_RUNOFF_TYPE[configs.runoff_type]
            if not isinstance(runoff, expected):
                raise TypeError(
                    f'runoff_type {configs.runoff_type!r} is read by {expected.__name__}, got {type(runoff).__name__}'
                )

        self.configs = configs
        self.network = network if network is not None else Network.from_configs(configs)
        self.runoff = runoff  # built from the configs the first time runoff is routed when None

        # configure logging - progress bar and info/debug logs are mutually exclusive
        self.logger = build_logger(self.configs, 'router')
        self.logger.debug('Logger initialized')

        # default discharge writer; overridable via set_discharge_writer
        self._discharge_writer = zarr_writer
        self.reach_indptr = np.zeros(0, dtype=np.int64)
        self.substeps = np.zeros(0, dtype=np.int64)
        return

    def __repr__(self) -> str:
        return f'{type(self).__name__}(params_file={self.configs.params_file!r})'

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
        state = pd.read_parquet(state_file, columns=['Q'])['Q'].to_numpy(dtype=np.float32)
        # a stabilized network may start from one value per sub-reach, which is checked once the layout is known
        if state.shape[0] != n_rivers and self.configs.network_type != 'stabilized':
            raise ValueError(
                f'channel_state_init_file has {state.shape[0]} values but {self.configs.params_file} has '
                f'{n_rivers} rivers. The state file must have one row per river in the same order.'
            )
        self.channel_state = state
        return

    def _check_channel_state(self) -> None:
        """
        Require channel_state to have one value per routed reach: one per river on a standard network, one per
        sub-reach on a stabilized one. Without an initial state file a stabilized network starts from zero in every
        sub-reach. A state file that does not match the layout raises; its values are never changed to fit.
        """
        if not self.reach_indptr.shape[0]:
            return
        n_reaches = int(self.reach_indptr[-1])
        if not self.configs.channel_state_init_file:
            self.channel_state = np.zeros(n_reaches, dtype=np.float32)
            return
        if self.channel_state.shape[0] != n_reaches:
            raise ValueError(
                f'channel_state_init_file has {self.channel_state.shape[0]} values, but the stabilized network has '
                f'{n_reaches} sub-reaches. The state file must have one row per sub-reach in the order a final '
                f'state file is written.'
            )
        return

    def _write_final_state(self) -> None:
        final_state_file = self.configs.channel_state_final_file
        if not final_state_file:
            return
        self.logger.debug('Writing Final State to Parquet')
        ids = self.network.river_ids
        ids = np.repeat(ids, np.diff(self.reach_indptr)) if self.reach_indptr.shape[0] else ids
        pd.DataFrame({'river_id': ids, 'Q': self.channel_state}).to_parquet(self.configs.channel_state_final_file)
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
        self.dt_runoff = self.dt_discharge  # Assigned to pass time validation. It is never used in channel routing
        self._validate_time_options()
        return

    def _set_forced_time_options(self, dates: DatetimeArray) -> None:
        self.dt_runoff = self.configs.dt_runoff or (dates[1] - dates[0]).astype('timedelta64[s]').astype(int)
        self.dt_discharge = self.configs.dt_discharge or self.dt_runoff
        self.dt_total = self.configs.dt_total or self.dt_runoff * dates.shape[0]
        if self.configs.dt_routing:
            self.dt_routing = self.configs.dt_routing
        elif self.configs.network_type == 'stabilized':
            self.dt_routing = self.network.largest_stable_dt(self.dt_runoff)
            self.logger.info(f'dt_routing was not provided; using {self.dt_routing} s, the largest stable divisor')
        else:
            self.logger.warning('dt_routing was not provided or is Null/False, defaulting to dt_runoff')
            self.dt_routing = self.dt_runoff
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
        if self.configs.network_type == 'stabilized':
            # a river too long for dt_routing is split into N sub-reaches that are each k/N long and take 1/N of its
            # lateral inflow, so one set of coefficients serves all of them. A river too short for it is sub-cycled
            # in m steps of dt_routing/m, so its coefficients are built for that step.
            self.subdivisions, substeps, _ = self.network.conditioning(self.dt_routing)
            self.reach_indptr = np.zeros(self.network.size + 1, dtype=np.int64)
            np.cumsum(self.subdivisions, out=self.reach_indptr[1:])
            self.substeps = substeps if np.any(substeps > 1) else np.zeros(0, dtype=np.int64)
            k = (self.network.k / self.subdivisions).astype(np.float32)
            dt_river = (self.dt_routing / substeps).astype(np.float32)  # float32 like the standard path's math
            self.logger.info(
                f'Stabilized network: {self.network.size} rivers routed as {int(self.reach_indptr[-1])} sub-reaches '
                f'at dt_routing={self.dt_routing} s, {int(np.count_nonzero(substeps > 1))} of them sub-cycled'
            )
        else:
            self.subdivisions = np.ones(self.network.size, dtype=np.int64)
            self.reach_indptr = np.zeros(0, dtype=np.int64)
            self.substeps = np.zeros(0, dtype=np.int64)
            k = self.network.k
            dt_river = self.dt_routing
        dt_div_k = dt_river / k
        denominator = dt_div_k + (2 * (1 - self.network.x))
        _2x = 2 * self.network.x
        # contiguous arrays iterate faster in kernels due to cpu and ram access patterns
        self.c1 = np.ascontiguousarray((dt_div_k - _2x) / denominator, dtype=np.float32)
        self.c2 = np.ascontiguousarray((dt_div_k + _2x) / denominator, dtype=np.float32)
        self.c3 = np.ascontiguousarray(((2 * (1 - self.network.x)) - dt_div_k) / denominator, dtype=np.float32)
        self.c4 = np.ascontiguousarray(self.c1 + self.c2, dtype=np.float32)
        self.c4_dt = np.ascontiguousarray(self.c4 / self.dt_runoff, dtype=np.float32)
        if self.reach_indptr.shape[0]:
            self.c4_dt = np.ascontiguousarray(self.c4_dt / self.subdivisions, dtype=np.float32)
        self.network.check_stability(
            self.dt_routing,
            action=self.configs.unstable_coefficients,
            subdivisions=self.subdivisions if self.reach_indptr.shape[0] else None,
            substeps=self.substeps if self.substeps.shape[0] else None,
        )
        if not np.allclose(self.c1 + self.c2 + self.c3, 1):
            self.logger.warning('Muskingum coefficients do not sum to 1')
            self.logger.debug(f'c1: {self.c1}')
            self.logger.debug(f'c2: {self.c2}')
            self.logger.debug(f'c3: {self.c3}')
            raise ValueError('Muskingum coefficients do not sum to 1, check routing parameters and time step')

        self._check_channel_state()
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
            threads: number of regions the network is split into when ``thread_pool`` is given. Gaussian grid
                runoff is aggregated inside those regions' routing passes. Also the concurrency limit of writers
                that follow it, like zarr_writer.

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
        # there are two types of loops, one for channel only, one if runoff forcing is provided.
        if self.configs.forcing == 'channel':
            self._execute_routing_channel(thread_pool)
        else:
            self._execute_routing_forced(thread_pool)
        return

    def _execute_routing_channel(self, thread_pool: ThreadPoolExecutor | None) -> None:
        self.logger.info('-' * 60)
        dispatcher = resolve_dispatcher(self.configs)

        # channel routing time options and (static) coefficients
        self._set_channel_time_options()
        self._set_static_muskingum_coefficients()

        self.logger.debug('Starting routing computation')
        q_t = self.channel_state.astype(np.float32, copy=True)
        discharge_array = self._new_discharge_buffer()
        dispatcher(self, q_t, discharge_array, thread_pool)
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

    def _execute_routing_forced(self, thread_pool: ThreadPoolExecutor | None) -> None:
        self._ensemble_member_states = []
        runoff_files = self.configs.runoff_files
        dispatcher = resolve_dispatcher(self.configs)  # raises before any runoff is read if there is no kernel
        if self.runoff is None:
            self.runoff = RUNOFF_CLASS_FOR_RUNOFF_TYPE[self.configs.runoff_type].from_configs(self.configs)
        if self.network.synthetic is not None:
            self.runoff.distribute(self.network)
        runoff_iter = self.runoff.generator(runoff_files)  # yields what the runoff_type's kernel reads

        total_files = len(runoff_files)
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
            self.logger.info(f'Routing catchment runoff: {runoff_file}')
            self._set_forced_time_options(dates)
            if self.num_runoff_steps > dates.shape[0]:
                raise ValueError(
                    f'dt_total={self.dt_total} s needs {self.num_runoff_steps} steps of catchment runoff but '
                    f'{runoff_file} provides {dates.shape[0]}. Lower dt_total or provide more input.'
                )
            if self.num_runoff_steps < dates.shape[0]:
                self.logger.debug(f'Using the first {self.num_runoff_steps} of {dates.shape[0]} steps for dt_total')
                dates = dates[: self.num_runoff_steps]
                if not isinstance(forcing, CellRunoff):  # the fused kernel reads only the steps it routes
                    forcing = forcing[:, : self.num_runoff_steps]
            # static coefficients depend only on dt_routing/dt_runoff, so rebuild them only when those change
            # across files (dynamic coefficients are rebuilt inside the kernel each substep)
            if self.configs.coefficients != 'dynamic' and (self.dt_routing, self.dt_runoff) != coeff_dt:
                self._set_static_muskingum_coefficients()
                coeff_dt = (self.dt_routing, self.dt_runoff)
            self.logger.debug('Starting routing computation')
            q_t = self.channel_state.astype(np.float32, copy=True)
            q_array = self._new_discharge_buffer()
            dispatcher(self, q_t, q_array, forcing, thread_pool)
            if self.configs.runoff_processing_mode == 'sequential':
                self.logger.debug('Updating Channel State for Next Sequential Computation')
                self.channel_state = q_t
            elif self.configs.runoff_processing_mode == 'ensemble':
                self.logger.debug('Recording Member State for Final State Aggregation')
                self._ensemble_member_states.append(q_t.copy())

            if self.dt_discharge > self.dt_runoff:
                # todo reduce to dt_discharge inside the kernel. That drops this resample pass and sizes the buffer to
                # the output: 0.44 GB rather than 10.6 GB for a year of hourly routing on the Amazon, where the
                # reshape-and-mean below costs 3.8 to 5.1 s on its own.
                self.logger.debug('Resampling dates and discharges to specified timestep')
                q_array = q_array.reshape(
                    (q_array.shape[0], int(self.dt_total / self.dt_discharge), int(self.dt_discharge / self.dt_runoff))
                ).mean(axis=2)
                dates = dates[:: self.num_runoff_steps_per_discharge]

            self.logger.debug('Writing Discharge Array to File')
            self._discharge_writer(self, dates, q_array, discharge_file, runoff_file)

        if self.configs.runoff_processing_mode == 'ensemble':
            self.channel_state = np.array(self._ensemble_member_states).mean(axis=0)
        self.logger.info('-' * 60)
        return

    def _new_discharge_buffer(self) -> FloatArray:
        """A zeroed C-order (river, time) array for the routed discharge of the original, non-synthetic rivers."""
        synthetic = self.network.synthetic
        n_rivers = self.network.river_ids.shape[0] if synthetic is None else np.count_nonzero(~synthetic)
        return np.zeros((n_rivers, self.num_runoff_steps), dtype=np.float32)

    ################################################
    # Dependency injection methods for users to overwrite default behaviors without subclassing
    ################################################

    def set_discharge_writer(self, func: WriteDischargesFn) -> Self:
        """Set how discharge results are saved to disc. See ._discharge_writer for function signature."""
        self._discharge_writer = func
        return self
