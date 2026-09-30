import logging
import time
from concurrent.futures import ThreadPoolExecutor
from types import ModuleType
from typing import Self

import numpy as np
import pandas as pd
from tqdm import tqdm

from .._logging import PROGRESS, build_logger
from ..configs import Configs, is_dev_null
from ..network import Network
from ..runoff import RUNOFF_CLASS_FOR_FORCING, Runoff
from ..types import DatetimeArray, FloatArray, WriteDischargesFn
from . import dynamic_muskingum, static_muskingum
from ._numba_kernels import Layout, route_region
from .writers import null_writer, zarr_writer

__all__ = ['Router', 'ROUTING_METHOD_FOR_COEFFICIENTS']

# the module of the routing method each coefficients config chooses; each routes one river with its own parameters
ROUTING_METHOD_FOR_COEFFICIENTS: dict[str, ModuleType] = {'static': static_muskingum, 'dynamic': dynamic_muskingum}


class Router:
    """
    Muskingum style river routing allowing
    - coefficients: static (muskingum) or dynamic (muskingum-cunge style)
    - forcing: channel-only or runoff on grids/meshes or already aggregated to catchments
    - network conditioning: standard, or "stabilized" making edits to put k and x within valid ranges for dt_routing
    - transformation: the uniform runoff transformation (e.g. runoff transform method)
    """

    configs: Configs
    network: Network  # ids, topology, k, x, the partition, and the stability analysis
    runoff: Runoff | None  # reads runoff_files, built from the configs the first time runoff is routed when None

    # what the routing method the coefficients config chooses prepared for the current time steps
    routing_parameters: static_muskingum.StaticMuskingum | dynamic_muskingum.DynamicMuskingum  # per river
    layout: Layout  # each river routed whole, or in substeps and subcycles on a stabilized network

    threads: int = 1  # the threads passed to the most recent route(), read by writers such as zarr_writer

    # State variables
    channel_state: FloatArray  # one value per river, or per sub-reach on a stabilized network
    _ensemble_member_states: list[FloatArray]  # for ensemble routing

    # Time options, in seconds
    dt_routing: int  # routing step, an integer divisor of dt_runoff
    dt_runoff: int  # time between catchment runoff steps, constant, and an integer divisor of dt_discharge
    dt_discharge: int  # time routed discharge is averaged over before it is written, an integer divisor of dt_total
    dt_total: int  # how long to simulate
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
            network: the network to route over. One is built from the network file here when none is given, so
                pass one to reuse a parsed and partitioned network across Routers.
            runoff: the Runoff that aggregates ``runoff_files`` to catchments. It must be the class for the
                ``forcing`` config. One is built from the configs the first time it is needed when none is given,
                so pass one to reuse a read weight table across Routers.
        """
        if not isinstance(configs, Configs):
            raise TypeError('Router must be given an rr.Configs object')
        if not isinstance(network, Network) and network is not None:
            raise TypeError(f'network must be a Network or None, got {type(network).__name__}')
        if not isinstance(runoff, Runoff) and runoff is not None:
            raise TypeError(f'runoff must be a Runoff or None, got {type(runoff).__name__}')
        if runoff is not None and configs.forcing != 'channel':
            expected = RUNOFF_CLASS_FOR_FORCING[configs.forcing]
            if not isinstance(runoff, expected):
                raise TypeError(
                    f'forcing {configs.forcing!r} is read by {expected.__name__}, got {type(runoff).__name__}'
                )

        self.configs = configs
        self.network = network if network is not None else Network.from_configs(configs)
        self.runoff = runoff  # built from the configs the first time runoff is routed when None

        # configure logging; Configs turns the progress bar off when logging is off
        self.logger = build_logger(self.configs, 'router')
        self.logger.debug('Logger initialized')

        # default discharge writer; overridable via set_discharge_writer
        self._discharge_writer = zarr_writer
        return

    def __repr__(self) -> str:
        return f'{type(self).__name__}(network_file={self.configs.network_file!r})'

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
                f'channel_state_init_file has {state.shape[0]} values but {self.configs.network_file} has '
                f'{n_rivers} rivers. The state file must have one row per river in the same order.'
            )
        self.channel_state = state
        return

    def _check_channel_state(self) -> None:
        """
        Require channel_state to have one value per routed reach: one per river on a standard network, one per
        sub-reach on a stabilized one. Without an initial state file a stabilized network starts from zero in every
        sub-reach. Any other state that does not match the layout raises, whether it was read from a state file or
        carried from the previous runoff file; its values are never changed to fit.
        """
        if not self.layout.reach_indptr.shape[0]:
            return
        n_reaches = int(self.layout.reach_indptr[-1])
        if self.channel_state.shape[0] == n_reaches:
            return
        if not self.configs.channel_state_init_file and not self.channel_state.any():
            self.channel_state = np.zeros(n_reaches, dtype=np.float32)  # a zero state is zero in any layout
            return
        raise ValueError(
            f'the channel state has {self.channel_state.shape[0]} values, but the stabilized network has {n_reaches} '
            f'sub-reaches at dt_routing={self.dt_routing} s. A channel_state_init_file must have one row per sub-reach '
            f'in the order a final state file is written, and the sub-reaches cannot change between the runoff files '
            f'of one sequential run: set dt_routing to route every file at the same step.'
        )

    def _write_final_state(self) -> None:
        final_state_file = self.configs.channel_state_final_file
        if not final_state_file:
            return
        self.logger.debug('Writing Final State to Parquet')
        ids = self.network.river_ids
        reach_indptr = self.layout.reach_indptr
        ids = np.repeat(ids, np.diff(reach_indptr)) if reach_indptr.shape[0] else ids
        pd.DataFrame({'river_id': ids, 'Q': self.channel_state}).to_parquet(self.configs.channel_state_final_file)
        return

    ################################################
    # Prepare arrays for routing
    ################################################

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
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)
        return

    def _prepare_routing_method(self) -> None:
        """Lay out each river and build the routing method's parameters for the current dt_routing and dt_runoff."""
        self.logger.debug('Preparing the routing method')
        method = ROUTING_METHOD_FOR_COEFFICIENTS[self.configs.coefficients]
        self.layout, self.routing_parameters = method.prepare_routing(
            self.network, self.configs.network_type, self.dt_routing, self.dt_runoff, self.configs.unstable_coefficients
        )
        reach_indptr, subcycles = self.layout
        if reach_indptr.shape[0]:
            self.logger.info(
                f'Stabilized network: {self.network.size} rivers routed in substeps as {int(reach_indptr[-1])} '
                f'sub-reaches at dt_routing={self.dt_routing} s, and {int(np.count_nonzero(subcycles > 1))} rivers '
                f'routed in subcycles'
            )
        self._check_channel_state()
        return

    ################################################
    # Methods to execute routing simulation
    ################################################

    def route(self, thread_pool: ThreadPoolExecutor | None = None, threads: int = 1) -> Self:
        """
        Execute the simulation described by the configs, over the network. The configs are validated for routing
        first. The thread pool and thread count are runtime resources, not configs, so they are given here.

        Args:
            thread_pool: optional pool to route and aggregate gridded runoff concurrently on. It is used as given and
                never shut down here, so it can be shared and closed by the caller's with block. Without one,
                routing is single-threaded regardless of ``threads``.
            threads: number of jobs the region's blocks are packed into when ``thread_pool`` is given. Grid
                runoff is aggregated inside those jobs as each river is routed. Also the concurrency limit of writers
                that follow it, like zarr_writer.

        Returns:
            Self: the class instance with updated channel_state and output files written to disk
        """
        if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
            raise ValueError(f'threads must be an integer >= 1, got {threads!r}')
        self.threads = threads
        # start timer
        self.logger.log(PROGRESS, 'Beginning routing')
        started = time.perf_counter()
        self.logger.debug('Validating configs')
        self.configs.validate_routing()
        self._select_discharge_writer()
        self.logger.debug(self)
        # read state, route, write state
        self._read_initial_state()
        if self.configs.forcing == 'channel':
            self._execute_routing_channel(thread_pool)
        else:
            self._execute_routing_forced(thread_pool)
        self._write_final_state()
        self.logger.log(PROGRESS, f'Routing completed in {time.perf_counter() - started:.3f} seconds')
        return self

    def _execute_routing_channel(self, thread_pool: ThreadPoolExecutor | None) -> None:
        self.logger.info('-' * 60)
        self._set_channel_time_options()
        self._prepare_routing_method()

        self.logger.debug('Starting routing computation')
        q_t = self.channel_state.astype(np.float32, copy=True)
        discharge_array = self._new_discharge_buffer()
        route_region(self, q_t, discharge_array, None, thread_pool)
        self.channel_state = q_t

        # generate dates since they cannot be copied from an external forcing file
        dates = pd.date_range(
            start=self.configs.start_datetime,
            periods=self.num_runoff_steps,
            freq=pd.to_timedelta(self.dt_discharge, unit='s'),
        ).to_numpy()

        # write outputs
        self.logger.debug('Writing Discharge Array to File')
        self._discharge_writer(self, dates, discharge_array, self.configs.discharge_files[0], '')  # no runoff file
        self.logger.info('-' * 60)
        return

    def _execute_routing_forced(self, thread_pool: ThreadPoolExecutor | None) -> None:
        self._ensemble_member_states = []
        runoff_files = self.configs.runoff_files
        if self.runoff is None:
            self.runoff = RUNOFF_CLASS_FOR_FORCING[self.configs.forcing].from_configs(self.configs)
        if self.network.synthetic is not None:
            self.runoff.distribute(self.network)
        runoff_iter = self.runoff.generator(runoff_files)  # yields the runoff routing reads, in its forcing's form

        total_files = len(runoff_files)
        file_iter = (
            (dates, runoff, runoff_file, discharge_file)
            for (dates, runoff, runoff_file), discharge_file in zip(
                runoff_iter, self.configs.discharge_files, strict=True
            )
        )
        if self.configs.progress_bar:
            file_iter = tqdm(file_iter, total=total_files, desc='Files Routed')

        prepared_for_time_steps: tuple[int, int] | None = None  # (dt_routing, dt_runoff) of the routing parameters
        for dates, runoff, runoff_file, discharge_file in file_iter:
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
                runoff = runoff.first_steps(self.num_runoff_steps)
            # the routing parameters depend only on dt_routing and dt_runoff, so they are rebuilt only when those change
            if (self.dt_routing, self.dt_runoff) != prepared_for_time_steps:
                self._prepare_routing_method()
                prepared_for_time_steps = (self.dt_routing, self.dt_runoff)
            self.logger.debug('Starting routing computation')
            q_t = self.channel_state.astype(np.float32, copy=True)
            discharge_array = self._new_discharge_buffer()
            route_region(self, q_t, discharge_array, runoff, thread_pool)
            if self.configs.runoff_processing_mode == 'sequential':
                self.logger.debug('Updating Channel State for Next Sequential Computation')
                self.channel_state = q_t
            elif self.configs.runoff_processing_mode == 'ensemble':
                self.logger.debug('Recording Member State for Final State Aggregation')
                self._ensemble_member_states.append(q_t.copy())

            if self.dt_discharge > self.dt_runoff:
                # the kernels write one value per runoff step, so coarser outputs are averaged from them here
                self.logger.debug('Resampling dates and discharges to specified timestep')
                discharge_array = discharge_array.reshape(
                    (
                        discharge_array.shape[0],
                        int(self.dt_total / self.dt_discharge),
                        int(self.dt_discharge / self.dt_runoff),
                    )
                ).mean(axis=2)
                dates = dates[:: self.num_runoff_steps_per_discharge]

            self.logger.debug('Writing Discharge Array to File')
            self._discharge_writer(self, dates, discharge_array, discharge_file, runoff_file)

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
    # Dependency injection methods
    ################################################

    def _select_discharge_writer(self) -> None:
        """Override the user's discharge writer when every output is the null device"""
        if not all(is_dev_null(f) for f in self.configs.discharge_files):
            return
        self.logger.warning('Discharge output is the null device: discharge will be routed and then discarded')
        self._discharge_writer = null_writer
        return

    def set_discharge_writer(self, func: WriteDischargesFn) -> Self:
        """Set how discharge results are saved to disk, with the function signature ``types.WriteDischargesFn``"""
        self._discharge_writer = func
        return self
