import logging
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Self

import numpy as np
import pandas as pd
from tqdm import tqdm

from .. import kernels
from .._logging import PROGRESS, configure_logging
from ..configs import Configs, check_channel_state_rivers, is_dev_null
from ..forcings import RUNOFF_CLASS_FOR_FORCING, Runoff
from ..network import Network
from ..types import DatetimeArray, Float64Array, FloatArray, Forcing, PathInput, WriteDischargesFn
from . import static_muskingum
from .static_muskingum import Coefficients
from .writers import null_writer, zarr_writer

__all__ = ['Router']

logger = logging.getLogger(__name__)

# the share of a file's discharge memory added for the files after it. Memory is paged in only where it is written, so
# it costs nothing until a longer file, such as a leap year after a common one, routes into it
DISCHARGE_MEMORY_HEADROOM = 0.125


class Router:
    """
    Muskingum river routing with static coefficients, each river's k and x fixed for the whole run
    - forcing: the channel alone, or runoff on grids or meshes or already aggregated to catchments, read from the
      forcing files
    - validity: dt_routing auto gives each river its own valid step, and conditioning routes rivers still invalid at
      their step in the sub-reaches or subcycles that make every coefficient non-negative
    """

    configs: Configs
    network: Network  # ids, topology, k, x, the partition, and the coefficient signs
    runoff: Runoff | None  # reads forcing_files, built from the configs the first time runoff is routed

    # what static_muskingum prepared for the current time steps
    coefficients: Coefficients  # each river's sub-reaches and the Muskingum coefficients of each

    # State variables
    channel_state: FloatArray  # one value per river, or with conditioning one per sub-reach

    # Time options, in seconds
    dt_routing: int | None  # every river's step when given, a divisor of dt_runoff; None for each river's own
    dt_runoff: int  # time between runoff steps, read from each forcing file, when rivers hand discharge downstream
    dt_discharge: int  # time between the routed discharges written, an integer divisor of dt_total
    dt_total: int  # how long to simulate
    num_runoff_steps: int
    num_runoff_steps_per_discharge: int

    # methods overridable via dependency injection
    _discharge_writer: WriteDischargesFn

    def __init__(self, configs: Configs) -> None:
        """
        Args:
            configs: options describing the simulation, from ``Configs(...)`` or ``Configs.from_json(path)``. They
                are validated for routing when ``route`` is called. A Router takes its options only from this
                object, and builds its network from the network file here and its runoff reader the first time
                runoff is routed, so build a new Configs and a Router from it to change them.
        """
        if not isinstance(configs, Configs):
            raise TypeError('Router must be given an rr.Configs object')
        self.configs = configs
        self.network = Network.from_configs(configs)
        self.runoff = None

        # point the river_route logger at these log options; Configs turns the progress bar off when logging is off
        configure_logging(self.configs)

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
            logger.warning('channel_state_init_file not provided. Defaulting to zero initial conditions')
            self.channel_state = np.zeros(n_rivers, dtype=np.float32)
            return
        logger.debug('Reading Initial State from Parquet')
        state = pd.read_parquet(state_file)
        missing = [column for column in ('riverId', 'Q') if column not in state.columns]
        if missing:
            raise ValueError(f'channel_state_init_file {state_file} is missing the column(s) {", ".join(missing)}')
        # conditioning starts from one value per sub-reach, whose count is checked once the coefficients are prepared
        check_channel_state_rivers(
            state['riverId'].to_numpy(), self.network.river_ids, self.configs.conditioning, state_file
        )
        self.channel_state = state['Q'].to_numpy(dtype=np.float32)
        return

    def _check_channel_state(self) -> None:
        """
        Require channel_state to have one value per routed reach: one per river, or with conditioning one per
        sub-reach. Without an initial state file conditioning starts from zero in every sub-reach. Any other state
        that does not match the sub-reaches raises, whether it was read from a state file or carried from the previous
        forcing file; its values are never changed to fit.
        """
        n_reaches = int(self.coefficients.sub_reach_starts[-1])
        if self.channel_state.shape[0] == n_reaches:
            return
        if not self.configs.channel_state_init_file and not self.channel_state.any():
            self.channel_state = np.zeros(n_reaches, dtype=np.float32)  # a zero state is zero for any sub-reaches
            return
        raise ValueError(
            f'the channel state has {self.channel_state.shape[0]} values, but conditioning routes {n_reaches} '
            f'sub-reaches at a runoff step of {self.dt_runoff} s. A channel_state_init_file must have one row per '
            f'sub-reach in the order a final state file is written, and the sub-reaches cannot change between the '
            f'forcing files of one sequential run, which they do when the files have different time steps.'
        )

    def _write_final_state(self) -> None:
        final_state_file = self.configs.channel_state_final_file
        if not final_state_file:
            return
        logger.debug('Writing Final State to Parquet')
        ids = np.repeat(self.network.river_ids, np.diff(self.coefficients.sub_reach_starts))  # a row per sub-reach
        pd.DataFrame({'riverId': ids, 'Q': self.channel_state}).to_parquet(self.configs.channel_state_final_file)
        return

    ################################################
    # Prepare arrays for routing
    ################################################

    def _set_channel_time_options(self) -> None:
        self.dt_routing = None if self.configs.dt_routing == 'auto' else self.configs.dt_routing
        self.dt_total = self.configs.dt_total
        self.dt_discharge = self.configs.dt_discharge or self.dt_routing or 0  # validate_routing requires one of them
        self.dt_runoff = self.dt_discharge  # with no runoff, rivers hand their discharge downstream every dt_discharge
        self._validate_time_options()
        return

    def _set_forced_time_options(self, dates: DatetimeArray, forcing_file: PathInput) -> None:
        if dates.shape[0] < 2:
            raise ValueError(f'{forcing_file} has a single time step, so its dt_runoff cannot be read from it')
        self.dt_runoff = (dates[1] - dates[0]).astype('timedelta64[s]').astype(int)
        if self.dt_runoff <= 0:
            raise ValueError(f'the time steps of {forcing_file} must increase, but its first is {self.dt_runoff} s')
        logger.debug(f'dt_runoff={self.dt_runoff} s, read from the time steps of {forcing_file}')
        # runoff with irregular time steps is for its provider to regularize; routing never resamples it
        steps = np.diff(dates)
        irregular = np.flatnonzero(steps != np.timedelta64(self.dt_runoff, 's'))
        if irregular.size:
            first = int(irregular[0])
            raise ValueError(
                f'{forcing_file} steps {steps[first].astype("timedelta64[s]")} from {dates[first]} to '
                f'{dates[first + 1]} after a first step of {self.dt_runoff} s, but runoff is routed in uniform '
                f'time steps'
            )
        self.dt_discharge = self.configs.dt_discharge or self.dt_runoff
        self.dt_total = self.configs.dt_total or self.dt_runoff * dates.shape[0]
        self.dt_routing = None if self.configs.dt_routing == 'auto' else self.configs.dt_routing
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
        if self.dt_routing is not None and self.dt_runoff < self.dt_routing:
            raise ValueError('dt_runoff must be >= dt_routing')
        # check that time options are evenly divisible
        if self.dt_total % self.dt_runoff != 0:
            raise ValueError('dt_total must be an integer multiple of dt_runoff')
        if self.dt_total % self.dt_discharge != 0:
            raise ValueError('dt_total must be an integer multiple of dt_discharge')
        if self.dt_discharge % self.dt_runoff != 0:
            raise ValueError('dt_discharge must be an integer multiple of dt_runoff')
        if self.dt_routing is not None and self.dt_runoff % self.dt_routing != 0:
            raise ValueError('dt_runoff must be an integer multiple of dt_routing')

        # Now that we know time parameters are valid, set time-derived parameters for computation cycles
        self.num_runoff_steps = int(self.dt_total / self.dt_runoff)
        self.num_runoff_steps_per_discharge = int(self.dt_discharge / self.dt_runoff)  # to resample/reshape results
        return

    def _prepare_coefficients(self) -> None:
        """Each river's sub-reaches and the coefficients of each for the current dt_runoff and dt_routing."""
        logger.debug('Preparing the Muskingum coefficients')
        self.coefficients = static_muskingum.prepare_routing(
            self.network, self.configs, self.dt_runoff, self.dt_routing
        )
        logger.info(
            f'{self.network.size} rivers routed as {int(self.coefficients.sub_reach_starts[-1])} sub-reaches, each '
            f'handing its discharge downstream every {self.dt_runoff} s'
        )
        self._check_channel_state()
        return

    ################################################
    # Methods to execute routing simulation
    ################################################

    def route(self, thread_pool: ThreadPoolExecutor | None = None, threads: int = 1, kernel: str = 'auto') -> Self:
        """
        Execute the simulation described by the configs, over the network. The configs are validated for routing
        first. The thread pool, thread count, and kernel are runtime resources, not configs, so they are given here.

        Args:
            thread_pool: optional pool to route and aggregate gridded runoff concurrently on. It is used as given and
                never shut down here, so it can be shared and closed by the caller's with block. Without one,
                routing is single-threaded regardless of ``threads``.
            threads: number of jobs the region's blocks are packed into when ``thread_pool`` is given. Grid
                runoff is aggregated inside those jobs as each river is routed, except cumulative or clipped runoff,
                which is aggregated first. The discharge writer is handed both, and zarr_writer writes its chunks on
                the pool.
            kernel: the routing kernel of river_route.kernels: 'auto' routes with the C kernel once ``rr compile`` has
                compiled it on this machine, and with numba until then; 'c' or 'numba' chooses one.

        Returns:
            Self: the class instance with updated channel_state and output files written to disk
        """
        if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
            raise ValueError(f'threads must be an integer >= 1, got {threads!r}')
        kernel = kernels.choose_kernel(kernel)
        # start timer
        logger.log(PROGRESS, 'Beginning routing')
        started = time.perf_counter()
        logger.debug('Validating configs')
        self.configs.validate_routing()
        self._select_discharge_writer()
        logger.debug(self)
        logger.info(f'Routing with the {kernel} kernel')
        # read state, route, write state
        self._read_initial_state()
        if self.configs.forcing == 'channel':
            self._execute_routing_channel(thread_pool, threads, kernel)
        else:
            self._execute_routing_forced(thread_pool, threads, kernel)
        self._write_final_state()
        logger.log(PROGRESS, f'Routing completed in {time.perf_counter() - started:.3f} seconds')
        return self

    def _execute_routing_channel(self, thread_pool: ThreadPoolExecutor | None, threads: int, kernel: str) -> None:
        logger.info('-' * 60)
        self._set_channel_time_options()
        self._prepare_coefficients()

        logger.debug('Starting routing computation')
        channel_state = self.channel_state.astype(np.float32, copy=True)
        memory = np.zeros(self._discharge_memory_size(), dtype=np.float32)
        discharge_array = self._route_region(channel_state, memory, None, thread_pool, threads, kernel)
        self.channel_state = channel_state

        # generate dates since they cannot be copied from an external forcing file
        dates = pd.date_range(
            start=self.configs.start_datetime,
            periods=self.num_runoff_steps,
            freq=pd.to_timedelta(self.dt_discharge, unit='s'),
        ).to_numpy()

        # write outputs
        logger.debug('Writing Discharge Array to File')
        discharge_file = self.configs.discharge_files[0]
        self._discharge_writer(
            self, dates, discharge_array, discharge_file, '', thread_pool=thread_pool, threads=threads
        )
        logger.info('-' * 60)
        return

    def _execute_routing_forced(self, thread_pool: ThreadPoolExecutor | None, threads: int, kernel: str) -> None:
        forcing_files = self.configs.forcing_files
        if self.runoff is None:
            self.runoff = RUNOFF_CLASS_FOR_FORCING[self.configs.forcing].from_configs(self.configs)
        if self.network.synthetic is not None:
            self.runoff.distribute(self.network)
        runoff_iter = self.runoff.generator(forcing_files)  # yields the runoff routing reads, in its forcing's form

        total_files = len(forcing_files)
        file_iter = (
            (dates, runoff, forcing_file, discharge_file)
            for (dates, runoff, forcing_file), discharge_file in zip(
                runoff_iter, self.configs.discharge_files, strict=True
            )
        )
        if self.configs.progress_bar:
            file_iter = tqdm(file_iter, total=total_files, desc='Files Routed')

        prepared_for_time_steps: tuple[int | None, int] | None = None  # (dt_routing, dt_runoff) of the coefficients
        member_state_sum: Float64Array | None = None  # the sum of the final states of the ensemble members
        discharge_memory = np.empty(0, dtype=np.float32)  # every file's routed discharge, grown only for a longer file
        for dates, runoff, forcing_file, discharge_file in file_iter:
            logger.info('-' * 60)
            logger.info(f'Routing forcing file: {forcing_file}')
            self._set_forced_time_options(dates, forcing_file)
            if self.num_runoff_steps > dates.shape[0]:
                raise ValueError(
                    f'dt_total={self.dt_total} s needs {self.num_runoff_steps} steps of catchment runoff but '
                    f'{forcing_file} provides {dates.shape[0]}. Lower dt_total or provide more input.'
                )
            if self.num_runoff_steps < dates.shape[0]:
                logger.debug(f'Using the first {self.num_runoff_steps} of {dates.shape[0]} steps for dt_total')
                dates = dates[: self.num_runoff_steps]
                runoff = runoff.first_steps(self.num_runoff_steps)
            # the coefficients depend only on dt_routing and dt_runoff, so they are rebuilt only when those change
            if (self.dt_routing, self.dt_runoff) != prepared_for_time_steps:
                self._prepare_coefficients()
                prepared_for_time_steps = (self.dt_routing, self.dt_runoff)
            logger.debug('Starting routing computation')
            channel_state = self.channel_state.astype(np.float32, copy=True)
            # the kernels write every value, so each file is routed into the memory of the files before it rather than
            # into a new array, which the operating system would page in again
            n_values = self._discharge_memory_size()
            if discharge_memory.shape[0] < n_values:
                discharge_memory = np.zeros(int(n_values * (1 + DISCHARGE_MEMORY_HEADROOM)), dtype=np.float32)
            routed_discharge = self._route_region(
                channel_state, discharge_memory[:n_values], runoff, thread_pool, threads, kernel
            )
            if self.configs.runoff_processing_mode == 'sequential':
                logger.debug('Updating Channel State for Next Sequential Computation')
                self.channel_state = channel_state
            else:  # every ensemble member starts from the same state, and the final state is the mean of theirs
                member_state_sum = (
                    channel_state.astype(np.float64) if member_state_sum is None else member_state_sum + channel_state
                )

            dates, discharge_array = self._discharge_at_end_of_discharge_steps(dates, routed_discharge)
            logger.debug('Writing Discharge Array to File')
            self._discharge_writer(
                self, dates, discharge_array, discharge_file, forcing_file, thread_pool=thread_pool, threads=threads
            )

        if member_state_sum is not None:
            logger.debug('Averaging the final states of the ensemble members')
            self.channel_state = (member_state_sum / total_files).astype(np.float32)
        logger.info('-' * 60)
        return

    def _discharge_memory_size(self) -> int:
        """The floats of the memory one file is routed into: a row per river of the network, each with the discharge at
        the start and at the end of every runoff step."""
        return self.network.size * (self.num_runoff_steps + 1)

    def _route_region(
        self,
        channel_state: FloatArray,
        memory: FloatArray,
        forcing: Forcing,
        thread_pool: ThreadPoolExecutor | None,
        threads: int,
        kernel: str,
    ) -> FloatArray:
        """
        Route one forcing file, or the channel alone, over the network into memory with the prepared coefficients.
        Returns the routed discharge of the original rivers at the end of each runoff step, a (river, time) array with
        each river's row contiguous.
        """
        discharge = memory.reshape(self.network.size, self.num_runoff_steps + 1)
        kernels.route_region(
            kernel, self.network, self.coefficients, channel_state, discharge, forcing, thread_pool, threads
        )
        synthetic = self.network.synthetic  # the sub-reaches stabilize added, whose discharge is never written
        return discharge[:, 1:] if synthetic is None else discharge[~synthetic, 1:]

    def _discharge_at_end_of_discharge_steps(
        self, dates: DatetimeArray, discharge_array: FloatArray
    ) -> tuple[DatetimeArray, FloatArray]:
        """The routed discharge at the end of each dt_discharge, that of its last runoff step, and the date each
        dt_discharge begins at, as a (river, time) array with each river's row contiguous."""
        if self.dt_discharge == self.dt_runoff:
            return dates, discharge_array
        logger.debug('Taking the discharge at the end of each dt_discharge')
        runoff_steps = self.num_runoff_steps_per_discharge
        return dates[::runoff_steps], np.ascontiguousarray(discharge_array[:, runoff_steps - 1 :: runoff_steps])

    ################################################
    # Dependency injection methods
    ################################################

    def _select_discharge_writer(self) -> None:
        """Override the user's discharge writer when every output is the null device"""
        if not all(is_dev_null(f) for f in self.configs.discharge_files):
            return
        logger.warning('Discharge output is the null device: discharge will be routed and then discarded')
        self._discharge_writer = null_writer
        return

    def set_discharge_writer(self, func: WriteDischargesFn) -> Self:
        """Set how discharge results are saved to disk, with the function signature ``types.WriteDischargesFn``"""
        self._discharge_writer = func
        return self
