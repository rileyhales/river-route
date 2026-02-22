import datetime
from typing import Any, Generator, Self

import numpy as np
import pandas as pd
import tqdm

from .AbstractRouter import AbstractRouter
from .types import DatetimeArray, FloatArray, PathInput
from ..transformers import AbstractBaseTransformer, Transformer, SCSUnitHydrograph

__all__ = ['UnitMuskingum']

GeneratorSignature = Generator[tuple[DatetimeArray, FloatArray, PathInput, PathInput], None, None]


class UnitMuskingum(AbstractRouter):
    """
    Muskingum routing with pluggable unit-hydrograph runoff transformation.

    Runoff depths/volumes are first convolved with a unit hydrograph kernel before being injected
    as lateral inflow into the Muskingum channel routing system.

    Unit hydrograph transformer config keys
    ----------------------------------------
    uh_type        : string name of the transformer to build (e.g. 'scs').  Required unless
                     a transformer is injected via set_transformer() or loaded via uh_kernel_file.
    uh_kernel_file : path to a parquet file (n_basins, n_time_steps) to warm-start the kernel.
                     When present, from_cached() is used instead of computing from network params.
    uh_state_file  : path to a parquet file (n_basins, n_time_steps) to warm-start the state.
                     Only used when uh_kernel_file is also provided.

    Input data config keys (choose one source)
    -------------------------------------------
    # option 1 - pre-computed per-catchment volumes
    runoff_volumes_files : list of 1 or more paths to netCDF files containing per-catchment volume
                           time series (m^3 per timestep) for each river segment. Intentionally
                           named differently from TeleportMuskingum's catchment_volumes_files to
                           keep the two routers' configs independent.
    # option 2 - raw runoff depths converted via a weight table
    runoff_depths_files  : list of 1 or more paths to netCDF files containing gridded runoff depth
                           time series. Requires weight_table_file. Intentionally named differently
                           from TeleportMuskingum's runoff_depths_files for the same reason.
    weight_table_file    : path to a netCDF weight table used to aggregate gridded depths to
                           per-catchment volumes. Required when using runoff_depths_files.

    discharge_files      : list of 1 or more paths to netCDF output files for routed discharge.
                           Count must match the number of input files.

    Required routing params columns: river_id, downstream_river_id, k, x, tc, area_sqm
    """

    _network_time_signature: tuple[Any, ...] | None = None
    _ensemble_member_states: list[FloatArray]

    tc: FloatArray
    area: FloatArray
    _runoff_transformer: AbstractBaseTransformer | None = None

    # Time options
    dt_runoff: float
    num_routing_steps: int
    num_routing_steps_per_runoff: int
    num_runoff_steps_per_discharge: int

    def __init__(self, config_file: PathInput | None = None, **kwargs: Any) -> None:
        super().__init__(config_file, **kwargs)

    # ------------------------------------------------------------------
    # Dependency injection
    # ------------------------------------------------------------------

    def set_transformer(self, transformer: AbstractBaseTransformer) -> Self:
        """
        Inject a user instantiated transformer object including custom subclasses of AbstractBaseTransformer allowing
        users to implement their own transformation logic in a custom transformer. The transformer must be a
        AbstractBaseTransformer subclass with a kernel and state set.

        You must call this before route() to have any affect.
        """
        if not isinstance(transformer, AbstractBaseTransformer):
            raise TypeError(f'transformer must be a AbstractBaseTransformer, got {type(transformer).__name__}')
        self._runoff_transformer = transformer
        return self

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate_router_configs(self) -> None:
        # todo: normalize runoff_volumes_files and runoff_depths_files from str/Path to list,
        #  following the same pattern as TeleportMuskingum._validate_router_configs
        # todo: raise ValueError if both runoff_volumes_files and runoff_depths_files are provided
        # todo: raise ValueError if len(discharge_files) != n_input_files
        # todo: raise FileNotFoundError for any input file that does not exist
        # todo: raise NotADirectoryError for any discharge_files output directory that does not exist
        # todo: raise ValueError if runoff_depths_files is set but weight_table_file is not provided
        # todo: raise FileNotFoundError if weight_table_file does not exist
        if (
                self._runoff_transformer is None
                and not self.conf.get('uh_type')
                and not self.conf.get('uh_kernel_file')
        ):
            raise RuntimeError(
                'No transformer configured. '
                'Set uh_type (and optionally uh_kernel_file/uh_state_file) in config, '
                'or call set_transformer() before routing.'
            )

    # ------------------------------------------------------------------
    # Network / time setup
    # ------------------------------------------------------------------

    def _set_network_dependent_vectors(self) -> None:
        super()._set_network_dependent_vectors()
        df = pd.read_parquet(self.conf['routing_params_file'])
        required = {'tc', 'area_sqm'}
        if not required.issubset(df.columns):
            raise ValueError(f'routing_params_file must also have columns: {sorted(required)}')
        self.tc = df['tc'].to_numpy(dtype=np.float64, copy=False)
        self.area = df['area_sqm'].to_numpy(dtype=np.float64, copy=False)
        if np.any(self.tc < 0):
            raise ValueError('tc (time of concentration) must be non-negative for all catchments')
        if np.any(self.area <= 0):
            raise ValueError('area_sqm must be > 0 for all catchments')

    def _set_network_and_time_dependent_vectors(self, dates: DatetimeArray) -> None:
        self.logger.debug('Setting and validating time parameters')
        self.dt_runoff = self.conf.get('dt_runoff', (dates[1] - dates[0]).astype('timedelta64[s]').astype(int))
        self.dt_discharge = self.conf.get('dt_discharge', self.dt_runoff)
        self.dt_total = self.conf.get('dt_total', self.dt_runoff * dates.shape[0])
        if not self.conf.get('dt_routing', 0):
            self.logger.warning('dt_routing was not provided or is Null/False, defaulting to dt_runoff')
        self.dt_routing = self.conf.get('dt_routing', self.dt_runoff)

        signature = (self.dt_total, self.dt_runoff, self.dt_discharge, self.dt_routing,)
        if self._network_time_signature == signature:
            return

        try:
            assert self.dt_total >= self.dt_runoff, 'dt_total must be >= dt_runoff'
            assert self.dt_total >= self.dt_discharge, 'dt_total must be >= dt_discharge'
            assert self.dt_discharge >= self.dt_runoff, 'dt_discharge must be >= dt_runoff'
            assert self.dt_runoff >= self.dt_routing, 'dt_runoff must be >= dt_routing'
            assert self.dt_total % self.dt_runoff == 0, 'dt_total must be an integer multiple of dt_runoff'
            assert self.dt_total % self.dt_discharge == 0, 'dt_total must be an integer multiple of dt_discharge'
            assert self.dt_discharge % self.dt_runoff == 0, 'dt_discharge must be an integer multiple of dt_runoff'
            assert self.dt_runoff % self.dt_routing == 0, 'dt_runoff must be an integer multiple of dt_routing'
        except AssertionError as e:
            self.logger.error(e)
            raise AssertionError('Time options are not valid')

        self.num_runoff_steps_per_discharge = int(self.dt_discharge / self.dt_runoff)
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)
        self.num_routing_steps = int(self.dt_total / self.dt_routing)

        self._set_muskingum_coefficients(self.dt_routing)
        self._network_time_signature = signature
        self._initialize_runoff_transformer()

    def _initialize_runoff_transformer(self) -> None:
        # If a transformer was already injected via set_transformer, leave it alone
        # todo should this be a setter, updater, initializer? i don't want people to change the transformer mid
        #  route. we would need functions to resample the kernel and state if the dt changes.
        if self._runoff_transformer is not None:
            return

        uh_kernel_file = self.conf.get('uh_kernel_file')
        if uh_kernel_file:
            # Kernel already computed: load directly — no uh_type needed
            self._runoff_transformer = Transformer.from_kernel(dt=self.dt_runoff, kernel=uh_kernel_file)
        else:
            # Build kernel from network params using the transformer named by uh_type
            uh_type = str(self.conf.get('uh_type', '')).lower()
            if uh_type == 'scs':
                self._runoff_transformer = SCSUnitHydrograph(dt=self.dt_runoff, tc=self.tc, area=self.area)
            else:
                raise ValueError(f'Unknown uh_type: {uh_type!r}. Valid options: "scs"')

        if self.conf.get('uh_state_file'):
            self._runoff_transformer.set_state(self.conf['uh_state_file'])

    # ------------------------------------------------------------------
    # Final state
    # ------------------------------------------------------------------

    def _write_final_state(self) -> None:
        final_state_file = self.conf.get('final_state_file', '')
        if final_state_file == '':
            return
        if self.conf['input_type'] == 'sequential':
            array = np.array(self.initial_state).T
        elif self.conf['input_type'] == 'ensemble':
            array = np.array(self._ensemble_member_states).mean(axis=0)
        else:
            raise ValueError('Unexpected input type when writing final states')
        self.logger.debug('Writing Final State to Parquet')
        pd.DataFrame(array, columns=['Q', 'R']).to_parquet(self.conf['final_state_file'])

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_runoff(self) -> GeneratorSignature:
        # todo: implement following the same two-option pattern as TeleportMuskingum._volumes_generator:
        #
        # option 1 — pre-computed volumes:
        #   iterate over zip(self.conf['runoff_volumes_files'], self.conf['discharge_files']),
        #   open each netCDF, read dates and the volume array, yield
        #   (dates, volumes_array, volume_file, discharge_file)
        #
        # option 2 — depths + weight table:
        #   iterate over zip(self.conf['runoff_depths_files'], self.conf['discharge_files']),
        #   call depth_to_volume(...) using self.conf['weight_table_file'] and the depth file,
        #   yield (dates, volumes_array, depth_file, discharge_file)
        #
        # note: config keys are intentionally different from TeleportMuskingum's so that the two
        #  routers can be configured independently and a config file cannot be accidentally used
        #  with the wrong router class.
        raise NotImplementedError('todo: implement _load_runoff')

    # ------------------------------------------------------------------
    # Routing
    # ------------------------------------------------------------------

    def route(self) -> Self:
        self.logger.info('Beginning routing')
        t1 = datetime.datetime.now()

        self._validate_configs()
        self._set_network_dependent_vectors()
        self.logger.debug(self)

        for dates, volumes_array, runoff_file, discharge_file in self._load_runoff():
            self.logger.info('-' * 80)
            self._set_network_and_time_dependent_vectors(dates)
            volumes_array = volumes_array.astype(np.float64, copy=False)
            # todo: confirm unit convention expected by _load_runoff — currently assumes m^3/timestep
            #  and converts to m^3/s here to match the Muskingum lateral inflow units
            np.divide(volumes_array, self.dt_runoff, out=volumes_array)
            discharge_array = self._router(dates, volumes_array)

            if self.dt_discharge > self.dt_runoff:
                self.logger.info('Resampling dates and discharges to specified timestep')
                discharge_array = (
                    discharge_array
                    .reshape((
                        int(self.dt_total / self.dt_discharge),
                        int(self.dt_discharge / self.dt_runoff),
                        self.A.shape[0],
                    ))
                    .mean(axis=1)
                )
                dates = dates[::self.num_runoff_steps_per_discharge]

            self.logger.info('Writing Discharge Array to File')
            np.round(discharge_array, decimals=2, out=discharge_array)
            discharge_array = discharge_array.astype(np.float32, copy=False)
            self._write_discharges(dates, discharge_array, discharge_file, runoff_file)

        self._write_final_state()

        t2 = datetime.datetime.now()
        self.logger.info('All runoff files routed')
        self.logger.info(f'Total job time: {(t2 - t1).total_seconds()}')
        return self

    def _router(self, dates: DatetimeArray, volumes: FloatArray) -> FloatArray:
        # todo why is dates not used?
        self.logger.debug('Getting initial state arrays')
        self._read_initial_state()
        q_init = self.initial_state[0]
        self._ensemble_member_states = []

        T, n = volumes.shape
        discharge_array = np.zeros((T, n), dtype=np.float64)
        q_t = q_init.astype(np.float64, copy=True)
        rhs = np.zeros(n, dtype=np.float64)
        work = np.zeros(n, dtype=np.float64)
        interval_sum = np.zeros(n, dtype=np.float64)

        if self._runoff_transformer is None:
            raise RuntimeError('Runoff transformer has not been initialized')
        transformer = self._runoff_transformer

        t1 = datetime.datetime.now()
        runoff_iter = range(T)
        if self.conf['progress_bar']:
            runoff_iter = tqdm.tqdm(runoff_iter, desc='Unit-Muskingum Routing', total=T)
        else:
            self.logger.info('Performing Unit-Muskingum routing iterations')

        for t in runoff_iter:
            ql_t = transformer.transform(volumes[t, :])
            interval_sum.fill(0.0)
            for _ in range(self.num_routing_steps_per_runoff):
                work[:] = self.A @ q_t
                np.multiply(self.c1, work, out=rhs)
                np.multiply(self.c3, q_t, out=work)
                np.add(rhs, work, out=rhs)
                np.add(rhs, ql_t, out=rhs)
                q_t[:] = self.lhs_factorized(rhs)
                interval_sum += q_t
            discharge_array[t, :] = interval_sum / self.num_routing_steps_per_runoff

        discharge_array[discharge_array < 0] = 0

        if self.conf['input_type'] == 'sequential':
            self.logger.debug('Updating Initial State for Next Sequential Computation')
            self.initial_state = (q_t, np.zeros_like(q_t))
        if self.conf['input_type'] == 'ensemble':
            self.logger.debug('Recording Member State for Final State Aggregation')
            self._ensemble_member_states.append(np.array([q_t, np.zeros_like(q_t)]).T)

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return discharge_array
