import datetime
from typing import Any

import numpy as np
import pandas as pd
import tqdm

from .TeleportMuskingum import TeleportMuskingum
from .typing import DatetimeArray, FloatArray, PathInput
from ..tools import adjacency_matrix
from ..transformers import AbstractBaseTransformer, Transformer, SCSUnitHydrograph

__all__ = ['UnitMuskingum']


class UnitMuskingum(TeleportMuskingum):
    """
    Muskingum routing with pluggable unit-hydrograph runoff transformation.

    Config keys
    -----------
    uh_type   : string name of the transformer to build (e.g. 'scs').  Required unless
                a transformer is injected via set_transformer() or loaded via uh_kernel.
    uh_kernel : path to a parquet file (n_basins, n_time_steps) to warm-start the kernel.
                When present, from_cached() is used instead of computing from network params.
    uh_state  : path to a parquet file (n_basins, n_time_steps) to warm-start the state.
                Only used when uh_kernel is also provided.

    Required routing params columns: river_id, downstream_river_id, k, x, tc, area_sqm
    """

    tc: FloatArray
    area: FloatArray
    _runoff_transformer: AbstractBaseTransformer | None = None

    def __init__(self, config_file: PathInput | None = None, **kwargs: Any) -> None:
        super().__init__(config_file, **kwargs)

    # ------------------------------------------------------------------
    # Dependency injection
    # ------------------------------------------------------------------

    def set_transformer(self, transformer: AbstractBaseTransformer) -> None:
        """
        Inject a user instantiated transformer object including custom subclasses of AbstractBaseTransformer allowing
        users to implement their own transformation logic in a custom transformer. The transformer must be a
        AbstractBaseTransformer subclass with a kernel and state set.

        Call this before route() to use custom transformation logic.
        """
        if not isinstance(transformer, AbstractBaseTransformer):
            raise TypeError(f'transformer must be a AbstractBaseTransformer, got {type(transformer).__name__}')
        self._runoff_transformer = transformer

    # ------------------------------------------------------------------
    # Network setup
    # ------------------------------------------------------------------

    def _set_network_dependent_vectors(self) -> None:
        self.logger.debug('Calculating network dependent vectors (UnitMuskingum)')
        df = pd.read_parquet(self.conf['routing_params_file'])

        required = {'river_id', 'downstream_river_id', 'k', 'x', 'tc', 'area_sqm'}
        if not required.issubset(df.columns):
            raise ValueError(
                f'routing_params_file must have columns: {sorted(required)}. '
                f'Found: {sorted(df.columns.tolist())}'
            )

        if df[self.conf['var_river_id']].duplicated().any():
            raise ValueError('routing_params_file contains duplicate river IDs.')

        self.river_ids = df[self.conf['var_river_id']].to_numpy(dtype=np.int64, copy=False)
        self.downstream_river_ids = df['downstream_river_id'].to_numpy(dtype=np.int64, copy=False)
        self.k = df['k'].to_numpy(dtype=np.float64, copy=False)
        self.x = df['x'].to_numpy(dtype=np.float64, copy=False)
        self.tc = df['tc'].to_numpy(dtype=np.float64, copy=False)
        self.area = df['area_sqm'].to_numpy(dtype=np.float64, copy=False)

        if np.any(self.tc < 0):
            raise ValueError('tc (time of concentration) must be non-negative for all catchments')
        if np.any(self.area <= 0):
            raise ValueError('area_sqm must be > 0 for all catchments')

        river_id_set = set(self.river_ids.tolist())
        downstream_ids = set(self.downstream_river_ids.tolist()) - {-1}
        unknown_downstream_ids = sorted(downstream_ids - river_id_set)
        if unknown_downstream_ids:
            raise ValueError(
                f'routing_params_file has downstream IDs not in river_id column: {unknown_downstream_ids[:10]}'
            )
        self.A = adjacency_matrix(self.river_ids, self.downstream_river_ids)

    def _set_network_and_time_dependent_vectors(self, dates: DatetimeArray) -> None:
        old_signature = self._network_time_signature
        super()._set_network_and_time_dependent_vectors(dates)
        if self._network_time_signature != old_signature:
            self._initialize_runoff_transformer()

    def _initialize_runoff_transformer(self) -> None:
        # If a transformer was already injected via set_transformer, leave it alone
        if self._runoff_transformer is not None:
            return

        uh_kernel = self.conf.get('uh_kernel')
        if uh_kernel:
            # Kernel already computed: load directly — no uh_type needed
            self._runoff_transformer = Transformer.from_kernel(dt=self.dt_runoff, kernel=uh_kernel)
        else:
            # Build kernel from network params using the transformer named by uh_type
            uh_type = str(self.conf.get('uh_type', '')).lower()
            if uh_type == 'scs':
                self._runoff_transformer = SCSUnitHydrograph(dt=self.dt_runoff, tc=self.tc, area=self.area)
            else:
                raise ValueError(f'Unknown uh_type: {uh_type!r}. Valid options: "scs"')

        if self.conf.get('uh_state'):
            self._runoff_transformer.set_state(self.conf['uh_state'])

    # ------------------------------------------------------------------
    # Routing
    # ------------------------------------------------------------------

    def route(self) -> 'UnitMuskingum':
        # Validate that a transformer will be available before entering the routing loop.
        # A transformer is valid if it was injected via set_transformer(), or the config
        # supplies uh_type (to build from params) or uh_kernel (to load from cache).
        if (
                self._runoff_transformer is None
                and not self.conf.get('uh_type')
                and not self.conf.get('uh_kernel')
        ):
            raise RuntimeError(
                'No transformer configured. '
                'Set uh_type (and optionally uh_kernel/uh_state) in config, '
                'or call set_transformer() before routing.'
            )
        return super().route()

    def _read_initial_state(self) -> None:
        super()._read_initial_state()

    def _write_final_state(self) -> None:
        super()._write_final_state()

    def _router(self, dates: DatetimeArray, volumes: FloatArray) -> FloatArray:
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
