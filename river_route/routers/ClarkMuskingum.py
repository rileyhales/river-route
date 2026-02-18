import datetime
import os

import netCDF4 as nc
import numpy as np
import pandas as pd
import tqdm
from numpy.typing import NDArray
from typing import Any

from .TeleportMuskingum import TeleportMuskingum
from ..tools import adjacency_matrix

__all__ = ['ClarkMuskingum', ]

FloatArray = NDArray[np.float64]
DatetimeArray = NDArray[np.datetime64]
PathInput = str


class ClarkMuskingum(TeleportMuskingum):
    # Clark-specific parameters from routing params
    tc: FloatArray  # n x 1 - time of concentration (seconds) per catchment
    R: FloatArray  # n x 1 - Clark linear reservoir storage coefficient (seconds) per catchment

    # Clark-specific state
    _time_area_curves: FloatArray | None = None # (n_points, n_catchments) cumulative area fractions
    _uh_matrix: FloatArray | None = None  # (max_uh_len, n_catchments) zero-padded UH ordinates
    _clark_buffer: FloatArray | None = None  # (max_uh_len, n_catchments) rolling convolution buffer

    def __init__(self, config_file: PathInput | None = None, **kwargs: Any) -> None:
        super().__init__(config_file, **kwargs)

    def set_time_area_curves(self, curves: FloatArray) -> 'ClarkMuskingum':
        """
        Set per-catchment cumulative time-area curves.

        Args:
            curves: 2D array of shape (n_points, n_catchments) with cumulative area fractions
                    at evenly spaced normalized times over [0, 1]. Values must be monotonically
                    non-decreasing in each column, bounded in [0, 1], starting near 0 and ending at 1.

        Returns:
            self for method chaining
        """
        curves = np.asarray(curves, dtype=np.float64)
        if curves.ndim != 2:
            raise ValueError('time_area_curves must be a 2D array with shape (n_points, n_catchments)')
        if np.any(curves < 0) or np.any(curves > 1):
            raise ValueError('time_area_curves values must be in [0, 1]')
        if np.any(np.diff(curves, axis=0) < -1e-12):
            raise ValueError('time_area_curves must be monotonically non-decreasing along axis 0')
        self._time_area_curves = curves
        return self

    @staticmethod
    def _default_scs_time_area_curve(n_points: int = 101) -> FloatArray:
        """
        Generate the SCS dimensionless cumulative time-area curve.

        A(t/tc) = 1.414 * (t/tc)^1.5           for t/tc <= 0.5
        A(t/tc) = 1 - 1.414 * (1 - t/tc)^1.5   for t/tc > 0.5

        Args:
            n_points: number of evenly spaced points over [0, 1]

        Returns:
            1D array of cumulative area fractions, shape (n_points,)
        """
        t_norm = np.linspace(0, 1, n_points)
        curve = np.where(
            t_norm <= 0.5,
            1.414 * t_norm ** 1.5,
            1 - 1.414 * (1 - t_norm) ** 1.5,
        )
        np.clip(curve, 0, 1, out=curve)
        return curve

    def _set_network_dependent_vectors(self) -> None:
        self.logger.debug('Calculating network dependent vectors (Clark)')
        df = pd.read_parquet(self.conf['routing_params_file'])

        required = {'river_id', 'downstream_river_id', 'k', 'x', 'tc', 'R'}
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
        self.R = df['R'].to_numpy(dtype=np.float64, copy=False)

        if np.any(self.tc < 0):
            raise ValueError('tc (time of concentration) must be non-negative for all catchments')
        if np.any(self.R <= 0):
            raise ValueError('R (reservoir storage coefficient) must be positive for all catchments')

        river_id_set = set(self.river_ids.tolist())
        downstream_ids = set(self.downstream_river_ids.tolist()) - {-1}
        unknown_downstream_ids = sorted(downstream_ids - river_id_set)
        if unknown_downstream_ids:
            raise ValueError(
                f'routing_params_file has downstream IDs not in river_id column: {unknown_downstream_ids[:10]}'
            )
        self.A = adjacency_matrix(self.river_ids, self.downstream_river_ids)

        # read time-area histograms from file if provided (and not already set programmatically)
        if self._time_area_curves is None and self.conf.get('time_area_file', ''):
            self._read_time_area_file()

    def _read_time_area_file(self) -> None:
        """
        Read per-catchment time-area histograms from a parquet file.

        Expected format:
            - Columns are river_id values (int)
            - Rows are evenly-spaced bins over normalized time [0, 1]
            - Values are incremental area fractions (non-negative, each column sums to ~1)

        The incremental histograms are converted to cumulative curves via cumsum
        with a leading zero row prepended.
        """
        ta_file = self.conf['time_area_file']
        self.logger.info(f'Reading time-area histograms from: {ta_file}')
        if not os.path.exists(ta_file):
            raise FileNotFoundError(f'time_area_file not found: {ta_file}')

        ta_df = pd.read_parquet(ta_file)
        ta_df.columns = ta_df.columns.astype(np.int64)

        # align columns to river_ids ordering
        missing = set(self.river_ids.tolist()) - set(ta_df.columns.tolist())
        if missing:
            raise ValueError(
                f'time_area_file is missing {len(missing)} river_ids, '
                f'first 10: {sorted(missing)[:10]}'
            )
        histograms = ta_df[self.river_ids].to_numpy(dtype=np.float64)  # (n_bins, n_catchments)

        if np.any(histograms < 0):
            raise ValueError('time_area_file contains negative values')

        # normalize each column to sum to 1
        col_sums = histograms.sum(axis=0)
        col_sums[col_sums == 0] = 1.0  # avoid division by zero
        histograms = histograms / col_sums

        # convert incremental histograms to cumulative curves with a leading zero row
        cumulative = np.cumsum(histograms, axis=0)
        cumulative = np.vstack([np.zeros((1, cumulative.shape[1])), cumulative])

        self._time_area_curves = cumulative

    def _set_network_and_time_dependent_vectors(self, dates: DatetimeArray) -> None:
        old_signature = self._network_time_signature
        super()._set_network_and_time_dependent_vectors(dates)
        if self._network_time_signature != old_signature:
            self.precompute_unit_hydrographs()

    def precompute_unit_hydrographs(self) -> 'ClarkMuskingum':
        """
        Precompute the Clark unit hydrograph ordinates for each catchment.

        For each catchment:
        1. Discretize the time-area curve at dt_runoff intervals over [0, tc]
        2. Compute incremental areas via np.diff (translation hydrograph)
        3. Route through a linear reservoir: O[j] = C1*I[j] + C2*O[j-1]
        4. Truncate the reservoir tail where C2^j < 1e-6
        5. Normalize so UH sums to 1
        6. Pack into a zero-padded matrix for vectorized convolution

        Returns:
            self for method chaining
        """
        self.logger.info('Precomputing Clark unit hydrographs')
        dt = self.dt_runoff
        n = len(self.tc)

        # get or build time-area curves
        if self._time_area_curves is not None:
            if self._time_area_curves.ndim == 1:
                ta_curves = np.tile(self._time_area_curves[:, np.newaxis], (1, n))
            else:
                if self._time_area_curves.shape[1] != n:
                    raise ValueError(
                        f'time_area_curves has {self._time_area_curves.shape[1]} catchments, expected {n}'
                    )
                ta_curves = self._time_area_curves
        else:
            scs = self._default_scs_time_area_curve()
            ta_curves = np.tile(scs[:, np.newaxis], (1, n))

        n_curve_points = ta_curves.shape[0]
        t_norm_curve = np.linspace(0, 1, n_curve_points)

        uh_list = []
        for i in range(n):
            tc_i = self.tc[i]
            R_i = self.R[i]

            if tc_i > 0:
                n_tc_steps = max(int(np.ceil(tc_i / dt)), 1)
                t_steps = np.arange(n_tc_steps + 1) * dt
                t_norm_steps = np.clip(t_steps / tc_i, 0, 1)
                cum_area = np.interp(t_norm_steps, t_norm_curve, ta_curves[:, i])
                translation = np.diff(cum_area)
            else:
                translation = np.array([1.0])

            C1 = dt / (R_i + 0.5 * dt)
            C2 = 1 - C1

            if C2 > 0:
                tail_len = int(np.ceil(np.log(1e-6) / np.log(C2)))
            else:
                tail_len = 0

            total_len = len(translation) + tail_len
            uh = np.zeros(total_len)

            for j in range(total_len):
                inflow = translation[j] if j < len(translation) else 0.0
                uh[j] = C1 * inflow + (C2 * uh[j - 1] if j > 0 else 0.0)

            uh_sum = uh.sum()
            if uh_sum > 0:
                uh /= uh_sum

            uh_list.append(uh)

        max_uh_len = max(len(uh) for uh in uh_list)
        self._uh_matrix = np.zeros((max_uh_len, n), dtype=np.float64)
        for i, uh in enumerate(uh_list):
            self._uh_matrix[:len(uh), i] = uh

        self.logger.info(f'UH matrix shape: {self._uh_matrix.shape}')
        return self

    def _read_initial_state(self) -> None:
        super()._read_initial_state()

        if self._clark_buffer is not None:
            return

        state_file = self.conf.get('initial_state_file', '')
        clark_file = self._clark_state_path(state_file) if state_file else ''

        n = self.A.shape[0]
        max_uh_len = self._uh_matrix.shape[0]

        if clark_file and os.path.exists(clark_file):
            self.logger.debug('Reading Clark buffer state from file')
            loaded = np.load(clark_file)
            # handle shape changes if UH was recomputed with different max_uh_len
            if loaded.shape[0] >= max_uh_len:
                self._clark_buffer = loaded[:max_uh_len, :n].astype(np.float64)
            else:
                self._clark_buffer = np.zeros((max_uh_len, n), dtype=np.float64)
                self._clark_buffer[:loaded.shape[0], :] = loaded[:, :n]
        else:
            self.logger.debug('Setting Clark buffer to zeros')
            self._clark_buffer = np.zeros((max_uh_len, n), dtype=np.float64)

    def _write_final_state(self) -> None:
        super()._write_final_state()

        final_state_file = self.conf.get('final_state_file', '')
        if not final_state_file:
            return
        clark_file = self._clark_state_path(final_state_file)
        self.logger.debug('Writing Clark buffer state to file')
        np.save(clark_file, self._clark_buffer)

    @staticmethod
    def _clark_state_path(state_file: str) -> str:
        base, _ = os.path.splitext(state_file)
        return base + '_clark.npy'

    def _router(self, dates: DatetimeArray, volumes: FloatArray) -> FloatArray:
        """
        Clark-Muskingum routing: UH-transformed lateral inflows added directly to the RHS.

        Muskingum channel routing equation (upstream only):
            (I - c2*A) @ Q(t+1) = c1*(A @ Q(t)) + c3*Q(t) + ql(t)

        Unlike parent Muskingum, lateral inflows are NOT weighted by c1/c2.
        The UH transformation already converts runoff volumes into discharge rates.

        Two modes controlled by conf['precomputed_lateral_inflows']:
            False (default): volumes are raw runoff rates, transformed on-the-fly via rolling buffer
            True: volumes are already UH-transformed lateral inflow rates, used directly
        """
        self.logger.debug('Getting initial state arrays')
        self._read_initial_state()
        q_init = self.initial_state[0]
        self._ensemble_member_states = []

        precomputed = self.conf.get('precomputed_lateral_inflows', False)
        T, n = volumes.shape

        discharge_array = np.zeros((T, n), dtype=np.float64)
        q_t = q_init.astype(np.float64, copy=True)
        rhs = np.zeros(n, dtype=np.float64)
        work = np.zeros(n, dtype=np.float64)
        interval_sum = np.zeros(n, dtype=np.float64)

        buf = self._clark_buffer  # None when precomputed_lateral_inflows=True

        t1 = datetime.datetime.now()
        runoff_iter = range(T)
        if self.conf['progress_bar']:
            runoff_iter = tqdm.tqdm(runoff_iter, desc='Clark Routing', total=T)
        else:
            self.logger.info('Performing Clark-Muskingum routing iterations')

        for t in runoff_iter:
            # --- get lateral inflow for this timestep ---
            if precomputed:
                ql_t = volumes[t, :]
            else:
                buf += self._uh_matrix * volumes[t, :]
                ql_t = buf[0, :]

            # --- Muskingum channel routing: only upstream discharge ---
            interval_sum.fill(0.0)
            for _ in range(self.num_routing_steps_per_runoff):
                # rhs = c1*(A @ q_t) + c3*q_t + ql_t
                work[:] = self.A @ q_t
                np.multiply(self.c1, work, out=rhs)
                np.multiply(self.c3, q_t, out=work)
                np.add(rhs, work, out=rhs)
                np.add(rhs, ql_t, out=rhs)
                q_t[:] = self.lhs_factorized(rhs)
                interval_sum += q_t
            discharge_array[t, :] = interval_sum / self.num_routing_steps_per_runoff

            # --- advance the rolling buffer after reading row 0 ---
            if not precomputed:
                buf[:-1, :] = buf[1:, :]
                buf[-1, :] = 0.0

        discharge_array[discharge_array < 0] = 0

        if self.conf['input_type'] == 'sequential':
            self.logger.debug('Updating Initial State for Next Sequential Computation')
            self.initial_state = (q_t, np.zeros_like(q_t))
            if not precomputed:
                self._clark_buffer = buf
        if self.conf['input_type'] == 'ensemble':
            self.logger.debug('Recording Member State for Final State Aggregation')
            self._ensemble_member_states.append(np.array([q_t, np.zeros_like(q_t)]).T)

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return discharge_array

    def compute_lateral_inflows(self, output_files: list[str]) -> 'ClarkMuskingum':
        """
        Pre-compute Clark UH-transformed lateral inflows and write to netCDF files.

        Reads the configured volume/runoff input files, applies the UH rolling buffer
        transformation, and writes the result as catchment volume files that can be used
        with precomputed_lateral_inflows=True for a subsequent route() call.

        The output files store lateral_rate * dt_runoff so the standard route() pipeline
        (which divides by dt_runoff) recovers the correct rates.

        Args:
            output_files: list of output netCDF file paths, one per input file

        Returns:
            self for method chaining
        """
        self.logger.info('Pre-computing Clark lateral inflows')
        self._validate_configs()
        self._set_network_dependent_vectors()

        input_files = self.conf.get('catchment_volumes_files', self.conf.get('runoff_depths_files', []))
        if len(output_files) != len(input_files):
            raise ValueError(
                f'Number of output_files ({len(output_files)}) must match '
                f'number of input files ({len(input_files)})'
            )

        file_idx = 0
        for dates, volumes_array, runoff_file, _ in self._volumes_generator():
            self._set_network_and_time_dependent_vectors(dates)
            volumes_array = volumes_array.astype(np.float64, copy=False)
            np.divide(volumes_array, self.dt_runoff, out=volumes_array)

            T, n = volumes_array.shape
            lateral_rates = np.empty_like(volumes_array)

            if self._clark_buffer is None:
                self._clark_buffer = np.zeros((self._uh_matrix.shape[0], n), dtype=np.float64)
            buf = self._clark_buffer

            for t in range(T):
                buf += self._uh_matrix * volumes_array[t, :]
                lateral_rates[t, :] = buf[0, :]
                buf[:-1, :] = buf[1:, :]
                buf[-1, :] = 0.0

            self._clark_buffer = buf

            # write as volumes (rate * dt) so route() pipeline divides back to rates
            lateral_volumes = lateral_rates * self.dt_runoff

            out_path = output_files[file_idx]
            self.logger.info(f'Writing lateral inflows to: {out_path}')
            with nc.Dataset(str(out_path), mode='w', format='NETCDF4') as ds:
                ds.createDimension('time', size=T)
                ds.createDimension(self.conf['var_river_id'], size=n)

                time_var = ds.createVariable('time', 'f8', ('time',))
                t0 = pd.Timestamp(dates[0]).strftime("%Y-%m-%d %H:%M:%S")
                time_var.units = f'seconds since {t0}'
                time_var[:] = (dates - dates[0]).astype('timedelta64[s]').astype(np.int64)

                id_var = ds.createVariable(self.conf['var_river_id'], 'i4', (self.conf['var_river_id'],))
                id_var[:] = self.river_ids

                vol_var = ds.createVariable(
                    self.conf['var_catchment_volume'], 'f4', ('time', self.conf['var_river_id']))
                vol_var[:] = lateral_volumes.astype(np.float32)
                vol_var.long_name = 'Clark UH-transformed lateral inflow volume'
                vol_var.units = 'm3'

            file_idx += 1

        self._write_final_state()
        self.logger.info('All lateral inflow files written')
        return self
