import datetime
from typing import Any

import netCDF4 as nc
import numpy as np
import pandas as pd
import tqdm
from scipy.sparse import csc_matrix
from scipy.sparse import diags
from scipy.sparse import eye
from scipy.sparse.linalg import factorized

from ..tools import adjacency_matrix
from .typing import FloatArray, IntArray, DatetimeArray, FactorizedSolveFn, PathInput, WriteDischargesFn
from .RoutingConfigs import RoutingConfigs

__all__ = ['Muskingum', ]


class Muskingum(RoutingConfigs):
    """
    Muskingum channel routing of discharge without lateral inflow.

    Required config keys:
        routing_params_file: parquet file with river_id, downstream_river_id, k, x columns
        initial_state_file: parquet file with Q and R columns for initial conditions
        dt_routing: routing computation timestep in seconds
        dt_discharge: output timestep in seconds (defaults to dt_routing)
        dt_total: total simulation duration in seconds
        discharge_file: output netCDF file path
        final_state_file: parquet file path to write final state

    Optional config keys:
        start_datetime: simulation start datetime string (defaults to '2000-01-01')
    """
    # Network dependent matrices and vectors from routing parameters file
    A: csc_matrix  # n x n - adjacency matrix
    river_ids: IntArray  # n x 1 - river ID for each segment
    downstream_river_ids: IntArray  # n x 1 - downstream river ID for each segment
    k: FloatArray  # n x 1 - K values for each segment
    x: FloatArray  # n x 1 - X values for each segment

    # Network and routing timestep dependent matrices and vectors
    c1: FloatArray  # n x 1 - C1 values for each segment => f(k, x, dt_routing)
    c2: FloatArray  # n x 1 - C2 values for each segment => f(k, x, dt_routing)
    c3: FloatArray  # n x 1 - C3 values for each segment => f(k, x, dt_routing)
    lhs: csc_matrix  # n x n - left hand side of matrix form of routing equation
    lhs_factorized: FactorizedSolveFn  # callable factorized left hand side matrix for direct

    # State variables
    initial_state: tuple[FloatArray, FloatArray]  # Q init, R init

    # Time options
    dt_total: float
    dt_discharge: float
    dt_routing: float

    # methods that are overridable via dependency injection
    _write_discharges: WriteDischargesFn

    def __init__(self, config_file: PathInput | None = None, **kwargs: Any):
        super().__init__(config_file, **kwargs)

    def _set_network_dependent_vectors(self) -> None:
        self.logger.debug('Calculating network dependent vectors')
        df = pd.read_parquet(self.conf['routing_params_file'])
        if not {'river_id', 'k', 'x', 'downstream_river_id'}.issubset(df.columns):
            raise ValueError(
                'routing_params_file must have 4 columns with names river_id, k, x, and downstream_river_id')

        if df[self.conf['var_river_id']].duplicated().any():
            raise ValueError('routing_params_file contains duplicate river IDs.')

        # todo how to recycle this so that other routing classes don't repeat but still get the same validation and
        #  setting of class properties? Not store as class properties? set in routing? Have a list of control vocab
        #  expected column names as a class property then this loops over those and sets set[key] = df[key].to_numpy(
        #  )? perhaps we only store column names and expected dtypes? Try/Except this block.
        self.river_ids = df[self.conf['var_river_id']].to_numpy(dtype=np.int64, copy=False)
        self.downstream_river_ids = df['downstream_river_id'].to_numpy(dtype=np.int64, copy=False)
        self.k = df['k'].to_numpy(dtype=np.float64, copy=False)
        self.x = df['x'].to_numpy(dtype=np.float64, copy=False)

        river_id_set = set(self.river_ids.tolist())
        downstream_ids = set(self.downstream_river_ids.tolist()) - {-1}
        unknown_downstream_ids = sorted(downstream_ids - river_id_set)
        if unknown_downstream_ids:
            raise ValueError(
                f'routing_params_file has downstream IDs not in river_id column: {unknown_downstream_ids[:10]}')
        self.A = adjacency_matrix(self.river_ids, self.downstream_river_ids)
        return

    def _set_muskingum_coefficients(self, dt_routing: float) -> None:
        self.logger.debug('Calculating Muskingum coefficients')
        dt_div_k = dt_routing / self.k
        denominator = dt_div_k + (2 * (1 - self.x))
        _2x = 2 * self.x
        self.c1 = (dt_div_k + _2x) / denominator
        self.c2 = (dt_div_k - _2x) / denominator
        self.c3 = ((2 * (1 - self.x)) - dt_div_k) / denominator

        if not np.allclose(self.c1 + self.c2 + self.c3, 1):
            self.logger.warning('Muskingum coefficients do not sum to 1')
            self.logger.debug(f'c1: {self.c1}')
            self.logger.debug(f'c2: {self.c2}')
            self.logger.debug(f'c3: {self.c3}')
            raise ValueError('Muskingum coefficients do not sum to 1, check routing parameters and time step')

        self.lhs = eye(self.A.shape[0]) - (diags(self.c2) @ self.A)
        self.lhs = self.lhs.tocsc()
        self.logger.info('Calculating factorized LHS matrix')
        self.lhs_factorized = factorized(self.lhs)
        return

    def _read_initial_state(self) -> None:
        if hasattr(self, 'initial_state'):
            return

        state_file = self.conf.get('initial_state_file', '')
        if state_file == '':
            # todo raise error on muskingum but warning on others? but avoid duplicating code?
            # raise RuntimeError('initial_state_file not provided. Cannot route without water in the channels')
            self.logger.warning('initial_state_file not provided. Defaulting to zero initial conditions')
            initial_state = (np.zeros(self.A.shape[0], dtype=np.float64), np.zeros(self.A.shape[0], dtype=np.float64))
            self.initial_state = initial_state
            return
        self.logger.debug('Reading Initial State from Parquet')
        initial_state = pd.read_parquet(state_file).values
        initial_state = (initial_state[:, 0].flatten(), initial_state[:, 1].flatten())
        self.initial_state = initial_state
        return

    def _write_final_state(self) -> None:
        final_state_file = self.conf.get('final_state_file', '')
        if final_state_file == '':
            return
        self.logger.debug('Writing Final State to Parquet')
        array = np.array(self.initial_state).T
        pd.DataFrame(array, columns=['Q', 'R']).to_parquet(self.conf['final_state_file'])
        return

    def route(self) -> 'Muskingum':
        self.logger.info('Beginning routing')
        t1 = datetime.datetime.now()

        self._validate_configs()
        self._set_network_dependent_vectors()
        self.logger.debug(self)

        # Time parameters from config
        self.dt_routing = self.conf['dt_routing']
        self.dt_total = self.conf['dt_total']
        self.dt_discharge = self.conf.get('dt_discharge', self.dt_routing)

        assert self.dt_total >= self.dt_discharge >= self.dt_routing, 'Must have dt_total >= dt_discharge >= dt_routing'
        assert self.dt_total % self.dt_discharge == 0, 'dt_total must be an integer multiple of dt_discharge'
        assert self.dt_discharge % self.dt_routing == 0, 'dt_discharge must be an integer multiple of dt_routing'

        num_output_steps = int(self.dt_total / self.dt_discharge)
        num_routing_per_output = int(self.dt_discharge / self.dt_routing)

        self._set_muskingum_coefficients(self.dt_routing)
        self._read_initial_state()

        discharge_array = self._router(num_output_steps, num_routing_per_output)

        # Generate date array for output
        start = np.datetime64(self.conf.get('start_datetime', '2000-01-01'))
        dt_sec = int(self.dt_discharge)
        dates = start + (np.arange(num_output_steps) * np.timedelta64(dt_sec, 's'))

        self.logger.info('Writing Discharge Array to File')
        np.round(discharge_array, decimals=2, out=discharge_array)
        discharge_array = discharge_array.astype(np.float32, copy=False)
        # todo update call signatures to be inclusive of variable amounts of other file paths that subclasses pass
        self._write_discharges(dates, discharge_array, self.conf['discharge_file'])
        self._write_final_state()

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return self

    def _router(self, num_output_steps: int, num_routing_per_output: int) -> FloatArray:
        """
        Route discharge without lateral inflow.

        Muskingum channel routing equation:
            (I - c2*A) @ Q(t+1) = c1*(A @ Q(t)) + c3*Q(t)
        """
        self.logger.debug('Getting initial state arrays')
        q_init, _ = self.initial_state

        n = self.A.shape[0]
        discharge_array = np.zeros((num_output_steps, n), dtype=np.float64)
        q_t = q_init.astype(np.float64, copy=True)
        rhs = np.zeros(n, dtype=np.float64)
        buffer = np.zeros(n, dtype=np.float64)
        interval_sum = np.zeros(n, dtype=np.float64)

        t1 = datetime.datetime.now()
        output_iter = range(num_output_steps)
        if self.conf['progress_bar']:
            output_iter = tqdm.tqdm(output_iter, desc='Channel Routing')
        else:
            self.logger.info('Performing routing computation iterations')

        for output_step in output_iter:
            interval_sum.fill(0.0)
            for _ in range(num_routing_per_output):
                # rhs = c1*(A @ q_t) + c3*q_t
                buffer[:] = self.A @ q_t
                np.multiply(self.c1, buffer, out=rhs)
                np.multiply(self.c3, q_t, out=buffer)
                np.add(rhs, buffer, out=rhs)
                q_t[:] = self.lhs_factorized(rhs)
                interval_sum += q_t
            discharge_array[output_step, :] = interval_sum / num_routing_per_output

        discharge_array[discharge_array < 0] = 0

        self.logger.debug('Updating Initial State')
        self.initial_state = (q_t, np.zeros_like(q_t))

        t2 = datetime.datetime.now()
        self.logger.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')
        return discharge_array

    def _write_discharges(self,
                          dates: DatetimeArray,
                          discharge_array: FloatArray,
                          discharge_file: PathInput,
                          runoff_file: PathInput, ) -> None:
        """
        Writes routed discharge from a routing simulation to a netcdf file.
        You can overwrite this method with a custom handler using set_write_discharges.

        Args:
            dates: datetime array corresponding to the discharge rows
            discharge_array: routed discharge values with shape (time, river)
            discharge_file: the file path to write the discharge data to
            runoff_file: the file path to the runoff file used to generate the discharge values

        Returns:
            None
        """
        # todo: resample before getting to this function.
        # todo: how can this function be made flexible to pass other file paths so that other classes don't have to
        #  repeat the same code just to pass more optional arguments??
        dates = dates[::self.num_runoff_steps_per_discharge].astype('datetime64[s]')
        with nc.Dataset(str(discharge_file), mode='w', format='NETCDF4') as ds:
            ds.createDimension('time', size=discharge_array.shape[0])
            ds.createDimension(self.conf['var_river_id'], size=discharge_array.shape[1])
            ds.runoff_file = str(runoff_file)
            time_var = ds.createVariable('time', 'f8', ('time',))
            time_var.units = f'seconds since {pd.Timestamp(dates[0]).strftime("%Y-%m-%d %H:%M:%S")}'
            time_var[:] = (dates - dates[0]).astype('timedelta64[s]').astype(np.int64)
            id_var = ds.createVariable(self.conf['var_river_id'], 'i4', (self.conf['var_river_id']), )
            id_var[:] = self.river_ids
            flow_var = ds.createVariable(self.conf['var_discharge'], 'f4', ('time', self.conf['var_river_id']))
            flow_var[:] = discharge_array
            flow_var.long_name = 'Discharge at catchment outlet'
            flow_var.standard_name = 'discharge'
            flow_var.aggregation_method = 'mean'
            flow_var.units = 'm3 s-1'
        return

    def set_write_discharges(self, func: WriteDischargesFn) -> 'Muskingum':
        """
        Overwrites the default write_discharges method to a custom function and returns the class instance so that you
        can chain the method with the constructor.

        Args:
            func (callable): function that takes dates, discharge_array, discharge_file, runoff_file and returns None

        Returns:
            river_route.LumpedMuskingum
        """
        self._write_discharges = func
        return self
