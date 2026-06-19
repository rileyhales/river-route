import datetime
import json
import logging
import sys
from typing import Any, Self

import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from tqdm import tqdm

from ..logging import PROGRESS
from ..runoff import runoff_to_qlateral
from ..types import DatetimeArray, FloatArray, IntArray, PathInput, VlateralGeneratorSignature, WriteDischargesFn
from ._Config import Configs
from ._kernel_registry import dispatch

__all__ = [
    "Router",
]


class Router:
    """
    Muskingum style river routing allowing
    - static or dynamic coefficients (e.g. muskingum vs muskingum-cunge style)
    - channel-only or lateral/external water inputs (forcing)
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
        if configs is not None and configs != "":
            configs = str(configs)
            if configs.endswith(".json"):
                with open(configs) as f:
                    raw = json.load(f)
            elif configs.endswith((".yml", ".yaml")):
                with open(configs) as f:
                    raw = yaml.load(f, Loader=yaml.FullLoader)
            else:
                raise RuntimeError("Unrecognized simulation config file type. Must be .json or .yaml")
        raw.update(kwargs)
        raw.pop("_router", None)
        self.cfg = Configs(**raw)

        # configure logging - progress bar and info/debug logs are mutually exclusive
        self.logger = logging.getLogger(f"river_route.{id(self):x}")
        self.logger.disabled = not self.cfg.log
        self.logger.setLevel(self.cfg.log_level)
        if self.cfg.log_stream == "stdout":
            self.logger.addHandler(logging.StreamHandler(sys.stdout))
        else:
            self.logger.addHandler(logging.FileHandler(self.cfg.log_stream))
        self.logger.handlers[0].setFormatter(logging.Formatter(self.cfg.log_format))
        self.logger.debug("Logger initialized")

        # default discharge writer; overridable via set_write_discharges
        self._write_discharges = self._default_write_discharges
        return

    def __repr__(self) -> str:
        return f"{type(self).__name__}(params_file={self.cfg.params_file!r})"

    ################################################
    # Initial and final state handling
    ################################################

    def _read_initial_state(self) -> None:
        if hasattr(self, "channel_state"):
            return

        state_file = self.cfg.channel_state_init_file
        if not state_file:
            self.logger.warning("channel_state_init_file not provided. Defaulting to zero initial conditions")
            self.channel_state = np.zeros(self.river_ids.shape[0], dtype=np.float32)
            return
        self.logger.debug("Reading Initial State from Parquet")
        self.channel_state = pd.read_parquet(state_file).values.flatten().astype(np.float32, copy=False)
        return

    def _write_final_state(self) -> None:
        final_state_file = self.cfg.channel_state_final_file
        if not final_state_file:
            return
        self.logger.debug("Writing Final State to Parquet")
        pd.DataFrame({"Q": self.channel_state}).to_parquet(self.cfg.channel_state_final_file)
        return

    ################################################
    # Prepare arrays for routing
    ################################################

    # todo parameterize all column names: missing downstream_river_id
    def _set_routing_vectors(self):
        """Build the time-independent network vectors once, before routing. Time options and the (static)
        Muskingum coefficients depend on the routing/runoff timestep, so they are set per routing pass in the
        _execute_routing_* methods, not here."""
        self._set_vectors_from_params()  # river_ids, next_river_ids, k, x, alpha/beta (dynamic only)
        self._set_connectivity_vectors()  # downstream_indices

    def _set_vectors_from_params(self):
        self.logger.debug("Calculating network dependent vectors")
        # todo use a cache for the params file to speed up possible repeated reads
        df = pd.read_parquet(self.cfg.params_file)
        # todo validate that the expected columns for river, downstream river, k, and x are here before assigning

        if df[self.cfg.var_river_id].duplicated().any():
            raise ValueError("params_file contains duplicate river IDs.")

        self.river_ids = np.ascontiguousarray(df[self.cfg.var_river_id].to_numpy(copy=False), dtype=np.int64)
        self.next_river_ids = np.ascontiguousarray(df["downstream_river_id"].to_numpy(copy=False), dtype=np.int64)
        self.k = np.ascontiguousarray(df["k"].to_numpy(copy=False), dtype=np.float32)
        self.x = np.ascontiguousarray(df["x"].to_numpy(copy=False), dtype=np.float32)

        # dynamic coefficients additionally need the alpha and beta columns
        if self.cfg.coeff == "dynamic":
            self.alpha = np.ascontiguousarray(df["alpha"].to_numpy(copy=False), dtype=np.float32)
            self.beta = np.ascontiguousarray(df["beta"].to_numpy(copy=False), dtype=np.float32)
            if np.any(self.alpha <= 0):
                raise ValueError("alpha column in params_file must be strictly positive")
        return

    def _validate_time_options(self):
        # check that time options have the correct relative sizes
        if self.dt_total < self.dt_runoff:
            raise ValueError("dt_total must be >= dt_runoff")
        if self.dt_total < self.dt_discharge:
            raise ValueError("dt_total must be >= dt_discharge")
        if self.dt_discharge < self.dt_runoff:
            raise ValueError("dt_discharge must be >= dt_runoff")
        if self.dt_runoff < self.dt_routing:
            raise ValueError("dt_runoff must be >= dt_routing")
        if not (self.dt_total >= self.dt_discharge >= self.dt_runoff >= self.dt_routing):
            raise ValueError("Need dt_total >= dt_discharge >= dt_runoff >= dt_routing")
        # check that time options are evenly divisible
        if self.dt_total % self.dt_runoff != 0:
            raise ValueError("dt_total must be an integer multiple of dt_runoff")
        if self.dt_total % self.dt_discharge != 0:
            raise ValueError("dt_total must be an integer multiple of dt_discharge")
        if self.dt_discharge % self.dt_runoff != 0:
            raise ValueError("dt_discharge must be an integer multiple of dt_runoff")
        if self.dt_runoff % self.dt_routing != 0:
            raise ValueError("dt_runoff must be an integer multiple of dt_routing")

        # Now that we know time parameters are valid, set time-derived parameters for computation cycles
        self.num_runoff_steps = int(self.dt_total / self.dt_runoff)
        self.num_runoff_steps_per_discharge = int(self.dt_discharge / self.dt_runoff)  # to resample/reshape results
        # todo for synthesized networks, this will have to become a vector
        self.num_routing_steps_per_runoff = int(self.dt_runoff / self.dt_routing)
        return

    def _set_forced_time_options(self, dates: DatetimeArray) -> None:
        self.dt_runoff = self.cfg.dt_runoff or (dates[1] - dates[0]).astype("timedelta64[s]").astype(int)
        self.dt_discharge = self.cfg.dt_discharge or self.dt_runoff
        self.dt_total = self.cfg.dt_total or self.dt_runoff * dates.shape[0]
        if not self.cfg.dt_routing:
            self.logger.warning("dt_routing was not provided or is Null/False, defaulting to dt_runoff")
        self.dt_routing = self.cfg.dt_routing or self.dt_runoff
        self._validate_time_options()
        return

    def _set_channel_time_options(self) -> None:
        self.dt_routing = self.cfg.dt_routing
        self.dt_total = self.cfg.dt_total
        self.dt_discharge = self.cfg.dt_discharge or self.dt_routing
        self.dt_runoff = self.dt_discharge  # Assign to pass time validation. It is never used in channel routing
        self._validate_time_options()
        return

    def _set_static_muskingum_coefficients(self):
        """
        implied dependency on having set time options and the network parameter vectors
        """
        self.logger.debug("Calculating Muskingum coefficients")
        dt_div_k = self.dt_routing / self.k
        denominator = dt_div_k + (2 * (1 - self.x))
        _2x = 2 * self.x
        # contiguous arrays iterate faster in kernels due to cpu and ram access patterns
        self.c1 = np.ascontiguousarray((dt_div_k - _2x) / denominator, dtype=np.float32)
        self.c2 = np.ascontiguousarray((dt_div_k + _2x) / denominator, dtype=np.float32)
        self.c3 = np.ascontiguousarray(((2 * (1 - self.x)) - dt_div_k) / denominator, dtype=np.float32)
        self.c4 = np.ascontiguousarray(self.c1 + self.c2, dtype=np.float32)
        self.c4_dt = np.ascontiguousarray(self.c4 / self.dt_runoff, dtype=np.float32)
        if not np.allclose(self.c1 + self.c2 + self.c3, 1):
            self.logger.warning("Muskingum coefficients do not sum to 1")
            self.logger.debug(f"c1: {self.c1}")
            self.logger.debug(f"c2: {self.c2}")
            self.logger.debug(f"c3: {self.c3}")
            raise ValueError("Muskingum coefficients do not sum to 1, check routing parameters and time step")

        # shuffling arrays to list coefficient of the downstream increases performance of kernel which can
        # read sequentially when solving each river, rather than essentially randomly throughout the array
        self.downstream_c1 = np.zeros(self.downstream_indices.shape[0], dtype=np.float32)
        self.downstream_c2 = np.zeros(self.downstream_indices.shape[0], dtype=np.float32)
        valid = self.downstream_indices >= 0
        self.downstream_c1[valid] = self.c1[self.downstream_indices[valid]]
        self.downstream_c2[valid] = self.c2[self.downstream_indices[valid]]
        return

    def _set_connectivity_vectors(self):
        """Build the index vectors describing network connectivity from river_ids and next_river_ids.

        downstream_indices: index of each segment's downstream segment in the parameter arrays, -1 at outlets.
        Requires the params file to be topologically sorted (upstream before downstream).
        """
        self.logger.debug("Calculating network connectivity vectors")
        n = self.river_ids.shape[0]
        river_index = {int(river_id): idx for idx, river_id in enumerate(self.river_ids.tolist())}

        # 1D array giving the index of the downstream river in the parameter arrays, -1 if none downstream
        self.downstream_indices = np.full(n, -1, dtype=np.int32)
        for upstream_idx, downstream_river_id in enumerate(self.next_river_ids.tolist()):
            if downstream_river_id < 0:
                continue
            downstream_idx = river_index.get(int(downstream_river_id))
            if downstream_idx is None:
                raise ValueError(f"params_file downstream_river_id {downstream_river_id} is not in the river_id column")
            if downstream_idx <= upstream_idx:
                raise ValueError("params_file must be topologically sorted upstream to downstream")
            self.downstream_indices[upstream_idx] = downstream_idx

        self.logger.log(logging.INFO, f"Network: {n} river segments")
        return

    #################################################
    # Generator for lateral inflow routing
    #################################################

    def _vlateral_generator(self) -> VlateralGeneratorSignature:
        if self.cfg.qlateral_files:
            for lateral_file, discharge_file in zip(self.cfg.qlateral_files, self.cfg.discharge_files, strict=True):
                self.logger.info("-" * 60)
                with xr.open_dataset(lateral_file) as ds:
                    dates = ds["time"].values.astype("datetime64[s]")
                    array = ds["qlateral"].values.astype(np.float32, copy=False)
                    yield dates, array, lateral_file, discharge_file
        elif self.cfg.grid_runoff_files and self.cfg.grid_weights_file:
            for runoff_file, discharge_file in zip(self.cfg.grid_runoff_files, self.cfg.discharge_files, strict=True):
                self.logger.info("-" * 60)
                self.logger.debug(f"Calculating qlateral: {runoff_file}")
                ds = runoff_to_qlateral(
                    runoff_file,
                    grid_weights_file=self.cfg.grid_weights_file,
                    var_runoff=self.cfg.var_grid_runoff,
                    var_x=self.cfg.var_x,
                    var_y=self.cfg.var_y,
                    var_t=self.cfg.var_t,
                    var_river_id=self.cfg.var_river_id,
                    cumulative=self.cfg.grid_accumulation_type == "cumulative",
                    as_volumes=True,
                )
                yield (
                    ds["time"].values.astype("datetime64[s]"),
                    ds["qlateral"].values.astype(np.float32, copy=False),
                    runoff_file,
                    discharge_file,
                )

    ################################################
    # Methods to execute routing simulation
    ################################################

    def route(self) -> Self:
        """
        Execute the simulation described by the provided configs and routing parameters. All configs, file paths,
        parameters, and options must be set when the object is initialized so that validation is performed before the
        simulation.

        Prep and validation of time, network, etc vectors is handled by tailored subroutines. Calls to set muskingum
        coefficients, time options, and so forth, do not appear as calls dispatched at the top level routing method

        Returns:
            Self: the class instance with updated channel_state and output files written to disk
        """
        # start timer
        self.logger.log(PROGRESS, "Beginning routing")
        t1 = datetime.datetime.now()
        # validate configuration options
        self.logger.debug("Validating configs")
        self.cfg.validate()
        self.logger.debug(self)
        # build network vectors, read state, route, write state
        self._set_routing_vectors()
        self._read_initial_state()
        self._execute_routing()
        self._write_final_state()
        # log total time
        t2 = datetime.datetime.now()
        self.logger.log(PROGRESS, f"Routing completed in {(t2 - t1).total_seconds()} seconds")
        return self

    def _execute_routing(self) -> None:
        # there are two types of loops, one for channel only, one if external forcing(s) are provided.
        if self.cfg.forcing == "channel":
            self._execute_routing_channel()
        else:
            self._execute_routing_lateral()
        return

    def _execute_routing_channel(self) -> None:
        self.logger.info("-" * 60)

        # channel routing time options and (static) coefficients
        self._set_channel_time_options()
        self._set_static_muskingum_coefficients()

        self.logger.debug("Starting routing computation")
        q_t = self.channel_state.astype(np.float32, copy=True)
        discharge_array = np.zeros((self.num_runoff_steps, self.river_ids.shape[0]), dtype=np.float32)
        dispatch(self, q_t=q_t, discharge_array=discharge_array)
        self.channel_state = q_t

        # Generate date array for output
        dates = pd.date_range(
            start=self.cfg.start_datetime,
            periods=self.num_runoff_steps,
            freq=pd.to_timedelta(self.dt_discharge, unit="s"),
        ).to_numpy()

        # write outputs
        self.logger.debug("Writing Discharge Array to File")
        discharge_array = discharge_array.astype(np.float32, copy=False)
        self._write_discharges(dates, discharge_array, self.cfg.discharge_files[0])
        self.logger.info("-" * 60)
        return

    def _execute_routing_lateral(self) -> None:
        """Lateral runoff forcing: loop over vlateral files"""
        self._ensemble_member_states = []

        total_files = len(self.cfg.qlateral_files or self.cfg.grid_runoff_files or [])
        file_iter = self._vlateral_generator()
        if self.cfg.progress_bar:
            file_iter = tqdm(file_iter, total=total_files, desc="Files Routed")

        coeff_dt: tuple[int, int] | None = None  # (dt_routing, dt_runoff) the static coefficients were built for
        for dates, qlateral, runoff_file, discharge_file in file_iter:
            self.logger.info(f"Routing qlateral: {runoff_file}")
            self._set_forced_time_options(dates)
            # static coefficients depend only on dt_routing/dt_runoff, so rebuild them only when those change
            # across files (dynamic coefficients are rebuilt inside the kernel each substep)
            if self.cfg.coeff != "dynamic" and (self.dt_routing, self.dt_runoff) != coeff_dt:
                self._set_static_muskingum_coefficients()
                coeff_dt = (self.dt_routing, self.dt_runoff)
            self.logger.debug("Starting routing computation")
            q_t = self.channel_state.astype(np.float32, copy=True)
            q_array = np.zeros((self.num_runoff_steps, self.river_ids.shape[0]), dtype=np.float32)
            dispatch(self, q_t=q_t, discharge_array=q_array, vlateral=np.ascontiguousarray(qlateral, dtype=np.float32))
            if self.cfg.runoff_processing_mode == "sequential":
                self.logger.debug("Updating Channel State for Next Sequential Computation")
                self.channel_state = q_t
            elif self.cfg.runoff_processing_mode == "ensemble":
                self.logger.debug("Recording Member State for Final State Aggregation")
                self._ensemble_member_states.append(q_t.copy())

            if self.dt_discharge > self.dt_runoff:
                self.logger.debug("Resampling dates and discharges to specified timestep")
                q_array = q_array.reshape(
                    (
                        int(self.dt_total / self.dt_discharge),
                        int(self.dt_discharge / self.dt_runoff),
                        self.river_ids.shape[0],
                    )
                ).mean(axis=1)
                dates = dates[:: self.num_runoff_steps_per_discharge]

            self.logger.debug("Writing Discharge Array to File")
            q_array = q_array.astype(np.float32, copy=False)
            self._write_discharges(dates, q_array, discharge_file, runoff_file)

        if self.cfg.runoff_processing_mode == "ensemble":
            self.channel_state = np.array(self._ensemble_member_states).mean(axis=0)
        self.logger.info("-" * 60)
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

    def _default_write_discharges(
        self,
        dates: DatetimeArray,
        q_array: FloatArray,
        q_file: PathInput,
        routed_file: PathInput = "",
    ) -> None:
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
        with nc.Dataset(str(q_file), mode="w", format="NETCDF4") as ds:
            ds.createDimension("time", size=q_array.shape[0])
            ds.createDimension(self.cfg.var_river_id, size=q_array.shape[1])
            ds.runoff_file = str(routed_file)
            time_var = ds.createVariable("time", "f8", ("time",))
            time_var.units = f"seconds since {pd.Timestamp(dates[0]).strftime('%Y-%m-%d %H:%M:%S')}"
            time_var[:] = (dates - dates[0]).astype("timedelta64[s]").astype(np.int64)
            id_var = ds.createVariable(
                self.cfg.var_river_id,
                "i4",
                self.cfg.var_river_id,
            )
            id_var[:] = self.river_ids
            flow_var = ds.createVariable(self.cfg.var_discharge, "f4", ("time", self.cfg.var_river_id))
            flow_var[:] = q_array
            flow_var.long_name = "Discharge at catchment outlet"
            flow_var.standard_name = "discharge"
            flow_var.aggregation_method = "mean"
            flow_var.units = "m3 s-1"
        return
