## Changelog

---

### [v3.0.1](https://github.com/rileyhales/river-route/tree/v3.0.1) — 2026-10-13

- The npm package's WebAssembly is compiled by `python river_route/kernels/build.py --wasm <folder>`, which needs
  nothing beyond Python's standard library. `rr compile --wasm <folder>` still compiles it too.
- `river_route.kernels.build.build(llvm_bin, checks, wasm)` compiles either the shared library or, given `wasm`, the
  WebAssembly modules. `build_wasm` is removed, and `build` no longer returns the library's path.
- npm releases are staged, then approved, before they go live.

---

### [v3.0.0](https://github.com/rileyhales/river-route/tree/v3.0.0) — 2026-10-13

The [v2 to v3 migration guide](migrating/v2-to-v3.md) gives the replacement for each removed class, config key, and
file format.

**Routing**

- One class, `rr.Router`, routes every forcing in place of the `Muskingum`, `RapidMuskingum`, and
  `UnitMuskingum` classes, with static Muskingum coefficients built from `muskingumK` and `muskingumX`. What it routes
  is chosen by the `forcing` config (`channel`, `catchment`, `grid`, or `ecmwf_grib`), and how each river computes by
  the compute options `dt_routing` (`'auto'` or a step in seconds), `use_subreaches`, and `use_subcycles` (`True`
  or `False`).
- Rivers are routed one at a time from upstream to downstream, each river's whole time series before the next river,
  instead of every river for one time step before the next step by sparse forward substitution. A river routed whole
  with static coefficients is routed in one loop over eight steps at a time, held in two SIMD vectors of four lanes,
  with the recurrence solved as prefix sums within each vector. It routed the Amazon 1.55 times faster than taking each
  river through its upstream terms, runoff, recurrence, downstream inflow, and discharge in loops of their own.
- Two routing kernels route the same way: a C kernel, `river_route/kernels/c`, and its numba twin. `rr compile`
  compiles the C for this machine with the clang of LLVM, and `--checks` also requires it to compile without a
  warning from clang's strictest set and to pass clang-tidy and clang's static analyzer. `Router.route(kernel='auto')`, the default, routes
  with the C once it is compiled and with numba until then, and `kernel='c'` or `kernel='numba'` chooses one.
- The C kernel is compiled to WebAssembly by `rr compile --wasm <folder>` and published to npm as `river-route`, with
  one module for runtimes with relaxed SIMD's fused multiply-adds and one for those without.
- Every river hands its discharge downstream, and records it, at the end of each runoff step, and computes within it
  in a plan of its own, settled in two jobs. `dt_routing` sets each river's step: `'auto'`, the default, gives each
  river its own valid step, the runoff step over the fewest subcycles that keep its coefficients non-negative, and a
  step in seconds, a divisor of the runoff step, is the step of every river. `use_subreaches` and `use_subcycles`, on
  by default, then make every coefficient non-negative at that step, and `'auto'` requires both. Without
  them every river is routed whole at its step, and a river with a negative coefficient at it is reported. The subcycles of each sub-reach are combined
  into one step of coefficients that give the discharge of taking them one at a time, built over the bits of the
  count of subcycles, so they cost one step however many there are. `Network.largest_valid_dt` is gone, since no step
  is searched for.
- `use_subreaches` and `use_subcycles`, the defaults, route each river with non-negative coefficients in its own
  step, the condition for valid Muskingum routing, where they can: with `use_subreaches` a river too long for its
  step is split into sub-reaches, the fewest
  sub-reaches in series inside the kernel, laid out by `Network.conditioning` to keep the river's mean travel
  time and travel time variance with x never negative (equal sub-reaches with `x = 1/2 - N*(1/2 - x)` where that is not
  negative, otherwise linear reservoirs of unequal travel times, 0.724 and 0.276 of k at `x = 0.2`), and with
  `use_subcycles` a river too short for it is routed in subcycles, shorter steps of its own. With both, a river that
  neither makes valid, which only a river with x above 1/3 can be, is split into the fewest sub-reaches of unequal
  travel times that each take subcycles of their own and keep the river's diffusion, `x_i = 1/2 - (1/2 - x) k / k_i`,
  found by `river_route.network.splitting`. A river at `x = 0.49` with a step of 0.6 of its k splits into three. They
  route at any step, with every forcing, and with `use_subreaches` state files hold one row per sub-reach.
- `Network.stabilize(dt)` rewrites a network in place with the sub-reaches `Network.conditioning(dt)` splits rivers
  into as rows of their own, laid out as `use_subreaches` and `use_subcycles` route them, and
  `Network.write_stabilized(dt)` writes it as a network table, in which added reaches are numbered
  down from -1,000,000, `synthetic` marks them, and `parentRiverId` maps each reach to the river it was split from.
  Such a network is routed with the `channel`, `grid`, or `ecmwf_grib` forcing, whose RunoffForcing classes share each
  river's runoff among its sub-reaches with `RunoffForcing.distribute`.
- Routing with `use_subreaches` or `use_subcycles`, and `Network.stabilize`, refuse a river with `muskingumX` above
  0.49, `river_route.network.MAX_STABILIZED_X`. A river needs up to `ceil(1/(1-2x))` sub-reaches, 50 at 0.49 and
  many more above it, and at x = 0.5 it is valid only at a step equal to its travel time, so that it can be
  impossible to stabilize. Rivers above it are routed only whole, with `dt_routing` in seconds and both options off.
- Rivers with a negative coefficient in their own step, outside `2*k*x <= step <= 2*k*(1-x)`, are reported as
  errors. The new `invalid_coefficients` config chooses `raise` (the default), which refuses to route them, `log`,
  which logs the error at the error level and routes on, or `ignore`.
- Negative discharge is written as routed. It was clamped to zero in the written discharge, which hid the rivers that
  produced it and broke the volume of their written discharge.
- The discharge written for each output interval is the discharge at its end, after routing it: that of the last
  subcycle of its last runoff step. It was the mean of the discharges computed at
  its routing steps. Channel routing writes the discharge at the end of each `dt_discharge` rather than their mean.
- A Network refuses a river with `muskingumK` not positive and finite, such as a river of zero length, which no step
  can route and which is resolved in the topology before routing; a river with `muskingumX` below 0, or above 0.5, at
  which the scheme amplifies its inflow; and a river that drains into itself.
- `Router.route(thread_pool=pool, threads=n)` routes concurrently on a `ThreadPoolExecutor` the caller passes. The
  region's rivers are divided into blocks, contiguous index ranges of rivers, and the blocks are packed into one job
  per thread. The pool is used as given and never shut down, and the package never creates one: without a pool,
  routing is single threaded. The discharge writer is handed the same pool and thread count, and `zarr_writer` writes
  its chunks on that pool.
- The unit hydrograph transform is postponed past v3: `UnitMuskingum`, `river_route.uhkernels`, and the
  `uh_kernel_file`, `uh_state_init_file`, and `uh_state_final_file` configs are removed.
- Repeated `route()` calls on one Router restart from `channel_state_init_file` instead of continuing from the
  previous run's final state.
- Every module logs to a child of the `river_route` logger, and a Router points that logger at its log options with
  one handler, which replaces the handler of the Router before it. Handlers no longer accumulate, and the runoff and
  weight table functions follow the log options too.

**Network**

- `rr.Network` owns the river network: the ids, topology, and parameters of the network file, the partition of the
  rivers into blocks, the validity analysis, and stabilization. A Router takes only a Configs and builds its
  Network from it, and reuses it, with its partition, for every `route()` call. A Network is built from a network
  file or a DataFrame, is stabilized once, since `stabilize` refuses a network that already is, and is written with
  `to_parquet(path)`. `Router.river_ids`, `k`, and `x` are read off the network, as `router.network.river_ids`.
- A Network holds `muskingumK` and `muskingumX` in float64, whatever the network file stores, and writes them in
  float64, so the travel times and x of a stabilized network are written exactly as they were split.
- The params file is renamed the network file, set with the `network_file` config in place of `params_file`: it
  describes the network the `Network` class is built from, not only its Muskingum parameters. Its columns are
  camelCase: `riverId`, `nextRiverId` (was `downstream_river_id`), `muskingumK` (was `k`), and `muskingumX` (was
  `x`), with `synthetic` and `parentRiverId` on a stabilized network.
- The rivers of the network file must be in depth first search (DFS) order, where v2 required only upstream before
  downstream: the rivers upstream of each river are the rows immediately before it, so every river's upstream
  watershed is one contiguous range of rows. Two new required columns describe that order: `riverIndex`, each river's
  position in it, and `upstreamCount`, the number of rivers upstream of each river. Every file with one entry per
  river lists the rivers in that order, and routing refuses catchment runoff or a grid weight table whose rivers are
  not the rivers of the network file in that order. `examples/migrate_v2_to_v3.py` converts v2 inputs, writing the
  params file as a network file in DFS order with both columns, and the grid weights, qlateral files, and channel
  state files in the same order.
- `Network` analyzes the validity of Muskingum coefficients with `validity_window`, `invalid_mask`, `sub_reaches`,
  `subcycles`, `conditioning`, which gives the `Layout` of each river's sub-reaches, and `check_validity`, refuses
  with `check_routable`, before any runoff is read, a
  network the kernels cannot route, such as a river split into more than 256 sub-reaches, and divides the rivers
  into the blocks threads route with `routing_blocks`, sized for the thread count. Each block is read from `upstreamCount`, so a
  `Network` refuses a table whose `upstreamCount` does not describe its DFS order.
- `river_route.tools` is removed. `subset_network_to_river` is in `river_route.network.streams`, with
  `assign_blocks` and `analyze_partitioning`, which partition networks, and cuts a basin as the rows its outlet's
  `riverIndex` and `upstreamCount` give, without building a graph. `connectivity_to_digraph` and `adjacency_matrix`
  are removed.

**Runoff**

- `forcing_files` replaces `qlateral_files` and `grid_runoff_files`, and is read in the form `forcing` names.
  `grid_weights_file` is required for the `grid` and `ecmwf_grib` forcings. Runoff options a forcing does not read
  are refused unless left at their defaults: `catchment` reads none, `grid_weights_file` included.
- `var_grid_runoff` is renamed `var_runoff`, and `grid_accumulation_type` is renamed `accumulation_type`.
- Catchment runoff files replace qlateral files. A file holds a `catchment_runoff` variable with dimensions
  `(riverId, time)`, in that order, and a `catchment_area` variable (m²) with dimension `riverId`, linked by the CF
  attribute `cell_measures = "area: catchment_area"`. The `units` attribute of `catchment_runoff` is required and
  marks it as volumes (`m3`) or depths (`m`, `mm`), which are multiplied by the catchment area when read.
  `RunoffForcing.to_netcdf` writes this schema.
- The runoff classes `CatchmentRunoffForcing`, `GridRunoff`, and `ECMWFGribReducedGrid` read the runoff of each forcing.
  A Router builds the one its `forcing` names the first time it routes runoff. `GridRunoff` reads and indexes its
  weight table once for every forcing file, and routing reads each river's grid cells directly into its forcing as the
  river is routed. `GridRunoff.to_dataset` and `GridRunoff.aggregate_to_file` replace `runoff.runoff_to_qlateral`.
- `ECMWFGribReducedGrid` routes GRIB files on a reduced gaussian grid, such as the octahedral O1280 grid of the ECMWF
  IFS, reading each message with eccodes and keeping only the weight table's cells. At the 18,721 cells of the
  Columbia in a 145 step O1280 forecast that reads in 0.72 s with memory for one message, where xarray and cfgrib
  took 5.8 s and 4 GB. `ReducedGaussianGrid.from_grib` reads a file's grid from its metadata,
  `ReducedGaussianGrid.cell_polygons` gives the area each cell represents, and `reduced_grid_weights` builds a weight
  table that locates each cell by its `cell_index`. IFS forecast runoff accumulates from the start of the forecast, so
  `ecmwf_grib` runoff is always cumulative and has no `accumulation_type` option.
- `runoff_depth_unit`, `force_positive_runoff`, and `as_volumes`, arguments of `runoff_to_qlateral` in v2, are
  configs.
- Runoff must come in uniform time steps. `force_uniform_timesteps`, which resampled irregular gridded runoff in v2,
  is removed: routing refuses a forcing file whose time steps are not all the same, and preparing runoff with a
  uniform time step is left to its provider. `dt_runoff` is no longer a config; it is read from the time steps of
  each forcing file.
- NaN runoff is set to zero before it is routed: a NaN cell contributes nothing while the other cells of its
  catchment still count.
- A runoff `units` attribute of `kg m-2` is read as millimeters of water; v2 read it as meters.
- The weight table functions `grid_weights`, `compute_voronoi_catchment_intersects`,
  `voronoi_diagram_from_regular_xy`, and `cell_xy_from_regular_grid` are in `river_route.forcings.weights`, and
  still importable from `river_route.forcings`, the package named after the `forcing` config that replaces
  `river_route.runoff`.

**Configs**

- `Configs` is frozen and keyword only. Build one from keyword arguments or with `Configs.from_json`, write one with
  `Configs.to_json`, and pass it to `Router(configs)`. Routers no longer take a config file path or options as
  keyword arguments.
- Config files are JSON only. YAML config files are no longer read, and `pyyaml` is no longer a dependency.
  `examples/config.json` lists every option.
- Unrecognized keys in a config file raise a `ValueError` naming them, and unknown keyword arguments a `TypeError`.
  `Router.route` validates the options and paths with `Configs.validate_routing`, and `GridRunoff.from_configs` with
  `Configs.validate_runoff`; passing one does not pass the other. `Configs.deep_validate()` reads the network file,
  grid weights, and initial state that are set and checks their contents and their consistency with each other;
  nothing calls it for you.
- Outputs named from `discharge_dir` are `discharge_<input name>.zarr`, taking the extension of the store the
  default writer writes instead of the extension of the forcing file. Input files with duplicate names are rejected
  rather than overwriting each other's output. `discharge_dir=os.devnull` routes and discards the discharge.
- `dt_discharge` and `dt_total` configure channel routing only. Routing with forcing routes every step of each forcing
  file and writes the discharge at the end of each, and refuses either option when it is set. Runoff is checked
  against the network file for its rivers, in order, before it is routed.
- `dt_discharge` and `dt_total` must be whole numbers of seconds, or 0 to leave them unset, and `dt_routing` must be
  `'auto'` or a positive whole number of seconds; any other value, and `dt_routing: auto` without `use_subreaches` and `use_subcycles`, is
  refused when the `Configs` is built. Runoff whose time steps do not increase is refused.

**Outputs**

- `river_route.router.writers` holds the discharge writers. `zarr_writer` is the default: a zarr store of `Q` with
  dimensions `(riverId, time)`, chunked so each chunk holds every time step of 500 rivers, rounded to
  `writers.ZARR_KEEPBITS` mantissa bits (a relative error of at most 2^-13) and compressed with Blosc lz4 and
  bitshuffle, 2.71x smaller than the raw array on a year of the Amazon. `netcdf_writer` writes an uncompressed
  netCDF file with dimensions `(riverId, time)`, where v2 wrote `(time, river_id)`. `null_writer` writes nothing.
  The discharge variable is always `Q`: the `var_discharge` config is removed, and a custom writer names it otherwise.
- Every file keyed by river names its river id `riverId`, as the network file does, where v2 named it `river_id`:
  channel state files, catchment runoff files, grid weight tables, catchment polygons read by the weight table
  functions, and discharge. The name is not configurable: the `var_river_id` config is removed, and a file that names
  it otherwise is renamed before it is read. The `var_` options that remain only name the variables of the gridded
  forcing files that are read; every file river-route writes uses the default names.
- `set_write_discharges` is renamed `set_discharge_writer`. A writer is called as
  `writer(router, dates, discharge_array, discharge_file, forcing_file, thread_pool=..., threads=...)`: it is handed
  the Router first, the discharge as a C-order `(river, time)` array, `forcing_file` `''` for channel routing, and the
  thread pool and thread count `Router.route` was given. `router.network.original_river_ids` is the river of each row
  of the discharge array, which on a stabilized network leaves out the synthetic sub-reaches.
- Channel state files hold a `riverId` column beside `Q`. Final state files are written with it, and an initial state
  file must have it: routing refuses one whose `riverId` does not list the rivers of the network file in order.
  `examples/migrate_v2_to_v3.py` adds it to v2 state files. River ids may be any integer. They are written as int64 in
  discharge, catchment runoff files, and network tables written by `Network`, and read as int64 from grid weight tables.

**Metrics**

- `metrics.kge2012` computes its variability ratio as the simulated over the observed coefficient of variation; v2
  inverted it.
- `metrics.mean_error` (`me`) is the simulated minus the observed values, positive when the simulation is too high,
  as in HydroErr and hydroGOF; v2 subtracted the simulated from the observed. Every metric returns a float.

**CLI**

- `rr route config.json` routes a config file, replacing `rr route --router <class> config` and the per-class
  commands `rr Muskingum`, `rr RapidMuskingum`, and `rr UnitMuskingum`.
- `rr subset <riverId> --network <in> --out-network <out>` cuts a network file, and optionally a grid weight table with
  `--weights` and `--out-weights`, to a river and every river upstream of it.
- `rr route --threads <n>` routes on a thread pool of n threads, and `rr route -` reads the config as JSON from stdin.
- `rr compile` compiles the C routing kernel, and `rr app` opens a window that sets up, saves, and runs routing
  configs.

**Packaging**

- Increased the minimum Python version to 3.14.
- Added lower version bounds to every dependency, and major version upper bounds to those that follow semantic
  versioning: `dask` and `xarray`, versioned by date, and `pyarrow` have none. Added `zarr`, `eccodes`, and
  `pyproj`, and removed `pyyaml` and `networkx`.
- river-route runs on macOS and Linux only. Importing it on any other platform, such as Windows, raises an
  ImportError.

---

### [v2.1.1](https://github.com/rileyhales/river-route/tree/v2.1.1) — 2026-04-11

- Small documentation edits.
- Update build system to hatchling and pyproject.toml

---

### [v2.1.0](https://github.com/rileyhales/river-route/tree/v2.1.0) — 2026-04-11

- Routing state, parameters (`k`, `x`), channel state, and all numba kernel scratch buffers now use `float32` instead of `float64`. Reduces memory by half, speeds up the inner routing loop.
- Fused the sparse matrix-vector multiply and forward-substitution into a single loop for modest speed gain.
- Better use of type aliases to reduce total number of annotations and improve readability.

---

### [v2.0.1](https://github.com/rileyhales/river-route/tree/v2.0.1) — 2026-03-10

- Adds pyarrow to dependencies which is included by conda installs as a pandas dependency but not in pip.

### [v2.0.0](https://github.com/rileyhales/river-route/tree/v2.0.0) — 2026-03-09

- Replaced `Muskingum` class with 3 separate classes `Muskingum`, `RapidMuskingum`, and `UnitMuskingum`.
- New class `Muskingum` is channel routing only with no runoff-transformation.
- New class `RapidMuskingum` is a reimplementation of the previous Muskingum class with routing
  and runoff transformation.
- New class `UnitMuskingum` is a new implementation of unit hydrograph runoff transformation and channel routing.
- Introduced `Configs` dataclass replacing the untyped config dictionary to centralize and validate configs.
- Merged `connectivity_file` into `params_file` (single-file network definition).
- Changed grid weights format from CSV to netCDF with `proportion` column.
- Simplified channel state files from two columns (Q, R) to one column (Q).
- Added topological sort validation on routing parameters.
- Added `types.py` and significantly improved type annotations and static checking coverage.
- Expanded `runoff.py` functions for Voronoi diagrams and grid weight computations.
- Renamed numerous config keys (see migration guide).
- Removed conversion utilities for RAPID inputs.
- Several new tutorials, references, and examples in the docs.
- Routing now uses numba JIT-compiled forward substitution replacing the scipy sparse linear solver.
- Headwater streams are excluded from the matrix solve in `UnitMuskingum`, reducing the system size roughly in half.
- Progress tracking uses tqdm at the file iteration level for `RapidMuskingum` and `UnitMuskingum`.
- Increased minimum Python version to 3.12.
- Added dependency on `numba`.

---

### [v1.3.0](https://github.com/rileyhales/river-route/tree/v1.3.0) — 2025-08-06

- Added support for routing runoff ensembles.

### [v1.2.3](https://github.com/rileyhales/river-route/tree/v1.2.3) — 2025-05-01

- Removed volume calculation keyword argument.

### [v1.2.2](https://github.com/rileyhales/river-route/tree/v1.2.2) — 2025-04-23

- Fixed catchment volume calculation to use the correct area column.

### [v1.2.1](https://github.com/rileyhales/river-route/tree/v1.2.1) — 2025-03-14

- Fixed bug where the time variable name argument was not applied correctly.

### [v1.2.0](https://github.com/rileyhales/river-route/tree/v1.2.0) — 2025-02-27

- Improved efficiency of catchment volume computations.

### [v1.1.0](https://github.com/rileyhales/river-route/tree/v1.1.0) — 2025-01-21

- Switched to a direct solver for the Muskingum linear system.

### [v1.0.3](https://github.com/rileyhales/river-route/tree/v1.0.3) — 2025-01-17

- Renamed `_MuskingumCunge` to `_Muskingum`.
- Added `metrics.py` module.
- Refactored `runoff.py`.
- Various bug fixes and documentation updates.

### [v1.0.2](https://github.com/rileyhales/river-route/tree/v1.0.2) — 2024-09-24

- Guaranteed consistent sort order for adjacency matrix construction.

### [v1.0.1](https://github.com/rileyhales/river-route/tree/v1.0.1) — 2024-08-19

- Fixed `input_type` config value not being set correctly.

### [v1.0.0](https://github.com/rileyhales/river-route/tree/v1.0.0) — 2024-08-17

- Initial stable release.
- Code cleanup and documentation restructuring.
