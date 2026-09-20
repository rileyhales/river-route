## Changelog

---

### Unreleased

- `RunoffGaussianGrid` is a dataclass. Its options are fields declared with their defaults and documented where
  they are declared, so the constructor is a `__post_init__` that reads the weight table, and the call signature is
  unchanged. Each field is named the same as the `Configs` option it comes from, so `from_configs` reads every one by
  its own field name instead of listing them. `var_runoff` and `cumulative`, the names the `Runoff` base class uses,
  are now properties derived from the `var_grid_runoff` and `grid_accumulation_type` fields rather than copies made in
  the constructor. `rivers_per_block` is a class attribute, not a constructor argument, as before.
- Added the `var_vlateral` config, the name of the lateral inflow variable in `vlateral_files`, default `vlateral`.
  The name was read literally before, so a file that called it something else had to be rewritten to be routed.
  `Runoff.to_netcdf` still writes `vlateral`, which is the name the file schema documents.
- `RunoffVlateral` is a dataclass and names what it reads: `var_vlateral` and `var_t`, each defaulting to the name
  the schema documents. It had no constructor before and read both names literally, which left the `var_runoff` and
  `var_t` the `Runoff` base class declares unset on it. `var_runoff` is now a property over the `var_vlateral` field,
  the same way `RunoffGaussianGrid` derives it from `var_grid_runoff`. `RunoffVlateral.from_configs` reads each field
  off the `Configs` by its own name, and `Router` builds the reader for `vlateral_files` with it, or with a
  `RunoffVlateral` given as `Router(configs, runoff=...)`, instead of a bare `RunoffVlateral()` that ignored both.
  It does not call `Configs.validate_runoff`, which requires the `grid_weights_file` that routing from
  `vlateral_files` does not use.
- Added the `routing_order` config, `'river'` (the default) or `'time'` (the only order before). River order routes
  each river's whole time series before the next river, solving that series eight steps per serial operation, and is
  6 to 8 times faster than time order on one thread and on 12. It has static channel, static vlateral, and dynamic
  vlateral kernels, and a fused kernel that aggregates gridded runoff and routes it in one pass without building a
  vlateral array. It routes concurrently on a `thread_pool` over the same regions as time order, packing the regions
  into one pass per thread. It reads vlateral in either `(time, river)` layout or, when handed the transpose of a
  `(river, time)` array, one contiguous row per river. See the new Routing Kernels reference.
- A discharge writer is handed the routed discharge as a C-order `(river, time)` array, which is the layout the
  kernels write in place, instead of the `(time, river)` array of earlier versions. `zarr_writer` and
  `netcdf_writer` write it without transposing anything; `parquet_writer`, whose columns are time steps, transposes
  it with `writers.to_time_major`, which a custom writer can call for the same reason. A custom writer that indexed
  `discharge_array[time_index]` or built a frame with `index=dates` needs to transpose or be reindexed.
- `zarr_writer` is now the default discharge writer, and both array formats store `(river_id, time)`: zarr as
  chunks holding every time step of a block of rivers, netCDF as a variable dimensioned `(river_id, time)`. The
  river dimension comes first because that is the layout the river order kernels write in place. On a year of the
  Amazon this took a netCDF run from 7.87 s to 4.81 s with no change to the values, because the kernel no longer
  transposes each block of rivers out into a `(time, river)` array. Anything reading these files by name, such as
  `ds['Q'].transpose('time', 'river_id')`, is unaffected; anything assuming the axis order is not. `parquet_writer`
  is unchanged: parquet is columnar, so its columns stay time steps and its rows stay rivers.
- `zarr_writer` rounds discharge to `writers.ZARR_KEEPBITS` mantissa bits and compresses each chunk with
  `writers.ZARR_COMPRESSOR`, Blosc lz4 with bitshuffle, where it wrote the values raw and uncompressed before. The
  store is 2.71x smaller on a year of the Amazon and faster to write than the uncompressed one, since less of it
  reaches the disk, at a relative error of at most `2^-13` per value. The rounding is done in numpy a chunk of
  rivers at a time rather than by a zarr filter, so the array handed to the writer is never modified, and a
  `float16` run is stored as it is because float16 holds fewer mantissa bits than the rounding keeps.
- Added the `discharge_dtype` config, `'float32'` (the default) or `'float16'`, which narrows the discharge buffer
  in memory. The routing math stays float32 and the channel state is never narrowed, so only the saved copy is
  rounded, bounded by `2^-11` relative to each value. It needs `routing_order='river'` and cannot be combined with
  a `dt_discharge` coarser than `dt_runoff`. zarr and parquet store float16 natively; netCDF has no half type, so
  `netcdf_writer` widens to float32. float16 only covers 6.1e-5 to 65,504: on a year of the Amazon 0.07% of routed
  values overflow to infinity on the main stem, so `route` warns whenever the option is used.
- NaN runoff is set to zero once, when the cell series are prepared, instead of after each river's area weighted
  sum. A NaN cell now contributes nothing while the other cells of its catchment still count, where before the whole
  river got zero for that step. vlateral read from files has NaN set to zero too, so no routing kernel sees NaN. The
  `replace_nan` argument of the aggregation kernels is removed.
- `streams.shreve_order` and `Network.shreve_order()` give each river's Shreve magnitude, and
  `streams.assign_regions(..., measure='shreve')` claims concurrent regions by it instead of by river count. On a
  network where confluences join two rivers the two measures give the same partition.
- Added `river_route.Network`, which owns the river network: the ids, topology, and Muskingum parameters read from
  the params file, the connectivity vectors, and the concurrent routing partition. A `Router` builds one from its
  `Configs` and reuses it, so the params file is parsed and the network partitioned once per `Network` instead of
  once per `route()`. Assign `router.network` to route over an already built `Network`, which is how one parsed
  network backs many simulations. The `Router.river_ids`, `next_river_ids`, `downstream_indices`, `k`, `x`,
  `alpha`, and `beta` attributes are removed: read them off the network, as `router.network.river_ids`. A custom
  discharge writer that read `router.river_ids` needs the new spelling. `Router._set_vectors_from_params`,
  `_set_connectivity_vectors`, `_set_region_schedule`, `_check_coefficient_stability`, and `_check_river_alignment`
  are removed; the equivalents are `Network.routing_schedule` and `Network.check_stability`.
- `Network.stability_report(dt)` returns a `StabilityReport` counting how many rivers are Muskingum-stable at a
  routing time step, how many are too long or too short for it, and how much bigger a fixed network would be.
  Reports at the same dt add together, so a sweep over many parameter files accumulates into one total.
- `Network.stabilize(dt)` builds, in memory only, the stabilized network: every reach too long for `dt` is
  replaced by sub-reaches in series that each route stably at it, returned as the flat CSR arrays a kernel
  consumes. The network gains reaches rather than being divided up, hence `StabilizedNetwork`.
  `mode='uniform'` gives every sub-reach of a river the same travel time; `mode='nonuniform'` packs pieces of the
  largest stable travel time and leaves the remainder last; `weights=` apportions each river's travel time over an
  explicit sequence of segment lengths. Rivers that are too short for `dt` are a known gap: subdivision cannot fix
  them, `Network.substeps_required(dt)` reports the temporal refinement they would need, and no routing kernel
  consumes it yet.
- `Router(configs, network=None, runoff=None)` takes its options from the `Configs` and nothing else.
  `Router(configs, **overrides)` is removed; build the `Configs` you want and pass it, so there is one way to
  set every option. `Configs.replace` is removed too: a Configs is set once when it is built and there is no
  copy-with-changes.
- `Network` and `RunoffGaussianGrid` take the options they need as ordinary arguments and each has a `from_configs`
  classmethod that reads those same values off a `Configs` and builds the identical object. `Router` builds both
  with `from_configs`. `RunoffGaussianGrid(configs, **overrides)` is removed: `RunoffGaussianGrid.from_configs(configs)` replaces it, and
  `RunoffGaussianGrid(grid_weights_file, var_x=..., ...)` builds one without a `Configs` at all. `Configs.validate_runoff` now
  runs in `RunoffGaussianGrid.from_configs` rather than in the constructor, since a directly built `RunoffGaussianGrid` has no `Configs`
  to validate. The docs home page maps every config option to the class it builds.

- Increased the minimum Python version to 3.14.
- Increased the minimum pandas version to 3.0.5. pandas 3 changed `DataFrame.to_numpy()` to return a
  read-only array, which broke the in-place unit conversion when preparing vlateral from runoff with
  irregular timesteps. The 3.0.5 floor also skips 3.0.4, which is yanked from PyPI for segfaults in
  datetime handling that this package hit on every CF encoded time axis it read.
- Increased the minimum numpy version to 2.5.
- Added the abstract `Runoff` base class of `RunoffGaussianGrid` and `RunoffVlateral`. Its `to_netcdf` writes a
  vlateral array in the format `RunoffVlateral` reads and `vlateral_files` routes.
- Added `RunoffVlateral`, which reads `vlateral_files` for routing. Gridded runoff is read with `router.runoff`, so a
  `RunoffGaussianGrid` built or subclassed by hand and passed to `Router(configs, runoff=...)` is what prepares the
  inflow. The `RunoffVlateral.reader` and `RunoffGaussianGrid.reader` methods yield one
  `(dates, vlateral, source_file)` tuple per input and take only their input files, never a router or output paths.
- Added `river_route.router.writers` with premade discharge writers for `Router.set_discharge_writer`:
  `netcdf_writer`, the
  default, and `zarr_writer`, which writes `(time, river_id)` discharge in chunks of 500 rivers that
  each span every time step, writing up to the `threads` given to `Router.route` chunks at once (optionally packed
  into shard files with `writers.ZARR_CHUNKS_PER_SHARD`), and `parquet_writer`, which writes one row per
  river and one column per time step, uncompressed and without dictionary encoding or statistics by default
  (`writers.PARQUET_WRITE_OPTIONS`). `null_writer` sends discharge to the operating system's null device
  (`os.devnull`) so a route keeps nothing on disk. zarr is now a dependency.
- Discharge writers now receive the `Router` as their first argument:
  `writer(router, dates, discharge_array, discharge_file, runoff_file)`, so they can read
  `router.network.river_ids` and `router.configs`. Custom writers written for the previous 4 argument signature need the new first argument.
  `Router._default_write_discharges` is removed; use `river_route.router.writers.netcdf_writer`.
- `Router.set_write_discharges` is renamed `Router.set_discharge_writer`.
- `Configs` is its own subpackage, `river_route.configs`, and is the one way options are given. Build a frozen
  `Configs` from keyword arguments or with `Configs.from_file`, `from_json`, or `from_yaml`, then pass it to
  `Router(configs)` or `RunoffGaussianGrid(configs)`. `Router` no longer takes a config file path or options as keyword
  arguments, and `RunoffGaussianGrid` no longer takes a weight table path and options. `Configs.to_json` and
  `Configs.to_yaml` write the options to a file that `from_file` reads back. `Router.route` calls
  `Configs.validate_routing` and
  `RunoffGaussianGrid` calls `Configs.validate_runoff`, which check the options and the paths once. `Configs.validate` is
  removed. The runoff options `runoff_depth_unit`,
  `force_positive_runoff`, `force_uniform_timesteps`, and `as_volumes` are now configs.
- Routing from gridded runoff no longer rebuilds, copies, or reallocates lateral inflow for every runoff file.
  The weight table is read and checked against the params file once per `route()`, and each file is aggregated
  in a single pass (area weighting, de-accumulation, clipping, NaN replacement, and the volume product) straight
  into a reused C-order buffer that the kernels read without a copy. The runoff read itself is unchanged. With
  `threads` above 1 the rivers are split into ranges of similar work that the same kernel aggregates
  concurrently on the `thread_pool`; single threaded, one range covers every river.
- Threading happens only on a thread pool the caller passes; the package never creates one. Threads are a runtime
  resource, not a config. `Router.route(thread_pool=pool, threads=n)` uses a `ThreadPoolExecutor` as given and never
  shuts it down, so it can be shared with `RunoffGaussianGrid.vlateral(..., thread_pool=pool, threads=n)` and closed by the
  caller's `with` block. `threads` sets how many regions the network is partitioned into when a pool is given;
  without one, routing is single-threaded. `RunoffGaussianGrid.aggregate` and `RunoffGaussianGrid.to_dataset` take the same arguments.
  `Router.thread_pool()` is removed.
- `river_route.router.writers` holds the discharge writers, beside the class whose io it is. It is not importable
  as `river_route.writers` any more. The runoff readers are methods of `RunoffVlateral` and `RunoffGaussianGrid`.
- `river_route.runoff` is now a subpackage following the layout of `river_route.router`. Runoff preparation is
  the `RunoffGaussianGrid` class (also exported as `river_route.RunoffGaussianGrid`), which reads the weight table once
  and provides `read_runoff`, `aggregate`, `vlateral` (into a reused buffer), `reader`, and `to_dataset`. The
  class was named `Runoff` earlier in this release cycle and is renamed `RunoffGaussianGrid`. It replaces the
  `runoff_to_vlateral` function: use `RunoffGaussianGrid(Configs(grid_weights_file=..., ...)).to_dataset(runoff_files)`. The
  aggregation kernel lives in `river_route.runoff._numba_kernels`, and the weight table functions
  (`grid_weights`, `compute_voronoi_catchment_intersects`, `voronoi_diagram_from_regular_xy`,
  `cell_xy_from_regular_grid`) in `river_route.runoff.weights`, still importable from `river_route.runoff`.
- Fixed vlateral preparation raising `output array is read-only` for irregular timesteps with
  `as_volumes=True`.
- Added tests for the gridded runoff path: weight table construction, area weighted aggregation, unit
  conversion, cumulative de-accumulation, irregular timestep resampling, and routing from grid files.
- Input arrays are validated before they reach the numba kernels. The kernels are compiled without bounds
  checking, so a lateral inflow file or channel state file that does not match the params file previously
  read and wrote past the end of its array instead of raising.
- Lateral inflow files are checked against the params file for river count and river id order. A file whose
  columns are ordered differently is rejected rather than routing each river's water down the wrong reach.
- Rivers whose parameters are not Muskingum-stable for `dt_routing` are now reported. Stability requires
  `2*k*x <= dt_routing <= 2*k*(1-x)`; outside that window the solution oscillates and clamping the negative
  discharges to zero does not conserve mass. Controlled by the new `unstable_coefficients` config
  (`warn` by default, or `raise` / `ignore`).
- `Configs.deep_validate()` reads every input file that is set and checks its contents. Nothing calls it for
  you: it repeats the read routing is about to do, so it is a method to run once on inputs you have not checked
  before rather than a cost every route pays. There is no `deep_validation` config. It also validates the
  `alpha` and `beta` columns when `coeff` is `dynamic`, and honors `var_river_id`.
- Fixed `dt_total` with `forcing: vlateral`. A value shorter than the input file now routes that portion,
  and one longer than the input file raises instead of reading past the end of the array.
- Unrecognized config keys raise a `ValueError` naming the key and suggesting the closest valid option,
  rather than a dataclass `TypeError`.
- Missing params file columns are reported by name instead of raising `KeyError`.
- `discharge_dir` rejects input files with duplicate basenames, which previously resolved to a single
  output file and silently overwrote each other.
- Repeated `route()` calls on one object restart from `channel_state_init_file` rather than silently
  continuing from the previous run's final state.
- Each `Router` gets its own log handler. Loggers were named from `id(self)`, which CPython reuses after
  garbage collection, so handlers accumulated and log lines were duplicated.
- Removed the invalid `var_vlateral` key from `examples/config.yaml`.
- Added synthetic network tests that run without the downloaded reference data, and CI now runs the test
  suite, ruff, and mypy on every push and pull request.
- Added lower and upper version bounds to all dependencies.

---

### [v3.0.0](https://github.com/rileyhales/river-route/tree/v3.0.0) — 2026-06-18

- Consolidated all routing into a single config-driven `rr.Router`.
- Routing procedure is now selected by config keys rather than class names
    - `coeff` (`static` | `dynamic`)
    - `forcing` (`channel` | `vlateral`)
    - `transform` (`uniform` | `unit_hydrograph`)
    - `network` (`standard`| `expanded`).
- **Removed** the `Muskingum`, `RapidMuskingum`, and `UnitMuskingum` classes
- Added a capability-keyed kernel registry to dispatch the routing kernel from the resolved selectors.
- Simplified the CLI to a single command: `rr route config.yaml`.

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
