## Changelog

---

### [v3.0.0](https://github.com/rileyhales/river-route/tree/v3.0.0) — 2026-09-29

The [v2 to v3 migration guide](migrating/v2-to-v3.md) gives the replacement for each removed class, config key, and
file format.

**Routing**

- One class, `rr.Router`, routes every combination in place of the `Muskingum`, `RapidMuskingum`, and
  `UnitMuskingum` classes. The routing procedure is chosen by configs: `coefficients` (`static` or `dynamic`),
  `forcing` (`channel`, `catchment`, `grid`, or `ecmwf_grib`), `transform` (`uniform`), and `network_type`
  (`standard` or `stabilized`). Options no routing method routes yet raise `NotImplementedError` before any runoff
  is read.
- Rivers are routed one at a time from upstream to downstream, each river's whole time series before the next river,
  by numba kernels, instead of every river for one time step before the next step by sparse forward substitution.
- `coefficients='dynamic'` routes with nonlinear Muskingum K = `dynamicAlpha` * Q ^ `dynamicBeta`, from those two
  columns of the params file, recomputed at every routing step. Dynamic coefficients route standard networks.
- `network_type='stabilized'` routes every river stably at `dt_routing`: a river too long for it is routed as
  sub-reaches in series, and a river too short for it is sub-cycled. `Network.stabilize(dt)` builds the stabilized
  network in memory and `Network.write_stabilized(dt)` writes it as a parameter table, in which added reaches are
  numbered down from -1,000,000, `synthetic` marks them, and `parent_river_id` maps each reach to the river it was
  split from. Stabilized networks are routed with static coefficients and the `channel`, `grid`, or `ecmwf_grib`
  forcing.
- Rivers whose parameters are not Muskingum stable for `dt_routing`, outside `2*k*x <= dt_routing <= 2*k*(1-x)`,
  are reported. The new `unstable_coefficients` config chooses `warn` (the default), `raise`, or `ignore`.
- `Router.route(thread_pool=pool, threads=n)` routes concurrently on a `ThreadPoolExecutor` the caller passes. The
  region's rivers are divided into blocks, contiguous index ranges of rivers, and the blocks are packed into one job
  per thread. The pool is used as given and never shut down, and the package never creates one: without a pool,
  routing is single threaded.
- The unit hydrograph transform is postponed past v3: `UnitMuskingum`, `river_route.uhkernels`, and the
  `uh_kernel_file`, `uh_state_init_file`, and `uh_state_final_file` configs are removed, and `transform` accepts
  only `uniform`.
- Repeated `route()` calls on one Router restart from `channel_state_init_file` instead of continuing from the
  previous run's final state.
- Each Router has its own log handler. Loggers were named from `id(self)`, which CPython reuses after garbage
  collection, so handlers accumulated and log lines were duplicated.

**Network**

- `rr.Network` owns the river network: the ids, topology, and parameters of the params file, the partition of the
  rivers into blocks, the stability analysis, and stabilization. A Router builds one from its Configs, or is given
  one as `Router(configs, network=network)`, so one parsed and partitioned network backs many simulations. A
  Network is built from a params file or a DataFrame. `Router.river_ids`, `k`, and `x` are read off the network, as
  `router.network.river_ids`.
- The params file column `downstream_river_id` is renamed `next_river_id`. Dynamic coefficients read the columns
  `dynamicAlpha` and `dynamicBeta`.
- The rivers of the params file must be in depth first search (DFS) order, where v2 required only upstream before
  downstream: the rivers upstream of each river are the rows immediately before it, so every river's upstream
  watershed is one contiguous range of rows. Every file with one entry per river lists the rivers in that order.
  `examples/migrate_v2_to_v3.py` converts v2 inputs, sorting the params file into DFS order and the grid weights,
  qlateral files, and channel state files into the same order.
- `river_route.tools` is removed. `subset_configs_to_river`, `connectivity_to_digraph`, and `adjacency_matrix` are in
  `river_route.network.streams` with the functions that analyze and partition networks.

**Runoff**

- `runoff_files` replaces `qlateral_files` and `grid_runoff_files`, and is read in the form `forcing` names.
  `grid_weights_file` is required for the `grid` and `ecmwf_grib` forcings and refused for `catchment`.
- Catchment runoff files replace qlateral files. A file holds a `catchment_runoff` variable with dimensions
  `(river_id, time)`, in that order, and a `catchment_area` variable (m²) with dimension `river_id`, linked by the CF
  attribute `cell_measures = "area: catchment_area"`. The `units` attribute of `catchment_runoff` is required and
  marks it as volumes (`m3`) or depths (`m`, `mm`), which are multiplied by the catchment area when read.
  `Runoff.to_netcdf` writes this schema.
- The runoff classes `CatchmentRunoff`, `GridRunoff`, and `ECMWFGribReducedGrid` read the runoff of each forcing.
  A Router builds the one its `forcing` names, or is given one as `Router(configs, runoff=runoff)` to reuse a read
  weight table. `GridRunoff` reads and indexes its weight table once for every runoff file, and routing reads each
  river's grid cells directly into its forcing as the river is routed. `GridRunoff.to_dataset` and
  `GridRunoff.aggregate_to_file` replace `runoff.runoff_to_qlateral`.
- `ECMWFGribReducedGrid` routes GRIB files on a reduced gaussian grid, such as the octahedral O1280 grid of the ECMWF
  IFS, reading each message with eccodes and keeping only the weight table's cells. At the 18,721 cells of the
  Columbia in a 145 step O1280 forecast that reads in 0.72 s with memory for one message, where xarray and cfgrib
  took 5.8 s and 4 GB. `ReducedGaussianGrid.from_grib` reads a file's grid from its metadata,
  `ReducedGaussianGrid.cell_polygons` gives the area each cell represents, and `reduced_grid_weights` builds a weight
  table that locates each cell by its `cell_index`. The `var_cell` config names the cell dimension.
- `runoff_depth_unit`, `force_positive_runoff`, `force_uniform_timesteps`, and `as_volumes`, arguments of
  `runoff_to_qlateral` in v2, are configs.
- NaN runoff is set to zero before it is routed: a NaN cell contributes nothing while the other cells of its
  catchment still count.
- The weight table functions `grid_weights`, `compute_voronoi_catchment_intersects`,
  `voronoi_diagram_from_regular_xy`, and `cell_xy_from_regular_grid` are in `river_route.runoff.weights`, and still
  importable from `river_route.runoff`.

**Configs**

- `Configs` is frozen and keyword only. Build one from keyword arguments or with `Configs.from_json`, write one with
  `Configs.to_json`, and pass it to `Router(configs)`. Routers no longer take a config file path or options as
  keyword arguments.
- Config files are JSON only. YAML config files are no longer read, and `pyyaml` is no longer a dependency.
  `examples/config.json` lists every option.
- Unrecognized config keys raise a `ValueError` naming them. `Router.route` validates the options and paths with
  `Configs.validate_routing`. `Configs.deep_validate()` reads every input file that is set and checks its contents
  and their consistency with each other; nothing calls it for you.
- Outputs named from `discharge_dir` are `discharge_<input name>.zarr`, taking the extension of the store the
  default writer writes instead of the extension of the runoff file. Input files with duplicate names are rejected
  rather than overwriting each other's output. `discharge_dir=os.devnull` routes and discards the discharge.
- `dt_total` shorter than a runoff file routes that portion, and one longer than the file raises instead of reading
  past the end of the array. Runoff is checked against the params file for its number of rivers before it is
  routed.

**Outputs**

- `river_route.router.writers` holds the discharge writers. `zarr_writer` is the default: a zarr store of `Q` with
  dimensions `(river_id, time)`, chunked so each chunk holds every time step of 500 rivers, rounded to
  `writers.ZARR_KEEPBITS` mantissa bits (a relative error of at most 2^-13) and compressed with Blosc lz4 and
  bitshuffle, 2.71x smaller than the raw array on a year of the Amazon. `netcdf_writer` writes an uncompressed
  netCDF file with dimensions `(river_id, time)`, where v2 wrote `(time, river_id)`. `null_writer` writes nothing.
- `set_write_discharges` is renamed `set_discharge_writer`. A writer is called as
  `writer(router, dates, discharge_array, discharge_file, runoff_file)`: it is handed the Router first, and the
  discharge as a C-order `(river, time)` array.

**CLI**

- `rr route config.json` routes a config file, replacing `rr route --router <class> config` and the per-class
  commands `rr Muskingum`, `rr RapidMuskingum`, and `rr UnitMuskingum`.
- `rr subset <river_id> --params <in> --out-params <out>` cuts a params file, and optionally a grid weight table with
  `--weights` and `--out-weights`, to a river and every river upstream of it.

**Packaging**

- Increased the minimum Python version to 3.14.
- Added lower and upper version bounds to every dependency. Added `zarr`, `cfgrib`, and `eccodes`, and removed
  `pyyaml`.

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
