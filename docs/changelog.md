## Changelog

---

### Unreleased

- Increased the minimum Python version to 3.14.
- Increased the minimum pandas version to 3.0.5. pandas 3 changed `DataFrame.to_numpy()` to return a
  read-only array, which broke the in-place unit conversion when preparing vlateral from runoff with
  irregular timesteps. The 3.0.5 floor also skips 3.0.4, which is yanked from PyPI for segfaults in
  datetime handling that this package hit on every CF encoded time axis it read.
- Increased the minimum numpy version to 2.5.
- Added `river_route.writers` with premade discharge writers for `Router.set_discharge_writer`: `netcdf_writer`, the
  default, and `zarr_writer`, which writes uncompressed `(time, river_id)` discharge in chunks of 500 rivers that
  each span every time step, writing up to the `threads` given to `Router.route` chunks at once (optionally packed
  into shard files with `writers.ZARR_CHUNKS_PER_SHARD`), and `parquet_writer`, which writes one row per
  river and one column per time step, uncompressed and without dictionary encoding or statistics by default
  (`writers.PARQUET_WRITE_OPTIONS`). `null_writer` sends discharge to the operating system's null device
  (`os.devnull`) so a route keeps nothing on disk. zarr is now a dependency.
- Discharge writers now receive the `Router` as their first argument:
  `writer(router, dates, discharge_array, discharge_file, runoff_file)`, so they can read `router.river_ids` and
  `router.cfg`. Custom writers written for the previous 4 argument signature need the new first argument.
  `Router._default_write_discharges` is removed; use `river_route.writers.netcdf_writer`.
- `Router.set_write_discharges` is renamed `Router.set_discharge_writer`.
- `Configs` is its own subpackage, `river_route.configs`, and is the one way options are given. Build a frozen
  `Configs` from keyword arguments or with `Configs.from_file`, `from_json`, or `from_yaml`, then pass it to
  `Router(configs, **overrides)` or `Runoff(configs, **overrides)`; the overrides change a copy. `Router` no longer
  takes a config file path or options as keyword arguments alone, and `Runoff` no longer takes a weight table path
  and options. `Configs.to_json` and `Configs.to_yaml` write the options to a file that `from_file` reads back, and
  `Configs.replace` returns a copy with options changed. `Router.route` calls `Configs.validate_routing` and
  `Runoff` calls `Configs.validate_runoff`, which check the options once and run `Configs.deep_validate` when
  `deep_validation` is True. `Configs.validate` is removed. The runoff options `runoff_depth_unit`,
  `force_positive_runoff`, `force_uniform_timesteps`, and `as_volumes` are now configs.
- Routing from gridded runoff no longer rebuilds, copies, or reallocates lateral inflow for every runoff file.
  The weight table is read and checked against the params file once per `route()`, and each file is aggregated
  in a single pass (area weighting, de-accumulation, clipping, NaN replacement, and the volume product) straight
  into a reused C-order buffer that the kernels read without a copy. The runoff read itself is unchanged. With
  `threads` above 1 the rivers are split into ranges of similar work that the same kernel aggregates
  concurrently on the `thread_pool`; single threaded, one range covers every river.
- Threading happens only on a thread pool the caller passes; the package never creates one. Threads are a runtime
  resource, not a config. `Router.route(thread_pool=pool, threads=n)` uses a `ThreadPoolExecutor` as given and never
  shuts it down, so it can be shared with `Runoff.vlateral(..., thread_pool=pool, threads=n)` and closed by the
  caller's `with` block. `threads` sets how many regions the network is partitioned into when a pool is given;
  without one, routing is single-threaded. `Runoff.aggregate` and `Runoff.to_dataset` take the same arguments.
  `Router.thread_pool()` is removed.
- `river_route.runoff` is now a subpackage following the layout of `river_route.routers`. Runoff preparation is
  the `Runoff` class (also exported as `river_route.Runoff`), which reads the weight table once
  and provides `read_runoff`, `aggregate`, `vlateral` (into a reused buffer), and `to_dataset`. It replaces the
  `runoff_to_vlateral` function: use `Runoff(Configs(grid_weights_file=..., ...)).to_dataset(runoff_files)`. The
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
- `Configs.deep_validate()` is now called from `route()`, controlled by the new `deep_validation` config.
  It also validates the `alpha` and `beta` columns when `coeff` is `dynamic`, and honors `var_river_id`.
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
