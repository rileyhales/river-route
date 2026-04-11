## Changelog

---

### v2.1.0 — 2026-04-11

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v2.1.0) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/2.1/)

- Routing state, parameters (`k`, `x`), channel state, and all numba kernel scratch buffers now use `float32` instead of `float64`. Reduces memory half, speeds up the inner routing loop.
- Fused the sparse matrix-vector multiply and forward-substitution into a single loop for modest speed gain.
- Better use of type aliases to reduce total number of annotations and improve readability.

---

### v2.0.1 — 2026-03-10

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v2.0.1) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/2.0/)

- Adds pyarrow to dependencies which is included by conda installs as a pandas dependency but not in pip.

### v2.0.0 — 2026-03-09

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v2.0.0) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/2.0/)

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

### v1.3.0 — 2025-08-06

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v1.3.0) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/1.3/)

- Added support for routing runoff ensembles.

### v1.2.3 — 2025-05-01

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v1.2.3) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/1.2/)

- Removed volume calculation keyword argument.

### v1.2.2 — 2025-04-23

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v1.2.2) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/1.2/)

- Fixed catchment volume calculation to use the correct area column.

### v1.2.1 — 2025-03-14

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v1.2.1) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/1.2/)

- Fixed bug where the time variable name argument was not applied correctly.

### v1.2.0 — 2025-02-27

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v1.2.0) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/1.2/)

- Improved efficiency of catchment volume computations.

### v1.1.0 — 2025-01-21

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v1.1.0) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/1.1/)

- Switched to a direct solver for the Muskingum linear system.

### v1.0.3 — 2025-01-17

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v1.0.3) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/1.0/)

- Renamed `_MuskingumCunge` to `_Muskingum`.
- Added `metrics.py` module.
- Refactored `runoff.py`.
- Various bug fixes and documentation updates.

### v1.0.2 — 2024-09-24

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v1.0.2) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/1.0/)

- Guaranteed consistent sort order for adjacency matrix construction.

### v1.0.1 — 2024-08-19

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v1.0.1) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/1.0/)

- Fixed `input_type` config value not being set correctly.

### v1.0.0 — 2024-08-17

[:octicons-mark-github-16: View Release](https://github.com/rileyhales/river-route/tree/v1.0.0) · [:octicons-book-16: View Docs](https://river-route.code.hales.app/1.0/)

- Initial stable release.
- Code cleanup and documentation restructuring.
