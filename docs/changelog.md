## Changelog

---

### [v2.0.0](https://github.com/rileyhales/river-route/tree/v2.0.0) — Unreleased

Major architectural overhaul. See the [Migration Guide](migrating/v1-to-v2.md) for upgrade
instructions.

- Replaced monolithic `Muskingum` class with `Muskingum`, `RapidMuskingum`, and `UnitMuskingum` router hierarchy.
- Introduced `Configs` dataclass replacing the untyped config dictionary.
- Merged `connectivity_file` into `params_file` (single-file network definition).
- Changed grid weights format from CSV to NetCDF with `proportion` column.
- Simplified channel state files from two columns (Q, R) to one column (Q).
- Added topological sort validation on routing parameters.
- Added `types.py` module centralizing type aliases.
- Expanded `runoff.py` with Voronoi-based grid weight computation pipeline.
- Renamed numerous config keys (see migration guide).
- Changed `tools.py` function signatures to accept arrays instead of file paths; removed RAPID conversion utilities.
- Added new documentation tutorials for channel routing, unit hydrograph routing, watershed preparation, and UH kernel generation.
- Increased minimum Python version to 3.12.

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
