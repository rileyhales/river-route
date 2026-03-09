## Migrating from v1 to v2

!!! tip
    For code examples to convert v1 inputs to v2 inputs, check `examples/migrate_v1_to_v2.py`.

---

## Router Class Changes

The `Muskingum` class has been renamed and is now supplemented with routing options.

| v1 Class    | v2 Replacement   | Use Case                                           |
|-------------|------------------|----------------------------------------------------|
|             | `Muskingum`      | Channel-only routing (no lateral inflow)           |
| `Muskingum` | `RapidMuskingum` | Routing with lateral runoff volumes (RAPID-style)  |
|             | `UnitMuskingum`  | Routing with unit hydrograph runoff transformation |

---

## Config Key Changes

The following config keys have been renamed or removed.

| v1 Key                    | v2 Key                     |
|---------------------------|----------------------------|
| `initial_state_file`      | `channel_state_init_file`  |
| `final_state_file`        | `channel_state_final_file` |
| `input_type`              | `runoff_processing_mode`   |
| `catchment_volumes_files` | `qlateral_files`           |
| `runoff_type`             | `grid_accumulation_type`   |
| `runoff_depths_files`     | `grid_runoff_files`        |
| `weight_table_file`       | `grid_weights_file`        |
| `var_catchment_volume`    | _(removed)_                |
| `connectivity_file`       | _(removed)_                |

---

## File Format Changes

### Routing Parameters

**v1** used two separate parquet files:

- `routing_params_file` — columns: `river_id`, `k`, `x`
- `connectivity_file` — columns: `river_id`, `ds_river_id`

**v2** merges these into a single `params_file` with columns:

| Column                | Type    | Description                                    |
|-----------------------|---------|------------------------------------------------|
| `river_id`            | int64   | Unique segment ID                              |
| `downstream_river_id` | int64   | Downstream segment ID; `-1` at outlets         |
| `k`                   | float64 | Muskingum K — wave travel time (seconds)       |
| `x`                   | float64 | Muskingum X — attenuation factor (0 ≤ x ≤ 0.5) |

!!! warning
    Rows **must be topologically sorted** (upstream before downstream).

### Grid Weights

**v1** used a CSV file with columns like `river_id`, `lon_index`, `lat_index`, `lon`, `lat`, `area_sqm`.

**v2** uses a NetCDF file with the following variables on an `index` dimension:

| Variable     | dtype   | Description                                                 |
|--------------|---------|-------------------------------------------------------------|
| `river_id`   | int64   | Catchment ID                                                |
| `x_index`    | int64   | Column index into the runoff grid                           |
| `y_index`    | int64   | Row index into the runoff grid                              |
| `x`          | float64 | X (usually longitude) value of grid cell center             |
| `y`          | float64 | Y (usually latitude) value of grid cell center              |
| `area_sqm`   | float64 | Cell/polygon intersection area in square meters             |
| `proportion` | float64 | **New** Fraction of catchment area (sums to 1 per river_id) |
