## Migrating from v1 to v2

!!! tip
    For code examples to convert v1 inputs to v2 inputs, check `examples/migrate_v1_to_v2.py`.

---

## Router Class Changes

!!! note
    This is the historical v1 → v2 guide. In **v3** the `Muskingum`, `RapidMuskingum`, and
    `UnitMuskingum` classes below were removed in favor of a single unified `Router` selected by
    config keys (`coeff`/`forcing`/`network`). See the [v2 → v3 migration guide](v2-to-v3.md).

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
| `k`                   | float32 | Muskingum K — wave travel time (seconds)       |
| `x`                   | float32 | Muskingum X — attenuation factor (0 ≤ x ≤ 0.5) |

!!! warning
    Rows **must be topologically sorted** (upstream before downstream).

### Grid Weights

**v1** used a CSV file with columns like `river_id`, `lon_index`, `lat_index`, `lon`, `lat`, `area_sqm`.

**v2** uses a NetCDF file with the following variables on an `index` dimension:

| Variable     | dtype | Description                                                 |
|--------------|-------|-------------------------------------------------------------|
| `river_id`   | int   | Catchment ID                                                |
| `x_index`    | int   | Column index into the runoff grid                           |
| `y_index`    | int   | Row index into the runoff grid                              |
| `x`          | float | X (usually longitude) value of grid cell center             |
| `y`          | float | Y (usually latitude) value of grid cell center              |
| `area_sqm`   | float | Cell/polygon intersection area in square meters             |
| `proportion` | float | **New** Fraction of catchment area (sums to 1 per river_id) |
