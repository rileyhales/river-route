## Migrating from v1 to v2

river-route v2 is a major release with breaking changes to the public API, config keys, and input file
formats. This guide covers everything you need to update.

!!! tip
    A migration script is provided at `examples/migrate_v1_to_v2.py` to automate the file format
    conversions described below. See [Using the Migration Script](#using-the-migration-script) for details.

---

## Router Class Changes

The monolithic `Muskingum` class has been replaced by a hierarchy of specialized routers:

| v1 Class    | v2 Replacement   | Use Case                                           |
|-------------|------------------|----------------------------------------------------|
| `Muskingum` | `Muskingum`      | Channel-only routing (no lateral inflow)           |
| `Muskingum` | `RapidMuskingum` | Routing with lateral runoff volumes (RAPID-style)  |
| `Muskingum` | `UnitMuskingum`  | Routing with unit hydrograph runoff transformation |

Update your imports:

```python
# v1
from river_route import Muskingum

m = Muskingum(config_file='config.yaml')

# v2 — choose the router that matches your workflow
from river_route import RapidMuskingum

m = RapidMuskingum(configs='config.yaml')
```

---

## Config Key Renames

The following config keys have been renamed. Update your YAML/JSON config files accordingly:

| v1 Key                    | v2 Key                          |
|---------------------------|---------------------------------|
| `initial_state_file`      | `channel_state_init_file`       |
| `final_state_file`        | `channel_state_final_file`      |
| `input_type`              | `runoff_processing_mode`        |
| `runoff_type`             | `runoff_accumulation_type`      |
| `catchment_volumes_files` | `catchment_runoff_files`        |
| `runoff_depths_files`     | `runoff_grid_files`             |
| `weight_table_file`       | `grid_weights_file`             |
| `var_catchment_volume`    | `var_catchment_runoff_variable` |

The `connectivity_file` key has been removed entirely — see
[Routing Parameters](#routing-parameters) below.

---

## File Format Changes

### Routing Parameters

**v1** used two separate parquet files:

- `routing_params_file` — columns: `river_id`, `k`, `x`
- `connectivity_file` — columns: `river_id`, `ds_river_id`

**v2** merges these into a single `params_file` with columns:

| Column                | Type    | Description                                          |
|-----------------------|---------|------------------------------------------------------|
| `river_id`            | int64   | Unique segment ID                                    |
| `downstream_river_id` | int64   | Downstream segment ID; `-1` (or negative) at outlets |
| `k`                   | float64 | Muskingum K — wave travel time (seconds)             |
| `x`                   | float64 | Muskingum X — attenuation factor (0 ≤ x ≤ 0.5)       |

!!! warning
    Rows **must be topologically sorted** (upstream before downstream). v2 validates this and raises
    a `ValueError` if the ordering is wrong. See
    [Preparing Watersheds](../tutorial/prepare-watersheds.md) for how to sort your data.

### Grid Weights

**v1** used a CSV file with columns like `river_id`, `lon_index`, `lat_index`, `lon`, `lat`, `area_sqm`.

**v2** uses a NetCDF file with the following variables on an `index` dimension:

| Variable     | dtype   | Description                                         |
|--------------|---------|-----------------------------------------------------|
| `river_id`   | int64   | Catchment ID                                        |
| `x_index`    | int64   | Column index into the runoff grid                   |
| `y_index`    | int64   | Row index into the runoff grid                      |
| `x`          | float64 | Longitude of grid cell center                       |
| `y`          | float64 | Latitude of grid cell center                        |
| `area_sqm`   | float64 | Intersection area in square meters                  |
| `proportion` | float64 | Fraction of catchment area (sums to 1 per river_id) |

The `proportion` column is new in v2 and is required.

### State Files

**v1** state files had two columns: `Q` and `R`.

**v2** state files have a single column: `Q`. The `R` state variable has been removed.

---

## API Changes in `tools.py`

### Removed Functions

- `routing_files_from_RAPID()` — no longer needed; RAPID CSV conversion is not supported.
- `inits_from_RAPID()` — removed along with the RAPID import workflow.

### Changed Signatures

`get_adjacency_matrix` has been renamed and now takes arrays instead of file paths:

```python
# v1
from river_route import tools
A = tools.get_adjacency_matrix('routing_params.parquet', 'connectivity.parquet')

# v2
import pandas as pd
from river_route import tools

df = pd.read_parquet('routing_params.parquet')
A = tools.adjacency_matrix(df['river_id'].values, df['downstream_river_id'].values)
```

`connectivity_to_digraph` also takes arrays instead of a file path:

```python
# v1
G = rr.tools.connectivity_to_digraph('connectivity.parquet')

# v2
G = rr.tools.connectivity_to_digraph(df['river_id'].values, df['downstream_river_id'].values)
```

---

## Using the Migration Script

The migration script at `examples/migrate_v1_to_v2.py` automates the three file format conversions.

### Merge Routing Parameters

Combine your v1 routing params and connectivity files into a single parquet:

```bash
python examples/migrate_v1_to_v2.py routing-params \
    routing_params.parquet connectivity.parquet \
    -o params_v2.parquet
```

This handles legacy column renames (`ds_river_id` → `downstream_river_id`) automatically.

!!! note
    The merged file must still be **topologically sorted**. If your v1 files were not sorted, you
    will need to sort them before using them with v2. See
    [Preparing Watersheds](../tutorial/prepare-watersheds.md) for a topological sort example.

### Convert Grid Weights

Convert a v1 CSV grid weights file to a v2 NetCDF file:

```bash
python examples/migrate_v1_to_v2.py grid-weights \
    weight_table.csv \
    -o grid_weights.nc
```

If the CSV has an `area_sqm` column but no `proportion` column, the script computes proportions
automatically.

### Convert State Files

Strip the `R` column from a v1 state file:

```bash
python examples/migrate_v1_to_v2.py state-file \
    initial_state.parquet \
    -o channel_state.parquet
```

---

## Manual Steps

The migration script handles file format conversions, but you will also need to:

1. **Rename config keys** in your YAML/JSON files (see [Config Key Renames](#config-key-renames)).
2. **Update your Python code** to import and use the new router classes (see
   [Router Class Changes](#router-class-changes)).
3. **Update any calls to `tools.py` functions** that changed signatures (see
   [API Changes](#api-changes-in-toolspy)).
4. **Topologically sort your routing parameters** if they are not already sorted.
