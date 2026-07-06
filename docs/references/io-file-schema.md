## Watershed Description Files

You can get example inputs from the GEOGLOWS River Forecast System available on AWS S3 at [s3://geoglows-v2/routing-test-data.zip/](https://geoglows-v2.s3.amazonaws.com/routing-test-data.zip).

### Routing Parameters

```yaml
params_file: '/path/to/params.parquet'
```

The routing parameters file is a parquet file. It has 1 row per river in the watershed.
Required for all routing:

| Column          | Data Type | Description                                                |
|-----------------|-----------|------------------------------------------------------------|
| `river_id`      | integer   | Unique ID of a river segment                               |
| `next_river_id` | integer   | ID of downstream river segment, or `-1` for outlet reaches |
| `k`             | float     | Muskingum `k` parameter (length / velocity)                |
| `x`             | float     | Muskingum `x` parameter, expected in `[0, 0.5]`            |

These routing parameters typically come from preprocessing and calibration workflows:

1. topology (`river_id`, `next_river_id`) from vector network processing
2. channel routing (`k`, `x`) from hydraulic assumptions and/or calibration

!!! warning "Topological Ordering Warning"
    Rows (rivers) ***must be sorted in topological order*** from upstream to downstream.

## Catchment Runoff Files

You need a time series of per-catchment runoff to be routed. There are 2 ways to provide it:

1. Pre-aggregated catchment files (`vlateral_files`)
2. Gridded runoff depths with a weight table (`grid_runoff_files` + `grid_weights_file`)

!!! warning "Runoff Depths Warning"
    There are many projections for grid cells, different names of variables, various file formats, and units of the
    runoff depths. You should be certain you can correctly calculate catchment volumes from runoff depth grids
    separately before using the calculations performed by `river-route`. Do not blindly trust this result!

### Pre-aggregated Catchment Files (recommended)

```yaml
vlateral_files:
  - '/path/to/catchment_runoff.nc'
```

!!! note "Ordering River IDs"
    The `river_id` values **must** be the same values and order as in the routing parameters

Catchment runoff is given as netcdf with 2 dimensions, `time` and `river_id`. The `river_id` dimension **must** contain
exactly the same IDs **and** be sorted in the same order as the `river_id` column of the routing parameters file. It
should have 1 data variable named `vlateral` which is an array of shape `(time, river_id)` of dtype float.
Lateral forcing (`forcing: vlateral`) expects runoff volumes (m³).

### Gridded Runoff Depths

```yaml
grid_runoff_files:
  - '/path/to/grid1.nc'
  - '/path/to/grid2.nc'
grid_weights_file: '/path/to/weight_table.nc'
```

!!! note "Ordering River IDs"
    The `river_id` values **must** be the same values and order as in the routing parameters

Runoff depths are given in a netCDF file with 3 dimensions: `time`, `y`, and `x`. The dimension names
can be overridden with `var_t`, `var_y`, and `var_x`. The runoff depth variable name can be overridden
with `var_grid_runoff` (default `'ro'`).

Weights need to be recomputed if the grid resolution, grid extent, or catchment boundaries change.
The grid weights netCDF has the following variables:

| Column       | Data Type | Description                                                                    |
|--------------|-----------|--------------------------------------------------------------------------------|
| `river_id`   | integer   | Unique ID of a river segment                                                   |
| `x_index`    | integer   | The x index of the runoff grid cell that overlaps with the catchment boundary  |
| `y_index`    | integer   | The y index of the runoff grid cell that overlaps with the catchment boundary  |
| `x`          | float     | The x coordinate of the runoff grid cell                                       |
| `y`          | float     | The y coordinate of the runoff grid cell                                       |
| `area_sqm`   | float     | Area of the grid cell–catchment overlap in square meters                       |
| `proportion` | float     | Fraction of catchment area covered by this grid cell, sums to 1.0 per river_id |

## Output Files

### Routed Discharge

Routed discharge outputs are given in a netCDF file with 2 dimensions: `time` and `river_id`. It will
have 1 variable named `Q` which is an array of shape `(time, river_id)` of dtype float.

You can change the structure of the output file by overriding the default write function.
See the [Advanced Uses](../tutorial/advanced.md) page for more information.

## Initial and Final States

```yaml
channel_state_init_file: '/path/to/initial.parquet'
channel_state_final_file: '/path/to/final.parquet'
```

State information is stored in parquet files. Muskingum routing solves for river discharge at time `t+1`
as a function of inflow at time `t` and `t+1`, and discharge at time `t`.

The parquet state file must contain 1 column in river order:

| Column | Description           |
|--------|-----------------------|
| `Q`    | River discharge state |
