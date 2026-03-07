## Watershed Description Files

You can view some example inputs by downloading the sample data used in the GEOGLOWS River Forecast System 
available on AWS S3 at [s3://geoglows-v2/routing-test-data.zip/](https://geoglows-v2.s3.amazonaws.com/routing-test-data.zip).

### Routing Parameters

```yaml
params_file: '/path/to/params.parquet'
```

The routing parameters file is a parquet file. It has 1 row per river in the watershed.
Required for all routers (`Muskingum`, `RapidMuskingum`, `UnitMuskingum`):

| Column                | Data Type | Description                                                |
|-----------------------|-----------|------------------------------------------------------------|
| `river_id`            | integer   | Unique ID of a river segment                               |
| `downstream_river_id` | integer   | ID of downstream river segment, or `-1` for outlet reaches |
| `k`                   | float     | Muskingum `k` parameter (length / velocity)                |
| `x`                   | float     | Muskingum `x` parameter, expected in `[0, 0.5]`            |

These routing parameters typically come from preprocessing and calibration workflows:

1. topology (`river_id`, `downstream_river_id`) from vector network processing
2. channel routing (`k`, `x`) from hydraulic assumptions and/or calibration

!!! warning "Topological Ordering Warning"
    Rows (rivers) ***must be sorted in topological order*** from upstream to downstream.

## Catchment Runoff Files

You need a time series of per-catchment runoff to be routed. There are 2 ways to provide it:

1. Pre-aggregated catchment files (`qlateral_files`)
2. Gridded runoff depths with a weight table (`grid_runoff_files` + `grid_weights_file`)

!!! warning "Runoff Depths Warning"
    There are many projections for grid cells, different names of variables, various file formats, and units of the
    runoff depths. You should be certain you can correctly calculate catchment volumes from runoff depth grids
    separately before using the calculations performed by `river-route`. Do not blindly trust this result!

### Pre-aggregated Catchment Files (recommended)

```yaml
qlateral_files:
  - '/path/to/catchment_runoff.nc'
```

!!! note "Ordering River IDs"
    The `river_id` values **must** be the same values and order as in the routing parameters

Catchment runoff is given as netcdf with 2 dimensions, `time` and `river_id`. The `river_id` dimension **must** contain
exactly the same IDs **and** be sorted in the same order as the `river_id` column of the routing parameters file. It
should have 1 data variable named `qlateral` which is an array of shape `(time, river_id)` of dtype float.
`RapidMuskingum` expects volumes (m³) and `UnitMuskingum` expects depths (m).

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

## UnitMuskingum UH State Files (Optional)

```yaml
transformer_kernel_file: '/path/to/kernel.npz'
uh_state_init_file: '/path/to/state.parquet'
uh_state_final_file: '/path/to/final_state.parquet'
```

`UnitMuskingum` reads a pre-computed convolution kernel and can optionally warm-start the UH
state from a previous run. The kernel is a scipy sparse npz file and the state files are parquet,
both with shape `(n_basins, n_time_steps)`, one row per basin.

- `transformer_kernel_file`: the unit hydrograph kernel (scipy sparse npz). Required for `UnitMuskingum`. Note
  that the kernel depends on `tc`, `area`, **and the routing timestep**.
- `uh_state_init_file`: warm-start the UH rolling state buffer from a prior run.
  Note, the **state depends on the routing timestep**.
- `uh_state_final_file`: path to write the final UH state after routing completes,
  for use as `uh_state_init_file` in a subsequent run.
