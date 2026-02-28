## Watershed Description Files

### Routing Parameters

```yaml
routing_params_file: '/path/to/params.parquet'
```

The routing parameters file is a parquet file. It has 1 row per river in the watershed. The index is ignored.
Rows (rivers) ***must be sorted in topological order*** from upstream to downstream.

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

`UnitMuskingum` does not require additional columns in `routing_params_file`. Catchment-specific parameters
(e.g. `tc`, `area_sqm` for the SCS unit hydrograph) are consumed when building the transformer kernel in
Python and are not read by the router itself.

Generally, use physically derived first guesses then calibrate against observed discharge where available.

## Catchment Runoff Files

You need a time series of per-catchment runoff to be routed. There are 2 ways to provide it:

1. Pre-aggregated catchment files (`catchment_runoff_files`)
2. Gridded runoff depths with a weight table (`runoff_grid_files` + `grid_weights_file`)

!!! warning "Runoff Depths Warning"
There are many projections for grid cells, different names of variables, various file formats, and units
of the runoff depths. You should be certain you can correctly calculate catchment volumes from runoff
depth grids separately before using the calculations performed by `river-route`.

### Pre-aggregated Catchment Files (recommended)

```yaml
catchment_runoff_files: '/path/to/catchment_runoff.nc'
```

Catchment runoff is given in a netCDF file. The file should have 2 dimensions: `time` and `river_id`.
The `river_id` dimension must contain exactly the same IDs *and* be sorted in the same order as the
`river_id` column of the routing parameters file.

The file should have 1 runoff variable (default name `runoff`) which is an array of shape
`(time, river_id)` of dtype float. The variable name can be overridden with `var_catchment_runoff_variable`.

!!! note "Calculating Catchment Volumes"
`river-route` is not a land surface modeling tool. Refer to the example case for guidance on
formatting these files.

### Gridded Runoff Depths

```yaml
runoff_grid_files: [ '/path/to/depths1.nc', '/path/to/depths2.nc', ... ]
grid_weights_file: '/path/to/weight_table.parquet'
```

Runoff depths are given in a netCDF file with 3 dimensions: `time`, `y`, and `x`. The dimension names
can be overridden with `var_t`, `var_y`, and `var_x`. The runoff depth variable name can be overridden
with `var_runoff_depth` (default `'ro'`).

When providing runoff depths, you must also provide a weight table parquet with the following columns:

| Column     | Data Type | Description                                                                   |
|------------|-----------|-------------------------------------------------------------------------------|
| `river_id` | integer   | Unique ID of a river segment                                                  |
| `x_index`  | integer   | The x index of the runoff grid cell that overlaps with the catchment boundary |
| `y_index`  | integer   | The y index of the runoff grid cell that overlaps with the catchment boundary |
| `x`        | float     | The x coordinate of the runoff grid cell                                      |
| `y`        | float     | The y coordinate of the runoff grid cell                                      |
| `area_sqm` | float     | Area of the grid cell–catchment overlap in square meters                      |

!!! note "Ordering Grid Weights"
The order of unique `river_id` values in the weight table should be the same as in the routing
parameters and topologically sorted from upstream to downstream.

Weights need to be recomputed if the grid resolution, grid extent, or catchment boundaries change.

## Output Files

### Routed Discharge

Routed discharge outputs are given in a netCDF file with 2 dimensions: `time` and `river_id`. It will
have 1 variable named `Q` which is an array of shape `(time, river_id)` of dtype float.

You can change the structure of the output file by overriding the default write function.
See the [Advanced Uses](../tutorial/advanced-tutorial.md) page for more information.

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

## UnitMuskingum Transformer State Files (Optional)

```yaml
transformer_kernel_file: '/path/to/kernel.parquet'
transformer_state_init_file: '/path/to/state.parquet'
transformer_state_final_file: '/path/to/final_state.parquet'
```

`UnitMuskingum` reads a pre-computed convolution kernel and can optionally warm-start the transformer
state from a previous run. Both are stored in shape `(n_basins, n_time_steps)`, one row per basin.

- `transformer_kernel_file`: the unit hydrograph kernel parquet. Required for `UnitMuskingum`. Note
  that the kernel depends on `tc`, `area`, **and the routing timestep**.
- `transformer_state_init_file`: warm-start the transformer's rolling state buffer from a prior run.
  Note, the **state depends on the routing timestep**.
- `transformer_state_final_file`: path to write the final transformer state after routing completes,
  for use as `transformer_state_init_file` in a subsequent run.
