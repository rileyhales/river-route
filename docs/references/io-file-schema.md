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

## Catchment Volumes or Runoff Depths

You need a time series of catchment volumes to be routed. The volumes should be in units of meters cubed and have a
uniform time step for all rivers with no missing values. There are 2 likely ways that you can get this information.

1. Generate it directly using a hydrological model.
2. Generate runoff depth grids and use zonal statistics to calculate catchment scale volumes.

!!! warning "Runoff Depths Warning"
There are many projections for grid cells, different names of variables, various file formats, and units of the runoff depths. You should be
certain you can correctly calculate catchment volumes from runoff depth grids separately before using the calculations performed by `river-route`.
The routing class may not correctly perform the conversions in all cases.

### Catchment Volumes (recommended)

```yaml
lateral_volume_files: '/path/to/volumes.nc'
```

Catchment volumes are given in a netCDF file.

The file should have 2 dimensions: time and river_id. The times can be given in any recognizable unit string. The river_id
dimension should have exactly the same IDs *AND* be sorted in the same order given in the river_id column of the routing parameters file.

The file should have 1 runoff volumes variables named "volume" which is an array of shape (time, river_id) of dtype
float.

!!! note "Calculating Catchment Volumes"
`river-route` is not a land surface modeling tool. It does have an example function illustrating how to perform
the calculations. It will not handle all file formats, land surface and/or hydrology models, etc. Refer to the
example case for guidance on formatting these files.

### Runoff Depths

```yaml
runoff_depth_grids: [ '/path/to/depths1.nc', '/path/to/depths2.nc', ... ]
grid_weights_file: '/path/to/weight_table.nc'
```

Runoff depths are given in a netCDF file.

The file should have 3 dimensions: time, y, and x. The times can be given in any unit with a recognizable units string.
The name of the x and y dimensions can vary and should be specified in the configuration file. The x and y dimensions
should be the same for all files.

The file should have 1 runoff depths variable. The name may vary so you should specify it in the configuration file. It
should contain array of shape (time, y, x) of dtype float.

When providing runoff depths, you must also provide a weight table. The weight table is a netCDF with 1 dimension, index, and 6 variables:
river_id, x_index, y_index, x, y, and area_sqm. There may be multiple rows with the same river_id but which reference difference
runoff grid cells with different x_index and y_index values. The area column should be the area of the grid cell which
overlaps with the catchment boundary and should be in units of meters squared.

!!! note "Ordering Grid Weights"
The order of unique river_id values in the weight table should be the same as in the routing parameters *AND*
should be topologically sorted from upstream to downstream.

| Column   | Data Type | Description                                                                           |
|----------|-----------|---------------------------------------------------------------------------------------|
| river_id | integer   | Unique ID of a river segment                                                          |
| x_index  | integer   | The x index of the runoff grid cell that overlaps with the catchment boundary         |
| y_index  | integer   | The y index of the runoff grid cell that overlaps with the catchment boundary         |
| x        | float     | The x coordinate of the runoff grid cell that overlaps with the catchment boundary    |
| y        | float     | The y coordinate of the runoff grid cell that overlaps with the catchment boundary    |
| area_sqm | float     | The area of the grid cell which overlaps with the catchment boundary in square meters |

## Output Files

### Routed Discharge

Routed discharge outputs are given in a netCDF file.

The file contains 2 dimensions: time and river_id. It will have 1 variable named "Q" which is an array of shape (time, river_id) of dtype float.

You can change the structure of the output file by overriding the default function to write outputs to disc. See the
[Saving Outputs](../tutorial/advanced-tutorial.md) page for more information.

## Initial and Final States

```yaml
channel_state_file: '/path/to/initial.parquet'
final_channel_state_file: '/path/to/final.parquet'
```

State information are stored in parquet files. Muskingum routing solves for river discharge at time t+1 as a
function of inflow volumes at time t and time t+1, and the discharge at time t.

The parquet state file must contain 1 column in river order:

| Column | Description           |
|--------|-----------------------|
| `Q`    | River discharge state |

## UnitMuskingum Transformer State Files (Optional)

```yaml
transformer_kernel_file: '/path/to/kernel.parquet'
transformer_state_file: '/path/to/state.parquet'
```

`UnitMuskingum` can optionally load a pre-computed transformer kernel and/or warm-start the transformer state
from parquet files. Both are stored in **tall format**: shape `(n_basins, n_time_steps)`, one row per basin.
The router transposes these to the internal `(n_time_steps, n_basins)` layout on load.

- `transformer_kernel_file`: skip the kernel computation step and use a previously saved kernel instead. Note, the kernel depends on tc, area, **and the routing timestep**
- `transformer_state_file`: warm-start the transformer's rolling state buffer from a prior run. Note, the **state depends on the routing timestep**!

Use `transformer.save_kernel(path)` to write a kernel produced by `SCSUnitHydrograph` or any other
`AbstractBaseTransformer` subclass.
