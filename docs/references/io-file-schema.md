## Watershed Description Files

### Routing Parameters

```yaml
routing_params_file: '/path/to/params.parquet'
```

The routing parameters file is a parquet file. It has 3 columns and 1 row per river in the watershed. The index is
ignored. The rows (rivers) ***must be sorted in topological order*** from upstream to downstream.

| Column   | Data Type | Description                                                         |
|----------|-----------|---------------------------------------------------------------------|
| river_id | integer   | Unique ID of a river segment                                        |
| k        | float     | the k parameter of the Muskingum routing equation length / velocity |
| x        | float     | the x parameter of the Muskingum routing equation. x : [0, 0.5]     |

### Connectivity File

```yaml
connectivity_file: '/path/to/connectivity.parquet'
```

The connectivity file is a csv with 2 columns and 1 row per river in the watershed. The index is ignored. This file
controls the topology of the rivers in the watershed. Each river segment must have at least 1 downstream segment. If the
river is an outlet then the downstream ID should be -1.

To specify the connectivity of braided rivers, a single river ID may have multiple rows with different IDs given as the
downstream segment. In this case, use the 3rd column to specify the percentage (decimal in the range (0, 1)) of
discharge from the river segment that flows to the downstream segment given on that row. All rivers that are not braided
should have a weight of 1.0. The weights column of rivers that are braided should sum to exactly 1.0 or else water will
be deleted or magically inserted into the rivers.

| Column              | Data Type | Description                                                                                         |
|---------------------|-----------|-----------------------------------------------------------------------------------------------------|
| river_id            | integer   | Unique ID of a river segment                                                                        |
| downstream_river_id | integer   | Unique ID of the downstream river segment                                                           |
| weight              | float     | Optional, the percentage of discharge from this river that should be routed to the downstream river |

## Catchment Volumes or Runoff Depths

You need a time series of catchment volumes to be routed. The volumes should be in units of meters cubed and have a
uniform time step for all rivers with no missing values. There are 2 likely ways that you can get this information.

1. Generate it directly using a hydrological model.
2. Generate runoff depth grids and use zonal statistics to calculate catchment scale volumes.

!!! warning "Runoff Depths Warning"
There is a wide variety of projections for the grid cells, different names of variables, various file formats, and
units of the runoff depths. For these reasons, you should be certain you can correctly calculate catchment volumes
from runoff depth grids separately before using the calculations performed by `river-route`. The routing class may
not correctly perform the conversions in all cases.

### Catchment Volumes (recommended)

```yaml
catchment_volumes_file: '/path/to/volumes.nc'
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
runoff_depths_files: [ '/path/to/depths1.nc', '/path/to/depths2.nc', ... ]
weight_table_file: '/path/to/weight_table.nc'
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
    The order of unique river_id values in the table made from the weight table should be the same as the order of the river_id column in the routing 
    parameters parquet *AND* should be topologically sorted from upstream to downstream.

| Column     | Data Type | Description                                                                           |
|------------|-----------|---------------------------------------------------------------------------------------|
| river_id   | integer   | Unique ID of a river segment                                                          |
| x_index    | integer   | The x index of the runoff grid cell that overlaps with the catchment boundary         |
| y_index    | integer   | The y index of the runoff grid cell that overlaps with the catchment boundary         |
| x          | float     | The x coordinate of the runoff grid cell that overlaps with the catchment boundary    |
| y          | float     | The y coordinate of the runoff grid cell that overlaps with the catchment boundary    |
| area_sqm   | float     | The area of the grid cell which overlaps with the catchment boundary in square meters |


## Output Files

### Routed Discharge

Routed discharge outputs are given in a netCDF file.

The file contains 2 dimensions: time and river_id.

It will have 1 variable named "Q" which is an array of shape (time, river_id) of dtype float.

You can change the structure of the output file by overriding the default function to write outputs to disc. See the
[Saving Outputs](../tutorial/advanced-tutorial.md) page for more information.

## Initial and Final State Files

```yaml
initial_state_file: '/path/to/initial.parquet'
final_state_file: '/path/to/final.parquet'
```

State information are stored in parquet files. Muskingum routing solves for river discharge at time t+1 as a
function of inflow volumes at time t and time t+1, and the discharge at time t. 
