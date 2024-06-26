## Input Files

### Routing Parameters

```yaml
routing_params_file: '/path/to/params.parquet'
```

The routing parameters file is a parquet file. It has 3 columns and 1 row per river in the watershed. The index is
ignored. If you use nonlinear routing, you can provide more columns with k and x values to use at different thresholds
of Q. Nonlinear routing columns should be named k_1, x_1, q_1, followed by k_2, x_2, q_2, for as many thresholds as you
chose. You may provide as many columns as you wish as long as you also provide a thresholds file. The rows (rivers)
***must be sorted in topological order*** from upstream to downstream.

| Column   | Data Type | Description                                                                    |
|----------|-----------|--------------------------------------------------------------------------------|
| river_id | integer   | Unique ID of a river segment                                                   |
| k        | float     | the k parameter of the MuskingumCunge Cunge routing equation length / velocity |
| x        | float     | the x parameter of the MuskingumCunge Cunge routing equation. x : [0, 0.5]     |
| k_1      | float     | Optional, the k parameter of the MuskingumCunge Cunge routing equation at Q_1  |
| x_1      | float     | Optional, the x parameter of the MuskingumCunge Cunge routing equation at Q_1  |
| q_1      | float     | Optional, the minimum value of Q at which to start using use k_1 and x_1       |

### Connectivity File

```yaml
connectivity_file: '/path/to/connectivity.parquet'
```

The connectivity files is a csv with 2 columns and 1 row per river in the watershed. The index is ignored. This file
controls the topology of the rivers in the watershed. Each river segment must have at least 1 downstream segment. If the
river is an outlet then it should have a downstream ID of -1.

To specify the connectivity of braided rivers, a single river ID may have multiple rows with difference IDs given as the
downstream segment. In this case, use the 3rd column to specify the percentage (decimal in the range (0, 1)) of
discharge from the river segment that flows to the downstream segment given on that row. All rivers that are not braided
should have a weight of 1.0. The weights column of rivers that are braided should sum to exactly 1.0 or else water will
be deleted or magically inserted into the rivers.

| Column              | Data Type | Description                                                                                         |
|---------------------|-----------|-----------------------------------------------------------------------------------------------------|
| river_id            | integer   | Unique ID of a river segment                                                                        |
| downstream_river_id | integer   | Unique ID of the downstream river segment                                                           |
| weight              | float     | Optional, the percentage of discharge from this river that should be routed to the downstream river |

### Catchment Volumes

```yaml
runoff_volumes_file: '/path/to/volumes.nc'
```

Runoff volumes are given in a netCDF file.

The file should have 2 dimensions: time and river_id. The times can be given in any unit with a recognizable units
string e.g. 'hours since 2000-01-01 00:00:00'. The river_id dimension should have exactly the same IDs *AND* be sorted
in the same order given in the river_id column of the routing parameters file.

The file should have 1 runoff volumes variables named "m3" which is an array of shape (time, river_id) of dtype float.

## Output Files

### Routed Discharge

Routed discharge outputs are given in a netCDF file.

The file contains 2 dimensions: time and river_id.

It will have 1 variable named "Q" which is an array of shape (time, river_id) of dtype float.

You can change the structure of the output file by overriding the default function to write outputs to disc. See the
[Advanced Skills](../advanced-skills.md)

## Initial and Final State Files

```yaml
initial_state_file: '/path/to/initial.parquet'
final_state_file: '/path/to/final.parquet'
```

State information are stored in parquet files. Muskingum Cunge routing solves for river discharge at time t+1 as a
function of inflow volumes at time t and time t+1, and the discharge at time t. 
