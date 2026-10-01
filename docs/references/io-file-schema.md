## Watershed Description Files

### Network File

```json
{
  "network_file": "/path/to/network.parquet"
}
```

The network file is a parquet file. It has 1 row per river in the watershed, and a network stabilized by
`Network.write_stabilized` also has 1 row per sub-reach it added.
Required for all routing:

| Column          | Data Type | Description                                                                     |
|-----------------|-----------|---------------------------------------------------------------------------------|
| `riverId`       | integer   | Unique ID of a river segment                                                    |
| `nextRiverId`   | integer   | ID of downstream river segment, or `-1` for outlet reaches                      |
| `muskingumK`    | float     | Muskingum `k` parameter (length / velocity), in seconds                         |
| `muskingumX`    | float     | Muskingum `x` parameter, expected in `[0, 0.5]`                                 |
| `riverIndex`    | integer   | Position of the river in depth first search order, one more than the row before |
| `upstreamCount` | integer   | Number of rivers upstream of the river, not counting itself                     |

Optional columns:

| Column            | Data Type | Description                                                                          |
|-------------------|-----------|--------------------------------------------------------------------------------------|
| `dynamicAlpha`    | float     | Required for `coefficients: dynamic`: K = dynamicAlpha * Q ^ dynamicBeta             |
| `dynamicBeta`     | float     | Required for `coefficients: dynamic`                                                 |
| `synthetic`       | boolean   | Written by `Network.write_stabilized`: True for an added sub-reach                   |
| `parentRiverId`   | integer   | Written by `Network.write_stabilized`: the river each sub-reach was split from       |

A network file with `synthetic` rows routes with the `channel`, `grid`, and `ecmwf_grib` forcings. The `catchment`
forcing refuses it, since a catchment runoff file has one row per river and none for the added sub-reaches.

These columns typically come from preprocessing and calibration workflows:

1. topology (`riverId`, `nextRiverId`) from vector network processing
2. depth first search order (`riverIndex`, `upstreamCount`) from a depth first search of that topology
3. channel routing (`muskingumK`, `muskingumX`) from hydraulic assumptions and/or calibration

!!! warning "River Ordering Requirement"
    Rows (rivers) ***must be sorted in depth first search (DFS) order***: each river comes after every river upstream
    of it, and the rivers upstream of a river are the rows immediately before it, so every river's whole upstream
    watershed is one contiguous range of rows ending at that river. `riverIndex` numbers the rows in that order and
    `upstreamCount` counts each river's upstream rivers, so a river's watershed is the rows from its `riverIndex`
    minus its `upstreamCount` to its `riverIndex`. That is what divides the rivers into blocks that route
    concurrently, and what `rr subset` cuts a basin with.

    river-route never reorders the network table, so every file with one entry per river (grid weights,
    catchment runoff, channel state) must list the rivers in this same order, and routing refuses runoff that does
    not. `Configs.deep_validate` checks the order and both columns against the topology:

    ```python
    import river_route as rr

    rr.Configs(network_file='/path/to/network.parquet').deep_validate()
    ```

    `examples/migrate_v2_to_v3.py` sorts a v2 params file into DFS order, writes both columns, and reorders the
    files that follow it.

## Catchment Runoff Files

You need a time series of per-catchment runoff to be routed. It is given as `runoff_files`, read in the form `forcing` names:

1. `catchment`: files already aggregated to catchments
2. `grid`: gridded runoff depths with x and y dimensions, aggregated with a weight table (`grid_weights_file`)
3. `ecmwf_grib`: ECMWF GRIB files of runoff depths on a reduced gaussian grid, aggregated with a weight table

!!! warning "Runoff Depths Warning"
    There are many projections for grid cells, different names of variables, various file formats, and units of the
    runoff depths. You should be certain you can correctly calculate catchment volumes from runoff depth grids
    separately before using the calculations performed by `river-route`. Do not blindly trust this result!

### Pre-aggregated Catchment Files (recommended)

```json
{
  "forcing": "catchment",
  "runoff_files": [
    "/path/to/catchment_runoff.nc"
  ]
}
```

!!! note "Ordering River IDs"
    The `riverId` values **must** be the same values and order as the `riverId` column of the network file

Catchment runoff is given as netcdf with 2 dimensions, `riverId` and `time`, in that order, so each river's series is
contiguous and is read straight into the river major arrays the router works in. The `riverId` dimension **must**
contain exactly the same IDs **and** be sorted in the same order as the `riverId` column of the network file. The names and the order are fixed and cannot be configured:

| Variable           | Dimensions           | Description                                                                 |
|--------------------|----------------------|-----------------------------------------------------------------------------|
| `catchment_runoff` | `(riverId, time)`    | Incremental runoff of each catchment per step, as a volume or a depth       |
| `catchment_area`   | `(riverId,)`         | Area of each catchment in m², the factor between depths and volumes         |

The `units` attribute of `catchment_runoff` is required and says which form it takes: `m3` for volumes, or a depth unit
(`m` or `mm`). Depths and volumes are equivalent: routing uses volumes, so depths are converted to meters and
multiplied by `catchment_area` when they are read. `catchment_runoff` names the area variable with the CF attribute
`cell_measures = "area: catchment_area"`. `Runoff.to_netcdf` writes this schema, and the grid runoff classes write it
from their grids with `aggregate_to_file`.

### Gridded Runoff Depths

```json
{
  "forcing": "grid",
  "runoff_files": [
    "/path/to/grid1.nc",
    "/path/to/grid2.nc"
  ],
  "grid_weights_file": "/path/to/weight_table.nc"
}
```

!!! note "Ordering River IDs"
    The `riverId` values **must** be the same values and order as the `riverId` column of the network file

Runoff depths are given in a netCDF file with 3 dimensions: `time`, `y`, and `x`. The dimension names
can be overridden with `var_t`, `var_y`, and `var_x`. The runoff depth variable name can be overridden
with `var_grid_runoff` (default `'ro'`).

Weights need to be recomputed if the grid resolution, grid extent, or catchment boundaries change.
The grid weights netCDF has the following variables, each with the one dimension `index`, a row per weight:

| Column       | Data Type | Description                                                                    |
|--------------|-----------|--------------------------------------------------------------------------------|
| `riverId`    | integer   | Unique ID of a river segment                                                   |
| `x_index`    | integer   | The x index of the runoff grid cell that overlaps with the catchment boundary  |
| `y_index`    | integer   | The y index of the runoff grid cell that overlaps with the catchment boundary  |
| `x`          | float     | The x coordinate of the runoff grid cell                                       |
| `y`          | float     | The y coordinate of the runoff grid cell                                       |
| `area_sqm`   | float     | Area of the grid cell–catchment overlap in square meters                       |
| `proportion` | float     | Fraction of catchment area covered by this grid cell, sums to 1.0 per riverId  |

The names are fixed. A weight table made before v3 names the river id `river_id`, and must be renamed `riverId` before
it is read. The weight table functions of `river_route.runoff.weights` read the catchments' `riverId` column and write
these names.

### Reduced Gaussian Grid Runoff Depths

```json
{
  "forcing": "ecmwf_grib",
  "runoff_files": [
    "/path/to/ro_20260927_00z_cf.grib"
  ],
  "grid_weights_file": "/path/to/gridweights_O1280.nc",
  "grid_accumulation_type": "cumulative"
}
```

Runoff depths are given as ECMWF GRIB files on a global reduced gaussian grid, such as the octahedral O1280 grid
of the IFS. This forcing is specialized to that format; runoff in any other form must be prepared as `catchment` or
`grid` forcing. Every message whose shortName is
`var_grid_runoff` (default `'ro'`) is read with eccodes as one time step, in the order of the files and of their
messages, at its validity date and time. Choosing files whose grid, ensemble member, and steps suit the weight table
and the routing is left to the caller. IFS forecast runoff accumulates from the start of the forecast, so it is
routed with `grid_accumulation_type` cumulative.

`river_route.runoff.ReducedGaussianGrid.from_grib` reads the grid of a file from its metadata: `N` and the number of
cells on each row. Its `cell_polygons` are the area each cell represents, one polygon per cell of the
world in cell order: a box of longitude and latitude halfway to the cells beside it and between latitude edges that
give each row of cells the area of its gaussian quadrature weight, in two parts on either edge of the map for the
cells centered on 180 degrees. `river_route.runoff.reduced_grid_weights` intersects those boxes with the catchments.
Its weight table has the columns of the table above with `cell_index`, the position of the cell in the values of a
GRIB message, in place of `x_index` and `y_index`.

## Channel State Files

A channel state file is a parquet file with two columns, one row per river, listing the rivers of the network file
in the same order. `channel_state_final_file` is written in this format, so a final state starts the next run as its
`channel_state_init_file`.

| Column     | Type    | Description                                            |
|------------|---------|--------------------------------------------------------|
| `riverId`  | integer | ID of the river, in the order of the network file      |
| `Q`        | float   | Discharge of the river in m³/s at the end of the run   |

With `network_type: stabilized` a state file has one row per sub-reach instead, each river's sub-reaches upstream to
downstream with its `riverId` repeated on each. Routing refuses a state file without `riverId` or whose rivers are
not the rivers of the network file in order.

## Output Files

### Routed Discharge

Routed discharge is written by `river_route.router.writers.zarr_writer` unless another writer is set, to a zarr
store with 2 dimensions: `riverId` and `time`. It has 1 variable named `Q` of shape `(riverId, time)` and dtype
float32, chunked so that each chunk holds every time step of a block of rivers.

The values are rounded to `writers.ZARR_KEEPBITS` mantissa bits, a relative error of at most `2^-13`, and each chunk
is compressed with `writers.ZARR_COMPRESSOR`, Blosc lz4 with bitshuffle. On a year of the Amazon that is 2.71x
smaller than the raw array and faster to write than storing it uncompressed, since less of it reaches the disk.

The river dimension comes first in every array format, because that is the layout the kernels write in place: each
river's whole series is contiguous. That is also the layout a writer is handed, as a C-order `(river, time)` array,
so nothing is transposed on the way to the file. Writing `(time, riverId)` instead costs about twice the kernel time
on a large network, since each river's series then has to be transposed out in blocks.

`river_route.router.writers.netcdf_writer` writes the same `(riverId, time)` layout to an uncompressed netCDF file.
