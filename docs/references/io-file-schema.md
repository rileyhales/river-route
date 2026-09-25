## Watershed Description Files

You can get example inputs from the GEOGLOWS River Forecast System available on AWS S3 at [s3://geoglows-v2/routing-test-data.zip/](https://geoglows-v2.s3.amazonaws.com/routing-test-data.zip).

### Routing Parameters

```json
{
  "params_file": "/path/to/params.parquet"
}
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

You need a time series of per-catchment runoff to be routed. It is given as `runoff_files`, read as the `runoff_type`:

1. `catchment`: files already aggregated to catchments
2. `gaussian_grid`: gridded runoff depths with x and y dimensions, aggregated with a weight table (`grid_weights_file`)
3. `reduced_gaussian_grid`: gridded runoff depths with one cell dimension (`var_cell`), not implemented yet

!!! warning "Runoff Depths Warning"
    There are many projections for grid cells, different names of variables, various file formats, and units of the
    runoff depths. You should be certain you can correctly calculate catchment volumes from runoff depth grids
    separately before using the calculations performed by `river-route`. Do not blindly trust this result!

### Pre-aggregated Catchment Files (recommended)

```json
{
  "runoff_type": "catchment",
  "runoff_files": [
    "/path/to/catchment_runoff.nc"
  ]
}
```

!!! note "Ordering River IDs"
    The `river_id` values **must** be the same values and order as in the routing parameters

Catchment runoff is given as netcdf with 2 dimensions, `river_id` and `time`, in that order, so each river's series is
contiguous and is read straight into the river major arrays the router works in. The `river_id` dimension **must**
contain exactly the same IDs **and** be sorted in the same order as the `river_id` column of the routing parameters
file. The names and the order are fixed and cannot be configured:

| Variable           | Dimensions           | Description                                                                 |
|--------------------|----------------------|-----------------------------------------------------------------------------|
| `catchment_runoff` | `(river_id, time)`   | Incremental runoff of each catchment per step, as a volume or a depth       |
| `catchment_area`   | `(river_id,)`        | Area of each catchment in m², the factor between depths and volumes         |

The `units` attribute of `catchment_runoff` is required and says which form it takes: `m3` for volumes, or a depth unit
(`m` or `mm`). Depths and volumes are equivalent: routing uses volumes, so depths are converted to meters and
multiplied by `catchment_area` when they are read. `catchment_runoff` names the area variable with the CF attribute
`cell_measures = "area: catchment_area"`. `Runoff.to_netcdf` writes this schema, and the grid runoff classes write it
from their grids with `aggregate_to_file`.

### Gridded Runoff Depths

```json
{
  "runoff_type": "gaussian_grid",
  "runoff_files": [
    "/path/to/grid1.nc",
    "/path/to/grid2.nc"
  ],
  "grid_weights_file": "/path/to/weight_table.nc"
}
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

Routed discharge is written by `river_route.router.writers.zarr_writer` unless another writer is set, to a zarr
store with 2 dimensions: `river_id` and `time`. It has 1 variable named `Q` of shape `(river_id, time)` and dtype
float32, chunked so that each chunk holds every time step of a block of rivers.

The values are rounded to `writers.ZARR_KEEPBITS` mantissa bits, a relative error of at most `2^-13`, and each chunk
is compressed with `writers.ZARR_COMPRESSOR`, Blosc lz4 with bitshuffle. On a year of the Amazon that is 2.71x
smaller than the raw array and faster to write than storing it uncompressed, since less of it reaches the disk.

The river dimension comes first in every array format, because that is the layout the kernels write in place: each
river's whole series is contiguous. That is also the layout a writer is handed, as a C-order `(river, time)` array,
so nothing is transposed on the way to the file. Writing `(time, river_id)` instead costs about twice the kernel time
on a large network, since each river's series then has to be transposed out in blocks.

`river_route.router.writers.netcdf_writer` writes the same `(river_id, time)` layout to an uncompressed netCDF file.
