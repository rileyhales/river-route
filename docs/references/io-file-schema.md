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
Lateral forcing (`forcing: vlateral`) expects runoff volumes (m³). The inflow variable name can be overridden with
`var_vlateral` and the time dimension name with `var_t`. `Runoff.to_netcdf` always writes the defaults.

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

Routed discharge is written by `river_route.router.writers.zarr_writer` unless another writer is set, to a zarr
store with 2 dimensions: `river_id` and `time`. It has 1 variable named `Q` of shape `(river_id, time)` and dtype
float32, chunked so that each chunk holds every time step of a block of rivers.

The values are rounded to `writers.ZARR_KEEPBITS` mantissa bits, a relative error of at most `2^-13`, and each chunk
is compressed with `writers.ZARR_COMPRESSOR`, Blosc lz4 with bitshuffle. On a year of the Amazon that is 2.71x
smaller than the raw array and faster to write than storing it uncompressed, since less of it reaches the disk. A
`float16` run is stored as it is, without rounding, because float16 holds fewer mantissa bits than the rounding
keeps.

The river dimension comes first in every array format, because that is the layout the kernels write in place: each
river's whole series is contiguous. That is also the layout a writer is handed, as a C-order `(river, time)` array,
so nothing is transposed on the way to the file. Writing `(time, river_id)` instead costs about twice the kernel time
on a large network, since each river's series then has to be transposed out in blocks.

`river_route.router.writers.netcdf_writer` writes the same `(river_id, time)` layout to an uncompressed netCDF file.
`river_route.router.writers.parquet_writer` writes a parquet file with a `river_id` column followed by one column per
time step, named `YYYY-MM-DDTHH:MM:SS`. Parquet is columnar, so a river major file would need one column per river,
which is hundreds of thousands of columns on a real network; its rows are rivers instead.

`Configs.discharge_dtype` may be set to `float16` to halve the memory the discharge buffer takes while routing. The
routing math is always float32 and the channel state is never narrowed, so this rounds the saved values only, bounded
by `2^-11` relative to each value.

float16 only covers 6.1e-5 to 65,504. Flows above that overflow to infinity and flows below it lose most of their
precision, so it suits smaller networks rather than the largest basins. On one year of the Amazon (303,097 rivers,
hourly) 0.07% of the routed values, 1.85 million of them, overflow to infinity on the main stem, and `route` logs a
warning whenever the option is used. zarr and parquet store float16 natively, so their files halve as well. netCDF has no half precision type, so
`netcdf_writer` widens to float32 and its file is the same size as a float32 run.
