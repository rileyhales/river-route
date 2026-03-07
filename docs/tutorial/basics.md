## Overview

`river-route` routes catchment-scale runoff through a vector river network. Three routers are available:

- **`Muskingum`**: pure channel routing with no lateral inflows. Routes an existing discharge state forward in
  time using only Muskingum channel equations. Requires an explicit initial state.
- **`RapidMuskingum`**: routes runoff volumes or depths directly into river channel inlets at each timestep.
  This is the most common starting point.
- **`UnitMuskingum`**: same as `RapidMuskingum` but convolves each timestep of runoff with a unit hydrograph
  kernel before adding it to the channel. See the
  [Channel Routing with Runoff Transformation](unit-hydrograph-routing.md) tutorial.

This tutorial uses `RapidMuskingum`.

## Vocabulary

- **VPU** (Vector Processing Unit): a named group of catchments and channels forming a complete routing domain.
- **Catchment**: a subunit of a watershed. Water enters at one upstream location and exits at exactly one outlet.
- **Topological order**: rivers sorted so that every upstream segment appears before all downstream segments.
  Required by `river-route` — the routing params file must be in topological order.

## Required Files

Three files are needed for a routing run:

1. **Routing parameters** (`params.parquet`) — river network topology and Muskingum coefficients.
2. **Lateral inflow** (`catchment_runoff.nc`) — per-catchment runoff time series.
3. **Routed discharge** (`discharge.nc`) — output path where results will be written.

See the [File Schemas reference](../references/io-file-schema.md) for field names and formats.

## Routing Parameters

The routing parameters parquet must contain at minimum these columns:

| Column                | Description                                                                |
|-----------------------|----------------------------------------------------------------------------|
| `river_id`            | Unique integer ID for each river segment                                   |
| `downstream_river_id` | ID of the downstream segment (`-1` or `<0` at outlets)                     |
| `k`                   | Muskingum K — travel time (seconds); typically channel length / wave speed |
| `x`                   | Muskingum X — attenuation factor (0 ≤ x ≤ 0.5)                             |

Rows must be in **topological order**: all upstream segments before their downstream neighbors.

See [Generating Config Files](prepare-watersheds.md) for guidance on building this file from GIS data.

## Config File

Config values can be passed as a YAML/JSON file, as keyword arguments, or both. Keyword arguments
override values from the config file.

```yaml
params_file: '/path/to/params.parquet'
qlateral_files: '/path/to/catchment_runoff.nc'
discharge_dir: '/path/to/output/'
```

## First Routing Run

```python
import river_route as rr

rr.RapidMuskingum('config.yaml').route()
```

Or pass arguments directly without a config file:

```python
import river_route as rr

(
    rr
    .RapidMuskingum(
        params_file='params.parquet',
        qlateral_files=['qlateral.nc', ],
        discharge_dir='./output/',
    )
    .route()
)
```

## Warm-Starting Channel State

By default, the channel starts at zero discharge. Provide a state file to initialize from a previous run:

```yaml
params_file: 'params.parquet'
qlateral_files: 'catchment_runoff.nc'
discharge_dir: 'output/'
channel_state_init_file: 'state.parquet'         # optional: initial channel state
channel_state_final_file: 'new_state.parquet'    # optional: save final state for next run
```

The state file is a parquet with a single column `Q` and one row per river segment, in the same order
as the routing params.

## Reading the Output

The routed discharge output is a netCDF file with dimensions `time` and `river_id`:

```python
import xarray as xr

river_of_interest = 123456789
ds = xr.open_dataset('discharge.nc')
series = ds['Q'].sel(river_id=river_of_interest).to_pandas()

# Save to CSV
series.to_csv('hydrograph.csv')

# Plot
series.plot()
```
