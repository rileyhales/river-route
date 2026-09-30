## Overview

`river-route` routes catchment-scale runoff through a vector river network. All routing runs through the
`Router` class, and the kind of routing it performs is chosen with the `forcing` config selector:

- **`forcing: channel`**: pure channel routing with no runoff entering the rivers. Routes an existing discharge state
  forward in time using only Muskingum channel equations. Requires an explicit initial state.
- **`forcing: catchment`**, **`grid`**, or **`ecmwf_grib`**: routes the runoff volumes or
  depths of `runoff_files`, in that form, directly into river channel inlets at each timestep. This is the most
  common starting point.

This tutorial routes catchment runoff (`forcing: catchment`).

## Vocabulary

- **VPU** (Vector Processing Unit): a named group of catchments and channels forming a complete routing domain.
- **Catchment**: a subunit of a watershed. Water enters at one upstream location and exits at exactly one outlet.
- **Depth first search (DFS) order**: rivers sorted so that each river comes after every river upstream of it and
  the rivers upstream of a river are the rows immediately before it, so every river's upstream watershed is one
  contiguous range of rows ending at that river. The network file must be in DFS order.

## Required Files

Three files are needed for a routing run:

1. **Network file** (`network.parquet`) — river network topology and Muskingum parameters.
2. **Catchment runoff** (`catchment_runoff.nc`) — per-catchment runoff time series.
3. **Routed discharge** (`discharge.zarr`) — output path where results will be written.

See the [File Schemas reference](../references/io-file-schema.md) for field names and formats.

## Network File

The network file must contain at minimum these columns:

| Column          | Description                                                                |
|-----------------|----------------------------------------------------------------------------|
| `riverId`       | Unique integer ID for each river segment                                   |
| `nextRiverId`   | ID of the downstream segment (`-1` at outlets)                             |
| `muskingumK`    | Muskingum K — travel time (seconds); typically channel length / wave speed |
| `muskingumX`    | Muskingum X — attenuation factor (0 ≤ x ≤ 0.5)                             |
| `riverIndex`    | Position of the river in DFS order                                         |
| `upstreamCount` | Number of rivers upstream of the river, not counting itself                |

Rows must be in **DFS order**, which `riverIndex` numbers: a river's upstream watershed is the rows from its
`riverIndex` minus its `upstreamCount` to its `riverIndex`. Every file with one entry per river, such as the catchment
runoff and channel state files, lists the rivers in this same order, and routing refuses one that does not. See the
[File Schemas reference](../references/io-file-schema.md#network-file) for how to check a table.

## Config File

Config values are held by a frozen `Configs` object. Build it from keyword arguments or read it from a JSON
file with `Configs.from_json`, then pass it to `Router`. A `Router` takes its options from a `Configs` and
nowhere else, and a `Configs` is set once when it is built, so an option is changed by building the `Configs` you
want.

```json
{
  "network_file": "/path/to/network.parquet",
  "forcing": "catchment",
  "runoff_files": "/path/to/catchment_runoff.nc",
  "discharge_dir": "/path/to/output/"
}
```

## First Routing Run

```python
import river_route as rr

configs = rr.Configs.from_json('config.json')
rr.Router(configs).route()
```

Or build the configs directly without a config file:

```python
import river_route as rr

configs = rr.Configs(
    network_file='network.parquet',
    runoff_files=['catchment_runoff.nc', ],
    discharge_dir='./output/',
    forcing='catchment',
)
rr.Router(configs).route()
```

A `Configs` is set once, when it is built, and is frozen afterward. There is no method to copy one with an
option changed: build the `Configs` you want. Use `configs.to_json(path)` to write the
options to a file that `Configs.from_json` reads back, e.g. to prepare many jobs for a scheduler.

## Warm-Starting Channel State

By default, the channel starts at zero discharge. Provide a state file to initialize from a previous run:

```json
{
  "network_file": "network.parquet",
  "forcing": "catchment",
  "runoff_files": "catchment_runoff.nc",
  "discharge_dir": "output/",
  "channel_state_init_file": "state.parquet",
  "channel_state_final_file": "new_state.parquet"
}
```

The state file is a parquet with a column `Q` and one row per river segment, in the same order as the routing
params. A final state file also holds each row's `river_id`.

## Reading the Output

The routed discharge output is a zarr store with dimensions `river_id` and `time`, named for its runoff file:

```python
import xarray as xr

river_of_interest = 123456789
ds = xr.open_zarr('output/discharge_catchment_runoff.zarr')
series = ds['Q'].sel(river_id=river_of_interest).to_pandas()

# Save to CSV
series.to_csv('hydrograph.csv')

# Plot
series.plot()
```
