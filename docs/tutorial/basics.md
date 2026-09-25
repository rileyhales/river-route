## Overview

`river-route` routes catchment-scale runoff through a vector river network. All routing runs through the
`Router` class, and the kind of routing it performs is chosen with the `forcing` config selector:

- **`forcing: channel`**: pure channel routing with no lateral inflows. Routes an existing discharge state
  forward in time using only Muskingum channel equations. Requires an explicit initial state.
- **`forcing: runoff`**: routes runoff volumes or depths directly into river channel inlets at each
  timestep. This is the most common starting point.

This tutorial uses lateral-runoff routing (`forcing: runoff`).

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

| Column          | Description                                                                |
|-----------------|----------------------------------------------------------------------------|
| `river_id`      | Unique integer ID for each river segment                                   |
| `next_river_id` | ID of the downstream segment (`-1` or `<0` at outlets)                     |
| `k`             | Muskingum K — travel time (seconds); typically channel length / wave speed |
| `x`             | Muskingum X — attenuation factor (0 ≤ x ≤ 0.5)                             |

Rows must be in **topological order**: all upstream segments before their downstream neighbors.

## Config File

Config values are held by a frozen `Configs` object. Build it from keyword arguments or read it from a JSON
file with `Configs.from_json`, then pass it to `Router`. A `Router` takes its options from a `Configs` and
nowhere else, and a `Configs` is set once when it is built, so an option is changed by building the `Configs` you
want.

```json
{
  "params_file": "/path/to/params.parquet",
  "forcing": "runoff",
  "runoff_type": "catchment",
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
    params_file='params.parquet',
    runoff_files=['catchment_runoff.nc', ],
    discharge_dir='./output/',
    forcing='runoff',
    runoff_type='catchment',
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
  "params_file": "params.parquet",
  "forcing": "runoff",
  "runoff_type": "catchment",
  "runoff_files": "catchment_runoff.nc",
  "discharge_dir": "output/",
  "channel_state_init_file": "state.parquet",
  "channel_state_final_file": "new_state.parquet"
}
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
