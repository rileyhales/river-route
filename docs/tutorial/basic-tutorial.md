!!! warning "Non-Comprehensive Tutorial"
    This is not a comprehensive hydrological modeling course. It will not teach watershed delineation,
    parameter calibration, or GIS skills. An approachable starting point is a GIS tutorial on watershed
    delineation and stream network extraction.

## Overview

`river-route` routes catchment-scale runoff through a vector river network using the Muskingum method.
Three routers are available:

- **`Router`**: pure channel routing with no lateral inflows. Routes an existing discharge state forward
  in time. Requires a non-zero initial state, `dt_routing`, and `dt_total`.
- **`RapidMuskingum`**: routes per-catchment runoff volumes or depths into channel inlets each timestep.
  The most common starting point.
- **`UnitMuskingum`**: same as `RapidMuskingum` but convolves each timestep of runoff with a precomputed
  unit hydrograph kernel before passing it to the channel equations.

## Quickstart

Three files are required for a `RapidMuskingum` run:

1. **Routing parameters** — parquet file describing the river network topology and Muskingum coefficients
2. **Catchment runoff** — netCDF file with per-catchment runoff volumes or depths over time
3. **Routed discharge** — output path where results will be written

```yaml
# config.yaml
routing_params_file: '/path/to/params.parquet'
catchment_runoff_files: '/path/to/catchment_runoff.nc'
discharge_files: '/path/to/discharge.nc'
dt_routing: 3600
```

```python
import river_route as rr

rr.RapidMuskingum('config.yaml').route()
```

Or pass arguments directly without a config file:

```python
import river_route as rr

rr.RapidMuskingum(
    routing_params_file='params.parquet',
    catchment_runoff_files='catchment_runoff.nc',
    discharge_files='discharge.nc',
    dt_routing=3600,
).route()
```

## Reading the Output

The routed discharge is written as a netCDF file with dimensions `time` and `river_id`:

```python
import xarray as xr

ds = xr.open_dataset('discharge.nc')
series = ds['Q'].sel(river_id=123456789).to_pandas()

series.to_csv('hydrograph.csv')
series.plot()
```

## Next Steps

- [Preparing Watersheds](prepare-watersheds.md) — build routing parameter files and config files from GIS data
- [Channel Routing](channel-routing.md) — sequential runs, warm-starting channel state, time step options
- [Channel Routing with Runoff Transformation](unit-hydrograph-routing.md) — unit hydrograph convolution with `UnitMuskingum`
- [Routing Runoff Ensembles](routing-ensembles.md) — parallel and sequential ensemble runs
- [Advanced Uses](advanced-tutorial.md) — customizing outputs, generating configs at runtime