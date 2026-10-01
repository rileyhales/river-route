## Migrating from v2 to v3

!!! tip
    See the [changelog](../changelog.md) for every v3 change and the [config file reference](../references/config-files.md) for the current keys.

---

## Major breaking changes

### Watershed Preparation

Watersheds prepared for routing in v3 must be sorted in a depth first search, DFS order, from the headwaters to the outlet, where v2 only
required upstream before downstream (topological order). This facilitates the greatest number of optimizations of the algorithm. This is
readily accomplished with code in many ways but likely not in the GIS software possibly being used to generate the watersheds.

In DFS order the rivers upstream of each river are the rows immediately before it, so every river's whole watershed is one contiguous range of
rows. The network file records that order in two columns: `riverIndex`, each river's position in DFS order, and `upstreamCount`, the number of
rivers upstream of it, so a river's watershed is the rows from its `riverIndex` minus its `upstreamCount` to its `riverIndex`. river-route
never reorders your files, so every file with one entry per river (grid weights, catchment runoff, channel state) must list the rivers in this
same order, and routing refuses runoff that does not.

```python
import river_route as rr

"""check that a network file is in DFS order and that its riverIndex and upstreamCount describe that order"""
rr.Configs(network_file='/path/to/network.parquet').deep_validate()
```

### Runoff and catchments

You no longer need to incorporate code to convert your runoff source to catchment level aggregates and then route those files. It was faster and 
more resource efficient to do this in v2. In v3 routing the grids directly is generally so much more efficient that it is discouraged to
prepare them in advance. The code 
will still accept precalculated catchment level runoff since there are still many reasons you may want your data this way. The difference is that 
this intermediate step should not be thought of as required or best or more efficient.


`examples/migrate_v2_to_v3.py` writes a v2 params file as a v3 network file in DFS order and rewrites the grid weights, qlateral files, and
channel state files in the same order.

## Router Class Changes

The `Muskingum`, `RapidMuskingum`, and `UnitMuskingum` classes are replaced by one class, `Router`. The routing
procedure is chosen by the `forcing` config instead of by the class.

| v2 Class                                  | v3 Replacement                      | Use Case                                           |
|-------------------------------------------|-------------------------------------|----------------------------------------------------|
| `Muskingum`                               | `Router` with `forcing='channel'`   | Channel-only routing (no runoff)                   |
| `RapidMuskingum` with `qlateral_files`    | `Router` with `forcing='catchment'` | Routing runoff already aggregated to catchments    |
| `RapidMuskingum` with `grid_runoff_files` | `Router` with `forcing='grid'`      | Routing gridded runoff depths with a weight table  |
| `UnitMuskingum`                           | _(Proper replacement postponed)_    | Routing with unit hydrograph runoff transformation |

## Configs class

The `rr.Router` class only accepts a `rr.Configs` object rather than all the arguments that were previously passed directly to the class.
The `Configs` object better validates all configs and ensures success better. Configs files are now only in JSON format and YAML files are
no longer supported. They could be easily in any format you prefer but the package reduced to JSON only rather than grow to support the
many format options being entertained.

A `Router` takes a `Configs` object, not a config file path or options as keyword arguments.

```python
import river_route as rr

# v2
"""Router accepted a config file path or kwargs"""
rr.RapidMuskingum('/path/to/config.yaml').route()
rr.RapidMuskingum('/path/to/config.yaml', dt_routing=900).route()

# v3
"""2 options to initialize a config object first then pass to router"""
conf1 = rr.Configs.from_json('/path/to/config.json')
conf2 = rr.Configs(
    network_file='/path/to/network.parquet',
    forcing='catchment',
    runoff_files=['/path/to/catchment_runoff.nc'],
    discharge_dir='/path/to/outputs/',
    dt_routing=900,
)

rr.Router(conf1).route()
```

## Config Key Changes

Config files are JSON only. Copy `examples/config.json`, which lists every key, and move the values of your YAML
file into it. The following config keys have been renamed, removed, or added.

| v2 Key                | v3 Key                                   |
|-----------------------|------------------------------------------|
| `params_file`         | `network_file`                           |
| `qlateral_files`      | `runoff_files` with `forcing: catchment` |
| `grid_runoff_files`   | `runoff_files` with `forcing: grid`      |
| `uh_kernel_file`      | _(removed)_                              |
| `uh_state_init_file`  | _(removed)_                              |
| `uh_state_final_file` | _(removed)_                              |
| `dt_runoff`           | _(removed)_, read from the runoff files  |
| `var_river_id`        | _(removed)_, every file uses `riverId`   |
| `var_discharge`       | _(removed)_, discharge is always `Q`     |

New keys, each with a default that matches v2 except `forcing`, which must be set to route runoff:

| Key                       | Default      | Description                                                                         |
|---------------------------|--------------|-------------------------------------------------------------------------------------|
| `coefficients`            | `'static'`   | `'static'` Muskingum K from `muskingumK`, or `'dynamic'` from `dynamicAlpha`, `dynamicBeta` |
| `forcing`                 | `'channel'`  | `'channel'`, or the form of `runoff_files`: `'catchment'`, `'grid'`, `'ecmwf_grib'` |
| `transform`               | `'uniform'`  | The only option                                                                     |
| `network_type`            | `'standard'` | `'standard'`, or `'stabilized'` to route every resolvable river stably               |
| `unstable_coefficients`   | `'warn'`     | `'warn'`, `'raise'`, or `'ignore'` rivers that are unstable at `dt_routing`         |
| `runoff_depth_unit`       | `None`       | Unit of gridded runoff depths; `None` reads the file attributes                     |
| `force_positive_runoff`   | `false`      | Clip negative runoff depths to zero                                                 |
| `as_volumes`              | `false`      | Prepare catchment runoff as volumes (m³) instead of depths (m)                      |

`forcing` must be set to route runoff: with the default, `channel`, any `runoff_files` are ignored.

---

## API Changes

Most functions kept their names but moved to the module that owns them. The network utilities from `tools` now live with the `Network` class
in `river_route.network.streams`, and the router reads river ids and parameters from its `Network` instead of copying them onto itself.

| v2                                          | v3                                                    |
|---------------------------------------------|-------------------------------------------------------|
| `router.set_write_discharges(func)`         | `router.set_discharge_writer(func)`                   |
| `router.river_ids`, `router.k`, `router.x`  | `router.network.river_ids`, `.k`, `.x`                |
| `river_route.tools.subset_configs_to_river` | `river_route.network.streams.subset_network_to_river` |
| `river_route.tools.connectivity_to_digraph` | _(removed)_                                           |
| `river_route.tools.adjacency_matrix`        | _(removed)_                                           |
| `river_route.runoff.runoff_to_qlateral`     | `river_route.GridRunoff(...).to_dataset`              |
| `river_route.uhkernels`                     | _(removed)_                                           |

---

## File Format Changes

!!! tip
    `examples/migrate_v2_to_v3.py` converts the params file, grid weights, qlateral files, channel state files, and config file from v2 to v3
    in one command. Run it with `--help` to see the arguments. Converting a YAML config file needs `pyyaml`, which v3 no longer depends on, so
    install it for the migration.

### Network File

The v2 params file is the v3 network file, configured with `network_file`. Its columns are camelCase, every outlet's `nextRiverId` must be
`-1`, and it has two new required columns, `riverIndex` and `upstreamCount`, that describe its DFS order; see
[Watershed Preparation](#watershed-preparation).

| v2 Column             | v3 Column                                                       |
|-----------------------|-----------------------------------------------------------------|
| `river_id`            | `riverId`                                                       |
| `downstream_river_id` | `nextRiverId`                                                   |
| `k`                   | `muskingumK`                                                    |
| `x`                   | `muskingumX`                                                    |
|                       | `riverIndex`, each river's position in DFS order                |
|                       | `upstreamCount`, the number of rivers upstream of each river    |

The migration script renames the columns, sorts the rivers into DFS order, and computes both new columns.

### Catchment Runoff (was qlateral)

v3 files are now all river-major C style arrays- meaning array dimensions are the shape `(riverId, time)`. This is a pivot or transpose of 
v2 files which were always of shape `(time, river_id)`. Every v3 file names the river id `riverId`, where v2 files named it `river_id`.

The kernels route one river's whole time series at a time, so keeping each river's series contiguous is what makes reading it fast. The
variable is renamed `catchment_runoff` and needs a `units` attribute, `m3` for volumes or `m` or `mm` for depths. The file also holds a
`catchment_area` variable (m²) which converts depths to volumes. v2 qlateral files were volumes. `Runoff.to_netcdf` writes this format.

```python
import xarray as xr

import river_route as rr

"""catchment areas are only needed for depths but the file stores them either way"""
with xr.open_dataset('/path/to/weights.nc') as weights:
    areas = weights[['river_id', 'area_sqm']].to_dataframe().groupby('river_id', sort=False)['area_sqm'].sum()

"""pivot the v2 qlateral array to (river, time) and write the v3 file, which names its river dimension riverId"""
with xr.open_dataset('/path/to/qlateral.nc') as ds:
    river_ids = ds['river_id'].values
    rr.CatchmentRunoff.to_netcdf(
        '/path/to/catchment_runoff.nc',
        dates=ds['time'].values,
        catchment_runoff=ds['qlateral'].transpose('river_id', 'time').values,
        river_ids=river_ids,
        catchment_area=areas.reindex(river_ids).to_numpy(),
        as_volumes=True,
    )
```

This keeps the rivers in the order of the qlateral file, so it only gives a valid file when that is the order of the network file. Otherwise use
the migration script, which writes the catchment runoff files in the same order as the sorted network file.

You do not need catchment runoff files for gridded runoff anymore. Routing with `forcing: grid` reads each river's grid cells while it routes,
which is faster than writing and then reading an intermediate file, except for cumulative or clipped runoff, which is aggregated on one thread
before it is routed. If you still want the files, `GridRunoff.aggregate_to_file` writes them.

### Grid Weights and Channel State

The grid weights name the river id `riverId` instead of `river_id`, so a v2 weight table must be renamed before it is read. A
channel state file holds a `Q` column and must now also hold each row's `riverId`. Both must list the rivers in the same order as the network
file. If your network file was reordered, these files need to be reordered with it. The migration script does this for you, and renames the
river id of both.

### Routed Discharge

The default writer is now zarr instead of netCDF and outputs in `discharge_dir` are named `discharge_<input name>.zarr`. Discharge is river-major
`(riverId, time)` like every other v3 file. The kernels produce discharge in that layout so writing it needs no transpose, and zarr compresses
it, which makes the files smaller and faster to write. Code that selects by dimension name must use the new name, `ds['Q'].sel(riverId=...)`,
and code that assumes the axis order must transpose. If you prefer netCDF, set `river_route.router.writers.netcdf_writer`, which also writes
`(riverId, time)`, and give its `.nc` output paths with `discharge_files`, since `discharge_dir` names `.zarr` outputs.

```python
import xarray as xr

"""zarr stores open with xarray just like the netCDF files did"""
discharge = xr.open_zarr('/path/to/discharge_runoff_2020.zarr')['Q']
```

---

## Custom Discharge Writers

Writers are now given the `Router` as their first argument so they can read anything they need from it, like the river of each row from
`router.network.original_river_ids` or the options from `router.configs`, and the thread pool and thread count `Router.route` was given as
keywords. The discharge array is river-major `(river, time)` instead of `(time, river)`.

```python
import river_route as rr

# v2
"""writers received the dates, a (time, river) array, and the file paths"""
def write(dates, discharge_array, discharge_file, runoff_file=''):
    ...


rr.RapidMuskingum('/path/to/config.yaml').set_write_discharges(write).route()

# v3
"""writers receive the router first and a (river, time) array"""
def write(router, dates, discharge_array, discharge_file, runoff_file='', *, thread_pool=None, threads=1):
    river_ids = router.network.original_river_ids  # the river of each row of discharge_array
    ...


rr.Router(rr.Configs.from_json('/path/to/config.json')).set_discharge_writer(write).route()
```

---

## CLI

There is only one routing command now since the config file chooses the routing procedure. The new `rr subset` command cuts a network file,
and optionally its grid weights, down to one river and every river upstream of it.

| v2                                             | v3                                                                      |
|------------------------------------------------|-------------------------------------------------------------------------|
| `rr route config.yaml --router RapidMuskingum` | `rr route config.json`                                                  |
| `rr Muskingum config.yaml`                     | `rr route config.json` with `forcing: channel`                          |
| `rr RapidMuskingum config.yaml`                | `rr route config.json` with `forcing` set                               |
| `rr UnitMuskingum config.yaml`                 | _(Proper replacement postponed)_                                        |
|                                                | `rr subset <riverId> --network <in.parquet> --out-network <out.parquet>` |
