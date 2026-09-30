# Runoff

The Runoff classes aggregate runoff to catchments. There is one for each `forcing` that routes runoff, and
`RUNOFF_CLASS_FOR_FORCING[configs.forcing]` is the class the `forcing` of a `Configs` names:

| `forcing`    | Class                  | Aggregates                                                                          |
|--------------|------------------------|-------------------------------------------------------------------------------------|
| `catchment`  | `CatchmentRunoff`      | nothing: its files are already aggregated to catchments, as volumes or depths       |
| `grid`       | `GridRunoff`           | grids with x and y dimensions, with a weight table                                  |
| `ecmwf_grib` | `ECMWFGribReducedGrid` | ECMWF GRIB files on a reduced gaussian grid, read with eccodes, with a weight table |

Gridded runoff is aggregated with a grid weight table built with the functions in `river_route.runoff.weights`.
Build a `GridRunoff` with `GridRunoff(grid_weights_file, ...)`, or with
`GridRunoff.from_configs(configs)` to read those same values off a `Configs`, which is how a `Router` builds one.
Only `from_configs` validates, since a directly built `GridRunoff` has no `Configs` to check.

Every class subclasses the abstract `Runoff`, whose `to_netcdf` writes a catchment runoff array in the format
`CatchmentRunoff` reads. The grid classes precompute that file from their grids with `aggregate_to_file`, so it can be
routed later with `forcing` catchment. Routing the grids directly is faster, because the aggregation then happens
inside the routing kernel.

Each class has a `generator` method that yields one `(dates, runoff, source_file)` tuple per input, with the runoff
in the form routing reads: `CatchmentRunoffVolumes`, or a `GridCellRunoff` whose grid cells
routing reads as it routes each river. It can be used on its own or by a `Router`. See [Customizing Runoff Inputs](../tutorial/advanced.md#customizing-runoff-inputs).

```python
import river_route as rr

configs = rr.Configs.from_json('config.json')
runoff = rr.GridRunoff.from_configs(configs)
router = rr.Router(configs, runoff=runoff)
```

::: river_route.runoff.Runoff
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.runoff.BaseGridRunoff
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.runoff.GridRunoff.GridRunoff
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
      show_root_full_path: false

::: river_route.runoff.CatchmentRunoff.CatchmentRunoff
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
      show_root_full_path: false

::: river_route.runoff.ECMWFGribReducedGrid.ECMWFGribReducedGrid
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
      show_root_full_path: false

::: river_route.runoff.ECMWFGribReducedGrid.ReducedGaussianGrid
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
      show_root_full_path: false

::: river_route.runoff.weights
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
