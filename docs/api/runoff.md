# Runoff

The Runoff classes aggregate runoff to catchments. There is one per `runoff_type`, and
`RUNOFF_CLASS_FOR_RUNOFF_TYPE[configs.runoff_type]` is the class the `runoff_type` of a `Configs` names:

| `runoff_type`           | Class                       | Aggregates                                                                    |
|-------------------------|-----------------------------|-------------------------------------------------------------------------------|
| `catchment`             | `CatchmentRunoff`           | nothing: its files are already aggregated to catchments, as volumes or depths |
| `gaussian_grid`         | `GaussianGridRunoff`        | grids with x and y dimensions, with a weight table                            |
| `reduced_gaussian_grid` | `ReducedGaussianGridRunoff` | grids with one cell dimension. A placeholder, not implemented yet             |

Gridded runoff is aggregated with a grid weight table built with the functions in `river_route.runoff.weights`.
Build a `GaussianGridRunoff` with `GaussianGridRunoff(grid_weights_file, ...)`, or with
`GaussianGridRunoff.from_configs(configs)` to read those same values off a `Configs`, which is how a `Router` builds one.
Only `from_configs` validates, since a directly built `GaussianGridRunoff` has no `Configs` to check.

Every class subclasses the abstract `Runoff`, whose `to_netcdf` writes a catchment runoff array in the format
`CatchmentRunoff` reads. The grid classes precompute that file from their grids with `aggregate_to_file`, so it can be
routed later with `runoff_type` catchment. Routing the grids directly is faster, because the aggregation then happens
inside the routing kernel.

Each class has a `reader` method that yields one `(dates, catchment_runoff, source_file)` tuple per input,
so it can be used on its own or by a `Router`. See [Customizing Runoff Inputs](../tutorial/advanced.md#customizing-runoff-inputs).

```python
import river_route as rr

configs = rr.Configs.from_json('config.json')
runoff = rr.GaussianGridRunoff.from_configs(configs)
router = rr.Router(configs, runoff=runoff)
```

::: river_route.runoff.Runoff
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.GaussianGridRunoff
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.CatchmentRunoff
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.ReducedGaussianGridRunoff
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.runoff.weights
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
