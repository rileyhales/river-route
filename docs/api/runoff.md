# Runoff

Gridded runoff is prepared as lateral inflow (vlateral) for routing by the `RunoffGaussianGrid` class, using a grid
weight table built with the functions in `river_route.runoff.weights`.

Build one with `RunoffGaussianGrid(grid_weights_file, ...)`, or with `RunoffGaussianGrid.from_configs(configs)` to read those same values
off a `Configs`, which is how a `Router` builds one when routing from `grid_runoff_files`. Only `from_configs`
validates, since a directly built `RunoffGaussianGrid` has no `Configs` to check.

Both classes subclass the abstract `Runoff`, whose `to_netcdf` writes a vlateral array in the format `RunoffVlateral`
reads, so aggregated grids can be saved once and routed later from `vlateral_files`.

Lateral inflow already prepared as volumes is read by the `RunoffVlateral` class, which a `Router` uses when routing
from `vlateral_files`.

Each class has a `reader` method that yields one `(dates, vlateral, source_file)` tuple per input,
so it can be used on its own or by a `Router`. See [Customizing Runoff Inputs](../tutorial/advanced.md#customizing-runoff-inputs).

```python
import river_route as rr

configs = rr.Configs.from_file('config.yaml')
runoff = rr.RunoffGaussianGrid.from_configs(configs)
router = rr.Router(configs, runoff=runoff)
```

::: river_route.runoff.Runoff
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.RunoffGaussianGrid
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.RunoffVlateral
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
