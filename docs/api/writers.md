# Writers

Premade discharge writers for `Router.set_discharge_writer`. `zarr_writer` is the default.

```python
import river_route as rr

router = rr.Router(rr.Configs.from_json('config.json')).set_discharge_writer(rr.router.writers.netcdf_writer)
```

::: river_route.router.writers
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
