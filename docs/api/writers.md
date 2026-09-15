# Writers

Premade discharge writers for `Router.set_discharge_writer`. `netcdf_writer` is the default.

```python
import river_route as rr

router = rr.Router(rr.Configs.from_file('config.yaml')).set_discharge_writer(rr.writers.zarr_writer)
```

::: river_route.writers
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
