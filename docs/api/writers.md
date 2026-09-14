# Writers

Premade discharge writers for `Router.set_write_discharges`. `netcdf_writer` is the default.

```python
import river_route as rr

router = rr.Router('config.yaml').set_write_discharges(rr.writers.zarr_writer)
```

::: river_route.writers
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
