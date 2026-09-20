# Configs

Every option for routing and for preparing gridded runoff is held by a frozen `Configs` object. Build one from keyword
arguments or read one from a YAML/JSON file with `Configs.from_file`, then pass it to `Router` or `RunoffGaussianGrid`.
`Router.route` validates it with `validate_routing` and `RunoffGaussianGrid` validates it with `validate_runoff`.

```python
import river_route as rr

configs = rr.Configs.from_file('config.yaml')
rr.Router(configs).route()
```

::: river_route.Configs
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
