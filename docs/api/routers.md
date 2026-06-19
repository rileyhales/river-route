# Router

All routing is performed by a single `Router` class. The routing procedure is described by the `coeff`,
`forcing`, and `network` configuration keys and resolved to a numerical kernel internally. See the
[v2 → v3 migration guide](../migrating/v2-to-v3.md) for the mapping from the former router classes.

::: river_route.Router
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.Configs
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
