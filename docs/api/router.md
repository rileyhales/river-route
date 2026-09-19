# Router

All routing is performed by the `Router` class. The routing procedure is described by the `coeff`, `forcing`,
`transform`, and `network` configuration keys and resolved to a numerical kernel internally.

A `Router` takes its options from a `Configs` and nothing else, plus optionally the two objects it would
otherwise build for itself: `Router(configs, network=..., runoff=...)`. It owns one simulation: the
coefficients, the time options, the channel state, and the routing loop. The network it routes over belongs to
[`Network`](network.md) and the weight table gridded runoff is aggregated with belongs to
[`RunoffGaussianGrid`](runoff.md). The network is built from the same `Configs` when one is not given. The `RunoffGaussianGrid` is only built if
routing from `grid_runoff_files` needs one, so a Router routing prepared vlateral never reads a weight table.
Pass either in to reuse one that already exists across many Routers.

::: river_route.Router
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
