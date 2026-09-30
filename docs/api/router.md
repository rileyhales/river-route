# Router

All routing is performed by the `Router` class. The routing procedure is described by the `coefficients`, `forcing`,
`transform`, and `network_type` configuration keys, from which the Router chooses the numerical kernel.

A `Router` takes its options from a `Configs` and nothing else, plus optionally the two objects it would
otherwise build for itself: `Router(configs, network=..., runoff=...)`. It owns one simulation: the
coefficients, the time options, the channel state, and the routing loop. The network it routes over belongs to
[`Network`](network.md) and the weight table gridded runoff is aggregated with belongs to
a [`Runoff`](runoff.md) class chosen by `forcing`. The network is built from the same `Configs` when one is not given.
The Runoff is only built when runoff is routed, so channel routing never reads a weight table. A Runoff passed in must be
the class for the `forcing`.
Pass either in to reuse one that already exists across many Routers.

::: river_route.router.Router.Router
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
      show_root_full_path: false
