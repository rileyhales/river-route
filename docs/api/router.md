# Router

All routing is performed by the `Router` class, Muskingum with static coefficients. What it routes is described by
the `forcing` configuration key, and how each river computes by the compute options `dt_routing`, `use_subreaches`, and
`use_subcycles`.

A `Router` takes a `Configs` and nothing else: `Router(configs)`. It owns one simulation: the coefficients, the time
options, the channel state, and the routing loop. The network it routes over belongs to [`Network`](network.md) and
the weight table gridded runoff is aggregated with belongs to a [`Runoff`](forcings.md) class chosen by `forcing`, and
the Router builds both from its `Configs`. The Runoff is only built when runoff is routed, so channel routing never
reads a weight table.

::: river_route.router.Router.Router
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
      show_root_full_path: false
