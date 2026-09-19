# Network

The river network a simulation routes over: identity, topology, and Muskingum parameters read from a parameter
table, plus the concurrent routing partition, the stability analysis derived from them, and the stabilized
network that analysis produces. A `Router` builds one
from its `Configs` and reuses it, so a parameter table is parsed and partitioned once no matter how many
simulations run over it.

::: river_route.Network
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.StabilityReport
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.StabilizedNetwork
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

## Streams

The static analysis and graph utilities a `Network` is built on, usable on a parameter table directly:
connectivity validation, stability and compute analysis, partitioning, subdivision, and subsetting.

::: river_route.network.streams
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
