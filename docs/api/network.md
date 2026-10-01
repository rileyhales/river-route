# Network

The river network a simulation routes over: identity, topology, and Muskingum parameters read from a network
file, plus the concurrent routing partition, the stability analysis derived from them, and the stabilized
network that analysis produces. A `Router` builds one
from its `Configs` and reuses it, so a network table is parsed and partitioned once no matter how many
simulations run over it.

::: river_route.network.Network.Network
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
      show_root_full_path: false

## Streams

The functions a `Network` partitions its rivers with, the analysis of how evenly a partition divides a network,
and subsetting a network file, and its grid weights, to one river and every river upstream of it.

::: river_route.network.streams
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
