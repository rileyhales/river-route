# Network

The river network a simulation routes over: identity, topology, and Muskingum parameters read from a network
file, plus the concurrent routing partition, the signs of the routing coefficients derived from them, the sub-reaches
and subcycles that `use_subreaches` and `use_subcycles` route to make every coefficient non-negative, and the stabilized network that writes
those sub-reaches as rows. A `Router` builds one from its `Configs` and reuses it, so a network table is parsed and
partitioned once no matter how many simulations run over it.

::: river_route.network.Network.Network
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
      show_root_full_path: false

::: river_route.network.Network.Layout
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
      show_root_full_path: false

## Splitting

The search for the fewest sub-reaches of unequal travel times, each in subcycles of its own, that
`Network.conditioning` splits a river into where equal sub-reaches and subcycles cannot make it valid alone.

::: river_route.network.splitting
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

## Streams

The functions a `Network` partitions its rivers with, the analysis of how evenly a partition divides a network,
and subsetting a network file, and its grid weights, to one river and every river upstream of it.

::: river_route.network.streams
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
