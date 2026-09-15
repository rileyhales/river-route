# Runoff

Gridded runoff is prepared as lateral inflow (vlateral) for routing by the `Runoff` class, using a grid
weight table built with the functions in `river_route.runoff.weights`.

::: river_route.Runoff
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false

::: river_route.runoff.weights
    handler: python
    options:
      members_order: source
      show_root_heading: true
      show_source: false
