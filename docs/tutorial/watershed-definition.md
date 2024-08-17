## Concepts and Vocabulary

`river-route` performs river routing on a vector watershed definition.

- Vector: In the sense of GIS data, vector data are points, lines, and polygons as opposed to gridded/raster data.
- Watershed: a watershed is a contiguous area that all drains to the same outlet. It is a boundary which water does not naturally flow across. No 
  water can flow into a watershed and has exactly 1 outlet.
- Catchment: a catchment is a subunit of a watershed. Water flows into it at the upstream side in exactly 1 location and leaves the catchment in 
  exactly 1 location.
- Subbasin: multiple hydraulically connected catchments forming a complete watershed which is part of a larger watershed.

### Components

A vector watershed definition **_must have_**:

1. Catchment boundary polygons
2. Stream polylines

Some definitions **_may optionally include_**:

1. Watershed boundaries (the dissolved, outermost boundary of all catchments)
2. Nexus points (points marking the location of confluences)
3. Outlet points (points at the downstream most location of a watershed)

### Requirements

These datsets should adhere to the following rules.

- Each catchment should follow 1 stream branch since routing occurs on a stream-to-stream connection, not catchment-to-catchment.
- No divergences or confluences (nexuses) occur **_within_** a catchment. It should occur exactly at the outlet of the catchment.
- Rivers are dendritic (do not diverge).

### Attributes

You need the following minimum attributes for each stream and catchment pair

- ID Number: some unique numeric identifier.
- Downstream ID: the ID of the downstream segment which the given stream flows in to.
- Length: the geodesic length of the river in meters.
- Area: the projected area of the catchment in square meters.
- Muskingum K: the k value to use for this river segment for routing. It can be calibrated later.
- Muskingum X: the x value to use for this river segment for routing. It can be calibrated later.
