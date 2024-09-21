!!! warning "Non-Comprehensive Tutorial"
    This is not a comprehensive hydrological modeling course which should teach you the theory of hydrological channel routing, calibration, and
    validation. It will not teach the prerequisite informatics and GIS skills to create a watershed representation. Hydraulics and hydrology software
    can be used to obtain the information needed for routing, or you can use a GIS software such as QGIS or ArcGIS. An approachable place to start
    learning these skills is a GIS tutorial demonstrating watershed delineation and stream network extraction as well as assigning and calculating
    attributes of the stream features.

## Overview

`river-route` is a Python package that routes catchment scale runoff volume time series on a vector watershed definition. It takes 3 input datasets
and writes 1 output dataset. Together, that makes a total of 4 files that need to be specified. The first 2 describe the river channel properties and
topology, the third is the input water being routed, and the fourth (the output) is the discharge time series calculated by the routing process. This
tutorial explain the process of preparing these files.

1. [Routing Parameters](../references/io-file-schema.md#routing-parameters) - parquet file
2. [Connectivity File](../references/io-file-schema.md#connectivity-file) - parquet file
3. [Catchment Volumes](../references/io-file-schema.md#catchment-volumes-or-runoff-depths) - netCDF file
4. [Routed Discharge](../references/io-file-schema.md#routed-discharge) - netCDF file

## Vocabulary

The following vocabulary is used in this tutorial.

- VPU: Vector Processing Unit. This refers to the collection of sub-basins or complete watersheds that you group together for computations. The term
  VPU is not required by river-route or a technical definition in hydrology or computational science. You can think of this as an ID or name of the
  watershed you are simulating.
- Vector: In the sense of GIS data, vector data are points, lines, and polygons as opposed to gridded/raster data.
- Watershed: a watershed is a contiguous area that all drains to the same outlet. It is a boundary which water does not naturally flow across. No
  water can flow into a watershed and has exactly 1 outlet.
- Catchment: a catchment is a subunit of a watershed. Water flows into it at the upstream side in exactly 1 location and leaves the catchment in
  exactly 1 location.
- Subbasin: multiple hydraulically connected catchments forming a complete watershed which is part of a larger watershed.

## Directory Structure

***You are not required to use a specific file structure or naming conventions***. If you do not already have an existing structure to match, you
could consider using the following pattern:

```
project_root_directory
├── configs
│   ├── <VPU Name 1>
│   │   ├── params.parquet
│   │   ├── connectivity.parquet
│   │   └── weight_table.csv
│   └── <VPU Name 2>
│       ├── params.parquet
│       ├── connectivity.parquet
│       └── weight_table.csv
├── volumes
│   ├── <VPU Name 1>
│   │   ├── volume_<vpu-name>_<first-time-step-or-ensemble-member>.nc
│   │   └── volume_<vpu-name>_<next-time-step-or-ensemble-member>.nc
│   └── <VPU Name 2>
│       ├── volume_<vpu-name>_<first-time-step-or-ensemble-member>.nc
│       └── volume_<vpu-name>_<next-time-step-or-ensemble-member>.nc
└── discharge
    ├── <VPU Name 1>
    │   ├── discharge_<vpu-name>_<first-time-step-or-ensemble-member>.nc
    │   └── discharge_<vpu-name>_<next-time-step-or-ensemble-member>.nc
    └── <VPU Name 2>
        ├── discharge_<vpu-name>_<first-time-step-or-ensemble-member>.nc
        └── discharge_<vpu-name>_<next-time-step-or-ensemble-member>.nc
```

## Preparing Vector Watersheds

A watershed definition is a series of GIS datasets that satisfy these criteria.

### Components

A vector watershed definition **_must have_**:

1. Catchment boundary polygons
2. Stream polylines

Some watershed definitions **_may optionally include_**:

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

## Runoff Depths vs Catchment Volumes

Routing requires a catchment level runoff volume time series. Runoff data is most commonly generated on regular grids from a land surface process
model. Runoff grid files may have a single file containing multiple time steps or multiple files which each contain a single time step. You might also
have multiple realizations of those data for each member in an ensemble simulation.

### Weight Table

You can use GIS methods to aggregate the distributed runoff depths to catchment level volumes. The steps are (1) intersect the grid cell boundaries
with the catchment boundaries to create polygons of each grid cell and catchment combination, then (2) multiply the runoff depth value of each polygon
times the area of that polygon (i.e. convert depth to volume). Repeat this for each time step. Repeat for each ensemble member, if applicable.

Intersections of grid cells and catchment boundaries are computationally expensive and do not change. One method to avoid repeating the expensive
intersection step (and therefore increase speed) is to cache the intersection results in a "weight table" that describes which grid cells contribute
to which catchments. The cached results are unique to the grid cell geometry (and catchment boundaries). You need to recreate the cached intersections
if you use runoff grids with different resolutions or extents or if you change your watershed definition (catchment boundaries).

## Routing

### Preparing a Config File

Your configuration file must contain at least the 4 essential file paths. There are many other options you can specify
that may be convenient for your use case or required for your specific datasets and file structures.

You may save the output file to any file path you select. Enter the correct file path to the other 3 existing files and
enter them into a YAML file that looks like this example.

```yaml
routing_params_file: '/path/to/params.parquet'
connectivity_file: '/path/to/connectivity.parquet'
catchment_volumes_file: '/path/to/volumes.nc'
outflow_file: '/path/to/outflows.nc'
```

### First Routing Run

You should now have the 4 essential files prepared: 2 watershed descriptors, 1 input catchment volumes time series,
and 1 output file path. With these files in place, you can now perform a routing computation. Use the following code
as a starting place.

```python
import river_route as rr

config_file_path = '/path/to/config.yaml'

m = (
    rr
    .MuskingumCunge(config_file_path)
    .route()
)
```

## Save and Plot Hydrographs

The default output format for the routed discharge is a netCDF file. You are free to write your own code or use any
compatible software to query data from that file. For quick access to a hydrograph for a single river, the routing class
has a `hydrograph` method to extract the hydrograph for a river number of interest and return it as a Pandas DataFrame.
From there, you can manipulate the data and plot is as normal.

!!! note
    The `hydrograph` method will only work if you are using the default output file format. If you have overridden the
    output file format, you will need to write your own function to extract the hydrograph.

```python
river_of_interest = 123456789
df = m.hydrograph(river_id=river_of_interest)

# Save the hydrograph to disc
df.to_csv('hydrograph.csv')

# Plot the hydrograph
df.plot()
```
