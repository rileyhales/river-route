!!! warning "Non Comprehensive Tutorial Content Disclaimer"
    This is not a comprehensive hydrological modeling course which should teach you the theory of hydrological channel routing, calibration, and 
    validation. It will not teach the prerequisite informatics and GIS skills to create a watershed representation. Hydraulics and hydrology software 
    can be used to obtain the information needed for routing or you can use a GIS software such as QGIS or ArcGIS. An approachable place to start 
    learning these skills is a GIS tutorial demonstrating watershed delineation and stream network extraction as well as assigning and calculating 
    attributes of the stream features.

## Essentials

`river-route` is a Python package that routes water volumes through a river network. It takes 3 input datasets and writes 1 output dataset. Together, 
that makes a total of 4 files that need to be specified. The first 2 describe the river channel properties and topology, the third is the input water 
being routed, and the fourth (the output) is the discharge time series calculated by the routing process. This tutorial explain the process of 
preparing these files.

1. [Routing Parameters](../references/io-file-schema.md#routing-parameters) - parquet file
2. [Connectivity File](../references/io-file-schema.md#connectivity-file) - parquet file
3. [Catchment Volumes](../references/io-file-schema.md#catchment-volumes-or-runoff-depths) - netCDF file
4. [Routed Discharge](../references/io-file-schema.md#routed-discharge) - netCDF file

## Recommended Directory Structure for Organizing Inputs and Outputs

***You are not required to use a specific file structure or naming conventions***. If you do
not already have a strong preference or existing structure to match, you could consider using the following pattern:

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

## Vocabulary

The following vocabulary is used in this tutorial.

**VPU**: Vector Processing Unit. This refers to the collection of sub-basins or complete watersheds that you group
together for computations. The term VPU is not required by river-route or a technical definition in hydrology or
computational science. You can think of this as an ID or name of the watershed you are simulating.
