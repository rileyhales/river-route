!!! warning "Non Comprehensive Tutorial Content Disclaimer"
    This is not a comprehensive hydrological modeling course which should teach you the theory of hydrological channel 
    routing, calibration, and validation. It will not teach the prerequisite informatics and GIS skills to create a 
    watershed representation. Hydraulics and hydrology software can be used to obtain the information needed for routing 
    or you can use a GIS software such as QGIS or ArcGIS. An approachable place to start learning these skills is a GIS 
    tutorial demonstrating watershed delineation and stream network extraction as well as assigning and calculating 
    attributes of the stream features (polylines).

## Vocabulary

**VPU**: Vector Processing Unit. This refers to the collection of sub-basins or complete watersheds that you group
together for computations. The term VPU is not required by river-route or a technical definition in hydrology or
computational science. You can think of this as an ID or name of the watershed you are simulating.

You may simulate a sub-basin (partial watershed), a single complete watershed, or multiple complete watersheds. You
might want to have multiple computation units for at least 2 common reasons.

1. Multiple divisions of watersheds enable parallelization of computations.
2. You may have multiple versions of input datasets for your watershed. For instance, different stream network
   definitions coming from different elevation sources or calibrated routing parameters for different events.

## Essentials

`river-route` is a Python package that routes water volumes through a river network. It takes 3 input datasets and
writes
1 output dataset. Together, that makes a total of 4 files that need to be specified. The first 2 describe the river
channel properties and topology, the third is the input water being routed, and the fourth (the output) is the discharge
time series calculated by the routing process. For more information on how to prepare these, see the following sections.

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
│   │   ├── volume_<vpu-name>_<first-time-steps-or-ensemble-member>.nc
│   │   └── volume_<vpu-name>_<next-time-steps-or-ensemble-member>.nc
│   └── <VPU Name 2>
│       ├── volume_<vpu-name>_<first-time-steps-or-ensemble-member>.nc
│       └── volume_<vpu-name>_<next-time-steps-or-ensemble-member>.nc
└── discharge
    ├── <VPU Name 1>
    │   ├── discharge_<vpu-name>_<first-time-steps-or-ensemble-member>.nc
    │   └── discharge_<vpu-name>_<next-time-steps-or-ensemble-member>.nc
    └── <VPU Name 2>
        ├── discharge_<vpu-name>_<first-time-steps-or-ensemble-member>.nc
        └── discharge_<vpu-name>_<next-time-steps-or-ensemble-member>.nc
```

## Preparing Watershed Inputs

## Preparing Volumes

...

## Preparing a Config File

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

## First Routing Run

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

## Saving and Plotting Hydrographs

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
