## Configuration File

`river-route` computations are controlled by several variables. You can pass these variables as keyword arguments to the
corresponding functions or provide a path to a configuration file. Supported file formats for config files are YAML or
JSON. Config files specify the following parameters:

1. Paths to input and output files
2. Options to control how the routing is performed
3. Options to control the numerical solver and tolerances
4. Information about the structure of the input and output files if non-standard
5. Logging parameters

## Minimum Required Inputs

Every `river-route` process needs at least the following 4 variables

- `routing_params_file` - path to the [routing parameters file](io-files.md#routing-parameters) (parquet)
- `connectivity_file` - path to the river network [connectivity file](io-files.md#connectivity-file) (parquet)
- `catchment_volumes_file` - path to the prepared [catchment volumes file](io-files.md#catchment-volumes) (netCDF)
- `outflow_file` - path where the [routed flows](io-files.md#routed-discharge) output file will be saved (netCDF)

## Config Options Table

The following table is a complete list of all configuration options and their purpose.

| Parameter Name         | Required | Type      | Group               | Description                                                           |                                                                                
|------------------------|----------|-----------|---------------------|-----------------------------------------------------------------------|
| routing_params_file    | True     | File Path | Required Input      | Path to the routing parameters parquet file.                          |                                                
| connectivity_file      | True     | File Path | Required Input      | Path to the network connectivity parquet file.                        |                                              
| catchment_volumes_file | True     | File Path | Required Input      | Path to the netCDF with catchment volume values to be routed.         |
| volume_file_vol_var    | False    | String    | File Management     | Name of the variable in the volume files.                             |
| runoff_depths_files    | True     | List      | Required Input      | List of paths to netCDF files with runoff depths to be routed.        |
| weight_table_file      | False    | File Path | Required Input      | Path to the weight table file.                                        |
| depth_file_x_var       | False    | String    | File Management     | Name of the x variable in the depth files.                            |
| depth_file_y_var       | False    | String    | File Management     | Name of the y variable in the depth files.                            |
| depth_file_t_var       | False    | String    | File Management     | Name of the time variable in the depth files.                         |
| depth_file_runoff_var  | False    | String    | File Management     | Name of the runoff variable in the depth files.                       |
| outflow_file           | True     | File Path | File Management     | Path where the outflows netCDF file should be saved.                  |                                       
| dt_routing             | True     | Integer   | Compute Options     | Time interval in seconds between routing computations.                |                              
| dt_outflows            | False    | Integer   | Compute Options     | Time interval in seconds between writing flows to disc.               |
| routing                | False    | String    | Compute Options     | The routing method to use: "linear" or "nonlinear".                   |
| runoff_type            | False    | String    | Compute Options     | Specify if runoff files are "sequential" time steps or an "ensemble"  |
| initial_state_file     | False    | File Path | Initialization Data | Path to the initial state file.                                       |                                                     
| final_state_file       | False    | File Path | Initialization Data | Path where the final state file should be saved.                      |                                    
| log                    | False    | Boolean   | Logging Options     | Whether to enable logging.                                            |                                       
| progress_bar           | False    | Boolean   | Logging Options     | Whether to display a progress bar when routing                        |
| log_level              | False    | String    | Logging Options     | The logging level to print: DEBUG, INFO, CRITICAL, WARNING            |
| log_stream             | False    | String    | Logging Options     | The destination for log messages: 'stdout', 'stderr', or a file path. |
| job_name               | False    | String    | Logging Options     | A name for this job printed in logs and debug statements.             |                           
| var_runoff_volume      | False    | String    | File Management     | Name of the variable in files containing runoff volumes               |
| var_river_id           | False    | String    | File Management     | Name of the variable in all files that contains the river IDs.        |
| var_outflow            | False    | String    | File Management     | Name of the variable in the outflows file that contains the outflows. |
| solver_atol            | False    | Float     | Solver Options      | Absolute tolerance for the solver.                                    |

## Example YAML File

An example yaml file is given below with the default values prepopulated where possible.

# specify the title and the language of the code block is YAML

```yaml title="Example Config YAML"
# Required Watershed Files
routing_params_file: ''
connectivity_file: ''
# Volume Inputs - (Option 1)
catchment_volumes_file: ''
var_runoff_volume: 'volume'
# Depth Inputs - (Option 2)
runoff_depths_files: ''
weight_table_file: ''
depth_file_x_var: 'lon'
depth_file_y_var: 'lat'
depth_file_t_var: 'time'
depth_file_runoff_var: 'ro'
# Output file
outflow_file: ''
# Input and Output file structure - Optional
var_river_id: 'river_id'
var_outflow: 'Q'
# Compute Options - Optional
routing: 'linear'  # linear or nonlinear
runoff_type: 'sequential'  # sequential or ensemble
dt_routing: 0  # defaults to time step of volume inputs
dt_outflows: 0  # defaults to time step of volume inputs
# Solver Options - Optional
solver_atol: 0.00001
# initial and final state files - Optional
initial_state_file: ''
final_state_file: ''
# simulation management and debugging - Optional
log: False
progress_bar: False
log_level: 'DEBUG'
log_stream: ''
job_name: ''
```
