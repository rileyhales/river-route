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

- `routing_params_file` - path to the [routing parameters file](io-file-schema.md#routing-parameters) (parquet)
- `connectivity_file` - path to the river network [connectivity file](io-file-schema.md#connectivity-file) (parquet)
- `catchment_volumes_file` - path to the prepared [catchment volumes file](io-file-schema.md#catchment-volumes-or-runoff-depths) (netCDF)
- `outflow_file` - path where the [routed flows](io-file-schema.md#routed-discharge) output file will be saved (netCDF)

## Example YAML File

An example yaml file is given below with the default values prepopulated where possible.

```yaml title="Config File Example river-route v0.15.0"
{% include-markdown "../../examples/config.yaml" %}
```

## Config Options Table

The following table is a complete list of all configuration options and their purpose.

| Parameter Name         | Required | Type      | Group               | Description                                                            |                                                                                
|------------------------|----------|-----------|---------------------|------------------------------------------------------------------------|
| routing_params_file    | True     | File Path | Required Input      | Path to the routing parameters parquet file.                           |                                                
| connectivity_file      | True     | File Path | Required Input      | Path to the network connectivity parquet file.                         |                                              
| catchment_volumes_file | True     | File Path | Required Input      | Path to the netCDF with catchment volume values to be routed.          |
| runoff_depths_files    | True     | File Path | Required Input      | List of paths to netCDF files with runoff depths to be routed.         |
| weight_table_file      | True     | File Path | Required Input      | Path to the weight table file.                                         |
| outflow_file           | True     | File Path | Required Input      | Path where the outflows netCDF file should be saved.                   |
| routing                | False    | String    | Compute Options     | The routing method to use: "linear" or "nonlinear".                    |
| input_type             | False    | String    | Compute Options     | Specify if runoff files are "sequential" time steps or an "ensemble"   |
| runoff_type            | False    | String    | Compute Options     | Specify if runoff files are "incremental" or "cumulative"              |
| dt_routing             | True     | Integer   | Compute Options     | Time interval in seconds between routing computations.                 |                              
| dt_outflows            | False    | Integer   | Compute Options     | Time interval in seconds between writing flows to disc.                |
| solver_atol            | False    | Float     | Compute Options     | Absolute tolerance for the solver.                                     |
| initial_state_file     | False    | File Path | Initialization Data | Path to the initial state file.                                        |                                                     
| final_state_file       | False    | File Path | Initialization Data | Path where the final state file should be saved.                       |                                    
| log                    | False    | Boolean   | Logging Options     | Whether to enable logging.                                             |                                       
| progress_bar           | False    | Boolean   | Logging Options     | Whether to display a progress bar when routing                         |
| log_level              | False    | String    | Logging Options     | The logging level to print: DEBUG, INFO, CRITICAL, WARNING             |
| log_stream             | False    | String    | Logging Options     | The destination for log messages: 'stdout', 'stderr', or a file path.  |
| var_x                  | False    | String    | File Management     | Name of the variable in all files that contains the x coordinates.     |
| var_y                  | False    | String    | File Management     | Name of the variable in all files that contains the y coordinates.     |
| var_t                  | False    | String    | File Management     | Name of the variable in all files that contains the time values.       |
| var_runoff_depths      | False    | String    | File Management     | Name of the variable in files containing runoff depths                 |
| var_catchment_volumes  | False    | String    | File Management     | Name of the variable in the catchment volumes file containing volumes. |
| var_river_id           | False    | String    | File Management     | Name of the variable in all files that contains the river IDs.         |
| var_outflow            | False    | String    | File Management     | Name of the variable in the outflows file that contains the outflows.  |
