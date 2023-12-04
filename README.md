# River Route

`river-route` is a Python package for routing runoff through a river network.
The routing calculations are vectorized and use numpy and scipy which keeps the array computation times on par with
faster compiled languages.

The implemented routing methods are:

- Matrix Form Muskingum Cunge - Analytical solution
- Matrix Form Muskingum Cunge - Numerical solution

## Quick Start Guide
You will need to prepare a configuration file for the routing.

```python
import river_route as rr

(
    rr
    .Muskingum('/path/to/config.yml')
    .route()
)
```
## Configuration File

The minimum required inputs in the configuration file are:

- `routing_params_file` - path to the routing parameters file (parquet)
- `connectivity_file` - path to the river network connectivity file (parquet)
- `runoff_file` - path to the prepared runoff data file (netCDF)
- `outflow_file` - path where the routed flows will be saved (netCDF)

It is recommended that you also provide these parameters:

- `lhs_file` - path where the LHS matrix will be cached (npz file). This parameter can be provided to future computations for a speed gain.
- `lhsinv_file` - path where the LHS inverse matrix will be cached (npz file). This parameter can be provided to future computations for a substantial speed gain.
- `adj_file` - path where the adjacency matrix will be cached (npz file). This parameter can be given instead of the connectivity file in future computations for a speed gain. 
- `dt_routing` - an integer time in seconds for the routing time step. Will default to 300s or 5 minutes
- `dt_outflows` - an integer time in seconds for the outflow time step. Will default to the inflow time step

You can provide initial conditions/state and save final conditions/state with these parameters:

- `qinit_file` - path to the initial flows file (parquet). Defaults to 0.0 for all rivers.
- `rinit_file` - path to the initial runoff file (netCDF). Defaults to 0.0 for all rivers.
- `qfinal_file` - path where the final flows file will be saved (parquet). It will not be saved if a path is not provided.
- `rfinal_file` - path where the final runoff file will be saved (netCDF). It will not be saved if a path is not provided.

You can provide logging options with these parameters:

- `log_file` - path to a file where the log will be written. Defaults to stdout (logs not saved to file).
- `job_name` - a name for this job printed in logs and debug statements.
- `progress_bar` - display a computations progress bar in logs: True or False. Defaults to True.

You can modify how the routing computations are performed with these parameters:

- `routing_method` - method for solving routing equations: 'numerical' or 'analytical'. Defaults to 'numerical'.
- `petsc_ksp_type` - name of the PETSC solver used when using the numerical solution. Recommend 'preonly' if you have a well-prepared network and 'richardson' otherwise.


A diagram of the possible configuration file parameters and their role in the routing computations is shown below.

```mermaid
graph LR
subgraph "Required Inputs"
Inputs["
Routing Parameters
Connectivity or Adjacency
Runoff Data
Outflow File
"]
end

subgraph "Time Options"
TimeParameters["
Routing Timestep
Outflow Timestep
"]
end

subgraph "Init Files"
Initialization["
Initial Flows
Initial Runoff
"]
end

Routing[Perform Routing]

subgraph "Main Output"
Result[Routed Discharge]
end

subgraph "Cached Files"
CachedFiles["
Final Flows State
Final Runoff State
LHS Matrix
LHS Inverse Matrix
Adjacency Matrix
"]
end

Inputs & TimeParameters & Initialization ==> Routing
Routing ==> Result & CachedFiles
```

| Parameter Name      | Required | Type      | Group                  | Help                                                               |                                                                                
|---------------------|----------|-----------|------------------------|--------------------------------------------------------------------|
| routing_params_file | True     | File Path | Stream Network Details | Path to the routing parameters parquet file.                       |                                                
| connectivity_file   | True     | File Path | Stream Network Details | Path to the network connectivity parquet file.                     |                                              
| runoff_file         | True     | File Path | Input Runoff           | Path to the netCDF with runoff values to be routed.                |                                   
| outflow_file        | True     | File Path | Output Discharge       | Path where the outflows netCDF file should be saved.               |                                       
| dt_routing          | True     | Integer   | Time Parameter         | Time interval in seconds between routing computations.             |                              
| dt_outflows         | False    | Integer   | Time Parameter         | Time interval in seconds between writing flows to disc.            |                             
| qinit_file          | False    | File Path | Initialization Data    | Path to the initial flows file.                                    |                                                     
| rinit_file          | False    | File Path | Initialization Data    | Path to the initial runoff file.                                   |                                                    
| qfinal_file         | False    | File Path | Output File Path       | Path where the final flows file should be saved.                   |                                    
| rfinal_file         | False    | File Path | Output File Path       | Path where the final runoff file should be saved.                  |                                   
| lhs_file            | False    | File Path | Cachable Network File  | Path where the LHS matrix should be cached.                        |                                         
| lhsinv_file         | False    | File Path | Cachable Network File  | Path where the LHS inverse matrix should be cached.                |                                 
| adj_file            | False    | File Path | Cachable Network File  | Path where the adjacency matrix should be cached.                  |                                   
| log_file            | False    | File Path | Logging Options        | Path to a file where the log will be written.                      |                                       
| job_name            | False    | String    | Logging Options        | A name for this job printed in logs and debug statements.          |                           
| progress_bar        | False    | Boolean   | Logging Options        | Display a computations progress bar in logs: true or false.        | 
| routing_method      | False    | String    | Solving Options        | Method for solving routing equations: 'numerical' or 'analytical'. |       
| petsc_ksp_type      | False    | String    | Solving Options        | Name of the PETSC solver used when using the numerical solution.   | 

## Time Parameters

### Parameter Descriptions

Only 1 time option is a required input in the configuration file:

- `dt_routing` - the time interval, in seconds, between routing calculation steps. It must be constant across all rivers
  and for the full simulation.

3 other time parameters are optional. They may be provided in the configuration file or they will be derived from the
runoff data.

- `dt_runoff` - the time interval, in seconds, between runoff values. It must be constant between all time steps of
  runoff.
- `dt_outflow` - the time interval, in seconds, between outflow values which get written to disc. It must be constant
  between all time steps of outflow.
- `dt_total` - the total time, in seconds, of the runoff data. It is equal to the number of time steps multiplied
  by `dt_runoff`.

### Parameter Relationships

```
dt_total >= dt_outflow >= dt_runoff >= dt_routing
```

- The total time of the runoff data must be greater than or equal to the time interval between outflow values.
- The time interval between outflow values must be greater than or equal to the time interval between runoff values.
- The time interval between runoff values must be greater than or equal to the time interval between routing
  calculations.

```
dt_total % dt_outflow == 0
dt_outflow % dt_runoff == 0
dt_runoff % dt_routing == 0
```

- Each time interval must be evenly divisible by the next smallest time interval so that the loops of calculations can
  be automatically constructed.

```
dt_total === dt_runoff * number_of_time_steps
```

- The length of the routing computations must be equal to the length of the runoff data.

### Limitations

1. You may not write to disc more frequently than the runoff interval.
2. You must route each interval of runoff data at least 1 time before the next interval of runoff data. You should
   resample the runoff data to larger or smaller time intervals using your own methods before routing.
3. You cannot run the routing calculations for longer than runoff data are available. If the runoff data cover 10 days,
   the routing computations will only last 10 days. Routing longer can be accomplished by padding the runoff data with
   zero values for more time steps.

## Tips for efficient computations

Results will vary based on your system specifications, the size of the river network, and the length of the simulation.
These tips may help you achieve faster results.

1. **Use fewer inflow files**: File IO operations can be relatively slow and are a probable bottleneck on HPC systems
   when I/O operations depend on networked drives. You may achieve faster results by doing a single computation
   covering 2 weeks instead of 14 computations covering 1 day each.
2. **Cache routing Matrices**: The adjacency matrix and inverted I-C2@A matrix can be time consuming to compute. Provide
   paths to store them in the config file to cache them between simulations
3. **Adjust the time step**: Using a longer time step will reduce the number of computations which takes less time to
   compute. It also requires storing fewer intermediate results in memory yielding a modest reduction in memory usage. A
   longer time step can increase performance if it does not induce numerical instability in the outflows.

