# River Route

`river-route` is a Python package for routing runoff through a river network.
The routing calculations are vectorized and use numpy and scipy which keeps the array computation times on par with faster compiled languages. 

The implemented routing methods are:
- Muskingum Cunge - Analytical solution
- Muskingum Cunge - Numerical solution

## Quick Start Guide
You will need to prepare a configuration file for the routing.

```python
import river_route as rr

rm = rr.Muskingum('/path/to/config.yml')
rm.route()

river_id_to_inspect = 12345
rm.plot(river_id_to_inspect)
```

The minimum required inputs in the configuration file are:
 - `runoff_file` - path to the prepared runoff data file (netCDF)
 - `routing_params_file` - path to the routing parameters file (parquet)
 - `connectivity_file` - path to the river network connectivity file (parquet)
 - `outflow_file` - path where the routed flows will be saved (netCDF)
 - `dt_routing` - an integer time in seconds for the routing time step (e.g. 300 for 5 minutes)


## Time Parameters

### Parameter Descriptions

Only 1 time option is a required input in the configuration file:
 - `dt_routing` - the time interval, in seconds, between routing calculation steps. It must be constant across all rivers and for the full simulation.

3 other time parameters are optional. They may be provided in the configuration file or they will be derived from the runoff data.
- `dt_runoff` - the time interval, in seconds, between runoff values. It must be constant between all time steps of runoff.
- `dt_outflow` - the time interval, in seconds, between outflow values which get written to disc. It must be constant between all time steps of outflow.
- `dt_total` - the total time, in seconds, of the runoff data. It is equal to the number of time steps multiplied by `dt_runoff`.

### Parameter Relationships

```
dt_total >= dt_outflow >= dt_runoff >= dt_routing
```
- The total time of the runoff data must be greater than or equal to the time interval between outflow values.
- The time interval between outflow values must be greater than or equal to the time interval between runoff values.
- The time interval between runoff values must be greater than or equal to the time interval between routing calculations.

```
dt_total % dt_outflow == 0
dt_outflow % dt_runoff == 0
dt_runoff % dt_routing == 0
```
- Each time interval must be evenly divisible by the next smallest time interval so that the loops of calculations can be automatically constructed.

```
dt_total === dt_runoff * number_of_time_steps
```
- The length of the routing computations must be equal to the length of the runoff data.

### Limitations
1. You may not write to disc more frequently than the runoff interval.
2. You must route each interval of runoff data at least 1 time before the next interval of runoff data. You should resample the runoff data to larger or smaller time intervals using your own methods before routing.
3. You cannot run the routing calculations for longer than runoff data are available. If the runoff data cover 10 days, the routing computations will only last 10 days. Routing longer can be accomplished by padding the runoff data with zero values for more time steps.

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

