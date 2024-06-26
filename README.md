# River Route

`river-route` is a Python package for routing runoff through a river network. Routing calculations are vectorized with
numpy and scipy which keeps the array computation times on par with faster compiled languages.

## Output File Schemas

### Discharge Data

Discharge data are given in a netCDF file. The file has 2 dimensions: time and river_id.
The times given will exactly match the times given in the runoff data unless a different time interval was specified.
The river_id dimension will exactly match the IDs given in the river_id column of the routing parameters file. Discharge
values will be in a variable named "Q", in an array of shape (time, river_id), and dtype float.



## Tips for efficient computations

Results will vary based on your system specifications, the size of the river network, and the length of the simulation.
These tips may help you achieve faster results.

1. **Use fewer inflow files**: File IO operations can be relatively slow and are a probable bottleneck on HPC systems
   when I/O operations depend on networked drives. You may achieve faster results by doing a single computation
   covering 2 weeks instead of 14 computations covering 1 day each. The trade-off is that larger arrays consume more
   memory which could eventually become too large if your river network is very large.
2. **Adjust the time step**: Using a longer time step will reduce the number of computations and therefore take less
   time to compute. It also requires storing fewer intermediate results in memory yielding a modest reduction in memory
   usage. A longer time step can increase performance if it does not induce numerical instability in the outflows.
