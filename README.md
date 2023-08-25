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

