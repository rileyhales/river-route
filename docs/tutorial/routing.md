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
