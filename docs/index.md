# River-Route

The `river-route` Python package routes runoff volumes and discharge through river networks using
numba-accelerated solvers. It is designed for large networks of rivers and catchments.

## Routers

You have these options for how to do hydrological routing

### Channel routing only

| Router      | Description                                               |
|-------------|-----------------------------------------------------------|
| `Muskingum` | Matrix muskingum channel routing with no lateral inflows. |

### Channel routing with runoff lateral inflows

| Router           | Description                                                                                                                    |
|------------------|--------------------------------------------------------------------------------------------------------------------------------|
| `RapidMuskingum` | Runoff is treated as a point source in the catchment which all enters the channel during the interval the runoff occurs.       | 
| `UnitMuskingum`  | Runoff transformed via convolution with user provided unit hydrographs at each catchment then routed in the downstream segment |

## Start Here

1. [Muskingum Channel Routing](tutorial/basics.md)
2. [Channel Routing with Runoff Transformation](tutorial/unit-hydrograph-routing.md)
3. [Routing Ensembles](tutorial/routing-ensembles.md)
4. [Advanced Uses](tutorial/advanced.md)

```commandline
pip install river-route
```

```python
import river_route as rr

# Most common: route runoff volumes through the network
(
    rr
    .RapidMuskingum('/path/to/config.yaml')
    .route()
)
```

## Usage Examples

Configuration options are passed as a config file path, keyword arguments, or both. Keyword arguments
override any value in the config file.

```python
import river_route as rr

# Option 1 - Config file only
(
    rr
    .RapidMuskingum('/path/to/config.yaml')
    .route()
)

# Option 2 - Keyword arguments only
(
    rr
    .RapidMuskingum(**{
        'params_file': '/path/to/routing_params.parquet',
        'qlateral_files': '/path/to/catchment_runoff.nc',
        'discharge_dir': '/path/to/output/',
    })
    .route()
)

# Option 3 - Config file with keyword argument overrides
(
    rr
    .RapidMuskingum(
        '/path/to/config.yaml',
        **{
            'qlateral_files': '/path/to/catchment_runoff.nc',
            'discharge_dir': '/path/to/output/',
        }
    )
    .route()
)
```
