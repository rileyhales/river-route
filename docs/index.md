# River-Route

`river-route` routes runoff and discharge through large river networks using numba-accelerated Muskingum-family solvers.

## Choose a Router

| Router           | Description                                                      |
|------------------|------------------------------------------------------------------|
| `Muskingum`      | Channel routing only (no lateral inflow).                        |
| `RapidMuskingum` | Routes lateral runoff directly to channels each runoff interval. |
| `UnitMuskingum`  | Convolves runoff with a unit hydrograph, then routes in-channel. |

## Quick Start

```bash
pip install river-route
```

```python
import river_route as rr

(
    rr
    .RapidMuskingum("/path/to/config.yaml")
    .route()
)
```

Config can be passed as:

1. A YAML/JSON file path.
2. Keyword arguments.
3. Both (keyword arguments override config file values).

```python
import river_route as rr

(
    rr
    .RapidMuskingum(
        "/path/to/config.yaml",
        qlateral_files=["/path/to/catchment_runoff.nc"],
        discharge_dir="/path/to/output/",
    )
    .route()
)
```

## Start Here

1. [Basics](tutorial/basics.md)
2. [Unit Hydrographs with Routing](tutorial/unit-hydrograph-routing.md)
3. [Routing Ensembles](tutorial/routing-ensembles.md)
4. [Advanced Uses](tutorial/advanced.md)

## Core References

1. [Configuration File](references/config-files.md)
2. [Input/Output File Schemas](references/io-file-schema.md)
3. [Time Variables](references/time-options.md)
4. [Math Derivations](references/math.md)
5. [Forward Substitution](references/forward-substitution.md)
