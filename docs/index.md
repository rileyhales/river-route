# River-Route

`river-route` routes runoff and discharge through large river networks using numba-accelerated Muskingum-family solvers.

## Describe your routing

There is one router, `rr.Router`. You describe the routing procedure with three config selector keys.

| Selector  | Options                                                                                            | Default |
|-----------|----------------------------------------------------------------------------------------------------|---------|
| `coeff`   | `'static'` (constant Muskingum K from columns `k`, `x`) or `'dynamic'` (nonlinear K = alpha*Q^beta from columns `alpha`, `beta`, `x`). | `'static'` |
| `forcing` | `'channel'` (channel routing only), `'lateral'` (route lateral inflow), or `'external'` (planned). One value only. | `'channel'` |
| `network` | `'standard'` (one reach per river). `'expanded'` (auto subdivide/substep unstable reaches) is planned, not yet available. | `'standard'` |

```python
import river_route as rr

rr.Router("/path/to/config.yaml", coeff="dynamic", forcing="lateral", network="standard").route()
```

`forcing='external'`, specifying multiple forcings, and Unit Hydrograph routing are planned for a later v3 release and are not yet available in v3.

## Quick Start

```bash
pip install river-route
```

```python
import river_route as rr

rr.Router("/path/to/config.yaml").route()
```

Config can be passed as:

1. A YAML/JSON file path.
2. Keyword arguments.
3. Both (keyword arguments override config file values).

```python
import river_route as rr

(
    rr
    .Router(
        "/path/to/config.yaml",
        forcing="lateral",
        qlateral_files=["/path/to/catchment_runoff.nc"],
        discharge_dir="/path/to/output/",
    )
    .route()
)
```

## Start Here

1. [Basics](tutorial/basics.md)
2. [Unit Hydrographs with Routing](tutorial/unit-hydrograph-routing.md) (Unit Hydrograph routing is planned for a later v3 release)
3. [Routing Ensembles](tutorial/routing-ensembles.md)
4. [Advanced Uses](tutorial/advanced.md)

## Core References

1. [Configuration File](references/config-files.md)
2. [Input/Output File Schemas](references/io-file-schema.md)
3. [Time Variables](references/time-options.md)
4. [Math Derivations](references/math.md)
5. [Parallelism](references/parallelism.md)
