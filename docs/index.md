# River-Route

The `river-route` Python package routes runoff volumes and discharge through river networks using vectorized,
sparse-matrix solvers. It is designed for large networks of rivers and catchments.

## Routers

Three routers are available:

| Router              | Description                                                                                                   |
|---------------------|---------------------------------------------------------------------------------------------------------------|
| `Muskingum`         | Channel routing with no lateral inflows.                                                                      |
| `TeleportMuskingum` | Overland runoff placed at catchment inlets uniformly over the runoff timestep. Muskingum channel routing.     |
| `UnitMuskingum`     | Overland runoff transformed by a pluggable unit hydrograph before channel routing. Muskingum channel routing. |

## Start Here

1. [Basic Walkthrough](tutorial/basic-tutorial.md)
2. [Advanced Concepts](tutorial/advanced-tutorial.md)
3. [Routing Ensembles](tutorial/routing-ensembles.md)
4. [Unit Hydrograph Routing](tutorial/unit-hydrograph-routing.md)

```commandline
pip install river-route
```

```python
import river_route as rr

# Most common: route runoff volumes through the network
(
    rr
    .TeleportMuskingum('/path/to/config.yaml')
    .route()
)
```

## Computation Process

```mermaid
---
title: River Route Process Diagram
---
graph LR
    subgraph "Required-Inputs"
        Inputs["Routing Parameters\nCatchment Volumes (or Runoff Depths)\nDischarge Files"]
    end

    subgraph "Compute-Options"
        co1["Routing Timestep\nDischarge Timestep\nRunoff Type\nInput Type"]
    end

    subgraph "Initialization"
        i1["Initial State File"]
    end

    subgraph "Logging-Options"
        management["Log File Path\nProgress Bar\nLog Level"]
    end

    subgraph "Computations"
    direction TB
    a[Build Adjacency Matrix] --> b
    b[Factorize LHS Matrix] --> c
    c[Read Volume/Depth Array\nTeleportMuskingum & UnitMuskingum] --> d
    d[Apply UH Transform\nUnitMuskingum only] --> e
    e[Iterate Routing Timesteps] --> f
    f[Solve Muskingum System] --> g
    g[Enforce Non-negative Flows] --> h & e
    h[Write Discharge to Disk] --> i
    i[Cache Final State]
    end

    subgraph "Main-Output"
    Result[Routed Discharge]
    end

    subgraph "Cachable-Files"
    CachedFiles["Final State\nUH Kernel & State\n(UnitMuskingum)"]
    end

    Required-Inputs & Compute-Options & Initialization & Logging-Options ==> Computations
    Computations ==> Main-Output & Cachable-Files
```

## UnitMuskingum Transformer

`UnitMuskingum` delegates runoff transformation to a pluggable `AbstractBaseTransformer` subclass.
Two initialization paths are available:

```mermaid
---
title: UnitMuskingum Transformer Initialization
---
graph TD
    A{uh_kernel\nprovided?} -->|Yes| B[Transformer.from_kernel]
    A -->|No| C{uh_type}
    C -->|'scs'| D[SCSUnitHydrograph.__init__]
    C -->|custom| E[set_transformer injection]
    B --> F{uh_state\nprovided?}
    D --> F
    E --> F
    F -->|Yes| G[set_state from parquet]
    F -->|No| H[zero initial state]
    G & H --> I[AbstractBaseTransformer ready]
```

```mermaid
---
title: Transformer Class Hierarchy
---
classDiagram
    class AbstractBaseTransformer {
        <<abstract>>
        +float dt
        +FloatArray kernel
        +FloatArray state
        +from_kernel(dt, kernel_path)$
        +set_state(state_path)
        +transform(runoff_vector) FloatArray
        +_build_kernel()* FloatArray
    }
    class Transformer {
        +_build_kernel() NotImplementedError
    }
    class SCSUnitHydrograph {
        +FloatArray tc
        +FloatArray area
        +_build_kernel() FloatArray
    }
    AbstractBaseTransformer <|-- Transformer : load from parquet
    AbstractBaseTransformer <|-- SCSUnitHydrograph : compute from params
```

## Usage Examples

Configuration options are passed as a config file path, keyword arguments, or both. Keyword arguments
override any value in the config file.

```python
import river_route as rr

# Option 1 - Config file only
(
    rr
    .TeleportMuskingum('/path/to/config.yaml')
    .route()
)

# Option 2 - Keyword arguments only
(
    rr
    .TeleportMuskingum(**{
        'routing_params_file': '/path/to/routing_params.parquet',
        'catchment_volumes_files': '/path/to/volumes.nc',
        'discharge_files': '/path/to/discharge.nc',
    })
    .route()
)

# Option 3 - Config file with keyword argument overrides
(
    rr
    .TeleportMuskingum(
        '/path/to/config.yaml',
        **{
            'catchment_volumes_files': '/path/to/volumes.nc',
            'discharge_files': '/path/to/discharge.nc',
        }
    )
    .route()
)

# Unit-Muskingum (with SCS triangular unit hydrograph)
(
    rr
    .UnitMuskingum('/path/to/config.yaml')
    .route()
)
```
