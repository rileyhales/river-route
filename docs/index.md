# River-Route

`river-route` routes runoff and discharge through large river networks using numba-accelerated Muskingum-family solvers.

## Describe your routing

All routing runs through `rr.Router`. You describe the routing procedure with config selector keys,
which together choose the kernel.

| Selector      | Options                                                                                                                                | Default      |
|---------------|----------------------------------------------------------------------------------------------------------------------------------------|--------------|
| `coeff`       | `'static'` (constant Muskingum K from columns `k`, `x`) or `'dynamic'` (nonlinear K = alpha*Q^beta from columns `alpha`, `beta`, `x`). | `'static'`   |
| `forcing`     | `'channel'` (channel routing only) or `'runoff'` (runoff enters the rivers in addition to routing). One value only.                    | `'channel'`  |
| `transform`   | `'uniform'` or `'unit_hydrograph'`. Only read when `forcing` is `'runoff'`.                                                            | `'uniform'`  |
| `runoff_type` | `'catchment'`, `'gaussian_grid'`, or `'reduced_gaussian_grid'`, the form of `runoff_files`. Required when `forcing` is `'runoff'`.     | none         |
| `network`     | `'standard'` (one reach per river).                                                                                                    | `'standard'` |

A combination with no kernel yet raises `NotImplementedError` naming it and listing those that exist.

```python
import river_route as rr

configs = rr.Configs.from_json("/path/to/config.json")  # sets coeff, forcing, network, and the other options
rr.Router(configs).route()
```

## Quick Start

```bash
pip install river-route
```

```python
import river_route as rr

configs = rr.Configs.from_json("/path/to/config.json")
rr.Router(configs).route()
```

Options are held by a frozen `Configs` object, built from keyword arguments or read from a JSON file with
`Configs.from_json`, and validated when they are used. `Configs` is frozen and set once when it is built, so an
option is changed by building the `Configs` you want. `to_json` writes the options to a file that
`Configs.from_json` reads back.

```python
import river_route as rr

configs = rr.Configs(
    params_file="/path/to/params.parquet",
    forcing="runoff",
    runoff_type="catchment",
    runoff_files=["/path/to/catchment_runoff.nc"],
    discharge_dir="/path/to/output/",
)
configs.to_json("/path/to/config.json")
rr.Router(configs).route()
```

## Classes

| Class     | Owns                                                                                                   |
|-----------|--------------------------------------------------------------------------------------------------------|
| `Configs` | Every option, frozen. Built once and given to the classes below.                                        |
| `Network` | The river network: ids, topology, `k` and `x`, the concurrent partition, stability analysis, subdivision. |
| `Router`  | One simulation over a `Network`: coefficients, time options, channel state, the routing loop.            |
| `GaussianGridRunoff`  | Gridded runoff to lateral inflow, reusing one weight table across many runoff files.                     |

`Router` takes a `Configs` and nothing else, and builds its `Network` and `GaussianGridRunoff` from it. `Network` and `GaussianGridRunoff`
take the options they need as ordinary arguments and each has a `from_configs` classmethod that reads those same
values off a `Configs`. `examples/config.json` lists every option.

A `Network` reads and partitions a parameter table once and every simulation over it reuses the result, so it can
be built up front and handed to a `Router` when many runs share one network:

```python
import river_route as rr

configs = rr.Configs.from_json('/path/to/config.json')

network = rr.Network.from_configs(configs)
stabilized = network.stabilize(3600)         # in-memory network with reaches added so every one routes stably

router = rr.Router(configs)
router.network = network                     # optional; the Router builds its own when not given one
router.route()
```

## Start Here

1. [Basics](tutorial/basics.md)
2. [Routing Ensembles](tutorial/routing-ensembles.md)
3. [Advanced Uses](tutorial/advanced.md)

## Core References

1. [Configuration File](references/config-files.md)
2. [Input/Output File Schemas](references/io-file-schema.md)
3. [Time Variables](references/time-options.md)
4. [Math Derivations](references/math.md)
5. [Parallelism](references/parallelism.md)
