# River-Route

`river-route` routes runoff and discharge through large river networks using numba-accelerated Muskingum-family solvers.

## Describe your routing

All routing runs through `rr.Router`. You describe the routing procedure with config selector keys,
which together choose the kernel.

| Selector       | Options                                                                                                                                                            | Default      |
|----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------|
| `coefficients` | `'static'` (constant Muskingum K from columns `k`, `x`) or `'dynamic'` (nonlinear K = dynamicAlpha*Q^dynamicBeta from columns `dynamicAlpha`, `dynamicBeta`, `x`). | `'static'`   |
| `forcing`      | `'channel'` (channel routing only), or the form of the `runoff_files` routed into the rivers: `'catchment'`, `'grid'`, or `'ecmwf_grib'`.                          | `'channel'`  |
| `transform`    | `'uniform'`, the only option. Only read when `forcing` is not `'channel'`.                                                                                         | `'uniform'`  |
| `network_type` | `'standard'` (one reach per river) or `'stabilized'` (every river routed stably at `dt_routing` with sub-reaches or substeps).                                     | `'standard'` |

A combination with no kernel yet raises `NotImplementedError` naming it and listing those that exist.

```python
import river_route as rr

configs = rr.Configs.from_json("/path/to/config.json")  # sets coefficients, forcing, network_type, and the other options
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
    forcing="catchment",
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
| `GridRunoff` | Gridded runoff to lateral inflow, reusing one weight table across many runoff files.                 |

`Router` takes a `Configs` and nothing else, and builds its `Network` and `GridRunoff` from it. `Network` and `GridRunoff`
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
