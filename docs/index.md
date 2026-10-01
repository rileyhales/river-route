# River-Route

`river-route` routes runoff and discharge through large river networks using numba-accelerated Muskingum-family
solvers that route each river's whole time series in turn, from upstream to downstream.

## Describe your routing

All routing runs through `rr.Router`. You describe the routing procedure with config selector keys,
which together choose the kernel.

| Selector       | Options                                                                                                                                                            | Default      |
|----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------|
| `coefficients` | `'static'` (constant Muskingum K from columns `muskingumK`, `muskingumX`) or `'dynamic'` (nonlinear K = dynamicAlpha*Q^dynamicBeta from columns `dynamicAlpha`, `dynamicBeta`, `muskingumX`). | `'static'`   |
| `forcing`      | `'channel'` (channel routing only), or the form of the `runoff_files` routed into the rivers: `'catchment'`, `'grid'`, or `'ecmwf_grib'`.                          | `'channel'`  |
| `transform`    | `'uniform'`, the only option: each step's catchment runoff enters its river at a constant rate over the step.                                                      | `'uniform'`  |
| `network_type` | `'standard'` (one reach per river) or `'stabilized'` (each river routed in the substeps or subcycles that make it stable at `dt_routing`, where any do).           | `'standard'` |

The one combination no routing method routes yet, `'dynamic'` coefficients on a `'stabilized'` network, raises
`NotImplementedError` before any runoff is read.

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
    network_file="/path/to/network.parquet",
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
| `Network` | The river network: ids, topology, `k` and `x`, the concurrent partition, stability analysis, substeps and subcycles. |
| `Router`  | One simulation over a `Network`: coefficients, time options, channel state, the routing loop.            |
| `GridRunoff` | Gridded runoff to catchment runoff, reusing one weight table across many runoff files.              |

`Router` takes only a `Configs`, and builds its `Network` and the `Runoff` class its `forcing` names from it.
`Network` and `GridRunoff` take the options they need as ordinary arguments and each has a `from_configs` classmethod
that reads those same values off a `Configs`. `examples/config.json` lists every option.

A `Router` reads and partitions its network table once, and every `route()` call on it reuses the result.

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
