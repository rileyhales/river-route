# River-Route

`river-route` routes runoff and discharge through large river networks using numba-accelerated Muskingum-family solvers.

## Describe your routing

All routing runs through `rr.Router`. You describe the routing procedure with four config selector keys,
which together resolve to one compiled kernel.

| Selector    | Options                                                                                                                                | Default      |
|-------------|----------------------------------------------------------------------------------------------------------------------------------------|--------------|
| `coeff`     | `'static'` (constant Muskingum K from columns `k`, `x`) or `'dynamic'` (nonlinear K = alpha*Q^beta from columns `alpha`, `beta`, `x`). | `'static'`   |
| `forcing`   | `'channel'` (channel routing only) or `'vlateral'` (route lateral inflow). One value only.                                             | `'channel'`  |
| `transform` | `'uniform'` or `'unit_hydrograph'`. Only read when `forcing` is `'vlateral'`.                                                          | `'uniform'`  |
| `network`   | `'standard'` (one reach per river).                                                                                                    | `'standard'` |

An unimplemented combination raises `NotImplementedError` listing the ones that exist.

```python
import river_route as rr

configs = rr.Configs.from_file("/path/to/config.yaml")  # sets coeff, forcing, network, and the other options
rr.Router(configs).route()
```

## Quick Start

```bash
pip install river-route
```

```python
import river_route as rr

configs = rr.Configs.from_file("/path/to/config.yaml")
rr.Router(configs).route()
```

Options are held by a frozen `Configs` object, built from keyword arguments or read from a YAML/JSON file with
`Configs.from_file`, and validated when they are used. `Configs` is frozen and set once when it is built, so an
option is changed by building the `Configs` you want. `to_yaml` and `to_json` write the options to a file that
`Configs.from_file` reads back.

```python
import river_route as rr

configs = rr.Configs(
    params_file="/path/to/params.parquet",
    forcing="vlateral",
    vlateral_files=["/path/to/catchment_runoff.nc"],
    discharge_dir="/path/to/output/",
)
configs.to_yaml("/path/to/config.yaml")
rr.Router(configs).route()
```

## Classes

| Class     | Owns                                                                                                   |
|-----------|--------------------------------------------------------------------------------------------------------|
| `Configs` | Every option, frozen. Built once and given to the classes below.                                        |
| `Network` | The river network: ids, topology, `k` and `x`, the concurrent partition, stability analysis, subdivision. |
| `Router`  | One simulation over a `Network`: coefficients, time options, channel state, the routing loop.            |
| `RunoffGaussianGrid`  | Gridded runoff to lateral inflow, reusing one weight table across many runoff files.                     |

`Router` takes a `Configs` and nothing else, and builds its `Network` and `RunoffGaussianGrid` from it. `Network` and `RunoffGaussianGrid`
take the options they need as ordinary arguments and each has a `from_configs` classmethod that reads those same
values off a `Configs`. `examples/config.yaml` groups every option under the class that reads it.

A `Network` reads and partitions a parameter table once and every simulation over it reuses the result, so it can
be built up front and handed to a `Router` when many runs share one network:

```python
import river_route as rr

configs = rr.Configs.from_file('/path/to/config.yaml')

network = rr.Network.from_configs(configs)
print(network.stability_report(3600))        # which rivers are Muskingum-stable at this dt
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
