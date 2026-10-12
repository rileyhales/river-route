# River-Route

`river-route` routes runoff and discharge through large river networks using numba-accelerated Muskingum-family
solvers that route each river's whole time series in turn, from upstream to downstream.

## Describe your routing

All routing runs through `rr.Router`, which routes Muskingum with static coefficients: each river's k and x, from the
`muskingumK` and `muskingumX` columns of the network file, are fixed for the whole run. You describe what is routed with
the `forcing` selector.

| Selector  | Options                                                                                                                                                                                                          | Default     |
|-----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------|
| `forcing` | `'channel'` (channel routing only), or the form of the `forcing_files` routed into the rivers: `'catchment'`, `'grid'`, or `'ecmwf_grib'`. Each step's runoff enters its river at a constant rate over the step. | `'channel'` |

Two compute options set how each river computes, in two jobs: `dt_routing` sets each river's step, and `conditioning`
then makes every coefficient non-negative at it, the condition for valid Muskingum routing.

| Compute option | Options                                                                                                                                                                                                       | Default  |
|----------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------|
| `dt_routing`   | `'auto'` gives each river its own valid step, the runoff step over the fewest subcycles that keep its coefficients non-negative; a step in seconds, a divisor of the runoff step, is the step of every river. | `'auto'` |
| `conditioning` | `True` splits each river too long for its step into sub-reaches and routes each river too short for it in subcycles of its step; `False` routes every river whole at its step. `'auto'` requires it.          | `True`   |

Every forcing routes, and `dt_routing: auto` makes every coefficient of every river non-negative.

```python
import river_route as rr

configs = rr.Configs.from_json("/path/to/config.json")  # sets the forcing and the other options
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
    forcing_files=["/path/to/catchment_runoff.nc"],
    discharge_dir="/path/to/output/",
)
configs.to_json("/path/to/config.json")
rr.Router(configs).route()
```

## Classes

| Class     | Owns                                                                                                   |
|-----------|--------------------------------------------------------------------------------------------------------|
| `Configs` | Every option, frozen. Built once and given to the classes below.                                        |
| `Network` | The river network: ids, topology, `k` and `x`, the concurrent partition, coefficient signs, sub-reaches, subcycles. |
| `Router`  | One simulation over a `Network`: coefficients, time options, channel state, the routing loop.            |
| `GridRunoff` | Gridded runoff to catchment runoff, reusing one weight table across many forcing files.             |

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
