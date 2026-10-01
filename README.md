# River Route

[![PyPI version](https://badge.fury.io/py/river-route.svg)](https://pypi.org/project/river-route/)

`river-route` is a Python package for routing runoff and discharge through large river
networks. Its numba-compiled kernels route each river's whole time series in turn, from upstream to downstream, for
efficient Muskingum-family routing at watershed scale.

## Router Options

All routing runs through `Router`. The routing procedure is described by config selector keys,
which together choose the kernel:

| Key            | Options                                                                                                                                                                                  |
|----------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `coefficients` | `static` (constant Muskingum K from `muskingumK`,`muskingumX`) or `dynamic` (nonlinear K = dynamicAlpha*Q^dynamicBeta from `dynamicAlpha`,`dynamicBeta`,`muskingumX`). Default `static`. |
| `forcing`      | `channel` (channel routing only, the default), or the form of the `runoff_files` routed into the rivers: `catchment`, `grid`, or `ecmwf_grib`.                                           |
| `transform`    | `uniform` (the default and only option): each step's catchment runoff enters its river at a constant rate over the step.                                                                 |
| `network_type` | `standard` (one reach per river, the default) or `stabilized` (each river routed in the substeps or subcycles that make it stable at `dt_routing`, where any do).                        |

The one combination no routing method routes yet, `dynamic` coefficients on a `stabilized` network, raises
`NotImplementedError` before any runoff is read.

## Installation

```bash
pip install river-route
```

Or from source with [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/rileyhales/river-route.git
cd river-route
uv sync                 # create the environment and install river-route
uv sync --group dev     # ...or include the test and docs tooling
```

## Quick Start

```python
import river_route as rr

configs = rr.Configs.from_json("./path/to/configs.json")
rr.Router(configs).route()
```

Configuration is held by a frozen `Configs` object, built from:

1. A JSON config file with `Configs.from_json`.
2. Keyword arguments, `Configs(network_file=..., ...)`.

`Configs.to_json` writes the options to a file that `Configs.from_json` reads back. A
`Configs` is set once when it is built and is frozen afterward, with no copy-with-changes: build the one you
want. `Router` takes only a `Configs`, and builds its `Network` and the `Runoff` class its `forcing` names from it.
`Network` and `GridRunoff` take the options they need directly and each has a `from_configs` classmethod that reads
those same values off a `Configs`. `examples/config.json` lists every option.

## Classes

| Class        | Owns                                                                                                             |
|--------------|------------------------------------------------------------------------------------------------------------------|
| `Configs`    | Every option, frozen. Built once and passed to the classes below.                                                |
| `Network`    | The river network: ids, topology, k and x, the concurrent partition, stability analysis, substeps and subcycles. |
| `Router`     | One simulation over a `Network`: coefficients, time options, channel state, the routing loop.                    |
| `GridRunoff` | Gridded runoff to catchment runoff, reusing one weight table across many runoff files.                           |

A `Router` parses and partitions its network table once, and every `route()` call on it reuses the result.

Core required inputs are:

- `network_file` (network topology and Muskingum parameters)
- `runoff_files` when `forcing` is not `channel`, plus `grid_weights_file` for the grid forcings
- `discharge_dir` (or explicit `discharge_files`)

## CLI

```bash
rr --help
rr route /path/to/config.json
```

Copy `examples/config.json` to start; it lists every config key.

Subset a network table, and its grid weight table, to one river and everything upstream of it. The target
river becomes the outlet of the subset.

```bash
rr subset 12345 --network network.parquet --out-network network_subset.parquet --weights weights.nc --out-weights weights_subset.nc
```
