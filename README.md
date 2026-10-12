# River Route

is there an argument to have the kernels in numba and c moved to river-route-kernels which is published separately?
[![PyPI version](https://badge.fury.io/py/river-route.svg)](https://pypi.org/project/river-route/)

`river-route` is a Python package for routing runoff and discharge through large river networks. The routing kernel is Muskingum with static
parameters with valid configuration options enforced everywhere possible.

The kernel is primarily developed in Python with Numba but a duplicate C kernel is also maintained. You can run ``rr compile`` to compile the C kernel
which may speed up performance. The compile CLI command expects LLVM on the machine to build and pass the checks.

The C code is also used to publish a WebAssembly version of the kernel with the most minimal harness to NPM, also named `river-route`. To use the web
version, use the primarily bot generated companion repo river-route-web

## Router Options

All routing runs through `Router`, with or without a forcing. The `forcing` selector describes what is routed:

| Key       | Options                                                |
|-----------|--------------------------------------------------------|
| `forcing` | `channel` (default) `catchment`, `grid`, `ecmwf_grib`. |

Compute options set the step each river computes in, and whether it is conditioned to keep its coefficients valid:

| Key            | Options                                                                                          |
|----------------|--------------------------------------------------------------------------------------------------|
| `dt_routing`   | `auto` (default), each river's own valid step, or a positive integer for every river             |
| `conditioning` | `true` (default): sub-reaches and subcycles for rivers invalid at their step; `auto` requires it |

Every forcing routes, and `dt_routing: auto` makes every coefficient of every river non-negative.

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

### The C kernel

River-route routes with numba until its C kernel is compiled on your machine, which needs the clang of
[LLVM](https://llvm.org) (on macOS, the Xcode Command Line Tools' clang will do). Compile it once per environment:

```bash
rr compile              # finds clang on the PATH, or give its folder: rr compile --llvm-bin /opt/homebrew/opt/llvm/bin
```

`Router.route` then routes with the C kernel, its default `kernel='auto'`; `kernel='numba'` or `kernel='c'` chooses one. Both kernels are the same
routines and give the same answers.

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
`Configs` is set once when it is built and is frozen afterward, with no copy-with-changes: build the one you want. `Router` takes only a `Configs`,
and builds its `Network` and the `Runoff` class its `forcing` names from it.
`Network` and `GridRunoff` take the options they need directly and each has a `from_configs` classmethod that reads those same values off a `Configs`.
`examples/config.json` lists every option.

## Classes

| Class        | Owns                                                                                                            |
|--------------|-----------------------------------------------------------------------------------------------------------------|
| `Configs`    | Every option, frozen. Built once and passed to the classes below.                                               |
| `Network`    | The river network: ids, topology, k and x, the concurrent partition, coefficient signs, sub-reaches, subcycles. |
| `Router`     | One simulation over a `Network`: coefficients, time options, channel state, the routing loop.                   |
| `GridRunoff` | Gridded runoff to catchment runoff, reusing one weight table across many forcing files.                         |

A `Router` parses and partitions its network table once, and every `route()` call on it reuses the result.

Core required inputs are:

- `network_file` (network topology and Muskingum parameters)
- `forcing_files` when `forcing` is not `channel`, plus `grid_weights_file` for the grid forcings
- `discharge_dir` (or explicit `discharge_files`)

## CLI

```bash
rr --help
rr route /path/to/config.json
```

Copy `examples/config.json` to start; it lists every config key.

Subset a network table, and its grid weight table, to one river and everything upstream of it. The target river becomes the outlet of the subset.

```bash
rr subset 12345 --network network.parquet --out-network network_subset.parquet --weights weights.nc --out-weights weights_subset.nc
```
