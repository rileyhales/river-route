# River Route

is there an argument to have the kernels in numba and c moved to river-route-kernels which is published separately?
[![PyPI version](https://badge.fury.io/py/river-route.svg)](https://pypi.org/project/river-route/)

`river-route` is a Python package for routing runoff and discharge through large river networks. The routing kernel is Muskingum with static
parameters with valid configuration options enforced everywhere possible.

The kernel is primarily developed in Python with Numba but a duplicate C kernel is also maintained. You can run ``rr compile`` to compile the C kernel
which may speed up performance. The compile CLI command expects LLVM on the machine to access clang and clang-tidy.

The C code is also used to publish a WebAssembly version of the kernel with the most minimal harness to NPM, also named `river-route`. To use the web
version, use the primarily bot generated companion repo [river-route-web](https://github.com/rileyhales/river-route-web)

## Router Options

All routing runs through `Router`, with or without a forcing. The `forcing` selector describes what is routed:

| Key       | Options                                                |
|-----------|--------------------------------------------------------|
| `forcing` | `channel` (default) `catchment`, `grid`, `ecmwf_grib`. |

Compute options set the step each river computes in, and whether it is split or subcycled to keep its coefficients valid:

| Key              | Options                                                                                     |
|------------------|---------------------------------------------------------------------------------------------|
| `dt_routing`     | `auto` (default), each river's own valid step, or a positive integer for every river        |
| `use_subreaches` | `true` (default): split rivers too long for their step into sub-reaches; `auto` requires it |
| `use_subcycles`  | `true` (default): route rivers too short for their step in subcycles; `auto` requires it    |

Some rivers cannot

Every forcing routes, and `dt_routing: auto` with `use_subreaches` makes every coefficient of every river non-negative.

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
`Configs` is set once when it is built and is frozen afterward, with no copy-with-changes: build the one you want. `Router` takes only a `Configs`,
and builds its `Network` and the `RunoffForcing` class its `forcing` names from it.
`Network` and `GridRunoff` take the options they need directly and each has a `from_configs` classmethod that reads those same values off a `Configs`.
`examples/config.json` lists every option.

Core required inputs are:

- `network_file` (network topology and Muskingum parameters)
- `forcing_files` when `forcing` is not `channel`, plus `grid_weights_file` for the grid forcings
- `discharge_dir` (or explicit `discharge_files`)
