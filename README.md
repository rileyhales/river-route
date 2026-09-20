# River Route

[![PyPI version](https://badge.fury.io/py/river-route.svg)](https://pypi.org/project/river-route/)

`river-route` is a Python package for routing runoff and discharge through large river
networks. It uses numba-compiled kernels and sparse matrix operations for efficient
Muskingum-family routing at watershed scale.

## Router Options

All routing runs through `Router`. The routing procedure is described by four config selector keys,
which together resolve to one compiled kernel:

| Key         | Options                                                                                                                           |
|-------------|-----------------------------------------------------------------------------------------------------------------------------------|
| `coeff`     | `static` (constant Muskingum K from `k`,`x`) or `dynamic` (nonlinear K = alpha*Q^beta from `alpha`,`beta`,`x`). Default `static`. |
| `forcing`   | `channel` (channel routing only, the default) or `vlateral` (lateral runoff inflow). One value only.                              |
| `transform` | `uniform` (the default) or `unit_hydrograph`. Only read when `forcing` is `vlateral`.                                             |
| `network`   | `standard` (one reach per river, the default).                                                                                    |

An unimplemented combination raises `NotImplementedError` listing the ones that exist.

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

configs = rr.Configs.from_file("./path/to/configs.yaml")
rr.Router(configs).route()
```

Configuration is held by a frozen `Configs` object, built from:

1. A YAML/JSON config file with `Configs.from_file`.
2. Keyword arguments, `Configs(params_file=..., ...)`.

`Configs.to_yaml` and `Configs.to_json` write the options to a file that `Configs.from_file` reads back. A
`Configs` is set once when it is built and is frozen afterward, with no copy-with-changes: build the one you
want. `Router` takes a `Configs`, and optionally a `Network` and a `RunoffGaussianGrid`. `Network` and `RunoffGaussianGrid` take the
options they need directly and each has a `from_configs` classmethod that reads those same values off a
`Configs`. `examples/config.yaml` groups every option under the class that reads it.

## Classes

| Class                | Owns                                                                                                  |
|----------------------|-------------------------------------------------------------------------------------------------------|
| `Configs`            | Every option, frozen. Built once and passed to the classes below.                                     |
| `Network`            | The river network: ids, topology, k and x, the concurrent partition, stability analysis, subdivision. |
| `Router`             | One simulation over a `Network`: coefficients, time options, channel state, the routing loop.         |
| `RunoffGaussianGrid` | Gridded runoff to lateral inflow, reusing one weight table across many runoff files.                  |

A `Network` parses and partitions a parameter table once and every simulation over it reuses the result, so it can
be built directly and handed to a `Router` when many runs share one network:

```python
network = rr.Network.from_configs(configs)
print(network.stability_report(3600))  # which rivers are Muskingum-stable at this dt
stabilized = network.stabilize(3600)  # in-memory network with reaches added so every one routes stably

router = rr.Router(configs)
router.network = network
router.route()
```

Core required inputs are:

- `params_file` (network topology and Muskingum parameters)
- One runoff source (`vlateral_files` or `grid_runoff_files` + `grid_weights_file`) when using `vlateral` forcing
- `discharge_dir` (or explicit `discharge_files`)

## CLI

```bash
rr --help
rr route /path/to/config.yaml
```

Copy `examples/config.yaml` to start; it lists every config key.

Subset a parameter table, and its grid weight table, to one river and everything upstream of it. The target
river becomes the outlet of the subset.

```bash
rr subset 12345 params.parquet params_subset.parquet --weights weights.nc --out-weights weights_subset.nc
```

## Testing

`pytest` is not a required dependency. You need to install `pytest` separately to run tests.

```bash
./tests/download_test_data.sh
pytest tests -v -s
```
