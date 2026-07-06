# River Route

[![PyPI version](https://badge.fury.io/py/river-route.svg)](https://pypi.org/project/river-route/)

`river-route` is a Python package for routing runoff and discharge through large river
networks. It uses numba-compiled kernels and sparse matrix operations for efficient
Muskingum-family routing at watershed scale.

## Router Options

The public API is a single `Router` class. The routing procedure is described by three
config selector keys:

| Key       | Options                                                                                                                           |
|-----------|-----------------------------------------------------------------------------------------------------------------------------------|
| `coeff`   | `static` (constant Muskingum K from `k`,`x`) or `dynamic` (nonlinear K = alpha*Q^beta from `alpha`,`beta`,`x`). Default `static`. |
| `forcing` | `channel` (channel routing only, the default) or `vlateral` (lateral runoff inflow). One value only.                              |
| `network` | `standard` (one reach per river, the default).                                                                                    |

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

rr.Router("./path/to/configs.yaml").route()

# selectors can also be passed (or overridden) as keyword arguments
rr.Router("./path/to/configs.yaml", forcing="vlateral").route()
```

Configuration can be provided by:

1. A YAML/JSON config file path.
2. Keyword arguments.
3. Both (kwargs override file values).

Core required inputs are:

- `params_file` (network topology and Muskingum parameters)
- One runoff source (`vlateral_files` or `grid_runoff_files` + `grid_weights_file`) when using `vlateral` forcing
- `discharge_dir` (or explicit `discharge_files`)

## CLI

```bash
rr --help
rr route examples/config.yaml
```

## Testing

`pytest` is not a required dependency. You need to install `pytest` separately to run tests.

```bash
./tests/download_test_data.sh
pytest tests -v -s
```
