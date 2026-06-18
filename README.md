# River Route

[![PyPI version](https://badge.fury.io/py/river-route.svg)](https://pypi.org/project/river-route/)

`river-route` is a Python package for routing runoff and discharge through large river
networks. It uses numba-compiled kernels and sparse matrix operations for efficient
Muskingum-family routing at watershed scale.

## Router Options

| Router           | Use case                                                        |
|------------------|-----------------------------------------------------------------|
| `Muskingum`      | Channel routing only (no lateral runoff input).                 |
| `RapidMuskingum` | Route runoff directly to channels at each timestep.             |
| `UnitMuskingum`  | Transform runoff with a unit hydrograph before channel routing. |

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

(
    rr
    .RapidMuskingum("examples/config_rapid_muskingum.yaml")
    .route()
)
```

Configuration can be provided by:

1. A YAML/JSON config file path.
2. Keyword arguments.
3. Both (kwargs override file values).

Core required inputs are:

- `params_file` (network topology and Muskingum parameters)
- One runoff source (`qlateral_files` or `grid_runoff_files` + `grid_weights_file`) for transform routers
- `discharge_dir` (or explicit `discharge_files`)

## CLI

```bash
rr --help
rr RapidMuskingum examples/config_rapid_muskingum.yaml
rr UnitMuskingum examples/config_unit_muskingum.yaml
```

## Testing

`pytest` is not a required dependency. You need to install `pytest` separately to run tests.

```bash
./tests/download_test_data.sh
pytest tests -v -s
```
