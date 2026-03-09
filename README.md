# River Route

[![Documentation Status](https://readthedocs.org/projects/river-route/badge/?version=latest)](https://river-route.hales.app/en/latest/)
[![PyPI version](https://badge.fury.io/py/river-route.svg)](https://pypi.org/project/river-route/)
[![GitHub repo size](https://img.shields.io/github/repo-size/rileyhales/river-route)](https://github.com/rileyhales/river-route)
![License](https://img.shields.io/github/license/rileyhales/river-route)

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

For local development:

```bash
git clone https://github.com/rileyhales/river-route.git
cd river-route
conda env create -f environment.yaml
conda activate rr
python -m pip install -e ".[all]"
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

## Documentation

- Hosted docs: https://river-route.hales.app
- Migration guide: [`docs/migrating/v1-to-v2.md`](docs/migrating/v1-to-v2.md)

Build docs locally:

```bash
mkdocs serve
```

## Testing

```bash
pytest tests -v
```

Integration tests use external sample data:

```bash
./tests/download_test_data.sh
pytest tests -v -s
```
