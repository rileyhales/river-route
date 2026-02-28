# River Route

[![Documentation Status](https://readthedocs.org/projects/river-route/badge/?version=latest)](https://river-route.hales.app/en/latest/)
[![PyPI version](https://badge.fury.io/py/river-route.svg)](https://pypi.org/project/river-route/)
[![GitHub repo size](https://img.shields.io/github/repo-size/rileyhales/river-route)](https://github.com/rileyhales/river-route)
![License](https://img.shields.io/github/license/rileyhales/river-route)

`river-route` is a Python package for routing catchment runoff volumes through a river network. It provides three
routers: `Muskingum`, `RapidMuskingum`, and `UnitMuskingum`. Routing calculations are vectorized using numpy
and scipy sparse matrices, keeping computation times competitive with compiled languages.

Please visit https://river-route.hales.app for documentation

```bash
# Install from PyPI
pip install river-route
```

```bash
# Install from source
git clone https://github.com/rileyhales/river-route.git
cd river-route
# Install dependencies using conda
conda env create -f environment.yml
conda activate rr
# Install package in editable mode
python -m pip install -e .
```

```python
import river_route as rr

# Muskingum: channel routing only without lateral inflows. Requires initial discharge state and explicit time parameters.
(
    rr
    .Muskingum('examples/config_muskingum.yaml')
    .route()
)
# RapidMuskingum: places runoff volumes at the catchment inlets ignoring overland flow travel times
(
    rr
    .RapidMuskingum('examples/config.yaml')
    .route()
)
# UnitMuskingum: applies a pluggable unit-hydrograph runoff transform before Muskingum channel routing
(
    rr
    .UnitMuskingum('examples/config_unit_muskingum.yaml')
    .route()
)
```
