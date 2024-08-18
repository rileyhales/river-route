# River Route

[![Documentation Status](https://readthedocs.org/projects/river-route/badge/?version=latest)](https://river-route.hales.app/en/latest/)
[![PyPI version](https://badge.fury.io/py/river-route.svg)](https://pypi.org/project/river-route/)
[![GitHub repo size](https://img.shields.io/github/repo-size/rileyhales/river-route)](https://github.com/rileyhales/river-route)
![License](https://img.shields.io/github/license/rileyhales/river-route)

`river-route` is a Python package for routing catchment volumes through a river network. Routing calculations are vectorized and use numpy and scipy
which keeps the array computation times on par with faster compiled languages.

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
