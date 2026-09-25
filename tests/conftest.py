"""
Fixtures every test file shares: the reference solution package and its manifest, and the Willamette River, a basin of
1,390 rivers cut from the package for tests that need to run in seconds.
"""

import json
import os
from pathlib import Path

import pytest
from reference_solution import GRID_NAMES, PACKAGE, WILLAMETTE, Basin

import river_route as rr
from river_route.network import streams


@pytest.fixture(scope='session')
def package() -> Path:
    path = Path(os.environ.get('RIVER_ROUTE_REFERENCE_SOLUTION', PACKAGE))
    if not (path / 'manifest.json').exists():
        pytest.fail(f'no reference solution at {path}; build it with tests/reference_solution.py')
    return path


@pytest.fixture(scope='session')
def manifest(package: Path) -> dict:
    return json.loads((package / 'manifest.json').read_text())


@pytest.fixture(scope='session')
def willamette(package: Path, manifest: dict, tmp_path_factory: pytest.TempPathFactory) -> Basin:
    configs = manifest['runs']['static_standard']['configs']
    directory = tmp_path_factory.mktemp('willamette')
    params_file, weights_file = directory / 'routing.parquet', directory / 'gridweights.nc'
    streams.subset_configs_to_river(
        WILLAMETTE, package / configs['params_file'], params_file, package / configs['grid_weights_file'], weights_file
    )
    runoff_files = [package / path for path in configs['runoff_files']]
    grid = rr.GaussianGridRunoff(weights_file, **GRID_NAMES)
    months = [(dates, volumes.copy()) for dates, volumes, _ in grid.catchment_reader(runoff_files)]
    return Basin(params_file, weights_file, runoff_files, months, rr.Network(params_file))
