"""
Fixtures every test file shares: the reference solution package and its manifest, and the Willamette River, a basin of
1,390 rivers cut from the package for tests that need to run in seconds. ``pytest --kernel numba`` or ``--kernel c``
routes every test with that kernel; by default the C kernel routes once ``rr compile`` has compiled it.
"""

import json
import os
from pathlib import Path

import pytest
import reference_solution
from reference_solution import GRID_NAMES, PACKAGE, WILLAMETTE, Basin

import river_route as rr
from river_route.network import streams


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        '--kernel',
        default='auto',
        choices=rr.kernels.KERNELS,
        help='the routing kernel the tests route with; the reference runs and test_kernels.py route with both',
    )


def pytest_configure(config: pytest.Config) -> None:
    reference_solution.KERNEL = config.getoption('kernel')


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
    network_file, weights_file = directory / 'network.parquet', directory / 'gridweights.nc'
    streams.subset_network_to_river(
        WILLAMETTE,
        package / configs['network_file'],
        network_file,
        package / configs['grid_weights_file'],
        weights_file,
    )
    forcing_files = [package / path for path in configs['forcing_files']]
    grid = rr.GridRunoff(weights_file, **GRID_NAMES)
    months = [(dates, volumes.copy()) for dates, volumes, _ in grid.catchment_reader(forcing_files)]
    return Basin(network_file, weights_file, forcing_files, months, rr.Network(network_file))
