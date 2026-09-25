"""
Route every run of the reference solution from its packaged inputs, and require the catchment runoff, the discharge,
and the final state to match the packaged outputs exactly. How the package is built is described in
reference_solution.py. Set RIVER_ROUTE_REFERENCE_SOLUTION to use a package stored somewhere other than tests/data.
"""

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from reference_solution import PACKAGE, RUN_NAMES, route_run, write_catchment_runoff


@pytest.fixture(scope='session')
def package() -> Path:
    path = Path(os.environ.get('RIVER_ROUTE_REFERENCE_SOLUTION', PACKAGE))
    if not (path / 'manifest.json').exists():
        pytest.fail(f'no reference solution at {path}; build it with tests/reference_solution.py')
    return path


@pytest.fixture(scope='session')
def manifest(package: Path) -> dict:
    return json.loads((package / 'manifest.json').read_text())


def test_catchment_runoff_is_reproduced(package: Path, manifest: dict, tmp_path: Path) -> None:
    catchment_runoff = manifest['catchment_runoff']
    for runoff_file, reference_file in zip(catchment_runoff['runoff_files'], catchment_runoff['files'], strict=True):
        regenerated_file = tmp_path / Path(reference_file).name
        write_catchment_runoff(package, catchment_runoff, runoff_file, regenerated_file)
        with xr.open_dataset(package / reference_file) as reference, xr.open_dataset(regenerated_file) as regenerated:
            xr.testing.assert_equal(regenerated, reference)  # every value and coordinate; attributes hold build dates


@pytest.mark.parametrize('name', RUN_NAMES)
def test_run_is_reproduced(name: str, package: Path, manifest: dict, tmp_path: Path) -> None:
    run = manifest['runs'][name]
    routed = {}  # each discharge file's name -> what the router handed its writer

    def keep_discharge(router, dates, discharge_array, discharge_file, runoff_file=''):
        routed[Path(discharge_file).name] = (dates.copy(), discharge_array.copy(), router.network.river_ids)

    discharge_files = [tmp_path / Path(reference_file).name for reference_file in run['discharge']]
    route_run(package, run, discharge_files, tmp_path / 'final_state.parquet', keep_discharge)
    for reference_file in run['discharge']:
        dates, discharge, river_ids = routed[Path(reference_file).name]
        with xr.open_zarr(package / reference_file) as reference:
            np.testing.assert_array_equal(discharge, reference['Q'].to_numpy())
            np.testing.assert_array_equal(dates.astype('datetime64[ns]'), reference['time'].to_numpy())
            np.testing.assert_array_equal(river_ids, reference['river_id'].to_numpy())
    pd.testing.assert_frame_equal(
        pd.read_parquet(tmp_path / 'final_state.parquet'),
        pd.read_parquet(package / run['final_state']),
        check_exact=True,
    )
