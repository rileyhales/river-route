"""
Route every run of the reference solution from its packaged inputs, and require the catchment runoff, the discharge,
and the final state to match the packaged outputs within TOLERANCE, with dates and river ids matching exactly. How the
package is built is described in reference_solution.py.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from reference_solution import assert_same, reference_runs, route_run, write_catchment_runoff


def test_catchment_runoff_is_reproduced(package: Path, manifest: dict, tmp_path: Path) -> None:
    catchment_runoff = manifest['catchment_runoff']
    for runoff_file, reference_file in zip(catchment_runoff['runoff_files'], catchment_runoff['files'], strict=True):
        regenerated_file = tmp_path / Path(reference_file).name
        write_catchment_runoff(package, catchment_runoff, runoff_file, regenerated_file)
        with xr.open_dataset(package / reference_file) as reference, xr.open_dataset(regenerated_file) as regenerated:
            for name in ('riverId', 'time'):
                np.testing.assert_array_equal(regenerated[name].to_numpy(), reference[name].to_numpy())
            for name in ('catchment_runoff', 'catchment_area'):
                assert_same(regenerated[name].to_numpy(), reference[name].to_numpy())


@pytest.mark.parametrize('name', reference_runs('', '', [], []))
def test_run_is_reproduced(name: str, package: Path, manifest: dict, tmp_path: Path) -> None:
    run = manifest['runs'][name]
    routed = {}  # each discharge file's name -> what the router handed its writer

    def keep_discharge(router, dates, discharge_array, discharge_file, runoff_file='', *, thread_pool=None, threads=1):
        routed[Path(discharge_file).name] = (dates.copy(), discharge_array.copy(), router.network.river_ids)

    discharge_files = [tmp_path / Path(reference_file).name for reference_file in run['discharge']]
    route_run(package, run, discharge_files, tmp_path / 'final_state.parquet', keep_discharge)
    for reference_file in run['discharge']:
        dates, discharge, river_ids = routed[Path(reference_file).name]
        with xr.open_zarr(package / reference_file) as reference:
            assert_same(discharge, reference['Q'].to_numpy())
            np.testing.assert_array_equal(dates.astype('datetime64[ns]'), reference['time'].to_numpy())
            np.testing.assert_array_equal(river_ids, reference['riverId'].to_numpy())
    final_state = pd.read_parquet(tmp_path / 'final_state.parquet')
    reference_state = pd.read_parquet(package / run['final_state'])
    np.testing.assert_array_equal(final_state['riverId'].to_numpy(), reference_state['riverId'].to_numpy())
    assert_same(final_state['Q'].to_numpy(), reference_state['Q'].to_numpy())
