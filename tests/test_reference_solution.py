"""
Route every run of the reference solution from its packaged inputs, and require the catchment runoff, the discharge,
and the final state to match the packaged outputs within TOLERANCE, with dates and river ids matching exactly. How the
package is built is described in reference_solution.py.

The benchmark tests are opt-in (pytest -m benchmark). They time the single-threaded runs and fail when one is more than
BENCHMARK_SLACK slower than the baseline stored in the package, which RIVER_ROUTE_RECORD_BENCHMARK=1 records.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from reference_solution import RUN_NAMES, assert_same, route_run, write_catchment_runoff

from river_route.router import writers

BENCHMARK_SLACK = 1.15  # run to run noise on one machine is about 5 to 10 percent
BENCHMARK_REPEATS = 2  # after one untimed run that compiles, averaged


def test_catchment_runoff_is_reproduced(package: Path, manifest: dict, tmp_path: Path) -> None:
    catchment_runoff = manifest['catchment_runoff']
    for runoff_file, reference_file in zip(catchment_runoff['runoff_files'], catchment_runoff['files'], strict=True):
        regenerated_file = tmp_path / Path(reference_file).name
        write_catchment_runoff(package, catchment_runoff, runoff_file, regenerated_file)
        with xr.open_dataset(package / reference_file) as reference, xr.open_dataset(regenerated_file) as regenerated:
            for name in ('river_id', 'time'):
                np.testing.assert_array_equal(regenerated[name].to_numpy(), reference[name].to_numpy())
            for name in ('catchment_runoff', 'catchment_area'):
                assert_same(regenerated[name].to_numpy(), reference[name].to_numpy())


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
            assert_same(discharge, reference['Q'].to_numpy())
            np.testing.assert_array_equal(dates.astype('datetime64[ns]'), reference['time'].to_numpy())
            np.testing.assert_array_equal(river_ids, reference['river_id'].to_numpy())
    final_state = pd.read_parquet(tmp_path / 'final_state.parquet')
    reference_state = pd.read_parquet(package / run['final_state'])
    np.testing.assert_array_equal(final_state['river_id'].to_numpy(), reference_state['river_id'].to_numpy())
    assert_same(final_state['Q'].to_numpy(), reference_state['Q'].to_numpy())


@pytest.mark.benchmark
@pytest.mark.parametrize('name', ['static_standard', 'static_stabilized'])
def test_routing_is_not_slower(name: str, package: Path, manifest: dict, tmp_path: Path) -> None:
    run = manifest['runs'][name]
    discharge_files = [tmp_path / Path(path).name for path in run['discharge']]
    seconds = []
    for repeat in range(BENCHMARK_REPEATS + 1):
        start = time.perf_counter()
        route_run(package, run, discharge_files, tmp_path / 'final_state.parquet', writers.null_writer)
        if repeat:
            seconds.append(time.perf_counter() - start)
    mean = float(np.mean(seconds))
    baseline_file = package / 'benchmark_baseline.json'
    baseline = json.loads(baseline_file.read_text()) if baseline_file.exists() else {}
    if os.environ.get('RIVER_ROUTE_RECORD_BENCHMARK'):
        baseline_file.write_text(json.dumps(baseline | {name: mean}, indent=2))
        return
    if name not in baseline:
        pytest.skip(f'no baseline for {name}; record one with RIVER_ROUTE_RECORD_BENCHMARK=1 pytest -m benchmark')
    assert mean <= baseline[name] * BENCHMARK_SLACK, (
        f'{name} took {mean:.3f} s against a baseline of {baseline[name]:.3f} s'
    )
