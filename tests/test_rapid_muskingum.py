"""Tests for RapidMuskingum routing — ERA5 grid-to-discharge integration tests."""
import os
import shutil
import tempfile

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import river_route as rr
from conftest import ERA5_FILES, ERA5_KWARGS, RFSv2ConfigsData


def test_rapid_muskingum_from_depths(vpu: RFSv2ConfigsData):
    """Route 1 month of ERA5 through grid weights and compare output against known discharge."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    tmpdir = tempfile.mkdtemp()
    try:
        discharge_file = os.path.join(tmpdir, 'q.nc')

        rr.RapidMuskingum(
            params_file=str(vpu.rr2_params_file),
            grid_weights_file=str(vpu.grid_weights_file),
            grid_runoff_files=[ERA5_FILES[0]],
            discharge_files=[discharge_file],
            log=False,
            progress_bar=False,
            **ERA5_KWARGS,
        ).route()

        assert os.path.exists(discharge_file)
        with xr.open_dataset(discharge_file) as ds_new, xr.open_dataset(vpu.discharge_files[0]) as ds_known:
            assert 'Q' in ds_new
            assert ds_new['Q'].shape == ds_known['Q'].shape, \
                f'Shape mismatch: {ds_new["Q"].shape} vs {ds_known["Q"].shape}'
            np.testing.assert_allclose(
                ds_new['Q'].values, ds_known['Q'].values, atol=0.01,
                err_msg='Routed discharge does not match known-good output',
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_initial_state_used(vpu: RFSv2ConfigsData):
    """Route same month with vs without initial state; first timestep should differ."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    params = pd.read_parquet(vpu.rr2_params_file)
    n_rivers = len(params)

    tmpdir = tempfile.mkdtemp()
    try:
        init_state_file = os.path.join(tmpdir, 'init_state.parquet')
        pd.DataFrame({'Q': np.full(n_rivers, 50.0)}).to_parquet(init_state_file)

        q_with_state = os.path.join(tmpdir, 'q_with.nc')
        q_without_state = os.path.join(tmpdir, 'q_without.nc')

        rr.RapidMuskingum(
            params_file=str(vpu.rr2_params_file),
            grid_weights_file=str(vpu.grid_weights_file),
            grid_runoff_files=[ERA5_FILES[0]],
            discharge_files=[q_with_state],
            channel_state_init_file=init_state_file,
            log=False, progress_bar=False,
            **ERA5_KWARGS,
        ).route()

        rr.RapidMuskingum(
            params_file=str(vpu.rr2_params_file),
            grid_weights_file=str(vpu.grid_weights_file),
            grid_runoff_files=[ERA5_FILES[0]],
            discharge_files=[q_without_state],
            log=False, progress_bar=False,
            **ERA5_KWARGS,
        ).route()

        with xr.open_dataset(q_with_state) as ds1, xr.open_dataset(q_without_state) as ds2:
            q1_first = ds1['Q'].values[0, :]
            q2_first = ds2['Q'].values[0, :]
            diff_first = np.abs(q1_first - q2_first).sum()
            assert diff_first > 0, 'First timestep should differ with different initial states'

            q1_last = ds1['Q'].values[-1, :]
            q2_last = ds2['Q'].values[-1, :]
            diff_last = np.abs(q1_last - q2_last).sum()
            assert diff_last < diff_first, 'Initial state influence should decay over time'
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_final_state_roundtrip(vpu: RFSv2ConfigsData):
    """Route months 1+2 at once vs month1 -> save state -> month2; month2 outputs must match."""
    if len(ERA5_FILES) < 2:
        pytest.skip('Need at least 2 ERA5 files')

    tmpdir = tempfile.mkdtemp()
    try:
        q_all = [os.path.join(tmpdir, f'q_all_{i}.nc') for i in range(2)]
        rr.RapidMuskingum(
            params_file=str(vpu.rr2_params_file),
            grid_weights_file=str(vpu.grid_weights_file),
            grid_runoff_files=ERA5_FILES[:2],
            discharge_files=q_all,
            log=False, progress_bar=False,
            **ERA5_KWARGS,
        ).route()

        q_m1 = os.path.join(tmpdir, 'q_m1.nc')
        state_after_m1 = os.path.join(tmpdir, 'state_m1.parquet')
        rr.RapidMuskingum(
            params_file=str(vpu.rr2_params_file),
            grid_weights_file=str(vpu.grid_weights_file),
            grid_runoff_files=[ERA5_FILES[0]],
            discharge_files=[q_m1],
            channel_state_final_file=state_after_m1,
            log=False, progress_bar=False,
            **ERA5_KWARGS,
        ).route()

        q_m2 = os.path.join(tmpdir, 'q_m2.nc')
        rr.RapidMuskingum(
            params_file=str(vpu.rr2_params_file),
            grid_weights_file=str(vpu.grid_weights_file),
            grid_runoff_files=[ERA5_FILES[1]],
            discharge_files=[q_m2],
            channel_state_init_file=state_after_m1,
            log=False, progress_bar=False,
            **ERA5_KWARGS,
        ).route()

        with xr.open_dataset(q_all[1]) as ds_all, xr.open_dataset(q_m2) as ds_split:
            assert ds_all['Q'].shape == ds_split['Q'].shape
            np.testing.assert_allclose(
                ds_all['Q'].values, ds_split['Q'].values,
                rtol=1e-4, atol=0.01,
                err_msg='State roundtrip: sequential does not match all-at-once',
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_ensemble_routing(vpu: RFSv2ConfigsData):
    """Route 2 months in ensemble mode; verify outputs written and final state is reasonable."""
    if len(ERA5_FILES) < 2:
        pytest.skip('Need at least 2 ERA5 files')

    tmpdir = tempfile.mkdtemp()
    try:
        discharge_files = [os.path.join(tmpdir, f'q_ens_{i}.nc') for i in range(2)]
        final_state_file = os.path.join(tmpdir, 'final_state_ensemble.parquet')

        rr.RapidMuskingum(
            params_file=str(vpu.rr2_params_file),
            grid_weights_file=str(vpu.grid_weights_file),
            grid_runoff_files=ERA5_FILES[:2],
            discharge_files=discharge_files,
            channel_state_final_file=final_state_file,
            runoff_processing_mode='ensemble',
            log=False, progress_bar=False,
            **ERA5_KWARGS,
        ).route()

        for qf in discharge_files:
            assert os.path.exists(qf), f'Missing output: {qf}'

        # Small negatives are expected from raw Muskingum math (discharge is clipped, state is not)
        assert os.path.exists(final_state_file)
        final_state = pd.read_parquet(final_state_file)
        assert 'Q' in final_state.columns
        assert np.all(final_state['Q'].values > -5), 'Excessively negative values in ensemble mean state'
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_rapid_muskingum_from_qlateral(vpu: RFSv2ConfigsData):
    """Route from pre-computed qlateral and compare against known discharge."""
    tmpdir = tempfile.mkdtemp()
    try:
        discharge_file = os.path.join(tmpdir, 'q.nc')

        rr.RapidMuskingum(
            params_file=str(vpu.rr2_params_file),
            qlateral_files=[vpu.qlateral_files[0]],
            discharge_files=[discharge_file],
            log=False,
            progress_bar=False,
        ).route()

        assert os.path.exists(discharge_file)
        with xr.open_dataset(discharge_file) as ds_new, xr.open_dataset(vpu.discharge_files[0]) as ds_known:
            assert 'Q' in ds_new
            assert ds_new['Q'].shape == ds_known['Q'].shape, \
                f'Shape mismatch: {ds_new["Q"].shape} vs {ds_known["Q"].shape}'
            np.testing.assert_allclose(
                ds_new['Q'].values, ds_known['Q'].values,
                rtol=1e-4, atol=0.01,
                err_msg='Discharge from qlateral does not match known-good output',
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
