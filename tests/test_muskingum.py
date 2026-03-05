"""Tests for Muskingum channel-only routing (no lateral inflow)."""
import os
import shutil
import tempfile

import numpy as np
import pandas as pd
import xarray as xr

import river_route as rr
from conftest import RFSv2ConfigsData


def test_muskingum_channel_only(vpu: RFSv2ConfigsData):
    """Route from a synthetic initial state with no lateral inflow; verify decay and final state."""
    params = pd.read_parquet(vpu.rr1_params_file)
    n_rivers = len(params)

    tmpdir = tempfile.mkdtemp()
    try:
        # Create a synthetic initial state (uniform 10 m3/s)
        init_state_file = os.path.join(tmpdir, 'init_state.parquet')
        pd.DataFrame({'Q': np.full(n_rivers, 10.0)}).to_parquet(init_state_file)

        discharge_file = os.path.join(tmpdir, 'q_channel_only.nc')
        final_state_file = os.path.join(tmpdir, 'final_state.parquet')

        rr.Muskingum(
            params_file=str(vpu.rr1_params_file),
            discharge_files=[discharge_file],
            channel_state_init_file=init_state_file,
            channel_state_final_file=final_state_file,
            dt_routing=3600,
            dt_total=86400,
            log=False,
        ).route()

        with xr.open_dataset(discharge_file) as ds:
            assert 'Q' in ds
            q_vals = ds['Q'].values
            assert q_vals.shape[0] == 24, f'Expected 24 time steps, got {q_vals.shape[0]}'
            assert np.all(q_vals >= 0), 'Negative discharge values found'
            # With no lateral inflow, total discharge should decrease over time
            total_per_step = q_vals.sum(axis=1)
            assert total_per_step[-1] <= total_per_step[0], \
                'Total discharge should not increase without lateral inflow'

        assert os.path.exists(final_state_file)
        final_state = pd.read_parquet(final_state_file)
        assert 'Q' in final_state.columns
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return


def test_muskingum_zero_initial_state(vpu: RFSv2ConfigsData):
    """Routing with zero initial state and no lateral inflow should produce all-zero discharge."""
    params = pd.read_parquet(vpu.rr1_params_file)
    n_rivers = len(params)

    tmpdir = tempfile.mkdtemp()
    try:
        # Muskingum requires channel_state_init_file — provide explicit zeros
        init_state_file = os.path.join(tmpdir, 'zero_state.parquet')
        pd.DataFrame({'Q': np.zeros(n_rivers)}).to_parquet(init_state_file)

        discharge_file = os.path.join(tmpdir, 'q.nc')

        rr.Muskingum(
            params_file=str(vpu.rr1_params_file),
            discharge_files=[discharge_file],
            channel_state_init_file=init_state_file,
            dt_routing=3600,
            dt_total=3600 * 6,
            log=False,
        ).route()

        with xr.open_dataset(discharge_file) as ds:
            assert np.all(ds['Q'].values == 0)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
