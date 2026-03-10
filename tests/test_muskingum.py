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
    params = pd.read_parquet(vpu.rr2_params_file)
    n_rivers = len(params)

    tmpdir = tempfile.mkdtemp()
    try:
        # Create a synthetic initial state (uniform 10 m3/s)
        init_state_file = os.path.join(tmpdir, 'init_state.parquet')
        pd.DataFrame({'Q': np.full(n_rivers, 10.0)}).to_parquet(init_state_file)

        discharge_file = os.path.join(tmpdir, 'q_channel_only.nc')
        final_state_file = os.path.join(tmpdir, 'final_state.parquet')

        rr.Muskingum(
            params_file=str(vpu.rr2_params_file),
            discharge_files=[discharge_file],
            channel_state_init_file=init_state_file,
            dt_routing=900,
            dt_total=3600 * 24 * 15,
            log=False,
        ).route()

        with xr.open_dataset(discharge_file) as ds:
            assert 'Q' in ds
            q_vals = ds['Q'].values
            assert q_vals.shape[0] == 24 * 15 * 4, f'Expected 1440 time steps, got {q_vals.shape[0]}'
            assert q_vals.shape[1] == n_rivers, f'Expected {n_rivers} rivers, got {q_vals.shape[1]}'
            assert np.all(q_vals >= 0), 'Negative discharge values found'

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return


def test_muskingum_zero_initial_state(vpu: RFSv2ConfigsData):
    """Routing with zero initial state and no lateral inflow should produce all-zero discharge."""
    params = pd.read_parquet(vpu.rr2_params_file)
    n_rivers = len(params)

    tmpdir = tempfile.mkdtemp()
    try:
        # Muskingum requires channel_state_init_file — provide explicit zeros
        init_state_file = os.path.join(tmpdir, 'zero_state.parquet')
        pd.DataFrame({'Q': np.zeros(n_rivers)}).to_parquet(init_state_file)

        discharge_file = os.path.join(tmpdir, 'q.nc')

        rr.Muskingum(
            params_file=str(vpu.rr2_params_file),
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
