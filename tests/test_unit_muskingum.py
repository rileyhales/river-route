"""Tests for UnitMuskingum routing — sparse UH convolution + Muskingum."""
import os
import shutil
import tempfile

import numpy as np
import pandas as pd
import scipy.sparse
import xarray as xr

import river_route as rr
from conftest import RFSv2ConfigsData


def _make_sparse_kernel(n_rivers, n_kernel_steps, tmpdir):
    """Build a trivial instant-response sparse kernel and save as npz."""
    kernel = np.zeros((n_kernel_steps, n_rivers), dtype=np.float64)
    kernel[0, :] = 1.0
    path = os.path.join(tmpdir, 'kernel.npz')
    scipy.sparse.save_npz(path, scipy.sparse.csr_matrix(kernel))
    return path


def test_unit_muskingum_synthetic(vpu: RFSv2ConfigsData):
    """Build synthetic runoff + trivial sparse kernel; verify UnitMuskingum runs and produces valid output."""
    params = pd.read_parquet(vpu.rr2_params_file)
    river_ids = params['river_id'].values
    n_rivers = len(river_ids)

    tmpdir = tempfile.mkdtemp()
    try:
        n_timesteps = 24
        n_kernel_steps = 3
        dates = pd.date_range('2020-01-01', periods=n_timesteps, freq='h')
        np.random.seed(42)
        depths = np.random.uniform(0, 0.001, (n_timesteps, n_rivers)).astype(np.float32)
        runoff_ds = xr.Dataset(
            {'runoff': xr.DataArray(depths, dims=('time', 'river_id'))},
            coords={'time': dates, 'river_id': river_ids},
        )
        runoff_file = os.path.join(tmpdir, 'synthetic_depths.nc')
        runoff_ds.to_netcdf(runoff_file)

        kernel_file = _make_sparse_kernel(n_rivers, n_kernel_steps, tmpdir)
        discharge_file = os.path.join(tmpdir, 'q_uh.nc')
        final_state_file = os.path.join(tmpdir, 'final_state.parquet')
        transformer_state_file = os.path.join(tmpdir, 'transformer_state.parquet')

        rr.UnitMuskingum(
            params_file=str(vpu.rr2_params_file),
            transformer_kernel_file=kernel_file,
            qlateral_files=[runoff_file],
            discharge_files=[discharge_file],
            channel_state_final_file=final_state_file,
            transformer_state_final_file=transformer_state_file,
            log=False,
            progress_bar=False,
        ).route()

        assert os.path.exists(discharge_file)
        with xr.open_dataset(discharge_file) as ds:
            assert 'Q' in ds
            assert np.all(ds['Q'].values >= 0), 'Negative discharge found'

        assert os.path.exists(final_state_file)
        assert os.path.exists(transformer_state_file)

        # State saved transposed: (n_rivers, n_kernel_steps)
        state_df = pd.read_parquet(transformer_state_file)
        assert state_df.shape == (n_rivers, n_kernel_steps), \
            f'Transformer state shape {state_df.shape} != expected ({n_rivers}, {n_kernel_steps})'
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_unit_muskingum_transformer_state_roundtrip(vpu: RFSv2ConfigsData):
    """Route 2 files at once vs file1 -> save state -> file2; file2 outputs must match."""
    params = pd.read_parquet(vpu.rr1_params_file)
    river_ids = params['river_id'].values
    n_rivers = len(river_ids)

    tmpdir = tempfile.mkdtemp()
    try:
        n_timesteps = 24
        n_kernel_steps = 3
        np.random.seed(42)

        runoff_files = []
        for i in range(2):
            dates = pd.date_range(f'2020-01-0{i + 1}', periods=n_timesteps, freq='h')
            depths = np.random.uniform(0, 0.001, (n_timesteps, n_rivers)).astype(np.float32)
            ds = xr.Dataset(
                {'runoff': xr.DataArray(depths, dims=('time', 'river_id'))},
                coords={'time': dates, 'river_id': river_ids},
            )
            path = os.path.join(tmpdir, f'runoff_{i}.nc')
            ds.to_netcdf(path)
            runoff_files.append(path)

        kernel_file = _make_sparse_kernel(n_rivers, n_kernel_steps, tmpdir)

        # Route both files at once
        q_all = [os.path.join(tmpdir, f'q_all_{i}.nc') for i in range(2)]
        rr.UnitMuskingum(
            params_file=str(vpu.rr2_params_file),
            transformer_kernel_file=kernel_file,
            qlateral_files=runoff_files,
            discharge_files=q_all,
            log=False, progress_bar=False,
        ).route()

        # Route file 1 and save both channel + transformer state
        q_f1 = os.path.join(tmpdir, 'q_f1.nc')
        channel_state = os.path.join(tmpdir, 'channel_state.parquet')
        transformer_state = os.path.join(tmpdir, 'transformer_state.parquet')
        rr.UnitMuskingum(
            params_file=str(vpu.rr2_params_file),
            transformer_kernel_file=kernel_file,
            qlateral_files=[runoff_files[0]],
            discharge_files=[q_f1],
            channel_state_final_file=channel_state,
            transformer_state_final_file=transformer_state,
            log=False, progress_bar=False,
        ).route()

        # Route file 2 from saved state
        q_f2 = os.path.join(tmpdir, 'q_f2.nc')
        rr.UnitMuskingum(
            params_file=str(vpu.rr2_params_file),
            transformer_kernel_file=kernel_file,
            qlateral_files=[runoff_files[1]],
            discharge_files=[q_f2],
            channel_state_init_file=channel_state,
            transformer_state_init_file=transformer_state,
            log=False, progress_bar=False,
        ).route()

        # File 2 output from split routing should match file 2 from all-at-once
        with xr.open_dataset(q_all[1]) as ds_all, xr.open_dataset(q_f2) as ds_split:
            assert ds_all['Q'].shape == ds_split['Q'].shape
            np.testing.assert_allclose(
                ds_all['Q'].values, ds_split['Q'].values,
                rtol=1e-4, atol=0.01,
                err_msg='UnitMuskingum state roundtrip: sequential does not match all-at-once',
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
