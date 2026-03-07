import os
import shutil
import tempfile

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import xarray as xr

from conftest import DATA_DIR, ERA5_FILES, RFSv2ConfigsData
from river_route.runoff import grid_weights, grid_to_qlateral


def test_grid_weights(vpu: RFSv2ConfigsData):
    """Generate a weight table and compare against existing known-good weight table."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    if not vpu.catchments.exists():
        pytest.skip(f'Missing {vpu.catchments}')

    # Detect the river_id column name in the catchments file (LINKNO vs linkno)
    catchment_cols = gpd.read_parquet(str(vpu.catchments)).columns
    river_id_col = next((c for c in catchment_cols if c.lower() == 'linkno'), None)
    if river_id_col is None:
        pytest.fail(f'No linkno/LINKNO column in {vpu.catchments}; found {list(catchment_cols)}')

    result = grid_weights(ERA5_FILES[0], str(vpu.catchments), var_x='longitude', var_y='latitude',
                          var_river_id=river_id_col)

    expected_cols = {river_id_col, 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion'}
    assert expected_cols.issubset(set(result.columns)), f'Missing columns: {expected_cols - set(result.columns)}'

    # Compare against known-good weight table by sorting both on (river_id, x_index, y_index)
    sort_keys = ['river_id', 'x_index', 'y_index', 'area_sqm']
    compare_cols = ['river_id', 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion']
    known = (
        xr.open_dataset(str(vpu.grid_weights_file))
        [compare_cols]
        .to_dataframe()
        .reset_index(drop=True)
        .sort_values(sort_keys)
        .reset_index(drop=True)
    )
    generated = (
        result
        .rename(columns={river_id_col: 'river_id'})
        [compare_cols]
        .sort_values(sort_keys)
        .reset_index(drop=True)
    )
    pd.testing.assert_frame_equal(generated, known, check_exact=False, check_dtype=False, rtol=1e-5)


def test_grid_to_qlateral(vpu: RFSv2ConfigsData):
    """Aggregate ERA5 gridded runoff to qlateral depths using the weight table."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    ds = grid_to_qlateral(ERA5_FILES[0], grid_weights_file=str(vpu.grid_weights_file), var_runoff='ro',
                          var_x='longitude', var_y='latitude', var_t='valid_time', as_volumes=True)
    assert 'qlateral' in ds
    assert 'time' in ds.dims
    assert 'river_id' in ds.dims

    # compare for exact match against known-good qlateral file for this month
    known_file = vpu.qlateral_files[0]
    ds_known = xr.open_dataset(known_file)
    np.testing.assert_allclose(
        ds['qlateral'].values, ds_known['qlateral'].values, atol=0.1,
        err_msg='Aggregated qlateral does not match known-good output',
    )


def test_grid_to_qlateral_cumulative_input(vpu: RFSv2ConfigsData):
    """Cumulative runoff input should produce the same volumes as incremental input."""
    incremental_file = DATA_DIR / 'era5' / 'era5_194001.nc'
    if not incremental_file.exists():
        pytest.skip('Missing ERA5 file')

    tmpdir = tempfile.mkdtemp()
    try:
        # Generate cumulative version on the fly
        cumulative_file = os.path.join(tmpdir, 'era5_194001_cumulative.nc')
        with xr.open_dataset(incremental_file) as ds:
            ro_cum = np.cumsum(ds['ro'].values.astype(np.float64), axis=0).astype(np.float32)
            ds_out = ds.copy()
            ds_out['ro'] = (ds['ro'].dims, ro_cum)
            ds_out.to_netcdf(cumulative_file)

        kwargs = dict(
            grid_weights_file=str(vpu.grid_weights_file),
            var_runoff='ro', var_x='longitude', var_y='latitude', var_t='valid_time',
        )

        ds_inc = grid_to_qlateral(str(incremental_file), cumulative=False, **kwargs)
        ds_cum = grid_to_qlateral(cumulative_file, cumulative=True, **kwargs)

        # Tolerance is loose because cumsum -> float32 storage -> diff loses precision
        # relative to direct incremental aggregation
        np.testing.assert_allclose(
            ds_inc['qlateral'].values, ds_cum['qlateral'].values,
            rtol=0.02, atol=0.15,
            err_msg='Cumulative input does not produce same qlateral as incremental',
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
