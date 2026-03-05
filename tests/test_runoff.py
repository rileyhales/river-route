"""Tests for river_route.runoff — grid coordinate extraction, voronoi diagrams, weight tables, qlateral aggregation."""
import os
import shutil
import tempfile

import numpy as np
import pytest
import xarray as xr

from conftest import DATA_DIR, ERA5_FILES, RFSv2ConfigsData
from river_route.runoff import (
    grid_weights,
    grid_to_qlateral,
)


# ═════════════════════════════════════════════════════════════════════════════
# grid_weights (end-to-end weight table generation)
# ═════════════════════════════════════════════════════════════════════════════

def test_grid_weights(vpu: RFSv2ConfigsData):
    """Generate a weight table and compare against existing known-good weight table."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    if not vpu.catchments.exists():
        pytest.skip(f'Missing {vpu.catchments}')

    result = grid_weights(
        ERA5_FILES[0], str(vpu.catchments),
        x_var='longitude', y_var='latitude', river_id_var='linkno',
    )

    expected_cols = {'linkno', 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion'}
    assert expected_cols.issubset(set(result.columns)), \
        f'Missing columns: {expected_cols - set(result.columns)}'
    # todo check for exact equality by sorting in a predictable way and then comparing the columns

    # Compare against known-good weight table
    known = xr.open_dataset(str(vpu.grid_weights_file))
    known_rivers = set(known['river_id'].values.tolist())
    result_rivers = set(result['linkno'].values.tolist())
    # Generated weights should cover the same rivers
    assert known_rivers == result_rivers, \
        f'{len(known_rivers - result_rivers)} rivers in known but not generated, ' \
        f'{len(result_rivers - known_rivers)} rivers generated but not in known'


# ═════════════════════════════════════════════════════════════════════════════
# grid_to_qlateral
# ═════════════════════════════════════════════════════════════════════════════

def test_grid_to_qlateral(vpu: RFSv2ConfigsData):
    """Aggregate ERA5 gridded runoff to qlateral depths and volumes using the weight table."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    ds = grid_to_qlateral(ERA5_FILES[0], grid_weights_file=str(vpu.grid_weights_file), var_runoff='ro',
                          var_x='longitude', var_y='latitude', var_t='valid_time')
    assert 'depth' in ds
    assert 'volume' in ds
    assert 'time' in ds.dims
    assert 'river_id' in ds.dims
    # Volumes should be non-negative (m³)
    assert np.all(ds['volume'].values >= 0), 'Negative volumes found'


def test_grid_to_qlateral_volumes_vs_depths_ratio(vpu: RFSv2ConfigsData):
    """Verify that volumes / depths ≈ catchment area for each river."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    ds = grid_to_qlateral(ERA5_FILES[0], grid_weights_file=str(vpu.grid_weights_file), var_runoff='ro',
                          var_x='longitude', var_y='latitude', var_t='valid_time')

    # Compute total catchment area per river from the weight table
    wt_df = xr.open_dataset(str(vpu.grid_weights_file))[['river_id', 'area_sqm']].to_dataframe()
    catchment_area = wt_df.groupby('river_id')['area_sqm'].sum()
    # Align to the river_id order in the output
    area = catchment_area.reindex(ds['river_id'].values).values

    # volumes should equal depths * area (for non-zero depths)
    depths = ds['depth'].values
    volumes = ds['volume'].values
    expected_volumes = depths * area[np.newaxis, :].astype(np.float32)

    # Only check where depths are non-trivial to avoid division noise
    mask = np.abs(depths) > 1e-12
    np.testing.assert_allclose(
        volumes[mask], expected_volumes[mask],
        rtol=1e-6,
        err_msg='volumes != depths * catchment_area',
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
            ds_inc['volume'].values, ds_cum['volume'].values,
            rtol=0.02, atol=0.15,
            err_msg='Cumulative input does not produce same volumes as incremental',
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
