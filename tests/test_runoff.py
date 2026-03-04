"""Tests for river_route.runoff — grid coordinate extraction, voronoi diagrams, weight tables, catchment aggregation."""
import os
import shutil
import tempfile

import numpy as np
import pytest
import xarray as xr

from conftest import DATA_DIR, ERA5_FILES, VPUData
from river_route.runoff import (
    cell_xy_from_regular_grid,
    voronoi_diagram_from_regular_xy,
    compute_voronoi_catchment_intersects,
    grid_weights,
    grid_to_catchment,
)


# ═════════════════════════════════════════════════════════════════════════════
# cell_xy_from_regular_grid
# ═════════════════════════════════════════════════════════════════════════════

def test_cell_xy_from_era5():
    """Extract cell coordinates from a real ERA5 file."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    x, y = cell_xy_from_regular_grid(ERA5_FILES[0], x_var='longitude', y_var='latitude')
    assert x.ndim == 1
    assert y.ndim == 1
    assert len(x) > 0
    assert len(y) > 0
    # ERA5 global grid should span roughly -180 to 180 / -90 to 90
    assert x.min() >= -180 and x.max() <= 360
    assert y.min() >= -90 and y.max() <= 90


def test_cell_xy_synthetic():
    """Extract cell coordinates from a synthetic regular grid."""
    tmpdir = tempfile.mkdtemp()
    try:
        ds = xr.Dataset(
            {'data': xr.DataArray(np.zeros((3, 5)), dims=('lat', 'lon'))},
            coords={'lon': np.arange(5) * 0.25, 'lat': np.arange(3) * 0.25},
        )
        path = os.path.join(tmpdir, 'grid.nc')
        ds.to_netcdf(path)

        x, y = cell_xy_from_regular_grid(path, x_var='lon', y_var='lat')
        assert len(x) == 5
        assert len(y) == 3
        np.testing.assert_array_equal(x, np.arange(5) * 0.25)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ═════════════════════════════════════════════════════════════════════════════
# voronoi_diagram_from_regular_grid_cell_xy
# ═════════════════════════════════════════════════════════════════════════════

def test_voronoi_diagram_synthetic():
    """Build Voronoi polygons from a small synthetic grid."""
    x = np.array([0.0, 1.0, 2.0])
    y = np.array([0.0, 1.0])
    gdf = voronoi_diagram_from_regular_xy(x, y)
    assert 'geometry' in gdf.columns
    assert 'x_index' in gdf.columns
    assert 'y_index' in gdf.columns
    # Should have one polygon per grid cell
    assert len(gdf) == len(x) * len(y)


def test_voronoi_diagram_from_era5():
    """Build Voronoi polygons from ERA5 grid coordinates."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    x, y = cell_xy_from_regular_grid(ERA5_FILES[0], x_var='longitude', y_var='latitude')
    gdf = voronoi_diagram_from_regular_xy(x, y)
    assert 'geometry' in gdf.columns
    assert 'x_index' in gdf.columns
    assert 'y_index' in gdf.columns
    # One polygon per grid cell
    assert len(gdf) == len(x) * len(y)
    # CRS should default to EPSG:4326
    assert gdf.crs is not None
    assert gdf.crs.to_epsg() == 4326


# ═════════════════════════════════════════════════════════════════════════════
# compute_voronoi_catchment_intersects
# ═════════════════════════════════════════════════════════════════════════════

def test_compute_voronoi_catchment_intersects(vpu: VPUData):
    """Intersect Voronoi polygons with catchment boundaries."""
    import geopandas as gpd

    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    catchments_file = vpu.hydrography_dir / 'catchments_718.parquet'
    if not catchments_file.exists():
        pytest.skip(f'Missing {catchments_file}')

    x, y = cell_xy_from_regular_grid(ERA5_FILES[0], x_var='longitude', y_var='latitude')
    voronoi_gdf = voronoi_diagram_from_regular_xy(x, y)
    catchments_gdf = gpd.read_parquet(catchments_file)

    result = compute_voronoi_catchment_intersects(voronoi_gdf, catchments_gdf, river_id_variable='linkno')

    expected_cols = {'linkno', 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion'}
    assert expected_cols.issubset(set(result.columns)), \
        f'Missing columns: {expected_cols - set(result.columns)}'
    # Proportions should sum to ~1.0 per catchment
    prop_sums = result.groupby('linkno')['proportion'].sum()
    np.testing.assert_allclose(
        prop_sums.values, 1.0, atol=0.02,
        err_msg='Proportions do not sum to ~1.0 per catchment',
    )


# ═════════════════════════════════════════════════════════════════════════════
# grid_weights (end-to-end weight table generation)
# ═════════════════════════════════════════════════════════════════════════════

def test_grid_weights(vpu: VPUData):
    """Generate a weight table and compare against existing known-good weight table."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    catchments_file = vpu.hydrography_dir / 'catchments_718.parquet'
    if not catchments_file.exists():
        pytest.skip(f'Missing {catchments_file}')

    result = grid_weights(
        ERA5_FILES[0], str(catchments_file),
        x_var='longitude', y_var='latitude', river_id_var='linkno',
    )

    expected_cols = {'linkno', 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion'}
    assert expected_cols.issubset(set(result.columns)), \
        f'Missing columns: {expected_cols - set(result.columns)}'

    # Compare against known-good weight table
    known = xr.open_dataset(str(vpu.grid_weights_file))
    known_rivers = set(known['river_id'].values.tolist())
    result_rivers = set(result['linkno'].values.tolist())
    # Generated weights should cover the same rivers
    assert known_rivers == result_rivers, \
        f'{len(known_rivers - result_rivers)} rivers in known but not generated, ' \
        f'{len(result_rivers - known_rivers)} rivers generated but not in known'


# ═════════════════════════════════════════════════════════════════════════════
# grid_to_catchment
# ═════════════════════════════════════════════════════════════════════════════

def test_grid_to_catchment_as_volumes(vpu: VPUData):
    """Aggregate ERA5 gridded runoff to catchment volumes using the weight table."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    ds = grid_to_catchment(
        ERA5_FILES[0],
        weight_table=str(vpu.grid_weights_file),
        runoff_var='ro',
        x_var='longitude',
        y_var='latitude',
        time_var='valid_time',
        as_volumes=True,
    )
    assert 'runoff' in ds
    assert 'time' in ds.dims
    assert 'river_id' in ds.dims
    # Volumes should be non-negative (m³)
    assert np.all(ds['runoff'].values >= 0), 'Negative volumes found'


def test_grid_to_catchment_as_depths(vpu: VPUData):
    """Aggregate ERA5 gridded runoff to catchment depths (for UnitMuskingum)."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    ds = grid_to_catchment(
        ERA5_FILES[0],
        weight_table=str(vpu.grid_weights_file),
        runoff_var='ro',
        x_var='longitude',
        y_var='latitude',
        time_var='valid_time',
        as_volumes=False,
    )
    assert 'runoff' in ds
    assert 'time' in ds.dims
    assert 'river_id' in ds.dims


def test_grid_to_catchment_volumes_vs_depths_ratio(vpu: VPUData):
    """Verify that volumes / depths ≈ catchment area for each river."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    kwargs = dict(
        runoff_var='ro', x_var='longitude', y_var='latitude', time_var='valid_time',
    )
    ds_vol = grid_to_catchment(
        ERA5_FILES[0], weight_table=str(vpu.grid_weights_file), as_volumes=True, **kwargs,
    )
    ds_dep = grid_to_catchment(
        ERA5_FILES[0], weight_table=str(vpu.grid_weights_file), as_volumes=False, **kwargs,
    )

    # Compute total catchment area per river from the weight table
    wt_df = xr.open_dataset(str(vpu.grid_weights_file))[['river_id', 'area_sqm']].to_dataframe()
    catchment_area = wt_df.groupby('river_id')['area_sqm'].sum()
    # Align to the river_id order in the output
    area = catchment_area.reindex(ds_vol['river_id'].values).values

    # volumes should equal depths * area (for non-zero depths)
    depths = ds_dep['runoff'].values
    volumes = ds_vol['runoff'].values
    expected_volumes = depths * area[np.newaxis, :]

    # Only check where depths are non-trivial to avoid division noise
    mask = np.abs(depths) > 1e-12
    np.testing.assert_allclose(
        volumes[mask], expected_volumes[mask],
        rtol=1e-6,
        err_msg='volumes != depths * catchment_area',
    )


def test_grid_to_catchment_cumulative_input(vpu: VPUData):
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
            weight_table=str(vpu.grid_weights_file),
            runoff_var='ro', x_var='longitude', y_var='latitude', time_var='valid_time',
            as_volumes=True,
        )

        ds_inc = grid_to_catchment(str(incremental_file), cumulative=False, **kwargs)
        ds_cum = grid_to_catchment(cumulative_file, cumulative=True, **kwargs)

        # Tolerance is loose because cumsum -> float32 storage -> diff loses precision
        # relative to direct incremental aggregation
        np.testing.assert_allclose(
            ds_inc['runoff'].values, ds_cum['runoff'].values,
            rtol=0.02, atol=0.15,
            err_msg='Cumulative input does not produce same volumes as incremental',
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
