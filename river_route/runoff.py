import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import scipy.sparse
import shapely.geometry
import shapely.ops
import xarray as xr

from ._metadata import __version__
from .types import PathInput

__all__ = [
    'cell_xy_from_regular_grid',
    'voronoi_diagram_from_regular_xy',
    'compute_voronoi_catchment_intersects',
    'grid_weights',
    'runoff_to_qlateral',
]

logger = logging.getLogger(__name__)


def cell_xy_from_regular_grid(
        dataset: PathInput, x_var: str = 'lon', y_var: str = 'lat',
) -> tuple[np.ndarray, np.ndarray]:
    """Get cell center x and y coordinates from a regular grid common dataset structure."""
    with xr.open_dataset(dataset) as ds:
        if x_var not in ds.variables:
            raise KeyError(f'{x_var} must be a variable in {dataset}')
        if y_var not in ds.variables:
            raise KeyError(f'{y_var} must be a variable in {dataset}')
        x = ds[x_var].values
        y = ds[y_var].values

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError('Regular grid requires 1D x/y coordinate arrays')
    return x, y


def voronoi_diagram_from_regular_xy(x: np.ndarray, y: np.ndarray, crs: int = 4326) -> gpd.GeoDataFrame:
    """Create a GeoDataFrame of Voronoi polygons around the center of each cell in a grid."""
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError('x and y must be 1D arrays')
    if x.shape[0] == 0 or y.shape[0] == 0:
        raise ValueError('x and y cannot be empty')

    x_grid, y_grid = np.meshgrid(x, y)
    x_grid = x_grid.flatten()
    y_grid = y_grid.flatten()

    if x_grid.shape[0] != y_grid.shape[0]:
        raise ValueError('x and y must have the same number of points')

    logger.info('Creating Voronoi polygons')
    regions = shapely.ops.voronoi_diagram(
        shapely.geometry.MultiPoint([shapely.geometry.Point(xi, yi) for xi, yi in zip(x_grid, y_grid)])
    )

    logger.info('Adding attributes to voronoi polygons')
    voronoi_gdf = gpd.GeoDataFrame(geometry=[region for region in regions.geoms], crs=crs)
    voronoi_gdf['x'] = voronoi_gdf.geometry.apply(lambda geom: geom.centroid.x).astype(float)
    voronoi_gdf['y'] = voronoi_gdf.geometry.apply(lambda geom: geom.centroid.y).astype(float)
    voronoi_gdf['x_index'] = voronoi_gdf['x'].apply(lambda value: np.argmin(np.abs(x - value))).astype(int)
    voronoi_gdf['y_index'] = voronoi_gdf['y'].apply(lambda value: np.argmin(np.abs(y - value))).astype(int)
    return voronoi_gdf.sort_values(by=['x', 'y']).reset_index(drop=True)


def compute_voronoi_catchment_intersects(voronoi_gdf: gpd.GeoDataFrame, catchments_gdf: gpd.GeoDataFrame,
                                         save_path: PathInput | None = None, attributes: dict | None = None,
                                         river_id_variable: str = 'river_id') -> pd.DataFrame:
    """
    Create a table of intersections between Voronoi polygons and catchments.
    """
    if river_id_variable not in catchments_gdf.columns:
        raise KeyError(f'catchments_gdf must contain a {river_id_variable} column')
    if not {'x_index', 'y_index', 'x', 'y'}.issubset(voronoi_gdf.columns):
        raise KeyError('voronoi_gdf must include x_index, y_index, x, and y columns')

    logger.info('Performing overlay operation')
    intersections = gpd.overlay(voronoi_gdf, catchments_gdf, how='intersection')
    logger.info('Calculating area of intersections')
    intersections['area_sqm'] = intersections.geometry.to_crs({'proj': 'cea'}).area

    df = (
        intersections[[river_id_variable, 'x_index', 'y_index', 'x', 'y', 'area_sqm']]
        .groupby([river_id_variable, 'x_index', 'y_index', 'x', 'y'], as_index=False)
        .agg({'area_sqm': 'sum'})
        .sort_values([river_id_variable, 'area_sqm'], ascending=[True, False])
        .reset_index(drop=True)
    )
    total_area = (
        df
        [[river_id_variable, 'area_sqm']]
        .groupby(river_id_variable)
        .sum()
        .rename(columns={'area_sqm': 'area_sqm_total'})
    )
    df = df.merge(total_area, left_on=river_id_variable, right_index=True, how='left')
    df['proportion'] = df['area_sqm'] / df['area_sqm_total']

    if save_path:
        (
            df
            .to_xarray()
            .assign_attrs({
                'description': 'proportions of runoff cells that intersect catchments for use with river',
                'voronoi_gdf_crs': voronoi_gdf.crs.to_string(),
                'catchments_gdf_crs': catchments_gdf.crs.to_string(),
                'river_route_version': __version__,
                **(attributes or {}),
            })
            .to_netcdf(save_path)
        )
    return df


def grid_weights(grid_path: PathInput, catchments_path: PathInput, *,
                 var_x: str = 'lon', var_y: str = 'lat', var_river_id: str = 'river_id', crs: int = 4326,
                 save_voronoi_path: PathInput | None = None,
                 save_weights_path: PathInput | None = None,
                 routing_params_path: PathInput | None = None) -> pd.DataFrame:
    """
    Compute the grid weights for a given grid and catchments.

    Args:
        grid_path: path to a NetCDF file containing the grid information (must include 'lon' and 'lat' variables)
        catchments_path: path to a GeoParquet file containing the catchment geometries (must include 'river_id' column)
        var_x: x-coordinate variable name in the grid file
        var_y: y-coordinate variable name in the grid file
        var_river_id: variable name for river ID in the catchments file
        crs: EPSG code for the grid coordinate reference system (default: 4326)
        save_voronoi_path: optional path to save the Voronoi polygons as a GeoParquet file
        save_weights_path: optional path to save the grid weights as a NetCDF file
        routing_params_path: optional path to a routing params parquet file whose river_id column order is used to
            topologically sort the weight table rows. When omitted the row order is spatial (not topological).

    Returns:
        pd.DataFrame: a DataFrame containing the grid weights with columns
            ['river_id', 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion']
    """
    x, y = cell_xy_from_regular_grid(grid_path, x_var=var_x, y_var=var_y)

    # Convert 0-360 longitudes to -180..180 for overlay with catchments
    x_geo = x.copy()
    x_geo[x_geo > 180] -= 360
    sort_order = np.argsort(x_geo)
    x_geo = x_geo[sort_order]

    voronoi_gdf = voronoi_diagram_from_regular_xy(x_geo, y, crs=crs)

    # Map x_index back to original grid indices (for use by grid_to_qlateral)
    voronoi_gdf['x_index'] = sort_order[voronoi_gdf['x_index'].values]
    if save_voronoi_path:
        voronoi_gdf.to_parquet(save_voronoi_path)
    catchments_gdf = gpd.read_parquet(catchments_path)
    df = compute_voronoi_catchment_intersects(
        voronoi_gdf, catchments_gdf, save_path=None,
        attributes=dict(grid_path=str(grid_path), catchments_path=str(catchments_path)),
        river_id_variable=var_river_id
    )

    if routing_params_path is not None:
        ordered_ids = pd.read_parquet(routing_params_path)[var_river_id].to_numpy()
        id_to_order = {int(rid): i for i, rid in enumerate(ordered_ids)}
        df = (
            df
            .assign(_sort_key=df[var_river_id].map(id_to_order))
            .sort_values(['_sort_key', 'area_sqm'], ascending=[True, False])
            .drop(columns='_sort_key')
            .reset_index(drop=True)
        )
    else:
        logger.warning('routing_params_path not provided; weight table row order may not match routing network order')

    if save_weights_path:
        (
            df
            [[var_river_id, 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion']]
            .to_xarray()
            .assign_attrs({
                'description': 'proportions of runoff cells that intersect river catchments',
                'grid_path': str(grid_path),
                'catchments_path': str(catchments_path),
                'river_route_version': __version__,
            })
            .to_netcdf(save_weights_path)
        )

    return df


def _cumulative_to_incremental(df) -> pd.DataFrame:
    return pd.DataFrame(
        np.vstack([df.values[0, :], np.diff(df.values, axis=0)]),
        index=df.index,
        columns=df.columns
    )


def _incremental_to_cumulative(df) -> pd.DataFrame:
    return df.cumsum()


def _get_conversion_factor(unit: str) -> int | float:
    if unit is None:
        logger.warning("No units attribute found. Assuming meters")
        return 1
    if unit in ('m', 'meters', 'kg m-2'):
        return 1
    elif unit in ('mm', 'millimeters'):
        return .001
    else:
        raise ValueError(f"Unknown units: {unit}")


def runoff_to_qlateral(
        runoff_data: PathInput | list[PathInput],
        grid_weights_file: PathInput,
        *,
        var_runoff: str = 'ro',
        var_x: str = 'lon',
        var_y: str = 'lat',
        var_t: str = 'time',
        var_river_id: str = 'river_id',
        runoff_depth_unit: str | None = None,
        cumulative: bool = False,
        force_positive_runoff: bool = False,
        force_uniform_timesteps: bool = True,
        as_volumes: bool = False,
) -> xr.Dataset:
    """
    Aggregates gridded runoff depths to catchment level qlateral as depths or volumes.
    The core computation is an area weighted average of the cell's runoff value.

    Args:
        runoff_data (str | list[str]): path(s) to runoff files
        grid_weights_file (str): path to the weight table netCDF produced by ``grid_weights()``
        var_runoff (str): runoff variable name in the LSM files
        var_x (str): x-coordinate variable name in the LSM files
        var_y (str): y-coordinate variable name in the LSM files
        var_t (str): time variable name in the LSM files
        var_river_id (str): river ID variable name in the weight table and parameters file
        runoff_depth_unit (str): unit of the depth values; checked for in file attributes, defaulting to meters
        cumulative (bool): whether the runoff data is cumulative; converted to incremental if True
        force_positive_runoff (bool): clip negative runoff values to zero
        force_uniform_timesteps (bool): resample to a uniform timestep if the input is irregular
        as_volumes (bool): if True, return volumes (m³) instead of depths (m)

    Returns:
        xr.Dataset: qlateral with dimensions ``time`` and ``river_id``.
            Contains a single variable ``qlateral`` in metres or m³.
    """
    with xr.open_dataset(grid_weights_file) as ds:
        weight_df = ds[[var_river_id, 'x_index', 'y_index', 'proportion', 'area_sqm']].to_dataframe()
    unique_indexes = (
        weight_df
        [['x_index', 'y_index']]
        .drop_duplicates()
        .reset_index(drop=True)
        .reset_index()
        .astype(int)
    )
    unique_sorted_rivers = weight_df[[var_river_id, ]].drop_duplicates().sort_index()  # index already topo sorted

    with xr.open_mfdataset(runoff_data) as ds:
        runoff_depth_unit = runoff_depth_unit or ds[var_runoff].attrs.get('units', 'm')
        conversion_factor = _get_conversion_factor(runoff_depth_unit)
        runoff_raw = (
            ds
            [var_runoff]
            .isel({
                var_x: xr.DataArray(unique_indexes['x_index'].values, dims="points"),
                var_y: xr.DataArray(unique_indexes['y_index'].values, dims="points")
            })
            .transpose(var_t, "points")
            .values
        )
        time_index = ds[var_t].to_numpy()

    # Build sparse weight matrix W: (n_rivers, n_unique_points)
    point_idx = (
        weight_df[['x_index', 'y_index']]
        .merge(unique_indexes, on=['x_index', 'y_index'], how='left')
        ['index'].values
    )
    river_ids_ordered = unique_sorted_rivers[var_river_id].values
    river_id_to_row = pd.Series(np.arange(len(river_ids_ordered)), index=river_ids_ordered)
    river_idx = river_id_to_row.loc[weight_df[var_river_id].values].values

    W = scipy.sparse.csr_matrix(
        (weight_df['proportion'].values * conversion_factor, (river_idx, point_idx)),
        shape=(len(river_ids_ordered), len(unique_indexes)),
    )

    # Area-weighted aggregation via sparse matrix multiply, then free the grid data
    qlateral = np.asarray(W @ runoff_raw.T).T  # (time, n_rivers)
    del runoff_raw

    catchment_area = (
        weight_df
        .groupby(var_river_id)['area_sqm']
        .sum()
        .reindex(unique_sorted_rivers[var_river_id].values)
        .to_numpy()
    )

    # In-place conversions to avoid allocating new arrays
    if cumulative:
        for i in range(qlateral.shape[0] - 1, 0, -1):
            qlateral[i] -= qlateral[i - 1]
    if force_positive_runoff:
        np.clip(qlateral, 0, None, out=qlateral)

    time_diff = np.diff(time_index)
    if not np.all(time_diff == time_index[1] - time_index[0]) and force_uniform_timesteps:
        timestep = int((time_index[1] - time_index[0]) / np.timedelta64(1, 's'))
        logger.warning(f'Time steps are not uniform, resampling to the first timestep: {timestep} seconds')
        df = pd.DataFrame(qlateral, index=time_index, columns=river_ids_ordered)
        df = (
            _incremental_to_cumulative(df)
            .resample(rule=f'{timestep}s')
            .interpolate(method='linear')
        )
        df = _cumulative_to_incremental(df)
        time_index = df.index.values
        qlateral = df.to_numpy(dtype=np.float64)
        del df

    mask = np.isnan(qlateral)
    if mask.any():
        qlateral[mask] = 0.0

    if as_volumes:
        qlateral *= catchment_area[np.newaxis, :]
        units = 'm3'
        long_name = 'Incremental qlateral volumes'
    else:
        units = 'm'
        long_name = 'Incremental qlateral depths'

    start_date = pd.Timestamp(time_index[0]).strftime('%Y%m%d%H')
    end_date = pd.Timestamp(time_index[-1]).strftime('%Y%m%d%H')
    timestep = int((time_index[1] - time_index[0]) / np.timedelta64(1, 's')) if len(time_index) > 1 else 0
    return xr.Dataset(
        {
            'qlateral': xr.DataArray(
                qlateral,
                dims=('time', 'river_id'),
                attrs={'long_name': long_name, 'units': units},
            ),
        },
        coords={
            'river_id': xr.DataArray(
                river_ids_ordered.astype(np.int64, copy=False),
                dims=('river_id',),
                attrs={'long_name': 'unique ID number for each river'},
            ),
            'time': xr.DataArray(
                time_index,
                dims=('time',),
                attrs={
                    'long_name': 'time',
                    'standard_name': 'time',
                    'axis': 'T',
                    'time_step': f'{timestep}',
                },
            ),
        },
        attrs={
            'title': f'Incremental qlateral {long_name.split()[-1]}',
            'description': f'Incremental qlateral ({units}) for each river',
            'source': f'river-route v{__version__}',
            'history': f'Created on {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}',
            'suggested_file_name': f'qlateral_{start_date}_{end_date}.nc',
        },
    )
