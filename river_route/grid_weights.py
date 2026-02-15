import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely.geometry
import shapely.ops
import xarray as xr
from .__metadata__ import __version__

logger = logging.getLogger(__name__)

__all__ = [
    'get_cell_xy_from_regular_grid',
    'get_cell_xy_from_reduced_grid',
    'voroni_diagram_from_cell_xy',
    'grid_weights_table',
]


def get_cell_xy_from_regular_grid(
        dataset_path: str,
        x_var: str = 'lon',
        y_var: str = 'lat',
) -> tuple[np.ndarray, np.ndarray]:
    """
    Get flattened x and y coordinates of cell centers from a regular grid.
    """
    with xr.open_dataset(dataset_path) as ds:
        if x_var not in ds.variables:
            raise KeyError(f'{x_var} must be a variable in {dataset_path}')
        if y_var not in ds.variables:
            raise KeyError(f'{y_var} must be a variable in {dataset_path}')
        x = ds[x_var].values
        y = ds[y_var].values

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError('Regular grid requires 1D x/y coordinate arrays')

    x_grid, y_grid = np.meshgrid(x, y)
    return x_grid.flatten(), y_grid.flatten()


def get_cell_xy_from_reduced_grid(
        dataset_path: str,
        x_var: str = 'lon',
        y_var: str = 'lat',
) -> tuple[np.ndarray, np.ndarray]:
    """
    Get flattened x and y coordinates of cell centers from a reduced/unstructured grid.
    """
    with xr.open_dataset(dataset_path) as ds:
        if x_var not in ds.variables or y_var not in ds.variables:
            raise KeyError(f'{x_var} and {y_var} must be variables in {dataset_path}')
        x = ds[x_var].values
        y = ds[y_var].values

    if x.ndim == 1 and y.ndim == 1:
        if x.shape != y.shape:
            raise ValueError('Reduced grid requires x and y arrays with the same shape')
        return x.flatten(), y.flatten()

    if x.ndim == 2 and y.ndim == 2:
        if x.shape != y.shape:
            raise ValueError('Reduced grid requires x and y arrays with the same shape')
        return x.flatten(), y.flatten()

    raise ValueError('Reduced grid expects either both 1D arrays of equal length or both 2D arrays of equal shape')


def voroni_diagram_from_cell_xy(x: np.ndarray, y: np.ndarray, crs: int = 4326) -> gpd.GeoDataFrame:
    """
    Create a GeoDataFrame of Voroni polygons around the center of each cell in a grid.
    """
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError('x and y must be 1D arrays')
    if x.shape[0] != y.shape[0]:
        raise ValueError('x and y must have the same number of points')
    if x.shape[0] == 0:
        raise ValueError('x and y cannot be empty')

    logger.info('Creating Voroni polygons')
    regions = shapely.ops.voronoi_diagram(
        shapely.geometry.MultiPoint([shapely.geometry.Point(xi, yi) for xi, yi in zip(x, y)])
    )

    logger.info('Adding attributes to voroni polygons')
    voroni_gdf = gpd.GeoDataFrame(geometry=[region for region in regions.geoms], crs=crs)
    voroni_gdf['x'] = voroni_gdf.geometry.apply(lambda geom: geom.centroid.x).astype(float)
    voroni_gdf['y'] = voroni_gdf.geometry.apply(lambda geom: geom.centroid.y).astype(float)
    voroni_gdf['x_index'] = voroni_gdf['x'].apply(lambda value: np.argmin(np.abs(x - value))).astype(int)
    voroni_gdf['y_index'] = voroni_gdf['y'].apply(lambda value: np.argmin(np.abs(y - value))).astype(int)
    return voroni_gdf


def grid_weights_table(voroni_gdf: gpd.GeoDataFrame, catchments_gdf: gpd.GeoDataFrame,
                       save_path: str | None = None, save_attributes: dict | None = None,
                       river_id_variable: str = 'river_id') -> pd.DataFrame:
    """
    Create a table of intersections between Voroni polygons and catchments.
    """
    if river_id_variable not in catchments_gdf.columns:
        raise KeyError('catchments_gdf must contain a river_id column')
    if not {'x_index', 'y_index', 'x', 'y'}.issubset(voroni_gdf.columns):
        raise KeyError('voroni_gdf must include x_index, y_index, x, and y columns')

    logger.info('Performing overlay operation')
    intersections = gpd.overlay(voroni_gdf, catchments_gdf, how='intersection')
    logger.info('Calculating area of intersections')
    intersections['area_sqm'] = intersections.geometry.to_crs({'proj': 'cea'}).area

    df = (
        intersections[['river_id', 'x_index', 'y_index', 'x', 'y', 'area_sqm']]
        .groupby(['river_id', 'x_index', 'y_index', 'x', 'y'], as_index=False)
        .agg({'area_sqm': 'sum'})
        .sort_values(['river_id', 'area_sqm'], ascending=[True, False])
        .reset_index(drop=True)
    )
    total_area = df[['river_id', 'area_sqm']].groupby('river_id').sum().rename(columns={'area_sqm': 'area_sqm_total'})
    df = df.merge(total_area, left_on='river_id', right_index=True, how='left')
    df['proportion'] = df['area_sqm'] / df['area_sqm_total']

    if save_path:
        save_attributes = {
            'description': 'proportions of runoff cells that intersect river catchments to be used with river-route',
            'voroni_gdf_crs': voroni_gdf.crs.to_string(),
            'catchments_gdf_crs': catchments_gdf.crs.to_string(),
            'river_route_version': __version__,
            **(save_attributes or {}),
        }
        ds = xr.Dataset({
            'river_id': ('index', df['river_id'].to_numpy(dtype=np.int64, copy=False)),
            'x_index': ('index', df['x_index'].to_numpy(dtype=np.int64, copy=False)),
            'y_index': ('index', df['y_index'].to_numpy(dtype=np.int64, copy=False)),
            'x': ('index', df['x'].to_numpy(dtype=np.float64, copy=False)),
            'y': ('index', df['y'].to_numpy(dtype=np.float64, copy=False)),
            'area_sqm': ('index', df['area_sqm'].to_numpy(dtype=np.float64, copy=False)),
            'proportion': ('index', df['proportion'].to_numpy(dtype=np.float64, copy=False)),
        })
        ds.attrs.update(save_attributes)
        ds.to_netcdf(save_path)
    return df
