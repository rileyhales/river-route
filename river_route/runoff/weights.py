import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely.geometry
import shapely.ops
import xarray as xr

from .._metadata import __version__
from ..types import PathInput

__all__ = [
    'cell_xy_from_regular_grid',
    'voronoi_diagram_from_regular_xy',
    'compute_voronoi_catchment_intersects',
    'grid_weights',
]

logger = logging.getLogger(__name__)


def cell_xy_from_regular_grid(
    dataset: PathInput, x_var: str = 'lon', y_var: str = 'lat'
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
        shapely.geometry.MultiPoint([shapely.geometry.Point(xi, yi) for xi, yi in zip(x_grid, y_grid, strict=True)])
    )

    logger.info('Adding attributes to voronoi polygons')
    voronoi_gdf = gpd.GeoDataFrame(geometry=[region for region in regions.geoms], crs=crs)
    voronoi_gdf['x'] = voronoi_gdf.geometry.apply(lambda geom: geom.centroid.x).astype(float)
    voronoi_gdf['y'] = voronoi_gdf.geometry.apply(lambda geom: geom.centroid.y).astype(float)
    voronoi_gdf['x_index'] = voronoi_gdf['x'].apply(lambda value: np.argmin(np.abs(x - value))).astype(int)
    voronoi_gdf['y_index'] = voronoi_gdf['y'].apply(lambda value: np.argmin(np.abs(y - value))).astype(int)
    return voronoi_gdf.sort_values(by=['x', 'y']).reset_index(drop=True)


def compute_voronoi_catchment_intersects(
    voronoi_gdf: gpd.GeoDataFrame,
    catchments_gdf: gpd.GeoDataFrame,
    save_path: PathInput | None = None,
    attributes: dict | None = None,
    river_id_variable: str = 'river_id',
) -> pd.DataFrame:
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
        df[[river_id_variable, 'area_sqm']]
        .groupby(river_id_variable)
        .sum()
        .rename(columns={'area_sqm': 'area_sqm_total'})
    )
    df = df.merge(total_area, left_on=river_id_variable, right_index=True, how='left')
    df['proportion'] = df['area_sqm'] / df['area_sqm_total']
    # computed in float64 so proportions come from exact areas, then stored as float32 to halve the table and every
    # catchment runoff array aggregated from it
    df = df.astype({column: np.float32 for column in ('x', 'y', 'area_sqm', 'area_sqm_total', 'proportion')})

    if save_path:
        (
            df.to_xarray()
            .assign_attrs(
                {
                    'description': 'proportions of runoff cells that intersect catchments for use with river',
                    'voronoi_gdf_crs': voronoi_gdf.crs.to_string(),
                    'catchments_gdf_crs': catchments_gdf.crs.to_string(),
                    'river_route_version': __version__,
                    **(attributes or {}),
                }
            )
            .to_netcdf(save_path)
        )
    return df


def grid_weights(
    grid_path: PathInput,
    catchments_path: PathInput,
    *,
    var_x: str = 'lon',
    var_y: str = 'lat',
    var_river_id: str = 'river_id',
    crs: int = 4326,
    save_voronoi_path: PathInput | None = None,
    save_weights_path: PathInput | None = None,
    routing_params_path: PathInput | None = None,
) -> pd.DataFrame:
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

    # Map x_index back to original grid indices (for use by GaussianGridRunoff)
    voronoi_gdf['x_index'] = sort_order[voronoi_gdf['x_index'].values]
    if save_voronoi_path:
        voronoi_gdf.to_parquet(save_voronoi_path)
    catchments_gdf = gpd.read_parquet(catchments_path)
    df = compute_voronoi_catchment_intersects(
        voronoi_gdf,
        catchments_gdf,
        save_path=None,
        attributes=dict(grid_path=str(grid_path), catchments_path=str(catchments_path)),
        river_id_variable=var_river_id,
    )

    if routing_params_path is not None:
        ordered_ids = pd.read_parquet(routing_params_path)[var_river_id].to_numpy()
        id_to_order = {int(rid): i for i, rid in enumerate(ordered_ids)}
        df = (
            df.assign(_sort_key=df[var_river_id].map(id_to_order))
            .sort_values(['_sort_key', 'area_sqm'], ascending=[True, False])
            .drop(columns='_sort_key')
            .reset_index(drop=True)
        )
    else:
        logger.warning('routing_params_path not provided; weight table row order may not match routing network order')

    if save_weights_path:
        (
            df[[var_river_id, 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion']]
            .to_xarray()
            .assign_attrs(
                {
                    'description': 'proportions of runoff cells that intersect river catchments',
                    'grid_path': str(grid_path),
                    'catchments_path': str(catchments_path),
                    'river_route_version': __version__,
                }
            )
            .to_netcdf(save_weights_path)
        )

    return df
