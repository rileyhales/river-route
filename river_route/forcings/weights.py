import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import shapely.ops
import xarray as xr
from pyproj import Transformer

from .._metadata import __version__
from ..types import PathInput
from .ECMWFGribReducedGrid import ReducedGaussianGrid

__all__ = [
    'cell_xy_from_regular_grid',
    'voronoi_diagram_from_regular_xy',
    'compute_voronoi_catchment_intersects',
    'grid_weights',
    'reduced_grid_weights',
]

logger = logging.getLogger(__name__)


def cell_xy_from_regular_grid(
    dataset: PathInput, var_x: str = 'lon', var_y: str = 'lat'
) -> tuple[np.ndarray, np.ndarray]:
    """
    The cell center x and y coordinates of a regular grid.

    Args:
        dataset: netCDF, or any file xarray opens, holding the grid's 1D x and y coordinate variables
        var_x: name of the x coordinate variable
        var_y: name of the y coordinate variable

    Returns:
        tuple: (x, y), the 1D coordinate arrays

    Raises:
        KeyError: if either variable is missing
        ValueError: if either coordinate is not 1D
    """
    with xr.open_dataset(dataset) as ds:
        if var_x not in ds.variables:
            raise KeyError(f'{var_x} must be a variable in {dataset}')
        if var_y not in ds.variables:
            raise KeyError(f'{var_y} must be a variable in {dataset}')
        x = ds[var_x].to_numpy()
        y = ds[var_y].to_numpy()

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError('Regular grid requires 1D x/y coordinate arrays')
    return x, y


def voronoi_diagram_from_regular_xy(x: np.ndarray, y: np.ndarray, crs: int = 4326) -> gpd.GeoDataFrame:
    """
    The Voronoi polygon around the center of each cell of a regular grid, with the cell's center and its position in
    the x and y coordinates.

    Args:
        x: (n_x,) cell center x coordinates
        y: (n_y,) cell center y coordinates
        crs: EPSG code of the coordinates

    Returns:
        gpd.GeoDataFrame: one polygon per cell with columns x, y, x_index, and y_index, sorted by x and y

    Raises:
        ValueError: if x or y is not a non-empty 1D array
    """
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
    regions = shapely.ops.voronoi_diagram(shapely.multipoints(shapely.points(x_grid, y_grid)))

    logger.info('Adding attributes to voronoi polygons')
    voronoi_gdf = gpd.GeoDataFrame(geometry=list(regions.geoms), crs=crs)
    centroids = shapely.centroid(voronoi_gdf.geometry.array)
    voronoi_gdf['x'] = shapely.get_x(centroids)
    voronoi_gdf['y'] = shapely.get_y(centroids)
    voronoi_gdf['x_index'] = voronoi_gdf['x'].apply(lambda value: np.argmin(np.abs(x - value))).astype(int)
    voronoi_gdf['y_index'] = voronoi_gdf['y'].apply(lambda value: np.argmin(np.abs(y - value))).astype(int)
    return voronoi_gdf.sort_values(by=['x', 'y']).reset_index(drop=True)


def compute_voronoi_catchment_intersects(
    voronoi_gdf: gpd.GeoDataFrame,
    catchments_gdf: gpd.GeoDataFrame,
    save_path: PathInput | None = None,
    attributes: dict | None = None,
) -> pd.DataFrame:
    """
    The weight table of the intersections between grid cell Voronoi polygons and catchments: the area of each
    catchment inside each cell and that area's proportion of the catchment.

    Args:
        voronoi_gdf: cell polygons with columns x_index, y_index, x, and y, as from voronoi_diagram_from_regular_xy
        catchments_gdf: catchment polygons with a riverId column
        save_path: optional netCDF file to save the weight table to
        attributes: optional attributes added to the saved file

    Returns:
        pd.DataFrame: columns riverId, x_index, y_index, x, y, area_sqm, area_sqm_total, and proportion

    Raises:
        KeyError: if a required column is missing
    """
    if 'riverId' not in catchments_gdf.columns:
        raise KeyError('catchments_gdf must contain a riverId column')
    if not {'x_index', 'y_index', 'x', 'y'}.issubset(voronoi_gdf.columns):
        raise KeyError('voronoi_gdf must include x_index, y_index, x, and y columns')

    logger.info('Performing overlay operation')
    intersections = gpd.overlay(voronoi_gdf, catchments_gdf, how='intersection')
    logger.info('Calculating area of intersections')
    intersections['area_sqm'] = intersections.geometry.to_crs({'proj': 'cea'}).area
    cell_columns = ('x_index', 'y_index', 'x', 'y')
    df = _proportions_of_catchment_areas(intersections, cell_columns)

    if save_path:
        (
            df.to_xarray()
            .assign_attrs(
                {
                    'description': 'proportions of runoff cells that intersect catchments for use with river-route',
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
    crs: int = 4326,
    save_voronoi_path: PathInput | None = None,
    save_weights_path: PathInput | None = None,
    network_path: PathInput | None = None,
) -> pd.DataFrame:
    """
    Compute the grid weights for a given grid and catchments.

    Args:
        grid_path: path to a NetCDF file containing the grid's 1D x and y coordinate variables
        catchments_path: path to a GeoParquet file containing the catchment geometries and a riverId column
        var_x: x-coordinate variable name in the grid file
        var_y: y-coordinate variable name in the grid file
        crs: EPSG code for the grid coordinate reference system (default: 4326)
        save_voronoi_path: optional path to save the Voronoi polygons as a GeoParquet file
        save_weights_path: optional path to save the grid weights as a NetCDF file
        network_path: optional path to a network file whose riverId column order is used to
            topologically sort the weight table rows. When omitted the row order is spatial (not topological).

    Returns:
        pd.DataFrame: a DataFrame containing the grid weights with columns
            ['riverId', 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'area_sqm_total', 'proportion']. The saved
            file omits area_sqm_total.
    """
    x, y = cell_xy_from_regular_grid(grid_path, var_x=var_x, var_y=var_y)

    # Convert 0-360 longitudes to -180..180 for overlay with catchments
    x_geo = x.copy()
    x_geo[x_geo > 180] -= 360
    sort_order = np.argsort(x_geo)
    x_geo = x_geo[sort_order]

    voronoi_gdf = voronoi_diagram_from_regular_xy(x_geo, y, crs=crs)

    # Map x_index back to original grid indices (for use by GridRunoff)
    voronoi_gdf['x_index'] = sort_order[voronoi_gdf['x_index'].to_numpy()]
    if save_voronoi_path:
        voronoi_gdf.to_parquet(save_voronoi_path)
    catchments_gdf = gpd.read_parquet(catchments_path)
    df = compute_voronoi_catchment_intersects(
        voronoi_gdf,
        catchments_gdf,
        save_path=None,
        attributes={'grid_path': str(grid_path), 'catchments_path': str(catchments_path)},
    )
    return _order_and_save_weight_table(
        df,
        cell_columns=('x_index', 'y_index', 'x', 'y'),
        network_path=network_path,
        save_weights_path=save_weights_path,
        attributes={'grid_path': str(grid_path), 'catchments_path': str(catchments_path)},
    )


def reduced_grid_weights(
    grib_path: PathInput,
    catchments_path: PathInput,
    *,
    save_weights_path: PathInput | None = None,
    network_path: PathInput | None = None,
) -> pd.DataFrame:
    """
    Compute the grid weights of the cells of a reduced gaussian grid GRIB file, read by ``ECMWFGribReducedGrid``: the
    area of each catchment inside each cell polygon of ``ReducedGaussianGrid.cell_polygons`` that overlaps it.

    Every part of a cell polygon is a rectangle of longitude and latitude, and stays a rectangle in the cylindrical
    equal-area projection, so the catchments are projected to it once and cut to each cell with the rectangle clipping
    of GEOS, whose pieces are measured where they are. A general polygon overlay followed by projecting every piece
    took 40 s for the 58,358 catchments of the Columbia (39.9 million vertices) on O1280, and this takes 6 s.

    Args:
        grib_path: GRIB file whose first message's grid is used
        catchments_path: GeoParquet file of the catchment polygons with a riverId column, in any CRS
        save_weights_path: optional path to save the grid weights as a netCDF file
        network_path: optional path to a network file whose riverId column order is used to
            topologically sort the weight table rows. When omitted the row order is spatial (not topological).

    Returns:
        pd.DataFrame: the grid weights with columns [riverId, cell_index, x, y, area_sqm, area_sqm_total, proportion]
    """
    grid = ReducedGaussianGrid.from_grib(grib_path)
    catchments_gdf = gpd.read_parquet(catchments_path, columns=['riverId', 'geometry']).to_crs({'proj': 'cea'})
    to_lonlat = Transformer.from_crs(catchments_gdf.crs, 4326, always_xy=True)
    min_x, min_y, max_x, max_y = catchments_gdf.total_bounds
    cells_gdf = grid.cell_polygons(bounds=(*to_lonlat.transform(min_x, min_y), *to_lonlat.transform(max_x, max_y)))

    parts, cell_of_part = shapely.get_parts(cells_gdf.geometry.array, return_index=True)
    west, south, east, north = shapely.bounds(parts).T
    to_equal_area = Transformer.from_crs(4326, catchments_gdf.crs, always_xy=True)
    west, south = to_equal_area.transform(west, south)
    east, north = to_equal_area.transform(east, north)
    catchment_geometry = catchments_gdf.geometry.array
    catchment, part = shapely.STRtree(shapely.box(west, south, east, north)).query(catchment_geometry)
    # grouped by cell part, so each part clips all of its catchments in one call
    by_part = np.argsort(part, kind='stable')
    catchment, part = catchment[by_part], part[by_part]
    starts = np.flatnonzero(np.diff(part, prepend=-1))
    stops = np.append(starts[1:], part.shape[0])
    area_sqm = np.empty(part.shape[0])
    for start, stop in zip(starts, stops, strict=True):
        p = part[start]
        pieces = shapely.clip_by_rect(catchment_geometry[catchment[start:stop]], west[p], south[p], east[p], north[p])
        area_sqm[start:stop] = shapely.area(pieces)

    cell_row = cell_of_part[part]
    pieces_df = pd.DataFrame(
        {
            'riverId': catchments_gdf['riverId'].to_numpy()[catchment],
            'cell_index': cells_gdf['cell_index'].to_numpy()[cell_row],
            'x': cells_gdf['x'].to_numpy()[cell_row],
            'y': cells_gdf['y'].to_numpy()[cell_row],
            'area_sqm': area_sqm,
        }
    )
    df = _proportions_of_catchment_areas(pieces_df[area_sqm > 0], ('cell_index', 'x', 'y'))
    return _order_and_save_weight_table(
        df,
        cell_columns=('cell_index', 'x', 'y'),
        network_path=network_path,
        save_weights_path=save_weights_path,
        attributes={'grid_path': str(grib_path), 'catchments_path': str(catchments_path)},
    )


def _proportions_of_catchment_areas(pieces: pd.DataFrame, cell_columns: tuple[str, ...]) -> pd.DataFrame:
    """Sum the area_sqm of the pieces of each catchment in each cell, then each cell's proportion of the catchment."""
    df = (
        pieces[['riverId', *cell_columns, 'area_sqm']]
        .groupby(['riverId', *cell_columns], as_index=False)
        .agg({'area_sqm': 'sum'})
        .sort_values(['riverId', 'area_sqm'], ascending=[True, False])
        .reset_index(drop=True)
    )
    total_area = df[['riverId', 'area_sqm']].groupby('riverId').sum().rename(columns={'area_sqm': 'area_sqm_total'})
    df = df.merge(total_area, left_on='riverId', right_index=True, how='left')
    df['proportion'] = df['area_sqm'] / df['area_sqm_total']
    # computed in float64 so proportions come from exact areas, then stored as float32 to halve the table and every
    # catchment runoff array aggregated from it
    return df.astype(dict.fromkeys(('x', 'y', 'area_sqm', 'area_sqm_total', 'proportion'), np.float32))


def _order_and_save_weight_table(
    df: pd.DataFrame,
    *,
    cell_columns: tuple[str, ...],
    network_path: PathInput | None,
    save_weights_path: PathInput | None,
    attributes: dict,
) -> pd.DataFrame:
    """Sort the weight table's rivers into the order of the network file, then save it when a path is given."""
    if network_path is not None:
        ordered_ids = pd.read_parquet(network_path, columns=['riverId'])['riverId'].to_numpy()
        id_to_order = {int(rid): i for i, rid in enumerate(ordered_ids)}
        df = (
            df.assign(_sort_key=df['riverId'].map(id_to_order))
            .sort_values(['_sort_key', 'area_sqm'], ascending=[True, False])
            .drop(columns='_sort_key')
            .reset_index(drop=True)
        )
    else:
        logger.warning('network_path not provided; weight table row order may not match routing network order')

    if save_weights_path:
        (
            df[['riverId', *cell_columns, 'area_sqm', 'proportion']]
            .to_xarray()
            .assign_attrs(
                {
                    'description': 'proportions of runoff cells that intersect river catchments',
                    **attributes,
                    'river_route_version': __version__,
                }
            )
            .to_netcdf(save_weights_path)
        )
    return df
