import logging
from typing import NamedTuple, Self

import eccodes
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import shapely.ops
import xarray as xr
from pyproj import Transformer

from .._metadata import __version__
from ..types import IntArray, PathInput

__all__ = [
    'cell_xy_from_regular_grid',
    'voronoi_diagram_from_regular_xy',
    'compute_voronoi_catchment_intersects',
    'grid_weights',
    'ReducedGaussianGrid',
    'reduced_grid_weights',
]

logger = logging.getLogger(__name__)

EARTH_RADIUS_M = 6371229.0  # radius of the sphere the ECMWF IFS models the earth as


class ReducedGaussianGrid(NamedTuple):
    """
    A global reduced gaussian grid: 2N rows of cells at the gaussian latitudes, north to south, each row's cells evenly
    spaced in longitude eastward from 0 degrees. Cells are numbered row after row, the order of the values in a GRIB
    message.
    """

    N: int  # rows of cells between a pole and the equator
    cells_per_row: IntArray  # (2N,) number of cells on each row, north to south, the GRIB pl array

    @classmethod
    def from_grib(cls, grib_path: PathInput) -> Self:
        """Read the grid from the metadata of the first message in a GRIB file."""
        with open(grib_path, 'rb') as file:
            handle = eccodes.codes_grib_new_from_file(file)
        if handle is None:
            raise ValueError(f'{grib_path} holds no GRIB messages')
        try:
            if eccodes.codes_get(handle, 'gridType') != 'reduced_gg':
                raise ValueError(f'{grib_path} is on a {eccodes.codes_get(handle, "gridType")} grid, not reduced_gg')
            grid = cls(eccodes.codes_get(handle, 'N'), eccodes.codes_get_array(handle, 'pl'))
            numbered_as_described = (
                eccodes.codes_get(handle, 'global') == 1
                and eccodes.codes_get(handle, 'iScansNegatively') == 0
                and eccodes.codes_get(handle, 'jScansPositively') == 0
                and eccodes.codes_get(handle, 'numberOfDataPoints') == grid.cells_per_row.sum()
            )
        finally:
            eccodes.codes_release(handle)
        if not numbered_as_described:
            raise ValueError(f'{grib_path} is not a global grid numbered north to south and west to east')
        if grid.cells_per_row.shape[0] != 2 * grid.N:
            raise ValueError(f'{grib_path} has {grid.cells_per_row.shape[0]} rows of cells, expected {2 * grid.N}')
        return grid

    @property
    def n_cells(self) -> int:
        return int(self.cells_per_row.sum())

    @property
    def spacing_km(self) -> float:
        """Side of a square with the mean cell area, the resolution ECMWF quotes, e.g. 9 km for O1280."""
        return float(np.sqrt(4 * np.pi * EARTH_RADIUS_M**2 / self.n_cells) / 1000)

    def cell_polygons(self, bounds: tuple[float, float, float, float] | None = None) -> gpd.GeoDataFrame:
        """
        The area each cell represents, as a box of longitude and latitude in EPSG:4326 with longitudes from -180 to
        180: halfway in longitude to the cells beside it, and in latitude between the edges of its row. The edges split
        the sphere into bands whose areas are the gaussian quadrature weights of the rows, so each cell covers the share
        of the sphere the grid's own quadrature weights its value by.

        Args:
            bounds: optional (min x, min y, max x, max y) in degrees; only the cells overlapping it are made. Without
                it every cell of the world is made.

        Returns:
            gpd.GeoDataFrame: one row per cell, in cell order: cell_index (position in a GRIB message's values), x and
                y (the cell center longitude and latitude), and geometry
        """
        sine_latitude, quadrature_weight = np.polynomial.legendre.leggauss(2 * self.N)  # south to north
        row_latitude = np.degrees(np.arcsin(sine_latitude[::-1]))
        sine_row_edge = 1 - np.concatenate(([0], np.cumsum(quadrature_weight[::-1])))
        sine_row_edge[-1] = -1  # the running sum of the weights reaches 2, the south pole, only to within rounding
        row_edge = np.degrees(np.arcsin(sine_row_edge))
        row = np.repeat(np.arange(2 * self.N), self.cells_per_row)
        cells_before_row = np.concatenate(([0], np.cumsum(self.cells_per_row)[:-1]))
        cells_in_row = self.cells_per_row[row]
        x = (360 * (np.arange(row.shape[0]) - cells_before_row[row]) / cells_in_row + 180) % 360 - 180
        west, east = x - 180 / cells_in_row, x + 180 / cells_in_row
        south, north = row_edge[row + 1], row_edge[row]
        straddles = west < -180  # the cells centered on -180 degrees, which reach past it onto the east edge of the map
        cell_index = np.arange(row.shape[0])
        if bounds is not None:
            min_x, min_y, max_x, max_y = bounds
            overlaps_x = ((east > min_x) & (west < max_x)) | (straddles & (west + 360 < max_x))
            cell_index = np.flatnonzero(overlaps_x & (north > min_y) & (south < max_y))
        west, east = west[cell_index], east[cell_index]
        south, north = south[cell_index], north[cell_index]
        straddles = straddles[cell_index]
        geometry = shapely.box(np.maximum(west, -180), south, east, north)
        # each cell reaching past -180 degrees is the two parts of its box, one on either edge of the map
        east_edge_part = shapely.box(west[straddles] + 360, south[straddles], 180, north[straddles])
        geometry[straddles] = shapely.multipolygons(np.stack([geometry[straddles], east_edge_part], axis=1))
        columns = {'cell_index': cell_index, 'x': x[cell_index], 'y': row_latitude[row[cell_index]]}
        return gpd.GeoDataFrame(columns, geometry=geometry, crs=4326)


def cell_xy_from_regular_grid(
    dataset: PathInput, x_var: str = 'lon', y_var: str = 'lat'
) -> tuple[np.ndarray, np.ndarray]:
    """Get cell center x and y coordinates from a regular grid common dataset structure."""
    with xr.open_dataset(dataset) as ds:
        if x_var not in ds.variables:
            raise KeyError(f'{x_var} must be a variable in {dataset}')
        if y_var not in ds.variables:
            raise KeyError(f'{y_var} must be a variable in {dataset}')
        x = ds[x_var].to_numpy()
        y = ds[y_var].to_numpy()

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
    cell_columns = ('x_index', 'y_index', 'x', 'y')
    df = _proportions_of_catchment_areas(intersections, river_id_variable, cell_columns)

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
        river_id_variable=var_river_id,
    )
    return _order_and_save_weight_table(
        df,
        cell_columns=('x_index', 'y_index', 'x', 'y'),
        var_river_id=var_river_id,
        routing_params_path=routing_params_path,
        save_weights_path=save_weights_path,
        attributes={'grid_path': str(grid_path), 'catchments_path': str(catchments_path)},
    )


def reduced_grid_weights(
    grib_path: PathInput,
    catchments_path: PathInput,
    *,
    var_river_id: str = 'river_id',
    var_catchment_id: str = 'river_id',
    save_weights_path: PathInput | None = None,
    routing_params_path: PathInput | None = None,
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
        catchments_path: GeoParquet file of the catchment polygons, in any CRS
        var_river_id: name of the river id column in the routing params and the weight table
        var_catchment_id: name of the river id column in the catchments file
        save_weights_path: optional path to save the grid weights as a netCDF file
        routing_params_path: optional path to a routing params parquet file whose river_id column order is used to
            topologically sort the weight table rows. When omitted the row order is spatial (not topological).

    Returns:
        pd.DataFrame: the grid weights with columns [river_id, cell_index, x, y, area_sqm, area_sqm_total, proportion]
    """
    grid = ReducedGaussianGrid.from_grib(grib_path)
    catchments_gdf = gpd.read_parquet(catchments_path, columns=[var_catchment_id, 'geometry']).to_crs({'proj': 'cea'})
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
            var_river_id: catchments_gdf[var_catchment_id].to_numpy()[catchment],
            'cell_index': cells_gdf['cell_index'].to_numpy()[cell_row],
            'x': cells_gdf['x'].to_numpy()[cell_row],
            'y': cells_gdf['y'].to_numpy()[cell_row],
            'area_sqm': area_sqm,
        }
    )
    df = _proportions_of_catchment_areas(pieces_df[area_sqm > 0], var_river_id, ('cell_index', 'x', 'y'))
    return _order_and_save_weight_table(
        df,
        cell_columns=('cell_index', 'x', 'y'),
        var_river_id=var_river_id,
        routing_params_path=routing_params_path,
        save_weights_path=save_weights_path,
        attributes={'grid_path': str(grib_path), 'catchments_path': str(catchments_path)},
    )


def _proportions_of_catchment_areas(
    pieces: pd.DataFrame, river_id_variable: str, cell_columns: tuple[str, ...]
) -> pd.DataFrame:
    """Sum the area_sqm of the pieces of each catchment in each cell, then each cell's proportion of the catchment."""
    df = (
        pieces[[river_id_variable, *cell_columns, 'area_sqm']]
        .groupby([river_id_variable, *cell_columns], as_index=False)
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
    return df.astype(dict.fromkeys(('x', 'y', 'area_sqm', 'area_sqm_total', 'proportion'), np.float32))


def _order_and_save_weight_table(
    df: pd.DataFrame,
    *,
    cell_columns: tuple[str, ...],
    var_river_id: str,
    routing_params_path: PathInput | None,
    save_weights_path: PathInput | None,
    attributes: dict,
) -> pd.DataFrame:
    """Sort the weight table's rivers into the order of the routing params, then save it when a path is given."""
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
            df[[var_river_id, *cell_columns, 'area_sqm', 'proportion']]
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
