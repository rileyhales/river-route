# Author: Riley Hales, PhD
# Date: 2024-08-18
# Copyright: 2024 Riley Hales, all rights reserved
# Description: Reference code for creating cached intersections of catchments and LSM grid cells
# mamba install -c conda-forge geopandas numpy pandas shapely xarray natsort

import glob
import logging
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely.geometry
import shapely.ops
import xarray as xr
from natsort import natsorted

# Get the logger
logger = logging.getLogger(__name__)
# Set up the logger formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# Create a console handler
console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)
# Add the handler to the logger
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)


def cell_polygons_from_xy_spacing(x: np.ndarray, y: np.ndarray, crs: int = 4326) -> gpd.GeoDataFrame:
    """
    Create a GeoDataFrame of Voroni polygons around the center of each cell in a 2D grid

    Args:
        y: a 1 dimensional np.ndarray of y coordinate centers
        x: a 1 dimensional np.ndarray of x coordinate centers
        crs: the coordinate reference system number of the grid to pass to as GeoDataFrame(crs=crs)

    Returns:
        gpd.GeoDataFrame: a GeoDataFrame of polygons around the center of each cell
    """
    x_grid, y_grid = np.meshgrid(x, y)
    x_grid = x_grid.flatten()
    y_grid = y_grid.flatten()

    # create cell polygon from the cell center points
    # the order of polygons in the voronoi diagram is **guaranteed not** the same as the order of the input points
    logger.info('Creating Voroni polygons')
    regions = shapely.ops.voronoi_diagram(
        shapely.geometry.MultiPoint([shapely.geometry.Point(x, y) for x, y in zip(x_grid, y_grid)])
    )

    logger.info('Adding attributes to voronoi polygons')
    voroni_gdf = gpd.GeoDataFrame(geometry=[region for region in regions.geoms], crs=crs)
    voroni_gdf['x'] = voroni_gdf.geometry.apply(lambda f: f.centroid.x).astype(float)
    voroni_gdf['y'] = voroni_gdf.geometry.apply(lambda f: f.centroid.y).astype(float)
    voroni_gdf['x_index'] = voroni_gdf['x'].apply(lambda v: np.argmin(np.abs(x - v)))
    voroni_gdf['y_index'] = voroni_gdf['y'].apply(lambda v: np.argmin(np.abs(y - v)))
    return voroni_gdf


def make_overlay_table(voroni_gdf: gpd.GeoDataFrame, catchments_gdf: gpd.GeoDataFrame,
                       catchment_id_field: str) -> pd.DataFrame:
    """
    Create a table of intersections between Voroni polygons and catchments

    Args:
        voroni_gdf: a GeoDataFrame of Voroni polygons around the center of each cell in a 2D grid
        catchments_gdf: a GeoDataFrame of catchments polygons
        catchment_id_field: the field in catchments_gdf to use as the id field

    Returns:
        pd.DataFrame: a table of intersections between Voroni polygons and catchments
    """
    logger.info('Performing overlay operation')
    intersections = gpd.overlay(voroni_gdf, catchments_gdf, how='intersection')
    logger.info('Calculating area of intersections')
    intersections['area_sqm'] = intersections.geometry.to_crs({'proj': 'cea'}).area
    # merge rows with the same catchment id, lon_index, lat_index, lon, and lat, summing the area_sqm
    return (
        intersections
        [[catchment_id_field, 'lon_index', 'lat_index', 'lon', 'lat', 'area_sqm']]
        .groupby([catchment_id_field, 'lon_index', 'lat_index', 'lon', 'lat'])
        .agg({'area_sqm': 'sum'})
        .reset_index()
        .sort_values([catchment_id_field, 'area_sqm'], ascending=[True, False])
    )


def make_voronoi_polygons(ds: str, save_dir: str, x_var: str = 'lon', y_var: str = 'lat', crs: int = 4326) -> None:
    """
    Create a GeoDataFrame of Voroni polygons around the center of each cell in a 2D grid given in a sample land surface
    output dataset

    Args:
        ds: Path to the dataset
        save_dir: Path to save the Voroni polygons
        x_var: The x variable name in the dataset
        y_var: The y variable name in the dataset
        crs: The coordinate reference system number of the grid to pass to as GeoDataFrame(crs=crs)

    Returns:
        None
    """
    with xr.open_dataset(ds) as ds:
        x = ds[x_var].values
        y = ds[y_var].values

    x_origin = x[0]
    y_origin = y[0]
    x_spacing = x[1] - x[0]
    y_spacing = y[1] - y[0]
    file_name = f'voronipolygons_x0={x_origin}_y0={y_origin}_dx={x_spacing}_dy={y_spacing}.geoparquet'

    voroni_gdf = cell_polygons_from_xy_spacing(x, y, crs)
    voroni_gdf.to_parquet(os.path.join(save_dir, file_name))
    return


if __name__ == '__main__':
    grids = natsorted(glob.glob('./sample_datasets/*.nc'))

    for grid in grids:
        make_voronoi_polygons(grid, './voroni_polygons')

    grids = natsorted(glob.glob('./voroni_polygons/*.geoparquet'))
    catchments = natsorted(glob.glob('./catchments/*.geoparquet'))
    id_field = 'LINKNO'

    for catchment in catchments:
        catchment_label = ''
        catchment_save_dir = os.path.join('./inputs', catchment_label)
        c_gdf = gpd.read_parquet(catchment)
        c_bbox = c_gdf.total_bounds
        for grid in grids:
            # determine the table name for the given catchment and grid combination
            table_name = str(
                os.path.basename(grid)
                .replace('voronipolygons', f'intersections_{catchment_label}')
                .replace('.geoparquet', '.parquet')
            )
            # read the voroni polygons and filter for slightly faster and more memory efficient processing
            v_gdf = gpd.read_parquet(grid)
            v_gdf = v_gdf.cx[c_bbox[0]:c_bbox[2], c_bbox[1]:c_bbox[3]]
            # make and save the overlay table
            (
                make_overlay_table(v_gdf, c_gdf, id_field)
                .to_parquet(os.path.join(catchment_save_dir, table_name), index=False)
            )
