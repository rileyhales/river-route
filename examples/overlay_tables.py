# Author: Riley Hales, PhD
# Date: 2024-08-18
# Copyright: 2024 Riley Hales, all rights reserved
# Description: Reference code for creating cached intersections of catchments and LSM grid cells
# mamba install -c conda-forge geopandas numpy pandas shapely xarray natsort

import glob
import os

import geopandas as gpd
from natsort import natsorted

from river_route.grid_weights import make_voronoi_diagram, overlay_table


if __name__ == '__main__':
    grids = natsorted(glob.glob('./sample_datasets/*.nc'))

    for grid in grids:
        make_voronoi_diagram(grid, './voroni_polygons')

    grids = natsorted(glob.glob('./voroni_polygons/*.geoparquet'))
    catchments = natsorted(glob.glob('./catchments/*.geoparquet'))

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
                overlay_table(v_gdf, c_gdf)
                .to_parquet(os.path.join(catchment_save_dir, table_name), index=False)
            )
