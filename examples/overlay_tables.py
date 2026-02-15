import geopandas as gpd

from river_route.runoff import get_cell_xy_from_regular_grid, voroni_diagram_from_cell_xy, compute_voroni_catchment_intersects

if __name__ == '__main__':
    sample_grid = './sample_grid.nc'
    catchments = './catchments/*.parquet'
    voroni_polygons_save_path = './voroni_polygons_sample_grid.parquet'
    grid_weights_save_path = './grid_weights_sample_grid.parquet'

    x, y = get_cell_xy_from_regular_grid(sample_grid, x_var='lon', y_var='lat')
    voroni_gdf = voroni_diagram_from_cell_xy(x, y, crs=4326)
    voroni_gdf.to_parquet(voroni_polygons_save_path)

    catchments_gdf = gpd.read_parquet(catchments)
    catchments_bounds = catchments_gdf.total_bounds
    voroni_gdf = voroni_gdf.cx[catchments_bounds[0]:catchments_bounds[2], catchments_bounds[1]:catchments_bounds[3]]
    compute_voroni_catchment_intersects(voroni_gdf, catchments_gdf, grid_weights_save_path)
