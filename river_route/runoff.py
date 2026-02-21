import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely.geometry
import shapely.ops
import xarray as xr

from .__metadata__ import __version__
from .typing import PathInput

__all__ = [
    # for making grid_weights
    'cell_xy_from_regular_grid',
    'voronoi_diagram_from_regular_grid_cell_xy',
    'compute_voronoi_catchment_intersects',
    'grid_weights',
    # for making catchment volumes from grid weights and depth grids
    'depth_to_volume',
]

logger = logging.getLogger(__name__)


def cell_xy_from_regular_grid(dataset: PathInput, x_var: str = 'lon', y_var: str = 'lat', ) -> tuple[
    np.ndarray, np.ndarray]:
    """
    Get cell center x and y coordinates from a regular grid in expected, common dataset structure.
    """
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


def voronoi_diagram_from_regular_grid_cell_xy(x: np.ndarray, y: np.ndarray, crs: int = 4326) -> gpd.GeoDataFrame:
    """
    Create a GeoDataFrame of Voronoi polygons around the center of each cell in a grid.
    The voronoi
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
                                         save_path: PathInput | None = None, save_attributes: dict | None = None,
                                         river_id_variable: str = 'river_id') -> pd.DataFrame:
    """
    Create a table of intersections between Voronoi polygons and catchments.
    """
    if river_id_variable not in catchments_gdf.columns:
        raise KeyError('catchments_gdf must contain a river_id column')
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
        save_attributes = {
            'description': 'proportions of runoff cells that intersect river catchments to be used with river-route',
            'voronoi_gdf_crs': voronoi_gdf.crs.to_string(),
            'catchments_gdf_crs': catchments_gdf.crs.to_string(),
            'river_route_version': __version__,
            **(save_attributes or {}),
        }
        ds = xr.Dataset({
            'river_id': ('index', df[river_id_variable].to_numpy(dtype=np.int64, copy=False)),
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


def grid_weights(grid_path: PathInput, catchments_path: PathInput, *,
                 x_var: str = 'lon', y_var: str = 'lat', river_id_var: str = 'river_id',
                 save_voronoi_path: PathInput | None = None,
                 save_weights_path: PathInput | None = None,
                 routing_params_path: PathInput | None = None) -> pd.DataFrame:
    """
    Compute the grid weights for a given grid and catchments.

    Args:
        grid_path: path to a NetCDF file containing the grid information (must include 'lon' and 'lat' variables)
        catchments_path: path to a GeoParquet file containing the catchment geometries (must include 'river_id' column)
        x_var: x-coordinate variable name in the grid file
        y_var: y-coordinate variable name in the grid file
        river_id_var: variable name for river ID in the catchments file
        save_voronoi_path: optional path to save the Voronoi polygons as a GeoParquet file
        save_weights_path: optional path to save the grid weights as a NetCDF file
        routing_params_path: optional path to a routing params parquet file whose river_id column order is used to
            topologically sort the weight table rows. When omitted the row order is spatial (not topological).

    Returns:
        pd.DataFrame: a DataFrame containing the grid weights with columns
            ['river_id', 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion']
    """
    x, y = cell_xy_from_regular_grid(grid_path, x_var=x_var, y_var=y_var)
    voronoi_gdf = voronoi_diagram_from_regular_grid_cell_xy(x, y, crs=4326)
    if save_voronoi_path:
        voronoi_gdf.to_parquet(save_voronoi_path)
    catchments_gdf = gpd.read_parquet(catchments_path)
    df = compute_voronoi_catchment_intersects(
        voronoi_gdf, catchments_gdf, save_path=None,
        save_attributes={'grid_path': str(grid_path), 'catchments_path': str(catchments_path)},
        river_id_variable=river_id_var
    )

    if routing_params_path is not None:
        ordered_ids = pd.read_parquet(routing_params_path)[river_id_var].to_numpy()
        id_to_order = {int(rid): i for i, rid in enumerate(ordered_ids)}
        df = (
            df
            .assign(_sort_key=df[river_id_var].map(id_to_order))
            .sort_values(['_sort_key', 'area_sqm'], ascending=[True, False])
            .drop(columns='_sort_key')
            .reset_index(drop=True)
        )
    else:
        logger.warning('routing_params_path not provided; weight table row order may not match routing network order')

    if save_weights_path:
        save_attributes = {
            'description': 'proportions of runoff cells that intersect river catchments to be used with river-route',
            'grid_path': str(grid_path),
            'catchments_path': str(catchments_path),
            'river_route_version': __version__,
        }
        ds = xr.Dataset({
            river_id_var: ('index', df[river_id_var].to_numpy(dtype=np.int64, copy=False)),
            'x_index': ('index', df['x_index'].to_numpy(dtype=np.int64, copy=False)),
            'y_index': ('index', df['y_index'].to_numpy(dtype=np.int64, copy=False)),
            'x': ('index', df['x'].to_numpy(dtype=np.float64, copy=False)),
            'y': ('index', df['y'].to_numpy(dtype=np.float64, copy=False)),
            'area_sqm': ('index', df['area_sqm'].to_numpy(dtype=np.float64, copy=False)),
            'proportion': ('index', df['proportion'].to_numpy(dtype=np.float64, copy=False)),
        })
        ds.attrs.update(save_attributes)
        ds.to_netcdf(save_weights_path)

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
        print("No units attribute found. Assuming meters")
        return 1
    if unit in ('m', 'meters', 'kg m-2'):
        return 1
    elif unit in ('mm', 'millimeters'):
        return .001
    else:
        raise ValueError(f"Unknown units: {unit}")


def depth_to_volume(
        runoff_data: PathInput | list[PathInput],
        weight_table: PathInput,
        *,
        runoff_var: str = 'ro',
        x_var: str = 'lon',
        y_var: str = 'lat',
        time_var: str = 'time',
        river_id_var: str = 'river_id',
        runoff_depth_unit: str = None,
        cumulative: bool = False,
        force_positive_runoff: bool = False,
        force_uniform_timesteps: bool = True,
) -> xr.Dataset:
    """
    Calculates the catchment runoff volumes from a given runoff dataset and a directory of VPU configs.

    Args:
        runoff_data (str | list[str]): path(s) to runoff files
        weight_table (str): a string path to the weight table
        runoff_var (str): the name of the runoff variable in the LSM files
        x_var (str): the name of the x variable in the LSM files
        y_var (str): the name of the y variable in the LSM files
        time_var (str): the name of the time variable in the LSM files
        river_id_var (str): the name of the river ID variable in the parameters file
        runoff_depth_unit (str): the unit of the depths in the runoff data
        cumulative (bool): whether the runoff data is cumulative or incremental
        force_positive_runoff (bool): whether to force all runoff values to be >= 0
        force_uniform_timesteps (bool): whether to force all timesteps to be uniform

    Returns:
        xr.Dataset: a Dataset of catchment runoff volumes with dimensions ``time`` and ``river_id``
    """
    with xr.open_dataset(weight_table) as ds:
        weight_df = ds[[river_id_var, 'x_index', 'y_index', 'area_sqm']].to_dataframe()
    unique_indexes = (
        weight_df
        [['x_index', 'y_index']]
        .drop_duplicates()
        .reset_index(drop=True)
        .reset_index()
        .astype(int)
    )
    unique_sorted_rivers = weight_df[[river_id_var, ]].drop_duplicates().sort_index()  # index -> already topo sorted

    with xr.open_mfdataset(runoff_data) as ds:
        runoff_depth_unit = runoff_depth_unit if runoff_depth_unit else ds[runoff_var].attrs.get('units', 'm')
        conversion_factor = _get_conversion_factor(runoff_depth_unit)
        df = pd.DataFrame(
            ds
            [runoff_var]
            .isel({
                x_var: xr.DataArray(unique_indexes['x_index'].values, dims="points"),
                y_var: xr.DataArray(unique_indexes['y_index'].values, dims="points")
            })
            .transpose(time_var, "points")
            .values,
            columns=unique_indexes[['x_index', 'y_index']].astype(str).apply('_'.join, axis=1),
            index=ds[time_var].to_numpy()
        )
    point_labels = weight_df[['x_index', 'y_index']].astype(str).apply('_'.join, axis=1)
    df = (
        df
        .loc[:, point_labels]
        .set_axis(weight_df[river_id_var], axis=1)
        .mul(weight_df['area_sqm'].values * conversion_factor, axis=1)
        .T.groupby(level=0).sum().T
        .loc[:, unique_sorted_rivers[river_id_var].values]
    )

    # conversions and restructuring
    if cumulative:
        df = _cumulative_to_incremental(df)
    if force_positive_runoff:
        df = df.clip(lower=0)
    time_diff = np.diff(df.index)
    if not np.all(time_diff == df.index[1] - df.index[0]) and force_uniform_timesteps:
        timestep = (df.index[1] - df.index[0]).astype('timedelta64[s]').astype(int)
        logger.warning(f'Time steps are not uniform, resampling to the first timestep: {timestep} seconds')
        df = (
            _incremental_to_cumulative(df)
            .resample(rule=f'{timestep}S')
            .interpolate(method='linear')
        )
        df = _cumulative_to_incremental(df)
    df = df.fillna(0)

    # create xr.Dataset to return
    start_date = df.index[0].strftime('%Y%m%d%H')
    end_date = df.index[-1].strftime('%Y%m%d%H')
    suggested_file_name = f'volumes_{start_date}_{end_date}.nc'
    timestep = int((df.index[1] - df.index[0]).total_seconds()) if df.shape[0] > 1 else 0
    return xr.Dataset(
        {
            'volume': xr.DataArray(
                df.to_numpy(dtype=np.float32, copy=False),
                dims=('time', 'river_id'),
                attrs={
                    'long_name': 'Incremental catchment runoff volume',
                    'units': 'm3',
                },
            ),
        },
        coords={
            'river_id': xr.DataArray(
                df.columns.to_numpy(dtype=np.int64, copy=False),
                dims=('river_id',),
                attrs={
                    'long_name': 'unique ID number for each river',
                },
            ),
            'time': xr.DataArray(
                (df.index - df.index[0]).astype('timedelta64[s]').astype(np.int64),
                dims=('time',),
                attrs={
                    'long_name': 'time',
                    'standard_name': 'time',
                    'units': f'seconds since {df.index[0].strftime("%Y-%m-%d %H:%M:%S")}',
                    'axis': 'T',
                    'time_step': f'{timestep}',
                },
            ),
        },
        attrs={
            'title': 'Catchment runoff volumes',
            'description': 'Incremental catchment runoff volumes in m3 for each river',
            'source': f'river-route v{__version__}',
            'history': f'Created on {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}',
            'suggested_file_name': suggested_file_name,
        },
    )
