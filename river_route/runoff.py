import logging
from concurrent.futures import Executor
from dataclasses import dataclass
from typing import Self

import geopandas as gpd
import numba
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
    'GridWeights',
    'read_grid_runoff',
    'aggregate_grid_runoff',
    'runoff_to_vlateral',
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
    # vlateral array aggregated from it
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

    # Map x_index back to original grid indices (for use by grid_to_vlateral)
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


def _cumulative_to_incremental(df) -> pd.DataFrame:
    return pd.DataFrame(np.vstack([df.values[0, :], np.diff(df.values, axis=0)]), index=df.index, columns=df.columns)


def _incremental_to_cumulative(df) -> pd.DataFrame:
    return df.cumsum()


def _get_conversion_factor(unit: str) -> int | float:
    if unit is None:
        logger.warning('No units attribute found. Assuming meters')
        return 1
    if unit in ('m', 'meters', 'kg m-2'):
        return 1
    elif unit in ('mm', 'millimeters'):
        return 0.001
    else:
        raise ValueError(f'Unknown units: {unit}')


@dataclass(frozen=True)
class GridWeights:
    """
    A grid weight table prepared for aggregating runoff onto rivers. Reading and indexing the table is the same work
    for every runoff file, so it is built once and reused for all of them.

    Attributes:
        river_ids: (n_rivers,) river ids in weight table order, which follows the routing params order
        x_index: (n_cells,) grid x index of each unique cell any catchment touches
        y_index: (n_cells,) grid y index of each unique cell
        indptr: (n_rivers + 1,) sparse row pointers, the weights of river r are indptr[r]:indptr[r + 1]
        cell: (n_weights,) position in x_index and y_index of the cell each weight applies to
        proportion: (n_weights,) share of the river's catchment area inside the cell
        catchment_area: (n_rivers,) total catchment area of each river in m²
    """

    river_ids: np.ndarray
    x_index: np.ndarray
    y_index: np.ndarray
    indptr: np.ndarray
    cell: np.ndarray
    proportion: np.ndarray
    catchment_area: np.ndarray

    @classmethod
    def from_file(cls, grid_weights_file: PathInput, var_river_id: str = 'river_id') -> Self:
        """Read a weight table netCDF produced by ``grid_weights()``."""
        with xr.open_dataset(grid_weights_file) as ds:
            weight_df = ds[[var_river_id, 'x_index', 'y_index', 'proportion', 'area_sqm']].to_dataframe()
        unique_indexes = (
            weight_df[['x_index', 'y_index']].drop_duplicates().reset_index(drop=True).reset_index().astype(int)
        )
        # index already topo sorted
        river_ids = weight_df[[var_river_id]].drop_duplicates().sort_index()[var_river_id].to_numpy()

        cells = weight_df[['x_index', 'y_index']].merge(unique_indexes, on=['x_index', 'y_index'], how='left')
        point_idx = cells['index'].values
        river_id_to_row = pd.Series(np.arange(len(river_ids)), index=river_ids)
        river_idx = river_id_to_row.loc[weight_df[var_river_id].values].values
        matrix = scipy.sparse.csr_matrix(
            (weight_df['proportion'].values, (river_idx, point_idx)), shape=(len(river_ids), len(unique_indexes))
        )
        catchment_area = weight_df.groupby(var_river_id)['area_sqm'].sum().reindex(river_ids).to_numpy()
        return cls(
            river_ids=river_ids,
            x_index=unique_indexes['x_index'].to_numpy(),
            y_index=unique_indexes['y_index'].to_numpy(),
            indptr=matrix.indptr,
            cell=matrix.indices,
            proportion=matrix.data,
            catchment_area=catchment_area,
        )


def read_grid_runoff(
    runoff_data: PathInput | list[PathInput],
    weights: GridWeights,
    *,
    var_runoff: str = 'ro',
    var_x: str = 'lon',
    var_y: str = 'lat',
    var_t: str = 'time',
    runoff_depth_unit: str | None = None,
) -> tuple[np.ndarray, np.ndarray, int | float]:
    """
    Read runoff at the grid cells a weight table touches.

    Args:
        runoff_data: path(s) to runoff files
        weights: the prepared grid weight table
        var_runoff: runoff variable name in the LSM files
        var_x: x-coordinate variable name in the LSM files
        var_y: y-coordinate variable name in the LSM files
        var_t: time variable name in the LSM files
        runoff_depth_unit: unit of the depth values; checked for in file attributes, defaulting to meters

    Returns:
        tuple: (runoff as a C-order (time, n_cells) array, time values, factor converting the depth unit to meters)
    """
    with xr.open_mfdataset(runoff_data, chunks=None) as ds:
        runoff_depth_unit = runoff_depth_unit or ds[var_runoff].attrs.get('units', 'm')
        conversion_factor = _get_conversion_factor(runoff_depth_unit)
        runoff = (
            ds[var_runoff]
            .isel(
                {
                    var_x: xr.DataArray(weights.x_index, dims='points'),
                    var_y: xr.DataArray(weights.y_index, dims='points'),
                }
            )
            .transpose(var_t, 'points')
            .values
        )
        time_index = ds[var_t].to_numpy()
    return np.ascontiguousarray(runoff), time_index, conversion_factor


# rivers aggregated together in one scratch block. Small enough that the block's time series stay in cache: of 8 to
# 4096 on a month of Amazon runoff, 32 and 64 were fastest and within 2% of each other.
_RIVERS_PER_BLOCK = 64


# no fastmath: it lets the compiler assume values are never NaN, which would remove the NaN replacement
@numba.njit(cache=True, nogil=True)
def _aggregate_to_rivers(
    runoff_by_cell,
    indptr,
    cell,
    weight,
    scale,
    zero,
    cumulative,
    force_positive,
    replace_nan,
    r_start,
    r_stop,
    scratch,
    out,
):
    """
    Area weighted sum of each river's cells, with every per-value conversion applied in the same pass:
    de-accumulation, clipping, NaN replacement, and the per-river scale (catchment area for volumes).

    Each river's whole time series is summed from contiguous (cell, time) series, which the compiler vectorizes,
    into a scratch block of rivers small enough to stay in cache. Each finished block is then written into the
    C-order (time, river) output one short contiguous run per row. ``scale`` is empty to skip scaling.

    Only rivers ``r_start`` to ``r_stop`` are aggregated. Disjoint ranges can run concurrently because each writes
    only its own output columns, provided each is given its own ``scratch``.
    """
    n_steps = runoff_by_cell.shape[1]
    block = scratch.shape[0]
    for r0 in range(r_start, r_stop, block):
        r1 = min(r0 + block, r_stop)
        for r in range(r0, r1):
            row = scratch[r - r0]
            for t in range(n_steps):
                row[t] = zero
            for k in range(indptr[r], indptr[r + 1]):
                w = weight[k]
                series = runoff_by_cell[cell[k]]
                for t in range(n_steps):
                    row[t] += w * series[t]
            if cumulative:
                # descending, so each step subtracts the previous step's cumulative total before it is replaced
                for t in range(n_steps - 1, 0, -1):
                    row[t] -= row[t - 1]
            if force_positive:
                for t in range(n_steps):
                    if row[t] < zero:
                        row[t] = zero
            if replace_nan:
                for t in range(n_steps):
                    if row[t] != row[t]:
                        row[t] = zero
            if scale.shape[0]:
                s = scale[r]
                for t in range(n_steps):
                    row[t] *= s
        for t in range(n_steps):
            destination = out[t]
            for b in range(r1 - r0):
                destination[r0 + b] = scratch[b, t]
    return


def _river_ranges(indptr: np.ndarray, n_ranges: int) -> np.ndarray:
    """
    Split the rivers into at most ``n_ranges`` contiguous ranges of similar aggregation work.

    A river costs one pass over its time series per weight plus one for its conversions and output, so ranges are
    balanced on weights plus rivers rather than on river count alone.

    Returns:
        np.ndarray: range boundaries ``[0, ..., n_rivers]``; range i is rivers ``bounds[i]`` to ``bounds[i + 1]``
    """
    n_rivers = indptr.shape[0] - 1
    work = indptr + np.arange(n_rivers + 1)
    targets = np.linspace(0, work[-1], max(1, n_ranges) + 1)
    return np.unique(np.concatenate(([0], np.searchsorted(work, targets), [n_rivers])))


def _aggregate_over_ranges(
    runoff_by_cell: np.ndarray,
    weights: GridWeights,
    weight: np.ndarray,
    scale: np.ndarray,
    dtype: np.dtype,
    cumulative: bool,
    force_positive: bool,
    replace_nan: bool,
    out: np.ndarray,
    pool: Executor | None,
    threads: int,
) -> None:
    """
    Run the aggregation kernel over every river: as one range covering the network when single threaded, or as
    concurrent ranges of similar work when a pool is given. The kernel is the same either way, so there is no separate
    serial code path to keep in sync.
    """
    n_steps = runoff_by_cell.shape[1]
    bounds = _river_ranges(weights.indptr, threads if pool is not None else 1)

    def run(i: int) -> None:
        r_start, r_stop = int(bounds[i]), int(bounds[i + 1])
        scratch = np.empty((min(_RIVERS_PER_BLOCK, r_stop - r_start), n_steps), dtype=dtype)
        _aggregate_to_rivers(
            runoff_by_cell,
            weights.indptr,
            weights.cell,
            weight,
            scale,
            dtype.type(0),
            cumulative,
            force_positive,
            replace_nan,
            r_start,
            r_stop,
            scratch,
            out,
        )

    n_ranges = bounds.shape[0] - 1
    if pool is None or n_ranges < 2:
        for i in range(n_ranges):
            run(i)
    else:
        list(pool.map(run, range(n_ranges)))  # list() so a worker exception propagates
    return


def aggregate_grid_runoff(
    runoff: np.ndarray,
    time_index: np.ndarray,
    weights: GridWeights,
    conversion_factor: int | float = 1,
    *,
    cumulative: bool = False,
    force_positive_runoff: bool = False,
    force_uniform_timesteps: bool = True,
    as_volumes: bool = False,
    out: np.ndarray | None = None,
    pool: Executor | None = None,
    threads: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Aggregate gridded runoff depths onto rivers as area weighted depths or volumes in a single pass, written straight
    into a C-order (time, river) array that the routing kernels can read without a copy.

    Args:
        runoff: (time, n_cells) runoff depths at the weight table's cells, as returned by ``read_grid_runoff``
        time_index: (time,) datetime64 values of the runoff steps
        weights: the prepared grid weight table
        conversion_factor: multiplier converting the runoff depth unit to meters
        cumulative: whether the runoff data is cumulative; converted to incremental if True
        force_positive_runoff: clip negative runoff values to zero
        force_uniform_timesteps: resample to a uniform timestep if the input is irregular
        as_volumes: if True, return volumes (m³) instead of depths (m)
        out: optional C-order array with at least as many rows as ``runoff`` and one column per river, reused as the
            output buffer so repeated calls allocate nothing. Not used when irregular timesteps are resampled.
        pool: optional thread pool. The rivers are split into ``threads`` ranges of similar work that the same
            kernel aggregates concurrently. Without a pool one range covers every river.
        threads: number of river ranges to aggregate concurrently when ``pool`` is given

    Returns:
        tuple: (vlateral as a C-order (time, n_rivers) array, its time values). When ``out`` is used the array is a
            view of its first rows.
    """
    n_steps, n_rivers = runoff.shape[0], weights.river_ids.shape[0]
    runoff_by_cell = np.ascontiguousarray(runoff.T)  # each cell's time series contiguous for the vectorized sums
    weight = weights.proportion * conversion_factor if conversion_factor != 1 else weights.proportion
    dtype = np.result_type(runoff.dtype, weight.dtype)
    no_scale = weights.catchment_area[:0]

    time_diff = np.diff(time_index)
    if force_uniform_timesteps and time_diff.size and not np.all(time_diff == time_diff[0]):
        vlateral = np.empty((n_steps, n_rivers), dtype=dtype)
        _aggregate_over_ranges(
            runoff_by_cell,
            weights,
            weight,
            no_scale,
            dtype,
            cumulative,
            force_positive_runoff,
            False,
            vlateral,
            pool,
            threads,
        )
        timestep = int((time_index[1] - time_index[0]) / np.timedelta64(1, 's'))
        logger.warning(f'Time steps are not uniform, resampling to the first timestep: {timestep} seconds')
        df = pd.DataFrame(vlateral, index=time_index, columns=weights.river_ids)
        df = _incremental_to_cumulative(df).resample(rule=f'{timestep}s').interpolate(method='linear')
        df = _cumulative_to_incremental(df)
        time_index = df.index.values
        # pandas 3 returns a read-only array from to_numpy, and the conversions below write in place
        vlateral = df.to_numpy(dtype=np.float32, copy=True)
        del df
        vlateral[np.isnan(vlateral)] = 0.0
        if as_volumes:
            vlateral *= weights.catchment_area[np.newaxis, :]
        return vlateral, time_index

    if out is None:
        out = np.empty((n_steps, n_rivers), dtype=dtype)
    elif out.ndim != 2 or out.shape[0] < n_steps or out.shape[1] != n_rivers or not out.flags.c_contiguous:
        raise ValueError(f'out must be a C-order array of at least ({n_steps}, {n_rivers}), got shape {out.shape}')
    vlateral = out[:n_steps]
    _aggregate_over_ranges(
        runoff_by_cell,
        weights,
        weight,
        weights.catchment_area if as_volumes else no_scale,
        dtype,
        cumulative,
        force_positive_runoff,
        True,
        vlateral,
        pool,
        threads,
    )
    return vlateral, time_index


def runoff_to_vlateral(
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
    Aggregates gridded runoff depths to catchment level vlateral as depths or volumes.
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
        xr.Dataset: vlateral with dimensions ``time`` and ``river_id``.
            Contains a single variable ``vlateral`` in metres or m³.
    """
    weights = GridWeights.from_file(grid_weights_file, var_river_id=var_river_id)
    runoff, time_index, conversion_factor = read_grid_runoff(
        runoff_data,
        weights,
        var_runoff=var_runoff,
        var_x=var_x,
        var_y=var_y,
        var_t=var_t,
        runoff_depth_unit=runoff_depth_unit,
    )
    vlateral, time_index = aggregate_grid_runoff(
        runoff,
        time_index,
        weights,
        conversion_factor,
        cumulative=cumulative,
        force_positive_runoff=force_positive_runoff,
        force_uniform_timesteps=force_uniform_timesteps,
        as_volumes=as_volumes,
    )
    del runoff

    units = 'm3' if as_volumes else 'm'
    long_name = 'Incremental vlateral volumes' if as_volumes else 'Incremental vlateral depths'
    start_date = pd.Timestamp(time_index[0]).strftime('%Y%m%d%H')
    end_date = pd.Timestamp(time_index[-1]).strftime('%Y%m%d%H')
    timestep = int((time_index[1] - time_index[0]) / np.timedelta64(1, 's')) if len(time_index) > 1 else 0
    return xr.Dataset(
        {'vlateral': xr.DataArray(vlateral, dims=('time', 'river_id'), attrs={'long_name': long_name, 'units': units})},
        coords={
            'river_id': xr.DataArray(
                weights.river_ids.astype(np.int64, copy=False),
                dims=('river_id',),
                attrs={'long_name': 'unique ID number for each river'},
            ),
            'time': xr.DataArray(
                time_index,
                dims=('time',),
                attrs={'long_name': 'time', 'standard_name': 'time', 'axis': 'T', 'time_step': f'{timestep}'},
            ),
        },
        attrs={
            'title': f'Incremental vlateral {long_name.split()[-1]}',
            'description': f'Incremental vlateral ({units}) for each river',
            'source': f'river-route v{__version__}',
            'history': f'Created on {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}',
            'suggested_file_name': f'vlateral_{start_date}_{end_date}.nc',
        },
    )
