import numpy as np
import pandas as pd
import xarray as xr


def merge_routing_params(routing_params_path: str, connectivity_path: str, output_path: str) -> None:
    """
    routing_params_path: Path to the v1 routing parameters parquet (river_id, k, x).
    connectivity_path: Path to the v1 connectivity parquet (river_id, ds_river_id).
    output_path: Where to write the merged v2 parquet (river_id, downstream_river_id, k, x).
    """
    params = pd.read_parquet(routing_params_path)
    connectivity = pd.read_parquet(connectivity_path).rename(columns={'ds_river_id': 'downstream_river_id'})
    merged = params.merge(connectivity, on='river_id', how='inner')
    required = ['river_id', 'downstream_river_id', 'k', 'x']
    missing = set(required) - set(merged.columns)
    if missing:
        raise ValueError(f'Merged file is missing required columns: {missing}')
    merged[list(required)].to_parquet(output_path, index=False)


def convert_grid_weights(csv_path: str, output_path: str) -> None:
    """
    csv_path: Path to the v1 grid weights CSV.
    output_path: Where to write the v2 netCDF file.
    """
    df = (
        pd
        .read_csv(csv_path)
        .rename(columns={
            'lon_index': 'x_index',
            'lat_index': 'y_index',
            'lon': 'x',
            'lat': 'y', }
        )
    )

    if 'proportion' not in df.columns:
        if 'area_sqm' not in df.columns:
            raise ValueError(
                'Cannot compute proportion: grid weights file has no "area_sqm" column. '
                'Add an "area_sqm" or "proportion" column manually.'
            )
        total_area = (
            df[['river_id', 'area_sqm']]
            .groupby('river_id')
            .sum()
            .rename(columns={'area_sqm': 'area_sqm_total'})
        )
        df = df.merge(total_area, left_on='river_id', right_index=True, how='left')
        df['proportion'] = df['area_sqm'] / df['area_sqm_total']
        df.drop(columns='area_sqm_total', inplace=True)

    ds = xr.Dataset(
        data_vars={
            'river_id': ('index', df['river_id'].to_numpy(dtype=np.int64)),
            'x_index': ('index', df['x_index'].to_numpy(dtype=np.int64)),
            'y_index': ('index', df['y_index'].to_numpy(dtype=np.int64)),
            'x': ('index', df['x'].to_numpy(dtype=np.float64)),
            'y': ('index', df['y'].to_numpy(dtype=np.float64)),
            'area_sqm': ('index', df['area_sqm'].to_numpy(dtype=np.float64)),
            'proportion': ('index', df['proportion'].to_numpy(dtype=np.float64)),
        },
        attrs={
            'description': 'proportions of runoff cells that intersect catchments for use with river-route',
        },
    )
    ds.to_netcdf(output_path)
