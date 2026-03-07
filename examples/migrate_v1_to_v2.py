"""
example code to convert river-route 1 files to river-route 2 format.
"""

import numpy as np
import pandas as pd
import xarray as xr

# Legacy column names that may appear in v1 files
_ROUTING_COLUMN_RENAMES = {
    'ds_river_id': 'downstream_river_id',
}

_GRID_WEIGHTS_COLUMN_RENAMES = {
    'lon_index': 'x_index',
    'lat_index': 'y_index',
    'lon': 'x',
    'lat': 'y',
}


def merge_routing_params(routing_params_path: str, connectivity_path: str, output_path: str) -> None:
    """
    Merge separate routing-params and connectivity parquets into one v2 file.

    Args:
        routing_params_path: Path to the v1 routing parameters parquet (river_id, k, x).
        connectivity_path: Path to the v1 connectivity parquet (river_id, ds_river_id / downstream_river_id).
        output_path: Where to write the merged v2 parquet.
    """
    params = pd.read_parquet(routing_params_path)
    connectivity = pd.read_parquet(connectivity_path)

    params.rename(columns=_ROUTING_COLUMN_RENAMES, inplace=True)
    connectivity.rename(columns=_ROUTING_COLUMN_RENAMES, inplace=True)

    merged = params.merge(connectivity, on='river_id', how='inner')

    required = {'river_id', 'downstream_river_id', 'k', 'x'}
    missing = required - set(merged.columns)
    if missing:
        raise ValueError(f'Merged file is missing required columns: {missing}')

    merged.to_parquet(output_path, index=False)
    print(f'Wrote merged routing params ({len(merged)} rows) to {output_path}')


def convert_grid_weights(csv_path: str, output_path: str) -> None:
    """Convert a v1 CSV grid-weights file to a v2 NetCDF file.

    Args:
        csv_path: Path to the v1 grid weights CSV.
        output_path: Where to write the v2 NetCDF file.
    """
    df = pd.read_csv(csv_path)
    df.rename(columns=_GRID_WEIGHTS_COLUMN_RENAMES, inplace=True)

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
    print(f'Wrote grid weights ({len(df)} rows) to {output_path}')


def convert_state_file(state_path: str, output_path: str) -> None:
    """Convert a v1 two-column state file (Q, R) to a v2 single-column file (Q).

    Args:
        state_path: Path to the v1 state parquet.
        output_path: Where to write the v2 state parquet.
    """
    df = pd.read_parquet(state_path)

    if 'Q' not in df.columns:
        raise ValueError(f'State file has no "Q" column. Found columns: {list(df.columns)}')

    df = df[['Q']]
    df.to_parquet(output_path, index=False)
    print(f'Wrote state file ({len(df)} rows) to {output_path}')
