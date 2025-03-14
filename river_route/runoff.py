import os

import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr

__all__ = [
    'calc_catchment_volumes',
    'write_catchment_volumes',
]


def _cumulative_to_incremental(df) -> pd.DataFrame:
    return pd.DataFrame(
        np.vstack([df.values[0, :], np.diff(df.values, axis=0)]),
        index=df.index,
        columns=df.columns
    )


def _incremental_to_cumulative(df) -> pd.DataFrame:
    return df.cumsum()


def _get_conversion_factor(unit: str) -> int or float:
    if unit is None:
        print("No units attribute found. Assuming meters")
        return 1
    if unit in ('m', 'meters', 'kg m-2'):
        return 1
    elif unit in ('mm', 'millimeters'):
        return .001
    else:
        raise ValueError(f"Unknown units: {unit}")


def calc_catchment_volumes(
        runoff_data: str,
        weight_table: str,
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
) -> pd.DataFrame:
    """
    Calculates the catchment runoff volumes from a given runoff dataset and a directory of VPU configs.

    Args:
        runoff_data (str): a string or list of strings of paths to LSM files
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
        pd.DataFrame: a DataFrame of the catchment runoff volumes with stream IDs as columns and a datetime index
    """
    with xr.open_dataset(weight_table) as ds:
        weight_df = ds[['river_id', 'x_index', 'y_index', 'proportion', 'area_sqm_total']].to_dataframe()
    unique_idxs = weight_df[['x_index', 'y_index']].drop_duplicates().reset_index(drop=True).reset_index().astype(int)
    unique_sorted_rivers = weight_df[['river_id', 'area_sqm_total']].drop_duplicates().sort_index()

    with xr.open_mfdataset(runoff_data) as ds:
        runoff_depth_unit = runoff_depth_unit if runoff_depth_unit else ds[runoff_var].attrs.get('units', 'm')
        conversion_factor = _get_conversion_factor(runoff_depth_unit)
        df = pd.DataFrame(
            ds
            [runoff_var]
            .isel({
                x_var: xr.DataArray(unique_idxs['x_index'].values, dims="points"),
                y_var: xr.DataArray(unique_idxs['y_index'].values, dims="points")
            })
            .transpose(time_var, "points")
            .values,
            columns=unique_idxs[['x_index', 'y_index']].astype(str).apply('_'.join, axis=1),
            index=ds[time_var].to_numpy()
        )
    df = df[weight_df[['x_index', 'y_index']].astype(str).apply('_'.join, axis=1)]
    df.columns = weight_df['river_id']
    df = df * weight_df['area_sqm_total'].values * conversion_factor
    df = df.T.groupby(by=df.columns).sum().T
    df = df[unique_sorted_rivers['river_id'].values]

    # conversions and restructuring
    if cumulative:
        df = _cumulative_to_incremental(df)
    if force_positive_runoff:
        df = df.clip(lower=0)
    time_diff = np.diff(df.index)
    if not np.all(time_diff == df.index[1] - df.index[0]) and force_uniform_timesteps:
        timestep = (df.index[1] - df.index[0]).astype('timedelta64[s]').astype(int)
        log.warning(f'Time steps are not uniform, resampling to the first timestep: {timestep} seconds')
        df = (
            _incremental_to_cumulative(df)
            .resample(rule=f'{timestep}S')
            .interpolate(method='linear')
        )
        df = _cumulative_to_incremental(df)
    df = df.fillna(0)
    return df


def write_catchment_volumes(df: pd.DataFrame, output_dir: str, label: str = None) -> None:
    """
    Write the catchment runoff volumes to file in the river-route expected format.

    Args:
        df: pd.DataFrame
            The catchment volumes timeseries dataframe with stream IDs as columns and a datetime index
        output_dir: str
            The directory to write the file to
        label: str
            An optional label to include in the file name for organization purposes
    """
    # Create output inflow netcdf data
    os.makedirs(output_dir, exist_ok=True)
    start_date = df.index[0].strftime('%Y%m%d%H')
    end_date = df.index[-1].strftime('%Y%m%d%H')
    file_name = f'volumes{f"_{label}" if label else ""}_{start_date}_{end_date}.nc'
    inflow_file_path = os.path.join(output_dir, file_name)

    with nc.Dataset(inflow_file_path, "w", format="NETCDF4") as ds:
        ds.createDimension('time', df.shape[0])
        ds.createDimension('river_id', df.shape[1])

        ro_vol_var = ds.createVariable('volume', 'f4', ('time', 'river_id'), zlib=True, complevel=5)
        ro_vol_var[:] = df.to_numpy()
        ro_vol_var.long_name = 'Incremental catchment runoff volume'
        ro_vol_var.units = 'm3'

        id_var = ds.createVariable('river_id', 'i4', ('river_id',), zlib=True, complevel=5)
        id_var[:] = df.columns.astype(int)
        id_var.long_name = 'unique ID number for each river'

        timestep = (df.index[1] - df.index[0]).seconds
        time_var = ds.createVariable('time', 'i4', ('time',), zlib=True, complevel=5)
        time_var[:] = (df.index - df.index[0]).astype('timedelta64[s]').astype(int)
        time_var.long_name = 'time'
        time_var.standard_name = 'time'
        time_var.units = f'seconds since {df.index[0].strftime("%Y-%m-%d %H:%M:%S")}'
        time_var.axis = 'T'
        time_var.time_step = f'{timestep}'
    return
