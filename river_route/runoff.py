import glob
import logging
import os
import re

import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr

__all__ = [
    'calc_catchment_volumes',
    'write_runoff_volumes',
]


def _cumulative_to_incremental(df) -> pd.DataFrame:
    return pd.DataFrame(
        np.vstack([df.values[0, :], np.diff(df.values, axis=0)]),
        index=df.index,
        columns=df.columns
    )


def _incremental_to_cumulative(df) -> pd.DataFrame:
    return df.cumsum()


def calc_catchment_volumes(
        runoff_data: str,
        weight_table: str,
        runoff_var: str = 'ro',
        x_var: str = 'lon',
        y_var: str = 'lat',
        river_id_var: str = 'river_id',
        time_var: str = 'time',
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
        river_id_var (str): the name of the river ID variable in the parameters file
        time_var (str): the name of the time variable in the LSM files
        cumulative (bool): whether the runoff data is cumulative or incremental
        force_positive_runoff (bool): whether to force all runoff values to be >= 0
        force_uniform_timesteps (bool): whether to force all timesteps to be uniform

    Returns:
        pd.DataFrame: a DataFrame of the catchment runoff volumes with stream IDs as columns and a datetime index
    """
    assert os.path.exists(runoff_data), f"Runoff data not found: {runoff_data}"

    weight_df = xr.open_dataset(weight_table)[['river_id', 'lat_index', 'lon_index', 'weight']].to_dataframe()
    sorted_rivid_array = xr.open_dataset(weight_table)['sorted_river_ids'].values.flatten()

    with xr.open_mfdataset(runoff_data) as ds:
        # Check units and conversion factors
        conversion_factor = 1
        units = ds.attrs.get('units', False)
        units = units if units else ds[runoff_var].attrs.get('units', False)
        if not units:
            log.warning("No units attribute found. Assuming meters")
            units = 'm'
        units = units.lower()
        if units in ('m', 'meters', 'kg m-2'):
            conversion_factor = 1
        elif units in ('mm', 'millimeters'):
            conversion_factor = .001
        else:
            raise ValueError(f"Unknown units: {ds.attrs['units']}")
        # get the time array from the dataset
        datetime_array = ds[time_var].to_numpy()

        if isinstance(runoff_var, list):
            ds = sum(*(ds[v] for v in runoff_var))
        else:
            ds = ds[runoff_var]
        if ds.ndim != 3:
            raise ValueError(f"Unknown number of dimensions: {ds.ndim}")

        vol_df = (
            ds
            .sel(**{
                time_var: slice(None),
                y_var: weight_df['lat_index'].values,
                x_var: weight_df['lon_index'].values,
            })
            .to_dataframe()
            .reset_index()
            .pivot(columns=river_id_var, index=time_var, values=runoff_var)
            .replace(np.nan, 0)
        )
    # This order of operations is important
    if cumulative:
        vol_df = _cumulative_to_incremental(vol_df)
    if force_positive_runoff:
        vol_df = vol_df.clip(lower=0)
    vol_df = vol_df * weight_df['area_sqm'].values * conversion_factor
    vol_df = vol_df.T.groupby(by=weight_df['river_id'].values).sum().T
    vol_df = vol_df[sorted_rivid_array]

    # Check that all time steps are the same
    time_diff = np.diff(datetime_array)
    if not np.all(time_diff == datetime_array[1] - datetime_array[0]) and force_uniform_timesteps:
        timestep = (datetime_array[1] - datetime_array[0]).astype('timedelta64[s]').astype(int)
        log.warning(f'Time steps are not uniform, resampling to the first timestep: {timestep} seconds')
        # forced to incremental before this step, get cumulative values then linear resample
        vol_df = (
            _incremental_to_cumulative(vol_df)
            .resample(rule=f'{timestep}S')
            .interpolate(method='linear')
        )
        vol_df = _cumulative_to_incremental(vol_df)

    return vol_df


def write_catchment_volumes(vol_df: pd.DataFrame, output_dir: str, label: str = None) -> None:
    """
    Write the catchment runoff volumes to file in the river-route expected format.

    Args:
        vol_df: pd.DataFrame
            The inflow data with stream IDs as columns and a datetime index
        output_dir: str
            The directory to write the file to
        label: str
            An optional label to include in the file name for organization purposes
    """
    # Create output inflow netcdf data
    os.makedirs(output_dir, exist_ok=True)
    start_date = vol_df.index[0].strftime('%Y%m%d%H')
    end_date = vol_df.index[-1].strftime('%Y%m%d%H')
    file_name = f'volumes{f"_{label}" if label else ""}_{start_date}_{end_date}.nc'
    inflow_file_path = os.path.join(output_dir, file_name)

    with nc.Dataset(inflow_file_path, "w", format="NETCDF4") as ds:
        ds.createDimension('time', vol_df.shape[0])
        ds.createDimension('river_id', vol_df.shape[1])

        ro_vol_var = ds.createVariable('volume', 'f4', ('time', 'river_id'), zlib=True, complevel=9)
        ro_vol_var[:] = vol_df.to_numpy()
        ro_vol_var.long_name = 'Incremental catchment runoff volume'
        ro_vol_var.units = 'm3'

        id_var = ds.createVariable('river_id', 'i4', ('river_id',), zlib=True, complevel=9)
        id_var[:] = vol_df.columns.astype(int)
        id_var.long_name = 'unique ID number for each river'

        timestep = (vol_df.index[1] - vol_df.index[0]).seconds
        time_var = ds.createVariable('time', 'i4', ('time',), zlib=True, complevel=9)
        time_var[:] = (vol_df.index - vol_df.index[0]).astype('timedelta64[s]').astype(int)
        time_var.long_name = 'time'
        time_var.standard_name = 'time'
        time_var.units = f'seconds since {vol_df.index[0].strftime("%Y-%m-%d %H:%M:%S")}'
        time_var.axis = 'T'
        time_var.time_step = f'{timestep}'
    return
