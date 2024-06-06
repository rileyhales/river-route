import glob
import logging
import os
import re

import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr

__all__ = [
    'calc_runoff_volumes',
    'write_runoff_volumes',
    'create_runoff_volumes_files',
]


def _cumulative_to_incremental(df) -> pd.DataFrame:
    return pd.DataFrame(
        np.vstack([df.values[0, :], np.diff(df.values, axis=0)]),
        index=df.index,
        columns=df.columns
    )


def _incremental_to_cumulative(df) -> pd.DataFrame:
    return df.cumsum()


def _guess_variable_name(var_name: str, possible_matches: list) -> str:
    if len(possible_matches) == 0:
        raise ValueError(f"No {var_name} variable found in LSM data. Check dataset or specify {var_name}_var")
    if len(possible_matches) == 1:
        return possible_matches[0]
    elif len(possible_matches) > 1:
        raise ValueError(f"Multiple {var_name} variables found. Specify with {var_name}_var: {possible_matches}")
    else:
        raise ValueError(f"Unexpected error finding {var_name} variable. Check dataset or specify {var_name}_var")


def _search_for_weight_table(input_dir: str, x0: int or float, y0: int or float, dx: int or float,
                             dy: int or float) -> str:
    # weight table name includes x0=<number>, y0=<number>, dx=<number>, dy=<number>
    tables = glob.glob(os.path.join(input_dir, f'weight_*'))
    tables = [x for x in tables if re.match(fr'weight_.*x0={x0}.*y0={y0}.*dx={dx}.*dy={dy}.*', x)]
    if not len(tables):
        logging.error(f'No weight table found in {input_dir} which matches the dataset')
        raise FileNotFoundError(f'Could not find a weight table in {input_dir} which matches the dataset')
    if len(tables) > 1:
        logging.error(f'Multiple weight tables found in {input_dir}. Please specify weight_table')
        raise ValueError(f'Multiple weight tables found in {input_dir}. Please specify weight_table')
    return tables[0]


def calc_runoff_volumes(
        runoff_data: str,
        vpu_dir: str,
        weight_table: str = None,
        cumulative: bool = False,
        force_positive_runoff: bool = False,
        force_uniform_timesteps: bool = True,
        runoff_var: str = None,
        x_var: str = None,
        y_var: str = None,
        time_var: str = None,
) -> pd.DataFrame:
    """
    Calculates the catchment runoff volumes from a given runoff dataset and a directory of VPU configs.

    Args:
        runoff_data: a string or list of strings of paths to LSM files
        vpu_dir: a string path to the directory of VPU files which contains weight tables
        weight_table: a string path to the weight table to override the automatic search
        cumulative: whether the provided runoff data is cumulative (else, incremental)
        force_positive_runoff: whether to replace negative runoff with zero
        force_uniform_timesteps: whether to linearly resample the data to uniform timesteps
        runoff_var: the name of the runoff variable in the LSM data
        x_var: the name of the x variable in the LSM data
        y_var: the name of the y variable in the LSM data
        time_var: the name of the time variable in the LSM data

    Returns:
        pd.DataFrame: a DataFrame of the catchment runoff volumes with stream IDs as columns and a datetime index
    """
    # open all the ncs and select only the area within the weight table
    if type(runoff_data) is list:
        # this is correct, a list of files is allowed
        assert all([os.path.exists(x) for x in runoff_data]), 'Not all files in the list exist'
    elif os.path.isdir(runoff_data):
        logging.warning(f'{runoff_data} is a directory. Guessing which files to use.')
        runoff_data = os.path.join(runoff_data, '*.nc*')
    elif os.path.isfile(runoff_data):
        ...  # this is correct, a single file is allowed
    elif '*' in runoff_data:
        ...  # this is correct, xarray will interpret the glob sting independently
    else:
        raise FileNotFoundError(f'{runoff_data} does not exist and is not a glob pattern')

    with xr.open_mfdataset(runoff_data) as ds:
        # Select the variable names
        if not runoff_var:
            logging.warning('Runoff variable not given. Guessing from default names')
            runoff_var = [x for x in ['ro', 'RO', 'runoff', 'RUNOFF'] if x in ds.variables]
            runoff_var = _guess_variable_name('runoff', runoff_var)
        if not x_var:
            logging.warning('X variable not given. Guessing from default names.')
            x_var = [x for x in ['x', 'lon', 'longitude', 'LONGITUDE', 'LON'] if x in ds.variables]
            x_var = _guess_variable_name('x', x_var)
        if not y_var:
            logging.warning('Y variable not given. Guessing from default names.')
            y_var = [x for x in ['y', 'lat', 'latitude', 'LATITUDE', 'LAT'] if x in ds.variables]
            y_var = _guess_variable_name('y', y_var)
        if not time_var:
            logging.warning('Time variable not given. Guessing from default names.')
            time_var = [x for x in ['time', 'TIME', ] if x in ds.variables]
            time_var = _guess_variable_name('time', time_var)

        # Check units and conversion factors
        conversion_factor = 1
        units = ds.attrs.get('units', False)
        if not units:
            logging.warning("No units attribute found. Assuming meters")
        elif ds.attrs['units'] == 'm':
            conversion_factor = 1
        elif ds.attrs['units'] == 'mm':
            conversion_factor = .001
        else:
            raise ValueError(f"Unknown units: {ds.attrs['units']}")

        # get the array shape descriptors
        dx = ds[x_var].values[1] - ds[x_var].values[0]
        x0 = ds[x_var].values[0]
        dy = ds[y_var].values[1] - ds[y_var].values[0]
        y0 = ds[y_var].values[0]

        # Find the right weight table if not specified
        if not weight_table:
            weight_table = _search_for_weight_table(vpu_dir, x0, y0, dx, dy)

        # load in weight table and get some information
        weight_df = pd.read_csv(weight_table)

        # todo replace with routing params file or the order of the entries in weight_df
        comid_df = pd.read_csv(os.path.join(vpu_dir, 'comid_lat_lon_z.csv'))
        sorted_rivid_array = comid_df.iloc[:, 0].to_numpy()

        # for readability, select certain cols from the weight table
        stream_ids = weight_df.iloc[:, 0].to_numpy()
        ds = ds[runoff_var]

        # get the time array from the dataset
        datetime_array = ds[time_var].to_numpy()

        # todo use Dataset.sel(**kwargs) using a dictionary of time, y_var, x_var key/values
        if ds.ndim == 3:
            vol_df = ds.values[:, weight_df['lat_index'].values, weight_df['lon_index'].values]
        elif ds.ndim == 4:
            vol_df = ds.values[:, :, weight_df['lat_index'].values, weight_df['lon_index'].values]
            vol_df = np.where(np.isnan(vol_df[:, 0, :]), vol_df[:, 1, :], vol_df[:, 0, :]),
        else:
            raise ValueError(f"Unknown number of dimensions: {ds.ndim}")

    # This order of operations is important
    vol_df = pd.DataFrame(vol_df, columns=stream_ids, index=datetime_array)
    vol_df = vol_df.replace(np.nan, 0)
    if cumulative:
        vol_df = _cumulative_to_incremental(vol_df)
    if force_positive_runoff:
        vol_df = vol_df.clip(lower=0)
    vol_df = vol_df * weight_df['area_sqm'].values * conversion_factor
    vol_df = vol_df.T.groupby(by=stream_ids).sum().T
    vol_df = vol_df[sorted_rivid_array]

    # Check that all time steps are the same
    time_diff = np.diff(datetime_array)
    if not np.all(time_diff == datetime_array[1] - datetime_array[0]) and force_uniform_timesteps:
        timestep = (datetime_array[1] - datetime_array[0]).astype('timedelta64[s]').astype(int)
        logging.warning(f'Time steps are not uniform, resampling to the first timestep: {timestep} seconds')
        # everything is forced to be incremental before this step so we can use cumsum to get the cumulative values
        vol_df = (
            _incremental_to_cumulative(vol_df)
            .resample(rule=f'{timestep}S')
            .interpolate(method='linear')
        )
        vol_df = _cumulative_to_incremental(vol_df)

    return vol_df


def write_runoff_volumes(vol_df: pd.DataFrame, output_dir: str, vpu_name: str, file_label: str = None) -> None:
    """
    Write the catchment runoff volumes to file in the river-route expected format.

    Args:
        vol_df: pd.DataFrame
            The inflow data with stream IDs as columns and a datetime index
        output_dir: str
            The directory to write the file to
        vpu_name: str
            The name of the VPU
        file_label: str
            A label to include in the file name for organization purposes
    """
    # Create output inflow netcdf data
    os.makedirs(output_dir, exist_ok=True)
    start_date = vol_df.index[0].strftime('%Y%m%d')
    end_date = vol_df.index[-1].strftime('%Y%m%d')
    file_name = f'm3_{vpu_name}_{start_date}_{end_date}.nc'
    if file_label is not None:
        file_name = f'm3_{vpu_name}_{start_date}_{end_date}_{file_label}.nc'
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


def create_runoff_volumes_files(runoff_data: str,
                                vpu_dir: str,
                                output_dir: str,
                                vpu_name: str = None,
                                weight_table: str = None,
                                cumulative: bool = False,
                                file_label: str = None,
                                force_positive_runoff: bool = False,
                                runoff_var: str = None,
                                x_var: str = None,
                                y_var: str = None,
                                time_var: str = None, ) -> None:
    """
    Create the inflow files for a given VPU, from a given runoff dataset, and write them to disc.
    Args:
        runoff_data: a string or list of strings of paths to LSM files
        vpu_dir: a string path to the directory of VPU files which contains weight tables
        output_dir: a string path to the directory to write the inflow files to
        vpu_name: a string name of the VPU
        weight_table: a string path to the weight table to override the automatic search
        cumulative: whether the provided runoff data is cumulative (else, incremental)
        file_label: a string label to include in the file name for organization purposes
        force_positive_runoff: whether to replace negative runoff with zero
        runoff_var: the name of the runoff variable in the LSM data
        x_var: the name of the x variable in the LSM data
        y_var: the name of the y variable in the LSM data
        time_var: the name of the time variable in the LSM data

    Returns:
        None
    """
    vol_df = calc_runoff_volumes(runoff_data=runoff_data, vpu_dir=vpu_dir, weight_table=weight_table,
                                 cumulative=cumulative, force_positive_runoff=force_positive_runoff,
                                 force_uniform_timesteps=True, runoff_var=runoff_var, x_var=x_var, y_var=y_var,
                                 time_var=time_var)

    vpu_name = vpu_name if vpu_name is not None else os.path.basename(vpu_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Write the inflow file
    write_runoff_volumes(vol_df, output_dir, vpu_name, file_label)
    return
