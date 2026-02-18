import logging
from pathlib import Path

import geopandas as gpd
import numpy as np
import xarray as xr

import river_route as rr

data_root = Path(__file__).resolve().parent / 'data'

# log to stdout with level DEBUG
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
log.addHandler(ch)


def compare_parquets(file1: Path, file2: Path, **kwargs) -> bool:
    df1 = gpd.read_parquet(file1)
    df2 = gpd.read_parquet(file2)
    return df1.equals(df2)


def compare_netcdfs(file1: Path, file2: Path, variable: str | list[str] | None = None) -> bool:
    ds1 = xr.open_dataset(file1)
    ds2 = xr.open_dataset(file2)
    ds1vars = set(ds1.data_vars)
    ds2vars = set(ds2.data_vars)

    if variable is None:
        log.warning('No variable specified, comparing only variables in common between the two datasets')
        log.warning(f'variables in {file1}: {ds1vars}')
        log.warning(f'variables in {file2}: {ds2vars}')
        variable = list(ds1vars.intersection(ds2vars))
    if isinstance(variable, str):
        if not variable in ds1vars or not variable in ds2vars:
            log.error(f'variable {variable} not found in both datasets. Cannot perform comparison.')
            return False
        variable = [variable]

    is_identical = True
    for var in variable:
        # check the shapes and for differences less than 0.1
        if ds1[var].equals(ds2[var]):
            continue
        is_identical = False
        log.error(f'variable {var} is not a perfect match between {file1} and {file2}')
        if ds1[var].shape != ds2[var].shape:
            log.error(f'variable {var} has different shapes: {ds1[var].shape} vs {ds2[var].shape}')
        if np.any((ds1[var].values - ds2[var].values) > 0.1):
            log.error(f'variable {var} has values that differ by more than 0.1')
    return is_identical


def compare_with_logs(function, file1, file2):
    log.info(f'Comparing {file1} and {file2}')
    log.debug(f'New file: {file1}')
    log.debug(f'Expected file: {file2}')
    if function(file1, file2):
        log.info(f'[OK] {file1} and {file2} are identical')
    else:
        log.error(f'[FAIL] {file1} and {file2} are not exactly identical')


def core_muskingum_feature_set() -> None:
    # todo write final state files and compare final states
    # todo use an initial state file and test that the initial values are used correctly

    # inputs
    routing_params_file = data_root / 'sample_watershed' / 'prepared_streams.parquet'
    catchments_file = data_root / 'sample_watershed' / 'prepared_catchments.parquet'
    runoff_depths_file = data_root / 'sample_watershed' / 'prepared_runoff_depths.nc'
    # files created
    voroni_diagram_file = data_root / 'generated' / 'voronoi.parquet'
    grid_weights_file = data_root / 'generated' / 'grid_weights.nc'
    catchment_volumes_file = data_root / 'generated' / 'volumes.nc'
    discharge_from_volumes_file = data_root / 'generated' / 'discharge.nc'
    discharge_from_depths_file = data_root / 'generated' / 'discharge.nc'
    # solution
    expected_grid_weights_file = data_root / 'solutions' / 'grid_weights.nc'
    expected_volumes_file = data_root / 'solutions' / 'volumes.nc'
    expected_voroni_diagram_file = data_root / 'solutions' / 'voronoi.parquet'
    expected_discharge_file = data_root / 'solutions' / 'discharge.nc'

    runoff_grid_variable_name_params = {
        'x_var': 'longitude',
        'y_var': 'latitude',
        'time_var': 'valid_time'
    }

    ############# 1 create a voroni diagram and weight table using rr.grid_weights
    log.info('Computing voroni diagram and grid weights using rr.runoff.grid_weights')
    rr.runoff.grid_weights(
        runoff_depths_file, catchments_file,
        x_var=runoff_grid_variable_name_params['x_var'],
        y_var=runoff_grid_variable_name_params['y_var'],
        save_voroni_path=voroni_diagram_file,
        save_weights_path=grid_weights_file,
    )
    compare_with_logs(compare_parquets, voroni_diagram_file, expected_voroni_diagram_file)
    compare_with_logs(compare_netcdfs, grid_weights_file, expected_grid_weights_file)

    ############# 2 compute catchment volumes from gridded depths and grid weights
    log.info('Computing catchment volumes from gridded depths and grid weights')
    (
        rr
        .runoff
        .depth_to_volume(
            runoff_depths_file,
            grid_weights_file,
            x_var=runoff_grid_variable_name_params['x_var'],
            y_var=runoff_grid_variable_name_params['y_var'],
            time_var=runoff_grid_variable_name_params['time_var'],
        )
        .to_netcdf(catchment_volumes_file)
    )
    compare_with_logs(compare_netcdfs, catchment_volumes_file, expected_volumes_file)

    ############# 3 route catchment volumes to discharge using rr.LumpedMuskingum
    log.info('Routing catchment volumes to discharge using rr.LumpedMuskingum')
    (
        rr
        .Muskingum(**{
            'routing_params_file': routing_params_file,
            'catchment_volumes_files': catchment_volumes_file,
            'discharge_files': discharge_from_volumes_file,
            'log': True,
        })
        .route()
    )
    # 4 route depths and weights to discharge using rr.LumpedMuskingum
    (
        rr
        .Muskingum(**{
            'routing_params_file': routing_params_file,
            'weight_table_file': grid_weights_file,
            'runoff_depths_files': runoff_depths_file,
            'var_x': runoff_grid_variable_name_params['x_var'],
            'var_y': runoff_grid_variable_name_params['y_var'],
            'var_t': runoff_grid_variable_name_params['time_var'],
            'discharge_files': discharge_from_depths_file,
            'log': True,
        })
        .route()
    )
    compare_with_logs(compare_netcdfs, discharge_from_volumes_file, expected_discharge_file)
    compare_with_logs(compare_netcdfs, discharge_from_depths_file, expected_discharge_file)


if __name__ == '__main__':
    core_muskingum_feature_set()
