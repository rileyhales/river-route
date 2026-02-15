import river_route as rr


def test_route_volumes_matches_known_output() -> None:
    # inputs
    routing_params_file = './data/sample_watershed/prepared_streams.parquet'
    catchments_file = './data/sample_watershed/prepared_catchments.parquet'
    runoff_depths_file = './data/sample_watershed/prepared_runoff_depths.nc'
    # files created
    voroni_grid_diagram_file = './data/generated/voroni_grid_diagram.parquet'
    grid_weights_file = './data/generated/grid_weights.nc'
    catchment_volumes_file = './data/generated/catchment_volumes.nc'
    discharge_from_volumes_file = './data/generated/discharge_from_volumes.nc'
    discharge_from_depths_file = './data/generated/discharge_from_depths.nc'
    # solution
    expected_grid_weights_file = './data/solutions/grid_weights.nc'
    expected_discharge_file = './data/solutions/discharge.nc'

    # 1 create a weight table using rr.grid_weights
    rr.runoff.grid_weights(
        runoff_depths_file, catchments_file,
        save_voroni_path=voroni_grid_diagram_file,
        save_weights_path=grid_weights_file,
    )
    # todo compare weight table to expected weight table
    # todo compare the voroni diagram to expected voroni diagram

    # 2 compute catchment volumes from gridded depths and grid weights
    rr.runoff.depth_to_volume(runoff_depths_file, grid_weights_file, time_var='valid_time')
    # todo compare catchment volumes to expected catchment volumes

    # 3 route catchment volumes to discharge using rr.Muskingum
    (
        rr
        .Muskingum(**{
            'routing_params_file': routing_params_file,
            'catchment_volumes_files': catchment_volumes_file,
            'discharge_files': discharge_from_volumes_file,
        })
        .route()
    )
    # 4 route depths and weights to discharge using rr.Muskingum
    (
        rr
        .Muskingum(**{
            'routing_params_file': routing_params_file,
            'grid_weights_file': grid_weights_file,
            'runoff_depths_files': runoff_depths_file,
            'discharge_files': discharge_from_depths_file,
        })
        .route()
    )
    # todo compare discharge from volumes to expected discharge
    # todo compare discharge from depths to expected discharge
