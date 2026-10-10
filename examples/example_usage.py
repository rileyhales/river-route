import argparse
import contextlib
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import river_route as rr

if __name__ == '__main__':
    # Amazon: 6020006540, Caribbean: 7020065090
    parser = argparse.ArgumentParser()
    parser.add_argument('--bucket-root', default=Path.home() / 'data' / 'rfsv3')
    parser.add_argument('--region', default='6020006540')
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()

    # find the routing files for the region
    routing_root = Path(args.bucket_root) / 'routing'
    network_file = routing_root / f'region={args.region}' / 'network.parquet'
    grid_weights_file = routing_root / f'region={args.region}' / f'gridweights_ERA5_{args.region}.nc'
    # find the runoff forcing files
    era5_root = Path(args.bucket_root) / 'forcings' / 'era5'
    runoff_files = sorted(era5_root.glob('194[0-9].zarr'))
    # provide an output directory
    discharge_dir = Path(os.devnull)  # devnull will pass validation and results will be discarded

    # ── 1. Configs: every option, frozen, in one object ──────────────────────
    conf = rr.Configs(
        network_file=network_file,
        # primary modeling choices
        parameters='static',   # 'static' or 'dynamic'
        forcing='grid',        # 'channel', or how the forcing is provided, see docs
        transform='uniform',   # 'uniform', pending update will permit 'unit-hydrograph'
        dt_routing='auto',     # 'auto', or a positive, integer number of seconds
        # Calculation quality control options
        restrict_x_range=True,
        invalid_coefficients='warn',
        # model state files
        channel_state_init_file=None,
        channel_state_final_file=None,
        # where to place the results
        discharge_dir=discharge_dir,
        # the runoff forcing data and how to read it
        runoff_files=runoff_files,
        grid_weights_file=grid_weights_file,
        var_grid_runoff='ro',
        var_x='longitude',
        var_y='latitude',
        var_t='valid_time',
        # logging and validation options
        progress_bar=True,
        log_level='PROGRESS',
    )

    pool_context = ThreadPoolExecutor(args.threads) if args.threads > 1 else contextlib.nullcontext()
    with pool_context as pool:
        rr.Router(conf).route(thread_pool=pool, threads=args.threads)
