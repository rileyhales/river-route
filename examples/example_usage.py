import argparse
import contextlib
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import river_route as rr

if __name__ == '__main__':
    # Amazon: 6020006540, Caribbean: 7020065090
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--bucket-root',
        required=True,
        help='folder holding routing/region=<region>/ with its network and ERA5 grid weights, and forcings/era5/',
    )
    parser.add_argument('--region', default='6020006540', help='the region to route')
    parser.add_argument('--threads', type=int, default=3, help='threads to route with')
    args = parser.parse_args()

    # find the routing files for the region
    routing_root = Path(args.bucket_root) / 'routing'
    network_file = routing_root / f'region={args.region}' / 'network.parquet'
    grid_weights_file = routing_root / f'region={args.region}' / f'gridweights_ERA5_{args.region}.nc'
    # find the runoff forcing files
    era5_root = Path(args.bucket_root) / 'forcings' / 'era5'
    forcing_files = sorted(era5_root.glob('194[0-9].zarr'))
    # provide an output directory
    discharge_dir = Path(os.devnull)  # devnull will pass validation and results will be discarded

    conf = rr.Configs(
        network_file=network_file,
        # primary modeling choices
        forcing='grid',        # 'channel', or how the forcing is provided, see docs
        # compute options
        dt_routing='auto',     # 'auto', or a positive, integer number of seconds
        use_subreaches=True,   # split rivers too long for their step into sub-reaches; 'auto' requires it
        use_subcycles=True,    # route rivers too short for their step in subcycles; 'auto' requires it
        invalid_coefficients='raise',
        # model state files
        channel_state_init_file=None,
        channel_state_final_file=None,
        # outputs save location
        discharge_dir=discharge_dir,
        # the forcing files and how to read them
        forcing_files=forcing_files,
        grid_weights_file=grid_weights_file,
        var_runoff='ro',
        var_x='longitude',
        var_y='latitude',
        var_t='valid_time',
        # logging and validation options
        progress_bar=True,
        log_level='PROGRESS',
    )

    pool_context = ThreadPoolExecutor(args.threads) if args.threads > 1 else contextlib.nullcontext()
    with pool_context as pool:
        # 'auto' routes with the C kernel once rr compile has compiled it, and with numba until then
        rr.Router(conf).route(thread_pool=pool, threads=args.threads, kernel='auto')
