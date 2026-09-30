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

    dt_routing = 3600  # all times in seconds, 1 hour

    # find the routing files for the region
    routing_root = Path(args.bucket_root) / 'routing'
    network_file = routing_root / f'region={args.region}' / 'routing.parquet'
    grid_weights_file = routing_root / f'region={args.region}' / f'gridweights_ERA5_{args.region}.nc'
    # find the runoff forcing files
    era5_root = Path(args.bucket_root) / 'forcings' / 'era5'
    runoff_files = sorted(era5_root.glob('194[0-9].zarr'))
    # provide an output directory
    discharge_dir = Path(os.devnull)  # devnull will pass validation and results will be discarded

    # ── 1. Configs: every option, frozen, in one object ──────────────────────
    conf = rr.Configs(
        # the network/network file
        network_file=network_file,
        # primary modeling choices
        coefficients='static',
        forcing='grid',
        transform='uniform',
        network_type='standard',
        dt_routing=dt_routing,
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
        log_level='ERROR',
        unstable_coefficients='ignore',
    )

    pool_context = ThreadPoolExecutor(args.threads) if args.threads > 1 else contextlib.nullcontext()
    with pool_context as pool:
        rr.Router(conf).route(thread_pool=pool, threads=args.threads)
