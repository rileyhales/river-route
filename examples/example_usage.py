import argparse
import contextlib
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import river_route as rr

if __name__ == '__main__':
    # Amazon: 6020006540, Caribbean: 7020065090, Mississippi:
    parser = argparse.ArgumentParser()
    parser.add_argument('--region', default='6020006540')
    parser.add_argument('--regions-root', default=Path.home() / 'data' / 'rfsv3' / 'hydrography')
    parser.add_argument('--runoff-root', default=Path.home() / 'data' / 'era5_zarr_16x16_12month')
    parser.add_argument('--dt-routing', type=int, default=3600)
    parser.add_argument('--network-conditioning', choices=('standard', 'stabilized'), default='standard')
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()

    runoff_files = sorted(list(Path(args.runoff_root).glob('year=*/*.zarr')))
    discharge_dir = Path(os.devnull)         # must pass something to pass config validation
    region_dir = Path(args.regions_root)
    params_file = region_dir / f'region={args.region}' / 'routing.parquet'
    grid_weights_file = region_dir / f'region={args.region}' / f'gridweights_ERA5_{args.region}.nc'

    # ── 1. Configs: every option, frozen, in one object ──────────────────────
    conf = rr.Configs(
        # the network/routing parameters file
        params_file=params_file,
        network_type=args.network_conditioning,
        # primary modeling choices
        coefficients='static',
        forcing='runoff',
        transform='uniform',
        dt_routing=args.dt_routing,
        # model state files
        channel_state_init_file=None,
        channel_state_final_file=None,
        # where to place the results
        discharge_dir=discharge_dir,
        # the runoff forcing data and how to read it
        runoff_type='gaussian_grid',
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
        (
            rr.Router(conf)
            .route(thread_pool=pool, threads=args.threads)
        )
