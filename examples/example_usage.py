import argparse
import contextlib
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import river_route as rr

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--region', default='6020006540')
    parser.add_argument('--regions-root', default=Path.home() / 'data' / 'v3TestData' / 'hydrography')
    parser.add_argument('--runoff-root', default=Path.home() / 'data' / 'era5_zarr_8x8_3month')
    parser.add_argument('--dt-routing', type=int, default=3600)
    parser.add_argument('--threads', type=int, default=12)
    args = parser.parse_args()

    runoff_files = sorted(list(Path(args.runoff_root).glob('year=*/*.zarr')))
    null_outputs = [Path(os.devnull) / f.name for f in runoff_files]  # must set something to pass config validation
    region_dir = Path(args.regions_root)
    params_file = region_dir / f'region={args.region}' / 'routing.parquet'
    grid_weights_file = region_dir / f'region={args.region}' / f'gridweights_ERA5_{args.region}.nc'

    # ── 1. Configs: every option, frozen, in one object ──────────────────────
    conf = rr.Configs(
        params_file=params_file,
        grid_runoff_files=runoff_files,
        grid_weights_file=grid_weights_file,
        discharge_files=null_outputs,
        forcing='vlateral',
        coeff='static',
        routing_order='river',
        dt_routing=args.dt_routing,
        runoff_processing_mode='sequential',
        var_grid_runoff='ro',
        var_x='longitude',
        var_y='latitude',
        var_t='valid_time',
        progress_bar=True,
        log_level='ERROR',
        unstable_coefficients='ignore',
    )

    conf.validate_routing()
    conf.validate_runoff()
    # deep_validate() reads every input file and checks its contents. It is a one off check to run on inputs you
    # have not used before, not part of a run: routing never calls it, and it costs more than the routing does.

    pool_context = ThreadPoolExecutor(args.threads) if args.threads > 1 else contextlib.nullcontext()
    with pool_context as pool:
        (
            rr.Router(conf)
            .set_discharge_writer(rr.router.writers.null_writer)
            .route(thread_pool=pool, threads=args.threads)
        )
