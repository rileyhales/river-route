import argparse
import contextlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

from .configs import Configs
from .kernels.build import build
from .network.streams import subset_network_to_river
from .router import Router

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='rr', description='river-route: configurable Muskingum river routing')
    subparsers = parser.add_subparsers(dest='command')

    route = subparsers.add_parser('route', help='Run routing from a config file')
    route.add_argument('config', type=str, help='Path to routing configuration file (JSON), or - to read it from stdin')
    route.add_argument('--threads', type=int, default=1, help='Threads to route with')

    subset = subparsers.add_parser('subset', help='Subset a network file and optionally its grid weights')
    subset.add_argument('river', type=int, help='riverId to subset to; it becomes the outlet')
    subset.add_argument('--network', type=str, required=True, help='Path to the full network file parquet')
    subset.add_argument('--out-network', type=str, required=True, help='Path to write the subset network file parquet')
    subset.add_argument('--weights', type=str, default=None, help='Path to the full grid weights netCDF')
    subset.add_argument('--out-weights', type=str, default=None, help='Path to write the subset grid weights netCDF')

    compile_ = subparsers.add_parser('compile', help='Compile the C routing kernel for this machine with LLVM')
    compile_.add_argument('--llvm-bin', type=str, default=None, help='Path to LLVM_BIN')
    compile_.add_argument('--checks', action='store_true', help='Require no warnings and pass clang-tidy')
    compile_.add_argument('--wasm', type=str, default=None, metavar='FOLDER', help='Compile WebAssembly into FOLDER')

    app = subparsers.add_parser('app', help='Open a window that sets up and runs routing')
    app.add_argument('--home', type=str, default=None, help='Folder every file dialog opens in (default: the current)')
    app.add_argument('--threads', type=int, default=1, help='Threads each run routes with, until changed (default: 1)')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        exit(0)

    if args.command == 'route':
        if args.threads < 1:
            route.error('--threads must be 1 or more')
        configs = Configs(**json.load(sys.stdin)) if args.config == '-' else Configs.from_json(args.config)
        with ThreadPoolExecutor(args.threads) if args.threads > 1 else contextlib.nullcontext() as pool:
            Router(configs).route(thread_pool=pool, threads=args.threads)

    if args.command == 'subset':
        if (args.weights is None) != (args.out_weights is None):
            subset.error('--weights and --out-weights must be given together')
        subset_network_to_river(args.river, args.network, args.out_network, args.weights, args.out_weights)

    if args.command == 'compile':
        try:
            build(args.llvm_bin, args.checks, args.wasm)
        except (FileNotFoundError, NotADirectoryError, subprocess.CalledProcessError) as error:
            compile_.error(str(error))

    if args.command == 'app':
        try:
            from .app import open_window

            open_window(args.home, args.threads)
        except (NotADirectoryError, ValueError, ImportError) as error:
            app.error(str(error))
