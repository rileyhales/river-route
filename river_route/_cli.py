import argparse

from .configs import Configs
from .network.streams import subset_configs_to_river
from .router import Router


def main():
    parser = argparse.ArgumentParser(prog='rr', description='river-route: configurable Muskingum river routing')
    subparsers = parser.add_subparsers(dest='command')

    route = subparsers.add_parser('route', help='Run routing from a config file')
    route.add_argument('config', type=str, help='Path to routing configuration file (JSON)')

    subset = subparsers.add_parser('subset', help='Subset a parameter table and optionally its grid weights')
    subset.add_argument('river', type=int, help='river_id to subset to; it becomes the outlet')
    subset.add_argument('--params', type=str, default=None, help='Path to the full parameters parquet')
    subset.add_argument('--out-params', type=str, default=None, help='Path to write the subset parameters parquet')
    subset.add_argument('--weights', type=str, default=None, help='Path to the full grid weights netCDF')
    subset.add_argument('--out-weights', type=str, default=None, help='Path to write the subset grid weights netCDF')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    if args.command == 'route':
        if args.config is None:
            route.error('Missing required argument: config')
        Router(Configs.from_json(args.config)).route()

    if args.command == 'subset':
        if (args.params is None) != (args.out_params is None):
            subset.error('--params and --out-params are given together: one is useless without the other')
        if (args.weights is None) != (args.out_weights is None):
            subset.error('--weights and --out-weights are given together: one is useless without the other')
        if (args.params is None) and (args.weights is None):
            subset.error('At least one of --params or --weights must be given to subset')
        subset_configs_to_river(args.river, args.params, args.out_params, args.weights, args.out_weights)


if __name__ == '__main__':
    main()
