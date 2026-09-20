import argparse

from .configs import Configs
from .network.streams import subset_configs_to_river
from .router import Router


def main():
    parser = argparse.ArgumentParser(prog='rr', description='river-route: configurable Muskingum river routing')
    subparsers = parser.add_subparsers(dest='command')

    route = subparsers.add_parser(
        'route',
        help='Run routing from a config file. The procedure is selected by the '
        'coeff/forcing/network keys in the config.',
    )
    route.add_argument('config', type=str, help='Path to routing configuration file (YAML or JSON)')

    subset = subparsers.add_parser(
        'subset',
        help='Subset a routing parameter table, and optionally its grid weight table, to one river and '
        'everything upstream of it. The target river becomes the outlet of the subset.',
    )
    subset.add_argument('river', type=int, help='river_id to subset to; it becomes the outlet')
    subset.add_argument('params', type=str, help='Path to the full routing parameters parquet file')
    subset.add_argument('out_params', type=str, help='Path to write the subsetted parameters parquet file')
    subset.add_argument('--weights', type=str, default=None, help='Path to the full grid weights netCDF file')
    subset.add_argument(
        '--out-weights', type=str, default=None, help='Path to write the subsetted grid weights netCDF file'
    )

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    if args.command == 'route':
        Router(Configs.from_file(args.config)).route()

    if args.command == 'subset':
        if (args.weights is None) != (args.out_weights is None):
            subset.error('--weights and --out-weights are given together: one is useless without the other')
        subset_configs_to_river(args.river, args.params, args.out_params, args.weights, args.out_weights)


if __name__ == '__main__':
    main()
