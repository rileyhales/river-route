import argparse

from .configs import Configs
from .network.streams import subset_network_to_river
from .router import Router


def main() -> None:
    """Parse the command line and run the rr subcommand it names."""
    parser = argparse.ArgumentParser(prog='rr', description='river-route: configurable Muskingum river routing')
    subparsers = parser.add_subparsers(dest='command')

    route = subparsers.add_parser('route', help='Run routing from a config file')
    route.add_argument('config', type=str, help='Path to routing configuration file (JSON)')

    subset = subparsers.add_parser('subset', help='Subset a network file and optionally its grid weights')
    subset.add_argument('river', type=int, help='riverId to subset to; it becomes the outlet')
    subset.add_argument('--network', type=str, required=True, help='Path to the full network file parquet')
    subset.add_argument('--out-network', type=str, required=True, help='Path to write the subset network file parquet')
    subset.add_argument('--weights', type=str, default=None, help='Path to the full grid weights netCDF')
    subset.add_argument('--out-weights', type=str, default=None, help='Path to write the subset grid weights netCDF')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    if args.command == 'route':
        Router(Configs.from_json(args.config)).route()

    if args.command == 'subset':
        if (args.weights is None) != (args.out_weights is None):
            subset.error('--weights and --out-weights must be given together')
        subset_network_to_river(args.river, args.network, args.out_network, args.weights, args.out_weights)


if __name__ == '__main__':
    main()
