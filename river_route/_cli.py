import argparse

import river_route as rr


def main():
    parser = argparse.ArgumentParser(
        description='CLI for the river_route Python package for performing river routing on large river networks',
        add_help=True,
    )
    parser.add_argument(
        '-c', '--config',
        type=str,
        help='Path to routing configuration file',
        default=None,
        required=False,
    )
    args = parser.parse_args()

    if args.config:
        return rr.MuskingumCunge(args.config).route()
    return parser.print_help()
