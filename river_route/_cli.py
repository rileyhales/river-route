import argparse

import river_route as rr


def main():
    p = argparse.ArgumentParser(
        description='CLI for the river_route Python package for performing river routing on large river networks',
        add_help=True,
    )
    p.add_argument(
        '-c', '--config',
        type=str,
        help='Path to routing configuration file',
        default=None,
        required=False,
    )

    args = p.parse_args()

    if not args.config:
        return p.print_help()

    rr.MuskingumCunge(args.config).route()
    return
