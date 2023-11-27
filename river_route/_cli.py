import argparse

import river_route as rr


def main():
    parser = argparse.ArgumentParser(
        description='CLI for the Python river_route package for performing river routing on large river networks',
        add_help=True,
    )
    parser.add_argument(
        '-c', '--config',
        type=str,
        help='Path to routing configuration file',
        default=None,
    )
    args = parser.parse_args()

    if args.config is not None:
        rr.Muskingum(args.config).route()
    else:
        parser.print_help()

    return
