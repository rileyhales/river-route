import argparse

import river_route as rr


def main():
    """
    CLI for the river-route Python package.

    Parameters:
        -h --help: Show this help message and exit
        -c --config: Path to routing configuration file

    Examples:
        ```bash
        rr -c path/to/config.yml
        ```
    """
    parser = argparse.ArgumentParser(
        description='CLI for the river_route Python package',
        add_help=True,
    )
    parser.add_argument(
        '-c', '--config',
        type=str,
        help='path to routing configuration file',
        default=None,
        required=False,
    )

    args = parser.parse_args()

    if not args.config:
        return parser.print_help()

    rr.MuskingumCunge(args.config).route()
    return
