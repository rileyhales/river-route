import argparse

import river_route as rr


def main():
    parser = argparse.ArgumentParser(
        description='CLI for the Python river_route package for performing river routing on large river networks',
    )

    # parse subcommands for route and uh
    commands_subparser = parser.add_subparsers(dest='subcommand', required=True)

    # route subcommand
    parser_route = commands_subparser.add_parser('muskingum', help='Route inflows through a river network')
    parser_route.add_argument(
        '-c', '--config',
        type=str,
        help='Path to routing configuration file',
        default=None,
    )
    args = parser.parse_args()

    if args.subcommand == 'muskingum':
        rr.Muskingum(args.config).route()
    else:
        parser.print_usage()

    return


if __name__ == '__main__':
    main()
