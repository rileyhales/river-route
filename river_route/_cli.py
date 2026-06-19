import argparse

from .routers import Router


def main():
    parser = argparse.ArgumentParser(
        prog='rr',
        description='river-route: configurable Muskingum river routing',
    )
    subparsers = parser.add_subparsers(dest='command')

    route = subparsers.add_parser(
        'route',
        help='Run routing from a config file. The procedure is selected by the '
             'coeff/forcing/network keys in the config.',
    )
    route.add_argument('config', type=str, help='Path to routing configuration file (YAML or JSON)')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    if args.command == 'route':
        Router(args.config).route()


if __name__ == '__main__':
    main()
