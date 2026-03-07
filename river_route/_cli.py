import argparse
import sys

from .routers import Muskingum, RapidMuskingum, UnitMuskingum

ROUTERS = {
    'Muskingum': Muskingum,
    'RapidMuskingum': RapidMuskingum,
    'UnitMuskingum': UnitMuskingum,
}


def _add_config_arg(subparser: argparse.ArgumentParser) -> None:
    subparser.add_argument('config', type=str, help='Path to routing configuration file')


def _route(args):
    """Run routing from a config file using the router specified by --router."""
    router_cls = ROUTERS.get(args.router)
    if router_cls is None:
        print(f'Unknown router: {args.router!r}. Must be one of: {", ".join(ROUTERS)}')
        sys.exit(1)

    router_cls(args.config).route()


def main():
    parser = argparse.ArgumentParser(
        prog='rr',
        description='river-route: vectorized Muskingum river routing',
    )
    subparsers = parser.add_subparsers(dest='command')

    route = subparsers.add_parser('route', help='Run routing from a config file with a specified router')
    _add_config_arg(route)
    route.add_argument('--router', type=str, required=True,
                       choices=list(ROUTERS.keys()),
                       help='Router class to use (Muskingum, RapidMuskingum, or UnitMuskingum)')

    channel = subparsers.add_parser('Muskingum', help='Channel-only Muskingum routing (no lateral inflow)')
    _add_config_arg(channel)

    rapid = subparsers.add_parser('RapidMuskingum', help='RAPID-style Muskingum routing with lateral runoff')
    _add_config_arg(rapid)

    unit = subparsers.add_parser('UnitMuskingum', help='Unit hydrograph transform then Muskingum routing')
    _add_config_arg(unit)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    if args.command == 'route':
        _route(args)
    else:
        ROUTERS[args.command](args.config).route()
