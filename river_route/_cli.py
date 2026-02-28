import argparse

from .routers import Router, RapidMuskingum, UnitMuskingum

ROUTERS = {
    'rapid': RapidMuskingum,
    'unit': UnitMuskingum,
    'channel': Router,
}


def _add_config_arg(subparser: argparse.ArgumentParser) -> None:
    subparser.add_argument('config', type=str, help='Path to routing configuration file')


def main():
    parser = argparse.ArgumentParser(
        prog='rr',
        description='river-route: vectorized Muskingum river routing',
    )
    subparsers = parser.add_subparsers(dest='command')

    rapid = subparsers.add_parser('rapid', help='RAPID-style Muskingum routing with lateral runoff')
    _add_config_arg(rapid)

    unit = subparsers.add_parser('unit', help='Unit hydrograph transform then Muskingum routing')
    _add_config_arg(unit)

    channel = subparsers.add_parser('channel', help='Channel-only Muskingum routing (no lateral inflow)')
    _add_config_arg(channel)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    ROUTERS[args.command](args.config).route()
