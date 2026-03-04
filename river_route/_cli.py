import argparse
import sys

from .routers import Muskingum, RapidMuskingum, UnitMuskingum

try:
    import river_route_app
    APP_AVAILABLE = True
except ImportError:
    APP_AVAILABLE = False

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


def _run_app(args):
    if not APP_AVAILABLE:
        print('river-route-app is not installed.')
        print('Install it with:  pip install river-route[app]')
        sys.exit(1)
    import uvicorn
    port = args.port
    uvicorn.run('river_route_app.server:app', host='127.0.0.1', port=port)


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

    rapid = subparsers.add_parser('RapidMuskingum', help='RAPID-style Muskingum routing with lateral runoff')
    _add_config_arg(rapid)

    unit = subparsers.add_parser('UnitMuskingum', help='Unit hydrograph transform then Muskingum routing')
    _add_config_arg(unit)

    channel = subparsers.add_parser('Muskingum', help='Channel-only Muskingum routing (no lateral inflow)')
    _add_config_arg(channel)

    app = subparsers.add_parser('app', help='Launch the browser-based app (requires pip install river-route[app])')
    app.add_argument('--port', type=int, default=8000, help='Port to listen on (default: 8000)')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    if args.command == 'app':
        _run_app(args)
    elif args.command == 'route':
        _route(args)
    else:
        ROUTERS[args.command](args.config).route()
