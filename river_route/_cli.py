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


def _launch_app(args):
    """Launch river-route-app from the rr namespace if installed."""
    try:
        import uvicorn
        from river_route_app.server import app as rr_app
    except Exception as e:
        print(
            'rr app requires river-route-app to be installed in this environment.\n'
            'Install it with: pip install river-route-app'
        )
        print(f'Details: {e}')
        sys.exit(1)

    uvicorn.run(rr_app, host=args.host, port=args.port)


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

    app = subparsers.add_parser('app', help='Launch river-route-app web UI')
    app.add_argument('--host', type=str, default='127.0.0.1', help='Host interface for the app server')
    app.add_argument('--port', type=int, default=8000, help='Port for the app server')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    if args.command == 'route':
        _route(args)
    elif args.command == 'app':
        _launch_app(args)
    else:
        ROUTERS[args.command](args.config).route()
