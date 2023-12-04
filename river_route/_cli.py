import argparse
import json
import os

import river_route as rr


def _json_to_python_type(json_type):
    if json_type == 'number':
        return int
    elif json_type == 'string':
        return str
    elif json_type == 'boolean':
        return bool
    else:
        raise ValueError(f'Unknown json type: {json_type}')


def _add_args_from_file(parser, params):
    for k, v in params.items():
        v['type'] = _json_to_python_type(v['type'])
        kwargs = {
            'help': v['help'],
            'type': v['type'],
            'action': 'store' if v['type'] != bool else 'store_true',
        }
        if v['type'] is bool:
            kwargs.pop('type')
        parser.add_argument(f'--{k}', **kwargs, )
    return


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
        required=True,
    )

    description_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config_files', 'descriptions.json')
    with open(description_file, 'r') as f:
        _add_args_from_file(parser, json.load(f))

    args = parser.parse_args()

    if not args.config:
        return parser.print_help()

    kwargs = vars(args)
    kwargs.pop('config')
    rr.Muskingum(args.config, **kwargs).route()
    return
