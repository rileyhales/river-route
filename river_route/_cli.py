import argparse
import os

import pandas as pd

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
    for param in params:
        param['Type'] = _json_to_python_type(param['Type'])
        kwargs = {
            'help': param['Description'],
            'type': param['Type'],
            'action': 'store' if param['Type'] != bool else 'store_true',
        }
        if param['Type'] is bool:
            kwargs.pop('type')
        parser.add_argument(f'--{param["Parameter Name"]}', **kwargs, )
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
        required=False,
    )

    description_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config_files', 'descriptions.csv')
    _add_args_from_file(parser, pd.read_csv(description_file).to_dict('records'))

    args = parser.parse_args()

    if not args.config:
        return parser.print_help()

    kwargs = vars(args)
    kwargs.pop('config')
    rr.Muskingum(args.config, **kwargs).route()
    return
