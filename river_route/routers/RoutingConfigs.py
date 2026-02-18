import json
import logging
import os
import random
import sys
from pathlib import Path
from typing import Any

import yaml

from .typing import ConfigDict, PathInput
from ..__metadata__ import __version__

__all__ = ['RoutingConfigs', ]


class RoutingConfigs:
    conf: ConfigDict
    logger: logging.Logger

    def __init__(self, config_file: PathInput | None = None, **kwargs: Any) -> None:
        self.logger = logging.getLogger(f'river_route.{''.join([str(random.randint(0, 9)) for _ in range(8)])}')
        self._set_configs(config_file, **kwargs)
        return

    def __repr__(self):
        messages = ['Configs:', ] + [f'\t{k}: {v}' for k, v in self.conf.items()]
        return '\n'.join(messages)

    def _set_configs(self, config_file: PathInput | None, **kwargs: Any) -> None:
        if config_file is None or config_file == '':
            self.conf = {}
        elif str(config_file).endswith('.json'):
            with open(config_file, 'r') as f:
                self.conf = json.load(f)
        elif str(config_file).endswith('.yml') or str(config_file).endswith('.yaml'):
            with open(config_file, 'r') as f:
                self.conf = yaml.load(f, Loader=yaml.FullLoader)
        else:
            raise RuntimeError('Unrecognized simulation config file type. Must be .json or .yaml')

        self.conf.update(kwargs)
        self.conf['river-route-version'] = __version__
        self.conf['log'] = bool(self.conf.get('log', True))
        self.conf['progress_bar'] = self.conf.get('progress_bar', self.conf['log'])
        self.conf['log_level'] = self.conf.get('log_level', 'INFO')

        # compute and routing options (time is validated separately at compute step)
        self.conf['input_type'] = self.conf.get('input_type', 'sequential')
        self.conf['runoff_type'] = self.conf.get('runoff_type', 'incremental')

        # expected variable names in input/output files
        self.conf['var_river_id'] = self.conf.get('var_river_id', 'river_id')
        self.conf['var_discharge'] = self.conf.get('var_discharge', 'Q')
        self.conf['var_x'] = self.conf.get('var_x', 'x')
        self.conf['var_y'] = self.conf.get('var_y', 'y')
        self.conf['var_t'] = self.conf.get('var_t', 'time')
        self.conf['var_catchment_volume'] = self.conf.get('var_catchment_volume', 'volume')
        self.conf['var_runoff_depth'] = self.conf.get('var_runoff_depth', 'ro')

        # configure logging
        self.logger.disabled = not self.conf.get('log', True)
        self.logger.setLevel(self.conf.get('log_level', 'INFO'))
        log_destination = self.conf.get('log_stream', 'stdout')
        log_format = self.conf.get('log_format', '%(asctime)s - %(levelname)s - %(message)s')
        if log_destination == 'stdout':
            self.logger.addHandler(logging.StreamHandler(sys.stdout))
        elif isinstance(log_destination, str):
            self.logger.addHandler(logging.FileHandler(log_destination))
        self.logger.handlers[0].setFormatter(logging.Formatter(log_format))
        self.logger.debug('Logger initialized')
        return

    def _validate_configs(self) -> None:
        self.logger.debug('Validating configs file')
        if not os.path.exists(self.conf['routing_params_file']):
            raise FileNotFoundError('Routing params file not found')

        # check for valid options
        if self.conf['input_type'] not in ['sequential', 'ensemble']:
            raise ValueError('Input type not recognized')
        if self.conf['runoff_type'] not in ['incremental', 'cumulative']:
            raise ValueError('Runoff type not recognized')

        # fill the files with empty lists so checks are consistent, we'll delete later
        self.conf['catchment_volumes_files'] = self.conf.get('catchment_volumes_files', [])
        self.conf['runoff_depths_files'] = self.conf.get('runoff_depths_files', [])
        self.conf['discharge_files'] = self.conf.get('discharge_files', [])

        # if only 1 file given as str, wrap in a list so future code that anticipates iterables behave properly
        if isinstance(self.conf.get('catchment_volumes_files', []), (str, Path)):
            self.conf['catchment_volumes_files'] = [self.conf['catchment_volumes_files'], ]
        if isinstance(self.conf.get('runoff_depths_files', []), (str, Path)):
            self.conf['runoff_depths_files'] = [self.conf['runoff_depths_files'], ]
        if isinstance(self.conf.get('discharge_files', []), (str, Path)):
            self.conf['discharge_files'] = [self.conf['discharge_files'], ]
        # provided volumes or runoffs must exist and output directory must exist
        n_files = len(self.conf['catchment_volumes_files']) + len(self.conf['runoff_depths_files'])
        if self.conf['catchment_volumes_files'] and self.conf['runoff_depths_files']:
            raise ValueError('Provide either catchment volumes files or runoff depths files, not both')
        if not len(self.conf['discharge_files']) == n_files:
            raise ValueError('Number of discharge files must match the number of input files (volumes or depths)')
        for file in self.conf['catchment_volumes_files']:
            if not os.path.exists(file):
                raise FileNotFoundError(f'Catchment volumes file not found at: {file}')
        for file in self.conf['runoff_depths_files']:
            if not os.path.exists(file):
                raise FileNotFoundError(f'Runoff depths file not found at: {file}')
        for directory in set(os.path.dirname(file) for file in self.conf['discharge_files']):
            if not os.path.exists(directory):
                raise NotADirectoryError(f'Output file directory not found at: {directory}')

        # if we're using runoff depths, a weight table must be provided.
        if self.conf['runoff_depths_files']:
            if not os.path.exists(self.conf['weight_table_file']):
                raise FileNotFoundError('Weight table file not found')

        # now delete empty arrays
        if not self.conf['catchment_volumes_files']:
            del self.conf['catchment_volumes_files']
        if not self.conf['runoff_depths_files']:
            del self.conf['runoff_depths_files']

        # if the initial state was provided but is falsey then remove it behaves with future checks
        if 'initial_state_file' in self.conf and not self.conf['initial_state_file']:
            del self.conf['initial_state_file']
        # if the initial state file is provided, it must be valid type and exist
        if self.conf.get('initial_state_file', ''):
            if not os.path.exists(self.conf['initial_state_file']):
                raise FileNotFoundError('Initial state file not found')

        # convert all relative paths to absolute paths
        for arg in [k for k in self.conf.keys() if 'file' in k]:
            if isinstance(self.conf[arg], list):
                self.conf[arg] = [os.path.abspath(path) for path in self.conf[arg]]
            elif isinstance(self.conf[arg], (str, Path)):
                self.conf[arg] = os.path.abspath(self.conf[arg])
        return
