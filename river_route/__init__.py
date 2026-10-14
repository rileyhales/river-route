import sys

if sys.platform not in ('darwin', 'linux'):  # its C kernel is compiled and its GRIB library installed on these only
    raise ImportError(f'river-route runs on macOS and Linux, not {sys.platform}')

from . import forcings, kernels, metrics, network, router
from ._metadata import __author__, __url__, __version__
from .configs import Configs
from .forcings import CatchmentRunoffForcing, ECMWFGribReducedGrid, GridRunoff
from .network import Network
from .router import Router

__all__ = [
    # configuration, the object every other class is built from
    'Configs',
    # the river network: topology, parameters, partitioning, coefficient signs, sub-reaches and subcycles
    'Network',
    # routing
    'Router',
    # forcings: the runoff classes that read and prepare forcing_files
    'CatchmentRunoffForcing',
    'GridRunoff',
    'ECMWFGribReducedGrid',
    # modules
    'network',
    'router',
    'kernels',
    'forcings',
    'metrics',
    # metadata
    '__version__',
    '__author__',
    '__url__',
]
