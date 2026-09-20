"""
Tests for the runoff source the configs name, run on the synthetic network from conftest.
"""

import pytest
from conftest import SyntheticNetwork

import river_route as rr


def test_routing_requires_a_runoff_source(network: SyntheticNetwork):
    configs = rr.Configs(
        forcing='vlateral',
        params_file=str(network.params_file),
        discharge_files=[network.path('q.nc')],
        log=False,
        progress_bar=False,
    )
    with pytest.raises(ValueError, match='Provide vlateral_files or grid_runoff_files'):
        rr.Router(configs).route()
