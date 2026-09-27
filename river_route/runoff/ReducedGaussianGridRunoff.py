from typing import Self

from ..configs import Configs
from ..types import PathInput, PathList, RunoffGenerator
from .Runoff import Runoff

__all__ = ['ReducedGaussianGridRunoff']


class ReducedGaussianGridRunoff(Runoff):
    """
    Placeholder for aggregating runoff on a reduced gaussian grid to catchments. A reduced grid indexes its cells along
    one dimension, named by the ``var_cell`` config, instead of the x and y dimensions of a ``GaussianGridRunoff``.
    Not implemented yet.
    """

    def __init__(self, *args, **kwargs) -> None:
        raise NotImplementedError('ReducedGaussianGridRunoff is not implemented yet')

    def generator(self, runoff_files: PathList) -> RunoffGenerator:
        raise NotImplementedError('ReducedGaussianGridRunoff is not implemented yet')

    def aggregate_to_file(self, runoff_data: PathInput | list[PathInput], path: PathInput) -> None:
        raise NotImplementedError('ReducedGaussianGridRunoff is not implemented yet')

    @classmethod
    def from_configs(cls, configs: Configs) -> Self:
        raise NotImplementedError('ReducedGaussianGridRunoff is not implemented yet')
