"""
``GridRunoff`` reads runoff depths on a grid with x and y dimensions with xarray, the forcing grid. It is a
``BaseGridRunoffForcing`` (bases.py), which yields a ``GridCellRunoff`` whose cells routing reads as it routes each
river, or ``CatchmentRunoffVolumes`` when a file must be aggregated first. Both are routed with the overloads registered
in bases.py, so this module defines none of its own.
"""

from dataclasses import KW_ONLY, dataclass
from typing import Literal

from .bases import BaseGridRunoffForcing

__all__ = ['GridRunoff']


@dataclass(eq=False, repr=False)
class GridRunoff(BaseGridRunoffForcing):
    """
    Runoff on a grid with x and y dimensions, such as a regular latitude and longitude grid, whose weight table
    locates each cell with ``x_index`` and ``y_index``.
    """

    _: KW_ONLY
    var_x: str = 'x'  # name of the grid x coordinate variable
    var_y: str = 'y'  # name of the grid y coordinate variable
    accumulation_type: Literal['incremental', 'cumulative'] = 'incremental'

    @property
    def cell_dimensions(self) -> dict[str, str]:
        """The weight table columns x_index and y_index, and the runoff dimensions each indexes."""
        return {'x_index': self.var_x, 'y_index': self.var_y}

    @property
    def cumulative(self) -> bool:
        """Whether the runoff values are totals accumulated from the start of each file, as accumulation_type says."""
        return self.accumulation_type == 'cumulative'
