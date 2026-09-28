"""
``GridRunoff`` reads runoff depths on a grid with x and y dimensions with xarray, the forcing grid. It is a
``BaseGridRunoff`` (bases.py), which yields a ``GridCellRunoff`` whose cells routing reads as it routes each river, or
``CatchmentRunoffVolumes`` when a file must be aggregated first. Both are routed with the overloads registered in
bases.py, so this module defines none of its own.
"""

from dataclasses import KW_ONLY, dataclass

from .bases import BaseGridRunoff

__all__ = ['GridRunoff']


@dataclass(eq=False, repr=False)
class GridRunoff(BaseGridRunoff):
    """
    Runoff on a grid with x and y dimensions, such as a regular latitude and longitude grid, whose weight table
    locates each cell with ``x_index`` and ``y_index``.
    """

    _: KW_ONLY
    var_x: str = 'x'  # name of the grid x coordinate variable
    var_y: str = 'y'  # name of the grid y coordinate variable

    @property
    def cell_dimensions(self) -> dict[str, str]:
        return {'x_index': self.var_x, 'y_index': self.var_y}
