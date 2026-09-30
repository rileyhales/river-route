from .bases import BaseGridRunoff, CatchmentRunoffVolumes, GridCellRunoff, Runoff
from .CatchmentRunoff import CatchmentRunoff
from .ECMWFGribReducedGrid import ECMWFGribReducedGrid, ReducedGaussianGrid
from .GridRunoff import GridRunoff
from .weights import (
    cell_xy_from_regular_grid,
    compute_voronoi_catchment_intersects,
    grid_weights,
    reduced_grid_weights,
    voronoi_diagram_from_regular_xy,
)

# the Runoff class that reads the runoff_files of each Configs forcing but channel, which routes no runoff
RUNOFF_CLASS_FOR_FORCING: dict[str, type[Runoff]] = {
    'catchment': CatchmentRunoff,
    'grid': GridRunoff,
    'ecmwf_grib': ECMWFGribReducedGrid,
}

__all__ = [
    'RUNOFF_CLASS_FOR_FORCING',
    'CatchmentRunoffVolumes',
    'GridCellRunoff',
    'Runoff',
    'BaseGridRunoff',
    'CatchmentRunoff',
    'GridRunoff',
    'ECMWFGribReducedGrid',
    'cell_xy_from_regular_grid',
    'voronoi_diagram_from_regular_xy',
    'compute_voronoi_catchment_intersects',
    'grid_weights',
    'ReducedGaussianGrid',
    'reduced_grid_weights',
]
