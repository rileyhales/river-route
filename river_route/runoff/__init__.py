from .CatchmentRunoff import CatchmentRunoff, CatchmentRunoffVolumes
from .GaussianGridRunoff import GaussianGridRunoff, GridCellRunoff
from .ReducedGaussianGridRunoff import ReducedGaussianGridRunoff
from .Runoff import Runoff
from .weights import (
    cell_xy_from_regular_grid,
    compute_voronoi_catchment_intersects,
    grid_weights,
    voronoi_diagram_from_regular_xy,
)

# the Runoff class that reads the runoff_files of each Configs forcing but channel, which routes no runoff
RUNOFF_CLASS_FOR_FORCING: dict[str, type[Runoff]] = {
    'catchment': CatchmentRunoff,
    'gaussian_grid': GaussianGridRunoff,
    'reduced_gaussian_grid': ReducedGaussianGridRunoff,
}

__all__ = [
    'RUNOFF_CLASS_FOR_FORCING',
    'CatchmentRunoffVolumes',
    'GridCellRunoff',
    'Runoff',
    'CatchmentRunoff',
    'GaussianGridRunoff',
    'ReducedGaussianGridRunoff',
    'cell_xy_from_regular_grid',
    'voronoi_diagram_from_regular_xy',
    'compute_voronoi_catchment_intersects',
    'grid_weights',
]
