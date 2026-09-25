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

# the Runoff class that reads each Configs runoff_type, which Configs requires whenever forcing is runoff
RUNOFF_CLASS_FOR_RUNOFF_TYPE: dict[str | None, type[Runoff]] = {
    'catchment': CatchmentRunoff,
    'gaussian_grid': GaussianGridRunoff,
    'reduced_gaussian_grid': ReducedGaussianGridRunoff,
}

__all__ = [
    'RUNOFF_CLASS_FOR_RUNOFF_TYPE',
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
