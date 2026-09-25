from .CatchmentRunoff import CatchmentRunoff
from .GaussianGridRunoff import CellRunoff, GaussianGridRunoff
from .ReducedGaussianGridRunoff import ReducedGaussianGridRunoff
from .Runoff import Runoff
from .weights import (
    cell_xy_from_regular_grid,
    compute_voronoi_catchment_intersects,
    grid_weights,
    voronoi_diagram_from_regular_xy,
)

# the Runoff class that reads each Configs runoff_type
RUNOFF_CLASS_FOR_RUNOFF_TYPE: dict[str, type[Runoff]] = {
    'catchment': CatchmentRunoff,
    'gaussian_grid': GaussianGridRunoff,
    'reduced_gaussian_grid': ReducedGaussianGridRunoff,
}

__all__ = [
    'RUNOFF_CLASS_FOR_RUNOFF_TYPE',
    'CellRunoff',
    'Runoff',
    'CatchmentRunoff',
    'GaussianGridRunoff',
    'ReducedGaussianGridRunoff',
    'cell_xy_from_regular_grid',
    'voronoi_diagram_from_regular_xy',
    'compute_voronoi_catchment_intersects',
    'grid_weights',
]
