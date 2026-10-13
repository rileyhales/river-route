from .bases import BaseGridRunoffForcing, CatchmentRunoffVolumes, GridCellRunoff, RunoffForcing
from .CatchmentRunoff import CatchmentRunoffForcing
from .ECMWFGribReducedGrid import ECMWFGribReducedGrid, ReducedGaussianGrid
from .GridRunoff import GridRunoff
from .weights import (
    cell_xy_from_regular_grid,
    compute_voronoi_catchment_intersects,
    grid_weights,
    reduced_grid_weights,
    voronoi_diagram_from_regular_xy,
)

# the RunoffForcing class that reads the forcing_files of each Configs forcing but channel, which routes no runoff
RUNOFF_CLASS_FOR_FORCING: dict[str, type[RunoffForcing]] = {
    'catchment': CatchmentRunoffForcing,
    'grid': GridRunoff,
    'ecmwf_grib': ECMWFGribReducedGrid,
}

__all__ = [
    'RUNOFF_CLASS_FOR_FORCING',
    'CatchmentRunoffVolumes',
    'GridCellRunoff',
    'RunoffForcing',
    'BaseGridRunoffForcing',
    'CatchmentRunoffForcing',
    'GridRunoff',
    'ECMWFGribReducedGrid',
    'cell_xy_from_regular_grid',
    'voronoi_diagram_from_regular_xy',
    'compute_voronoi_catchment_intersects',
    'grid_weights',
    'ReducedGaussianGrid',
    'reduced_grid_weights',
]
