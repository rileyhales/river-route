from .TransformMuskingum import TransformMuskingum
from ..types import FloatArray

__all__ = ['RapidMuskingum', ]


class RapidMuskingum(TransformMuskingum):
    """
    Muskingum router that accepts pre-aggregated catchment volumes as lateral inflow.

    Runoff volumes (m³) are divided by the runoff time step to convert to flow rates (m³/s)
    before being routed through the channel network using the Muskingum method. This matches
    the RAPID model convention where overland flow is assumed to enter the channel uniformly
    over the runoff time step.

    Required configs:
    - params_file: path to routing parameters parquet file.
    - discharge_dir: directory where output netCDF files are written (named after input files).
    Lateral input — one of:
    - catchment_runoff_files: list of paths to netCDF files containing per-catchment volumes (m³).
    - runoff_grid_files + grid_weights_file: gridded runoff depth inputs remapped to catchments.
    """

    @property
    def _catchment_runoff_as_volume(self) -> bool:
        return True  # True -> means as volume

    def transform_runoff(self, r_t: FloatArray) -> FloatArray:
        """Convert a runoff volume (m³) to the lateral term for the routing equation."""
        return r_t / self.dt_runoff * (self.c1 + self.c2)  # note c1 + c2 = c4 in muskingum cunge
