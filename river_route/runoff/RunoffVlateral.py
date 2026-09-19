import numpy as np
import xarray as xr

from ..types import PathList, VlateralGenerator
from .Runoff import Runoff

__all__ = ['RunoffVlateral']


class RunoffVlateral(Runoff):
    """
    Reads lateral inflow (vlateral) that is already prepared as volumes in netCDF files, for routing.
    """

    def __repr__(self) -> str:
        return f'{type(self).__name__}()'

    def reader(self, vlateral_files: PathList) -> VlateralGenerator:
        """
        Read lateral inflow straight from netCDF files, one at a time.

        Args:
            vlateral_files: netCDF files each holding a ``vlateral`` variable with dimensions (time, river)

        Yields:
            tuple: (dates, vlateral, source_file) per input file
        """
        for lateral_file in vlateral_files:
            with xr.open_dataset(lateral_file) as ds:
                dates = ds['time'].values.astype('datetime64[s]')
                array = ds['vlateral'].values.astype(np.float32, copy=False)
                yield dates, array, lateral_file
