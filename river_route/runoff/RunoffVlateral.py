from dataclasses import dataclass, fields
from typing import Self

import numpy as np
import xarray as xr

from ..configs import Configs
from ..types import PathList, VlateralGenerator
from .Runoff import Runoff

__all__ = ['RunoffVlateral']


@dataclass(kw_only=True, eq=False)
class RunoffVlateral(Runoff):
    """
    Reads lateral inflow (vlateral) that is already prepared as volumes in netCDF files, for routing.

    The names it reads the inflow and its time steps by are fields, so a file that names them differently can be
    routed without rewriting it. ``Runoff.to_netcdf`` always writes the defaults, which are the names documented for
    ``vlateral_files``, so files this package wrote are read with a plain ``RunoffVlateral()``.

    Every option is named the same here as it is on ``Configs``, so ``RunoffVlateral.from_configs`` reads each one off
    a ``Configs`` by its own name, which is how a ``Router`` builds one.
    """

    var_vlateral: str = 'vlateral'  # name of the lateral inflow variable in the vlateral files
    var_t: str = 'time'  # name of the time coordinate variable

    @property
    def var_runoff(self) -> str:
        """The ``Runoff`` interface name for ``var_vlateral``, the inflow variable of the vlateral files."""
        return self.var_vlateral

    def reader(self, vlateral_files: PathList) -> VlateralGenerator:
        """
        Read lateral inflow straight from netCDF files, one at a time.

        Args:
            vlateral_files: netCDF files each holding a ``vlateral`` variable with dimensions (time, river)

        Yields:
            tuple: (dates, vlateral, source_file) per input file, with NaN volumes replaced by zero
        """
        for lateral_file in vlateral_files:
            with xr.open_dataset(lateral_file) as ds:
                dates = ds[self.var_t].values.astype('datetime64[s]')
                array = ds[self.var_vlateral].values.astype(np.float32, copy=False)
            np.nan_to_num(array, copy=False, nan=0.0)  # the routing kernels never see NaN
            yield dates, array, lateral_file

    @classmethod
    def from_configs(cls, configs: Configs) -> Self:
        """
        Build a RunoffVlateral from the vlateral options on a ``Configs``. Each option is read off the Configs by the
        field's own name, so a field with no matching option raises AttributeError rather than silently keeping its
        default. The configs are not validated for runoff, since ``validate_runoff`` requires the
        ``grid_weights_file`` that routing from ``vlateral_files`` does not use.
        """
        if not isinstance(configs, Configs):
            raise TypeError(
                f'from_configs takes a Configs, got {type(configs).__name__}. '
                f'Use Configs(...) or Configs.from_file(path).'
            )
        return cls(**{f.name: getattr(configs, f.name) for f in fields(cls) if f.init})
