from typing import Self

import numpy as np
import xarray as xr

from ..configs import Configs
from ..types import PathList, RunoffGenerator
from .Runoff import CATCHMENT_AREA, CATCHMENT_RUNOFF, VOLUME_UNITS, Runoff

__all__ = ['CatchmentRunoff']


class CatchmentRunoff(Runoff):
    """
    A degenerate Runoff aggregator class which only allows for iterating over files already converted to
    catchment level volumes or depths.
    """

    as_volumes = True  # catchment runoff is always read as volumes, the quantity routing uses

    def generator(self, runoff_files: PathList) -> RunoffGenerator:
        """
        Read catchment runoff straight from netCDF files, one at a time.

        Args:
            runoff_files: netCDF files each holding ``catchment_runoff`` with dimensions (river_id, time)
                and, when it holds depths, ``catchment_area`` with dimension river_id

        Yields:
            tuple: (dates, catchment_runoff, source_file) per input file, as a C-order (n_rivers, time) array of
                volumes (m³) with NaN replaced by zero
        """
        for runoff_file in runoff_files:
            with xr.open_dataset(runoff_file) as ds:
                dates = ds['time'].values.astype('datetime64[s]')
                runoff = ds[CATCHMENT_RUNOFF]
                if runoff.dims != ('river_id', 'time'):
                    raise ValueError(
                        f'{CATCHMENT_RUNOFF} in {runoff_file} has dimensions {runoff.dims}, expected (river_id, time)'
                    )
                units = runoff.attrs.get('units')
                if units is None:
                    raise ValueError(f'{CATCHMENT_RUNOFF} in {runoff_file} has no units attribute, m3 or a depth unit')
                array = runoff.values.astype(np.float32, copy=False)
                if units not in VOLUME_UNITS:
                    area = ds[CATCHMENT_AREA].values.astype(np.float32) * np.float32(self._get_conversion_factor(units))
                    array *= area[:, np.newaxis]
            np.nan_to_num(array, copy=False, nan=0.0)  # the routing kernels never see NaN
            yield dates, array, runoff_file

    @classmethod
    def from_configs(cls, configs: Configs) -> Self:
        """
        Build a CatchmentRunoff. The catchment runoff file schema is fixed, so no option of the ``Configs`` is read.
        """
        if not isinstance(configs, Configs):
            raise TypeError(
                f'from_configs takes a Configs, got {type(configs).__name__}. '
                f'Use Configs(...) or Configs.from_json(path).'
            )
        return cls()
