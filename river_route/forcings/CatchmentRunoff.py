from typing import Self

import numpy as np
import xarray as xr

from .bases import CATCHMENT_AREA, CATCHMENT_RUNOFF, VOLUME_UNITS, CatchmentRunoffVolumes, RunoffForcing
from ..configs import Configs
from ..types import PathList, RunoffGenerator

__all__ = ['CatchmentRunoffForcing']


class CatchmentRunoffForcing(RunoffForcing):
    """
    Reads runoff already aggregated to catchments and yields each as ``CatchmentRunoffVolumes`` (bases.py),
    whose rows routing reads in place, a cell per river. It defines no form of runoff of its own.
    A degenerate Runoff aggregator class which only allows for iterating over files already converted to
    catchment level volumes or depths.
    """

    def generator(self, forcing_files: PathList) -> RunoffGenerator:
        """
        Read catchment runoff straight from netCDF files, one at a time.

        Args:
            forcing_files: netCDF files each holding ``catchment_runoff`` with dimensions (riverId, time)
                and, when it holds depths, ``catchment_area`` with dimension riverId

        Yields:
            tuple: (dates, catchment_runoff, source_file) per input file, the runoff as CatchmentRunoffVolumes:
                C-order (n_rivers, time) volumes (m³) with NaN replaced by zero
        """
        for forcing_file in forcing_files:
            with xr.open_dataset(forcing_file) as ds:
                dates = ds['time'].values.astype('datetime64[s]')
                runoff = ds[CATCHMENT_RUNOFF]
                if runoff.dims != ('riverId', 'time'):
                    raise ValueError(
                        f'{CATCHMENT_RUNOFF} in {forcing_file} has dimensions {runoff.dims}, expected (riverId, time)'
                    )
                if 'riverId' not in ds.variables:
                    raise ValueError(f'{forcing_file} has no riverId coordinate naming the river of each row')
                units = runoff.attrs.get('units')
                if units is None:
                    raise ValueError(f'{CATCHMENT_RUNOFF} in {forcing_file} has no units attribute, m3 or a depth unit')
                river_ids = ds['riverId'].values
                array = runoff.values.astype(np.float32, copy=False)
                if units not in VOLUME_UNITS:
                    area = ds[CATCHMENT_AREA].values.astype(np.float32) * np.float32(self._get_conversion_factor(units))
                    array *= area[:, np.newaxis]
            np.nan_to_num(array, copy=False, nan=0.0)  # the routing kernels never see NaN
            yield dates, CatchmentRunoffVolumes(array, river_ids), forcing_file

    @classmethod
    def from_configs(cls, configs: Configs) -> Self:
        """
        Build a CatchmentRunoffForcing. The catchment runoff file schema is fixed, so no option of the ``Configs`` is
        read.
        """
        if not isinstance(configs, Configs):
            raise TypeError(
                f'from_configs takes a Configs, got {type(configs).__name__}. '
                f'Use Configs(...) or Configs.from_json(path).'
            )
        return cls()
