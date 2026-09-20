import logging
from abc import ABC, abstractclassmethod, abstractmethod

import numpy as np
import pandas as pd
import xarray as xr

from .._metadata import __version__
from ..configs import Configs
from ..types import DatetimeArray, FloatArray, IntArray, PathInput, VlateralGenerator

__all__ = ['Runoff']

logger = logging.getLogger(__name__)


class Runoff(ABC):
    """
    Base class for the sources lateral inflow (vlateral) is routed from. Each subclass provides a ``reader`` that
    yields its vlateral arrays, and every subclass writes them to disk with ``to_netcdf`` in the one format
    ``RunoffVlateral`` reads.
    """

    var_runoff: str
    var_x: str
    var_y: str
    var_t: str
    runoff_depth_unit: str | None
    cumulative: bool
    force_positive_runoff: bool
    force_uniform_timesteps: bool
    as_volumes: bool = True  # whether the arrays are volumes (m³) rather than depths (m)

    @abstractmethod
    def reader(self, *args, **kwargs) -> VlateralGenerator:
        """
        Yield one (dates, vlateral, source_file) tuple per input.
        """

    @abstractclassmethod
    def from_configs(cls, configs: Configs) -> Runoff:
        """
        Build a ``Runoff`` subclass from a ``Configs`` object.
        """

    def to_netcdf(self, path: PathInput, dates: DatetimeArray, vlateral: FloatArray, river_ids: IntArray) -> None:
        """
        Write a vlateral array to a netCDF file that ``RunoffVlateral`` reads and routes with ``vlateral_files``.

        Args:
            path: netCDF file to write
            dates: (time,) datetime64 values of the steps
            vlateral: (time, n_rivers) array ordered like ``river_ids``
            river_ids: (n_rivers,) river id of each column
        """
        self._vlateral_dataset(dates, vlateral, river_ids).to_netcdf(path)
        return

    def _vlateral_dataset(self, dates: DatetimeArray, vlateral: FloatArray, river_ids: IntArray) -> xr.Dataset:
        """Build the vlateral dataset with dimensions ``time`` and ``river_id`` that ``to_netcdf`` writes."""
        units = 'm3' if self.as_volumes else 'm'
        long_name = 'Incremental vlateral volumes' if self.as_volumes else 'Incremental vlateral depths'
        start_date = pd.Timestamp(dates[0]).strftime('%Y%m%d%H')
        end_date = pd.Timestamp(dates[-1]).strftime('%Y%m%d%H')
        timestep = int((dates[1] - dates[0]) / np.timedelta64(1, 's')) if len(dates) > 1 else 0
        return xr.Dataset(
            {
                'vlateral': xr.DataArray(
                    vlateral, dims=('time', 'river_id'), attrs={'long_name': long_name, 'units': units}
                )
            },
            coords={
                'river_id': xr.DataArray(
                    np.asarray(river_ids).astype(np.int64, copy=False),
                    dims=('river_id',),
                    attrs={'long_name': 'unique ID number for each river'},
                ),
                'time': xr.DataArray(
                    dates,
                    dims=('time',),
                    attrs={'long_name': 'time', 'standard_name': 'time', 'axis': 'T', 'time_step': f'{timestep}'},
                ),
            },
            attrs={
                'title': f'Incremental vlateral {long_name.split()[-1]}',
                'description': f'Incremental vlateral ({units}) for each river',
                'source': f'river-route v{__version__}',
                'history': f'Created on {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}',
                'suggested_file_name': f'vlateral_{start_date}_{end_date}.nc',
            },
        )

    @staticmethod
    def _cumulative_to_incremental(df: pd.DataFrame) -> pd.DataFrame:
        return pd.DataFrame(
            np.vstack([df.values[0, :], np.diff(df.values, axis=0)]), index=df.index, columns=df.columns
        )

    @staticmethod
    def _get_conversion_factor(unit: str | None) -> int | float:
        if unit is None:
            logger.warning('No units attribute found. Assuming meters')
            return 1
        if unit in ('m', 'meters', 'kg m-2'):
            return 1
        elif unit in ('mm', 'millimeters'):
            return 0.001
        else:
            raise ValueError(f'Unknown units: {unit}')
