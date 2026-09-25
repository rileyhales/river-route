import logging
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Self

import numpy as np
import pandas as pd
import xarray as xr

from .._metadata import __version__
from ..types import DatetimeArray, FloatArray, IntArray, PathInput, PathList, RunoffGenerator

if TYPE_CHECKING:
    from ..configs.Configs import Configs
    from ..network.Network import Network

__all__ = ['Runoff', 'CATCHMENT_RUNOFF', 'CATCHMENT_AREA', 'VOLUME_UNITS']

# the fixed schema of a catchment runoff file: one depth or volume series per catchment and the area of each catchment
CATCHMENT_RUNOFF = 'catchment_runoff'
CATCHMENT_AREA = 'catchment_area'
VOLUME_UNITS = ('m3', 'm^3', 'm³')  # units attributes that mark catchment runoff as volumes, anything else is a depth

logger = logging.getLogger(__name__)


class Runoff(ABC):
    """
    Base class for the sources of catchment runoff, the runoff volume of each catchment before it is transformed into
    lateral inflow to its river. Each subclass is built from a Configs with ``from_configs`` and provides a
    ``generator`` that yields its runoff in the form routing reads, and every subclass writes catchment runoff to disk
    with ``to_netcdf`` in the one format ``CatchmentRunoff`` reads.
    """

    as_volumes: bool = True  # whether the arrays are volumes (m³) rather than depths (m)

    @classmethod
    @abstractmethod
    def from_configs(cls, configs: Configs) -> Self:
        """Build this Runoff from the options on a Configs, as a Router does when it is not given one."""

    @abstractmethod
    def generator(self, runoff_files: PathList) -> RunoffGenerator:
        """
        Yield one (dates, runoff, source_file) tuple per input, where runoff is what routing reads: its
        CatchmentRunoffVolumes, or a GridCellRunoff that routing aggregates as it routes. Each checks its own arrays
        with ``check`` and gives the runoff of its first steps with ``first_steps``.
        """

    def distribute(self, network: Network) -> None:
        """
        Match the runoff to a stabilized network, whose synthetic sub-reaches each need a share of their parent
        river's runoff. Classes that can split their runoff override this; the others raise when the network has rows
        they have no runoff for.
        """
        if network.synthetic is not None:
            raise NotImplementedError(f'{type(self).__name__} cannot be routed on a stabilized network yet')
        return

    def to_netcdf(
        self,
        path: PathInput,
        dates: DatetimeArray,
        catchment_runoff: FloatArray,
        river_ids: IntArray,
        catchment_area: FloatArray,
    ) -> None:
        """
        Write a catchment runoff array to a netCDF file that ``CatchmentRunoff`` reads, routed as ``runoff_files``
        with ``runoff_type`` catchment.

        Args:
            path: netCDF file to write
            dates: (time,) datetime64 values of the steps
            catchment_runoff: (n_rivers, time) depths (m) or volumes (m³), rows ordered like ``river_ids``
            river_ids: (n_rivers,) river id of each column
            catchment_area: (n_rivers,) area of each catchment in m², the factor between depths and volumes
        """
        self._catchment_runoff_dataset(dates, catchment_runoff, river_ids, catchment_area).to_netcdf(path)
        return

    def _catchment_runoff_dataset(
        self, dates: DatetimeArray, catchment_runoff: FloatArray, river_ids: IntArray, catchment_area: FloatArray
    ) -> xr.Dataset:
        """Build the catchment runoff dataset with dimensions (``river_id``, ``time``) that ``to_netcdf`` writes."""
        units = 'm3' if self.as_volumes else 'm'
        kind = 'volumes' if self.as_volumes else 'depths'
        start_date = pd.Timestamp(dates[0]).strftime('%Y%m%d%H')
        end_date = pd.Timestamp(dates[-1]).strftime('%Y%m%d%H')
        timestep = int((dates[1] - dates[0]) / np.timedelta64(1, 's')) if len(dates) > 1 else 0
        return xr.Dataset(
            {
                CATCHMENT_RUNOFF: xr.DataArray(
                    catchment_runoff,
                    dims=('river_id', 'time'),
                    attrs={
                        'long_name': f'Incremental catchment runoff {kind}',
                        'units': units,
                        'cell_measures': f'area: {CATCHMENT_AREA}',
                    },
                ),
                CATCHMENT_AREA: xr.DataArray(
                    np.asarray(catchment_area).astype(np.float32, copy=False),
                    dims=('river_id',),
                    attrs={'long_name': 'catchment area of each river', 'units': 'm2'},
                ),
            },
            coords={
                'river_id': xr.DataArray(
                    np.asarray(river_ids).astype(np.int32, copy=False),
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
                'title': f'Incremental catchment runoff {kind}',
                'description': f'Incremental catchment runoff ({units}) for each river',
                'source': f'river-route v{__version__}',
                'history': f'Created on {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}',
                'suggested_file_name': f'catchment_runoff_{start_date}_{end_date}.nc',
            },
        )

    @staticmethod
    def _cumulative_to_incremental(df: pd.DataFrame) -> pd.DataFrame:
        values = df.to_numpy()
        return pd.DataFrame(np.vstack([values[:1], np.diff(values, axis=0)]), index=df.index, columns=df.columns)

    @staticmethod
    def _get_conversion_factor(unit: str | None) -> float:
        if unit is None:
            logger.warning('No units attribute found. Assuming meters')
            return 1
        if unit in ('m', 'meters', 'kg m-2'):
            return 1
        if unit in ('mm', 'millimeters'):
            return 0.001
        raise ValueError(f'Unknown units: {unit}')
