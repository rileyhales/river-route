from typing import NamedTuple, Self

import numpy as np
import xarray as xr
from numba.extending import overload

from ..configs import Configs
from ..router._routing_passes import (
    count_rivers_prepared_together,
    get_river_catchment_runoff,
    is_argument_type,
    prepare_runoff_of_rivers,
)
from ..types import FloatArray, PathList, RunoffGenerator
from .Runoff import CATCHMENT_AREA, CATCHMENT_RUNOFF, VOLUME_UNITS, Runoff

__all__ = ['CatchmentRunoffVolumes', 'CatchmentRunoff']


class CatchmentRunoffVolumes(NamedTuple):
    """Each river's catchment runoff volume (m³) in every runoff step, as a C-order (river, time) array whose rows
    routing reads in place."""

    runoff: FloatArray  # (n_rivers, n_steps) float32 volumes (m³)

    def check(self, n_rivers: int, n_steps: int) -> None:
        """Raise ValueError unless it is (n_rivers, n_steps) with each river's row contiguous."""
        if self.runoff.ndim != 2 or self.runoff.strides[1] != self.runoff.itemsize:
            raise ValueError('catchment runoff must be a (river, time) array with each river row contiguous')
        if self.runoff.shape != (n_rivers, n_steps):
            raise ValueError(f'catchment runoff has shape {self.runoff.shape}, expected {(n_rivers, n_steps)}')

    def first_steps(self, n_steps: int) -> CatchmentRunoffVolumes:
        """The runoff of the first n_steps steps, a view whose rows stay contiguous."""
        return CatchmentRunoffVolumes(self.runoff[:, :n_steps])


@overload(count_rivers_prepared_together)
def _catchment_runoff_is_read_in_place(runoff):
    if is_argument_type(runoff, CatchmentRunoffVolumes):
        return lambda runoff: 0


@overload(prepare_runoff_of_rivers)
def _catchment_runoff_needs_no_preparing(runoff, first_river, stop_river, scratch):
    if is_argument_type(runoff, CatchmentRunoffVolumes):
        return lambda runoff, first_river, stop_river, scratch: None


@overload(get_river_catchment_runoff)
def _catchment_runoff_of_river(runoff, r, first_river, scratch):
    if is_argument_type(runoff, CatchmentRunoffVolumes):
        return lambda runoff, r, first_river, scratch: runoff.runoff[r]


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
            tuple: (dates, catchment_runoff, source_file) per input file, the runoff as CatchmentRunoffVolumes:
                C-order (n_rivers, time) volumes (m³) with NaN replaced by zero
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
            yield dates, CatchmentRunoffVolumes(array), runoff_file

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
