"""
``ECMWFGribReducedGrid`` reads runoff depths from ECMWF GRIB files on a reduced gaussian grid with eccodes, the
forcing ecmwf_grib. It is a ``BaseGridRunoff`` (bases.py) that overrides only ``read_runoff``, so it yields the same
``GridCellRunoff``, or ``CatchmentRunoffVolumes`` when a file must be aggregated first, routed with the overloads
registered in bases.py. It defines no form of runoff or overload of its own.
"""

import os
from dataclasses import KW_ONLY, dataclass

import eccodes
import numpy as np

from ..types import DatetimeArray, FloatArray, IntArray, PathInput
from .bases import BaseGridRunoff

__all__ = ['ECMWFGribReducedGrid']


@dataclass(eq=False, repr=False)
class ECMWFGribReducedGrid(BaseGridRunoff):
    """
    Runoff from ECMWF GRIB files on a reduced gaussian grid, such as the octahedral O1280 grid of the IFS, with a
    weight table made by ``reduced_grid_weights`` that locates each cell by its ``cell_index``, its position in the
    values of a GRIB message. This reader is specialized to that format; runoff in any other form must be prepared as
    ``catchment`` or ``grid`` forcing.

    GRIB files are read with eccodes one message at a time: each message is decoded whole and only the weight table's
    cells are kept. At the 18,721 cells of the Columbia River basin in a 145 step O1280 forecast that took 0.72 s with
    memory for one message, where xarray and cfgrib took 5.8 s and built the whole (step, cell) array in 4 GB.

    Every message whose shortName is ``var_grid_runoff`` is read as one time step, in the order of the files and of
    their messages, at its validity date and time. Choosing files whose grid, ensemble member, and steps suit the
    weight table and the routing is left to the caller. A GRIB message has no named dimensions, so ``var_cell`` and
    ``var_t`` are not read: each cell is taken by its position in the message's values. IFS forecast runoff
    accumulates from the start of the forecast, so it is routed with ``grid_accumulation_type`` cumulative.
    """

    _: KW_ONLY
    var_cell: str = 'cell'  # name of the grid cell dimension

    @property
    def cell_dimensions(self) -> dict[str, str]:
        return {'cell_index': self.var_cell}

    def read_runoff(self, runoff_data: PathInput | list[PathInput]) -> tuple[FloatArray, DatetimeArray, float]:
        """
        Read runoff at the grid cells the weight table touches from GRIB files, in the order of their messages.

        Args:
            runoff_data: path(s) to GRIB files

        Returns:
            tuple: (runoff as a C-order (time, n_cells) array, validity times, factor converting the depth unit to
                meters)
        """
        runoff_files = [runoff_data] if isinstance(runoff_data, str | os.PathLike) else list(runoff_data)
        messages_per_file = []
        for runoff_file in runoff_files:
            with open(runoff_file, 'rb') as file:
                messages_per_file.append(eccodes.codes_count_in_file(file))
        cell_index = self.cell_indexes['cell_index']
        runoff = np.empty((sum(messages_per_file), cell_index.shape[0]), dtype=np.float32)
        dates = np.empty(runoff.shape[0], dtype='datetime64[s]')
        units = set()
        n_steps = 0
        for runoff_file, n_messages in zip(runoff_files, messages_per_file, strict=True):
            with open(runoff_file, 'rb') as file:
                for _ in range(n_messages):
                    handle = eccodes.codes_grib_new_from_file(file)
                    try:
                        if eccodes.codes_get(handle, 'shortName') != self.var_grid_runoff:
                            continue
                        dates[n_steps] = self._read_message(handle, cell_index, runoff[n_steps], runoff_file)
                        units.add(eccodes.codes_get(handle, 'units'))
                        n_steps += 1
                    finally:
                        eccodes.codes_release(handle)
        if n_steps == 0:
            raise ValueError(f'no GRIB message in {runoff_data} has shortName {self.var_grid_runoff}')
        if len(units) != 1:
            raise ValueError(f'the {self.var_grid_runoff} messages of {runoff_data} have units {sorted(units)}')
        conversion_factor = self._get_conversion_factor(self.runoff_depth_unit or units.pop())
        return runoff[:n_steps], dates[:n_steps], conversion_factor

    def _read_message(
        self, handle: int, cell_index: IntArray, row: FloatArray, runoff_file: PathInput
    ) -> np.datetime64:
        """Write a message's runoff at the weight table's cells into row, missing values as NaN, and return its
        validity time."""
        values = eccodes.codes_get_array(handle, 'values', np.float32)
        if cell_index.size and cell_index.max() >= values.shape[0]:
            raise ValueError(f'the weight table references cells beyond the {values.shape[0]} in {runoff_file}')
        np.take(values, cell_index, out=row)
        if eccodes.codes_get(handle, 'bitmapPresent'):
            row[row == np.float32(eccodes.codes_get(handle, 'missingValue'))] = np.nan
        date, time = eccodes.codes_get(handle, 'validityDate'), eccodes.codes_get(handle, 'validityTime')
        return np.datetime64(
            f'{date // 10000:04d}-{date // 100 % 100:02d}-{date % 100:02d}T{time // 100:02d}:{time % 100:02d}', 's'
        )
