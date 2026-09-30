"""
``ECMWFGribReducedGrid`` reads runoff depths from ECMWF GRIB files on a reduced gaussian grid with eccodes, the
forcing ecmwf_grib. It is a ``BaseGridRunoff`` (bases.py) that overrides only ``read_runoff``, so it yields the same
``GridCellRunoff``, or ``CatchmentRunoffVolumes`` when a file must be aggregated first, routed with the overloads
registered in bases.py. It defines no form of runoff or overload of its own. ``ReducedGaussianGrid`` is the grid of
those files, read from their metadata, on which ``reduced_grid_weights`` (weights.py) builds a weight table.
"""

import os
from dataclasses import KW_ONLY, dataclass
from typing import NamedTuple, Self

import eccodes
import geopandas as gpd
import numpy as np
import shapely

from ..types import DatetimeArray, FloatArray, IntArray, PathInput
from .bases import BaseGridRunoff

__all__ = ['ECMWFGribReducedGrid', 'ReducedGaussianGrid']

EARTH_RADIUS_M = 6371229.0  # radius of the sphere the ECMWF IFS models the earth as


class ReducedGaussianGrid(NamedTuple):
    """
    A global reduced gaussian grid: 2N rows of cells at the gaussian latitudes, north to south, each row's cells evenly
    spaced in longitude eastward from 0 degrees. Cells are numbered row after row, the order of the values in a GRIB
    message.
    """

    N: int  # rows of cells between a pole and the equator
    cells_per_row: IntArray  # (2N,) number of cells on each row, north to south, the GRIB pl array

    @classmethod
    def from_grib(cls, grib_path: PathInput) -> Self:
        """Read the grid from the metadata of the first message in a GRIB file."""
        with open(grib_path, 'rb') as file:
            handle = eccodes.codes_grib_new_from_file(file)
        if handle is None:
            raise ValueError(f'{grib_path} holds no GRIB messages')
        try:
            if eccodes.codes_get(handle, 'gridType') != 'reduced_gg':
                raise ValueError(f'{grib_path} is on a {eccodes.codes_get(handle, "gridType")} grid, not reduced_gg')
            grid = cls(eccodes.codes_get(handle, 'N'), eccodes.codes_get_array(handle, 'pl'))
            numbered_as_described = (
                eccodes.codes_get(handle, 'global') == 1
                and eccodes.codes_get(handle, 'iScansNegatively') == 0
                and eccodes.codes_get(handle, 'jScansPositively') == 0
                and eccodes.codes_get(handle, 'numberOfDataPoints') == grid.cells_per_row.sum()
            )
        finally:
            eccodes.codes_release(handle)
        if not numbered_as_described:
            raise ValueError(f'{grib_path} is not a global grid numbered north to south and west to east')
        if grid.cells_per_row.shape[0] != 2 * grid.N:
            raise ValueError(f'{grib_path} has {grid.cells_per_row.shape[0]} rows of cells, expected {2 * grid.N}')
        return grid

    @property
    def n_cells(self) -> int:
        """The number of cells of the whole grid, the number of values in a GRIB message on it."""
        return int(self.cells_per_row.sum())

    @property
    def spacing_km(self) -> float:
        """Side of a square with the mean cell area, the resolution ECMWF quotes, e.g. 9 km for O1280."""
        return float(np.sqrt(4 * np.pi * EARTH_RADIUS_M**2 / self.n_cells) / 1000)

    def cell_polygons(self, bounds: tuple[float, float, float, float] | None = None) -> gpd.GeoDataFrame:
        """
        The area each cell represents, as a box of longitude and latitude in EPSG:4326 with longitudes from -180 to
        180: halfway in longitude to the cells beside it, and in latitude between the edges of its row. The edges split
        the sphere into bands whose areas are the gaussian quadrature weights of the rows, so each cell covers the share
        of the sphere the grid's own quadrature weights its value by.

        Args:
            bounds: optional (min x, min y, max x, max y) in degrees; only the cells overlapping it are made. Without
                it every cell of the world is made.

        Returns:
            gpd.GeoDataFrame: one row per cell, in cell order: cell_index (position in a GRIB message's values), x and
                y (the cell center longitude and latitude), and geometry
        """
        sine_latitude, quadrature_weight = np.polynomial.legendre.leggauss(2 * self.N)  # south to north
        row_latitude = np.degrees(np.arcsin(sine_latitude[::-1]))
        sine_row_edge = 1 - np.concatenate(([0], np.cumsum(quadrature_weight[::-1])))
        sine_row_edge[-1] = -1  # the running sum of the weights reaches 2, the south pole, only to within rounding
        row_edge = np.degrees(np.arcsin(sine_row_edge))
        row = np.repeat(np.arange(2 * self.N), self.cells_per_row)
        cells_before_row = np.concatenate(([0], np.cumsum(self.cells_per_row)[:-1]))
        cells_in_row = self.cells_per_row[row]
        x = (360 * (np.arange(row.shape[0]) - cells_before_row[row]) / cells_in_row + 180) % 360 - 180
        west, east = x - 180 / cells_in_row, x + 180 / cells_in_row
        south, north = row_edge[row + 1], row_edge[row]
        straddles = west < -180  # the cells centered on -180 degrees, which reach past it onto the east edge of the map
        cell_index = np.arange(row.shape[0])
        if bounds is not None:
            min_x, min_y, max_x, max_y = bounds
            overlaps_x = ((east > min_x) & (west < max_x)) | (straddles & (west + 360 < max_x))
            cell_index = np.flatnonzero(overlaps_x & (north > min_y) & (south < max_y))
        west, east = west[cell_index], east[cell_index]
        south, north = south[cell_index], north[cell_index]
        straddles = straddles[cell_index]
        geometry = shapely.box(np.maximum(west, -180), south, east, north)
        # each cell reaching past -180 degrees is the two parts of its box, one on either edge of the map
        east_edge_part = shapely.box(west[straddles] + 360, south[straddles], 180, north[straddles])
        geometry[straddles] = shapely.multipolygons(np.stack([geometry[straddles], east_edge_part], axis=1))
        columns = {'cell_index': cell_index, 'x': x[cell_index], 'y': row_latitude[row[cell_index]]}
        return gpd.GeoDataFrame(columns, geometry=geometry, crs=4326)


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
        """The weight table column cell_index, the position of each cell in the values of a GRIB message."""
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
