"""
The abstract bases of the runoff classes and the forms of runoff they hand to routing.

``Runoff`` is the base of every runoff class. Each is built from a Configs with ``from_configs`` and yields its runoff
one file at a time from ``generator``, in one of the NamedTuples defined here. ``BaseGridRunoff`` is the base of the
gridded runoff classes, which read runoff depths at the grid cells of a weight table and aggregate them onto rivers. A
subclass names the weight table columns that locate a cell with ``cell_dimensions`` and may read its files its own way
by overriding ``read_runoff``.

A routing job takes each river through the stages without a body in router/_numba_kernels.py. Each form of runoff
implements two of them with numba overloads registered here, next to its NamedTuple:

    get_river_forcing   river r's catchment runoff, from the runoff a file was read into
    add_river_forcing   adds that catchment runoff into the river's forcing

    CatchmentRunoffVolumes  (river, time) volumes whose row r is river r's catchment runoff series. A series, a 1D
                            array, is also the form of the catchment runoff a transform gives.
    GridCellRunoff          runoff depths at grid cells and the sparse weights of each river. River r's catchment
                            runoff is a RiverGridCells, its cells and weights, added into its forcing one cell at a time
                            so no catchment runoff series is built.

A new form of runoff needs a NamedTuple with ``check`` and ``first_steps``, and an overload of both stages for it.
"""

import logging
from abc import ABC, abstractmethod
from collections.abc import Iterator
from concurrent.futures import Executor
from dataclasses import KW_ONLY, dataclass, field, fields
from typing import TYPE_CHECKING, Literal, NamedTuple, Self

import numpy as np
import pandas as pd
import scipy.sparse
import xarray as xr
from numba import types
from numba.extending import overload

from .._metadata import __version__
from ..configs import Configs
from ..router._numba_kernels import add_river_forcing, get_river_forcing, is_argument_type
from ..types import DatetimeArray, FloatArray, IntArray, PathInput, PathList, RunoffGenerator
from . import _numba_kernels as kernels

if TYPE_CHECKING:
    from ..network.Network import Network

__all__ = [
    'CATCHMENT_RUNOFF',
    'CATCHMENT_AREA',
    'VOLUME_UNITS',
    'CatchmentRunoffVolumes',
    'GridCellRunoff',
    'Runoff',
    'BaseGridRunoff',
]

# the fixed schema of a catchment runoff file: one depth or volume series per catchment and the area of each catchment
CATCHMENT_RUNOFF = 'catchment_runoff'
CATCHMENT_AREA = 'catchment_area'
VOLUME_UNITS = ('m3', 'm^3', 'm³')  # units attributes that mark catchment runoff as volumes, anything else is a depth

logger = logging.getLogger(__name__)


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


@overload(get_river_forcing)
def _catchment_runoff_of_river(runoff, r):
    if is_argument_type(runoff, CatchmentRunoffVolumes):
        return lambda runoff, r: runoff.runoff[r]


def _add_catchment_runoff_series(catchment_runoff, multiplier, n_steps, n_per_step, work):
    """Add multiplier times each step of a river's catchment runoff series, read in place, into its steps of work."""
    if n_per_step == 1:
        for t in range(n_steps):
            work[t] += multiplier * np.float32(catchment_runoff[t])
    else:
        for t in range(n_steps):
            external = multiplier * np.float32(catchment_runoff[t])
            for h in range(t * n_per_step, (t + 1) * n_per_step):
                work[h] += external


# a series is the form of a CatchmentRunoffVolumes row, and of any catchment runoff a transform gives
@overload(add_river_forcing, jit_options={'nogil': True, 'fastmath': {'contract'}})
def _catchment_runoff_series_is_added(catchment_runoff, multiplier, n_steps, n_per_step, work):
    if isinstance(catchment_runoff, types.Array) and catchment_runoff.ndim == 1:
        return _add_catchment_runoff_series


class GridCellRunoff(NamedTuple):
    """
    One file of gridded runoff depths read at the weight table's grid cells, with the weights that turn them into each
    river's catchment runoff volume as routing routes it.
    """

    runoff_by_cell: FloatArray  # (n_cells, time) C-order runoff depths, each cell's time series contiguous
    indptr: IntArray  # (n_rivers + 1,) sparse row pointers, the weights of river r are indptr[r]:indptr[r + 1]
    cell: IntArray  # (n_weights,) row of runoff_by_cell each weight applies to
    # (n_weights,) float32 volume (m³) a unit of the cell's runoff depth gives the river: the weight's proportion of the
    # river's catchment area, times that area, times the factor converting the depth unit to meters
    weight: FloatArray

    def check(self, n_rivers: int, n_steps: int) -> None:
        """Raise ValueError unless the weight table describes n_rivers rivers and the runoff has n_steps steps."""
        if self.runoff_by_cell.ndim != 2 or self.runoff_by_cell.shape[1] < n_steps:
            raise ValueError(f'runoff_by_cell has shape {self.runoff_by_cell.shape}, expected (n_cells, >= {n_steps})')
        if not self.runoff_by_cell.flags.c_contiguous:
            raise ValueError('runoff_by_cell must be C-contiguous, each cell series contiguous')
        if self.indptr.shape[0] != n_rivers + 1 or self.weight.shape != self.cell.shape:
            raise ValueError(f'the weight table does not describe the {n_rivers} rivers being routed')
        if self.cell.size and self.cell.max() >= self.runoff_by_cell.shape[0]:
            raise ValueError(f'the weight table references cells beyond the {self.runoff_by_cell.shape[0]} read')

    def first_steps(self, n_steps: int) -> GridCellRunoff:
        """Itself, since routing reads only the steps it routes."""
        return self


class RiverGridCells(NamedTuple):
    """The grid cells of one river's catchment: the runoff of every cell, and the cell and weight of each of its own."""

    runoff_by_cell: FloatArray  # (n_cells, time) C-order runoff depths
    cell: IntArray  # (n_river_weights,) row of runoff_by_cell each of the river's weights applies to
    weight: FloatArray  # (n_river_weights,) volume a unit of the cell's runoff depth gives the river


@overload(get_river_forcing)
def _grid_cells_of_river(runoff, r):
    if is_argument_type(runoff, GridCellRunoff):
        return lambda runoff, r: RiverGridCells(
            runoff.runoff_by_cell,
            runoff.cell[runoff.indptr[r] : runoff.indptr[r + 1]],
            runoff.weight[runoff.indptr[r] : runoff.indptr[r + 1]],
        )


def _add_runoff_of_river_grid_cells(catchment_runoff, multiplier, n_steps, n_per_step, work):
    """
    Add multiplier times the river's catchment runoff volume into its steps of work a cell at a time, each cell's
    runoff depths read straight from its row, so no catchment runoff series is built first. The C kernel of jsrr, the
    browser port of river-route, reads the cells this way, and doing the same routed the Columbia about 1.4 times faster
    than aggregating the runoff of 64 rivers at a time into scratch rows and reading those as each river was routed.
    """
    for k in range(catchment_runoff.cell.shape[0]):
        series = catchment_runoff.runoff_by_cell[catchment_runoff.cell[k]]
        forcing = multiplier * catchment_runoff.weight[k]
        if n_per_step == 1:
            for t in range(n_steps):
                work[t] += forcing * series[t]
        else:
            for t in range(n_steps):
                external = forcing * series[t]
                for h in range(t * n_per_step, (t + 1) * n_per_step):
                    work[h] += external


@overload(add_river_forcing, jit_options={'nogil': True, 'fastmath': {'contract'}})
def _river_grid_cells_are_added(catchment_runoff, multiplier, n_steps, n_per_step, work):
    if is_argument_type(catchment_runoff, RiverGridCells):
        return _add_runoff_of_river_grid_cells


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
        CatchmentRunoffVolumes, or a GridCellRunoff whose grid cells routing reads as it routes. Each checks its arrays
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
        with ``forcing`` catchment.

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


@dataclass(eq=False, repr=False)
class BaseGridRunoff(Runoff):
    """
    Prepares catchment runoff for routing from gridded runoff depths and a grid weight table. Each subclass names the
    weight table columns that locate a grid cell and the runoff dimension each one indexes, with ``cell_dimensions``.

    Only ``grid_weights_file`` is required. Every option is named the same here as it is on ``Configs``, so
    ``from_configs`` reads each one off a ``Configs`` by its own name, which is how a ``Router`` builds one.

    Reading and indexing the weight table is the same work for every runoff file, so it is done once when the
    runoff is created and reused for all of them. Each runoff file is read at the grid cells the table touches.
    Routing reads each river's cells directly as it routes the river, and ``aggregate`` sums them onto rivers as area
    weighted depths or volumes in a single pass, each river's series written straight into its row of a C-order
    (river, time) array.
    """

    grid_weights_file: PathInput  # weight table netCDF produced by the functions in ``river_route.runoff.weights``
    _: KW_ONLY
    var_river_id: str = 'river_id'  # name of the river id variable in the weight table
    var_grid_runoff: str = 'ro'  # name of the runoff variable in the gridded runoff files
    var_t: str = 'time'  # name of the time coordinate variable
    # ``incremental`` runoff per step, or ``cumulative`` totals to difference
    grid_accumulation_type: Literal['incremental', 'cumulative'] = 'incremental'
    runoff_depth_unit: str | None = None  # unit of the runoff depths; None reads the file attributes, else meters
    force_positive_runoff: bool = False  # clip negative runoff depths to zero
    force_uniform_timesteps: bool = True  # resample runoff with irregular timesteps to the first timestep
    as_volumes: bool = False  # prepare volumes (m3) instead of depths (m); routing always uses volumes

    # Weight table read once by __post_init__, in row order of the table which follows the routing params order
    river_ids: IntArray = field(init=False)  # (n_rivers,) river id of each catchment runoff column
    # each cell_dimensions column's (n_cells,) index of every unique cell any catchment touches
    cell_indexes: dict[str, IntArray] = field(init=False)
    # (n_rivers + 1,) sparse row pointers, the weights of river r are indptr[r]:indptr[r + 1]
    indptr: IntArray = field(init=False)
    cell: IntArray = field(init=False)  # (n_weights,) position in x_index and y_index of the cell of each weight
    proportion: FloatArray = field(init=False)  # (n_weights,) share of the river's catchment area inside the cell
    catchment_area: FloatArray = field(init=False)  # (n_rivers,) total catchment area of each river in m²
    _buffer: FloatArray = field(init=False)  # flat reused catchment_runoff() output, sized to the longest file

    def __post_init__(self) -> None:
        if self.grid_weights_file is None:
            raise ValueError(f'grid_weights_file is required to build a {type(self).__name__}')
        self._read_weights()
        self._buffer = np.empty(0, dtype=np.float32)
        return

    @property
    @abstractmethod
    def cell_dimensions(self) -> dict[str, str]:
        """Each weight table column that locates a grid cell, and the runoff dimension it indexes."""

    @property
    def cumulative(self) -> bool:
        """Whether the runoff values are cumulative totals to difference rather than incremental per step."""
        return self.grid_accumulation_type == 'cumulative'

    @classmethod
    def from_configs(cls, configs: Configs) -> Self:
        """Build this grid runoff from the runoff options on a ``Configs``. The configs are validated for runoff
        before the weight table is read. Each option is read off the Configs by the field's own name, so a field with
        no matching option raises AttributeError rather than silently keeping its default."""
        if not isinstance(configs, Configs):
            raise TypeError(
                f'from_configs takes a Configs, got {type(configs).__name__}. '
                f'Use Configs(...) or Configs.from_json(path).'
            )
        configs.validate_runoff()
        return cls(**{f.name: getattr(configs, f.name) for f in fields(cls) if f.init})

    def __repr__(self) -> str:
        n_cells = next(iter(self.cell_indexes.values())).shape[0]
        return f'{type(self).__name__}(n_rivers={self.river_ids.shape[0]}, n_cells={n_cells})'

    ################################################
    # Weight table
    ################################################

    def _read_weights(self) -> None:
        """Read the weight table netCDF into the sparse arrays the aggregation kernel consumes."""
        var_river_id = self.var_river_id
        cell_columns = list(self.cell_dimensions)
        with xr.open_dataset(self.grid_weights_file) as ds:
            weight_df = ds[[var_river_id, *cell_columns, 'proportion', 'area_sqm']].to_dataframe()
        unique_indexes = weight_df[cell_columns].drop_duplicates().reset_index(drop=True).reset_index().astype(int)
        # index already topo sorted
        river_ids = weight_df[[var_river_id]].drop_duplicates().sort_index()[var_river_id].to_numpy(dtype=np.int32)

        cells = weight_df[cell_columns].merge(unique_indexes, on=cell_columns, how='left')
        point_idx = cells['index'].to_numpy()
        river_id_to_row = pd.Series(np.arange(len(river_ids)), index=river_ids)
        river_idx = river_id_to_row.loc[weight_df[var_river_id].to_numpy()].to_numpy()
        matrix = scipy.sparse.csr_matrix(
            (weight_df['proportion'].to_numpy(), (river_idx, point_idx)), shape=(len(river_ids), len(unique_indexes))
        )
        self.river_ids = river_ids
        self.cell_indexes = {column: unique_indexes[column].to_numpy() for column in cell_columns}
        self.indptr = matrix.indptr
        self.cell = matrix.indices
        self.proportion = matrix.data
        self.catchment_area = weight_df.groupby(var_river_id)['area_sqm'].sum().reindex(river_ids).to_numpy()
        return

    def distribute(self, network: Network) -> None:
        """
        Match the weight table to a stabilized network: each river's runoff volume is split evenly between every row
        with it as parent_river_id, which is the river and the synthetic sub-reaches injected upstream of it. Does
        nothing when the table already has a row per network row.
        """
        n = network.river_ids.shape[0]
        if self.river_ids.shape[0] == n:
            return
        parent = pd.Index(self.river_ids).get_indexer(network.parent_river_ids)  # weight table row of each parent
        if np.any(parent < 0):
            raise ValueError(f'{self.grid_weights_file} is missing parent_river_id values of the network')
        counts = np.diff(self.indptr)[parent]
        indptr = np.concatenate(([0], np.cumsum(counts))).astype(self.indptr.dtype)
        take = np.repeat(self.indptr[parent] - indptr[:-1], counts) + np.arange(indptr[-1])
        self.indptr, self.cell, self.proportion = indptr, self.cell[take], self.proportion[take]
        self.catchment_area = self.catchment_area[parent] / np.bincount(parent, minlength=len(self.river_ids))[parent]
        self.river_ids = network.river_ids
        self._buffer = np.empty(0, dtype=np.float32)
        return

    ################################################
    # Prepare catchment runoff from gridded runoff
    ################################################

    def catchment_runoff(
        self, runoff_data: PathInput | list[PathInput], thread_pool: Executor | None = None, threads: int = 1
    ) -> tuple[FloatArray, DatetimeArray]:
        """
        Read and aggregate runoff into one reused buffer so that repeated calls allocate nothing.

        Args:
            runoff_data: path(s) to runoff files
            thread_pool: optional thread pool to aggregate on, used as given and never shut down here
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given

        Returns:
            tuple: (catchment runoff as a C-order (n_rivers, time) array, its time values). The array is a view of
                the reused buffer and is overwritten by the next call.
        """
        runoff, time_index, conversion_factor = self.read_runoff(runoff_data)
        if self._buffer.shape[0] < runoff.shape[0] * self.river_ids.shape[0]:
            self._buffer = np.empty(runoff.shape[0] * self.river_ids.shape[0], dtype=np.float32)
        return self.aggregate(
            runoff, time_index, conversion_factor, out=self._buffer, thread_pool=thread_pool, threads=threads
        )

    def catchment_reader(
        self, runoff_files: PathList, thread_pool: Executor | None = None, threads: int = 1
    ) -> Iterator[tuple[DatetimeArray, FloatArray, PathInput]]:
        """
        Aggregate gridded runoff files into catchment runoff volumes, one file at a time, with the weight table this
        instance already holds. Routing does not use this: ``generator`` hands routing the grid cells to aggregate as
        it routes.

        Every file is aggregated into one reused C-order buffer that the kernels read directly: no per-file allocation
        or copy. The yielded array is overwritten by the next file, so a consumer that needs to keep it must copy it.

        Args:
            runoff_files: gridded runoff files to aggregate
            thread_pool: optional thread pool to aggregate on, used as given and never shut down here
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given

        Yields:
            tuple: (dates, catchment_runoff, source_file) per input file
        """
        self.as_volumes = True  # routing uses catchment runoff volumes
        for runoff_file in runoff_files:
            catchment_runoff, dates = self.catchment_runoff(runoff_file, thread_pool=thread_pool, threads=threads)
            yield dates.astype('datetime64[s]'), catchment_runoff.astype(np.float32, copy=False), runoff_file

    def generator(self, runoff_files: PathList) -> RunoffGenerator:
        """
        Read gridded runoff files for routing, which reads each river's grid cells as it routes the river, so no
        catchment runoff array is built. A file whose timesteps must be resampled, or whose catchment runoff must be
        de-accumulated or clipped at zero, needs each river's whole catchment runoff series first, so it is aggregated
        here and yielded as catchment runoff instead.

        Args:
            runoff_files: gridded runoff files to read

        Yields:
            tuple: (dates, runoff, source_file) per input file, where runoff is a GridCellRunoff, or
                CatchmentRunoffVolumes when the file was aggregated here
        """
        self.as_volumes = True  # routing uses catchment runoff volumes
        for runoff_file in runoff_files:
            runoff, time_index, conversion_factor = self.read_runoff(runoff_file)
            if self._needs_resampling(time_index) or self.cumulative or self.force_positive_runoff:
                catchment_runoff, time_index = self.aggregate(runoff, time_index, conversion_factor)
                aggregated = CatchmentRunoffVolumes(catchment_runoff.astype(np.float32, copy=False))
                yield time_index.astype('datetime64[s]'), aggregated, runoff_file
                continue
            volume_of_depth = self._weight(conversion_factor) * np.repeat(self.catchment_area, np.diff(self.indptr))
            forcing = GridCellRunoff(
                runoff_by_cell=self._by_cell(runoff),
                indptr=self.indptr,
                cell=self.cell,
                weight=volume_of_depth.astype(np.float32),
            )
            del runoff
            yield time_index.astype('datetime64[s]'), forcing, runoff_file

    def to_dataset(
        self, runoff_data: PathInput | list[PathInput], thread_pool: Executor | None = None, threads: int = 1
    ) -> xr.Dataset:
        """
        Read and aggregate runoff into a catchment runoff dataset that can be saved and routed as ``runoff_files``
        with ``forcing`` catchment.

        Args:
            runoff_data: path(s) to runoff files
            thread_pool: optional thread pool to aggregate on, used as given and never shut down here
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given

        Returns:
            xr.Dataset: catchment runoff with dimensions ``time`` and ``river_id``. Contains ``catchment_runoff`` in
                meters or m³ and ``catchment_area`` in m².
        """
        runoff, time_index, conversion_factor = self.read_runoff(runoff_data)
        catchment_runoff, time_index = self.aggregate(
            runoff, time_index, conversion_factor, thread_pool=thread_pool, threads=threads
        )
        del runoff
        return self._catchment_runoff_dataset(time_index, catchment_runoff, self.river_ids, self.catchment_area)

    def aggregate_to_file(
        self,
        runoff_data: PathInput | list[PathInput],
        path: PathInput,
        thread_pool: Executor | None = None,
        threads: int = 1,
    ) -> None:
        """
        Precompute the catchment runoff of gridded runoff and write it to a netCDF file that is routed as
        ``runoff_files`` with ``forcing`` catchment. Routing the grids directly aggregates inside the routing
        kernel instead, which is faster than routing a precomputed file.

        Args:
            runoff_data: path(s) to runoff files
            path: netCDF file to write
            thread_pool: optional thread pool to aggregate on, used as given and never shut down here
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given
        """
        self.to_dataset(runoff_data, thread_pool=thread_pool, threads=threads).to_netcdf(path)
        return

    def read_runoff(self, runoff_data: PathInput | list[PathInput]) -> tuple[FloatArray, DatetimeArray, float]:
        """
        Read runoff at the grid cells the weight table touches.

        Args:
            runoff_data: path(s) to runoff files

        Returns:
            tuple: (runoff as a C-order (time, n_cells) array, time values, factor converting the depth unit to meters)
        """
        with xr.open_mfdataset(runoff_data, chunks=None) as ds:
            runoff_depth_unit = self.runoff_depth_unit or ds[self.var_grid_runoff].attrs.get('units', 'm')
            conversion_factor = self._get_conversion_factor(runoff_depth_unit)
            runoff = (
                ds[self.var_grid_runoff]
                .isel(
                    {
                        dimension: xr.DataArray(self.cell_indexes[column], dims='points')
                        for column, dimension in self.cell_dimensions.items()
                    }
                )
                .transpose(self.var_t, 'points')
                .to_numpy()
            )
            time_index = ds[self.var_t].to_numpy()
        return np.ascontiguousarray(runoff), time_index, conversion_factor

    def aggregate(
        self,
        runoff: FloatArray,
        time_index: DatetimeArray,
        conversion_factor: float = 1,
        out: FloatArray | None = None,
        thread_pool: Executor | None = None,
        threads: int = 1,
    ) -> tuple[FloatArray, DatetimeArray]:
        """
        Aggregate gridded runoff depths onto rivers as area weighted depths or volumes in a single pass.

        Args:
            runoff: (time, n_cells) runoff depths at the weight table's cells, as returned by ``read_runoff``
            time_index: (time,) datetime64 values of the runoff steps
            conversion_factor: multiplier converting the runoff depth unit to meters
            out: optional flat buffer of at least n_rivers * time values, used as the output. Not used when irregular
                timesteps are resampled.
            thread_pool: optional thread pool, used as given and never shut down here. The rivers are split into
                ``threads`` ranges of similar work that the same kernel aggregates concurrently. Without a pool one
                range covers every river.
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given

        Returns:
            tuple: (catchment runoff as a C-order (n_rivers, time) array, its time values). When ``out`` is used the
                array is a view of its start.
        """
        n_steps, n_rivers = runoff.shape[0], self.river_ids.shape[0]
        runoff_by_cell = self._by_cell(runoff)
        weight = self._weight(conversion_factor)
        dtype = np.result_type(runoff.dtype, weight.dtype)
        no_scale = self.catchment_area[:0]

        if self._needs_resampling(time_index):
            catchment_runoff = np.empty((n_rivers, n_steps), dtype=dtype)
            self._aggregate_over_ranges(runoff_by_cell, weight, no_scale, catchment_runoff, thread_pool, threads)
            timestep = int((time_index[1] - time_index[0]) / np.timedelta64(1, 's'))
            logger.warning(f'Time steps are not uniform, resampling to the first timestep: {timestep} seconds')
            df = pd.DataFrame(catchment_runoff.T, index=time_index, columns=self.river_ids)
            df = df.cumsum().resample(rule=f'{timestep}s').interpolate(method='linear')
            df = self._cumulative_to_incremental(df)
            time_index = df.index.to_numpy()
            # resampling works on (time, river) columns, so this rare path transposes back to (river, time) once. It is
            # copied explicitly: pandas hands out a read-only view of its storage, whose transpose is already contiguous
            catchment_runoff = np.array(df.to_numpy(dtype=np.float32).T, order='C')
            del df
            catchment_runoff[np.isnan(catchment_runoff)] = 0.0
            if self.as_volumes:
                catchment_runoff *= self.catchment_area[:, np.newaxis]
            return catchment_runoff, time_index

        if out is None:
            catchment_runoff = np.empty((n_rivers, n_steps), dtype=dtype)
        elif out.ndim != 1 or out.shape[0] < n_rivers * n_steps:
            raise ValueError(f'out must be a flat buffer of at least {n_rivers * n_steps} values')
        else:
            catchment_runoff = out[: n_rivers * n_steps].reshape(n_rivers, n_steps)
        scale = self.catchment_area if self.as_volumes else no_scale
        self._aggregate_over_ranges(runoff_by_cell, weight, scale, catchment_runoff, thread_pool, threads)
        return catchment_runoff, time_index

    @staticmethod
    def _by_cell(runoff: FloatArray) -> FloatArray:
        """
        The (time, n_cells) runoff as C-order (n_cells, time), each cell's series contiguous for the vectorized sums,
        with NaN replaced by zero. This is the one place NaN is handled, so a missing cell contributes nothing while
        the other cells of its catchments still count, and no kernel downstream ever sees a NaN.
        """
        runoff = np.ascontiguousarray(runoff)
        runoff_by_cell = np.empty(runoff.shape[::-1], dtype=runoff.dtype)
        kernels.cells_by_time(runoff, runoff_by_cell)
        return runoff_by_cell

    def _weight(self, conversion_factor: float) -> FloatArray:
        """The weight table proportions multiplied by the factor converting the runoff depth unit to meters."""
        return self.proportion * conversion_factor if conversion_factor != 1 else self.proportion

    def _needs_resampling(self, time_index: DatetimeArray) -> bool:
        """Whether irregular timesteps must be resampled onto the first timestep before routing."""
        time_diff = np.diff(time_index)
        return bool(self.force_uniform_timesteps and time_diff.size and not np.all(time_diff == time_diff[0]))

    def _aggregate_over_ranges(
        self,
        runoff_by_cell: FloatArray,
        weight: FloatArray,
        scale: FloatArray,
        out: FloatArray,
        thread_pool: Executor | None,
        threads: int,
    ) -> None:
        """
        Run the aggregation kernel over every river: as one range covering the network when single threaded, or as
        concurrent ranges of similar work when a pool is given. The kernel is the same either way, so there is no
        separate serial code path to keep in sync.
        """
        dtype = np.result_type(runoff_by_cell.dtype, weight.dtype)
        bounds = self._river_ranges(self.indptr, threads if thread_pool is not None else 1)

        def run(i: int) -> None:
            r_start, r_stop = int(bounds[i]), int(bounds[i + 1])
            kernels.aggregate_to_rivers(
                runoff_by_cell=runoff_by_cell,
                indptr=self.indptr,
                cell=self.cell,
                weight=weight,
                scale=scale,
                zero=dtype.type(0),
                cumulative=self.cumulative,
                force_positive=self.force_positive_runoff,
                r_start=r_start,
                r_stop=r_stop,
                out=out,
            )

        n_ranges = bounds.shape[0] - 1
        if thread_pool is None or n_ranges < 2:
            for river_range in range(n_ranges):
                run(river_range)
        else:
            list(thread_pool.map(run, range(n_ranges)))  # list() so a worker exception propagates
        return

    @staticmethod
    def _river_ranges(indptr: IntArray, n_ranges: int) -> IntArray:
        """
        Split the rivers into at most ``n_ranges`` contiguous ranges of similar aggregation work.

        A river costs one pass over its time series per weight plus one for its conversions and output, so ranges are
        balanced on weights plus rivers rather than on river count alone.

        Returns:
            np.ndarray: range boundaries ``[0, ..., n_rivers]``; range i is rivers ``bounds[i]`` to ``bounds[i + 1]``
        """
        n_rivers = indptr.shape[0] - 1
        work = indptr + np.arange(n_rivers + 1)
        targets = np.linspace(0, work[-1], max(1, n_ranges) + 1)
        return np.unique(np.concatenate(([0], np.searchsorted(work, targets), [n_rivers])))
