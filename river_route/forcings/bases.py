"""
The abstract bases of the runoff classes and the forms of runoff they hand to routing.

``RunoffForcing`` is the base of every runoff class. Each is built from a Configs with ``from_configs`` and yields its
runoff one file at a time from ``generator``, in one of the NamedTuples defined here. ``BaseGridRunoffForcing`` is the
base of the gridded runoff classes, which read runoff depths at the grid cells of a weight table and aggregate them onto
rivers. A subclass names the weight table columns that locate a cell with ``cell_dimensions`` and may read its files its
own way by overriding ``read_runoff``.

The forms of runoff routing reads:

    CatchmentRunoffVolumes  (river, time) volumes, each row the runoff of one river's catchment.
    GridCellRunoff          runoff depths at grid cells and the sparse weights of each river, read as each river is
                            routed so no series of its runoff volumes is built.

A new form of runoff needs a NamedTuple with ``check`` and ``cells_and_weights``.

Both routing kernels of river_route/kernels read every form the same way, from ``cells_and_weights``: runoff at cells
and each river's weights of them. Runoff aggregated to catchments is a cell per river weighted by one.
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

from .._metadata import __version__
from ..configs import Configs
from ..types import DatetimeArray, FloatArray, Int32Array, IntArray, PathInput, PathList, RunoffGenerator
from . import _numba_kernels as kernels

if TYPE_CHECKING:
    from ..network.Network import Network

__all__ = [
    'CATCHMENT_RUNOFF',
    'CATCHMENT_AREA',
    'VOLUME_UNITS',
    'CatchmentRunoffVolumes',
    'GridCellRunoff',
    'RunoffForcing',
    'BaseGridRunoffForcing',
]

# the fixed schema of a catchment runoff file: one depth or volume series per catchment and the area of each catchment
CATCHMENT_RUNOFF = 'catchment_runoff'
CATCHMENT_AREA = 'catchment_area'
VOLUME_UNITS = ('m3', 'm^3', 'm³')  # units attributes that mark catchment runoff as volumes, anything else is a depth

logger = logging.getLogger(__name__)


def _check_river_order(runoff_river_ids: IntArray, river_ids: IntArray) -> None:
    """Raise ValueError unless the runoff is for the rivers of the network in the same order, since routing reads the
    runoff of the river at each river_index from the row at that index."""
    if runoff_river_ids.shape != river_ids.shape:
        raise ValueError(f'the runoff is for {runoff_river_ids.shape[0]} rivers, the network has {river_ids.shape[0]}')
    if np.array_equal(runoff_river_ids, river_ids):
        return
    row = int(np.argmax(runoff_river_ids != river_ids))
    raise ValueError(
        f'Row {row} of the runoff is riverId {runoff_river_ids[row]} where the network has {river_ids[row]}. '
        f'The runoff must list the rivers of the network file in the same order.'
    )


class CatchmentRunoffVolumes(NamedTuple):
    """Each river's catchment runoff volume (m³) in every runoff step, as a C-order (river, time) array whose rows
    routing reads in place."""

    runoff: FloatArray  # (n_rivers, n_steps) float32 volumes (m³)
    river_ids: IntArray  # (n_rivers,) river id of each row

    def check(self, river_ids: IntArray, n_steps: int) -> None:
        """Raise ValueError unless it holds n_steps steps for ``river_ids`` in that order, each river's row
        contiguous."""
        if self.runoff.ndim != 2 or self.runoff.strides[1] != self.runoff.itemsize:
            raise ValueError('catchment runoff must be a (river, time) array with each river row contiguous')
        if self.runoff.shape != (river_ids.shape[0], n_steps):
            raise ValueError(
                f'catchment runoff has shape {self.runoff.shape}, expected {(river_ids.shape[0], n_steps)}'
            )
        _check_river_order(self.river_ids, river_ids)

    def cells_and_weights(self) -> tuple[FloatArray, Int32Array, Int32Array, FloatArray]:
        """The runoff as the routing kernels read it: a cell per river, its own volumes, each weighted by one."""
        n_rivers = self.runoff.shape[0]
        weight_starts = np.arange(n_rivers + 1, dtype=np.int32)
        return self.runoff, weight_starts, weight_starts[:-1], np.ones(n_rivers, dtype=np.float32)


class GridCellRunoff(NamedTuple):
    """
    One file of gridded runoff depths read at the weight table's grid cells, with the weights that turn them into each
    river's catchment runoff volume as routing routes it.
    """
    runoff_by_cell: FloatArray  # (n_cells, time) C-order runoff depths, each cell's time series contiguous
    indptr: IntArray  # (n_rivers + 1,) sparse row pointers to the weights of each river_index
    cell: IntArray  # (n_weights,) row of runoff_by_cell each weight applies to
    # (n_weights,) float32 volume (m³) a unit of the cell's runoff depth gives the river: the weight's proportion of the
    # river's catchment area, times that area, times the factor converting the depth unit to meters
    weight: FloatArray
    river_ids: IntArray  # (n_rivers,) river id of each river the weights describe, in the order of indptr

    def check(self, river_ids: IntArray, n_steps: int) -> None:
        """Raise ValueError unless the weight table describes ``river_ids`` in that order and the runoff has n_steps
        steps."""
        n_rivers = river_ids.shape[0]
        if self.runoff_by_cell.ndim != 2 or self.runoff_by_cell.shape[1] < n_steps:
            raise ValueError(f'runoff_by_cell has shape {self.runoff_by_cell.shape}, expected (n_cells, >= {n_steps})')
        if not self.runoff_by_cell.flags.c_contiguous:
            raise ValueError('runoff_by_cell must be C-contiguous, each cell series contiguous')
        if self.indptr.shape[0] != n_rivers + 1 or self.weight.shape != self.cell.shape:
            raise ValueError(f'the weight table does not describe the {n_rivers} rivers being routed')
        if self.cell.size and self.cell.max() >= self.runoff_by_cell.shape[0]:
            raise ValueError(f'the weight table references cells beyond the {self.runoff_by_cell.shape[0]} read')
        _check_river_order(self.river_ids, river_ids)

    def cells_and_weights(self) -> tuple[FloatArray, Int32Array, Int32Array, FloatArray]:
        """The runoff as the routing kernels read it: the depths of the cells, and each river's weights of them."""
        return (
            self.runoff_by_cell,
            self.indptr.astype(np.int32),
            self.cell.astype(np.int32),
            self.weight.astype(np.float32, copy=False),
        )


class RunoffForcing(ABC):
    """
    Base class for the sources of catchment runoff, the runoff volume each catchment adds to its river in each step.
    Each subclass is built from a Configs with ``from_configs`` and provides a
    ``generator`` that yields its runoff in the form routing reads, and every subclass writes catchment runoff to disk
    with ``to_netcdf`` in the one format ``CatchmentRunoffForcing`` reads.
    """

    @classmethod
    @abstractmethod
    def from_configs(cls, configs: Configs) -> Self:
        """Build this RunoffForcing from the options on a Configs, as a Router does."""

    @abstractmethod
    def generator(self, forcing_files: PathList) -> RunoffGenerator:
        """
        Yield one (dates, runoff, source_file) tuple per input, where runoff is what routing reads: its
        CatchmentRunoffVolumes, or a GridCellRunoff whose grid cells routing reads as it routes. Each checks its arrays
        with ``check``.
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

    @staticmethod
    def to_netcdf(
        path: PathInput,
        dates: DatetimeArray,
        catchment_runoff: FloatArray,
        river_ids: IntArray,
        catchment_area: FloatArray,
        *,
        as_volumes: bool,
    ) -> None:
        """
        Write a catchment runoff array to a netCDF file that ``CatchmentRunoffForcing`` reads, routed as
        ``forcing_files`` with ``forcing`` catchment.

        Args:
            path: netCDF file to write
            dates: (time,) datetime64 values of the steps
            catchment_runoff: (n_rivers, time) depths (m) or volumes (m³), rows ordered like ``river_ids``
            river_ids: (n_rivers,) river id of each row
            catchment_area: (n_rivers,) area of each catchment in m², the factor between depths and volumes
            as_volumes: whether catchment_runoff holds volumes (m³) rather than depths (m), which sets the units
                written. It is given with the array rather than read from a RunoffForcing's options, since
                ``catchment_reader`` yields volumes whatever ``as_volumes`` option a grid runoff was built with.
        """
        RunoffForcing._catchment_runoff_dataset(
            dates, catchment_runoff, river_ids, catchment_area, as_volumes
        ).to_netcdf(path)
        return

    @staticmethod
    def _catchment_runoff_dataset(
        dates: DatetimeArray,
        catchment_runoff: FloatArray,
        river_ids: IntArray,
        catchment_area: FloatArray,
        as_volumes: bool,
    ) -> xr.Dataset:
        """Build the catchment runoff dataset with dimensions (``riverId``, ``time``) that ``to_netcdf`` writes."""
        units = 'm3' if as_volumes else 'm'
        kind = 'volumes' if as_volumes else 'depths'
        start_date = pd.Timestamp(dates[0]).strftime('%Y%m%d%H')
        end_date = pd.Timestamp(dates[-1]).strftime('%Y%m%d%H')
        timestep = int((dates[1] - dates[0]) / np.timedelta64(1, 's')) if len(dates) > 1 else 0
        return xr.Dataset(
            {
                CATCHMENT_RUNOFF: xr.DataArray(
                    catchment_runoff,
                    dims=('riverId', 'time'),
                    attrs={
                        'long_name': f'Incremental catchment runoff {kind}',
                        'units': units,
                        'cell_measures': f'area: {CATCHMENT_AREA}',
                    },
                ),
                CATCHMENT_AREA: xr.DataArray(
                    np.asarray(catchment_area).astype(np.float32, copy=False),
                    dims=('riverId',),
                    attrs={'long_name': 'catchment area of each river', 'units': 'm2'},
                ),
            },
            coords={
                'riverId': xr.DataArray(
                    np.asarray(river_ids).astype(np.int64, copy=False),
                    dims=('riverId',),
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
    def _get_conversion_factor(unit: str | None) -> float:
        if unit is None:
            logger.warning('No units attribute found. Assuming meters')
            return 1
        if unit in ('m', 'meters'):
            return 1
        if unit in ('mm', 'millimeters', 'kg m-2'):  # a kilogram of water on a square meter is a millimeter deep
            return 0.001
        raise ValueError(f'Unknown units: {unit}')


@dataclass(eq=False, repr=False)
class BaseGridRunoffForcing(RunoffForcing):
    """
    Prepares catchment runoff for routing from gridded runoff depths and a grid weight table. Each subclass names the
    weight table columns that locate a grid cell and the runoff dimension each one indexes, with ``cell_dimensions``.

    Only ``grid_weights_file`` is required. Every option is named the same here as it is on ``Configs``, so
    ``from_configs`` reads each one off a ``Configs`` by its own name, which is how a ``Router`` builds one.

    Reading and indexing the weight table is the same work for every forcing file, so it is done once when the
    runoff is created and reused for all of them. Each forcing file is read at the grid cells the table touches.
    Routing reads each river's cells directly as it routes the river, and ``aggregate`` sums them onto rivers as area
    weighted depths or volumes one river at a time, each river's series written straight into its row of a C-order
    (river, time) array.
    """

    grid_weights_file: PathInput  # weight table netCDF produced by the functions in ``river_route.forcings.weights``
    _: KW_ONLY
    var_runoff: str = 'ro'  # name of the runoff variable in the gridded forcing files
    var_t: str = 'time'  # name of the time coordinate variable
    accumulation_type: Literal['incremental', 'cumulative'] = 'incremental'
    force_positive_runoff: bool = False  # clip negative runoff depths to zero
    runoff_depth_unit: str | None = None  # unit of the runoff depths; None reads the file attributes, else meters
    as_volumes: bool = False  # prepare volumes (m3) instead of depths (m); routing always uses volumes

    # Weight table read once by __post_init__, in row order of the table which follows the network file order
    river_ids: IntArray = field(init=False)  # (n_rivers,) river id of each catchment runoff column
    # each cell_dimensions column's (n_cells,) index of every unique cell any catchment touches
    cell_indexes: dict[str, IntArray] = field(init=False)
    # (n_rivers + 1,) sparse row pointers to the weights of each river_index
    indptr: IntArray = field(init=False)
    cell: IntArray = field(init=False)  # (n_weights,) position in cell_indexes of the cell of each weight
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
        return self.accumulation_type == 'cumulative'

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

    def _read_weights(self) -> None:
        """Read the weight table netCDF into the sparse arrays the aggregation kernel consumes."""
        cell_columns = list(self.cell_dimensions)
        with xr.open_dataset(self.grid_weights_file) as ds:
            weight_df = ds[['riverId', *cell_columns, 'proportion', 'area_sqm']].to_dataframe()
        unique_indexes = weight_df[cell_columns].drop_duplicates().reset_index(drop=True).reset_index().astype(int)
        # index already topo sorted
        river_ids = weight_df[['riverId']].drop_duplicates().sort_index()['riverId'].to_numpy(dtype=np.int64)

        cells = weight_df[cell_columns].merge(unique_indexes, on=cell_columns, how='left')
        point_idx = cells['index'].to_numpy()
        river_id_to_row = pd.Series(np.arange(len(river_ids)), index=river_ids)
        river_idx = river_id_to_row.loc[weight_df['riverId'].to_numpy()].to_numpy()
        matrix = scipy.sparse.csr_matrix(
            (weight_df['proportion'].to_numpy(), (river_idx, point_idx)), shape=(len(river_ids), len(unique_indexes))
        )
        self.river_ids = river_ids
        self.cell_indexes = {column: unique_indexes[column].to_numpy() for column in cell_columns}
        self.indptr = matrix.indptr
        self.cell = matrix.indices
        self.proportion = matrix.data
        self.catchment_area = weight_df.groupby('riverId')['area_sqm'].sum().reindex(river_ids).to_numpy()
        return

    def distribute(self, network: Network) -> None:
        """
        Match the weight table to a stabilized network: each river's runoff volume is split between every row with it
        as parentRiverId, which is the river and the synthetic sub-reaches injected upstream of it, in proportion to
        their travel times, the share of the river each one is. Does nothing when the table already has a row per
        network row.
        """
        n = network.river_ids.shape[0]
        if self.river_ids.shape[0] == n:
            return
        parent = pd.Index(self.river_ids).get_indexer(network.parent_river_ids)  # weight table row of each parent
        if np.any(parent < 0):
            raise ValueError(f'{self.grid_weights_file} is missing parentRiverId values of the network')
        counts = np.diff(self.indptr)[parent]
        indptr = np.concatenate(([0], np.cumsum(counts))).astype(self.indptr.dtype)
        take = np.repeat(self.indptr[parent] - indptr[:-1], counts) + np.arange(indptr[-1])
        self.indptr, self.cell, self.proportion = indptr, self.cell[take], self.proportion[take]
        k = network.k.astype(np.float64)
        share = k / np.bincount(parent, weights=k, minlength=len(self.river_ids))[parent]
        self.catchment_area = (self.catchment_area[parent] * share).astype(self.catchment_area.dtype)
        self.river_ids = network.river_ids
        self._buffer = np.empty(0, dtype=np.float32)
        return

    def catchment_runoff(
        self,
        runoff_data: PathInput | list[PathInput],
        thread_pool: Executor | None = None,
        threads: int = 1,
        *,
        as_volumes: bool | None = None,
    ) -> tuple[FloatArray, DatetimeArray]:
        """
        Read and aggregate runoff into one reused output buffer, so repeated calls allocate no new output. Reading the
        runoff still allocates its arrays on every call.

        Args:
            runoff_data: path(s) to gridded forcing files
            thread_pool: optional thread pool to aggregate on, used as given and never shut down here
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given
            as_volumes: volumes (m³) rather than depths (m); None follows the ``as_volumes`` option

        Returns:
            tuple: (catchment runoff as a C-order (n_rivers, time) array, its time values). The array is a view of
                the reused buffer and is overwritten by the next call.
        """
        runoff, time_index, conversion_factor = self.read_runoff(runoff_data)
        if self._buffer.shape[0] < runoff.shape[0] * self.river_ids.shape[0]:
            self._buffer = np.empty(runoff.shape[0] * self.river_ids.shape[0], dtype=np.float32)
        catchment_runoff = self.aggregate(
            runoff, conversion_factor, out=self._buffer, thread_pool=thread_pool, threads=threads, as_volumes=as_volumes
        )
        return catchment_runoff, time_index

    def catchment_reader(
        self, forcing_files: PathList, thread_pool: Executor | None = None, threads: int = 1
    ) -> Iterator[tuple[DatetimeArray, FloatArray, PathInput]]:
        """
        Aggregate gridded forcing files into catchment runoff volumes, one file at a time, with the weight table this
        instance already holds. Routing does not use this: ``generator`` hands routing the grid cells to aggregate as
        it routes.

        Every file is aggregated into one reused C-order buffer: no per-file allocation of the output or copy of it.
        The yielded array is overwritten by the next file, so a consumer that needs to keep it must copy it.

        Args:
            forcing_files: gridded forcing files to aggregate
            thread_pool: optional thread pool to aggregate on, used as given and never shut down here
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given

        Yields:
            tuple: (dates, catchment_runoff, source_file) per input file
        """
        for forcing_file in forcing_files:
            catchment_runoff, dates = self.catchment_runoff(
                forcing_file, thread_pool=thread_pool, threads=threads, as_volumes=True
            )
            yield dates.astype('datetime64[s]'), catchment_runoff.astype(np.float32, copy=False), forcing_file

    def generator(self, forcing_files: PathList) -> RunoffGenerator:
        """
        Read gridded forcing files for routing, which reads each river's grid cells as it routes the river, so no
        catchment runoff array is built. A file whose catchment runoff must be de-accumulated or clipped at zero needs
        each river's whole catchment runoff series first, so it is aggregated here and yielded as catchment runoff
        instead.

        Args:
            forcing_files: gridded forcing files to read

        Yields:
            tuple: (dates, runoff, source_file) per input file, where runoff is a GridCellRunoff, or
                CatchmentRunoffVolumes when the file was aggregated here
        """
        for forcing_file in forcing_files:
            runoff, time_index, conversion_factor = self.read_runoff(forcing_file)
            if self.cumulative or self.force_positive_runoff:
                catchment_runoff = self.aggregate(runoff, conversion_factor, as_volumes=True)
                aggregated = CatchmentRunoffVolumes(catchment_runoff.astype(np.float32, copy=False), self.river_ids)
                yield time_index.astype('datetime64[s]'), aggregated, forcing_file
                continue
            volume_of_depth = self._weight(conversion_factor) * np.repeat(self.catchment_area, np.diff(self.indptr))
            forcing = GridCellRunoff(
                runoff_by_cell=self._by_cell(runoff),
                indptr=self.indptr,
                cell=self.cell,
                weight=volume_of_depth.astype(np.float32),
                river_ids=self.river_ids,
            )
            del runoff
            yield time_index.astype('datetime64[s]'), forcing, forcing_file

    def to_dataset(
        self, runoff_data: PathInput | list[PathInput], thread_pool: Executor | None = None, threads: int = 1
    ) -> xr.Dataset:
        """
        Read and aggregate runoff into a catchment runoff dataset that can be saved and routed as ``forcing_files``
        with ``forcing`` catchment.

        Args:
            runoff_data: path(s) to gridded forcing files
            thread_pool: optional thread pool to aggregate on, used as given and never shut down here
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given

        Returns:
            xr.Dataset: catchment runoff with dimensions (``riverId``, ``time``). Contains ``catchment_runoff`` in
                meters or m³ and ``catchment_area`` in m².
        """
        runoff, time_index, conversion_factor = self.read_runoff(runoff_data)
        catchment_runoff = self.aggregate(runoff, conversion_factor, thread_pool=thread_pool, threads=threads)
        del runoff
        return self._catchment_runoff_dataset(
            time_index, catchment_runoff, self.river_ids, self.catchment_area, self.as_volumes
        )

    def aggregate_to_file(
        self,
        runoff_data: PathInput | list[PathInput],
        path: PathInput,
        thread_pool: Executor | None = None,
        threads: int = 1,
    ) -> None:
        """
        Precompute the catchment runoff of gridded runoff and write it to a netCDF file that is routed as
        ``forcing_files`` with ``forcing`` catchment. Routing the grids directly aggregates inside the routing
        kernel instead, which is faster than routing a precomputed file.

        Args:
            runoff_data: path(s) to gridded forcing files
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
            runoff_data: path(s) to gridded forcing files

        Returns:
            tuple: (runoff as a C-order (time, n_cells) array, time values, factor converting the depth unit to meters)
        """
        with xr.open_mfdataset(runoff_data) as ds:
            runoff_depth_unit = self.runoff_depth_unit or ds[self.var_runoff].attrs.get('units', 'm')
            conversion_factor = self._get_conversion_factor(runoff_depth_unit)
            runoff = (
                ds[self.var_runoff]
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
        conversion_factor: float = 1,
        out: FloatArray | None = None,
        thread_pool: Executor | None = None,
        threads: int = 1,
        *,
        as_volumes: bool | None = None,
    ) -> FloatArray:
        """
        Aggregate gridded runoff depths onto rivers as area weighted depths or volumes, one river at a time.

        Args:
            runoff: (time, n_cells) runoff depths at the weight table's cells, as returned by ``read_runoff``
            conversion_factor: multiplier converting the runoff depth unit to meters
            out: optional flat buffer of at least n_rivers * time values, used as the output
            thread_pool: optional thread pool, used as given and never shut down here. The rivers are split into
                ``threads`` ranges of similar work that the same kernel aggregates concurrently. Without a pool one
                range covers every river.
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given
            as_volumes: volumes (m³) rather than depths (m); None follows the ``as_volumes`` option

        Returns:
            FloatArray: catchment runoff as a C-order (n_rivers, time) array, a view of the start of ``out`` when it is
                given
        """
        volumes = self.as_volumes if as_volumes is None else as_volumes
        n_steps, n_rivers = runoff.shape[0], self.river_ids.shape[0]
        runoff_by_cell = self._by_cell(runoff)
        weight = self._weight(conversion_factor)
        dtype = np.result_type(runoff.dtype, weight.dtype)
        if out is None:
            catchment_runoff = np.empty((n_rivers, n_steps), dtype=dtype)
        elif out.ndim != 1 or out.shape[0] < n_rivers * n_steps:
            raise ValueError(f'out must be a flat buffer of at least {n_rivers * n_steps} values')
        else:
            catchment_runoff = out[: n_rivers * n_steps].reshape(n_rivers, n_steps)
        scale = self.catchment_area if volumes else self.catchment_area[:0]  # empty to skip scaling
        self._aggregate_over_ranges(runoff_by_cell, weight, scale, catchment_runoff, thread_pool, threads)
        return catchment_runoff

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
            kernels.aggregate_grid_cells_to_rivers(
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

        A river costs one loop over its time series per weight plus one for its conversions and output, so ranges are
        balanced on weights plus rivers rather than on river count alone.

        Returns:
            np.ndarray: range boundaries ``[0, ..., n_rivers]``; range i is rivers ``bounds[i]`` to ``bounds[i + 1]``
        """
        n_rivers = indptr.shape[0] - 1
        work = indptr + np.arange(n_rivers + 1)
        targets = np.linspace(0, work[-1], max(1, n_ranges) + 1)
        return np.unique(np.concatenate(([0], np.searchsorted(work, targets), [n_rivers])))
