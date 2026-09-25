import logging
from collections.abc import Iterator
from concurrent.futures import Executor
from dataclasses import KW_ONLY, dataclass, field, fields
from typing import TYPE_CHECKING, Literal, NamedTuple, Self

import numpy as np
import pandas as pd
import scipy.sparse
import xarray as xr

from ..configs import Configs
from ..types import DatetimeArray, FloatArray, IntArray, PathInput, PathList, RunoffGenerator
from . import _numba_kernels as kernels
from .Runoff import Runoff

if TYPE_CHECKING:
    from ..network import Network

__all__ = ['CellRunoff', 'GaussianGridRunoff']

logger = logging.getLogger(__name__)


class CellRunoff(NamedTuple):
    """
    One runoff file read at the weight table's cells, with everything the gaussian grid kernel needs to aggregate it
    onto the rivers as it routes them. The fields are in the order ``aggregate_river`` takes them.
    """

    runoff_by_cell: FloatArray  # (n_cells, time) C-order runoff depths, each cell's time series contiguous
    indptr: IntArray  # (n_rivers + 1,) sparse row pointers, the weights of river r are indptr[r]:indptr[r + 1]
    cell: IntArray  # (n_weights,) row of runoff_by_cell each weight applies to
    weight: FloatArray  # (n_weights,) proportions already multiplied by the depth unit conversion factor
    scale: FloatArray  # (n_rivers,) catchment areas that turn depths into volumes; EMPTY to keep depths
    cumulative: bool  # de-accumulate each river's series from cumulative to incremental
    force_positive: bool  # clip negative runoff to zero


@dataclass(eq=False, repr=False)
class GaussianGridRunoff(Runoff):
    """
    Prepares catchment runoff for routing from gridded runoff depths and a grid weight table.

    Only ``grid_weights_file`` is required. Every option is named the same here as it is on ``Configs``, so
    ``GaussianGridRunoff.from_configs`` reads each one off a ``Configs`` by its own name, which is how a ``Router``
    builds one.

    Reading and indexing the weight table is the same work for every runoff file, so it is done once when the
    GaussianGridRunoff is created and reused for all of them. Each runoff file is read at the grid cells the table
    touches and aggregated onto rivers as area weighted depths or volumes in a single pass, each river's series
    written straight into its row of a C-order (river, time) array that the routing kernels read in place.
    """

    grid_weights_file: PathInput  # weight table netCDF produced by ``river_route.runoff.grid_weights``
    _: KW_ONLY
    var_river_id: str = 'river_id'  # name of the river id variable in the weight table
    var_grid_runoff: str = 'ro'  # name of the runoff variable in the gridded runoff files
    var_x: str = 'x'  # name of the grid x coordinate variable
    var_y: str = 'y'  # name of the grid y coordinate variable
    var_t: str = 'time'  # name of the time coordinate variable
    # ``incremental`` runoff per step, or ``cumulative`` totals to difference
    grid_accumulation_type: Literal['incremental', 'cumulative'] = 'incremental'
    runoff_depth_unit: str | None = None  # unit of the runoff depths; None reads the file attributes, else meters
    force_positive_runoff: bool = False  # clip negative runoff depths to zero
    force_uniform_timesteps: bool = True  # resample runoff with irregular timesteps to the first timestep
    as_volumes: bool = False  # prepare volumes (m3) instead of depths (m); routing always uses volumes

    # Weight table read once by __post_init__, in row order of the table which follows the routing params order
    river_ids: IntArray = field(init=False)  # (n_rivers,) river id of each catchment runoff column
    x_index: IntArray = field(init=False)  # (n_cells,) grid x index of each unique cell any catchment touches
    y_index: IntArray = field(init=False)  # (n_cells,) grid y index of each unique cell
    # (n_rivers + 1,) sparse row pointers, the weights of river r are indptr[r]:indptr[r + 1]
    indptr: IntArray = field(init=False)
    cell: IntArray = field(init=False)  # (n_weights,) position in x_index and y_index of the cell of each weight
    proportion: FloatArray = field(init=False)  # (n_weights,) share of the river's catchment area inside the cell
    catchment_area: FloatArray = field(init=False)  # (n_rivers,) total catchment area of each river in m²
    _buffer: FloatArray = field(init=False)  # flat reused catchment_runoff() output, sized to the longest file

    def __post_init__(self) -> None:
        if self.grid_weights_file is None:
            raise ValueError('grid_weights_file is required to build a GaussianGridRunoff')
        self._read_weights()
        self._buffer = np.empty(0, dtype=np.float32)
        return

    @property
    def var_runoff(self) -> str:
        """The ``Runoff`` interface name for ``var_grid_runoff``, the runoff variable of the gridded files."""
        return self.var_grid_runoff

    @property
    def cumulative(self) -> bool:
        """Whether the runoff values are cumulative totals to difference rather than incremental per step."""
        return self.grid_accumulation_type == 'cumulative'

    @classmethod
    def from_configs(cls, configs: Configs) -> Self:
        """Build a GaussianGridRunoff from the runoff options on a ``Configs``. The configs are validated for runoff
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
        return f'{type(self).__name__}(n_rivers={self.river_ids.shape[0]}, n_cells={self.x_index.shape[0]})'

    ################################################
    # Weight table
    ################################################

    def _read_weights(self) -> None:
        """Read the weight table netCDF into the sparse arrays the aggregation kernel consumes."""
        var_river_id = self.var_river_id
        with xr.open_dataset(self.grid_weights_file) as ds:
            weight_df = ds[[var_river_id, 'x_index', 'y_index', 'proportion', 'area_sqm']].to_dataframe()
        unique_indexes = (
            weight_df[['x_index', 'y_index']].drop_duplicates().reset_index(drop=True).reset_index().astype(int)
        )
        # index already topo sorted
        river_ids = weight_df[[var_river_id]].drop_duplicates().sort_index()[var_river_id].to_numpy(dtype=np.int32)

        cells = weight_df[['x_index', 'y_index']].merge(unique_indexes, on=['x_index', 'y_index'], how='left')
        point_idx = cells['index'].values
        river_id_to_row = pd.Series(np.arange(len(river_ids)), index=river_ids)
        river_idx = river_id_to_row.loc[weight_df[var_river_id].values].values
        matrix = scipy.sparse.csr_matrix(
            (weight_df['proportion'].values, (river_idx, point_idx)), shape=(len(river_ids), len(unique_indexes))
        )
        self.river_ids = river_ids
        self.x_index = unique_indexes['x_index'].to_numpy()
        self.y_index = unique_indexes['y_index'].to_numpy()
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
    ) -> RunoffGenerator:
        """
        Aggregate gridded runoff files into catchment runoff volumes, one file at a time, with the weight table this
        instance already holds. Routing does not use this: ``reader`` hands the gaussian grid kernel the grid cells to
        aggregate as it routes.

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

    def generator(self, runoff_files: PathList) -> Iterator[tuple[DatetimeArray, CellRunoff | FloatArray, PathInput]]:
        """
        Read gridded runoff files for the gaussian grid kernel, which aggregates and routes in one pass, so no
        catchment runoff array is built. A file whose timesteps must be resampled cannot be aggregated inside the
        routing sweep, so it is aggregated here and yielded as a catchment runoff array instead.

        Args:
            runoff_files: gridded runoff files to read

        Yields:
            tuple: (dates, forcing, source_file) per input file, where forcing is a CellRunoff, or a C-order
                (n_rivers, time) catchment runoff array when the file was resampled
        """
        self.as_volumes = True  # routing uses catchment runoff volumes
        for runoff_file in runoff_files:
            runoff, time_index, conversion_factor = self.read_runoff(runoff_file)
            if self._needs_resampling(time_index):
                catchment_runoff, time_index = self.aggregate(runoff, time_index, conversion_factor)
                yield time_index.astype('datetime64[s]'), catchment_runoff.astype(np.float32, copy=False), runoff_file
                continue
            forcing = CellRunoff(
                runoff_by_cell=self._by_cell(runoff),
                indptr=self.indptr,
                cell=self.cell,
                weight=self._weight(conversion_factor),
                scale=self.catchment_area if self.as_volumes else self.catchment_area[:0],
                cumulative=self.cumulative,
                force_positive=self.force_positive_runoff,
            )
            del runoff
            yield time_index.astype('datetime64[s]'), forcing, runoff_file

    def to_dataset(
        self, runoff_data: PathInput | list[PathInput], thread_pool: Executor | None = None, threads: int = 1
    ) -> xr.Dataset:
        """
        Read and aggregate runoff into a catchment runoff dataset that can be saved and routed as ``runoff_files``
        with ``runoff_type`` catchment.

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
        ``runoff_files`` with ``runoff_type`` catchment. Routing the grids directly aggregates inside the routing
        kernel instead, which is faster than routing a precomputed file.

        Args:
            runoff_data: path(s) to runoff files
            path: netCDF file to write
            thread_pool: optional thread pool to aggregate on, used as given and never shut down here
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given
        """
        self.to_dataset(runoff_data, thread_pool=thread_pool, threads=threads).to_netcdf(path)
        return

    def read_runoff(self, runoff_data: PathInput | list[PathInput]) -> tuple[FloatArray, DatetimeArray, int | float]:
        """
        Read runoff at the grid cells the weight table touches.

        Args:
            runoff_data: path(s) to runoff files

        Returns:
            tuple: (runoff as a C-order (time, n_cells) array, time values, factor converting the depth unit to meters)
        """
        with xr.open_mfdataset(runoff_data, chunks=None) as ds:
            runoff_depth_unit = self.runoff_depth_unit or ds[self.var_runoff].attrs.get('units', 'm')
            conversion_factor = self._get_conversion_factor(runoff_depth_unit)
            runoff = (
                ds[self.var_runoff]
                .isel(
                    {
                        self.var_x: xr.DataArray(self.x_index, dims='points'),
                        self.var_y: xr.DataArray(self.y_index, dims='points'),
                    }
                )
                .transpose(self.var_t, 'points')
                .values
            )
            time_index = ds[self.var_t].to_numpy()
        return np.ascontiguousarray(runoff), time_index, conversion_factor

    def aggregate(
        self,
        runoff: FloatArray,
        time_index: DatetimeArray,
        conversion_factor: int | float = 1,
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
            time_index = df.index.values
            # resampling works on (time, river) columns, so this rare path transposes back to (river, time) once
            catchment_runoff = np.ascontiguousarray(df.to_numpy(dtype=np.float32).T)
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

    def _weight(self, conversion_factor: int | float) -> FloatArray:
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
            for i in range(n_ranges):
                run(i)
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
