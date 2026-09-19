import logging
from collections.abc import Iterator
from concurrent.futures import Executor
from typing import Literal, NamedTuple, Self

import numpy as np
import pandas as pd
import scipy.sparse
import xarray as xr

from ..configs import Configs
from ..types import DatetimeArray, FloatArray, IntArray, PathInput, PathList, VlateralGenerator
from . import _numba_kernels as kernels
from .Runoff import Runoff

__all__ = ['CellRunoff', 'RunoffGaussianGrid']

logger = logging.getLogger(__name__)


class CellRunoff(NamedTuple):
    """One runoff file read at the weight table's cells, ready for the fused aggregate-and-route kernel."""

    runoff_by_cell: FloatArray  # (n_cells, time) C-order runoff depths, each cell's time series contiguous
    weight: FloatArray  # (n_weights,) proportions already multiplied by the depth unit conversion factor
    scale: FloatArray  # (n_rivers,) catchment areas that turn depths into volumes; EMPTY to keep depths


class RunoffGaussianGrid(Runoff):
    """
    Prepares lateral inflow (vlateral) for routing from gridded runoff depths and a grid weight table.

    Reading and indexing the weight table is the same work for every runoff file, so it is done once when the
    RunoffGaussianGrid is created and reused for all of them. Each runoff file is read at the grid cells the table
    touches and aggregated onto rivers as area weighted depths or volumes in a single pass, written straight into a
    C-order (time, river) array that the routing kernels read without a copy.
    """

    # Weight table, in row order of the table which follows the routing params order
    river_ids: IntArray  # (n_rivers,) river id of each vlateral column
    x_index: IntArray  # (n_cells,) grid x index of each unique cell any catchment touches
    y_index: IntArray  # (n_cells,) grid y index of each unique cell
    indptr: IntArray  # (n_rivers + 1,) sparse row pointers, the weights of river r are indptr[r]:indptr[r + 1]
    cell: IntArray  # (n_weights,) position in x_index and y_index of the cell each weight applies to
    proportion: FloatArray  # (n_weights,) share of the river's catchment area inside the cell
    catchment_area: FloatArray  # (n_rivers,) total catchment area of each river in m²

    # Customizable/optimizable: 32 and 64 were benchmarked as near equivalents.
    rivers_per_block: int = 64
    _buffer: FloatArray  # (time, n_rivers) reused output of vlateral(), grown to the longest file seen

    def __init__(
        self,
        grid_weights_file: PathInput,
        *,
        var_river_id: str = 'river_id',
        var_grid_runoff: str = 'ro',
        var_x: str = 'x',
        var_y: str = 'y',
        var_t: str = 'time',
        grid_accumulation_type: Literal['incremental', 'cumulative'] = 'incremental',
        runoff_depth_unit: str | None = None,
        force_positive_runoff: bool = False,
        force_uniform_timesteps: bool = True,
        as_volumes: bool = False,
    ) -> None:
        """
        Build a RunoffGaussianGrid from the values it needs. ``RunoffGaussianGrid.from_configs`` builds the same thing
        by reading these values off a ``Configs``, which is how a ``Router`` builds one.

        Args:
            grid_weights_file: weight table netCDF produced by ``river_route.runoff.grid_weights``
            var_river_id: name of the river id variable in the weight table
            var_grid_runoff: name of the runoff variable in the gridded runoff files
            var_x: name of the grid x coordinate variable
            var_y: name of the grid y coordinate variable
            var_t: name of the time coordinate variable
            grid_accumulation_type: ``incremental`` runoff per step, or ``cumulative`` totals to difference
            runoff_depth_unit: unit of the runoff depths; None reads the file attributes, else meters
            force_positive_runoff: clip negative runoff depths to zero
            force_uniform_timesteps: resample runoff with irregular timesteps to the first timestep
            as_volumes: prepare volumes (m3) instead of depths (m); routing always uses volumes
        """
        if grid_weights_file is None:
            raise ValueError('grid_weights_file is required to build a RunoffGaussianGrid')
        self.var_runoff = var_grid_runoff
        self.var_x = var_x
        self.var_y = var_y
        self.var_t = var_t
        self.runoff_depth_unit = runoff_depth_unit
        self.cumulative = grid_accumulation_type == 'cumulative'
        self.force_positive_runoff = force_positive_runoff
        self.force_uniform_timesteps = force_uniform_timesteps
        self.as_volumes = as_volumes
        self._read_weights(grid_weights_file, var_river_id)
        self._buffer = np.empty((0, self.river_ids.shape[0]), dtype=np.float32)
        return

    @classmethod
    def from_configs(cls, configs: Configs) -> Self:
        """Build a RunoffGaussianGrid from the runoff options on a ``Configs``. The configs are validated for runoff
        before the weight table is read."""
        if not isinstance(configs, Configs):
            raise TypeError(
                f'from_configs takes a Configs, got {type(configs).__name__}. '
                f'Use Configs(...) or Configs.from_file(path).'
            )
        configs.validate_runoff()
        return cls(
            configs.grid_weights_file,
            var_river_id=configs.var_river_id,
            var_grid_runoff=configs.var_grid_runoff,
            var_x=configs.var_x,
            var_y=configs.var_y,
            var_t=configs.var_t,
            grid_accumulation_type=configs.grid_accumulation_type,
            runoff_depth_unit=configs.runoff_depth_unit,
            force_positive_runoff=configs.force_positive_runoff,
            force_uniform_timesteps=configs.force_uniform_timesteps,
            as_volumes=configs.as_volumes,
        )

    def __repr__(self) -> str:
        return f'{type(self).__name__}(n_rivers={self.river_ids.shape[0]}, n_cells={self.x_index.shape[0]})'

    ################################################
    # Weight table
    ################################################

    def _read_weights(self, grid_weights_file: PathInput, var_river_id: str) -> None:
        """Read a weight table netCDF into the sparse arrays the aggregation kernel consumes."""
        with xr.open_dataset(grid_weights_file) as ds:
            weight_df = ds[[var_river_id, 'x_index', 'y_index', 'proportion', 'area_sqm']].to_dataframe()
        unique_indexes = (
            weight_df[['x_index', 'y_index']].drop_duplicates().reset_index(drop=True).reset_index().astype(int)
        )
        # index already topo sorted
        river_ids = weight_df[[var_river_id]].drop_duplicates().sort_index()[var_river_id].to_numpy()

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

    ################################################
    # Prepare vlateral from runoff
    ################################################

    def vlateral(
        self, runoff_data: PathInput | list[PathInput], thread_pool: Executor | None = None, threads: int = 1
    ) -> tuple[FloatArray, DatetimeArray]:
        """
        Read and aggregate runoff into one reused buffer so that repeated calls allocate nothing.

        Args:
            runoff_data: path(s) to runoff files
            thread_pool: optional thread pool to aggregate on, used as given and never shut down here
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given

        Returns:
            tuple: (vlateral as a C-order (time, n_rivers) array, its time values). The array is a view of the
                reused buffer and is overwritten by the next call.
        """
        runoff, time_index, conversion_factor = self.read_runoff(runoff_data)
        if self._buffer.shape[0] < runoff.shape[0]:
            self._buffer = np.empty((runoff.shape[0], self.river_ids.shape[0]), dtype=np.float32)
        return self.aggregate(
            runoff, time_index, conversion_factor, out=self._buffer, thread_pool=thread_pool, threads=threads
        )

    def reader(
        self, grid_runoff_files: PathList, thread_pool: Executor | None = None, threads: int = 1
    ) -> VlateralGenerator:
        """
        Aggregate gridded runoff files into lateral inflow volumes, one file at a time, with the weight table this
        instance already holds.

        Every file is aggregated into one reused C-order buffer that the kernels read directly: no per-file allocation
        or copy. The yielded array is overwritten by the next file, so a consumer that needs to keep it must copy it.

        Args:
            grid_runoff_files: gridded runoff files to aggregate
            thread_pool: optional thread pool to aggregate on, used as given and never shut down here
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given

        Yields:
            tuple: (dates, vlateral, source_file) per input file
        """
        self.as_volumes = True  # lateral inflow is a volume
        for runoff_file in grid_runoff_files:
            vlateral, dates = self.vlateral(runoff_file, thread_pool=thread_pool, threads=threads)
            yield dates.astype('datetime64[s]'), vlateral.astype(np.float32, copy=False), runoff_file

    def cell_reader(
        self, grid_runoff_files: PathList
    ) -> Iterator[tuple[DatetimeArray, CellRunoff | FloatArray, PathInput]]:
        """
        Read gridded runoff files for the fused kernel, which aggregates and routes in one pass, so no vlateral array
        is built. A file whose timesteps must be resampled cannot be aggregated inside the routing sweep, so it is
        aggregated here and yielded as a vlateral array instead.

        Args:
            grid_runoff_files: gridded runoff files to read

        Yields:
            tuple: (dates, forcing, source_file) per input file, where forcing is a CellRunoff, or a C-order
                (time, n_rivers) vlateral array when the file was resampled
        """
        self.as_volumes = True  # lateral inflow is a volume
        for runoff_file in grid_runoff_files:
            runoff, time_index, conversion_factor = self.read_runoff(runoff_file)
            if self._needs_resampling(time_index):
                vlateral, time_index = self.aggregate(runoff, time_index, conversion_factor)
                yield time_index.astype('datetime64[s]'), vlateral.astype(np.float32, copy=False), runoff_file
                continue
            forcing = CellRunoff(
                runoff_by_cell=self._by_cell(runoff),
                weight=self._weight(conversion_factor),
                scale=self.catchment_area if self.as_volumes else self.catchment_area[:0],
            )
            del runoff
            yield time_index.astype('datetime64[s]'), forcing, runoff_file

    def to_dataset(
        self, runoff_data: PathInput | list[PathInput], thread_pool: Executor | None = None, threads: int = 1
    ) -> xr.Dataset:
        """
        Read and aggregate runoff into a vlateral dataset that can be saved and routed with ``vlateral_files``.

        Args:
            runoff_data: path(s) to runoff files
            thread_pool: optional thread pool to aggregate on, used as given and never shut down here
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given

        Returns:
            xr.Dataset: vlateral with dimensions ``time`` and ``river_id``.
                Contains a single variable ``vlateral`` in meters or m³.
        """
        runoff, time_index, conversion_factor = self.read_runoff(runoff_data)
        vlateral, time_index = self.aggregate(
            runoff, time_index, conversion_factor, thread_pool=thread_pool, threads=threads
        )
        del runoff
        return self._vlateral_dataset(time_index, vlateral, self.river_ids)

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
            out: optional C-order array with at least as many rows as ``runoff`` and one column per river, used as
                the output buffer. Not used when irregular timesteps are resampled.
            thread_pool: optional thread pool, used as given and never shut down here. The rivers are split into
                ``threads`` ranges of similar work that the same kernel aggregates concurrently. Without a pool one
                range covers every river.
            threads: number of river ranges to aggregate concurrently when ``thread_pool`` is given

        Returns:
            tuple: (vlateral as a C-order (time, n_rivers) array, its time values). When ``out`` is used the array is
                a view of its first rows.
        """
        n_steps, n_rivers = runoff.shape[0], self.river_ids.shape[0]
        runoff_by_cell = self._by_cell(runoff)
        weight = self._weight(conversion_factor)
        dtype = np.result_type(runoff.dtype, weight.dtype)
        no_scale = self.catchment_area[:0]

        if self._needs_resampling(time_index):
            vlateral = np.empty((n_steps, n_rivers), dtype=dtype)
            self._aggregate_over_ranges(runoff_by_cell, weight, no_scale, vlateral, thread_pool, threads)
            timestep = int((time_index[1] - time_index[0]) / np.timedelta64(1, 's'))
            logger.warning(f'Time steps are not uniform, resampling to the first timestep: {timestep} seconds')
            df = pd.DataFrame(vlateral, index=time_index, columns=self.river_ids)
            df = df.cumsum().resample(rule=f'{timestep}s').interpolate(method='linear')
            df = self._cumulative_to_incremental(df)
            time_index = df.index.values
            # pandas 3 returns a read-only array from to_numpy, and the conversions below write in place
            vlateral = df.to_numpy(dtype=np.float32, copy=True)
            del df
            vlateral[np.isnan(vlateral)] = 0.0
            if self.as_volumes:
                vlateral *= self.catchment_area[np.newaxis, :]
            return vlateral, time_index

        if out is None:
            out = np.empty((n_steps, n_rivers), dtype=dtype)
        elif out.ndim != 2 or out.shape[0] < n_steps or out.shape[1] != n_rivers or not out.flags.c_contiguous:
            raise ValueError(f'out must be a C-order array of at least ({n_steps}, {n_rivers}), got shape {out.shape}')
        vlateral = out[:n_steps]
        scale = self.catchment_area if self.as_volumes else no_scale
        self._aggregate_over_ranges(runoff_by_cell, weight, scale, vlateral, thread_pool, threads)
        return vlateral, time_index

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
        n_steps = runoff_by_cell.shape[1]
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
                scratch=np.empty((min(self.rivers_per_block, r_stop - r_start), n_steps), dtype=dtype),
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
