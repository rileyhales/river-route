import warnings
from pathlib import Path
from typing import Literal, Self

import numpy as np
import pandas as pd

from ..configs import Configs
from ..types import Float64Array, FloatArray, Int32Array, IntArray, JobBlocks, PathInput
from . import streams

__all__ = ['Network', 'NETWORK_DTYPES', 'REQUIRED_COLUMNS']

# the dtype of every column a network table may have, enforced whenever one is written or stabilized
NETWORK_DTYPES = {
    'riverId': np.int32,
    'nextRiverId': np.int32,
    'muskingumK': np.float32,
    'muskingumX': np.float32,
    'dynamicAlpha': np.float32,
    'dynamicBeta': np.float32,
    'riverIndex': np.int32,
    'upstreamCount': np.int32,
    'synthetic': bool,
    'parentRiverId': np.int32,
}
# the columns every network file has: its topology, its Muskingum parameters, and its DFS order
REQUIRED_COLUMNS = ('riverId', 'nextRiverId', 'muskingumK', 'muskingumX', 'riverIndex', 'upstreamCount')


class Network:
    """
    A wrapper around the network table, the dataframe that describes the number, topology, and DFS order of the rivers
    and their Muskingum parameters.

    When instantiating this class, the network file is eagerly read and the connectivity is set.
    This may make instantiation have a small but non-negligible time cost.
    """

    # inputs files
    network_file: PathInput | None  # the parquet file the network table was read from, None when given a DataFrame

    # the contents of the network table
    _df: pd.DataFrame
    downstream_indices: Int32Array  # (n,) index of each river's downstream river, -1 at a basin outlet
    _routing_blocks: dict[int, tuple[tuple[JobBlocks, ...], Int32Array]]  # threads -> (job_blocks, cut_target)

    @property
    def size(self) -> int:
        """The number of rivers, one per row of the network table."""
        return self._df.shape[0]

    @property
    def river_ids(self) -> IntArray:
        """(n,) id of each river, in the order of the network table."""
        return self._df['riverId'].to_numpy()

    @property
    def next_river_ids(self) -> IntArray:
        """(n,) id of the river each river drains into, -1 at a basin outlet."""
        return self._df['nextRiverId'].to_numpy()

    @property
    def k(self) -> FloatArray:
        """(n,) Muskingum k of each river, its travel time in seconds."""
        return self._df['muskingumK'].to_numpy()

    @property
    def x(self) -> FloatArray:
        """(n,) Muskingum x of each river, its weighting of inflow against outflow, in [0, 0.5]."""
        return self._df['muskingumX'].to_numpy()

    @property
    def upstream_counts(self) -> IntArray:
        """(n,) number of rivers upstream of each river, which in DFS order are the rows immediately before it."""
        return self._df['upstreamCount'].to_numpy()

    @property
    def dynamicAlpha(self) -> FloatArray | None:
        """(n,) alpha of K = dynamicAlpha * Q ** dynamicBeta for dynamic coefficients, or None without the column."""
        return self._df['dynamicAlpha'].to_numpy() if 'dynamicAlpha' in self._df.columns else None

    @property
    def dynamicBeta(self) -> FloatArray | None:
        """(n,) beta of K = dynamicAlpha * Q ** dynamicBeta for dynamic coefficients, or None without the column."""
        return self._df['dynamicBeta'].to_numpy() if 'dynamicBeta' in self._df.columns else None

    @property
    def synthetic(self) -> np.ndarray | None:
        """(n,) True for each sub-reach ``stabilize`` added, or None on a network that was never stabilized."""
        return self._df['synthetic'].to_numpy() if 'synthetic' in self._df.columns else None

    @property
    def parent_river_ids(self) -> IntArray | None:
        """(n,) id of the original river each row was split from on a stabilized network, or None."""
        return self._df['parentRiverId'].to_numpy() if 'parentRiverId' in self._df.columns else None

    @property
    def supports_dynamic_coefficients(self) -> bool:
        """Whether the network table has the dynamicAlpha and dynamicBeta columns dynamic coefficients need."""
        return 'dynamicAlpha' in self._df.columns and 'dynamicBeta' in self._df.columns

    def __init__(self, network_file: PathInput | pd.DataFrame) -> None:
        """
        Args:
            network_file: the network file parquet, or the network table itself as a DataFrame, which is copied
        """
        if network_file is None:
            raise ValueError('network_file is required to build a Network')
        self._routing_blocks = {}

        # read the parameters
        if isinstance(network_file, pd.DataFrame):
            self.network_file = None
            self._df = network_file.copy()
        else:
            self.network_file = network_file
            self.read_network_file()
        self.set_connectivity()
        return

    def read_network_file(self) -> Self:
        """Read the network table from its parquet file."""
        self._df = pd.read_parquet(self.network_file)
        return self

    @classmethod
    def from_configs(cls, configs: Configs) -> Self:
        """Build a Network from the network_file of a Configs object."""
        if not isinstance(configs, Configs):
            raise TypeError('provide configs is not of type rr.Configs')
        return cls(configs.network_file)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> Self:
        """Instantiate a Network class with the DataFrame in memory already."""
        if not isinstance(df, pd.DataFrame):
            raise TypeError('provide df is not of type pd.DataFrame')
        return cls(df)

    def __repr__(self) -> str:
        return f'{type(self).__name__}(n_rivers={self.size:,}, network_file={self.network_file!r})'

    def set_connectivity(self) -> None:
        """Build the index vectors describing network connectivity from the riverId and nextRiverId columns"""
        table = self.network_file or 'the network table'
        missing = [column for column in REQUIRED_COLUMNS if column not in self._df.columns]
        if missing:
            raise ValueError(f'{table} is missing the required column(s) {", ".join(missing)}')
        n = self.river_ids.shape[0]

        # The lookup is a hash join, not a row at a time dict lookup: at several million rivers the loop it
        # replaces is most of the cost of building a Network. get_indexer needs unique ids to resolve every id to
        # exactly one row, so duplicates are refused first.
        river_index = pd.Index(self.river_ids)
        if not river_index.is_unique:
            raise ValueError(f'{table} has duplicate riverId values')
        has_downstream = self.next_river_ids != -1
        upstream_idx = np.flatnonzero(has_downstream)
        downstream_idx = river_index.get_indexer(self.next_river_ids[has_downstream])

        missing = downstream_idx < 0  # get_indexer returns -1 for an id that is not in the riverId column
        if missing.any():
            unknown = self.next_river_ids[has_downstream][missing]
            raise ValueError(f'{table} nextRiverId {unknown[0]} is not in the riverId column')
        if np.any(downstream_idx <= upstream_idx):
            raise ValueError(f'{table} must be topologically sorted upstream to downstream')

        # the blocks threads route are read from upstreamCount alone, so it must count each river's upstream rivers
        # and those rivers must be the rows immediately before it: each upstream river's rows inside its downstream's
        upstream_counts = self.upstream_counts.astype(np.int64)
        counted = np.bincount(downstream_idx, weights=upstream_counts[upstream_idx] + 1, minlength=n)
        if not np.array_equal(upstream_counts, counted):
            raise ValueError(f'{table} upstreamCount must count every river upstream of each river')
        if np.any(upstream_idx - upstream_counts[upstream_idx] < downstream_idx - upstream_counts[downstream_idx]):
            raise ValueError(f'{table} is not in DFS order: the rivers upstream of a river must be the rows before it')

        # 1D array giving the index of the downstream river in the parameter arrays, -1 if none downstream
        self.downstream_indices = np.full(n, -1, dtype=np.int32)
        self.downstream_indices[upstream_idx] = downstream_idx.astype(np.int32, copy=False)
        return

    def routing_blocks(self, threads: int = 1) -> tuple[tuple[JobBlocks, ...], Int32Array]:
        """
        The blocks of rivers each job routes, packed for ``threads`` jobs by ``streams.assign_blocks``. They are
        derived once per thread count and cached, so repeated simulations over this network never rebuild them.

        Each block is a contiguous range of rows read from the upstreamCount column, so nothing is reordered: the
        forcing, the state, and the routed discharge stay in network file order from end to end. The best blocks
        depend on the thread count, so they are rebuilt for each one. On one thread, the main stem is the whole region.

        Args:
            threads: number of jobs the sub-watershed blocks are packed into

        Returns:
            tuple: (job_blocks, cut_target). Each entry of job_blocks is (block_starts, block_stops, block_outlet,
                block_number) for one job, and the last is the main stem. cut_target gives the river each block's
                outlet drains into, -1 at a basin outlet.
        """
        if threads not in self._routing_blocks:
            blocks = streams.assign_blocks(self.upstream_counts, self.downstream_indices, threads)
            self._routing_blocks[threads] = blocks
        return self._routing_blocks[threads]

    def stability_window(self, substeps: IntArray | None = None) -> tuple[Float64Array, Float64Array]:
        """
        The inclusive range of routing time steps each river is Muskingum-stable over, ``(2*k*x, 2*k*(1-x))``.
        A river is stable for dt exactly when ``dt`` falls inside its own window. With ``substeps`` the window is that
        of one of the river's equal sub-reaches, whose travel time is ``k / substeps``.
        """
        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        if substeps is not None:
            k /= substeps
        return 2 * k * x, 2 * k * (1 - x)

    def unstable_mask(
        self, dt: float, substeps: IntArray | None = None, subcycles: IntArray | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Per-river masks of the two ways a river fails Muskingum stability at ``dt``, for the whole river or, with
        ``substeps``, for each of its equal sub-reaches. With ``subcycles`` each river is checked at its own shorter
        routing step ``dt / subcycles``.

        Returns:
            tuple: (too_long, too_short). ``too_long`` is ``2*k*x > dt``, which makes c1 negative and is fixable
                by substeps. ``too_short`` is ``2*k*(1-x) < dt``, which makes c3 negative and is fixable by
                subcycles.
        """
        lower, upper = self.stability_window(substeps)
        step = dt if subcycles is None else dt / subcycles
        return lower > step, upper < step

    def substeps(self, dt: float) -> tuple[IntArray, np.ndarray]:
        """
        The substeps of each river at ``dt``: the fewest equal sub-reaches a river is divided into so that each routes
        stably at ``dt``, which is what ``Configs(network_type='stabilized')`` routes for rivers too long for dt. The
        count is ``ceil(2*k*x/dt)``: the smallest that brings every sub-reach's travel time ``k/N`` under the upper
        bound ``dt/(2*x)``. It is the same count ``stabilize`` uses, laid out the same way, so a per-reach state from
        one lines up with the other.

        A river that is too short for dt cannot be fixed by substeps, since its sub-reaches are shorter still, and one
        whose window holds no whole count is left alone too. Both keep one reach and are reported as unresolvable.

        Returns:
            tuple: (substeps, resolvable). ``substeps`` is int64 and at least 1 for every river. ``resolvable`` is
                True where the sub-reaches are all stable at dt.
        """
        if dt <= 0:
            raise ValueError(f'dt must be positive, got {dt}')
        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        substeps = np.maximum(1, np.ceil(2 * k * x / dt)).astype(np.int64)
        k_lo = dt / (2 * (1 - x))
        resolvable = k / substeps >= k_lo
        return np.where(resolvable, substeps, 1), resolvable

    def subcycles(self, dt: float) -> tuple[IntArray, np.ndarray]:
        """
        The subcycles of each river at ``dt``: the fewest equal steps a river is routed in within each step of ``dt``,
        overriding dt with the shorter routing step ``dt / subcycles`` of its own, which is what
        ``Configs(network_type='stabilized')`` does for rivers too short for dt. The count is ``ceil(dt/(2*k*(1-x)))``:
        the smallest that brings the river's own step under the upper bound ``2*k*(1-x)``. A river that is not too
        short takes 1.

        Its own step can also fall under the lower bound ``2*k*x``, when the window is narrower than one whole count of
        subcycles; such a river keeps one step and is reported as unresolvable. For ``x <= 1/3`` the window always
        holds one.

        Returns:
            tuple: (subcycles, resolvable). ``subcycles`` is int64 and at least 1 for every river. ``resolvable`` is
                True where the river is stable at its own step.
        """
        if dt <= 0:
            raise ValueError(f'dt must be positive, got {dt}')
        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        upper = 2 * k * (1 - x)
        with np.errstate(divide='ignore'):
            subcycles = np.ceil(np.where(upper > 0, dt / upper, np.inf))
        subcycles = np.maximum(1, np.nan_to_num(subcycles, posinf=np.iinfo(np.int32).max)).astype(np.int64)
        lower = 2 * k * x
        resolvable = dt / subcycles >= lower
        return np.where(resolvable, subcycles, 1), resolvable

    def conditioning(self, dt: float) -> tuple[IntArray, IntArray, np.ndarray]:
        """
        How ``Configs(network_type='stabilized')`` makes each river stable at ``dt``: a river too long for dt is
        divided into ``substeps`` equal sub-reaches, and a river too short for it is routed in ``subcycles`` equal
        steps of its own. No river needs both, since for ``x <= 0.5`` a river cannot be too long and too short at once.

        Returns:
            tuple: (substeps, subcycles, resolvable). ``resolvable`` is True where the river routes stably.
        """
        substeps, substeps_resolvable = self.substeps(dt)
        subcycles, subcycles_resolvable = self.subcycles(dt)
        too_long, too_short = self.unstable_mask(dt)
        resolvable = np.where(too_short, subcycles_resolvable, substeps_resolvable)
        return np.where(too_long, substeps, 1), np.where(too_short, subcycles, 1), resolvable

    def largest_stable_dt(self, period: int) -> int:
        """
        The largest routing time step that divides ``period`` evenly and at which every river can be made stable by
        ``conditioning``. The fewest substeps a river needs, ``ceil(2*k*x/dt)``, only shrinks as dt grows, so the
        largest such dt also gives the fewest reaches and the fewest steps. Rivers too short for dt take subcycles of
        their own, so they do not pull the whole network down to the dt of its shortest river. For ``x <= 1/3`` every
        river is resolvable at any dt, so this returns ``period`` itself.
        """
        period = int(period)
        if period <= 0:
            raise ValueError(f'period must be a positive integer, got {period}')
        for dt in streams.divisors_of(period)[::-1]:
            _, _, resolvable = self.conditioning(float(dt))
            if resolvable.all():
                return int(dt)
        return 1

    def check_stability(
        self,
        dt: float,
        action: Literal['warn', 'raise', 'ignore'] = 'warn',
        substeps: IntArray | None = None,
        subcycles: IntArray | None = None,
    ) -> None:
        """
        Report rivers that are not Muskingum-stable for ``dt`` and take the configured action. With ``substeps`` and
        ``subcycles`` a river counts as stable when its equal sub-reaches are stable at its own routing step.

        Stability requires ``2*k*x <= dt <= 2*k*(1-x)``. Outside that window the solution oscillates and the
        kernels clamp the resulting negative discharges to zero, which does not conserve mass. The c1+c2+c3 == 1
        identity holds for negative coefficients too, so it cannot detect this.

        Args:
            dt: the routing time step to check against
            action: ``warn`` issues a warning, ``raise`` raises ValueError, ``ignore`` does nothing
            substeps: optional sub-reaches per river, as from ``conditioning``
            subcycles: optional steps of its own per routing step for each river, as from ``conditioning``
        """
        if action == 'ignore':
            return
        too_long, too_short = self.unstable_mask(dt, substeps, subcycles)
        n_long, n_short = int(np.count_nonzero(too_long)), int(np.count_nonzero(too_short))
        if n_long == 0 and n_short == 0:
            return
        message = (
            f'{n_long + n_short} of {self.size} rivers are not Muskingum-stable for dt_routing={dt} s ({n_long} need '
            f'a larger dt_routing, {n_short} need a smaller one). Stability requires 2*k*x <= dt_routing <= '
            f'2*k*(1-x) for every river. Routed discharge for these rivers oscillates and negative values are clamped '
            f'to zero, which does not conserve mass. Use Network.unstable_mask to inspect the network or '
            f"Network.stabilize to build a stabilized network, route with network_type='stabilized' to route rivers "
            f'that are too long in substeps and rivers that are too short in subcycles, or set unstable_coefficients '
            f"to 'ignore' to silence this."
        )
        if action == 'raise':
            raise ValueError(message)
        warnings.warn(message, stacklevel=2)
        return

    def stabilize(
        self, dt: float, *, mode: Literal['uniform', 'nonuniform'] = 'uniform', weights: list | None = None
    ) -> Self:
        """
        Stabilize this network in place: every river too long for ``dt`` is replaced by sub-reaches in series that
        each route stably at that fixed ``dt``. The original river keeps its id and becomes the outlet sub-reach,
        and the added sub-reaches are injected directly upstream of it with ``synthetic`` True and ids counting down
        from -1,000,000. Deleting the synthetic rows gives back the original row order. A network that is already
        stabilized is returned unchanged.

        A reach of travel time k divided into N substeps k_1..k_N in series is stable when every piece satisfies
        ``dt/(2*(1-x)) <= k_i <= dt/(2*x)``. Since the pieces must sum to k, N substeps are feasible exactly when
        ``N*k_lo <= k <= N*k_hi`` -- which does not depend on how the pieces are distributed. Uniform and nonuniform
        substeps therefore produce the SAME reach count and fix the same set of rivers; they differ only in where the
        sub-reach boundaries fall.

        Args:
            dt: the fixed routing time step every sub-reach must be stable for
            mode: ``uniform`` gives every sub-reach of a river the same travel time ``k/N``, so boundaries are
                evenly spaced along the reach. ``nonuniform`` packs each river into the fewest pieces of the
                largest stable travel time ``dt/(2*x)`` and leaves the leftover as a shorter final piece, so
                sub-reaches are the same length across the whole network rather than the same count per river.
                A river whose leftover piece would itself be too short falls back to uniform pieces.
            weights: optional explicit substeps, one sequence of relative lengths per river in network table
                order. Each river's k is apportioned in proportion to its weights, which is what matches
                sub-reaches to real segment geometry. Overrides ``mode``. A river with a single weight is not
                divided. Unlike the automatic modes, substeps given here are always built as asked: a river
                whose pieces are not all stable is reported through ``resolvable`` rather than collapsed back to
                a single reach.

        Returns:
            Self: this network, stabilized, with a ``synthetic`` flag on every row
        """
        if self.synthetic is not None:
            return self
        if dt <= 0:
            raise ValueError(f'dt must be positive, got {dt}')
        if mode not in ('uniform', 'nonuniform'):
            raise ValueError(f"mode must be 'uniform' or 'nonuniform', got {mode!r}")

        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        n_rivers = k.shape[0]

        # the travel time window one sub-reach must land in; x == 0 puts no upper bound on a piece
        with np.errstate(divide='ignore'):
            k_hi = np.where(x > 0, dt / (2 * x), np.inf)
        k_lo = dt / (2 * (1 - x))

        if weights is not None:
            substeps, k_reach = self._pieces_from_weights(weights, k)
        else:
            # both modes need the same count: the fewest pieces whose travel time is at or under the upper bound
            ratio = np.where(np.isfinite(k_hi), k / k_hi, 1.0)
            substeps = np.maximum(1, np.ceil(ratio)).astype(np.int64)
            k_reach = (
                np.repeat(k / substeps, substeps)
                if mode == 'uniform'
                else self._pieces_at_target(k, k_hi, k_lo, substeps)
            )

        # a river is fixed only when every one of its pieces lands inside the window
        reach_lo = np.repeat(k_lo, substeps)
        reach_hi = np.repeat(k_hi, substeps)
        piece_ok = (k_reach >= reach_lo) & (k_reach <= reach_hi)
        indptr = np.zeros(n_rivers + 1, dtype=np.int64)
        np.cumsum(substeps, out=indptr[1:])
        resolvable = np.logical_and.reduceat(piece_ok, indptr[:-1]) if k_reach.size else np.ones(0, dtype=bool)
        resolvable &= substeps > 0

        # An unresolvable river is kept whole rather than split into pieces that are still unstable, since the
        # automatic modes exist to fix stability and extra reaches that do not fix it are only extra compute.
        # Explicit weights are a request to represent real geometry, so they are honored as given and the river
        # is only flagged. Both layouts are in river order, so the kept rivers' pieces drop straight into the
        # slots the final layout leaves them.
        if weights is None and not resolvable.all():
            kept_pieces = k_reach[np.repeat(resolvable, substeps)]
            substeps = np.where(resolvable, substeps, 1)
            keep = np.repeat(resolvable, substeps)
            k_reach = np.empty(int(substeps.sum()), dtype=np.float64)
            k_reach[keep] = kept_pieces
            k_reach[~keep] = k[~resolvable]

        # each river's block is its synthetic sub-reaches followed by the river itself as the outlet
        first_piece = np.concatenate(([0], np.cumsum(substeps)[:-1]))  # the row of each river's first piece
        df = self._df.iloc[np.repeat(np.arange(n_rivers), substeps)].reset_index(drop=True)
        synthetic = np.ones(len(df), dtype=bool)
        synthetic[first_piece + substeps - 1] = False
        df['synthetic'] = synthetic
        df['parentRiverId'] = df['riverId']
        df.loc[synthetic, 'riverId'] = -1_000_000 - np.arange(np.count_nonzero(synthetic), dtype=np.int32)
        # a sub-reach drains into the next row; an outlet drains into the head of its downstream river's block
        river_ids = df['riverId'].to_numpy()
        next_river_ids = np.append(river_ids[1:], -1).astype(np.int32)
        down = self.downstream_indices
        next_river_ids[~synthetic] = np.where(down >= 0, river_ids[first_piece][np.clip(down, 0, n_rivers - 1)], -1)
        df['nextRiverId'] = next_river_ids
        # the split travel times are fractions, so they take the table's float dtype even where k was integers
        if 'dynamicAlpha' in df.columns:
            alpha = df['dynamicAlpha'].to_numpy(dtype=np.float64) * k_reach / np.repeat(k, substeps)
            df['dynamicAlpha'] = alpha.astype(NETWORK_DTYPES['dynamicAlpha'])
        df['muskingumK'] = k_reach.astype(NETWORK_DTYPES['muskingumK'])
        # the added sub-reaches shift every row after them, so the DFS order columns are renumbered: the rows upstream
        # of each row run back to the first piece of the first river upstream of the river it belongs to
        first_row_upstream = first_piece[np.arange(n_rivers) - self.upstream_counts]
        rows = np.arange(len(df))
        df['riverIndex'] = (int(self._df['riverIndex'].iloc[0]) + rows).astype(NETWORK_DTYPES['riverIndex'])
        df['upstreamCount'] = (rows - np.repeat(first_row_upstream, substeps)).astype(NETWORK_DTYPES['upstreamCount'])
        self._df = df
        self._routing_blocks = {}
        self.set_connectivity()
        return self

    def write_stabilized(
        self,
        dt: float,
        path: PathInput | None = None,
        *,
        mode: Literal['uniform', 'nonuniform'] = 'uniform',
        weights: list | None = None,
    ) -> Path:
        """
        Stabilize this network for ``dt`` and write it as a network table with the ``synthetic`` column.

        Args:
            dt: the fixed routing time step every sub-reach must be stable for
            path: where to write the parquet. Defaults to ``<network file stem>_stabilized<dt>.parquet`` next to the
                network file this Network was read from.
            mode: passed to ``stabilize``
            weights: passed to ``stabilize``

        Returns:
            Path: the file written
        """
        if path is None:
            if self.network_file is None:
                raise ValueError('this Network was built from a DataFrame, so write_stabilized needs a path')
            network_file = Path(self.network_file)
            path = network_file.with_name(f'{network_file.stem}_stabilized{dt:g}.parquet')
        path = Path(path)
        self.stabilize(dt, mode=mode, weights=weights)
        _enforce_dtypes(self._df).to_parquet(path, index=False)
        return path

    @staticmethod
    def _pieces_at_target(k: Float64Array, k_hi: Float64Array, k_lo: Float64Array, substeps: IntArray) -> Float64Array:
        """
        Pieces of the largest stable travel time with the leftover as a shorter final piece.

        Where the leftover would itself fall below the lower bound the river is evened out into uniform pieces
        instead, which is the only other distribution guaranteed to be inside the window whenever any is.
        """
        total = int(substeps.sum())
        pos_in_block = np.arange(total) - np.repeat(np.concatenate(([0], np.cumsum(substeps)[:-1])), substeps)
        block_size = np.repeat(substeps, substeps)
        target = np.repeat(np.where(np.isfinite(k_hi), k_hi, k), substeps)
        remainder = np.repeat(k - (substeps - 1) * np.where(np.isfinite(k_hi), k_hi, 0.0), substeps)
        packed = np.where(pos_in_block == block_size - 1, remainder, target)
        # even the river out when its leftover piece is below the lower bound
        even_out = np.repeat((k - (substeps - 1) * np.where(np.isfinite(k_hi), k_hi, 0.0)) < k_lo, substeps)
        return np.where(even_out, np.repeat(k / substeps, substeps), packed)

    @staticmethod
    def _pieces_from_weights(weights: list, k: Float64Array) -> tuple[IntArray, Float64Array]:
        """Apportion each river's k over an explicit sequence of relative lengths."""
        n_rivers = k.shape[0]
        if len(weights) != n_rivers:
            raise ValueError(f'weights must have one sequence per river: got {len(weights)} for {n_rivers} rivers')
        substeps = np.array([len(w) for w in weights], dtype=np.int64)
        if np.any(substeps < 1):
            raise ValueError('every river needs at least one weight')
        flat = np.concatenate([np.asarray(w, dtype=np.float64).ravel() for w in weights])
        if np.any(flat <= 0):
            raise ValueError('weights must be strictly positive')
        totals = np.add.reduceat(flat, np.concatenate(([0], np.cumsum(substeps)[:-1])))
        return substeps, flat * np.repeat(k / totals, substeps)

    def to_parquet(self) -> None:
        """Write the network table to the parquet file it was read from, overwriting it."""
        if self.network_file is None:
            raise ValueError('this Network was built from a DataFrame, so it has no network file to write to')
        _enforce_dtypes(self._df).to_parquet(self.network_file, index=False)
        return


def _enforce_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """The network table with every column in NETWORK_DTYPES cast to its dtype."""
    return df.astype({column: dtype for column, dtype in NETWORK_DTYPES.items() if column in df.columns})
