import warnings
from pathlib import Path
from typing import Literal, Self

import numpy as np
import pandas as pd

from ..configs import Configs
from ..types import Float64Array, FloatArray, Int32Array, IntArray, PathInput
from . import streams

__all__ = ['Network', 'PARAMETER_DTYPES']

# the dtype of every column a parameter table may have, enforced whenever one is read or written
PARAMETER_DTYPES = {
    'river_id': np.int32,
    'next_river_id': np.int32,
    'k': np.float32,
    'x': np.float32,
    'dynamicAlpha': np.float32,
    'dynamicBeta': np.float32,
    'synthetic': bool,
    'parent_river_id': np.int32,
    'group': np.int32,
}

# todo: better ability to identify invalid reaches.
# todo: network should be able to be expanded then write the new network descriptor files to disc


class Network:
    """
    A wrapper around the dataframe of network parameters, or routing parameters, which describes the number and
    topology of the rivers as well as routing parameters.

    When instantiating this class, the network's parameters parquet file are eagerly read and the connectivity is set.
    This may make instantiation have a small but non-negligible time cost.
    """

    # inputs files
    params_file: PathInput | None  # the parquet file the vectors were read from, None when given a DataFrame

    # the contents of the parameter table
    _df: pd.DataFrame
    downstream_indices: Int32Array  # (n,) index of each river's downstream river, -1 at a basin outlet
    _schedules: dict[int, tuple[tuple, Int32Array]]  # threads -> (routing_jobs, cut_target)

    @property
    def size(self) -> int:
        return self._df.shape[0]

    @property
    def river_ids(self) -> IntArray:
        return self._df['river_id'].to_numpy()

    @property
    def next_river_ids(self) -> IntArray:
        return self._df['next_river_id'].to_numpy()

    @property
    def k(self) -> FloatArray:
        return self._df['k'].to_numpy()

    @property
    def x(self) -> FloatArray:
        return self._df['x'].to_numpy()

    @property
    def dynamicAlpha(self) -> FloatArray | None:
        return self._df['dynamicAlpha'].to_numpy() if 'dynamicAlpha' in self._df.columns else None

    @property
    def dynamicBeta(self) -> FloatArray | None:
        return self._df['dynamicBeta'].to_numpy() if 'dynamicBeta' in self._df.columns else None

    @property
    def synthetic(self) -> np.ndarray | None:
        return self._df['synthetic'].to_numpy() if 'synthetic' in self._df.columns else None

    @property
    def parent_river_ids(self) -> IntArray | None:
        return self._df['parent_river_id'].to_numpy() if 'parent_river_id' in self._df.columns else None

    @property
    def groups(self) -> Int32Array | None:
        return self._df['group'].to_numpy() if 'group' in self._df.columns else None

    @property
    def supports_dynamic_coefficients(self) -> bool:
        return 'dynamicAlpha' in self._df.columns and 'dynamicBeta' in self._df.columns

    def __init__(self, params_file: PathInput | pd.DataFrame) -> None:
        """
        Args:
            params_file: the parameter table parquet, or the parameter table itself as a DataFrame, which is copied
        """
        if params_file is None:
            raise ValueError('params_file is required to build a Network')
        self._schedules = {}

        # read the parameters
        if isinstance(params_file, pd.DataFrame):
            self.params_file = None
            self._df = params_file.copy()
        else:
            self.params_file = params_file
            self.read_parameters()
        self.set_connectivity()
        return

    def read_parameters(self) -> Self:
        """Read the parameter table from the parquet file and set the vectors"""
        self._df = pd.read_parquet(self.params_file)
        return self

    @classmethod
    def from_configs(cls, configs: Configs) -> Self:
        """Build a Network from a Configs object"""
        if not isinstance(configs, Configs):
            raise TypeError('provide configs is not of type rr.Configs')
        return cls(configs.params_file)

    def __repr__(self) -> str:
        return f'{type(self).__name__}(n_rivers={self.size:,}, params_file={self.params_file!r})'

    ################################################
    # Reading the parameter table
    ################################################

    def set_connectivity(self) -> None:
        """Build the index vectors describing network connectivity from river_ids and next_river_ids"""
        n = self.river_ids.shape[0]
        table = self.params_file or 'the parameter table'

        # The lookup is a hash join, not a row at a time dict lookup: at several million rivers the loop it
        # replaces is most of the cost of building a Network. get_indexer needs unique ids to resolve every id to
        # exactly one row, so duplicates are refused first.
        river_index = pd.Index(self.river_ids)
        if not river_index.is_unique:
            raise ValueError(f'{table} has duplicate river_id values')
        has_downstream = self.next_river_ids != -1
        upstream_idx = np.flatnonzero(has_downstream)
        downstream_idx = river_index.get_indexer(self.next_river_ids[has_downstream])

        missing = downstream_idx < 0  # get_indexer returns -1 for an id that is not in the river_id column
        if missing.any():
            unknown = self.next_river_ids[has_downstream][missing]
            raise ValueError(f'{table} next_river_id {unknown[0]} is not in the river_id column')
        if np.any(downstream_idx <= upstream_idx):
            raise ValueError(f'{table} must be topologically sorted upstream to downstream')

        # 1D array giving the index of the downstream river in the parameter arrays, -1 if none downstream
        self.downstream_indices = np.full(n, -1, dtype=np.int32)
        self.downstream_indices[upstream_idx] = downstream_idx.astype(np.int32, copy=False)
        return

    ################################################
    # Concurrent routing schedule
    ################################################

    def routing_schedule(self, threads: int = 1, concurrent: bool = True) -> tuple[tuple, Int32Array]:
        """
        Build the index ranges of the rivers each routing pass routes. Derived once per thread count and cached, so
        repeated simulations over this network never rebuild the partition.

        The parameter table is always used in the order it is given. Nothing in the routing path reorders a river,
        so the forcing, the state and the routed discharge stay in parameter file order from end to end. Threaded
        routing additionally requires that order to be DFS computation order -- a river following all of its own
        upstream rivers, which makes every subtree a contiguous block that a worker can be handed as a plain index
        range. That is checked against the input and reported if absent, never corrected here.

        The jobs are ordered longest first so the pool packs the very uneven region sizes a river network produces,
        with the main stem last. That orders the work queue only; the rivers inside each block keep their file
        positions. A single-threaded schedule is one job spanning the whole network.

        The regions come from the params file's ``group`` column when it has one, and otherwise from
        ``streams.assign_regions``, which sizes them for ``threads``: it searches a range of caps on the largest
        region and keeps whichever partition packs onto that many workers fastest. A thread-independent cut such as
        ``recommend_compute_groups`` puts a hard floor under the wall clock -- one oversized region no thread count
        can split -- so the partition has to be rebuilt per thread count rather than derived once from topology.

        Args:
            threads: worker count the partition is sized for
            concurrent: False forces the single-job schedule regardless of ``threads``, for a caller that has no
                thread pool to run the regions on

        Returns:
            tuple: (routing_jobs, cut_target). Each job is (block_starts, block_stops, outlet, region); cut_target
                gives the river each region's outlet drains into, -1 at a basin outlet.
        """
        key = threads if concurrent else 1
        if key in self._schedules:
            return self._schedules[key]

        n = self.river_ids.shape[0]
        whole = (np.array([0], dtype=np.int32), np.array([n], dtype=np.int32), -1, 0)
        single = ((whole,), np.zeros(0, dtype=np.int32))
        if not concurrent or threads < 2:
            self._schedules[key] = single
            return single

        downstream_index = self.downstream_indices.astype(np.int64)
        region = self.groups
        if region is None:
            # sized for this thread count; nothing is cached on the frame, so a later call for a different thread
            # count is free to cut the network differently
            region, _ = streams.assign_regions(downstream_index, threads=threads)
        layout = streams.regions_to_layout(np.ascontiguousarray(region, dtype=np.int32), downstream_index)

        n_regions = layout['n_regions']
        if not n_regions:
            self._schedules[key] = single
            return single

        starts, stops = layout['region_starts'], layout['region_stops']
        order = np.argsort(starts - stops)  # submission order for the pool only; river order is untouched
        jobs = [(starts[r : r + 1], stops[r : r + 1], int(layout['region_outlet'][r]), int(r)) for r in order.tolist()]
        jobs.append((layout['stem_starts'], layout['stem_stops'], -1, 0))
        schedule = (tuple(jobs), layout['cut_target'])
        self._schedules[key] = schedule
        return schedule

    def recommend_compute_groups(self) -> IntArray:
        """
        Recommend a group for every river from the topology alone: the tributary of a basin's main stem it belongs
        to, and -1 on the main stem, as ``streams.tributary_groups`` finds them. Offered for writing a ``group``
        column into a params file; routing does not call it, since these groups are the same however many threads
        route them and a network's largest tributary is far too big to be one worker's share.
        """
        return streams.tributary_groups(self.downstream_indices)

    ################################################
    # Muskingum stability
    ################################################

    def stability_window(self, subdivisions: IntArray | None = None) -> tuple[Float64Array, Float64Array]:
        """
        The inclusive range of routing time steps each river is Muskingum-stable over, ``(2*k*x, 2*k*(1-x))``.
        A river is stable for dt exactly when ``dt`` falls inside its own window. With ``subdivisions`` the window
        is that of one of the river's equal sub-reaches, whose travel time is ``k / subdivisions``.
        """
        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        if subdivisions is not None:
            k /= subdivisions
        return 2 * k * x, 2 * k * (1 - x)

    def unstable_mask(
        self, dt: float, subdivisions: IntArray | None = None, substeps: IntArray | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Per-river masks of the two ways a river fails Muskingum stability at ``dt``, for the whole river or, with
        ``subdivisions``, for each of its equal sub-reaches. With ``substeps`` each river is checked at its own step
        ``dt / substeps``.

        Returns:
            tuple: (too_long, too_short). ``too_long`` is ``2*k*x > dt``, which makes c1 negative and is fixable
                by subdivision. ``too_short`` is ``2*k*(1-x) < dt``, which makes c3 negative and is fixable by
                substeps.
        """
        lower, upper = self.stability_window(subdivisions)
        step = dt if substeps is None else dt / substeps
        return lower > step, upper < step

    def subdivisions(self, dt: float) -> tuple[IntArray, np.ndarray]:
        """
        The fewest equal sub-reaches each river is split into to route stably at ``dt``, which is what
        ``Configs(network_type='stabilized')`` routes. The count is ``ceil(2*k*x/dt)``: the smallest that
        brings every sub-reach's travel time ``k/N`` under the upper bound ``dt/(2*x)``. It is the same count
        ``stabilize`` uses, laid out the same way, so a per-reach state from one lines up with the other.

        A river that is too short for dt cannot be fixed by splitting, since splitting shortens it further, and one
        whose window holds no whole count is left alone too. Both keep one reach and are reported as unresolvable.

        Returns:
            tuple: (subdivisions, resolvable). ``subdivisions`` is int64 and at least 1 for every river.
                ``resolvable`` is True where the sub-reaches are all stable at dt.
        """
        if dt <= 0:
            raise ValueError(f'dt must be positive, got {dt}')
        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        n_sub = np.maximum(1, np.ceil(2 * k * x / dt)).astype(np.int64)
        k_lo = dt / (2 * (1 - x))
        resolvable = k / n_sub >= k_lo
        return np.where(resolvable, n_sub, 1), resolvable

    def substeps(self, dt: float) -> tuple[IntArray, np.ndarray]:
        """
        The fewest equal substeps each river is sub-cycled in within ``dt`` to route stably, which is what
        ``Configs(network_type='stabilized')`` does for rivers too short for dt. The count is
        ``ceil(dt/(2*k*(1-x)))``: the smallest that brings the river's own step ``dt/m`` under the upper bound
        ``2*k*(1-x)``. A river that is not too short takes 1.

        A step can also fall under the lower bound ``2*k*x``, when the window is narrower than one whole count of
        substeps; such a river keeps one step and is reported as unresolvable. For ``x <= 1/3`` the window always
        holds one.

        Returns:
            tuple: (substeps, resolvable). ``substeps`` is int64 and at least 1 for every river. ``resolvable`` is
                True where the river is stable at its own step.
        """
        if dt <= 0:
            raise ValueError(f'dt must be positive, got {dt}')
        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        upper = 2 * k * (1 - x)
        with np.errstate(divide='ignore'):
            m = np.ceil(np.where(upper > 0, dt / upper, np.inf))
        m = np.maximum(1, np.nan_to_num(m, posinf=np.iinfo(np.int32).max)).astype(np.int64)
        lower = 2 * k * x
        resolvable = dt / m >= lower
        return np.where(resolvable, m, 1), resolvable

    def conditioning(self, dt: float) -> tuple[IntArray, IntArray, np.ndarray]:
        """
        How ``Configs(network_type='stabilized')`` makes each river stable at ``dt``: a river too long for dt
        is split into ``subdivisions`` equal sub-reaches, and a river too short for it is sub-cycled in ``substeps``
        equal steps. No river needs both, since for ``x <= 0.5`` a river cannot be too long and too short at once.

        Returns:
            tuple: (subdivisions, substeps, resolvable). ``resolvable`` is True where the river routes stably.
        """
        subdivisions, split_ok = self.subdivisions(dt)
        substeps, cycle_ok = self.substeps(dt)
        too_long, too_short = self.unstable_mask(dt)
        substeps = np.where(too_short, substeps, 1)
        resolvable = np.where(too_short, cycle_ok, split_ok)
        return np.where(too_long, subdivisions, 1), substeps, resolvable

    def largest_stable_dt(self, period: int) -> int:
        """
        The largest routing time step that divides ``period`` evenly and at which every river can be made stable by
        ``conditioning``. The fewest sub-reaches a river needs, ``ceil(2*k*x/dt)``, only shrinks as dt grows, so the
        largest such dt also gives the fewest reaches and the fewest steps. Rivers too short for dt are sub-cycled
        in their own steps, so they do not pull the whole network down to the dt of its shortest river. For
        ``x <= 1/3`` every river is resolvable at any dt, so this returns ``period`` itself.
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
        subdivisions: IntArray | None = None,
        substeps: IntArray | None = None,
    ) -> None:
        """
        Report rivers that are not Muskingum-stable for ``dt`` and take the configured action. With
        ``subdivisions`` and ``substeps`` a river counts as stable when its equal sub-reaches are at its own step.

        Stability requires ``2*k*x <= dt <= 2*k*(1-x)``. Outside that window the solution oscillates and the
        kernels clamp the resulting negative discharges to zero, which does not conserve mass. The c1+c2+c3 == 1
        identity holds for negative coefficients too, so it cannot detect this.

        Args:
            dt: the routing time step to check against
            action: ``warn`` issues a warning, ``raise`` raises ValueError, ``ignore`` does nothing
            subdivisions: optional sub-reach count per river, as from ``conditioning``
            substeps: optional substep count per river, as from ``conditioning``
        """
        if action == 'ignore':
            return
        too_long, too_short = self.unstable_mask(dt, subdivisions, substeps)
        n_long, n_short = int(np.count_nonzero(too_long)), int(np.count_nonzero(too_short))
        if n_long == 0 and n_short == 0:
            return
        message = (
            f'{n_long + n_short} of {self.size} rivers are not Muskingum-stable for dt_routing={dt} s ({n_long} need '
            f'a larger dt_routing, {n_short} need a smaller one). Stability requires 2*k*x <= dt_routing <= '
            f'2*k*(1-x) for every river. Routed discharge for these rivers oscillates and negative values are clamped '
            f'to zero, which does not conserve mass. Use Network.unstable_mask to inspect the network or '
            f"Network.stabilize to build a stabilized network, route with network_type='stabilized' "
            f"to split rivers that are too long, or set unstable_coefficients to 'ignore' to silence this."
        )
        if action == 'raise':
            raise ValueError(message)
        warnings.warn(message, stacklevel=2)
        return

    ################################################
    # Synthetic stabilized network
    ################################################

    def stabilize(
        self, dt: float, *, mode: Literal['uniform', 'nonuniform'] = 'uniform', weights: list | None = None
    ) -> Self:
        """
        Stabilize this network in place: every reach too long for ``dt`` is replaced by sub-reaches in series that
        each route stably at that fixed ``dt``. The original river keeps its id and becomes the outlet sub-reach,
        and the added sub-reaches are injected directly upstream of it with ``synthetic`` True and ids counting down
        from -1,000,000. Deleting the synthetic rows gives back the original row order. A network that is already
        stabilized is returned unchanged.

        A reach of travel time k split into pieces k_1..k_N in series is stable when every piece satisfies
        ``dt/(2*(1-x)) <= k_i <= dt/(2*x)``. Since the pieces must sum to k, a split of size N is feasible exactly
        when ``N*k_lo <= k <= N*k_hi`` -- which does not depend on how the pieces are distributed. Uniform and
        nonuniform subdivision therefore produce the SAME reach count and fix the same set of rivers; they differ
        only in where the sub-reach boundaries fall.

        Args:
            dt: the fixed routing time step every sub-reach must be stable for
            mode: ``uniform`` gives every sub-reach of a river the same travel time ``k/N``, so boundaries are
                evenly spaced along the reach. ``nonuniform`` packs each river into the fewest pieces of the
                largest stable travel time ``dt/(2*x)`` and leaves the leftover as a shorter final piece, so
                sub-reaches are the same length across the whole network rather than the same count per river.
                A river whose leftover piece would itself be too short falls back to uniform pieces.
            weights: optional explicit subdivision, one sequence of relative lengths per river in parameter table
                order. Each river's k is apportioned in proportion to its weights, which is what matches
                sub-reaches to real segment geometry. Overrides ``mode``. A river with a single weight is not
                split. Unlike the automatic modes, a subdivision given here is always built as asked: a river
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
            n_sub, k_reach = self._pieces_from_weights(weights, k)
        else:
            # both modes need the same count: the fewest pieces whose travel time is at or under the upper bound
            ratio = np.where(np.isfinite(k_hi), k / k_hi, 1.0)
            n_sub = np.maximum(1, np.ceil(ratio)).astype(np.int64)
            k_reach = np.repeat(k / n_sub, n_sub) if mode == 'uniform' else self._pieces_at_target(k, k_hi, k_lo, n_sub)

        # a river is fixed only when every one of its pieces lands inside the window
        reach_lo = np.repeat(k_lo, n_sub)
        reach_hi = np.repeat(k_hi, n_sub)
        piece_ok = (k_reach >= reach_lo) & (k_reach <= reach_hi)
        indptr = np.zeros(n_rivers + 1, dtype=np.int64)
        np.cumsum(n_sub, out=indptr[1:])
        resolvable = np.logical_and.reduceat(piece_ok, indptr[:-1]) if k_reach.size else np.ones(0, dtype=bool)
        resolvable &= n_sub > 0

        # An unresolvable river is kept whole rather than split into pieces that are still unstable, since the
        # automatic modes exist to fix stability and extra reaches that do not fix it are only extra compute.
        # Explicit weights are a request to represent real geometry, so they are honored as given and the river
        # is only flagged. Both layouts are in river order, so the kept rivers' pieces drop straight into the
        # slots the final layout leaves them.
        if weights is None and not resolvable.all():
            kept_pieces = k_reach[np.repeat(resolvable, n_sub)]
            n_sub = np.where(resolvable, n_sub, 1)
            keep = np.repeat(resolvable, n_sub)
            k_reach = np.empty(int(n_sub.sum()), dtype=np.float64)
            k_reach[keep] = kept_pieces
            k_reach[~keep] = k[~resolvable]

        # each river's block is its synthetic sub-reaches followed by the river itself as the outlet
        group_start = np.concatenate(([0], np.cumsum(n_sub)[:-1]))
        df = self._df.iloc[np.repeat(np.arange(n_rivers), n_sub)].reset_index(drop=True)
        synthetic = np.ones(len(df), dtype=bool)
        synthetic[group_start + n_sub - 1] = False
        df['synthetic'] = synthetic
        df['parent_river_id'] = df['river_id']
        df.loc[synthetic, 'river_id'] = -1_000_000 - np.arange(np.count_nonzero(synthetic), dtype=np.int32)
        # a sub-reach drains into the next row; an outlet drains into the head of its downstream river's block
        river_ids = df['river_id'].to_numpy()
        next_river_ids = np.append(river_ids[1:], -1).astype(np.int32)
        down = self.downstream_indices
        next_river_ids[~synthetic] = np.where(down >= 0, river_ids[group_start][np.clip(down, 0, n_rivers - 1)], -1)
        df['next_river_id'] = next_river_ids
        if 'dynamicAlpha' in df.columns:
            df['dynamicAlpha'] *= (k_reach / np.repeat(k, n_sub)).astype(df['dynamicAlpha'].dtype)
        df['k'] = k_reach.astype(df['k'].dtype)
        self._df = df
        self._schedules = {}
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
        Stabilize this network for ``dt`` and write it as a parameter table with the ``synthetic`` column.

        Args:
            dt: the fixed routing time step every sub-reach must be stable for
            path: where to write the parquet. Defaults to ``<params stem>_stabilized<dt>.parquet`` next to the
                params file this Network was read from.
            mode: passed to ``stabilize``
            weights: passed to ``stabilize``

        Returns:
            Path: the file written
        """
        if path is None:
            if self.params_file is None:
                raise ValueError('this Network was built from a DataFrame, so write_stabilized needs a path')
            params_file = Path(self.params_file)
            path = params_file.with_name(f'{params_file.stem}_stabilized{dt:g}.parquet')
        path = Path(path)
        self.stabilize(dt, mode=mode, weights=weights)
        _enforce_dtypes(self._df).to_parquet(path, index=False)
        return path

    @staticmethod
    def _pieces_at_target(k: Float64Array, k_hi: Float64Array, k_lo: Float64Array, n_sub: IntArray) -> Float64Array:
        """
        Pieces of the largest stable travel time with the leftover as a shorter final piece.

        Where the leftover would itself fall below the lower bound the river is evened out into uniform pieces
        instead, which is the only other distribution guaranteed to be inside the window whenever any is.
        """
        total = int(n_sub.sum())
        pos_in_block = np.arange(total) - np.repeat(np.concatenate(([0], np.cumsum(n_sub)[:-1])), n_sub)
        block_size = np.repeat(n_sub, n_sub)
        target = np.repeat(np.where(np.isfinite(k_hi), k_hi, k), n_sub)
        remainder = np.repeat(k - (n_sub - 1) * np.where(np.isfinite(k_hi), k_hi, 0.0), n_sub)
        packed = np.where(pos_in_block == block_size - 1, remainder, target)
        # even the river out when its leftover piece is below the lower bound
        even_out = np.repeat((k - (n_sub - 1) * np.where(np.isfinite(k_hi), k_hi, 0.0)) < k_lo, n_sub)
        return np.where(even_out, np.repeat(k / n_sub, n_sub), packed)

    @staticmethod
    def _pieces_from_weights(weights: list, k: Float64Array) -> tuple[IntArray, Float64Array]:
        """Apportion each river's k over an explicit sequence of relative lengths."""
        n_rivers = k.shape[0]
        if len(weights) != n_rivers:
            raise ValueError(f'weights must have one sequence per river: got {len(weights)} for {n_rivers} rivers')
        n_sub = np.array([len(w) for w in weights], dtype=np.int64)
        if np.any(n_sub < 1):
            raise ValueError('every river needs at least one weight')
        flat = np.concatenate([np.asarray(w, dtype=np.float64).ravel() for w in weights])
        if np.any(flat <= 0):
            raise ValueError('weights must be strictly positive')
        totals = np.add.reduceat(flat, np.concatenate(([0], np.cumsum(n_sub)[:-1])))
        return n_sub, flat * np.repeat(k / totals, n_sub)

    def to_parquet(self) -> None:
        """Write the parameter table to the parquet file it was read from, overwriting it."""
        if self.params_file is None:
            raise ValueError('this Network was built from a DataFrame, so it has no params file to write to')
        _enforce_dtypes(self._df).to_parquet(self.params_file, index=False)
        return


def _enforce_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """The parameter table with every column in PARAMETER_DTYPES cast to its dtype."""
    return df.astype({column: dtype for column, dtype in PARAMETER_DTYPES.items() if column in df.columns})
