import logging
from dataclasses import dataclass
from typing import Any, Literal, Self

import numpy as np
import pandas as pd

from .._logging import build_logger
from ..configs import Configs
from ..types import FloatArray, IntArray, PathInput
from . import streams

__all__ = ['Network', 'StabilizedNetwork', 'StabilityReport',]

_NULL_LOGGER = logging.getLogger('river_route.network.null')
_NULL_LOGGER.addHandler(logging.NullHandler())

# k and x are stored float32, so a window bound derived from them carries about 1e-7 relative round-off. A reach
# whose travel time is an exact multiple of the bound would otherwise be read as just over it and split one time
# more than it needs. This slack is round-off only: it is far below any difference in k or x that means anything
# physically, and far above float32 epsilon.
_FLOAT32_SLACK = 1e-6


class Network:
    """
    The river network a simulation routes over: the identity, topology, and Muskingum parameters read from a
    parameter table, plus everything derived from them that does not depend on the forcing, the time options, or
    the coefficients.

    Reading and indexing a parameter table is the same work for every simulation over that network, so it is done
    once when the Network is created and reused. The concurrent-routing schedule is derived on demand and cached
    per thread count, so repeated routing, a parameter sweep, or a calibration loop over one network parses and
    partitions it exactly once.

    Nothing here ever reorders a parameter table. The order it is given in is the order it is routed and written
    in, so the forcing, the state, and the routed discharge always line up with the file the user provided.
    """

    # Identity and topology, in parameter table order
    river_ids: IntArray  # (n,) river id of each row
    next_river_ids: IntArray  # (n,) downstream river id, -1 where there is none
    downstream_indices: IntArray  # (n,) int32 position of the downstream river, -1 where there is none

    # Muskingum parameters
    k: FloatArray  # (n,) travel time
    x: FloatArray  # (n,) weighting factor
    alpha: FloatArray | None  # (n,) dynamic coefficients only
    beta: FloatArray | None  # (n,) dynamic coefficients only

    # Provenance, used in error messages
    source: PathInput  # the params file the vectors were read from, used in error messages
    var_river_id: str

    _stored_region: IntArray | None  # the params file 'region' column, when it has one
    _schedules: dict[int, tuple[tuple, IntArray]]  # threads -> (routing_jobs, cut_target)

    def __init__(
        self,
        params_file: PathInput,
        *,
        var_river_id: str = 'river_id',
        coeff: Literal['static', 'dynamic'] = 'static',
        logger: logging.Logger | None = None,
    ) -> None:
        """
        Build a Network from the values it needs. ``Network.from_configs`` builds the same thing by reading these
        values off a ``Configs``, which is how a ``Router`` builds one.

        Args:
            params_file: routing parameter parquet. Required columns are ``var_river_id``, ``next_river_id``,
                ``k``, and ``x``, plus ``alpha`` and ``beta`` when ``coeff='dynamic'``. An optional ``region``
                column from ``network.streams.partition_network`` is reused as the concurrent partition instead of
                deriving one.
            var_river_id: name of the river id column
            coeff: which coefficient scheme the parameters are read for; ``dynamic`` additionally requires the
                alpha and beta columns
            logger: where to send progress messages. Messages are dropped when none is given; ``from_configs``
                builds one from the log options on the Configs.
        """
        if params_file is None:
            raise ValueError('params_file is required to build a Network')
        self.logger = logger if logger is not None else _NULL_LOGGER
        self.var_river_id = var_river_id
        self.coeff = coeff
        self.source = params_file

        self.logger.debug(f'Reading network parameters: {self.source}')
        self._set_vectors(pd.read_parquet(self.source))
        self._set_connectivity()
        self._schedules = {}
        return

    @classmethod
    def from_configs(cls, configs: Configs) -> Self:
        """Build a Network from the options on a ``Configs``: params_file, var_river_id, coeff, and the log
        options the logger is built from."""
        if not isinstance(configs, Configs):
            raise TypeError(
                f'from_configs takes a Configs, got {type(configs).__name__}. '
                f'Use Configs(...) or Configs.from_file(path).'
            )
        return cls(
            configs.params_file,
            var_river_id=configs.var_river_id,
            coeff=configs.coeff,
            logger=build_logger(configs, 'network'),
        )

    def __len__(self) -> int:
        return self.river_ids.shape[0]

    def __repr__(self) -> str:
        return f'{type(self).__name__}(n_rivers={len(self):,}, source={self.source!r})'

    ################################################
    # Reading the parameter table
    ################################################

    def _set_vectors(self, df: pd.DataFrame) -> None:
        required = [self.var_river_id, 'next_river_id', 'k', 'x']
        if self.coeff == 'dynamic':
            required += ['alpha', 'beta']
        missing = [column for column in required if column not in df.columns]
        if missing:
            raise ValueError(f'{self.source} is missing required column(s): {", ".join(missing)}')

        if df[self.var_river_id].duplicated().any():
            raise ValueError(f'{self.source} contains duplicate river IDs.')

        # a stored partition (network.streams.partition_network) is reused as-is; absent, one is derived at setup
        self._stored_region = (
            np.ascontiguousarray(df['region'].to_numpy(copy=False), dtype=np.int64) if 'region' in df.columns else None
        )

        self.river_ids = np.ascontiguousarray(df[self.var_river_id].to_numpy(copy=False), dtype=np.int64)
        self.next_river_ids = np.ascontiguousarray(df['next_river_id'].to_numpy(copy=False), dtype=np.int64)
        self.k = np.ascontiguousarray(df['k'].to_numpy(copy=False), dtype=np.float32)
        self.x = np.ascontiguousarray(df['x'].to_numpy(copy=False), dtype=np.float32)

        self.alpha = None
        self.beta = None
        if self.coeff == 'dynamic':
            self.alpha = np.ascontiguousarray(df['alpha'].to_numpy(copy=False), dtype=np.float32)
            self.beta = np.ascontiguousarray(df['beta'].to_numpy(copy=False), dtype=np.float32)
            if np.any(self.alpha <= 0):
                raise ValueError(f'alpha column in {self.source} must be strictly positive')
        return

    def _set_connectivity(self) -> None:
        """Build the index vectors describing network connectivity from river_ids and next_river_ids"""
        self.logger.debug('Calculating network connectivity vectors')
        n = self.river_ids.shape[0]
        river_index = {int(river_id): idx for idx, river_id in enumerate(self.river_ids.tolist())}

        # 1D array giving the index of the downstream river in the parameter arrays, -1 if none downstream
        self.downstream_indices = np.full(n, -1, dtype=np.int32)
        for upstream_idx, next_river_id in enumerate(self.next_river_ids.tolist()):
            if next_river_id < 0:
                continue
            downstream_idx = river_index.get(int(next_river_id))
            if downstream_idx is None:
                raise ValueError(f'{self.source} next_river_id {next_river_id} is not in the river_id column')
            if downstream_idx <= upstream_idx:
                raise ValueError(f'{self.source} must be topologically sorted upstream to downstream')
            self.downstream_indices[upstream_idx] = downstream_idx

        self.logger.log(logging.INFO, f'Network: {n} river segments')
        return

    ################################################
    # Concurrent routing schedule
    ################################################

    def routing_schedule(self, threads: int = 1, concurrent: bool = True) -> tuple[tuple, IntArray]:
        """
        Build the list of index ranges the kernels sweep. Derived once per thread count and cached, so repeated
        simulations over this network never rebuild the partition.

        The parameter table is always used in the order it is given. Nothing in the routing path reorders a river,
        so the forcing, the state and the routed discharge stay in parameter file order from end to end. Threaded
        routing additionally requires that order to be DFS computation order -- a river following all of its own
        upstream rivers, which makes every subtree a contiguous block that a worker can be handed as a plain index
        range. That is checked against the input and reported if absent, never corrected here.

        The jobs are ordered longest first so the pool packs the very uneven region sizes a river network produces,
        with the main stem last. That orders the work queue only; the rivers inside each block keep their file
        positions. A single-threaded schedule is one job spanning the whole network.

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
        if self._stored_region is not None:
            self.logger.debug('Using the region column stored in the params file')
            region = self._stored_region
        else:
            self.logger.debug('Deriving a network partition for threaded routing')
            region, _ = streams.assign_regions(downstream_index, threads=threads)
        layout = streams.regions_to_layout(region, downstream_index)

        n_regions = layout['n_regions']
        if not n_regions:
            self.logger.warning('The network did not split into any concurrent regions; routing single-threaded')
            self._schedules[key] = single
            return single

        starts, stops = layout['region_starts'], layout['region_stops']
        order = np.argsort(starts - stops)  # submission order for the pool only; river order is untouched
        jobs = [(starts[r : r + 1], stops[r : r + 1], int(layout['region_outlet'][r]), int(r)) for r in order.tolist()]
        jobs.append((layout['stem_starts'], layout['stem_stops'], -1, 0))

        main_stem = int((layout['stem_stops'] - layout['stem_starts']).sum())
        self.logger.log(
            logging.INFO,
            f'Partition: {n_regions} concurrent regions on {threads} threads, '
            f'{main_stem} rivers ({main_stem / n:.2%}) routed sequentially as the main stem',
        )
        schedule = (tuple(jobs), layout['cut_target'])
        self._schedules[key] = schedule
        return schedule

    ################################################
    # Muskingum stability
    ################################################

    def stability_window(self) -> tuple[FloatArray, FloatArray]:
        """
        The inclusive range of routing time steps each river is Muskingum-stable over, ``(2*k*x, 2*k*(1-x))``.
        A river is stable for dt exactly when ``dt`` falls inside its own window.
        """
        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        return 2 * k * x, 2 * k * (1 - x)

    def unstable_mask(self, dt: float) -> tuple[np.ndarray, np.ndarray]:
        """
        Per-river masks of the two ways a river fails Muskingum stability at ``dt``.

        Returns:
            tuple: (too_long, too_short). ``too_long`` is ``2*k*x > dt``, which makes c1 negative and is fixable
                by subdivision. ``too_short`` is ``2*k*(1-x) < dt``, which makes c3 negative and is not.
        """
        lower, upper = self.stability_window()
        return lower > dt, upper < dt

    def substeps_required(self, dt: float) -> IntArray:
        """
        How many times each river would have to be sub-cycled within ``dt`` for its own step to fall inside its
        stability window. 1 means the river needs no temporal refinement. This is the size of the gap described on
        ``StabilizedNetwork``: no routing kernel consumes it yet.
        """
        _, upper = self.stability_window()
        with np.errstate(divide='ignore', invalid='ignore'):
            needed = np.ceil(np.where(upper > 0, dt / upper, np.inf))
        return np.maximum(1, np.nan_to_num(needed, nan=1.0, posinf=np.iinfo(np.int32).max)).astype(np.int64)

    def stability_report(self, dt: float) -> StabilityReport:
        """
        Count how this network fares at ``dt`` and how much bigger the stabilized network would be.

        Returns a ``StabilityReport``, which prints as a readable summary and adds to other reports at the same dt
        so a sweep over many parameter tables accumulates into one total.
        """
        too_long, too_short = self.unstable_mask(dt)
        n_subreaches, resolvable = streams.required_subreaches(self.k, self.x, dt)
        stable = ~(too_long | too_short)
        return StabilityReport(
            dt=dt,
            n_rivers=int(self.river_ids.shape[0]),
            n_stable=int(np.count_nonzero(stable)),
            n_too_long=int(np.count_nonzero(too_long)),
            n_too_short=int(np.count_nonzero(too_short)),
            n_resolved_by_split=int(np.count_nonzero(resolvable & (n_subreaches > 1))),
            n_unresolvable=int(np.count_nonzero(~resolvable)),
            total_subreaches=int(n_subreaches.sum()),
            max_subreaches=int(n_subreaches.max()) if n_subreaches.size else 0,
        )

    def check_stability(self, dt: float, action: Literal['warn', 'raise', 'ignore'] = 'warn') -> StabilityReport | None:
        """
        Report rivers that are not Muskingum-stable for ``dt`` and take the configured action.

        Stability requires ``2*k*x <= dt <= 2*k*(1-x)``. Outside that window the solution oscillates and the
        kernels clamp the resulting negative discharges to zero, which does not conserve mass. The c1+c2+c3 == 1
        identity holds for negative coefficients too, so it cannot detect this.

        Args:
            dt: the routing time step to check against
            action: ``warn`` logs the summary, ``raise`` raises ValueError, ``ignore`` does nothing and returns None

        Returns:
            The report, or None when ``action='ignore'`` or every river is stable.
        """
        if action == 'ignore':
            return None
        report = self.stability_report(dt)
        if report.n_too_long == 0 and report.n_too_short == 0:
            return None
        message = (
            f'{report.n_too_long + report.n_too_short} of {report.n_rivers} rivers are not Muskingum-stable for '
            f'dt_routing={dt} s ({report.n_too_long} need a larger dt_routing, {report.n_too_short} need a '
            f'smaller one). Stability requires 2*k*x <= dt_routing <= 2*k*(1-x) for every river. '
            f'Routed discharge for these rivers oscillates and negative values are clamped to zero, '
            f'which does not conserve mass. Use Network.stability_report to inspect the network or '
            f'Network.stabilize to build a stabilized network, '
            f"or set unstable_coefficients to 'ignore' to silence this."
        )
        if action == 'raise':
            raise ValueError(message)
        self.logger.warning(message)
        return report

    ################################################
    # Synthetic stabilized network
    ################################################

    def stabilize(
        self, dt: float, *, mode: Literal['uniform', 'nonuniform'] = 'uniform', weights: list | None = None
    ) -> StabilizedNetwork:
        """
        Build, in memory, the stabilized network: every reach too long for ``dt`` is replaced by sub-reaches in
        series that each route stably at that fixed ``dt``. The network ends up with more reaches than it started
        with; it is not divided into pieces. Nothing is written and this Network is not modified.

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
            StabilizedNetwork: the flat CSR arrays a routing kernel consumes, plus the per-river ``resolvable``,
                ``too_short``, and ``substeps_required`` flags describing what adding reaches could not stabilize.
        """
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
            # both modes need the same count: the fewest pieces whose travel time is at or under the upper bound.
            # The ratio is nudged down by the round-off slack so a reach that divides the bound exactly is not
            # counted as needing one extra piece.
            ratio = np.where(np.isfinite(k_hi), k / k_hi, 1.0)
            n_sub = np.maximum(1, np.ceil(ratio * (1 - _FLOAT32_SLACK))).astype(np.int64)
            k_reach = np.repeat(k / n_sub, n_sub) if mode == 'uniform' else self._pieces_at_target(k, k_hi, k_lo, n_sub)

        # a river is fixed only when every one of its pieces lands inside the window
        reach_lo = np.repeat(k_lo, n_sub)
        reach_hi = np.repeat(k_hi, n_sub)
        piece_ok = (k_reach >= reach_lo - _tolerance(reach_lo)) & (k_reach <= reach_hi + _tolerance(reach_hi))
        indptr = np.zeros(n_rivers + 1, dtype=np.int64)
        np.cumsum(n_sub, out=indptr[1:])
        resolvable = np.logical_and.reduceat(piece_ok, indptr[:-1]) if k_reach.size else np.ones(0, dtype=bool)
        resolvable &= n_sub > 0

        # An unresolvable river is kept whole rather than split into pieces that are still unstable, since the
        # automatic modes exist to fix stability and extra reaches that do not fix it are only extra compute.
        # Explicit weights are a request to represent real geometry, so they are honored as given and the river
        # is only flagged. Both layouts are in river order, so the kept rivers' pieces drop straight into the
        # slots the final layout leaves them.
        _, too_short = self.unstable_mask(dt)
        if weights is None and not resolvable.all():
            kept_pieces = k_reach[np.repeat(resolvable, n_sub)]
            n_sub = np.where(resolvable, n_sub, 1)
            keep = np.repeat(resolvable, n_sub)
            k_reach = np.empty(int(n_sub.sum()), dtype=np.float64)
            k_reach[keep] = kept_pieces
            k_reach[~keep] = k[~resolvable]

        return self._build_stabilized(
            dt=dt,
            mode='weights' if weights is not None else mode,
            n_sub=n_sub,
            k_reach=k_reach,
            resolvable=resolvable,
            too_short=too_short,
            substeps=self.substeps_required(dt),
        )

    @staticmethod
    def _pieces_at_target(k: FloatArray, k_hi: FloatArray, k_lo: FloatArray, n_sub: IntArray) -> FloatArray:
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

    def _pieces_from_weights(self, weights: list, k: FloatArray) -> tuple[IntArray, FloatArray]:
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

    def _build_stabilized(
        self,
        dt: float,
        mode: str,
        n_sub: IntArray,
        k_reach: FloatArray,
        resolvable: np.ndarray,
        too_short: np.ndarray,
        substeps: IntArray,
    ) -> StabilizedNetwork:
        """Lay the per-river pieces out as the flat, CSR-indexed, topologically sorted reach arrays."""
        n_rivers = self.river_ids.shape[0]
        total = int(n_sub.sum())

        reach_indptr = np.zeros(n_rivers + 1, dtype=np.int64)
        np.cumsum(n_sub, out=reach_indptr[1:])
        group_start = reach_indptr[:-1]
        outlet_index = reach_indptr[1:] - 1

        parent_index = np.repeat(np.arange(n_rivers, dtype=np.int64), n_sub)
        reach_river_id = np.repeat(self.river_ids, n_sub)
        x_reach = np.repeat(self.x.astype(np.float64), n_sub)
        lateral_scale = k_reach / np.repeat(self.k.astype(np.float64), n_sub)

        block_size = np.repeat(n_sub, n_sub)
        pos_in_block = np.arange(total) - np.repeat(group_start, n_sub)
        is_outlet = pos_in_block == (block_size - 1)
        # the outlet is subreach_number 0 and the upstream pieces count up, so re-expanding the same river yields
        # the same (reach_river_id, subreach_number) pairs and saved per-reach state can be matched back
        subreach_number = (block_size - 1 - pos_in_block).astype(np.int32)

        # non-outlet reaches feed the next reach; outlets feed the head of the downstream river's block
        downstream_index = np.empty(total, dtype=np.int64)
        non_outlet = np.nonzero(~is_outlet)[0]
        downstream_index[non_outlet] = non_outlet + 1
        down_river = self.downstream_indices.astype(np.int64)
        downstream_index[is_outlet] = np.where(down_river >= 0, group_start[np.clip(down_river, 0, n_rivers - 1)], -1)

        # the kernel sweeps in array order, so every reach must feed a later one or nothing at all
        if not np.all((downstream_index < 0) | (downstream_index > np.arange(total))):
            raise ValueError(
                'the stabilized network is not topologically sorted upstream-before-downstream; the parameter '
                'table must list every river before its downstream'
            )

        return StabilizedNetwork(
            n_reaches=total,
            n_rivers=n_rivers,
            dt=dt,
            mode=mode,
            k=k_reach.astype(np.float32),
            x=x_reach.astype(np.float32),
            lateral_scale=lateral_scale.astype(np.float32),
            downstream_index=downstream_index.astype(np.int32),
            parent_index=parent_index.astype(np.int32),
            reach_river_id=reach_river_id,
            subreach_number=subreach_number,
            reach_indptr=reach_indptr.astype(np.int32),
            outlet_index=outlet_index.astype(np.int32),
            subdivisions=n_sub.astype(np.int32),
            river_id=self.river_ids,
            resolvable=resolvable,
            too_short=too_short,
            substeps_required=substeps,
        )


@dataclass(frozen=True)
class StabilizedNetwork:
    """
    A stabilized network held in memory as flat, CSR-indexed arrays: the same rivers, with every reach too long
    for ``dt`` replaced by however many sub-reaches in series it takes for each one to route stably at that dt.
    The network gains reaches rather than being divided up, which is why it is a stabilized network and not a
    subdivided one. Nothing is written to disk and the original parameter table is never modified.

    Every river is replaced by a contiguous block of ``subdivisions`` reaches in series, laid out in the original
    topological order: river i owns ``reach_indptr[i]:reach_indptr[i+1]``, the first entry is its inlet and the
    last is its outlet, and the outlet is the only reach whose discharge is reported for river i. Within a block
    each reach feeds the next; an outlet feeds the inlet of its original downstream river's block. That keeps
    ``downstream_index`` a strictly-lower-triangular DAG, so a kernel can sweep the reaches in array order.

    Per-reach ``k`` is the reach's own share of the river's travel time and ``lateral_scale`` is that same share as
    a fraction, so a river's lateral inflow is apportioned by travel time and sums back to the original inflow.
    Uniform subdivision gives every reach ``k/N`` and ``1/N``; nonuniform gives unequal shares that still sum to
    the river's k and to 1.

    A reach is identified across rebuilds by the ``(reach_river_id, subreach_number)`` pair, with subreach_number 0
    at the outlet and counting upstream, so per-reach state can be matched back after the network is rebuilt.

    TODO: rivers whose travel time is TOO SHORT for dt (``2*k*(1-x) < dt``) are a known gap. Adding reaches cannot
      fix them -- splitting a reach shrinks k, which moves a too-short river further outside the stability window.
      The fix is temporal: sub-cycle those reaches ``substeps_required`` times at ``dt / substeps_required`` and
      average the substep outflows. ``substeps_required`` is computed and reported here so the size of the gap is
      visible, but no routing kernel consumes it: ``_kernel_registry`` registers only ``network='standard'``, so
      these rivers are routed as a single unstable reach and are flagged by ``unresolvable``/``too_short``.
    """

    n_reaches: int
    n_rivers: int
    dt: float
    mode: str
    # per reach (n_reaches,)
    k: FloatArray
    x: FloatArray
    lateral_scale: FloatArray
    downstream_index: IntArray
    parent_index: IntArray
    reach_river_id: IntArray
    subreach_number: IntArray
    # per original river (n_rivers,) and the CSR offsets (n_rivers + 1,)
    reach_indptr: IntArray
    outlet_index: IntArray
    subdivisions: IntArray
    river_id: IntArray
    resolvable: np.ndarray
    too_short: np.ndarray
    substeps_required: IntArray

    @property
    def reaches_added(self) -> int:
        """Reaches the stabilized network has beyond the one-per-river it started with."""
        return self.n_reaches - self.n_rivers

    def __repr__(self) -> str:
        return (
            f'StabilizedNetwork(mode={self.mode!r}, dt={self.dt:g}, n_rivers={self.n_rivers:,}, '
            f'n_reaches={self.n_reaches:,}, unresolvable={int((~self.resolvable).sum()):,})'
        )

    def broadcast_state(self, channel_state: FloatArray) -> FloatArray:
        """Seed a per-reach state array from a per-river channel state by copying each river's value into all of
        its sub-reaches, in the order they are laid out here."""
        channel_state = np.asarray(channel_state, dtype=np.float32)
        if channel_state.shape[0] != self.n_rivers:
            raise ValueError(f'channel_state has {channel_state.shape[0]} values for {self.n_rivers} rivers')
        return channel_state[self.parent_index]

    def collapse_state(self, reach_state: FloatArray) -> FloatArray:
        """Reduce a per-reach state array back to one value per river by taking each river's outlet reach."""
        reach_state = np.asarray(reach_state, dtype=np.float32)
        if reach_state.shape[0] != self.n_reaches:
            raise ValueError(f'reach_state has {reach_state.shape[0]} values for {self.n_reaches} reaches')
        return reach_state[self.outlet_index]


@dataclass(frozen=True)
class StabilityReport:
    """
    Counts from a static Muskingum stability check of one network at one routing time step.

    A reach is stable when ``2*k*x <= dt <= 2*k*(1-x)``. Outside that window a Muskingum coefficient is negative,
    the solution oscillates, and the kernels clamp the negative discharges to zero, which does not conserve mass.
    For ``x <= 0.5`` the window is non-empty so the two failure directions are mutually exclusive:

        too_long  (``2*k*x > dt``)      the reach travel time is too long for dt. Splitting the reach into
                                        sub-reaches in series shrinks each k and fixes it, when an integer split
                                        count also satisfies the lower bound. See ``Network.stabilize``.
        too_short (``2*k*(1-x) < dt``)  the reach travel time is too short for dt. Splitting makes this worse.
                                        The fix is temporal (sub-cycling at a smaller dt), which is NOT
                                        implemented; see the ``substeps`` gap noted on ``StabilizedNetwork``.

    Reports add together so a run over many parameter files can accumulate one total, provided every report used
    the same dt.
    """

    dt: float
    n_rivers: int
    n_stable: int
    n_too_long: int
    n_too_short: int
    n_resolved_by_split: int
    n_unresolvable: int
    total_subreaches: int
    max_subreaches: int

    @property
    def subreaches_to_create(self) -> int:
        """Synthetic reaches the fixed network adds on top of the rivers already in the parameter table."""
        return self.total_subreaches - self.n_rivers

    def __add__(self, other: StabilityReport) -> StabilityReport:
        if not isinstance(other, StabilityReport):
            return NotImplemented
        if self.dt != other.dt:
            raise ValueError(f'cannot add stability reports for different dt: {self.dt} and {other.dt}')
        return StabilityReport(
            dt=self.dt,
            n_rivers=self.n_rivers + other.n_rivers,
            n_stable=self.n_stable + other.n_stable,
            n_too_long=self.n_too_long + other.n_too_long,
            n_too_short=self.n_too_short + other.n_too_short,
            n_resolved_by_split=self.n_resolved_by_split + other.n_resolved_by_split,
            n_unresolvable=self.n_unresolvable + other.n_unresolvable,
            total_subreaches=self.total_subreaches + other.total_subreaches,
            max_subreaches=max(self.max_subreaches, other.max_subreaches),
        )

    def __radd__(self, other: Any) -> StabilityReport:
        # so sum() over an iterable of reports works without an explicit start value
        return self if other == 0 else NotImplemented

    def to_dict(self) -> dict[str, float]:
        return {
            'dt': self.dt,
            'n_rivers': self.n_rivers,
            'n_stable': self.n_stable,
            'n_too_long': self.n_too_long,
            'n_too_short': self.n_too_short,
            'n_resolved_by_split': self.n_resolved_by_split,
            'n_unresolvable': self.n_unresolvable,
            'total_subreaches': self.total_subreaches,
            'subreaches_to_create': self.subreaches_to_create,
            'max_subreaches': self.max_subreaches,
        }

    def __str__(self) -> str:
        return (
            f'Muskingum stability for dt={self.dt:g} s\n'
            f'  rivers in:              {self.n_rivers:,}\n'
            f'  already stable:         {self.n_stable:,}\n'
            f'  too long for dt:        {self.n_too_long:,} '
            f'({self.n_resolved_by_split:,} fixable by splitting, up to {self.max_subreaches:,} sub-reaches each)\n'
            f'  too short for dt:       {self.n_too_short:,} (needs substeps, not implemented)\n'
            f'  unresolvable by split:  {self.n_unresolvable:,} (kept as 1 reach, still an error)\n'
            f'  reaches in fixed network: {self.total_subreaches:,}\n'
            f'  sub-reaches to create:    {self.subreaches_to_create:,}'
        )


def _tolerance(bound: FloatArray) -> FloatArray:
    """Float slack for a window comparison, so a piece computed as exactly the bound is not rejected by rounding."""
    return np.where(np.isfinite(bound), np.abs(bound) * _FLOAT32_SLACK, 0.0)
