import warnings
from pathlib import Path
from typing import Literal, Self

import numpy as np
import pandas as pd

from ..configs.Configs import Configs
from ..types import Float64Array, FloatArray, Int32Array, IntArray, JobBlocks, PathInput
from . import streams

__all__ = [
    'Network',
    'NETWORK_DTYPES',
    'COLUMNS_FOR_PARAMETERS',
    'REQUIRED_COLUMNS',
    'RESTRICTED_LARGEST_X',
    'sub_reach_shares_and_x',
]

# the largest x routed with Configs(restrict_x_range=True). Every coefficient of a river is non-negative only when its
# step is within (1 - 2x) k of its travel time, so an x near 1/2 needs many sub-reaches and subcycles to find one.
RESTRICTED_LARGEST_X = 0.45

# the dtype of every column a network table may have, enforced whenever one is written or stabilized
NETWORK_DTYPES = {
    'riverId': np.int64,
    'nextRiverId': np.int64,
    'muskingumK': np.float32,
    'muskingumX': np.float32,
    'cungeKCoefficient': np.float32,
    'cungeKExponent': np.float32,
    'cungeXCoefficient': np.float32,
    'cungeXExponent': np.float32,
    'cungeMinQ': np.float32,
    'cungeMaxQ': np.float32,
    'riverIndex': np.int32,
    'upstreamCount': np.int32,
    'synthetic': bool,
    'parentRiverId': np.int64,
}
# the columns every network file has: its topology, its Muskingum parameters, and its DFS order
REQUIRED_COLUMNS = ('riverId', 'nextRiverId', 'muskingumK', 'muskingumX', 'riverIndex', 'upstreamCount')
# the columns each dynamic method's parameters are built from, in the order of the fields of its NamedTuple. static
# builds its coefficients from muskingumK and muskingumX
COLUMNS_FOR_PARAMETERS = {
    'static': (),
    'cunge': ('cungeKCoefficient', 'cungeKExponent', 'cungeXCoefficient', 'cungeXExponent', 'cungeMinQ', 'cungeMaxQ'),
}


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
    _routing_blocks: dict[int, tuple[tuple[JobBlocks, ...], Int32Array]]  # threads -> (job_blocks, boundary_target)

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
        """(n,) Muskingum x of each river, its weighting of inflow against outflow, from 0 to 0.5."""
        return self._df['muskingumX'].to_numpy()

    @property
    def upstream_counts(self) -> IntArray:
        """(n,) number of rivers upstream of each river, which in DFS order are the rows immediately before it."""
        return self._df['upstreamCount'].to_numpy()

    @property
    def synthetic(self) -> np.ndarray | None:
        """(n,) True for each sub-reach ``stabilize`` added, or None on a network that was never stabilized."""
        return self._df['synthetic'].to_numpy() if 'synthetic' in self._df.columns else None

    @property
    def parent_river_ids(self) -> IntArray | None:
        """(n,) id of the original river each row was split from on a stabilized network, or None."""
        return self._df['parentRiverId'].to_numpy() if 'parentRiverId' in self._df.columns else None

    @property
    def original_river_ids(self) -> IntArray:
        """Id of each river that is not a synthetic sub-reach, in order: the rows of every routed discharge array."""
        synthetic = self.synthetic
        return self.river_ids if synthetic is None else self.river_ids[~synthetic]

    def columns_for_parameters(self, parameters: str) -> tuple[FloatArray, ...]:
        """
        (n,) arrays of the columns the ``parameters`` are built from, in COLUMNS_FOR_PARAMETERS order.

        Raises:
            ValueError: if the parameters are unknown or the network table lacks one of their columns
        """
        if parameters not in COLUMNS_FOR_PARAMETERS:
            raise ValueError(f'parameters must be one of {sorted(COLUMNS_FOR_PARAMETERS)}, got {parameters!r}')
        columns = COLUMNS_FOR_PARAMETERS[parameters]
        missing = [column for column in columns if column not in self._df.columns]
        if missing:
            raise ValueError(
                f'{self.network_file or "the network table"} is missing the column(s) {", ".join(missing)} that the '
                f'{parameters} parameters read'
            )
        return tuple(self._df[column].to_numpy() for column in columns)

    def __init__(self, network_file: PathInput | pd.DataFrame) -> None:
        """
        Args:
            network_file: the network file parquet, or the network table itself as a DataFrame, which is copied
        """
        if network_file is None:
            raise ValueError('network_file is required to build a Network')
        self._routing_blocks = {}

        if isinstance(network_file, pd.DataFrame):
            self.network_file = None
            self._df = network_file.copy()
        else:
            self.network_file = network_file
            self._df = pd.read_parquet(network_file)
        self.set_connectivity()
        self.check_parameters()
        return

    @classmethod
    def from_configs(cls, configs: Configs) -> Self:
        """
        Build a Network from the network_file of a Configs object, with every x held at RESTRICTED_LARGEST_X or less
        when ``configs.restrict_x_range``.
        """
        if not isinstance(configs, Configs):
            raise TypeError(f'from_configs takes an rr.Configs, got {type(configs).__name__}')
        network = cls(configs.network_file)
        if configs.restrict_x_range:
            network._df['muskingumX'] = np.minimum(network.x, np.float32(RESTRICTED_LARGEST_X))
        return network

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
        draining_into_itself = downstream_idx == upstream_idx
        if draining_into_itself.any():
            river = self.river_ids[upstream_idx[draining_into_itself][0]]
            raise ValueError(
                f'{table} riverId {river} drains into itself, a loop no river network has: make it a basin outlet '
                f'with nextRiverId -1, or remove it, before routing'
            )
        if np.any(downstream_idx < upstream_idx):
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

    def check_parameters(self) -> None:
        """
        Refuse travel times and weights that no routing step can route. k must be positive and finite: a river of zero
        length is an artifact of delineation, with no travel time or storage, whose coefficients are undefined at every
        step, and it is resolved in the topology before routing by joining the rivers upstream of it to the river
        downstream of it. x must be finite and from 0 to 1/2: above 1/2 the scheme amplifies its inflow, and below 0
        the weight of the inflow at the start of each step, c2, is negative at every step shorter than ``-2*k*x``.
        """
        table = self.network_file or 'the network table'
        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        bad_k = ~(np.isfinite(k) & (k > 0))
        if bad_k.any():
            raise ValueError(
                f'{table} muskingumK must be positive and finite, but {int(np.count_nonzero(bad_k))} rivers have '
                f'another k, first riverId {self.river_ids[np.argmax(bad_k)]} with k = {k[np.argmax(bad_k)]}. A river '
                f'of zero length cannot be routed at any step: join the rivers upstream of it to the river downstream '
                f'of it before routing'
            )
        bad_x = ~(np.isfinite(x) & (x >= 0) & (x <= 0.5))
        if bad_x.any():
            raise ValueError(
                f'{table} muskingumX must be finite and from 0 to 0.5, but {int(np.count_nonzero(bad_x))} rivers have '
                f'another x, first riverId {self.river_ids[np.argmax(bad_x)]} with x = {x[np.argmax(bad_x)]}'
            )
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
            tuple: (job_blocks, boundary_target). Each entry of job_blocks is (block_starts, block_stops, block_outlet,
                block_number) for one job, and the last is the main stem. boundary_target gives the river each block's
                outlet drains into, -1 at a basin outlet.
        """
        if threads not in self._routing_blocks:
            blocks = streams.assign_blocks(self.upstream_counts, self.downstream_indices, threads)
            self._routing_blocks[threads] = blocks
        return self._routing_blocks[threads]

    def validity_window(self, sub_reaches: IntArray | None = None) -> tuple[Float64Array, Float64Array]:
        """
        The inclusive range of routing time steps over which every Muskingum coefficient of each river is non-negative,
        ``(2*k*x, 2*k*(1-x))``, the steps at which its routing is valid. Outside it the scheme is still stable and
        conserves mass, but its discharge can go negative. With ``sub_reaches`` the window is that of every one of the
        river's sub-reaches laid out by ``sub_reach_k_and_x``: from c1 of its longest to c3 of its shortest.
        """
        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        if sub_reaches is None:
            return 2 * k * x, 2 * k * (1 - x)
        larger_share, last_share, x = self._sub_reach_shares_and_x(sub_reaches)
        return 2 * larger_share * k * x, 2 * last_share * k * (1 - x)

    def invalid_mask(
        self, dt: float, sub_reaches: IntArray | None = None, subcycles: IntArray | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Per-river masks of the two ways a river's coefficients go negative at ``dt``, for the whole river or, with
        ``sub_reaches``, for any of its sub-reaches. With ``subcycles`` each river is checked at its own shorter
        routing step ``dt / subcycles``.

        Returns:
            tuple: (too_long, too_short). ``too_long`` is ``2*k*x > dt``, which makes c1 negative and is fixable
                by sub-reaches. ``too_short`` is ``2*k*(1-x) < dt``, which makes c3 negative and is fixable by
                subcycles.
        """
        lower, upper = self.validity_window(sub_reaches)
        step = dt if subcycles is None else dt / subcycles
        return lower > step, upper < step

    def sub_reach_k_and_x(self, sub_reaches: IntArray) -> tuple[Float64Array, Float64Array]:
        """
        The travel time and x of every sub-reach of each river split into N = ``sub_reaches`` sub-reaches in series,
        each river's sub-reaches upstream to downstream. Sub-reaches whose travel times are the shares w of k keep the
        river's mean travel time k whatever their x, and all taking ``x = 1/2 - (1/2 - x)/sum(w**2)`` they also keep its
        travel time variance ``k**2 * (1 - 2x)``. Equal shares take ``1/2 - N*(1/2 - x)``, the x Muskingum-Cunge gives
        a reach N times shorter, and they are used wherever it is not negative. Muskingum x is never negative, so where
        it would be, as for two or more sub-reaches of a river with ``x < 1/4``, the sub-reaches take x = 0, linear
        reservoirs, and the first N - 1 take equal shares of k larger than the last's, the shares whose squares sum to
        ``1 - 2x``. Two sub-reaches of a river with ``x = 0.2`` take 0.724 and 0.276 of its k.

        Returns:
            tuple: (k, x), each with a value per sub-reach
        """
        larger_share, last_share, x = self._sub_reach_shares_and_x(sub_reaches)
        share = np.repeat(larger_share, sub_reaches)
        share[np.cumsum(sub_reaches) - 1] = last_share
        return share * np.repeat(self.k.astype(np.float64), sub_reaches), np.repeat(x, sub_reaches)

    def _sub_reach_shares_and_x(self, sub_reaches: IntArray) -> tuple[Float64Array, Float64Array, Float64Array]:
        """(n,) share of k of each of a river's first N - 1 sub-reaches, of its last, and their x; see
        sub_reach_k_and_x."""
        return sub_reach_shares_and_x(sub_reaches, self.x.astype(np.float64))

    def sub_reaches(self, dt: float) -> tuple[IntArray, np.ndarray]:
        """
        The sub-reaches of each river at ``dt``: the fewest sub-reaches, laid out by ``sub_reach_k_and_x``, whose
        coefficients are all non-negative at ``dt``, which is what ``Configs(force_validity=True)`` routes for rivers
        too long for dt. An equal sub-reach of travel time ``k/N`` has a non-negative c1 when
        ``k/N <= dt + (1-2x)*k``, and a sub-reach whose x is 0 has one at any step, so the count is
        ``ceil(k / (dt + (1-2x)*k))``: 2 for every river too long at ``x = 0.2``, and never more than
        ``ceil(1/(1-2x))``. It is the same layout ``stabilize`` uses, so a per-reach state from one lines up with the
        other.

        The shortest sub-reach has a non-negative c3 when ``dt <= 2*k_N*(1 - x_N)``. A river too short for dt cannot
        meet it, since its sub-reaches are shorter still, and for ``x <= 3/8`` every river too long for dt does. A river
        that cannot keeps one reach and is reported as unresolvable.

        Returns:
            tuple: (sub_reaches, resolvable). ``sub_reaches`` is int64 and at least 1 for every river. ``resolvable``
                is True where every coefficient of the sub-reaches is non-negative at dt.
        """
        if dt <= 0:
            raise ValueError(f'dt must be positive, got {dt}')
        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        sub_reaches = np.maximum(1, np.ceil(k / (dt + (1 - 2 * x) * k))).astype(np.int64)
        _, shortest_upper_bound = self.validity_window(sub_reaches)
        resolvable = dt <= shortest_upper_bound
        return np.where(resolvable, sub_reaches, 1), resolvable

    def subcycles(self, dt: float) -> tuple[IntArray, np.ndarray]:
        """
        The subcycles of each river at ``dt``: the fewest equal steps a river is routed in within each step of ``dt``,
        overriding dt with the shorter routing step ``dt / subcycles`` of its own, which is what
        ``Configs(force_validity=True)`` does for rivers too short for dt. The count is ``ceil(dt/(2*k*(1-x)))``:
        the smallest that brings the river's own step under the upper bound ``2*k*(1-x)``. A river that is not too
        short takes 1.

        Its own step can also fall under the lower bound ``2*k*x``, when the window is narrower than one whole count of
        subcycles; such a river keeps one step and is reported as unresolvable. For ``x <= 1/3`` the window always
        holds one.

        Returns:
            tuple: (subcycles, resolvable). ``subcycles`` is int64 and at least 1 for every river. ``resolvable`` is
                True where every coefficient of the river is non-negative at its own step.
        """
        if dt <= 0:
            raise ValueError(f'dt must be positive, got {dt}')
        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        subcycles = np.maximum(1, np.ceil(dt / (2 * k * (1 - x)))).astype(np.int64)
        resolvable = dt / subcycles >= 2 * k * x
        return np.where(resolvable, subcycles, 1), resolvable

    def conditioning(self, dt: float) -> tuple[IntArray, IntArray, np.ndarray]:
        """
        How ``Configs(force_validity=True)`` makes every coefficient of each river non-negative at ``dt``: a
        river too long for dt is split into ``sub_reaches`` sub-reaches laid out by ``sub_reach_k_and_x``, and a river
        too short for it is routed in ``subcycles`` equal steps of its own. No river needs both, since for
        ``x <= 0.5`` a river cannot be too long and too short at once.

        Returns:
            tuple: (sub_reaches, subcycles, resolvable). ``resolvable`` is True where every coefficient is non-negative.
        """
        sub_reaches, sub_reaches_resolvable = self.sub_reaches(dt)
        subcycles, subcycles_resolvable = self.subcycles(dt)
        too_long, too_short = self.invalid_mask(dt)
        resolvable = np.where(too_short, subcycles_resolvable, sub_reaches_resolvable)
        return np.where(too_long, sub_reaches, 1), np.where(too_short, subcycles, 1), resolvable

    def largest_valid_dt(self, period: int) -> int:
        """
        The largest routing time step that divides ``period`` evenly and at which ``conditioning`` can make every
        coefficient of every river non-negative. The fewest sub-reaches a river needs, ``ceil(k / (dt + (1-2x)*k))``,
        only shrinks as dt grows, so the largest such dt also gives the fewest reaches and the fewest steps. Rivers too
        short for dt take subcycles of their own, so they do not pull the whole network down to the dt of its shortest
        river. For ``x <= 1/3`` every river is resolvable at any dt, so this returns ``period`` itself.
        """
        period = int(period)
        if period <= 0:
            raise ValueError(f'period must be a positive integer, got {period}')
        for dt in streams.divisors_of(period)[::-1]:
            _, _, resolvable = self.conditioning(float(dt))
            if resolvable.all():
                return int(dt)
        return 1

    def check_validity(
        self,
        dt: float,
        action: Literal['warn', 'raise', 'ignore'] = 'warn',
        sub_reaches: IntArray | None = None,
        subcycles: IntArray | None = None,
    ) -> None:
        """
        Report rivers with a negative Muskingum coefficient at ``dt`` and take the configured action. With
        ``sub_reaches`` and ``subcycles`` a river is checked as its sub-reaches at its own routing step.

        Every coefficient is non-negative when ``2*k*x <= dt <= 2*k*(1-x)``. Outside that window the scheme is still
        stable and still conserves mass, since c1+c2+c3 == 1 holds for negative coefficients too, but its discharge can
        dip below zero ahead of a rise (c1 < 0) or alternate from step to step and go negative (c3 < 0). The kernels
        write negative discharge as routed.

        Args:
            dt: the routing time step to check against
            action: ``warn`` issues a warning, ``raise`` raises ValueError, ``ignore`` does nothing
            sub_reaches: optional sub-reaches per river, as from ``conditioning``
            subcycles: optional steps of its own per routing step for each river, as from ``conditioning``
        """
        if action == 'ignore':
            return
        too_long, too_short = self.invalid_mask(dt, sub_reaches, subcycles)
        n_long, n_short = int(np.count_nonzero(too_long)), int(np.count_nonzero(too_short))
        if n_long == 0 and n_short == 0:
            return
        message = (
            f'{n_long + n_short} of {self.size} rivers have a negative Muskingum coefficient at dt_routing={dt} s '
            f'({n_long} are too long for it, with c1 < 0, and {n_short} too short, with c3 < 0), so their routing is '
            f'not valid. Every coefficient is non-negative when 2*k*x <= dt_routing <= 2*k*(1-x). Outside that window '
            f'routing is stable and conserves mass, but discharge can dip below zero ahead of a rise or alternate and '
            f'go negative, and it is written as routed. Use Network.invalid_mask to find these rivers, route with '
            f'force_validity=True to route rivers too long in sub-reaches and rivers too short in subcycles, or set '
            f"invalid_coefficients to 'ignore' to silence this."
        )
        if action == 'raise':
            raise ValueError(message)
        warnings.warn(message, stacklevel=2)
        return

    def stabilize(self, dt: float, *, weights: list | None = None) -> Self:
        """
        Stabilize this network in place: every river too long for ``dt`` is replaced by sub-reaches in series whose
        coefficients are all non-negative at that fixed ``dt``. The original river keeps its id and becomes the outlet
        sub-reach, and the added sub-reaches are injected directly upstream of it with ``synthetic`` True and ids
        counting down from -1,000,000. Deleting the synthetic rows gives back the original row order. A network is
        stabilized once: to stabilize it at another dt, build a new Network from the original network file.

        The sub-reaches are those of ``sub_reaches``, laid out by ``sub_reach_k_and_x``, which is what
        ``Configs(force_validity=True)`` routes inside the kernel, so they keep each river's mean travel time k and its
        travel time variance ``k**2 * (1 - 2x)``. A river they cannot make non-negative is kept whole. Rivers left
        with a negative coefficient are not reported here: ``invalid_mask`` finds them, and routing reports them as
        ``invalid_coefficients`` says.

        Args:
            dt: the fixed routing time step at which every coefficient of every sub-reach must be non-negative
            weights: optional explicit sub-reaches, one sequence of relative lengths per river in network table
                order. Each river's k is apportioned in proportion to its weights, which is what matches sub-reaches
                to real segment geometry, and its sub-reaches all take ``x = 1/2 - (1/2 - x)/sum(w**2)`` of their
                shares w of k, which keeps its variance, or 0 where that is negative. A river with a single weight is
                not split. Unlike the sub-reaches laid out for dt, which keep whole a river they cannot make
                non-negative, sub-reaches given here are always built as asked, even for a river with a sub-reach that
                has a negative coefficient.

        Returns:
            Self: this network, stabilized, with a ``synthetic`` flag on every row

        Raises:
            ValueError: if this network is already stabilized or dt is not positive
        """
        if self.synthetic is not None:
            raise ValueError(
                f'{self.network_file or "this network"} is already stabilized; build a new Network from the original '
                f'network file to stabilize it at another dt'
            )
        if dt <= 0:
            raise ValueError(f'dt must be positive, got {dt}')

        k = self.k.astype(np.float64)
        x = self.x.astype(np.float64)
        n_rivers = k.shape[0]
        if weights is None:
            sub_reaches, _ = self.sub_reaches(dt)
            k_reach, x_reach = self.sub_reach_k_and_x(sub_reaches)
        else:
            sub_reaches, k_reach = self._sub_reach_k_from_weights(weights, k)
        first_sub_reach = np.concatenate(([0], np.cumsum(sub_reaches)[:-1]))  # the row of each river's first sub-reach
        share = k_reach / np.repeat(k, sub_reaches)
        squared_shares = np.add.reduceat(share**2, first_sub_reach)
        if weights is not None:
            x_reach = np.repeat(np.maximum(0.0, 0.5 - (0.5 - x) / squared_shares), sub_reaches)

        # each river's block is its synthetic sub-reaches followed by the river itself as the outlet
        df = self._df.iloc[np.repeat(np.arange(n_rivers), sub_reaches)].reset_index(drop=True)
        synthetic = np.ones(len(df), dtype=bool)
        synthetic[first_sub_reach + sub_reaches - 1] = False
        df['synthetic'] = synthetic
        df['parentRiverId'] = df['riverId']
        river_ids = df['riverId'].to_numpy(dtype=np.int64, copy=True)
        river_ids[synthetic] = -1_000_000 - np.arange(np.count_nonzero(synthetic), dtype=np.int64)
        df['riverId'] = river_ids
        # a sub-reach drains into the next row; an outlet drains into the head of its downstream river's block
        next_river_ids = np.append(river_ids[1:], -1)
        down = self.downstream_indices
        next_river_ids[~synthetic] = np.where(down >= 0, river_ids[first_sub_reach][np.clip(down, 0, n_rivers - 1)], -1)
        df['nextRiverId'] = next_river_ids
        # the split travel times are fractions, so they take the table's float dtype even where k was integers. The
        # Muskingum-Cunge k of each sub-reach is its share of the river's, and the 1/2 - x that keeps the river's
        # variance at every discharge is the river's divided by the sum of the squared shares.
        cunge_scales = (('cungeKCoefficient', share), ('cungeXCoefficient', np.repeat(1 / squared_shares, sub_reaches)))
        for column, scale in cunge_scales:
            if column in df.columns:
                df[column] = (df[column].to_numpy(dtype=np.float64) * scale).astype(NETWORK_DTYPES[column])
        df['muskingumK'] = k_reach.astype(NETWORK_DTYPES['muskingumK'])
        df['muskingumX'] = x_reach.astype(NETWORK_DTYPES['muskingumX'])
        # the added sub-reaches shift every row after them, so the DFS order columns are renumbered: the rows upstream
        # of each row run back to the first sub-reach of the first river upstream of the river it belongs to
        first_row_upstream = first_sub_reach[np.arange(n_rivers) - self.upstream_counts]
        rows = np.arange(len(df))
        df['riverIndex'] = (int(self._df['riverIndex'].iloc[0]) + rows).astype(NETWORK_DTYPES['riverIndex'])
        df['upstreamCount'] = (rows - np.repeat(first_row_upstream, sub_reaches)).astype(
            NETWORK_DTYPES['upstreamCount']
        )
        self._df = df
        self._routing_blocks = {}
        self.set_connectivity()
        return self

    def write_stabilized(self, dt: float, path: PathInput | None = None, *, weights: list | None = None) -> Path:
        """
        Stabilize this network for ``dt`` and write it as a network table with the ``synthetic`` column.

        Args:
            dt: the fixed routing time step at which every coefficient of every sub-reach must be non-negative
            path: where to write the parquet. Defaults to ``<network file stem>_stabilized<dt>.parquet`` next to the
                network file this Network was read from.
            weights: passed to ``stabilize``

        Returns:
            Path: the file written
        """
        if path is None:
            if self.network_file is None:
                raise ValueError('this Network was built from a DataFrame, so write_stabilized needs a path')
            network_file = Path(self.network_file)
            path = network_file.with_name(f'{network_file.stem}_stabilized{dt:g}.parquet')
        self.stabilize(dt, weights=weights)
        return self.to_parquet(path)

    @staticmethod
    def _sub_reach_k_from_weights(weights: list, k: Float64Array) -> tuple[IntArray, Float64Array]:
        """The sub-reaches of each river and the k of each, apportioning its k over a sequence of relative lengths."""
        n_rivers = k.shape[0]
        if len(weights) != n_rivers:
            raise ValueError(f'weights must have one sequence per river: got {len(weights)} for {n_rivers} rivers')
        sub_reaches = np.array([len(w) for w in weights], dtype=np.int64)
        if np.any(sub_reaches < 1):
            raise ValueError('every river needs at least one weight')
        flat = np.concatenate([np.asarray(w, dtype=np.float64).ravel() for w in weights])
        if np.any(flat <= 0):
            raise ValueError('weights must be strictly positive')
        totals = np.add.reduceat(flat, np.concatenate(([0], np.cumsum(sub_reaches)[:-1])))
        return sub_reaches, flat * np.repeat(k / totals, sub_reaches)

    def to_parquet(self, path: PathInput) -> Path:
        """
        Write the network table, with every column in its dtype of NETWORK_DTYPES, to a parquet file.

        Args:
            path: the parquet file to write

        Returns:
            Path: the file written
        """
        path = Path(path)
        _enforce_dtypes(self._df).to_parquet(path, index=False)
        return path


def sub_reach_shares_and_x(sub_reaches: IntArray, x: Float64Array) -> tuple[Float64Array, Float64Array, Float64Array]:
    """
    (n,) share of k of each of a river's first N = ``sub_reaches`` - 1 sub-reaches, the share of its last, and the x
    they all take, for rivers of weight x: the layout Network.sub_reach_k_and_x describes, which keeps each river's
    mean travel time and its travel time variance with an x of 0 or more.
    """
    n = sub_reaches.astype(np.float64)
    sub_reach_x = np.maximum(0.0, 0.5 - n * (0.5 - x))
    # the shares w that keep the variance at that x have sum(w**2) = (1 - 2x) / (1 - 2 x_N), which is 1/N where x_N is
    # above 0, the equal shares, and is 1 - 2x where it is held at 0. At x = 1/2 every share keeps it.
    one_minus_two_sub_reach_x = 1 - 2 * sub_reach_x
    squared_shares = np.divide(1 - 2 * x, one_minus_two_sub_reach_x, out=1 / n, where=one_minus_two_sub_reach_x > 0)
    # the larger root p of (N - 1) p**2 + (1 - (N - 1) p)**2 = sum(w**2)
    others = n - 1
    discriminant = np.maximum(others**2 - others * n * (1 - squared_shares), 0.0)
    larger_share = np.divide(others + np.sqrt(discriminant), others * n, out=np.ones_like(n), where=others > 0)
    return larger_share, 1 - others * larger_share, sub_reach_x


def _enforce_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """The network table with every column in NETWORK_DTYPES cast to its dtype."""
    return df.astype({column: dtype for column, dtype in NETWORK_DTYPES.items() if column in df.columns})
