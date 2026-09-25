# Routing Kernels

Every kernel solves one river's whole time series before moving to the next river. The time series of one river is a
simple recurrence, which is solved eight steps per serial operation instead of one, and gridded runoff can be turned
into lateral inflow in the same pass so no catchment runoff array is built. The cost is that every time step of a file's
forcing must be available at once.

Each routing method is one module that routes a single river, chosen by the `coefficients` config. One numba pass
takes care of the network, the runoff, and the threads for every method, so a method module only says how one river
is routed. Options no method routes yet raise `NotImplementedError` before any runoff is read. Every combination
routes concurrently on a `thread_pool` by splitting the network into regions.

| `coefficients` | Module                     | `forcing`           | `runoff_type`                | `network_type`           |
|----------------|----------------------------|---------------------|------------------------------|--------------------------|
| `static`       | `router/static_muskingum`  | `channel`, `runoff` | `catchment`, `gaussian_grid` | `standard`, `stabilized` |
| `dynamic`      | `router/dynamic_muskingum` | `channel`, `runoff` | `catchment`, `gaussian_grid` | `standard`               |

There is no routing yet for `reduced_gaussian_grid` runoff or the `unit_hydrograph` transform. A stabilized network is
described by the `Layout` the pass hands each river, so a method routes one by following the layout rather than by
being a second method. Each Runoff class's `generator` yields the runoff in the form the pass reads.

## Measured

One year of hourly ERA5 runoff routed over region 6020006540 (303,097 rivers, 8,777 steps) on an Apple M3 Max,
through the fused kernel that aggregates gridded runoff and routes it in one pass. Seconds, minimum of repeated
passes, with the runoff read from a warm page cache. "Write" is a real file on disk, not fsynced.

| Writer | Kernel, 1 thread | Kernel, 12 threads | Write, 1 thread | Write, 12 threads | Total, 12 threads |     File |
|--------|-----------------:|-------------------:|----------------:|------------------:|------------------:|---------:|
| netCDF |             2.74 |               0.32 |            1.75 |              1.72 |              2.37 | 10.64 GB |
| zarr   |             2.78 |               0.32 |            2.15 |              1.93 |              2.63 | 10.64 GB |

The kernel scales about 8.7x from 1 thread to 12, which leaves the writer as most of the job.

Reading the year's runoff adds about 0.18 s warm, or 2 to 3 s from cold storage, and the transpose that prepares the
cell series adds 0.07 s. Both are outside the kernel and writer columns above.

## How routing works

What follows describes the implementation in `river_route/router/_routing_passes.py`, with the recurrence of each
river in `static_muskingum.py`.

### Why one river at a time

Routing every river for one time step before moving on to the next step is one chain of dependent loads and stores,
because in DFS order a river's downstream is usually the next index, so each iteration reads the value the previous
one just added to.
Muskingum couples a river only to its upstreams at the same step and to itself at the previous step, so as long as
every upstream index is below its downstream index, routing rivers one at a time gives the same result:

$$
Q_i (g) = c3_i\,Q_i (g-1) + c4dt_i\,R_i (t) + c1_i\,U_i (g) + c2_i\,U_i (g-1), \qquad U_i (g) = \sum_u Q_u (g)
$$

where $g$ counts routing steps, $t$ is the runoff step that $g$ falls in, $R_i$ is the catchment runoff volume, and
$U_i$ sums the upstream discharges.
Everything but $c3_i Q_i (g-1)$ is known before the step starts, so the serial chain is just that one term, and the
recurrence advances eight steps per serial multiply-add. The trade is that every step of a river's forcing must be
available when that river is routed.

### Passes and regions

The network is routed in passes over blocks of rivers. A region is a contiguous upstream-closed block; its outlet's whole unclamped series is written to its row of the boundary
buffer instead of to an inflow row, which is what lets regions route concurrently. One pass may hold many regions as
separate blocks, keeping per-pass overhead to once per thread rather than once per region. The main stem pass runs
last and adds each boundary row into the inflow of the river it drains into before routing that river. A single
pass over the whole network, with no outlet and no cuts, is the single threaded case.

### Inflow rows

$U_i$ is accumulated in a row of an inflow pool that river `i`'s upstreams add their unclamped series into as they are
routed. Element 0 holds $U_i (-1)$, the sum of the upstreams' initial states, and the rest hold $U_i (g)$. A row is live
only from when a river's first upstream is routed until the river itself is routed, so a small pool of rows is reused
instead of holding an `(n_rivers, n_routing_steps)` array. Two rows follow the pool: a row of zeros that headwaters
read as their inflow, and a sink that basin outlets write into.

### Layout

Everything is river major. Catchment runoff and discharge are C-order `(river, time)` arrays, so each river's series
is one contiguous row, and the kernels read and write each river's row in place. Nothing is transposed on the way in
or out. The runoff arrives one of three ways:

| Runoff                   | Defined in              | How the pass reads it                                                            |
|--------------------------|-------------------------|----------------------------------------------------------------------------------|
| `None`                   |                         | none, channel routing only                                                       |
| `CatchmentRunoffVolumes` | `CatchmentRunoff.py`    | `(river, time)` C-order catchment runoff volumes, each river's row read in place |
| `GridCellRunoff`         | `GaussianGridRunoff.py` | gridded runoff aggregated for 64 rivers at a time into scratch, then routed      |

The one numba pass, `route_scheduled_rivers`, takes each river through stages: finding its catchment runoff,
transforming it, and routing it with the routing method's parameters (`StaticMuskingum` or `DynamicMuskingum`). Each
stage is a function without a body, implemented by numba overloads registered next to the type of argument they read,
so numba compiles a version of the pass for each combination of argument types. A new kind of runoff is a new type
and its overloads of `count_rivers_prepared_together`, `prepare_runoff_of_rivers`, and `get_river_catchment_runoff`.
A new routing method is a new module with its parameters, `NETWORK_TYPES`, `prepare_routing`, and its overload of
`route_river`.

`route_scheduled_rivers` routes each schedule block in groups of rivers whose runoff is prepared together before they
are routed.
Gridded runoff is aggregated for groups of 64 consecutive rivers (`_RIVERS_AGGREGATED_TOGETHER`); any other runoff
needs no preparing, so its whole block is one group. Neighboring catchments share grid cells, so aggregating them back
to back reads each shared cell's series while it is still in cache; routing a river between each aggregation evicts
it. That measured 2-3% faster than aggregating one river at a time on the Amazon.

## Runoff

| You have                                          | What happens                                                      |
|---------------------------------------------------|-------------------------------------------------------------------|
| Gridded runoff and a weight table                 | aggregation is fused into routing, no catchment runoff array      |
| Catchment runoff files, or your own runoff reader | the `(river, time)` array is read in place, one row per river     |

A custom Runoff must yield catchment runoff as `CatchmentRunoffVolumes`: a `(river, time)` array of volumes whose rows
are contiguous. Discharge is always written `(river, time)`, and that C-order array is what a discharge writer
receives.

Inputs are NaN free by the time they reach a kernel: gridded runoff has NaN cells set to zero before it is aggregated,
so a missing cell contributes nothing while the other cells of its catchment still count, and catchment runoff read from
files has NaN values set to zero.
