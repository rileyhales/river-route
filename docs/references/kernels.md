# Routing Kernels

Every kernel solves one river's whole time series before moving to the next river. The time series of one river is a
simple recurrence, which is solved eight steps per serial operation instead of one, and gridded runoff can be turned
into lateral inflow in the same pass so no vlateral array is built. The cost is that every time step of a file's
forcing must be available at once.

A kernel is resolved from four selector keys, `coeff`, `forcing`, `transform`, and `network_conditioning`. All of
them route concurrently on a `thread_pool` by splitting the network into regions.

## Measured

One year of hourly ERA5 runoff routed over region 6020006540 (303,097 rivers, 8,777 steps) on an Apple M3 Max,
through the fused kernel that aggregates gridded runoff and routes it in one pass. Seconds, minimum of repeated
passes, with the runoff read from a warm page cache. "Write" is a real file on disk, not fsynced.

| Writer | `discharge_dtype` | Kernel, 1 thread | Kernel, 12 threads | Write, 1 thread | Write, 12 threads | Total, 12 threads | File |
|--------|-------------------|-----------------:|-------------------:|----------------:|------------------:|------------------:|-----:|
| zarr   | `float16`         |             2.45 |           **0.28** |            1.24 |          **0.98** |          **1.58** | 5.33 GB |
| netCDF | `float32`         |             2.74 |               0.32 |            1.75 |              1.72 |              2.37 | 10.64 GB |
| zarr   | `float32`         |             2.78 |               0.32 |            2.15 |              1.93 |              2.63 | 10.64 GB |
| netCDF | `float16`         |             2.44 |               0.25 |            3.39 |              3.57 |              4.11 | 10.64 GB |

The kernel scales about 8.7x from 1 thread to 12, which leaves the writer as most of the job. `float16` halves the
bytes written wherever the format stores it natively, which zarr and parquet do; netCDF has no half type, so it
widens back to float32 and pays for the pass without shrinking the file.

Reading the year's runoff adds about 0.18 s warm, or 2 to 3 s from cold storage, and the transpose that prepares the
cell series adds 0.07 s. Both are outside the kernel and writer columns above.

## How routing works

What follows describes the implementation in `river_route/router/_river_kernels.py`.

### Why one river at a time

Sweeping the whole network once per time step is one chain of dependent loads and stores, because in DFS order a
river's downstream is usually the next index, so each iteration reads the value the previous one just added to.
Muskingum couples a river only to its upstreams at the same step and to itself at the previous step, so as long as
every upstream index is below its downstream index, routing rivers one at a time gives the same result:

$$
Q_i(g) = c3_i\,Q_i(g-1) + c4dt_i\,vlateral_i(t) + c1_i\,U_i(g) + c2_i\,U_i(g-1), \qquad U_i(g) = \sum_u Q_u(g)
$$

where $g$ counts routing steps, $t$ is the runoff step that $g$ falls in, and $U_i$ sums the upstream discharges.
Everything but $c3_i Q_i(g-1)$ is known before the step starts, so the serial chain is just that one term, and the
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
routed. Element 0 holds $U_i(-1)$, the sum of the upstreams' initial states, and the rest hold $U_i(g)$. A row is live
only from when a river's first upstream is routed until the river itself is routed, so a small pool of rows is reused
instead of holding an `(n_rivers, n_routing_steps)` array. Two rows follow the pool: a row of zeros that headwaters
read as their inflow, and a sink that basin outlets write into.

### Blocks and layout

Rivers are handled in blocks of `BLOCK` rivers (64) so that `(time, river)` C-order arrays are read and written one
short contiguous run per time step, rather than one strided element per river per step. A block's scratch is always
`(rivers, time)`, so the routing itself always works river major; the blocking exists to bridge to and from
`(time, river)` arrays, and to keep a block of aggregated runoff in cache while it is routed. Lateral inflow arrives
one of three ways:

| Kernel                                              | Lateral inflow                                        |
|-----------------------------------------------------|-------------------------------------------------------|
| `static_vlateral` / `dynamic_vlateral`, `by_river=False` | `(time, river)` C-order, copied into a scratch block   |
| `static_vlateral` / `dynamic_vlateral`, `by_river=True`  | `(river, time)` C-order, each river's row read in place |
| `static_grid`                                        | gridded runoff aggregated per block, no vlateral array |


## Lateral inflow

| You have                                                     | What happens                                                      |
|--------------------------------------------------------------|-------------------------------------------------------------------|
| Gridded runoff and a weight table                            | aggregation is fused into routing, no vlateral array is built      |
| vlateral files, or your own runoff reader                    | the array is read as it is given                                   |
| A reader that builds each river's series contiguously        | hand the router the `(river, time)` array's `.T`, read with no copy |

A `(time, river)` vlateral array is used as is, one block of rivers copied into scratch at a time. A reader that
builds `(river, time)` arrays should pass their transpose, which keeps the `(time, river)` shape but is read one
river row at a time with no copy at all. Discharge is always written `(river, time)`, and that C-order array is what
a discharge writer receives; `writers.to_time_major` transposes it for a format that needs each time step's rivers
contiguous, as `parquet_writer` does.

Inputs are NaN free by the time they reach a kernel: gridded runoff has NaN cells set to zero before it is aggregated,
so a missing cell contributes nothing while the other cells of its catchment still count, and vlateral read from files
has NaN volumes set to zero.
