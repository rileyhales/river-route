# Routing Kernels

`routing_order` sets how a kernel walks the network. Both orders compute the same discharge to within float32
round-off, and both route concurrently on a `thread_pool` by splitting the network into regions.

- **`'river'`** solves one river's whole time series, then moves to the next river. The time series of one river
  is a simple recurrence, which is solved eight steps per serial operation instead of one. It can also convert gridded
  runoff to lateral inflow in the same pass, so no vlateral array is built. It needs every time step of a file's
  forcing at once.
- **`'time'`** solves every river for one time step, then moves to the next step. Each river waits on the one
  before it within a step, which makes it several times slower, but it advances the whole network one step at a
  time, which is the shape needed to exchange state with another model every step.

## Ranking

Region 6020006540 (303,097 rivers), 3 months of hourly ERA5 runoff per file, mean of 2 runs over 4 files (one year),
on an Apple M3 Max. Routing only, in seconds. Rows are fastest first within each problem.

| Problem         | Order | Lateral inflow it reads                   | Discharge written as      | 1 thread | 12 threads |
|-----------------|-------|-------------------------------------------|---------------------------|---------:|-----------:|
| static, lateral | river | vlateral `(river, time)`                  | `(river, time)` ¹         | **1.74** |   **0.21** |
| static, lateral | river | gridded runoff and a weight table (fused) ² | `(river, time)` ¹       |     2.09 |       0.23 |
| static, lateral | river | gridded runoff and a weight table (fused) ² | `(time, river)`         |     2.72 |       0.45 |
| static, lateral | river | vlateral `(time, river)`                  | `(river, time)` ¹         |     2.99 |       0.63 |
| static, lateral | time  | vlateral `(time, river)`                  | `(time, river)`           |    13.89 |       1.46 |
| static, channel | river | none                                      | `(river, time)` ¹         | **1.58** |   **0.18** |
| static, channel | time  | none                                      | `(time, river)`           |    13.13 |       1.25 |

¹ The Router writes river order discharge as `(river, time)` when the discharge writer declares
`discharge_layout = 'river'`, as `zarr_writer` and `null_writer` do, and as `(time, river)` for `netcdf_writer`,
`parquet_writer`, and custom writers, which read that layout faster.

² Includes turning runoff into lateral inflow. The other lateral rows start from vlateral, which takes 1.9 s on one
thread or 0.3 s on 12 threads to aggregate from the same runoff. From gridded runoff to discharge, river order takes
2.1 s on one thread and 0.23 s on 12; time order takes 15.8 s and 1.8 s.

## Which to use

| You have                                                     | Use                                                               |
|--------------------------------------------------------------|-------------------------------------------------------------------|
| Gridded runoff and a weight table                            | `routing_order: river`; aggregation is fused into routing          |
| vlateral files, or your own runoff reader                    | `routing_order: river`                                            |
| A reader that builds each river's series contiguously        | `routing_order: river`, and hand the router the `(river, time)` array's `.T` |
| A coupling that must exchange state with another model every step | `routing_order: time`                                         |

River order reads whichever vlateral layout it is given. A `(time, river)` array is used as is. A reader that builds
`(river, time)` arrays should pass their transpose, which keeps the `(time, river)` shape but is read one river at a
time with no copy. Time order always reads `(time, river)` and copies anything else into that layout first.

Inputs are NaN free by the time they reach a kernel: gridded runoff has NaN cells set to zero before it is aggregated,
so a missing cell contributes nothing while the other cells of its catchment still count, and vlateral read from files
has NaN volumes set to zero.
