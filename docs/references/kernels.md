# Routing Kernels

`routing_order` sets how a kernel walks the network. Both orders compute the same discharge to within float32
round-off.

- **`'time'`** solves every river for one time step, then moves to the next step. The network can be split into
  regions and routed on a thread pool, so it is the fastest order when you have cores to spare.
- **`'river'`** solves one river's whole time series, then moves to the next river. Each step waits on a single
  multiply-add instead of on the river before it, so it is 2 to 2.5 times faster on one core. It can also convert
  gridded runoff to lateral inflow inside the same pass, so no vlateral array is built. It needs every time step of a
  file's forcing at once, and it runs on one thread.

## Ranking

Region 6020006540 (303,097 rivers), 3 months of hourly ERA5 runoff per file, mean of 2 runs over 4 files (one year),
on an Apple M3 Max. Times are for routing only unless noted, in seconds.

Rows are fastest first within each problem.

| Problem             | Order | Lateral inflow it reads                     | Threads | 1 thread | 12 threads |
|---------------------|-------|---------------------------------------------|---------|---------:|-----------:|
| static, lateral     | time  | vlateral `(time, river)`                    | any     |     15.6 |    **1.4** |
| static, lateral     | river | vlateral `(river, time)`                    | 1       |  **6.4** |            |
| static, lateral     | river | gridded runoff and a weight table (fused) ¹ | 1       |      7.5 |            |
| static, lateral     | river | vlateral `(time, river)`                    | 1       |      8.7 |            |
| dynamic, lateral ²  | time  | vlateral `(time, river)`                    | any     |     22.1 |    **2.3** |
| dynamic, lateral ²  | river | vlateral `(time, river)` or `(river, time)` | 1       | **10.8** |            |
| static, channel     | time  | none                                        | any     |     13.2 |    **1.3** |
| static, channel     | river | none                                        | 1       |  **6.2** |            |

¹ Includes turning runoff into lateral inflow. The other rows start from vlateral, which takes 3.3 s on one thread or
0.5 s on 12 threads to aggregate from the same runoff. So going from gridded runoff to discharge takes 18.9 s in time
order on one thread, 7.5 s in river order, and 1.9 s in time order on 12 threads.

² Timed with alpha = k and beta = 0, because the region's parameters have no alpha or beta. The kernel still does all
of its nonlinear work.

## Which to use

| You have                                              | Use                                                                 |
|-------------------------------------------------------|---------------------------------------------------------------------|
| Several free cores                                    | `routing_order: time` with a `thread_pool`                           |
| One core, gridded runoff and a weight table           | `routing_order: river` (aggregation is fused into routing)          |
| One core, vlateral files or your own runoff reader    | `routing_order: river`                                               |
| A reader that builds each river's series contiguously | `routing_order: river`, and hand the router the `(river, time)` array's `.T` |
| Many independent simulations at once                  | `routing_order: river`, one simulation per core                      |

River order reads whichever layout it is given. A `(time, river)` vlateral array is used as is. A reader that builds
`(river, time)` arrays should pass their transpose, which keeps the `(time, river)` shape but is read one river at a
time with no copy. Time order always reads `(time, river)` and copies anything else into that layout first.
