# Parallelism in River Routing

Many scientific computations jump quickly to various methods of parallelism and GPU acceleration. River routing
math is very well suited for parallelization, but it is even more advantageous to worry first about the efficiency 
of the algorithm. There are more efficiency gains to be had through careful choice of solving procedures and 
preparing the inputs in more careful ways than there is through jumping immediately to parallelization. The 
effectiveness of parallelization depends on the strategy, the size and complexity of the network, and the 
hardware being used. 

Not all parallelization strategies are worth pursing in river routing cases. This page is a list of the strategies 
tested in `river-route` and recommendations based on using these methods to operate a global hydrological model 
generating a 5 trillion data point simulation product.

The following guide is the order of significance in how effective various methods are in speeding up and increasing 
resource efficiency of routing computation. This is shared both as justification for why `river-route` is designed 
in this way and also as an educational guide for other projects. 

## What cannot be parallelized?

The two fundamental constraints in river routing. First, it is a time stepping process. At any given river, you must 
solve for discharge at the current time `t` before solving for the next time `t+1`. Second, there is an 
order dependency that river segment upstream must be solved before the current river and the downstream river.

Within a time step, the solve is a forward substitution: the value at row $i$ depends on all previously solved
rows $1, \ldots, i-1$. Water only moves downstream, so a river cannot be computed before every river upstream of it.
A large river's main stem depends on its whole basin and is always computed on one thread.

## Better ways to prepare inputs

Much of the slowness in a routing scheme can be entirely avoided, not just sped up, by choosing the best ways to
prepare computations.

### Topological river sorting and Depth First Search (DFS)

TBD

### Splitting watershed subgraphs

[//]: # (todo: talk about identificaiton of subgraphs and balancing size with resources)

**Summary**: If your computations contains several independent watersheds, you can route them simultaneously in separate processes.

**Conclusion**: This is more beneficial as job sizes get larger. Watersheds have no dependencies on others. Separate watersheds and
process simultaneously or combine them into a single config file if compute times are small enough.

### Formats of inputs and outputs

[//]: # (todo)

after you carefully prepare inputs and the algorithm, a large, possibly the largest, portion of remaining time is spent reading and writing data.
You should reduce the number of times the code needs to read/write data and the number of total files it needs to read/write.
Faster storage formats and fewer, better sized files do more than parallelizing file I/O.

## The forward substitution solving algorithm

Forward substitution is inherently sequential: the value at row $i$ depends on all previously
solved rows $1, \ldots, i-1$. There is no way to compute row $i$ before its upstream dependencies
are known. This means the core solve cannot be split across threads or cores in a straightforward
way.

However, the routing kernels have been made about as minimal as possible: a single topological sweep per routing
step, sparse connectivity instead of matrices, and JIT compiled numba code. In my experience, this is preferable to
multiprocessing methods even though it uses an inherently sequential forward substitution algorithm. This approach is
the best method in my experience using it on a wide range of scales up to global computations of hourly resolution
discharge on millions of rivers and producing a 5 trillion data point simulation. It has the advantages that it:

1. needs only mainstream scientific python dependencies simply installed on a variety of hardware and Python versions
2. is the most memory efficient option
3. is the computationally fastest option because it does not iterate or do any matrix conditioning or pivoting
4. is the direct solution so there is no error due to solver convergence tolerances.

What you control is how much work you ask that kernel to do.

1. **Choose the simplest routing procedure your problem needs.** `coeff: static` computes Muskingum coefficients
   once and reuses them for every file with the same time steps. `coeff: dynamic` rebuilds them inside the kernel
   on every substep. Only pay for dynamic coefficients when the application needs them. See the
   [config file reference](config-files.md#routing-procedure-selectors).
2. **Use the largest stable routing time step.** Every routing substep is a full sweep of the network, so
   `dt_routing` directly sets the amount of work. Check stability with `Network.stability_report(dt)` rather than
   defaulting to a small step. See [Time Variables](time-options.md).
3. **Only produce the output you will use.** A coarser `dt_discharge` averages results before they are written, and
   a [custom writer](../tutorial/advanced.md#customizing-outputs) can save only the rivers you need. The premade
   `zarr_writer` is built to write as fast as possible when you need everything.
4. **Route what matters.** Rivers that do not contribute to the results you need do not have to be in the network.

## Options for parallelism

Once a job is efficient on a single core, these strategies can add to it. Measure each one against your
single-threaded result on your own hardware before keeping it.

### Multithreading matrix solvers

**Summary**: A network in DFS computation order can be split into independent upstream regions that are routed
concurrently, followed by the main stem on a single thread.

`river-route` never creates threads on its own. Threads are a runtime resource, not a config, so pass a
`ThreadPoolExecutor` and `threads`, the number of regions to split the network into, to `Router.route`. The same pool
is used to aggregate gridded runoff when routing from `grid_runoff_files`.

```python title="Threaded Routing"
from concurrent.futures import ThreadPoolExecutor

import river_route as rr

router = rr.Router(rr.Configs.from_file('config.yaml'))
with ThreadPoolExecutor(max_workers=8) as pool:
    router.route(thread_pool=pool, threads=8)
```

The speedup is limited by the rivers left in the sequential main stem and by memory bandwidth, which the routing
kernel is largely bound by. Use `river_route.network.streams.analyze_partitioning` to see how a network splits and the upper
bound on speedup before committing to it. The partition depends only on connectivity and `threads`, never on the
forcing, dt, or coefficients, so the `Network` derives it once and caches it per thread count: every simulation
over one `Network` reuses it. `river_route.network.streams.partition_network` can also store it as a `region` column in
the parameter file, which a `Network` reuses as-is instead of deriving one at all.

**Conclusion**: Meaningful speedup is possible with multiple threads but only if you sort the network into
independent but ordered subgraphs. This is an optional addition to a job that is already efficient single threaded.
It is not a substitute for a better algorithm or better prepared inputs.

### Concurrent jobs vs multiple threads in one job

**Summary**: Simulations of many inputs, perhaps from an ensemble of runoff projections, share
only the initial state. Multiple members can be processed concurrently in separate processes.

In production, you will usually add logic to:

1. Set the config variables and routing parameter files
2. Find all the catchment runoff files, 1 for each member
3. Determine unique output names for each so you don't overwrite or corrupt files
4. Submit each input/output file pair to a parallel processing framework

```python title="Parallel Routing Jobs for Ensemble Members"
from multiprocessing import Pool

import river_route as rr

params_file = 'routing_parameters.parquet'
runoff_files = ['catchment_runoff_member_1.nc',
                'catchment_runoff_member_2.nc', ]
output_files = ['discharges_member_1.nc',
                'discharges_member_2.nc', ]


def route(input_file: str, output_file: str) -> None:
    configs = rr.Configs(
        forcing='vlateral',
        params_file=params_file,
        vlateral_files=[input_file, ],
        discharge_files=[output_file, ],
    )
    rr.Router(configs).route()


if __name__ == '__main__':
    with Pool() as pool:
        pool.starmap(route, zip(runoff_files, output_files))
```

**Conclusion**: This is the best way to speed up ensemble simulations. 

### Separate pipelines for reading inputs, compute, writing outputs

**Summary**: A single routing process can have up to 3 meaningful tasks: 1 reads inputs from disk, 1
does computations, 1 writes results to disk. All 3 can be operating independently at the same time.

```mermaid
block-beta
    columns 7
    space:1 s1["Step 1"] s2["Step 2"] s3["Step 3"] s4["Step 4"] s5["Step 5"] s6["Step 6"]
    r["Read"]:1 r1["t=1"] r2["t=2"] r3["t=3"] r4["t=4"] space:2
    c["Compute"]:1 space:1 c1["t=1"] c2["t=2"] c3["t=3"] c4["t=4"] space:1
    w["Write"]:1 space:2 w1["t=1"] w2["t=2"] w3["t=3"] w4["t=4"]
```

If you have an unfavorable combination of slow I/O, slow CPU, and large computations, this solution
might help. Individual routing jobs get faster by making threads for portions that depend on
different hardware. However, using this method means you probably won't be able to use it in
combination with another parallelization strategy because you more quickly consume memory and disk
I/O bandwidth with one job. In my experience, this speedup is at most a few percent.

**Conclusion**: This speeds up individual jobs bottlenecked by I/O but not by much given modern hardware
capabilities. Faster storage formats and fewer, better sized files do more.

## Conclusions

The most important thing to work on is the efficiency of the algorithm and statically determining subgraphs.
