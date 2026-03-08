# Parallelism in River Routing

This page explains why river-route uses a serial execution model and documents the parallelism
strategies that were investigated and ultimately rejected. Each approach is described along with
the reasons it did not produce a net performance benefit. The intent is to save future developers
the effort of re-exploring these paths and to illustrate broader lessons about when parallelism
helps and when it does not.

---

## Runtime Profile

After profiling and optimizing the routing pipeline on a representative workload (12 input files,
162,000 river segments), the wall-clock time breaks down as:

| Component              | Time  | Share |
|------------------------|-------|-------|
| Numba routing kernel   | ~8.4s | 45%   |
| I/O reads              | ~4.4s | 24%   |
| Python overhead, dtype | ~2.1s | 11%   |
| I/O writes             | ~0.7s | 4%    |
| Everything else        | ~2.4s | 13%   |

The routing kernel dominates. It solves a sparse triangular system via forward substitution at
every time step. Understanding why that solve resists parallelism is the key to understanding why
the strategies below did not pan out.

---

## The Core Constraint: Forward Substitution Is Serial

The Muskingum routing equation at each time step has the form $L\, Q_{t+1} = b$, where $L$ is a
unit lower triangular matrix derived from the topologically-sorted river network adjacency matrix.
This system is solved by forward substitution (see the [Forward Substitution](forward-substitution.md)
page for details).

Forward substitution is inherently sequential: the value at row $i$ depends on all previously
solved rows $1, \ldots, i-1$. There is no way to compute row $i$ before its upstream dependencies
are known. This means the core solve cannot be split across threads or cores in a straightforward
way.

Every parallelism strategy must work around this constraint, and each one encounters a different
limiting factor.

---

## Strategy 1: Pipelined I/O with Thread Queues

### Concept

The routing loop processes input files one at a time: read a file, route it, write the output,
then move to the next file. Since disk I/O and CPU computation use different hardware resources,
they could in principle overlap. A background reader thread pre-reads the next input file while
the main thread routes the current one. A background writer thread writes the previous output
while the main thread has already moved on to the next file.

This is a classic producer-consumer pattern using Python's `queue.Queue` for coordination. The
Numba JIT-compiled routing kernel releases the GIL, so the I/O threads can genuinely run
concurrently with the computation.

### Why It Did Not Help

This strategy was implemented and benchmarked. On the largest computation groups running as a
single process, it produced a 10-15% wall-time reduction. However, this gain disappeared in
practice for two reasons:

1. **The I/O is already fast relative to compute.** Reads account for 24% of total time and writes
   only 4%. After other optimizations reduced the compute time and improved array handling, the
   I/O share shrank further. There is not enough I/O latency left to hide.

2. **Production runs use multiple processes.** In deployment, each computational group (VPU) runs
   as a separate OS process. When many processes run concurrently, the OS scheduler naturally
   interleaves their I/O and compute phases. Process A reads while process B routes and process C
   writes. This implicit pipelining achieves the same overlap without any threading code. Adding
   explicit per-process I/O threads on top of this provides no additional benefit and increases
   memory pressure, since each process must hold two input arrays in memory simultaneously instead
   of one.

### Lesson

Pipelined I/O helps when I/O is a large fraction of wall time and the process has exclusive access
to the disk. When the code is already compute-dominated or when multiple processes share the
system, the technique adds complexity without measurable benefit.

---

## Strategy 2: File-Level Parallelism for Ensemble Mode

### Concept

In ensemble mode, each input file starts from the same initial channel state and produces
independent results. Since there are no data dependencies between files, they could be routed in
parallel using a process pool (`concurrent.futures.ProcessPoolExecutor`). Each worker receives a
copy of the initial state and routing parameters, routes one file, and returns the result. The
final channel state is the mean of all member final states, same as the serial version.

### Why It Was Not Pursued

While theoretically sound, several practical issues make this unattractive:

1. **Memory cost.** Each worker needs its own copy of all routing arrays. For 162,000 rivers at
   731 time steps, each file's lateral inflow array is roughly 900 MB. With 8 workers, that is
   over 7 GB of additional memory just for the input arrays, plus duplicated coefficient arrays
   and output buffers.

2. **Serialization overhead.** Passing the router state to worker processes requires serializing
   large numpy arrays and sparse matrices across process boundaries. This overhead partially
   offsets the parallelism gain.

3. **Limited applicability.** This only works for ensemble mode. The more common sequential mode,
   where each file's final state becomes the next file's initial state, has an inherent serial
   dependency between files that cannot be broken.

4. **Diminishing returns from the routing kernel.** The routing kernel itself is already fast. The
   per-file compute time is on the order of one second. Process startup, data transfer, and
   synchronization overhead consume a meaningful fraction of that.

### Lesson

Embarrassingly parallel workloads are the easiest to parallelize, but the gain must outweigh the
overhead of distributing work. When each unit of work is small (a few seconds) and the data is
large (hundreds of megabytes), the coordination cost can dominate.

---

## Strategy 3: Independent Sub-Basin Parallelism

### Concept

The adjacency matrix represents a forest of trees, where each tree is a connected river basin
rooted at an outlet segment. Sub-basins that do not flow into each other have no data dependencies,
so their forward substitution solves could run in parallel. The approach would be to decompose the
network into connected components at initialization, reorder the adjacency matrix into
block-diagonal form, and use `numba.prange` to solve independent blocks concurrently.

### Why It Was Not Pursued

1. **Network topology is unfavorable.** Many VPUs are dominated by a single large connected basin.
   If 90% of the segments belong to one component, only 10% of the work can be parallelized. The
   speedup is bounded by the largest component (Amdahl's law).

2. **Implementation complexity.** Decomposing the network requires graph analysis at initialization,
   reordering all arrays to make each block contiguous, maintaining per-block CSC index arrays or
   offset schemes, and mapping the output discharge array back to the original segment ordering.
   This is substantial code complexity for an uncertain payoff.

3. **Thread overhead for small blocks.** River networks often have many tiny tributary components
   alongside one large trunk system. Launching parallel work for blocks of 5-10 segments is slower
   than just solving them serially. A block-size threshold would need to be tuned, adding another
   dimension of complexity.

### Lesson

Parallelism that depends on the structure of the input data is fragile. If the favorable case
(many large independent blocks) rarely occurs in practice, the engineering investment is wasted.
Analyze real-world data distributions before committing to a topology-dependent parallelization
strategy.

---

## Strategy 4: Numba prange on Element-Wise Loops

### Concept

Within each time step, the routing kernel performs three phases:

1. **RHS assembly** — compute `rhs[i] = c3[i] * q_t[i] + lateral[i]` for each element
   independently.
2. **Sparse matrix-vector product** — accumulate upstream contributions into the RHS vector.
3. **Forward solve** — solve the triangular system sequentially.

Phases 1 and 2 involve loops over all elements that are independent per element, making them
candidates for `numba.prange` (Numba's parallel range). This requires adding `parallel=True` to
the `@numba.njit` decorator.

### Why It Made Things Worse

This strategy was implemented and benchmarked. Performance was consistently worse than the serial
version, even on networks with 162,000 segments. The reasons:

1. **Memory bandwidth saturation.** The element-wise loops perform 2-3 floating-point operations
   per memory access (a multiply-add on array elements). This is far below the compute-to-memory
   ratio needed to benefit from multiple cores. A single core can already saturate the memory bus
   for sequential array traversal. Adding more cores does not help because they are all waiting on
   the same memory bus.

2. **Thread synchronization overhead.** `prange` partitions the loop across threads and
   synchronizes them at the end of each parallel region. With multiple parallel regions per time
   step and thousands of time steps, the cumulative synchronization cost is significant. Each
   barrier costs microseconds, but multiplied by hundreds of thousands of invocations, it adds up
   to measurable overhead.

3. **The bottleneck is elsewhere.** The forward solve (phase 3) is the most expensive phase and
   is strictly serial. Even if phases 1 and 2 were made infinitely fast, the total speedup would
   be modest because the serial phase dominates. This is Amdahl's law in action: parallelizing the
   non-bottleneck phases yields diminishing returns.

4. **Sparse matrix-vector products have write conflicts.** The SpMV loop iterates by column (CSC
   format) and accumulates into rows. Different columns can write to the same row, creating race
   conditions that prevent direct parallelization. A row-oriented (CSR) formulation would fix this,
   but the rest of the algorithm (forward solve) requires CSC access patterns. Maintaining both
   formats doubles the sparse matrix storage.

### Lesson

Not all independent loops benefit from parallelism. When the work per element is small and the
data access pattern is sequential, the overhead of thread management exceeds the benefit of
concurrent execution. Profile the memory bandwidth utilization before parallelizing array loops.
If a single core already saturates the bus, adding more cores will not help.

---

## Summary

| Strategy                   | Tested | Outcome                  | Primary Limiting Factor                    |
|----------------------------|--------|--------------------------|--------------------------------------------|
| Pipelined I/O              | Yes    | Small gain, not retained | I/O share too small; implicit OS pipelining|
| File-level ensemble        | No     | Not pursued              | Memory cost; limited to ensemble mode      |
| Sub-basin decomposition    | No     | Not pursued              | Unfavorable network topology               |
| prange on element loops    | Yes    | Negative gain            | Memory bandwidth bound; sync overhead      |

The serial implementation is already efficient. The routing kernel runs at $O(n + m)$ per time step
with JIT-compiled native code, the sparse matrix storage is compact, and the I/O layer uses
efficient binary formats. The total time for a 162,000-segment, 12-file computation group is under
20 seconds. The places where parallelism could theoretically help are either too small a fraction
of the total (I/O), limited by memory bandwidth (element loops), or structurally unfavorable
(network topology, file dependencies).

The broader takeaway is that parallelism is not free. Thread synchronization, memory duplication,
process serialization, and memory bus contention all have costs. When the serial code is already
close to the hardware limits for its workload pattern, those costs exceed the benefits. Profiling
and understanding the hardware bottleneck must come before reaching for concurrency.
