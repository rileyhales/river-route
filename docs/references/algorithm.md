# Routing Algorithm

[Math Derivations](math.md) derives the equation every river is routed with. This page shows how much work the x of
each river asks of the routing kernels, the arrays they read and write, how each is laid out in memory, and how a
kernel finds a river's values in them. The names are those of `river_route/kernels/c/routing.h`, which the numba
kernel, `numba_kernel.py`, shares.

---

## The order of routing

One call of the kernel, `rr_route_rivers`, routes one job (see [Parallelism](parallelism.md)): the job's blocks of
rivers in routing order, each river over every step before the next river, and each of a river's sub-reaches over every
step before the next sub-reach.

```python
for block in range(n_blocks):                                     # in routing order
    for river in range(block_starts[block], block_stops[block]):  # in routing order
        for sub_reach in range(sub_reach_starts[river], sub_reach_starts[river + 1]):  # upstream to downstream
            for step in range(n_steps):                           # eight steps at a time, then any left one at a time
                discharge[river, step + 1] = ...                  # from the inflow and runoff of the step
```

A river's first sub-reach takes its inflow from the discharge rows of the river's inflow rivers. Every sub-reach after
it takes its inflow from the river's own row, where the sub-reach before it left its discharge, and writes its own
discharge over that row, so the row ends with the discharge of the river's last sub-reach.

## How x sets the work of a river

The kernel routes each sub-reach of a river over every step, so a river split into $N$ sub-reaches takes $N$ times the
work of a river routed whole. Subcycles take no more work, since a river's subcycles are combined into one step of
coefficients before routing. The work of a network in each step is its count of sub-reaches, and with
`use_subreaches` that count is set by the x of each river.

### The sub-reaches a river needs

Every coefficient of a river is non-negative when its step is inside the window of
[Signs of the coefficients](math.md#signs-of-the-coefficients):

$$
2kx \le \Delta t \le 2k(1-x)
$$

A river too long for its step, $\Delta t < 2kx$, is split into $N$ equal sub-reaches of travel time $k/N$, each with
$x_N = 1/2 - N(1/2 - x)$, which keep the river's mean travel time $k$ and its travel time variance $k^2(1-2x)$. With
$2 (k/N) x_N = k/N - (1-2x)k$ and $2 (k/N)(1 - x_N) = k/N + (1-2x)k$, the window of each sub-reach is

$$
\frac{k}{N} - (1-2x)k \le \Delta t \le \frac{k}{N} + (1-2x)k
$$

It is centered on the travel time of the sub-reach, $k/N$, and is $(1-2x)k$ wide on either side, the diffusion of the
river. The fewest sub-reaches whose lower bound is under the step are

$$
N = \left\lceil \frac{k}{\Delta t + (1-2x)k} \right\rceil
$$

and since $\Delta t > 0$, the fraction is always less than $k / ((1-2x)k)$, so

$$
N \le \left\lceil \frac{1}{1-2x} \right\rceil
$$

for every river, whatever its travel time and whatever the step. That bound is what x costs: 3 at $x = 0.3$, 10 at
$x = 0.45$, 50 at $x = 0.49$, and 500 at $x = 0.499$. At $x = 1/2$ there is no bound. The window closes to the single
step $\Delta t = k/N$, so the river needs $\lceil k/\Delta t \rceil$ sub-reaches, each exactly one step long, and their
count grows with the river's travel time alone. Such a river can be made valid only when its $k$ is exactly a whole
multiple of $\Delta t$, so most rivers at $x = 1/2$ are kept whole and reported.

### Sub-reaches that take subcycles of their own

Equal sub-reaches can be too short for the step, with a negative $c_3$ in the shortest, when
$\Delta t > k/N + (1-2x)k$. That happens only for $x > 3/8$, and a river too short for the step can likewise fall
between two whole counts of subcycles only for $x > 1/3$. With `use_subreaches` and `use_subcycles` both, such a river
is split into sub-reaches of unequal travel times $k_i$, each taking $m_i$ subcycles of its own. Every one takes
$x_i = 1/2 - (1/2 - x)\,k/k_i$, the river's own diffusion, which keeps its travel time variance whatever the $k_i$,
and its window is

$$
\left| \frac{\Delta t}{m_i} - k_i \right| \le (1-2x)k
$$

as wide for every sub-reach as for the river's equal ones. The travel times a sub-reach may take are windows about
$\Delta t$, $\Delta t/2$, $\Delta t/3$, and so on, none shorter than $(1-2x)k$, where $x_i = 0$, and the river splits
into $n$ sub-reaches when $k$ is a sum of $n$ of them. `river_route.network.splitting` finds the fewest by growing the
sums of 1, 2, 3, ... windows as intervals, and each sub-reach takes the fewest subcycles whose step is in its window.
They are never fewer than the $N$ above, and in 200,000 rivers of random travel times, from a hundredth of the step to
fifty steps, and x from 1/3 to 0.4999, never more than $\lceil 1/(1-2x) \rceil$. A river of
$\Delta t = 0.6k$ at $x = 0.49$, which no count of equal sub-reaches makes valid, splits into three, of 0.601, 0.301,
and 0.098 of its $k$ in 1, 2, and 6 subcycles.

No sub-reach is shorter than 1/64 of the river or of the step, whichever is shorter. Near $x = 1/2$ the windows
narrow, and a split can need a sub-reach far shorter, taking thousands of subcycles only to make the travel times sum
to $k$; a river whose every split needs one is kept whole and reported, as is a river with no valid split into the 256
sub-reaches the kernels route one river as. Subcycles add no routing work, so the work of a river is still its count of
sub-reaches.

### An example

The sub-reaches of rivers routed at an hourly step, $\Delta t$ = 1 h, and the spread of each river's response, the
standard deviation of its travel time over its mean, $\sigma / k = \sqrt{1-2x}$:

| $x$   | Most any river needs | $k$ = 1 day | $k$ = 10 days | $\sigma / k$ |
|-------|----------------------|-------------|---------------|--------------|
| 0.3   | 3                    | 3           | 3             | 0.63         |
| 0.4   | 5                    | 5           | 5             | 0.45         |
| 0.45  | 10                   | 8           | 10            | 0.32         |
| 0.48  | 25                   | 13          | 23            | 0.20         |
| 0.49  | 50                   | 17          | 42            | 0.14         |
| 0.495 | 100                  | 20          | 71            | 0.10         |
| 0.499 | 500                  | 23          | 163           | 0.04         |
| 0.5   | no bound             | 24          | 240           | 0            |

The river of 10 days at $x = 0.45$ needs

$$
N = \left\lceil \frac{240}{1 + 0.1 \cdot 240} \right\rceil = \left\lceil \frac{240}{25} \right\rceil = \lceil 9.6 \rceil = 10
$$

and at $x = 0.499$ it needs

$$
N = \left\lceil \frac{240}{1 + 0.002 \cdot 240} \right\rceil = \left\lceil \frac{240}{1.48} \right\rceil = \lceil 162.2 \rceil = 163
$$

so it is routed 163 times over every step instead of 10. A river of 30 days at $x = 0.499$ needs
$\lceil 720 / 2.44 \rceil = 296$, more than the 256 sub-reaches the kernels route one river as, and
`Network.check_routable` refuses the network before any runoff is read.

### What holding x at 0.45 changes

A river's response to its upstream inflow has mean $k$ and variance $k^2(1-2x)$ at any step, so x sets how much it
attenuates a flood wave: from a linear reservoir at $x = 0$, whose response spreads with a standard deviation of $k$,
to translation without attenuation at $x = 1/2$. In Muskingum-Cunge terms, a river of length $\Delta x$ and wave
celerity $c$ represents a diffusion of $(1/2 - x)\, c\, \Delta x$. The rivers that need the most sub-reaches are those
near $x = 1/2$, which attenuate least, so holding x at 0.45 changes their routing little:

- It adds a variance of at most $0.1 k^2$ to any river, a standard deviation of at most $0.32 k$, and a diffusion of at
  most $0.05\, c\, \Delta x$.
- Along rivers in series, travel times add and so do variances, so along a path of $n$ rivers of equal $k$ the
  standard deviation it adds is at most $0.32 / \sqrt{n}$ of the path's travel time: 6% along 25 rivers.
- It holds every river to 10 sub-reaches or fewer, so a network takes at most 10 times the work of routing every river
  whole, and no river reaches the kernels' limit of 256.

### Holding x lower

`river-route` routes every river with the x of its network file, exactly as it is given. A modeler who would rather
bound the work may choose to hold x at 0.45 in the network file before routing:

```python
import numpy as np
import pandas as pd

network = pd.read_parquet('network.parquet')
network['muskingumX'] = np.minimum(network['muskingumX'], np.float32(0.45))
network.to_parquet('network_x_at_most_045.parquet')
```

Routing with `dt_routing` in seconds and `use_subreaches: false` instead routes every river too long for its step
whole, with no added work, and with `invalid_coefficients: 'log'`, since the negative $c_1$ of each is an error
that the default, `'raise'`, refuses to route.

## The arrays of a call

A call is given three structs, the network, the runoff, and the job, holding these arrays and the counts of rivers,
sub-reaches, cells, steps, and blocks:

| Struct       | Array                         | Length          | What it holds                                                                           |
|--------------|-------------------------------|-----------------|-----------------------------------------------------------------------------------------|
| `rr_network` | `sub_reach_starts`            | rivers + 1      | each river's first sub-reach, then the count of sub-reaches                             |
|              | `c1`, `c2`, `c3`, `c4`        | sub-reaches     | each sub-reach's coefficients                                                           |
|              | `inflow_river_starts`         | rivers + 1      | each river's first entry in `inflow_rivers`, then the count of entries                  |
|              | `inflow_rivers`               | entries         | the rivers draining into each river                                                     |
|              | `weight_starts`               | rivers + 1      | each river's first weight, then the count of weights                                    |
|              | `weight_cells`, `weights`     | weights         | the cell each weight reads, and the runoff volume (m³) a meter of runoff depth gives it |
| `rr_runoff`  | `depth`                       | a row per cell  | each cell's runoff depth (m) in each step                                               |
| `rr_job`     | `block_starts`, `block_stops` | blocks          | the first river of each block, and one past its last                                    |
|              | `channel_state`               | sub-reaches     | each sub-reach's discharge, at the start and then at the end                            |
|              | `discharge`                   | a row per river | each river's discharge at the start and at the end of each step                         |

Only `depth` and `discharge` grow with the number of steps. They are tables, a row for each cell or river, and the
next sections show how their rows are found. Every other array is a list for every river.

## A list for every river

Each river has a list of sub-reaches, a list of inflow rivers, and a list of weights, each of its own length. The
kernels keep the lists of every river one after another in one array, and a starts array of `n_rivers + 1` entries
gives where each river's list begins: river `r`'s list is the entries from `starts[r]` up to `starts[r + 1]`. The last
entry is the length of the array, where the last river's list ends. A river with an empty list, such as a headwater's
list of inflow rivers, begins where the next river's list does.

<div style="text-align: center;">
<div class="algorithm-figure">
--8<-- "docs/static/images/river-lists.svg"
</div>
<figcaption><em>Figure 1: The lists of three rivers. Rivers 0 and 1 drain into river 2, which is routed as two sub-reaches. Rivers 0 and 1 each read one runoff cell and river 2 reads both. A river's color marks its entries, and each arrow points to where a river's list begins. Rivers 0 and 1 have no inflow rivers, so their lists in <code>inflow_rivers</code> are empty and begin where river 2's does.</em></figcaption>
</div>

The coefficients and `channel_state` are indexed by sub-reach, so `sub_reach_starts` finds a river's entries in all
five. `weight_cells` and `weights` are indexed by weight, so `weight_starts` finds a river's entries in both.

## Rows and strides

Memory is one line of values, so a table is stored one row after another, and a kernel is given a pointer to its first
value. The stride of a table is how many values apart its rows are: the count of values from a row's value at one step
to the next row's value at the same step, which is the length of a row in memory. A kernel finds every value of a table
with it, cell `c` or river `r` at step `t`:

```c
depth[c * cell_stride + t]
discharge[r * river_stride + t]
```

<div style="text-align: center;">
<div class="algorithm-figure">
--8<-- "docs/static/images/depth-rows.svg"
</div>
<figcaption><em>Figure 2: The runoff depth of two cells over six steps, as a table and as it is stored. The numbers under the stored values are their positions in memory. Cell 1's depth in step 2 is <code>cell_stride</code>, 6 values, after cell 0's.</em></figcaption>
</div>

A discharge row holds one more value than the steps routed. Its first value, `t = 0`, is the river's discharge at the
start, which the kernel copies from the `channel_state` of the river's last sub-reach before routing, and `t = 1` to
`n_steps` are its discharge at the end of each step.

<div style="text-align: center;">
<div class="algorithm-figure">
--8<-- "docs/static/images/discharge-rows.svg"
</div>
<figcaption><em>Figure 3: The discharge of three rivers over four steps, as a table and as it is stored. The darker first value of each row, t = 0, is the river's discharge at the start, and t = 1 to 4 are its discharge at the end of each step. Each row is 5 values long, so <code>river_stride</code> is 5.</em></figcaption>
</div>

With both strides and the lists of Figure 1, step `t` of river `r` reads its inflow at the start and at the end of the
step from columns `t` and `t + 1` of its inflow rivers' rows, and its runoff volume from column `t` of its weights'
cells' rows:

```python
entries = range(inflow_river_starts[r], inflow_river_starts[r + 1])
inflow_before = sum(discharge[inflow_rivers[e] * river_stride + t] for e in entries)  # I at the start of step t
inflow_after = sum(discharge[inflow_rivers[e] * river_stride + t + 1] for e in entries)  # I at its end
runoff_volume = sum(
    weights[w] * depth[weight_cells[w] * cell_stride + t] for w in range(weight_starts[r], weight_starts[r + 1])
)  # V in step t
```

The two tables have strides of their own: their rows hold different counts of values, and either can be part of a
longer table, as in routing in windows. The kernels check that every row has room for the steps routed, `cell_stride`
at least `n_steps` and `river_stride` more than `n_steps`, and read only the first `n_steps` values of each depth row
and write only the first `n_steps + 1` values of each discharge row.

## Routing in windows

A caller can route a run a window of steps at a time, calling a kernel once for each window. `channel_state` carries
every sub-reach's discharge from the end of one window to the start of the next, so the windows route the same
discharge as routing every step at once. Each call can be given a pointer to the window's first step in the runoff
depth of the whole run. The rows are then the run's length, so `cell_stride` stays the run's count of steps while
`n_steps` is the window's. The discharge can be a table of its own for the window, a row of the window's steps plus
one for every river, which every window writes over. A last window shorter than the others writes only the first
`n_steps + 1` values of each row.

<div style="text-align: center;">
<div class="algorithm-figure">
--8<-- "docs/static/images/routing-windows.svg"
</div>
<figcaption><em>Figure 4: The second of two windows of 4 steps routing a run of 6 steps. The call routes steps 4 and 5, so <code>n_steps</code> is 2. <code>depth</code> points to cell 0's depth in step 4, and cell 1's is <code>cell_stride</code>, 6 values, after it. Each discharge row has room for a window of 4 steps, so <code>river_stride</code> is 5, and the call writes t = 0 to 2 of each row, where t = 0 is the discharge at the end of step 3, from <code>channel_state</code>.</em></figcaption>
</div>

## Where the strides are set

A kernel never computes a stride. The caller that made the tables knows how long their rows are and passes it.

| Caller                                                    | `cell_stride`                                   | `river_stride`                               |
|-----------------------------------------------------------|-------------------------------------------------|----------------------------------------------|
| `kernels.route`, for the C kernel                         | `depth.strides[0] // depth.itemsize`            | `discharge.strides[0] // discharge.itemsize` |
| `flat_view` in `numba_kernel.py`, for the numba kernel    | `depth.strides[0] // depth.itemsize`            | `discharge.strides[0] // discharge.itemsize` |
| JavaScript laying out windows in a WebAssembly kernel     | the window's count of steps                     | the window's count of steps + 1              |

numpy keeps an array's strides in bytes, so dividing by the size of a value, 4 bytes for float32, gives them in
values. A slice keeps the strides of the array it is cut from: `depth[:, 64:128]` is 64 steps of every row, and its
`cell_stride` is still the length of a whole row. The Router routes each forcing file into the rows of
`memory.reshape(network.size, num_runoff_steps + 1)`, so its `river_stride` is the steps routed plus one, and its
depth rows hold every step of the forcing file, so its `cell_stride` is the file's count of steps.
