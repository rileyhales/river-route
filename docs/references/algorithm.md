# Routing Algorithm

[Math Derivations](math.md) derives the equation every river is routed with. This page shows the arrays the routing
kernels read and write, how each is laid out in memory, and how a kernel finds a river's values in them. The names are
those of `river_route/kernels/c/routing.h`, which the numba kernel, `numba_kernel.py`, shares.

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
depth rows hold every step of the forcing file, so its `cell_stride` is the file's count of steps even when `dt_total`
routes fewer.
