# Mathematical Derivations

This page starts from the Muskingum-Cunge equation of one river and rewrites it, one form at a time, into the form the
routing kernels solve. Every form routes the same discharge. Each one needs fewer instructions per step, makes fewer of
those instructions wait on the step before, or stores less.

`river-route` computes in float32. Discharge and runoff are not known to more digits than a float32 holds, and a
float32 takes half the memory of a float64 and fits twice as many values in a vector register.

---

## Summary

1. Rivers are routed in depth first search order, so every river is routed after the rivers draining into it and reads
   their discharge rows.
2. Routing each river's whole time series before the next river leaves every inflow known before a river's first step.
   Each river is then one recurrence, $Q_{t+1} = c_3 Q_t + b_t$, whose terms $b_t$ are all known in advance.
3. Computing $b_t$ first and adding $c_3 Q_t$ last, in one fused multiply-add, leaves one instruction per step that
   waits on the step before.
4. A vector register holds four float32 lanes, and one instruction computes all four. The kernels put four consecutive
   steps of one river in the lanes, so everything that does not wait on the step before is computed four steps at a
   time.
5. Applied eight times, the recurrence is $Q_{t+j} = c_3^{\,j} Q_t + S_j$, where the partial sums $S_j$ do not depend
   on $Q$. One multiply-add per eight steps waits on the steps before, against eight one step at a time.
6. Constants fold into the coefficients once per run, and a river's subcycles fold into one step, so neither costs
   anything per step.
7. Only the discharge and the runoff depth grow with the number of steps. Inflow, runoff volume, terms, and sub-reach
   discharge are computed where they are used and never stored.

## The Muskingum-Cunge equation

Each river, or each sub-reach of a river, is routed by

$$
Q_{t+1} = c_1\, I_{t+1} + c_2\, I_t + c_3\, Q_t + c_4\, V_t \tag{1}
$$

- $Q_t$ is the river's discharge at time $t$, in m³/s.
- $I_t$ is its inflow at time $t$, the summed discharge of the rivers draining into it, in m³/s.
- $V_t$ is its runoff volume over the step from $t$ to $t + 1$, in m³.
- $\Delta t$ is the step, $k$ the travel time of a flood wave through the river, both in seconds, and $x$ a weight
  from 0 to 0.5.

The first three terms are Muskingum routing: continuity over the step,
$\tfrac{1}{2}(I_t + I_{t+1}) - \tfrac{1}{2}(Q_t + Q_{t+1}) = (S_{t+1} - S_t) / \Delta t$, with the storage
$S_t = k\,[x\, I_t + (1 - x)\, Q_t]$, solved for $Q_{t+1}$. The fourth adds runoff as RAPID does (David et al., 2011):
the runoff volume enters as a constant flow $V_t / \Delta t$ added to the inflow at both ends of the step, so $c_1$ and
$c_2$ each multiply it.

Written with $u = \Delta t / k$ and $v = 1 - 2x$, the coefficients are

$$
c_1 = \frac{u - 1 + v}{u + 1 + v}, \qquad
c_2 = \frac{u + 1 - v}{u + 1 + v}, \qquad
c_3 = \frac{1 + v - u}{u + 1 + v}, \qquad
c_4 = \frac{c_1 + c_2}{\Delta t} \tag{2}
$$

These are the usual $(\Delta t/k - 2x) / (\Delta t/k + 2(1 - x))$ and so on, with $u$ and $v$ substituted. In this form
the three share one denominator and their numerators differ only in signs, so computing them takes one division and a
few additions, and the sign of each coefficient can be read from its numerator (see
[Signs of the coefficients](#signs-of-the-coefficients)).

The coefficients sum to one, $c_1 + c_2 + c_3 = 1$, so mass is conserved and a river with a constant inflow and no
runoff settles at $Q = I$.

Written as it stands, each step of each river sums the rows of its inflow rivers twice, once at $t$ and once at
$t + 1$, sums its runoff over its weights, and computes four products and three sums. The rest of this page rewrites
that work.

## Every step of a river before the next river

Rivers drain only downstream and never in loops, so they can be ordered with every river after every river upstream of
it. `river-route` requires the order of a depth first search from each outlet, which also makes every river's
watershed one contiguous range of rows ending at that river (see the
[network file](io-file-schema.md#network-file)). That is what divides the rivers into blocks for threads (see
[Parallelism](parallelism.md)) and what keeps few rows needed at once (see [Memory](#memory)).

The kernels keep each river's inflow rivers in one list, `inflow_rivers`, with each river's first entry in it in
`inflow_river_starts`, and route the rivers of a job in routing order, each river's whole time series before the next
river:

```python
for river in range(n_rivers):   # routing order
    for t in range(n_steps):    # every step of the river before the next river
        discharge[river, t + 1] = ...   # Equation 1, its inflow summed from the rows of its inflow rivers
```

When a river is reached, the rows of its inflow rivers are complete. Its inflow is known at every step before its first
step is routed, and so is its runoff. Everything in Equation 1 except $c_3 Q_t$ is known in advance:

$$
Q_{t+1} = c_3\, Q_t + b_t, \qquad b_t = c_1\, I_{t+1} + c_2\, I_t + c_4\, V_t \tag{3}
$$

Each river is now one first order recurrence with known terms $b_t$, routed one river at a time. The kernels call
$b_t$ the river's *terms*.

A river reads each row from start to end, and the rows it reads were written by the rivers routed shortly before it,
so they are still in the processor's cache.

## Adding c3 Q last

A fused multiply-add computes $a \cdot b + c$ as one instruction, rounded once. Equation 3 in fused multiply-adds is

$$
b_t = \operatorname{fma}\bigl(c_1, I_{t+1}, \operatorname{fma}(c_2, I_t, c_4 V_t)\bigr), \qquad
Q_{t+1} = \operatorname{fma}(c_3, Q_t, b_t)
$$

Seven floating point operations, four products and three sums, become four instructions: one multiply and three
multiply-adds. Only the last of the four reads $Q_t$. The first three of the next step do not depend on this step, so
the processor computes them while the last multiply-add of this step is still running. One multiply-add per step waits
on the step before.

The order matters even though the algebra does not. The same equation written with $Q_t$ innermost,
$\operatorname{fma}(c_1, I_{t+1}, \operatorname{fma}(c_2, I_t, \operatorname{fma}(c_3, Q_t, c_4 V_t)))$, makes three
multiply-adds wait one after another at every step. A compiler keeps the written order, since reordering
floating point sums changes their rounding. In a Rust version of the kernel routing one step at a time, $Q_t$
innermost took 4.4 s and $c_3 Q_t$ added last 1.5 s. Every time on this page was measured while the kernels were
designed, on one thread of an Apple M3 Max routing the Columbia region's 58,358 rivers through a year of hourly ERA5
runoff (1940).

A multiply-add takes about four cycles before its result can be used. Routing one step at a time therefore takes at
least four cycles per step of every sub-reach, however few other instructions it has.

The inflow at $t$ of one step is the inflow at $t + 1$ of the step before. Keeping it, each step reads the rows of the
inflow rivers once instead of twice.

## SIMD: four lanes in one register

A processor's vector registers hold several values at once, and one instruction does the same operation on each of
them: single instruction, multiple data, or SIMD. Each value's place in the register is a *lane*. The registers every
target of `routing.c` has are 128 bits wide: SSE on x86 processors, NEON on ARM processors, and SIMD128 in
WebAssembly. A float32 is 32 bits, so one register holds four float32 lanes, and one vector multiply-add does four
multiply-adds in the time of one. Wider registers, such as the 256 bits of AVX on many x86 processors, would hold eight
lanes, but WebAssembly has only 128 bit vectors and the same C file is compiled for both, so the kernels use four.

### What a lane holds

The four lanes hold four consecutive steps of one sub-reach, not four rivers. A river's discharge row and a cell's
runoff depths are stored one step after another, so one load reads four consecutive steps of a row into the four lanes
of a register. Four rivers would need four reads from four rows for every vector, and one river may drain into another.

The lanes of one sub-reach do depend on one another, through $Q_t$: lane 1 needs the result of lane 0. Everything that
does not depend on $Q_t$ (the inflow, the runoff, the terms $b_t$, and the partial sums of
[Eight steps at a time](#eight-steps-at-a-time)) is computed in all four lanes at once. Only the multiply-add that
carries $Q_t$ from one group of steps to the next waits.

### How routing.c writes it

`routing.c` declares a vector of four floats with the vector type of GCC and clang,

```c
typedef float four_lanes __attribute__((vector_size(16)));
```

on which `+`, `-`, `*`, and `/` work lane by lane. clang compiles it to the vector instructions of each target, so no
processor's instructions are named in the C. A handful of functions, named the same in the numba kernel's `_simd.py`,
do everything the kernels need with lanes:

```text
broadcast(c)                      ( c,         c,         c,         c        )
load_four(row + t)                ( row[t],    row[t+1],  row[t+2],  row[t+3] )
store_four(row + t, q)            row[t] to row[t+3] = q0, q1, q2, q3
multiply_add(a, b, c)             ( a0 b0+c0,  a1 b1+c1,  a2 b2+c2,  a3 b3+c3 )
shift_in_one_lane(before, x)      ( before3,   x0,        x1,        x2       )
shift_in_two_lanes(before, x)     ( before2,   before3,   x0,        x1       )
last_lane_in_every_lane(x)        ( x3,        x3,        x3,        x3       )
```

- `broadcast` puts a coefficient in every lane, once per sub-reach, so it multiplies four steps at once.
- `load_four` and `store_four` read and write four steps from any step of a row. They do not require a row to start on
  a 16 byte boundary, so a group of steps can start at any step.
- `multiply_add` is fused where the target has a fused multiply-add: every 64 bit ARM processor, x86 processors with
  FMA3 (the native build compiles for the processor it is built on), and WebAssembly with relaxed SIMD. Plain SIMD128
  has none, so there it is a multiply and an add; that build routed 40 to 60% slower than the relaxed SIMD build.
- `shift_in_one_lane` and `shift_in_two_lanes` move each lane up one or two lanes and fill the first lanes from the
  last lanes of another vector. They give the inflow one step earlier without reading the row again, and they line up
  the lanes that the partial sums add together.
- `last_lane_in_every_lane` hands the result of the last of four steps to every lane of the next four.

A real vector type is required. The first C kernel kept its lanes in a `struct { float lane[4]; }`, which clang kept as
four separate floats wherever the struct was carried from one loop iteration to the next, and it took 0.58 s against
0.39 s with the vector type. numba has no vector type and turns off the LLVM optimization that would combine four
floats into one register, so its lanes as tuples took 4.3 s. `_simd.py` builds the same operations from LLVM vector
instructions through llvmlite (`multiply_add` is LLVM's `llvm.fma.v4f32`), and the numba kernel took 0.37 s, within 9%
of the C.

## Terms four steps at a time

Because $b_t$ does not depend on $Q$, the terms of four consecutive steps are computed in the four lanes of one
register. Writing an arrow over a symbol for its four consecutive steps,
$\vec{I}_\text{after} = (I_{t+1}, I_{t+2}, I_{t+3}, I_{t+4})$,
$\vec{I}_\text{before} = (I_t, I_{t+1}, I_{t+2}, I_{t+3})$, and $\vec{V} = (V_t, V_{t+1}, V_{t+2}, V_{t+3})$,

$$
\vec{b} = c_1\, \vec{I}_\text{after} + c_2\, \vec{I}_\text{before} + c_4\, \vec{V}
$$

- $\vec{I}_\text{after}$ is one `load_four` from the row of each inflow river, added together.
- $\vec{I}_\text{before}$ is $\vec{I}_\text{after}$ moved up one lane, with the last inflow of the four steps before in
  the first lane (`shift_in_one_lane`). No row is read twice.
- $\vec{V}$ is, for each weight, one `load_four` of the weight's cell and one `multiply_add` with the weight in
  every lane.
- $\vec{b}$ is then one multiply and two multiply-adds for four steps.

The recurrence itself still waits: four multiply-adds, one after another, carry $Q$ through the four steps. Routing
stays bound by the four cycles each of them waits.

## Eight steps at a time

Applying Equation 3 $j$ times from step $t$ gives

$$
Q_{t+j} = c_3^{\,j}\, Q_t + S_j, \qquad S_j = \sum_{i=0}^{j-1} c_3^{\,j-1-i}\, b_{t+i} \tag{4}
$$

which holds for $j = 1$ and carries from $j$ to $j + 1$, since
$Q_{t+j+1} = c_3 Q_{t+j} + b_{t+j} = c_3^{\,j+1} Q_t + (c_3 S_j + b_{t+j})$. The partial sums $S_j$ follow the same
recurrence from $S_0 = 0$, $S_{j+1} = c_3 S_j + b_{t+j}$, and hold only terms. They do not depend on $Q_t$, so they are
computed before $Q_t$ is known, in lanes.

### Partial sums in lanes

Writing $b_i$ for $b_{t+i}$, the partial sums of four lanes take two shifted multiply-adds (`scan_four_steps`):

$$
\begin{aligned}
\vec{s} &= \vec{b} + c_3\,(0,\ b_0,\ b_1,\ b_2) = (b_0,\ b_1 + c_3 b_0,\ b_2 + c_3 b_1,\ b_3 + c_3 b_2) \\
\vec{S} &= \vec{s} + c_3^2\,(0,\ 0,\ s_0,\ s_1)
= (b_0,\ b_1 + c_3 b_0,\ b_2 + c_3 b_1 + c_3^2 b_0,\ b_3 + c_3 b_2 + c_3^2 b_1 + c_3^3 b_0)
\end{aligned}
$$

so lane $j$ holds $S_{j+1}$. Each lane adds $c_3$ times the lane before it, then $c_3^2$ times the sum two lanes before
it: four lanes take $\log_2 4 = 2$ shifted multiply-adds (Hillis and Steele, 1986). The second four steps take their
own partial sums the same way, then add $S_4$ carried through them, $S_{4+j} = c_3^{\,j} S_4 + {}$ the $j$-th partial
sum of $b_4$ to $b_7$: one multiply-add of $(c_3, c_3^2, c_3^3, c_3^4)$ by $S_4$ in every lane.

### The discharge of eight steps

With the partial sums known, Equation 4 gives all eight discharges from $Q_t$ in two multiply-adds:

$$
\begin{aligned}
(Q_{t+1},\ Q_{t+2},\ Q_{t+3},\ Q_{t+4}) &= (c_3,\ c_3^2,\ c_3^3,\ c_3^4)\; Q_t + (S_1,\ S_2,\ S_3,\ S_4) \\
(Q_{t+5},\ Q_{t+6},\ Q_{t+7},\ Q_{t+8}) &= (c_3^5,\ c_3^6,\ c_3^7,\ c_3^8)\; Q_t + (S_5,\ S_6,\ S_7,\ S_8)
\end{aligned}
$$

with $Q_t$ in every lane. The last lane, $Q_{t+8}$, goes into every lane for the next eight steps. The powers of $c_3$
are computed once per sub-reach.

Only these two multiply-adds read $Q_t$, and the next eight steps need only the second. One multiply-add and one
shuffle per eight steps wait on the steps before, against eight multiply-adds one step at a time. The loads, terms, and
partial sums of the next eight steps do not depend on $Q$, so the processor computes them while that multiply-add runs.
Routing one step at a time took 1.40 s, and eight steps at a time 0.34 s, 4.1 times faster. Any steps left after the
last whole eight are routed one step at a time by Equation 3.

The eight-step form adds the same products in a different order, so its discharge differs from one step at a time only
by float32 rounding. With non-negative coefficients (see [Signs of the coefficients](#signs-of-the-coefficients)),
every power of $c_3$ is from 0 to 1, and with non-negative inflow and runoff every product is non-negative, so no sum
loses digits to cancellation.

## Operation counts

Counted from the source for one sub-reach over eight steps, with $N_u$ inflow rivers and $N_w$ weights. A vector
instruction counts as one. Reading the index of each inflow river and the cell and weight of each weight is left out;
it favors one step at a time, which reads them at every step rather than once per eight steps.

| Work, per sub-reach and eight steps | As written (Equation 1) | Equation 3, $c_3 Q$ last | Eight steps at a time |
|-------------------------------------|-------------------------|--------------------------|-----------------------|
| Read inflow                         | $16N_u$ loads, $16N_u$ adds | $8N_u$ loads, $8N_u$ adds | $2N_u$ loads, $2N_u$ adds |
| Read runoff                         | $8N_w$ loads, $8N_w$ multiplies, $8N_w$ adds | $8N_w$ loads, $8N_w$ multiply-adds | $2N_w$ loads, $2N_w$ multiply-adds |
| Terms and recurrence                | 32 multiplies, 24 adds  | 8 multiplies, 24 multiply-adds | 2 multiplies, 11 multiply-adds, 8 shuffles |
| Write discharge                     | 8 stores                | 8 stores                 | 2 stores              |
| Total, $N_u = N_w = 2$              | 176                     | 104                      | 39                    |
| Wait on the step before             | 24 (a multiply and two adds per step) | 8 multiply-adds | 1 multiply-add, 1 shuffle |

The eight-step column is two multiplies and four multiply-adds for the terms of two vectors, four multiply-adds for
their partial sums, one to carry $S_4$ into the second four, and two for the discharge, with two shuffles for the
inflow before each step, four in the partial sums, and two to put one lane in every lane.

The measured gain is larger than the 2.7 times fewer instructions from one step at a time to eight, because one step
at a time is bound by waiting, not by instructions. A current processor starts several instructions every cycle, so
its 104 instructions take fewer cycles than its eight multiply-adds that wait one after another, about 32 cycles.
Eight steps at a time waits about 6 cycles per eight steps and is bound by its 39 instructions instead.

## Folding constants into the coefficients

Any product of values that are the same at every step is computed once per run instead of once per step.

- **The runoff coefficient.** A sub-reach taking a share $w$ of its river's runoff has $c_4 = (c_1 + c_2)\, w / \Delta t$:
  the weight of the runoff in the Muskingum step, the river's share, and the conversion from a volume to a flow in one
  coefficient, one multiply per step instead of three.
- **Runoff weights.** A weight is the volume a meter of runoff depth in its cell gives the river, in m³, which is the
  area of the river's catchment inside the cell. The runoff volume is $V_t = \sum_c w_c\, d_{c,t}$, one multiply-add
  per weight per step, with no area or unit conversion left to do. Runoff already aggregated to catchments is a grid
  with one cell per river, and channel routing is a network with no weights.

## Folding subcycles into one step

A river too short for the runoff step is routed in $s$ subcycles of its own within each step (see
[Signs of the coefficients](#signs-of-the-coefficients)), with its inflow rising linearly from $I_t$ to $I_{t+1}$
across them and its runoff spread evenly over them. Subcycle $j$ of $s$ ends at the fraction $j / s$ of the step, where
the inflow is $I(j/s) = (1 - j/s)\, I_t + (j/s)\, I_{t+1}$:

$$
Q_j = c_3\, Q_{j-1} + c_1\, I(j/s) + c_2\, I\bigl((j-1)/s\bigr) + c_4\, V_t
$$

This is Equation 3 with $s$ steps, so Equation 4 gives $Q_s$ from $Q_0 = Q_t$. Collecting the products of $I_t$,
$I_{t+1}$, and $V_t$ gives one step of the whole runoff step:

$$
\begin{gathered}
Q_{t+1} = c_1'\, I_{t+1} + c_2'\, I_t + c_3'\, Q_t + c_4'\, V_t, \qquad c_3' = c_3^{\,s}, \\
c_1' = \sum_{j=1}^{s} c_3^{\,s-j}\, \frac{c_1\, j + c_2\, (j - 1)}{s}, \qquad
c_2' = \sum_{j=1}^{s} c_3^{\,s-j}\, \frac{c_1\, (s - j) + c_2\, (s - j + 1)}{s}, \qquad
c_4' = c_4 \sum_{j=1}^{s} c_3^{\,s-j} = c_4\, \frac{1 - c_3^{\,s}}{1 - c_3}
\end{gathered} \tag{5}
$$

The folded coefficients keep every property of a subcycle's:

- $c_1' + c_2' = (c_1 + c_2) \sum_j c_3^{\,s-j} = (1 - c_3)\, \frac{1 - c_3^{\,s}}{1 - c_3} = 1 - c_3'$, so they sum
  to one.
- $c_4' = (c_1' + c_2')\, w / \Delta t$, the same relation as before folding.
- None is negative where none of the subcycle's is.

`static_muskingum._combine_subcycles` computes them once per run in float64. The kernels never see subcycles: a river
of $s$ subcycles costs one step, not $s$, and Equation 4 applies to it unchanged. At an hourly runoff step with
`dt_routing: auto`, 830 of the Columbia region's 58,358 rivers take subcycles, one of them 210. Folded, the region
routes 84,007 sub-reach steps each hour instead of 89,855. In the Amazon region, 6,295 rivers take up to 141
subcycles, and folding routes 448,165 sub-reach steps each hour instead of 468,491.

## Sub-reaches in series

A river too long for its step is split into $N$ sub-reaches in series, and those do not fold. Eliminating the first
sub-reach's discharge from the second's recurrence leaves a recurrence that reaches back two steps; $N$ sub-reaches make
one that reaches back $N$ steps, with coefficients built from sums and products of theirs, and Equation 4 no longer
applies. A river of $N$ sub-reaches costs $N$ recurrences, so the kernels route its sub-reaches one after another, each
over every step. The Columbia region's 84,007 sub-reaches at an hourly step are 1.44 per river.

They are routed in place in the river's own discharge row. The first sub-reach reads the rows of the river's inflow
rivers and writes the river's row. Each later sub-reach reads the row, the discharge the sub-reach before it left, and
writes its own over it. That works because each group of eight steps loads its inflow at all eight steps before it
stores anything, and the inflow before its first step comes from the last lane of the group before, carried in a
register, since the row already holds the sub-reach's own discharge there. A river needs one row however many
sub-reaches it has, and `channel_state` holds each sub-reach's discharge at the end of the steps routed. Chaining the
sub-reaches in registers instead, each eight steps through every sub-reach in turn, took 0.52 s against 0.34 s, since
each eight steps then waited on the whole chain.

## Memory

### What is stored

Only two arrays grow with the number of steps $T$:

- the discharge, a row of $T + 1$ float32 values for every river: the output;
- the runoff depth, $T$ float32 values for every cell that a weight reads: the input.

Routing the Amazon region's 302,199 rivers through a year of hourly steps writes 10.6 GB of discharge. Everything else
is sized by rivers, sub-reaches, or weights: the inflow river lists, four coefficients and one `channel_state` value
per sub-reach, and a cell and a weight per weight.

### What is never stored

Each form above computes a quantity where it is used instead of storing it:

| Quantity                       | Stored, it would take | What the kernels keep instead                       |
|--------------------------------|-----------------------|-----------------------------------------------------|
| Inflow $I$                     | $n \times T$          | the inflow rivers' rows, read eight steps at a time |
| Runoff volume $V$              | $n \times T$          | weights times cell depths, eight steps at a time    |
| Terms $b$ and partial sums $S$ | $n \times T$          | two vectors of four lanes                           |
| Sub-reach discharge            | sub-reaches × $T$     | the river's own row, written over in place          |

Computing them in place also costs less time. Summing the inflow into a row of its own before routing took 11% longer
than reading the inflow rows in the eight-step loop. Aggregating the runoff to each river first made the kernel 12%
faster, 0.30 s against 0.34 s, but the aggregating took 0.55 s.

### Options that store less

**Route fewer steps at once.** `channel_state` carries everything a river needs from one set of steps to the next, so
the discharge needs rows only as long as the steps being routed: $n \times (W + 1)$ values for $W$ steps. The Router
already routes each runoff file as its own set of steps, into the memory of the file before. Shorter sets bound the
memory further at a cost per river for each set: routing every river over windows of 64 steps before the next window
took about twice as long as routing every step at once.

**Release each row once it has been read.** In routing order, a river's row is read only by the river it drains into,
and only when that river is routed. A river can write its own discharge over the row of one of its inflow rivers, the
way a sub-reach writes over the row it reads. The rows that must be held at once are only those written and not yet
read. In the depth first order of their network files, the Columbia region's 58,358 rivers never need more than 19 rows
at once, and the Amazon region's 302,199 rivers never more than 18. For the Amazon region's year of hourly steps, that
is 631 KB instead of 10.6 GB. Routing first the inflow river whose watershed needs the most rows lowers it to the least
any order can reach (Sethi and Ullman, 1970), which for confluences of two rivers is the network's Strahler order
(Ershov, 1958; Flajolet et al., 1979). These networks' Strahler orders are 8 and 9; their least are 17 and 16, set by confluences of 17 and 16
rivers. The kernels write every river's row because the rows are the output. Releasing rows pays when only some
rivers' discharge is kept, such as outlets or gauges, or when each row is written to the output file as soon as it is
complete.

**Keep only the rivers that are needed.** Rivers that contribute nothing to the results do not have to be in the
network, and with released rows, rivers that contribute but are not reported cost no memory beyond the rows held at
once.

## Signs of the coefficients

The coefficients share a denominator, $u + 1 + v$, positive for any positive step, so from Equation 2 each is
non-negative when its numerator is:

$$
\begin{aligned}
c_1 \ge 0 &\iff u \ge 1 - v \iff \Delta t \ge 2kx \\
c_2 \ge 0 &\iff u \ge v - 1, \text{ which holds for every } x \ge 0 \\
c_3 \ge 0 &\iff u \le 1 + v \iff \Delta t \le 2k(1 - x)
\end{aligned}
$$

Every coefficient is non-negative in the window $2kx \le \Delta t \le 2k(1-x)$: from 0 to $2k$ at $x = 0$, and only
$\Delta t = k$ at $x = 0.5$. The window spans a factor of $(1 - x) / x$ in $k$, 4 at $x = 0.2$, so a network whose
travel times span more than that has no one step at which every river is inside it. A negative coefficient gives a
negative response to a positive inflow. Whatever the step, a river's response to its inflow has mean $k$ and variance
$k^2 (1 - 2x)$, so a shorter step is not a more accurate one. Muskingum $x$ is never negative.

**A river too long for its step**, $\Delta t < 2kx$, has $c_1 < 0$ and is split into the fewest sub-reaches that keep
$c_1$ non-negative, $N = \lceil k / (\Delta t + (1 - 2x)k) \rceil$, whose travel times are the shares $w_i$ of $k$. They
keep the river's mean travel time $k$, and with the common

$$
x_N = \frac12 - \frac{\tfrac12 - x}{\sum_i w_i^2}
$$

they keep its variance, since $\sum_i (w_i k)^2 (1 - 2x_N) = k^2 (1 - 2x)$. Equal shares give
$x_N = \tfrac12 - N(\tfrac12 - x)$, used wherever it is not negative. Where it would be, the sub-reaches take $x_N = 0$
and unequal shares, the first $N - 1$ equal and larger than the last, whose squares sum to $1 - 2x$: 0.724 and 0.276 of
$k$ for two sub-reaches at $x = 0.2$. Each sub-reach takes its share of the river's runoff.

**A river too short for its step**, $\Delta t > 2k(1 - x)$, has $c_3 < 0$ and is routed in
$s = \lceil \Delta t / (2k(1 - x)) \rceil$ subcycles, which Equation 5 folds into one step.

Non-negative coefficients also keep the fast forms accurate: every power of $c_3$ in Equations 4 and 5 is from 0 to 1,
and every product they sum is non-negative. `river-route` writes negative discharge as routed; a river with negative
discharge is a river outside its window, which sub-reaches or subcycles treat (see [Time Variables](time-options.md)).

## References

- Cunge, J. A. (1969). On the subject of a flood propagation computation method (Muskingum method). *Journal of
  Hydraulic Research*, 7(2), 205–230. [doi:10.1080/00221686909500264](https://doi.org/10.1080/00221686909500264)
- David, C. H. (2011). River Network Routing on the NHDPlus Dataset. *Journal of Hydrometeorology*.
  [doi:10.1175/2011JHM1345.1](https://doi.org/10.1175/2011JHM1345.1)
- HEC-HMS Users Manual introduction to the Muskingum Model
  <https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/channel-flow/muskingum-model>
- HEC-HMS Users Manual introduction to the Muskingum Cunge Model
  <https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/channel-flow/muskingum-cunge-model>
- Hillis, W. D., and Steele, G. L. (1986). Data parallel algorithms. *Communications of the ACM*, 29(12), 1170–1183.
  [doi:10.1145/7902.7903](https://doi.org/10.1145/7902.7903)
- Blelloch, G. E. (1990). *Prefix Sums and Their Applications*. Technical Report CMU-CS-90-190, Carnegie Mellon
  University.
- Ershov, A. P. (1958). On programming of arithmetic operations. *Communications of the ACM*, 1(8), 3–6.
  [doi:10.1145/368892.368907](https://doi.org/10.1145/368892.368907)
- Flajolet, P., Raoult, J.-C., and Vuillemin, J. (1979). The number of registers required for evaluating arithmetic
  expressions. *Theoretical Computer Science*, 9(1), 99–125.
  [doi:10.1016/0304-3975(79)90009-4](https://doi.org/10.1016/0304-3975(79)90009-4)
- Sethi, R., and Ullman, J. D. (1970). The generation of optimal code for arithmetic expressions. *Journal of the
  ACM*, 17(4), 715–728. [doi:10.1145/321607.321620](https://doi.org/10.1145/321607.321620)
- Strahler, A. N. (1957). Quantitative analysis of watershed geomorphology. *Transactions, American Geophysical
  Union*, 38(6), 913–920. [doi:10.1029/TR038i006p00913](https://doi.org/10.1029/TR038i006p00913)
- Wikipedia: [Topological sorting](https://en.wikipedia.org/wiki/Topological_sorting).
