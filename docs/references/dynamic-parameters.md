# Dynamic Parameters

Static parameters fix each river's Muskingum $k$ and $x$ for the whole run. Dynamic parameters change them with the
river's discharge: every routing step, a dynamic method gives each river's $k$ and $x$ at its discharge, and its
coefficients are rebuilt from them. This is the variable parameter Muskingum method. `river-route` has one dynamic
method, Muskingum-Cunge, chosen with `parameters='cunge'`.

This page describes the equations of the method, where its parameters come from, how the kernel routes a river with
them and keeps every coefficient non-negative while they change, and what that costs in accuracy, mass, and time.

---

## Muskingum-Cunge k and x

Cunge (1969) showed that the Muskingum method solves the diffusion wave equation when its numerical diffusion is made
equal to the physical diffusion of the wave. For a river of length $L$ and bed slope $S_0$, carrying discharge $Q$ with
a top width $T$ and a kinematic wave celerity $c = dQ/dA$, that gives

$$
k = \frac{L}{c}, \qquad x = \frac{1}{2} - \frac{Q}{2\, T\, S_0\, c\, L}.
$$

Both change with discharge through $c$ and $T$. When the celerity and the top width of a river are power laws of
discharge, its hydraulic geometry (Leopold and Maddock, 1953), $k$ and $\tfrac12 - x$ are power laws of discharge too.
Every river is routed with that one form:

$$
k = \alpha\, Q^{\beta}, \qquad x = \frac{1}{2} - \gamma\, Q^{\delta},
$$

read from six columns of the network file:

| Column              | Symbol   | Meaning                                         |
|---------------------|----------|-------------------------------------------------|
| `cungeKCoefficient` | $\alpha$ | k in seconds at a discharge of 1 m³/s           |
| `cungeKExponent`    | $\beta$  | how k changes with discharge                    |
| `cungeXCoefficient` | $\gamma$ | $\tfrac12 - x$ at a discharge of 1 m³/s         |
| `cungeXExponent`    | $\delta$ | how $\tfrac12 - x$ changes with discharge       |
| `cungeMinQ`         |          | the least discharge k and x are evaluated at    |
| `cungeMaxQ`         |          | the greatest discharge k and x are evaluated at |

Three bounds hold the laws to values that can be routed:

- **Discharge** is held within `cungeMinQ` and `cungeMaxQ`, the discharges the regression was fit over, so neither law
  is extrapolated, and the power of a negative or zero discharge is never taken.
- **x is never negative.** Below 0 the weight of the inflow at the start of each step, $c_2$, is negative at every step
  shorter than $-2kx$ (see [Signs of the coefficients](math.md#signs-of-the-coefficients)).
- **x is at most 0.45** with `restrict_x_range=True`, the default, or 0.5 without it. See
  [Why x is restricted](#why-x-is-restricted).

Static parameters are the special case $\beta = \delta = 0$, with $\alpha = k$ and $\gamma = \tfrac12 - x$: a river
with exponents of 0 routes with the same coefficients at every step.

## Where the parameters come from

The regression is fit outside `river-route`, by whoever prepares the network file, which then only reads its columns.
The fit described here is the one made for the Columbia River of the test data from the GEOGLOWS flood maps, which hold
for every river a representative cross section built by the Automated Rating Curve (ARC) model and the discharge range
its flood maps were made for.

1. **The stage table.** At depth stages of 0.1 m, the representative cross section gives the discharge $Q$, the flow
   area $A$, and the top width $T$, and the reach its slope $S_0$. The length $L$ is the river's length in the network.
2. **k and x at each stage.** The celerity is the slope of the stage table, $c = dQ/dA$, from which
   $k = L / c$ and $x = \tfrac12 - Q / (2 T S_0 c L)$. Stages whose discharge or area does not rise above every stage
   below them are dropped, so $c$ is positive.
3. **The discharge range.** The regression is fit from the river's baseflow to the largest discharge of its flood maps,
   within the discharges of its stage table. Below the baseflow the stage table describes the narrow channel ARC
   estimates under the water surface of the elevation model, not measured terrain, so it is left out.
4. **The regression.** $\log k$ and $\log(\tfrac12 - x)$ are fit as straight lines in $\log Q$ by least squares over 50
   discharges spread evenly in $\log Q$ across the range. The line of $\log(\tfrac12 - x)$ is constrained to keep $x$
   from 0 to 0.45: a straight line in $\log Q$ is set by its values at the two ends of the range, so holding both ends
   within the bounds holds every discharge between, and the constrained least squares is a problem in those two values.
5. **Rivers without data** keep their static $k$ and $x$ as exponents of 0, with `cungeMinQ` and `cungeMaxQ` of 1.

On the Columbia, 26,352 of its 29,405 rivers were fit; the rest are in Canada, beyond the flood maps, or are left out
of the flood maps, including much of the main stems of the Columbia and Snake. Fit to the $k$ of each stage table, the
power law is off by a median of 4% (18% at the 95th percentile), against 16% for a constant $k$. A quadratic in
$\log Q$ would be off by 3%, but its extra parameters extrapolate poorly and the gain is small, so every river uses the
power law.

## How the kernel routes a river

Rivers are routed one at a time from upstream to downstream, each river's whole time series before the next (see
[Parallelism](parallelism.md)). So when a river is routed, the inflow from every river upstream of it is known for every
step of the run, and so is its own runoff. With `force_validity`, the default, each river is routed as the sub-reaches
in series that [Forcing validity](#forcing-validity) lays out, each sub-reach's whole series before the next, and each
routing step of a sub-reach in the subcycles its $k$ and $x$ need at that step. A river that needs neither is one
sub-reach of one subcycle.

### The discharge k and x are evaluated at

Each step's $k$ and $x$ are evaluated at the discharge entering the sub-reach over the step: the mean of its inflow at
the start and end of the step plus its share of the river's runoff,

$$
Q_t = \frac{I_t + I_{t+1}}{2} + w\,L_t,
$$

where $L_t$ is the river's runoff volume over the step divided by its length in seconds and $w$ is the sub-reach's
share of the river, 1 for a river routed whole. Every $Q_t$ is known before the sub-reach is routed, so $k$ and $x$ are
evaluated for its whole series at once.

### The coefficients

Written with $u = \Delta t / k$ and $v = 1 - 2x$, the Muskingum coefficients of each step, or of each of its $s$
subcycles with $u$ divided by $s$, are

$$
\begin{gathered}
c_1 = \frac{u - 1 + v}{u + 1 + v}, \qquad
c_2 = \frac{u + 1 - v}{u + 1 + v}, \qquad
c_3 = \frac{1 + v - u}{u + 1 + v}, \\
Q_{t+1} = c_{3,t}\, Q_t + c_{1,t}\, I_{t+1} + c_{2,t}\, I_t + (c_{1,t} + c_{2,t})\, w\, L_t.
\end{gathered}
$$

They sum to 1 at every step, as the static coefficients do, so a river whose inflow and runoff are constant is held at
a constant discharge. All three are non-negative when $1 - v \le u \le 1 + v$, which is
$2kx \le \Delta t \le 2k(1-x)$.

### The steps of the kernel

For each river, the kernel writes its runoff rate $L_t$ at every routing step from its forcing, and then for each of
its sub-reaches, upstream to downstream:

1. writes the discharge $Q_t$ of every step, held within `cungeMinQ` and `cungeMaxQ`;
2. evaluates the river's $u$ and $v$ at every step from the laws as base 2 logs,
   $u = 2^{\log_2(\Delta t / \alpha) - \beta \log_2 Q}$ and $v = 2^{\log_2(2\gamma) + \delta \log_2 Q}$, with $v$
   held from $1 - 2 x_{max}$ to 1, where $x_{max}$ is 0.45, or 0.5 without `restrict_x_range`. The logs of the
   parameters are taken once per run, when the routing step is known;
3. turns them into the sub-reach's: $u / w$, and $v$ scaled by its layout, held at 1 or less;
4. writes the subcycles, $c_1$, and $c_3$ of every step, with $c_2 = 1 - c_1 - c_3$;
5. solves the outflow of every step through its subcycles, with the inflow interpolated linearly across them, the
   only part where each step waits on the one before, and passes it on as the inflow of the next sub-reach;

and adds the outflow of its last sub-reach into the inflow of the river downstream and writes its discharge at the
end of each runoff step, as static routing does.

Steps 1 to 4 have no dependence between steps, and the compiler evaluates several steps at a time with vector
instructions. numba calls `np.log2` and `np.exp2` one value at a time, so the kernel computes them itself: $\log_2$ as
the exponent of a float32 plus a polynomial of its mantissa, within $9 \times 10^{-7}$, and $2^y$ as a polynomial of the
fraction of $y$ times 2 to its whole part, built in the exponent bits of a float32, within $3 \times 10^{-7}$ of its
value. Those errors are far smaller than the error of the regression itself.

## Forcing validity

With static parameters a river's window of valid steps is fixed, so a river too long for the routing step is split
into sub-reaches and a river too short for it takes subcycles, once (see
[Signs of the coefficients](math.md#signs-of-the-coefficients)). With dynamic parameters the window moves with
discharge, and one river can be too long at low discharge, where its $k$ is large, and too short at high discharge.

- **Sub-reaches are planned once, for every discharge.** Each sub-reach holds its own state, so a river's sub-reaches
  are fixed for the run. Discharge is held within `cungeMinQ` and `cungeMaxQ`, so every $k$ and $x$ the river can take
  is known before routing. The plan is the fewest sub-reaches for which some whole number of subcycles keeps every
  coefficient non-negative at each of 64 discharges spread evenly in $\log Q$ across the range.
- **The sub-reaches keep the river's variance.** Sub-reaches of shares $w_i$ of $k$ that all take
  $x_N = \tfrac12 - (\tfrac12 - x)/\sum_i w_i^2$ keep the river's mean travel time and its variance at every discharge.
  The shares are set by the river's smallest $x$ across its range: equal wherever $x_N$ stays at 0 or more there, and
  otherwise the unequal shares static sub-reaches take, whose squares sum to $1 - 2x$.
- **Subcycles are chosen at every step.** They hold no state, so each step of each sub-reach takes the fewest subcycles
  that keep its $c_3$ non-negative, $\lceil u / (1 + v) \rceil$, for its own $k$ and $x$.
- **Every step is valid.** Between the discharges the plan checks, a step's window can hold no whole number of
  subcycles, and those subcycles leave its $c_1$ negative. That step's $x$ is lowered just enough to make $c_1$ zero.
  The plan makes this rare and small, and with it no coefficient of any step is ever negative.
- **Rivers whose x changes widely.** Shares that keep a river's variance at its smallest $x$ can leave a sub-reach too
  long at its largest. A river with no valid split that keeps its variance is split into equal sub-reaches whose $x$ is
  held at 0 where keeping the variance would take it lower, and a river with no valid split at all is routed whole with
  its $x$ lowered at the steps that need it. Both are reported as `invalid_coefficients` says.

On the Columbia at a one hour step, with $x$ held at 0.45 or less, the plan splits its 29,405 rivers into 66,426
sub-reaches:

| Sub-reaches per river | 1     | 2      | 3     | 4     | 5   | 6 to 9 |
|-----------------------|-------|--------|-------|-------|-----|--------|
| Rivers                | 5,557 | 15,171 | 4,858 | 3,363 | 302 | 154    |

1,061 rivers take unequal shares, 135 rivers whose $x$ changes widely take equal shares with $x$ held at 0 at some
discharges, and none is routed whole. The most subcycles any discharge of a river needs is 3 at the median, 16 at the
99th percentile, and 451 at most. Routing March to June 2000, no discharge of any river is negative.

## Choosing the discharge k and x are evaluated at

Any discharge of the reach over the step could set its $k$ and $x$, and each is an approximation of where along the
reach and when during the step the wave's celerity applies. Ponce and Yevjevich (1978) and Ponce and Chaganti (1994)
evaluate them at averages of the inflow and outflow. Any choice that involves the outflow of the step is known only once
the step before it is solved, so every step would wait on the logs and powers of the one before. The inflow is known
for the whole series, which is what lets the kernel evaluate every step at once and leaves only the multiply-adds of
the outflow waiting on each other. On a sharp rise, inflow exceeds outflow, so the inflow gives a higher discharge and,
as $k$ falls with discharge, a slightly shorter travel time than the outflow would.

The routing step itself changes the answer a little, since $k$ and $x$ are evaluated once per routing step. Routing the
Columbia with forced validity at one hour and at shorter steps:

| Compared with a one hour step                               | 15 minutes | 5 minutes |
|-------------------------------------------------------------|------------|-----------|
| Outlet Nash-Sutcliffe efficiency, cunge                     | 0.999994   | 0.999993  |
| Outlet peak, cunge                                          | +0.08%     | +0.09%    |
| Nash-Sutcliffe efficiency per river, median, cunge          | 0.998      | 0.998     |
| Nash-Sutcliffe efficiency per river, 1st percentile, cunge  | 0.975      | 0.971     |
| Nash-Sutcliffe efficiency per river, 1st percentile, static | 0.9997     | 0.9997    |

## Conservation of mass

Static Muskingum routing conserves mass: with $k$ and $x$ fixed, the storage $S = k(xI + (1-x)Q)$ that one step ends
with is the storage the next begins with. When $k$ and $x$ change between steps, the storage a step's coefficients end
with is not the storage the next step's coefficients begin with, and the difference is water created or lost. This is a
property of every variable parameter Muskingum method; Todini (2007) derives a variant that conserves it.

Routing the Columbia's runoff of March to June 2000 from an empty network with forced validity, and then the channel
alone for 120 days until it holds no water:

| Parameters | Water out of the outlet less the runoff in |
|------------|--------------------------------------------|
| static     | +0.00003%, the rounding of float32         |
| cunge      | +0.018%                                    |

## Speed

On the Columbia, 29,405 rivers routed with forced validity for the 2,928 hours of March to June 2000 with grid runoff,
on an Apple M3 Max:

| Parameters | ns per river and hour | 1 thread | 4 threads | 8 threads |
|------------|-----------------------|----------|-----------|-----------|
| static     | 1.3                   | 0.108 s  | 0.035 s   | 0.025 s   |
| cunge      | 17                    | 1.48 s   | 0.42 s    | 0.23 s    |

These are the times of the routing kernels alone. Reading the grid runoff took about one second more in every case.
Each hour of a cunge river is 2.3 sub-reaches on average, each with its $k$ and $x$ evaluated and its subcycles solved.

## Why x is restricted

Every coefficient of a reach is non-negative only when its step is within $(1 - 2x)\,k$ of its travel time:
$2kx \le \Delta t \le 2k(1-x)$. The Muskingum-Cunge $x$ of a long, steep reach is close to $\tfrac12$: on the Columbia
the median $1 - 2x$ of the stage tables at their largest discharge is 0.012, so a step must match a reach's travel time
within about half a minute, while that travel time changes with discharge. Forcing validity at such an $x$ would take
about 160 times as many sub-reach steps as routing each river whole. Holding $x$ at 0.45 or less widens the window to
10% of $k$ either side and brings that to the 2.3 sub-reaches and few subcycles above. It adds diffusion to the reaches
whose Muskingum-Cunge $x$ is above 0.45: their travel time variance becomes at least a tenth of $k^2$.

## References

- Cunge, J. A. (1969). On the subject of a flood propagation computation method (Muskingum method). *Journal of
  Hydraulic Research*, 7(2), 205-230.
- Leopold, L. B., and Maddock, T. (1953). The hydraulic geometry of stream channels and some physiographic
  implications. *U.S. Geological Survey Professional Paper* 252.
- Ponce, V. M., and Yevjevich, V. (1978). Muskingum-Cunge method with variable parameters. *Journal of the Hydraulics
  Division*, 104(12), 1663-1667.
- Ponce, V. M., and Chaganti, P. V. (1994). Variable-parameter Muskingum-Cunge method revisited. *Journal of
  Hydrology*, 162(3-4), 433-439.
- Todini, E. (2007). A mass conservative and water storage consistent variable parameter Muskingum-Cunge approach.
  *Hydrology and Earth System Sciences*, 11(5), 1645-1659.
