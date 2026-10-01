# Mathematical Derivations

Derivations of the Muskingum routing equation and alternative methods for the most efficient methods for solving the linear system.

Note that the `river-route` implementation uses float32 precision for speed and efficiency since discharge and runoff measurements
are not known to a level needing 64-bit precision during the computation.

---

## Summary

Some key insights applying linear algebra and graph theory to river networks and the Muskingum equation are:

1. River networks are directed acyclic graphs (DAGs).
2. River segments can be topologically sorted so that upstream always comes before downstream.
3. Topological ordering makes the adjacency matrix $A$ strictly lower triangular.
4. River network adjacency matrices are extremely sparse with exactly one nonzero entry per column (except for outlets).
5. The Muskingum equation LHS $\mathbf{I} - c_1 A$ is unit lower triangular. The identity contributes ones on the
   diagonal; $c_1 A$ contributes entries only below.
6. Unit lower triangular systems are best solved with forward substitution.

## River Network Ordering

Rivers are often described as "networks" or "systems". When being modeled, river networks have a few properties that are useful
to take advantage of for mathematically more efficient algorithms.

1. They are "directed" -- meaning water only flows in one direction from upstream to downstream.
2. They are "acyclic" -- meaning there are no loops because water cannot flow upstream.
3. They are "dendritic" -- meaning they branch out when going upstream and merge when going downstream (ignoring braided rivers and deltas, for instance).

Rivers can be topologically sorted. Rather than sorting them from high to low by an attribute or an ID, topologically sorting means
sorting them in the order they are connected in the network. That is, from "upstream to downstream". The further upstream a river is
the earlier it should appear in the sorted list. A useful tool for conceptualizing and diagramming this is the Strahler stream order.
The Strahler order assigns the number 1 to the most upstream, or headwater, segments. When two segments of the same order merge, the
downstream segment is assigned an order 1 higher. If two different order merge, the downstream segment is assigned the higher of the
two inlet orders. Streams that are a headwater area have no upstream segments.

Some river datasets will have multiple segments in a row which have the same river order. In those cases, you could sort rivers of
the same order by increasing cumulative drainage area or another attribute that increases as you go downstream. There are multiple
valid ways to sort rivers which are all topologically sorted. It is not unique. The only requirement is that upstream segments
appear before downstream segments in the sorted list.

### Breadth First Search (BFS)

<div style="text-align: center;">

```mermaid
graph TD
    R1@{ shape: sm-circ } -->|" #1 · Order 1 "| R5@{ shape: sm-circ }
    R2@{ shape: sm-circ } -->|" #2 · Order 1 "|R5
    R3@{ shape: sm-circ } -->|" #3 · Order 1 "|R6@{ shape: sm-circ }
    R4@{ shape: sm-circ } -->|" #4 · Order 1 "|R6
    R5 -->|" #5 · Order 2 "|R7@{ shape: sm-circ }
    R6 -->|" #6 · Order 2 "|R7
    R8@{ shape: sm-circ } -->|" #8 · Order 1 "|R9@{ shape: sm-circ }
    R7 -->|" #7 · Order 3 "|R9
```

<figcaption><em>Figure 1: A topologically sorted river network labeled with Strahler stream orders.</em></figcaption>
</div>

In the diagram above, rivers 1 through 4 and 8 are headwaters with no upstream dependencies. Rivers 5 and 6 each receive two headwater
tributaries and are indexed after their upstream sources. River 7 merges two second-order streams and river 9, the outlet, appears last.
Another way to describe rivers that are topologically sorted is that they are sorted in order of independence. Segments at the top of
the list depend on no rivers and rivers further down the list depend on a greater number of upstream segments to get their inflow.
River 5's inflow depends on what is discharged from rivers 1 and 2. A river's discharge cannot be computed until all of upstream
contributors are known.

### Depth first search (DFS) order

A depth first search finds a topological order. It starts at each outlet and walks upstream, adding a river to the
order only after every river upstream of it has been added. Its order also places the rivers
upstream of each river in the rows immediately before it, so every river's whole upstream watershed is one contiguous
range of rows ending at that river. `river-route` requires this order (see the
[network file](io-file-schema.md#network-file) requirements).

<div style="text-align: center;">

```mermaid
---
config:
  themeVariables:
    fontSize: 12px
---
graph TD
    S0@{ shape: sm-circ } -->|" idx: 0<br>count: 0 "| S2@{ shape: sm-circ }
    S1@{ shape: sm-circ } -->|" idx: 1<br>count: 0 "| S2
    S3@{ shape: sm-circ } -->|" idx: 3<br>count: 0 "| S5@{ shape: sm-circ }
    S4@{ shape: sm-circ } -->|" idx: 4<br>count: 0 "| S5
    S2 -->|" idx: 2<br>count: 2 "| S6@{ shape: sm-circ }
    S5 -->|" idx: 5<br>count: 2 "| S6
    S6 -->|" idx: 6<br>count: 6 "| S8@{ shape: sm-circ }
    S7@{ shape: sm-circ } -->|" idx: 7<br>count: 0 "| S8
    S9@{ shape: sm-circ } -->|" idx: 9<br>count: 0 "| S12@{ shape: sm-circ }
    S10@{ shape: sm-circ } -->|" idx: 10<br>count: 0 "| S12
    S11@{ shape: sm-circ } -->|" idx: 11<br>count: 0 "| S12
    S8 -->|" idx: 8<br>count: 8 "| S13@{ shape: sm-circ }
    S12 -->|" idx: 12<br>count: 3 "| S13
    S13 -->|" idx: 13<br>count: 13 "| OUTLET@{ shape: sm-circ }
    linkStyle default stroke-width:3px
```

<figcaption><em>Figure 2: A river network in DFS order, each river labeled with its <code>riverIndex</code> and <code>upstreamCount</code>.</em></figcaption>
</div>

A river's watershed is the rows from its `riverIndex` minus its `upstreamCount` to its `riverIndex`. In Figure 2,
river 5's watershed is rivers 3 through 5, river 6's is rivers 0 through 6, river 12's is rivers 9 through 12, and river
13's, at the outlet, is all 14 rivers.

Figure 1 is numbered in a topological order that is not a DFS order: rivers 3 and 4 fall between river 5 and its
upstream rivers 1 and 2. A depth first search from the outlet orders the same network 1, 2, 5, 3, 4, 6, 7, 8, 9.
`examples/migrate_v2_to_v3.py` sorts a network into DFS order.

## Muskingum Routing

The Muskingum equation relates the outflow $Q_{t+1}$ to the inflow at the next step $I_{t+1}$, inflow at the current step $I_{t}$,
and the current discharge $Q_t$. Equivalent forms also use the notations $Q_{t}$ and $Q_{t-1}$. For a primer on the
derivation of the Muskingum equation as a relationship of storage, inflow, and outflow, try the HEC-HMS manual pages on the
[Muskingum Model](https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/channel-flow/muskingum-model) and the
[Muskingum-Cunge Model](https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/channel-flow/muskingum-cunge-model).

$$
Q_{t+1} = c_1\, I_{t+1} + c_2\, I_t + c_3\, Q_t
$$

Where the coefficients $c_1$, $c_2$, $c_3$ are given by:

$$
c_1 = \frac{\Delta t / k - 2x}{\Delta t / k + 2(1-x)}
\qquad
c_2 = \frac{\Delta t / k + 2x}{\Delta t / k + 2(1-x)}
\qquad
c_3 = \frac{2(1-x) - \Delta t / k}{\Delta t / k + 2(1-x)}
$$

Note that:

- Mass is conserved.
- $c_1 + c_2 + c_3 = 1$
- The $k$ parameters can be shown to be the flood wave travel time along the channel in seconds
- The $x$ parameter is a dimensionless "attenuation" factor between 0 (max attenuation) and 0.5 (no attenuation).
- Every time step depends on the step before. The solution must be found sequentially rather than parallelized across time steps.

The Muskingum Cunge equation adds the term $c_4$ to weight adding a lateral inflow term $Q_l$ to
each segment. In the RAPID assumption, lateral flow is the runoff volume divided by the runoff
time step, meaning all runoff enters the channel and exits the basin in the interval it is generated.
No overland flow time or attenuation occurs.

$$
Q_{t+1} = c_1\, I_{t+1} + c_2\, I_t + c_3\, Q_t + c_4\, Q_{l,t}
$$

where $c_4 = c_1 + c_2$.

## Derivation of Matrix Muskingum

### Adjacency matrix

An adjacency matrix $A$ encodes the connectivity of the river network into a square matrix of size $n_\text{segments} \times n_\text{segments}$.
The entry $A_{ij}$ is 1 if segment $j$ flows into segment $i$, and 0 otherwise. Because of the topological sorting, $A$ is strictly lower
triangular (all zeros on the diagonal and above). In each row, the 1 indicates that the river in that column is directly upstream. Conversely,
in each column, a 1 indicates that the river is directly downstream. Every column should have exactly 1 nonzero entry except for outlets which
have no values. Rows with no values are headwater segments. $A^T$ is also common and has an inverse interpretation of rows and columns.

<div class="matrix-grid">
<table>
  <caption><em>Table 1: Adjacency matrix for the river network in Figure 1. Entry A<sub>ij</sub> = 1 indicates segment j flows into segment i.</em></caption>
  <tr><th></th><th>R1</th><th>R2</th><th>R3</th><th>R4</th><th>R5</th><th>R6</th><th>R7</th><th>R8</th><th>R9</th></tr>
  <tr><th>R1</th><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr>
  <tr><th>R2</th><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr>
  <tr><th>R3</th><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr>
  <tr><th>R4</th><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr>
  <tr><th>R5</th><td>1</td><td>1</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr>
  <tr><th>R6</th><td></td><td></td><td>1</td><td>1</td><td></td><td></td><td></td><td></td><td></td></tr>
  <tr><th>R7</th><td></td><td></td><td></td><td></td><td>1</td><td>1</td><td></td><td></td><td></td></tr>
  <tr><th>R8</th><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr>
  <tr><th>R9</th><td></td><td></td><td></td><td></td><td></td><td></td><td>1</td><td>1</td><td></td></tr>
</table>
</div>

When you matrix-multiply $A$ (shape $(n_\text{segments}, n_\text{segments})$) by a column vector of discharges $Q_t$ (shape $(n_\text{segments}, 1)$),
the result is a column vector of inflows $I_t$ with shape $(n_\text{segments}, 1)$. The zero entries of $A$ drop the rows of $Q_t$ that are not directly
upstream. The upstream segments are summed for the total inflow.

$$
I_t = A\, Q_t
$$

Because the rivers are listed in a topological order, $A$ is strictly lower triangular and the diagonal is zero. $A$ is extremely sparse so computationally
it's more efficient to store it in a sparse format and do math only on the non-zero elements.

### Matrix Muskingum

Using the adjacency matrix and the definition that $I_t$ is the sum of upstream discharges, we can replace all $I$ with $A\, Q$.

Muskingum equation:

$$
Q_{t+1} = c_1\, I_{t+1} + c_2\, I_t + c_3\, Q_t
$$

Substitute all $I$ terms with $A\, Q$:

$$
Q_{t+1} = c_1\, \bigl(A\, Q_{t+1}\bigr) + c_2\, \bigl(A\, Q_t\bigr) + c_3\, Q_t
$$

Move all $Q_{t+1}$ terms to the left-hand side:

$$
Q_{t+1} - c_1\, \bigl(A\, Q_{t+1}\bigr) = c_2\, \bigl(A\, Q_t\bigr) + c_3\, Q_t
$$

Factor out $Q_{t+1}$ by right side matrix multiplication:

$$
\bigl(\mathbf{I} - c_1\, A\bigr)\; Q_{t+1} = c_2\, \bigl(A\, Q_t\bigr) + c_3\, Q_t
$$

With the lateral inflow term of Muskingum Cunge:

$$
\bigl(\mathbf{I} - c_1\, A\bigr)\; Q_{t+1} = c_2\, \bigl(A\, Q_t\bigr) + c_3\, Q_t + c_4\, Q_{l,t}
$$

Notes:

- The LHS is the same for every time step (only depends on $A$ and $c_1$).
- The RHS changes every sequential time step because it depends on the previous discharge.
- $c_2$ and $c_3$ can be expressed as 1D vectors of length $n_\text{segments}$ for vector multiplication.
- $c_1$ is an $NxN$ diagonal matrix so it can be subtracted from $I$ after multiplication with A.

## Unit Hydrograph Lateral Inflow (Planned)

### Derivation

The planned unit-hydrograph routing procedure would combine Muskingum channel routing with unit hydrograph lateral inflow. The unit hydrograph shape is developed
in a way that accounts for all the attenuation and travel time during the overland flow process. Thus, we cannot directly
add it to the equation using the $c_4$ term as in the Muskingum Cunge equation because additional attenuation and travel time
will be applied. Instead, a unique method for solving uses the superposition principle where the unit hydrograph convolution
discharge is superimposed on the routed discharge so the signal of the runoff transformation is preserved. In this form,
$Q_l$ represents the discharge generated from a unit hydrograph convolution during the current timestep, $t$.

Establish the relationship between total, routed, and lateral flow.

$$
Q_{\text{total},t+1} = Q_{\text{channel},t+1} + Q_{\text{lateral},t+1}
$$

$$
I_{\text{total},t+1} = I_{\text{channel},t+1} + I_{\text{lateral},t+1}
$$

Annotate the Matrix Muskingum equation with a subscript "ch" or "total".

$$
Q_{\text{channel},t+1} = c_1\, I_{\text{total},t+1} + c_2\, I_{\text{total},t} + c_3\, Q_{\text{channel},t}
$$

Substitute the relationship $I_t = A\, Q_t$

$$
Q_{\text{channel},t+1} = c_1\, \bigl(A\, Q_{\text{channel},t+1} + A\, Q_{\text{lateral},t+1}\bigr) + c_2\, \bigl(A\, Q_{\text{total},t}\bigr) + c_3\, Q_{\text{channel},t}
$$

Move all $Q_{\text{channel},t+1}$ terms to the left-hand side:

$$
Q_{\text{channel},t+1} - c_1\, \bigl(A\, Q_{\text{channel},t+1}\bigr) = c_1\, \bigl(A\, Q_{\text{lateral},t+1}\bigr) + c_2\, \bigl(A\, Q_{\text{total},t}\bigr) + c_3\, Q_{\text{channel},t}
$$

Factor out $Q_{\text{channel},t+1}$ by right side matrix multiplication:

$$
\bigl(\mathbf{I} - c_1\, A\bigr)\; Q_{\text{channel},t+1} = c_1\, \bigl(A\, Q_{\text{lateral},t+1}\bigr) + c_2\, \bigl(A\, Q_{\text{total},t}\bigr) + c_3\, Q_{\text{channel},t}
$$

Calculate the superposition of the channel and lateral flow:

$$
Q_{\text{total},t+1} = Q_{\text{channel},t+1} + Q_{\text{lateral},t+1}
$$

### Reduced Inner System

The headwater segments have no upstream dependencies and their discharge would be the unit hydrograph convolution output. In the planned 
procedure they could be excluded from the matrix solve and their outflow would enter the system as a known right-hand-side contribution 
during the superposition step.

### Kernel Structure

A unit hydrograph kernel has shape $(n_\text{steps},\; n_\text{basins})$. It is a 2D array of each basin's unit hydrograph, discretized 
to the routing time step, and concatenated into 1 array.

- Each column is the discretized unit hydrograph for one basin assuming a unit runoff depth ($R = 1\,\text{m}$).
- Each row value is the average flow ($\text{m}^2/\text{s}$) over the corresponding time step.
- Volume conservation requires: $\displaystyle\sum_{i} K_{i,j} \cdot \Delta t = A_j$ where $A_j$ is the basin area ($\text{m}^2$).

### Unit Hydrograph Convolution

Given a timeseries of runoff depths $r_t$ (meters per time step), the lateral inflow at time $t$ is found by convolving the runoff with the kernel:

$$
Q_{l,t} = \sum_{\tau=0}^{n_\text{steps}-1} K_\tau \cdot r_{t-\tau}
$$

The planned implementation would compute this as a fourier transform over the full timeseries using `scipy.signal.fftconvolve`.
That convolution helper is not part of the current release and is not wired into routing.

## Forward Substitution Algorithm

Because of the careful preparation of the routing matrices, the linear system can be solved for directly by forward substitution, 
one river at a time from upstream to downstream, without iterative methods, preconditioning, or factorization. This is a significant simplification and efficiency gain in terms of the 
complexity of the algorithm as well as how efficiently it can be implemented and compiled in code.

For a unit lower triangular system $L\, x = b$ of size $n$:

$$
x_i = b_i - \sum_{j < i} L_{ij}\, x_j \qquad \text{for } i = 1, 2, \ldots, n
$$

Because $L_{ii} = 1$, no division is needed. Each unknown $x_i$ depends only on previously
solved values $x_1, \ldots, x_{i-1}$, so the system is solved sequentially from the first
row to the last. `river-route` does not store the matrix $L$ at all, and it does not solve it one time step at a
time. Once every river upstream of a river is routed, the whole right-hand side of that river is known at every step,
so v3 solves the system one river at a time: each river's whole time series is routed before the next river, in
topological order, and pushed onto its single downstream river through a `downstream_indices` vector.

```
for j = 1, 2, ..., n:
    x[j] = b[j]                          # diagonal is 1, so x[j] = b[j] directly
    for each row i where L[i,j] != 0:    # only the nonzero entries below the diagonal
        b[i] -= L[i,j] * x[j]            # subtract the now-known contribution
```

*Listing 1: Conceptual column-oriented forward substitution pseudocode for a generic unit lower triangular system (not river-route's storage layout).*

In the v3 kernels (`route_job` in `river_route/router/_numba_kernels.py` and the `route_river` routing methods next
to it) there is no matrix. A river's inflow row already holds the summed discharge of its upstream rivers at every
routing step when the river is reached, so its whole series is the recurrence
$Q_{t+1} = c_1\, I_{t+1} + c_2\, I_t + c_3\, Q_t + c_4\, Q_{l,t}$ with every $I$ known, and once it is solved the
series is added into the inflow row of the river downstream:

```python
for i in range(n_rivers):                        # topological order: upstream before downstream
    inflow = inflow_rows[i]                      # summed discharge of the rivers upstream of i, at every level
    series[0] = q[i]
    for t in range(n_steps):
        q[i] = c1[i] * inflow[t + 1] + c2[i] * inflow[t] + c3[i] * q[i] + c4[i] * runoff_rate[i, t]
        series[t + 1] = q[i]
    if downstream_indices[i] >= 0:               # outlets have no downstream (-1)
        inflow_rows[downstream_indices[i]] += series
```

*Listing 2: River-at-a-time forward substitution over downstream indices, the order in which `route_job` routes.*

- **Time:** $O((n + m)\, T)$ where $n$ is the number of river segments, $m$ is the number of edges
  (upstream-downstream connections), and $T$ is the number of routing steps. For tree-structured river networks,
  $m = n - 1$.
- **Space:** No matrix is stored. Only the inflow rows of the rivers that some but not all of their upstream rivers
  have been routed into are held at once.

This is optimal — every edge is visited exactly once per routing step.

## Numerical stability

### Relationship between dt, k, x

Because the Muskingum equation is a valid solution to a partial differential equation, the equation will conserve mass and route water correctly
regardless of the choice of dt, k, and x. However, the choice of those parameters can cause physically impossible results causing either 1) negative
discharge or 2) oscillation from negative to positive discharge. These conditions happen when either $c_1$ or $c_3$ is negative. The coefficients
are all fractions with the same denominator which will always be positive for positive $dt$ and $k$. We can create inequalities describing when the
numerators are positive comparing dt (a subjective choice) to k and x (physically derived parameters):

this set of equations shows that for each coefficient, there is an inequality, and there is a range that dt must fall between depending on how large x is

$$
\begin{aligned}
c_1 > 0 &\implies \Delta t > 2kx \\
c_2 > 0 &\implies \Delta t > -2kx \\
c_3 > 0 &\implies \Delta t < 2k(1-x)
\end{aligned}
$$

Note that $c_2$ is always positive because $dt$, $k$, and $x$ are all positive.
The remaining 2 inequalities can be combined to find the range of valid $dt$ values:

$$
2kx < \Delta t < 2k(1-x) \\
x = 0 \implies 0 < \Delta t < 2k \\
x = 0.5 \implies k < \Delta t < k \implies \Delta t = k
$$

For a given $k$ and $x$, the valid range of $dt$ values is:

$$
2k(1-x) - 2kx \\
2k(1-2x)
$$

There are several noteworthy insights from these equations:

- When the is no attenuation ($x = 0.5$), the only valid $dt$ is $k$.
- When the is maximum attenuation ($x = 0$), the valid range of $dt$ is from 0 to $2k$.
- The lower bound of valid $dt$ values is 0 when $x = 0$ meaning maximum attenuation such as at a reservoir.
- The upper bound of valid $dt$ values approaches $k$ as $x$ approaches 0.5 meaning no attenuation.

### When k >>> dt

Long river segments have larger $k$ values. It can be hard to pick a large enough $dt$ to satisfy the $c1$ inequality $\Delta t > 2kx$ while
still keeping it small enough for the simulation result to be meaningful. Rivers in this case should be split into a series of smaller rivers 
with smaller k. When doing this, you can divide k proportional to the length of the smaller river segment. The `river-route` code refers to 
this as making **_"substeps"_** down the river.

### When dt >>> k

Short river segments have smaller $k$ values. It can be hard to pick a small enough $dt$ to satisfy the $c3$ inequality $\Delta t < 2k(1-x)$ 
while still keeping it large enough for the simulation to compute in a reasonable amount of time. Many higher density delineated rivers will 
place confluences within a few pixels of each other causing short segments to be made to fill gap between confluences. This happens in real 
rivers also but more often in DEM delineations. 

**Option 1:** You can handle these rivers by decreasing the time step for all rivers or for only the invalid rivers. The `river-route` code 
refers to this as making **_"subcycles"_** or multiple computations within what is normally only a single computation cycle. This is a good 
solution as long as it doesn't make you simulation take an unreasonable amount of additional time. Sometimes the segments are only a few meters 
long which theoretically need a $dt$ to a few seconds. That is burdensome for marginal accuracy gains. You might also consider:

**Option 2:** Editing your delineated streams to force confluences to overlap that are within a few pixels of each other. The error introduced
by this is likely quite small relative to the uncertainty you expect in hydrological reconstructions. This is more intense of an exercise 
because it requires search for rivers and editing the topology information (ID and downstream ID) at minimum but more likely also the geometry.

## References

- HEC-HMS Users Manual introduction to Muskingum Model
  <https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/channel-flow/muskingum-model>
- HEC-HMS Users Manual introduction to Muskingum Cunge Model
  <https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/channel-flow/muskingum-cunge-model>
- HEC-HMS Users Manual introduction to Unit Hydrographs
  <https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/transform/unit-hydrograph-basic-concepts>
- David, C. H. (2011) River Network Routing on the NHDPlus Dataset *Journal of Hydrometeorology* [doi:10.1175/2011JHM1345.1](https://doi.org/10.1175/2011JHM1345.1)
- NRCS (2010). *National Engineering Handbook*, Part 630: Hydrology, Chapter 16: Hydrographs. United States Department of Agriculture.
- Wikipedia: [Triangular matrix — Forward substitution](https://en.wikipedia.org/wiki/Triangular_matrix#Forward_substitution).
- Wikipedia: [Topological sorting](https://en.wikipedia.org/wiki/Topological_sorting).
