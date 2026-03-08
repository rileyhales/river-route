# Mathematical Derivations

This page collects the mathematical foundations behind the routing algorithms and unit hydrograph transformers in river-route.

---

## Muskingum Routing

The Muskingum method models the storage $S$ in a river segment as a linear combination of inflow $I$ and outflow $Q$:

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

## Matrix Formulation

### Network Adjacency

For a river network with $n$ segments whose rows are topologically sorted (upstream before
downstream), define the **adjacency matrix** $A$ where $A_{ij} = 1$ if segment $j$ flows
directly into segment $i$. Because of the topological ordering, $A$ is strictly lower triangular.

The inflow to each segment is the sum of outflow from its upstream neighbors: $I = A\, Q$.
Substituting into the Muskingum recursion and collecting terms yields the matrix equation:

$$
\bigl(\mathbf{I} - c_1\, A\bigr)\; Q_{t+1} = c_2\, \bigl(A\, Q_t\bigr) + c_3\, Q_t
$$

where $\mathbf{I}$ is the identity matrix and $c_1$, $c_2$, $c_3$ are diagonal matrices of
per-segment coefficients.

The left-hand-side matrix $L = \mathbf{I} - c_1\, A$ is **unit lower triangular** (ones on the
diagonal, nonzero entries only below the diagonal) because $A$ is strictly lower triangular.
This structure allows the system to be solved by
[forward substitution](forward-substitution.md) without any factorization.

### Muskingum Cunge Routing and the RAPID assumption

The RapidMuskingum adds a lateral inflow term $Q_l$ directly to each segment. The lateral flow is
the runoff volume divided by the runoff time step, meaning all runoff enters the channel in the
interval it is generated (no overland flow delay).

$$
Q_{t+1} = c_1\, I_{t+1} + c_2\, I_t + c_3\, Q_t + c_4\, Q_{l,t}
$$

where $c_4 = c_1 + c_2$. In matrix form:

$$
\bigl(\mathbf{I} - c_1\, A\bigr)\; Q_{t+1} = c_2\, \bigl(A\, Q_t\bigr) + c_3\, Q_t + c_4\, Q_{l,t}
$$

### UnitMuskingum — Unit Hydrograph Lateral Inflow

UnitMuskingum combines Muskingum channel routing with unit hydrograph lateral inflow. Discharge
at each segment is the superposition of channel-routed flow and overland flow transformed
through a unit hydrograph convolution. Let $Q_l$ denote the [UH convolution](#unit-hydrograph-convolution).

The channel-routed component $Q_\text{ch}$ is solved separately from $Q_l$:

$$
\bigl(\mathbf{I} - c_1\, A\bigr)\; Q_{\text{ch},t+1}
= c_1\, \bigl(A\, Q_{l,t+1}\bigr)
+ c_2\, \bigl(A\, Q_{\text{full},t}\bigr)
+ c_3\, Q_{\text{ch},t}
$$

$$
Q_{\text{full},t+1} = Q_{\text{ch},t+1} + Q_{l,t+1}
$$

Headwater segments (no upstream connections) have no channel inflow, so their discharge is entirely 
the UH convolution output. These segments are excluded from the matrix solve, and their outflow enters 
the reduced system as a known right-hand-side contribution. This typically cuts the sparse linear 
solve size roughly in half.

---

## Unit Hydrograph Theory

### Kernel Structure

A unit hydrograph kernel has shape $(n_\text{steps},\; n_\text{basins})$:

- Each column is the discretized unit hydrograph for one basin assuming a unit runoff depth ($R = 1\,\text{m}$).
- Each row value is the average flow ($\text{m}^2/\text{s}$) over the corresponding time step.
- Volume conservation requires: $\displaystyle\sum_{i} K_{i,j} \cdot \Delta t = A_j$ where $A_j$ is the basin area ($\text{m}^2$).

### Unit Hydrograph Convolution

Given a timeseries of runoff depths $r_t$ (meters per time step), the lateral inflow at time $t$
is obtained by convolving the runoff with the kernel:

$$
Q_{l,t} = \sum_{\tau=0}^{n_\text{steps}-1} K_\tau \cdot r_{t-\tau}
$$

river-route implements this as an FFT-accelerated convolution over the full timeseries using
`scipy.signal.fftconvolve`, with carryover state maintained between successive calls for
warm-starting.

---

## SCS Unit Hydrograph Parameterization

The SCS (NRCS) dimensionless unit hydrograph method generates kernels from basin properties.
Both the triangular and curvilinear variants share the same core parameterization.

**Source:** NRCS National Engineering Handbook (NEH) Part 630, Chapter 16, Table 16-1.

### Common Parameters

Given time-of-concentration $t_c$ (seconds) and runoff duration $t_r$ (seconds):

$$
t_l = 0.6\, t_c \qquad \text{(lag time)}
$$

$$
t_p = t_l + \frac{t_r}{2} \qquad \text{(time to peak)}
$$

### Triangular UH

$$
t_b = 2.67\, t_p \qquad \text{(base time)}
$$

$$
q_p = \frac{2\, A}{t_b} \qquad \text{(peak flow,}\; \text{m}^2/\text{s}\text{)}
$$

The dimensionless shape is a triangle: rising linearly from $(0, 0)$ to $(t_p,\, q_p)$, then
falling linearly to $(t_b,\, 0)$. The kernel is obtained by analytically integrating this
piecewise-linear shape over each $t_r$ interval.

### Curvilinear UH

$$
t_b = 5\, t_p \qquad \text{(base time, table extends to } t/t_p = 5.0\text{)}
$$

$$
q_p = \frac{A}{I_\text{dim} \cdot t_p} \qquad \text{(peak flow;}\; I_\text{dim} \approx 1.333 \text{ from table integral)}
$$

The dimensionless shape follows the standard NRCS tabulated curve (NEH 630, Ch. 16, Table 16-1)
with 33 ordinates from $t/t_p = 0$ to $t/t_p = 5$. The kernel is built by interpolating the
cumulative dimensionless curve and differencing over each $t_r$ interval.

### Volume Conservation

Both methods conserve volume exactly:

$$
\sum_i K_{i,j} \cdot t_r = A_j
$$
