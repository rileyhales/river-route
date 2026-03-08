# Forward Substitution for River Routing

This page explains the forward substitution algorithm that river-route uses to solve the Matrix
Muskingum routing equations at each time step. For the full mathematical derivation of the
routing equations, see the [Math Derivations](math.md) page.

---

## Why Forward Substitution?

The Muskingum matrix equation has the form:

$$
L\; Q_{t+1} = b
$$

where $L = \mathbf{I} - c_1\, A$ and $b$ is a known right-hand-side vector assembled from the
previous time step.

Because the river segments are **topologically sorted** (every segment appears after all of its
upstream contributors), the adjacency matrix $A$ is strictly lower triangular. Subtracting it
from the identity produces a **unit lower triangular** matrix — ones on the diagonal, nonzeros
only below. This special structure means the system can be solved in a single forward pass with
no factorization, pivoting, or iteration required.

This is significant because the routing equation must be solved at every time step. Avoiding
matrix factorization entirely makes the per-step cost proportional to the number of edges in the
river network rather than the number of segments cubed.

---

## The Algorithm

For a unit lower triangular system $L\, x = b$ of size $n$:

$$
x_i = b_i - \sum_{j < i} L_{ij}\, x_j \qquad \text{for } i = 1, 2, \ldots, n
$$

Because $L_{ii} = 1$, no division is needed. Each unknown $x_i$ depends only on previously
solved values $x_1, \ldots, x_{i-1}$, so the system is solved sequentially from the first
row to the last.

### Column-Oriented Variant (CSC)

river-route stores the sparse matrix in **compressed sparse column (CSC)** format and uses a
column-oriented forward substitution. Instead of computing one row at a time, it processes
one column at a time: once $x_j$ is known, its contribution is subtracted from all rows below.

```
for j = 1, 2, ..., n:
    x[j] = b[j]                          # diagonal is 1, so x[j] = b[j] directly
    for each row i where L[i,j] != 0:    # only the nonzero entries below the diagonal
        b[i] -= L[i,j] * x[j]            # subtract the now-known contribution
```

This maps directly to the CSC storage where `indptr` and `indices` give fast column-wise access
to nonzero entries. In the implementation (`_numba_kernels.py`), the loop looks like:

```python
for col in range(n):
    q[col] = rhs[col]
    for j in range(csc_indptr[col], csc_indptr[col + 1]):
        rhs[csc_indices[j]] -= lhs_off_data[j] * q[col]
```

### Computational Cost

- **Time:** $O(n + m)$ where $n$ is the number of river segments and $m$ is the number of edges
  (upstream-downstream connections). For tree-structured river networks, $m = n - 1$.
- **Space:** Only the sparse matrix entries are stored. No fill-in occurs because no
  factorization is performed.

This is optimal — every edge is visited exactly once per time step.

---

## Why It Works for River Networks

The key insight connecting linear algebra to hydrology:

1. **River networks are directed acyclic graphs (DAGs).** Water flows downstream; there are no
   cycles.
2. **Every DAG has a topological ordering.** Segments can be numbered so that upstream always
   comes before downstream.
3. **Topological ordering makes the adjacency matrix strictly lower triangular.** If segment $j$
   flows into segment $i$, then $j < i$, so $A_{ij}$ is below the diagonal.
4. **$\mathbf{I} - c_1 A$ is unit lower triangular.** The identity contributes ones on the
   diagonal; $c_1 A$ contributes entries only below.
5. **Unit lower triangular systems are solved by forward substitution.** No pivoting, no
   factorization, no iteration — just a single sweep through the segments in topological order.

The topological sort is performed once when the routing parameters are loaded. The forward
substitution then runs at every routing time step, compiled to native code via Numba JIT.

---

## Scipy Reference Implementations

The production routing kernels use Numba JIT-compiled forward substitution on raw CSC arrays
for performance. The equivalent math expressed with scipy sparse operations is shown below as
a more readable reference. These use `scipy.sparse.linalg.spsolve_triangular` which performs
the same forward substitution internally.

### Muskingum (Channel Only)

```python
from scipy.sparse import diags, eye
from scipy.sparse.linalg import spsolve_triangular

# LHS = I - diag(c1) @ A is unit lower triangular
LHS = eye(n, format='csc') - diags(c1) @ A

for output_step in range(num_output_steps):
    interval_sum = np.zeros(n)

    for _ in range(num_routing_per_output):
        # RHS from previous state: c2 * (A @ q_t) + c3 * q_t
        rhs = c2 * (A @ q_t) + c3 * q_t

        # Solve LHS @ q(t+1) = rhs
        q_t = spsolve_triangular(LHS, rhs, lower=True, unit_diagonal=True)
        interval_sum += q_t

    discharge_array[output_step] = interval_sum / num_routing_per_output
```

### RapidMuskingum (Direct Lateral Inflow)

Adds the lateral inflow term `ql_t` (with `c4 = c1 + c2` already applied by the caller)
directly to the right-hand side.

```python
LHS = eye(n, format='csc') - diags(c1) @ A

for _ in range(num_substeps):
    # RHS adds lateral inflow to the channel routing terms
    rhs = c2 * (A @ q_t) + c3 * q_t + ql_t

    q_t = spsolve_triangular(LHS, rhs, lower=True, unit_diagonal=True)
    interval_sum += q_t
```

### UnitMuskingum (Reduced Inner System)

Headwater segments are excluded from the matrix solve. Their discharge is the UH convolution
output and enters the inner system as a known RHS contribution via `A_hw_to_inner`.

```python
LHS = eye(n_inner, format='csc') - diags(c1_inner) @ A_inner

# Headwater contribution is constant within a runoff interval
hw_contrib = A_hw_to_inner @ ql_hw

for _ in range(num_substeps):
    # RHS combines lateral (c1), upstream full flow (c2), and channel state (c3)
    rhs = (
        c1_inner * (A_inner @ ql_inner + hw_contrib)
        + c2_inner * (A_inner @ q_full + hw_contrib)
        + c3_inner * q_ch
    )

    q_ch = spsolve_triangular(LHS, rhs, lower=True, unit_diagonal=True)
    q_full = q_ch + ql_inner
    interval_sum += q_full
```

---

## References

### Textbooks

- Golub, G. H. & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins
  University Press. — Chapter 3 covers triangular system solution including forward and back
  substitution.

- Trefethen, L. N. & Bau, D. (1997). *Numerical Linear Algebra*. SIAM. — Lecture 17 discusses
  triangular solvers and their stability properties.

### Muskingum Routing

- Cunge, J. A. (1969). On the subject of a flood propagation computation method (Muskingum
  method). *Journal of Hydraulic Research*, 7(2), 205–230.
  [doi:10.1080/00221686909500264](https://doi.org/10.1080/00221686909500264)

- Chow, V. T., Maidment, D. R., & Mays, L. W. (1988). *Applied Hydrology*. McGraw-Hill. —
  Chapter 8 derives the Muskingum method from the storage equation.

### Matrix Formulations for River Network Routing

- David, C. H., Maidment, D. R., Niu, G.-Y., Yang, Z.-L., Habets, F., & Eijkhout, V. (2011).
  River network routing on the NHDPlus dataset. *Journal of Hydrometeorology*, 12(5), 913–934.
  [doi:10.1175/2011JHM1345.1](https://doi.org/10.1175/2011JHM1345.1) — Introduces the matrix
  formulation of Muskingum routing used by RAPID and describes the connection between network
  topology and triangular matrix structure.

### Unit Hydrograph Methods

- NRCS (2010). *National Engineering Handbook*, Part 630: Hydrology, Chapter 16: Hydrographs.
  United States Department of Agriculture. — Source for the SCS dimensionless unit hydrograph
  tables and parameterization equations.

### Accessible Introductions

- Wikipedia: [Triangular matrix — Forward substitution](https://en.wikipedia.org/wiki/Triangular_matrix#Forward_substitution).
  Clear algorithmic description with pseudocode.

- Wikipedia: [Topological sorting](https://en.wikipedia.org/wiki/Topological_sorting).
  Background on why DAGs can always be linearly ordered.
