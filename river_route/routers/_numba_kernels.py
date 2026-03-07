"""
Numba-accelerated routing kernels using the forward substitution algorithm for unit lower triangular systems.
"""
import numba
import numpy as np


@numba.njit(cache=True, fastmath=True)
def muskingum_route(
        csc_indptr, csc_indices, lhs_off_data,
        c2, c3, q_t,
        discharge_array,
        num_output_steps, num_routing_per_output,
):
    """Full Muskingum channel-only routing loop."""
    n = len(q_t)
    rhs = np.empty(n, dtype=np.float64)
    interval_sum = np.empty(n, dtype=np.float64)
    inv_nrpo = 1.0 / num_routing_per_output

    for output_step in range(num_output_steps):
        for i in range(n):
            interval_sum[i] = 0.0

        for _ in range(num_routing_per_output):
            # RHS = c2 * (A @ q_t) + c3 * q_t
            for i in range(n):
                rhs[i] = c3[i] * q_t[i]
            for col in range(n):
                q_val = q_t[col]
                for j in range(csc_indptr[col], csc_indptr[col + 1]):
                    row = csc_indices[j]
                    rhs[row] += c2[row] * q_val

            # Forward solve: unit lower triangular
            for col in range(n):
                q_t[col] = rhs[col]
                for j in range(csc_indptr[col], csc_indptr[col + 1]):
                    rhs[csc_indices[j]] -= lhs_off_data[j] * q_t[col]

            for i in range(n):
                interval_sum[i] += q_t[i]

        for i in range(n):
            discharge_array[output_step, i] = interval_sum[i] * inv_nrpo


@numba.njit(cache=True, fastmath=True)
def rapid_substeps(
        csc_indptr, csc_indices, lhs_off_data,
        c2, c3, q_t, ql_t,
        rhs, interval_sum,
        num_substeps,
):
    """Inner substep loop for RapidMuskingum with lateral inflow."""
    n = len(q_t)
    for _ in range(num_substeps):
        for i in range(n):
            rhs[i] = c3[i] * q_t[i] + ql_t[i]
        for col in range(n):
            q_val = q_t[col]
            for j in range(csc_indptr[col], csc_indptr[col + 1]):
                row = csc_indices[j]
                rhs[row] += c2[row] * q_val
        for col in range(n):
            q_t[col] = rhs[col]
            for j in range(csc_indptr[col], csc_indptr[col + 1]):
                rhs[csc_indices[j]] -= lhs_off_data[j] * q_t[col]
        for i in range(n):
            interval_sum[i] += q_t[i]


# noinspection PyPep8Naming
@numba.njit(cache=True, fastmath=True)
def unit_substeps(
        csc_indptr, csc_indices, lhs_off_data,
        c2, c3,
        c1_A_ql, hw_contrib, ql_inner,
        q_ch, q_full,
        rhs, interval_sum,
        num_substeps,
):
    """Inner substep loop for UnitMuskingum on the reduced inner system."""
    n = len(q_ch)
    for _ in range(num_substeps):
        # RHS = c1*A*ql + c2*(hw_contrib + A_inner @ q_full) + c3*q_ch
        for i in range(n):
            rhs[i] = c1_A_ql[i] + c2[i] * hw_contrib[i] + c3[i] * q_ch[i]
        for col in range(n):
            q_val = q_full[col]
            for j in range(csc_indptr[col], csc_indptr[col + 1]):
                row = csc_indices[j]
                rhs[row] += c2[row] * q_val

        # Forward solve
        for col in range(n):
            q_ch[col] = rhs[col]
            for j in range(csc_indptr[col], csc_indptr[col + 1]):
                rhs[csc_indices[j]] -= lhs_off_data[j] * q_ch[col]

        # Post-solve
        for i in range(n):
            q_full[i] = q_ch[i] + ql_inner[i]
            interval_sum[i] += q_full[i]
