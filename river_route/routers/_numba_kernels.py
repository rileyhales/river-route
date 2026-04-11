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
    rhs = np.empty(n, dtype=np.float32)
    interval_sum = np.empty(n, dtype=np.float32)
    inv_nrpo = np.float32(1.0 / num_routing_per_output)

    for output_step in range(num_output_steps):
        for i in range(n):
            interval_sum[i] = 0.0

        for _ in range(num_routing_per_output):
            # RHS = c3 * Q(t)
            for i in range(n):
                rhs[i] = c3[i] * q_t[i]

            # Fused SpMV and forward solve: (I - c1*A) @ Q(t+1) = c2*(A @ Q(t)) + c3*Q(t)
            # SpMV scatter: rhs += c2 * (A @ Q(t))
            # Forward solve: Q(t+1) = rhs after back-substitution of (I - c1*A)
            for col in range(n):
                q_old = q_t[col]
                q_t[col] = rhs[col]
                for j in range(csc_indptr[col], csc_indptr[col + 1]):
                    row = csc_indices[j]
                    rhs[row] += c2[row] * q_old - lhs_off_data[j] * q_t[col]

            # Accumulate substep results for interval averaging
            for i in range(n):
                interval_sum[i] += q_t[i]

        # Average over substeps, clamp negative values
        for i in range(n):
            val = interval_sum[i] * inv_nrpo
            discharge_array[output_step, i] = val if val > 0.0 else 0.0


@numba.njit(cache=True, fastmath=True)
def rapid_route(
        csc_indptr, csc_indices, lhs_off_data,
        c2, c3, c4_dt, q_t, qlateral,
        discharge_array,
        num_substeps,
):
    """Full RapidMuskingum routing loop with lateral inflow."""
    n = len(q_t)
    rhs = np.empty(n, dtype=np.float32)
    interval_sum = np.empty(n, dtype=np.float32)
    inv_substeps = np.float32(1.0 / num_substeps)
    num_runoff_steps = qlateral.shape[0]

    for t in range(num_runoff_steps):
        for i in range(n):
            interval_sum[i] = 0.0

        for _ in range(num_substeps):
            # RHS = c3 * Q(t) + c4/dt * qlateral(t)
            for i in range(n):
                rhs[i] = c3[i] * q_t[i] + c4_dt[i] * qlateral[t, i]

            # Fused SpMV and forward solve: (I - c1*A) @ Q(t+1) = c2*(A @ Q(t)) + c3*Q(t) + c4/dt*ql(t)
            # SpMV scatter: rhs += c2 * (A @ Q(t))
            # Forward solve: Q(t+1) = rhs after back-substitution of (I - c1*A)
            for col in range(n):
                q_old = q_t[col]
                q_t[col] = rhs[col]
                for j in range(csc_indptr[col], csc_indptr[col + 1]):
                    row = csc_indices[j]
                    rhs[row] += c2[row] * q_old - lhs_off_data[j] * q_t[col]

            # Accumulate substep results for interval averaging
            for i in range(n):
                interval_sum[i] += q_t[i]

        # Average over substeps, clamp negative values
        for i in range(n):
            val = interval_sum[i] * inv_substeps
            discharge_array[t, i] = val if val > 0.0 else 0.0


# noinspection PyPep8Naming
@numba.njit(cache=True, fastmath=True)
def unit_route(
        lhs_indptr, lhs_indices, lhs_off_data,
        a_inner_indptr, a_inner_indices, a_inner_data,
        a_hw_indptr, a_hw_indices, a_hw_data,
        c1_inner, c2_inner, c3_inner,
        hw_idx, inner_idx,
        q_ch, q_full,
        convolved_lateral,
        discharge_array,
        num_substeps,
):
    """Full UnitMuskingum routing loop with unit hydrograph lateral inflow."""
    n_inner = len(inner_idx)
    n_hw = len(hw_idx)
    num_runoff_steps = convolved_lateral.shape[0]
    inv_substeps = np.float32(1.0 / num_substeps)

    rhs = np.empty(n_inner, dtype=np.float32)
    interval_sum = np.empty(n_inner, dtype=np.float32)
    ql_hw = np.empty(n_hw, dtype=np.float32)
    ql_inner = np.empty(n_inner, dtype=np.float32)
    a_inner_result = np.empty(n_inner, dtype=np.float32)
    a_hw_result = np.empty(n_inner, dtype=np.float32)
    c1_A_ql = np.empty(n_inner, dtype=np.float32)

    for t in range(num_runoff_steps):
        # Extract headwater and inner lateral inflows from convolved runoff
        for i in range(n_hw):
            ql_hw[i] = convolved_lateral[t, hw_idx[i]]
        for i in range(n_inner):
            ql_inner[i] = convolved_lateral[t, inner_idx[i]]

        # Headwater segments: Q = ql (no channel routing, UH convolution only)
        for i in range(n_hw):
            discharge_array[t, hw_idx[i]] = ql_hw[i]

        # SpMV: A_inner @ ql_inner (inner-to-inner lateral flow contributions)
        for i in range(n_inner):
            a_inner_result[i] = 0.0
        for col in range(n_inner):
            val = ql_inner[col]
            for j in range(a_inner_indptr[col], a_inner_indptr[col + 1]):
                a_inner_result[a_inner_indices[j]] += a_inner_data[j] * val

        # SpMV: A_hw_to_inner @ ql_hw (headwater-to-inner lateral flow contributions)
        for i in range(n_inner):
            a_hw_result[i] = 0.0
        for col in range(n_hw):
            val = ql_hw[col]
            for j in range(a_hw_indptr[col], a_hw_indptr[col + 1]):
                a_hw_result[a_hw_indices[j]] += a_hw_data[j] * val

        # c1 * A @ ql = c1 * (A_inner @ ql_inner + A_hw_to_inner @ ql_hw)
        for i in range(n_inner):
            c1_A_ql[i] = c1_inner[i] * (a_inner_result[i] + a_hw_result[i])

        for i in range(n_inner):
            interval_sum[i] = 0.0

        for _ in range(num_substeps):
            # RHS = c1*A*ql + c2*(A_hw_to_inner @ ql_hw) + c3*Q_ch(t)
            for i in range(n_inner):
                rhs[i] = c1_A_ql[i] + c2_inner[i] * a_hw_result[i] + c3_inner[i] * q_ch[i]

            # Fused SpMV and forward solve: (I - c1*A) @ Q_ch(t+1) = RHS + c2*(A_inner @ Q_full(t))
            # SpMV scatter: rhs += c2 * (A_inner @ Q_full(t))
            # Forward solve: Q_ch(t+1) = rhs after back-substitution of (I - c1*A)
            for col in range(n_inner):
                q_old = q_full[col]
                q_ch[col] = rhs[col]
                for j in range(lhs_indptr[col], lhs_indptr[col + 1]):
                    row = lhs_indices[j]
                    rhs[row] += c2_inner[row] * q_old - lhs_off_data[j] * q_ch[col]

            # Q_full = Q_ch + ql, accumulate for interval averaging
            for i in range(n_inner):
                q_full[i] = q_ch[i] + ql_inner[i]
                interval_sum[i] += q_full[i]

        # Average over substeps, clamp negative values
        for i in range(n_inner):
            val = interval_sum[i] * inv_substeps
            discharge_array[t, inner_idx[i]] = val if val > 0.0 else 0.0
