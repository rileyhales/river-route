"""Numba-accelerated kernels for Muskingum-family routing recurrences."""
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
            # Start each river i with the term that depends only on its own
            # previous discharge: rhs_i = c3_i * Q_i(t).
            for i in range(n):
                rhs[i] = c3[i] * q_t[i]

            # Move upstream to downstream. When river u=col is reached, rhs[u]
            # already contains every upstream contribution, so it is Q_u(t+dt).
            # Then add u's effect to each downstream river d:
            # rhs_d += c2_d * Q_u(t) + c1_d * Q_u(t+dt).
            for col in range(n):
                q_old = q_t[col]
                rhs_col = rhs[col]
                q_t[col] = rhs_col
                interval_sum[col] += rhs_col
                for j in range(csc_indptr[col], csc_indptr[col + 1]):
                    row = csc_indices[j]
                    rhs[row] += c2[row] * q_old - lhs_off_data[j] * rhs_col

        # Output is the nonnegative mean over routing substeps in this
        # discharge interval.
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
    """
    (I - c1*A) @ Q(t+1) = c2*(A @ Q(t)) + c3*Q(t) + c4/dt*ql(t)
    Full RapidMuskingum routing loop with lateral inflow.
    """
    n_rivers = len(q_t)
    n_runoff_steps = qlateral.shape[0]
    rhs = np.empty(n_rivers, dtype=np.float32)

    if num_substeps == 1:
        for t in range(n_runoff_steps):
            # solve the lateral contribution
            # rhs = c3 * Q(t) + (c4 / dt_runoff) * qlateral(t).
            for i in range(n_rivers):
                rhs[i] = c3[i] * q_t[i] + c4_dt[i] * qlateral[t, i]
            # for each river (col in discharge array), add its contribution to each downstream river (row):
            # rhs += c2 * Q_up(t) + c1 * Q_up(t+dt).
            for col in range(n_rivers):
                q_old = q_t[col]
                rhs_col = rhs[col]
                q_t[col] = rhs_col
                discharge_array[t, col] = rhs_col if rhs_col > 0.0 else 0.0
                # for each downstream segment (row), add the contribution from this upstream segment (col)
                # the range(indptr[col], indptr[col+1]) is length 1 if there is a downstream, 0 if outlet
                for j in range(csc_indptr[col], csc_indptr[col + 1]):
                    row = csc_indices[j]
                    rhs[row] += c2[row] * q_old - lhs_off_data[j] * rhs_col
        return

    interval_sum = np.empty(n_rivers, dtype=np.float32)
    inv_substeps = np.float32(1.0 / num_substeps)

    for t in range(n_runoff_steps):
        for i in range(n_rivers):
            interval_sum[i] = 0.0

        for _ in range(num_substeps):
            # Start each river i with local carryover plus lateral inflow:
            # rhs_i = c3_i * Q_i(t) + (c4_i / dt_runoff) * qlateral_i(t).
            for i in range(n_rivers):
                rhs[i] = c3[i] * q_t[i] + c4_dt[i] * qlateral[t, i]

            # Move upstream to downstream. Once Q_u(t+dt) is known, add u's
            # old and new discharge contributions to each downstream river d:
            # rhs_d += c2_d * Q_u(t) + c1_d * Q_u(t+dt).
            for col in range(n_rivers):
                q_old = q_t[col]
                rhs_col = rhs[col]
                q_t[col] = rhs_col
                interval_sum[col] += rhs_col
                for j in range(csc_indptr[col], csc_indptr[col + 1]):
                    row = csc_indices[j]
                    rhs[row] += c2[row] * q_old - lhs_off_data[j] * rhs_col

        # Output is the nonnegative mean over routing substeps in this runoff
        # interval.
        for i in range(n_rivers):
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

    rhs = np.empty(n_inner, dtype=np.float32)
    ql_hw = np.empty(n_hw, dtype=np.float32)
    ql_inner = np.empty(n_inner, dtype=np.float32)
    a_inner_result = np.empty(n_inner, dtype=np.float32)
    a_hw_result = np.empty(n_inner, dtype=np.float32)
    c1_A_ql = np.empty(n_inner, dtype=np.float32)

    if num_substeps == 1:
        # One routing step per runoff interval: Q_out(t) = max(Q_full(t+dt), 0).
        for t in range(num_runoff_steps):
            # Current unit-hydrograph lateral response split by river type.
            for i in range(n_hw):
                ql_hw[i] = convolved_lateral[t, hw_idx[i]]
            for i in range(n_inner):
                ql_inner[i] = convolved_lateral[t, inner_idx[i]]

            # Headwaters have no upstream channel contribution in this reduced
            # system, so their discharge is their lateral response.
            for i in range(n_hw):
                discharge_array[t, hw_idx[i]] = ql_hw[i]

            # For each inner river i, collect lateral water from inner upstream
            # rivers: sum_u A_inner[i, u] * ql_inner[u].
            for i in range(n_inner):
                a_inner_result[i] = 0.0
            for col in range(n_inner):
                val = ql_inner[col]
                for j in range(a_inner_indptr[col], a_inner_indptr[col + 1]):
                    a_inner_result[a_inner_indices[j]] += a_inner_data[j] * val

            # For each inner river i, collect lateral water from headwater
            # upstream rivers: sum_h A_hw[i, h] * ql_hw[h].
            for i in range(n_inner):
                a_hw_result[i] = 0.0
            for col in range(n_hw):
                val = ql_hw[col]
                for j in range(a_hw_indptr[col], a_hw_indptr[col + 1]):
                    a_hw_result[a_hw_indices[j]] += a_hw_data[j] * val

            # New-step lateral contribution for inner river i:
            # c1_i * sum_all_upstream_lateral_i(t).
            for i in range(n_inner):
                c1_A_ql[i] = c1_inner[i] * (a_inner_result[i] + a_hw_result[i])

            # Start each inner river i with terms known before inner upstream
            # channel states are added:
            # c1_i * upstream_lateral(t) + c2_i * headwater_lateral(t)
            # + c3_i * Q_ch_i(t).
            for i in range(n_inner):
                rhs[i] = c1_A_ql[i] + c2_inner[i] * a_hw_result[i] + c3_inner[i] * q_ch[i]

            # Move inner rivers upstream to downstream. Once Q_ch,u(t+dt) is
            # known, total discharge is Q_full,u(t+dt) = Q_ch,u(t+dt) + ql_u(t).
            # Add u's channel contribution to each downstream inner river d:
            # rhs_d += c2_d * Q_full,u(t) + c1_d * Q_ch,u(t+dt).
            for col in range(n_inner):
                q_old = q_full[col]
                rhs_col = rhs[col]
                q_ch[col] = rhs_col
                q_full_col = rhs_col + ql_inner[col]
                q_full[col] = q_full_col
                discharge_array[t, inner_idx[col]] = q_full_col if q_full_col > 0.0 else 0.0
                for j in range(lhs_indptr[col], lhs_indptr[col + 1]):
                    row = lhs_indices[j]
                    rhs[row] += c2_inner[row] * q_old - lhs_off_data[j] * rhs_col
        return

    interval_sum = np.empty(n_inner, dtype=np.float32)
    inv_substeps = np.float32(1.0 / num_substeps)

    for t in range(num_runoff_steps):
        # Current unit-hydrograph lateral response split by river type.
        for i in range(n_hw):
            ql_hw[i] = convolved_lateral[t, hw_idx[i]]
        for i in range(n_inner):
            ql_inner[i] = convolved_lateral[t, inner_idx[i]]

        # Headwaters have no upstream channel contribution in this reduced
        # system, so their discharge is their lateral response.
        for i in range(n_hw):
            discharge_array[t, hw_idx[i]] = ql_hw[i]

        # For each inner river i, collect lateral water from inner upstream
        # rivers: sum_u A_inner[i, u] * ql_inner[u].
        for i in range(n_inner):
            a_inner_result[i] = 0.0
        for col in range(n_inner):
            val = ql_inner[col]
            for j in range(a_inner_indptr[col], a_inner_indptr[col + 1]):
                a_inner_result[a_inner_indices[j]] += a_inner_data[j] * val

        # For each inner river i, collect lateral water from headwater upstream
        # rivers: sum_h A_hw[i, h] * ql_hw[h].
        for i in range(n_inner):
            a_hw_result[i] = 0.0
        for col in range(n_hw):
            val = ql_hw[col]
            for j in range(a_hw_indptr[col], a_hw_indptr[col + 1]):
                a_hw_result[a_hw_indices[j]] += a_hw_data[j] * val

        # New-step lateral contribution for inner river i:
        # c1_i * sum_all_upstream_lateral_i(t).
        for i in range(n_inner):
            c1_A_ql[i] = c1_inner[i] * (a_inner_result[i] + a_hw_result[i])

        for i in range(n_inner):
            interval_sum[i] = 0.0

        for _ in range(num_substeps):
            # Start each inner river i with terms known before inner upstream
            # channel states are added:
            # c1_i * upstream_lateral(t) + c2_i * headwater_lateral(t)
            # + c3_i * Q_ch_i(t).
            for i in range(n_inner):
                rhs[i] = c1_A_ql[i] + c2_inner[i] * a_hw_result[i] + c3_inner[i] * q_ch[i]

            # Move inner rivers upstream to downstream. Once Q_ch,u(t+dt) is
            # known, total discharge is Q_full,u(t+dt) = Q_ch,u(t+dt) + ql_u(t).
            # Add u's channel contribution to each downstream inner river d:
            # rhs_d += c2_d * Q_full,u(t) + c1_d * Q_ch,u(t+dt).
            for col in range(n_inner):
                q_old = q_full[col]
                rhs_col = rhs[col]
                q_ch[col] = rhs_col
                q_full_col = rhs_col + ql_inner[col]
                q_full[col] = q_full_col
                interval_sum[col] += q_full_col
                for j in range(lhs_indptr[col], lhs_indptr[col + 1]):
                    row = lhs_indices[j]
                    rhs[row] += c2_inner[row] * q_old - lhs_off_data[j] * rhs_col

        # Output is the nonnegative mean of total discharge over routing
        # substeps in this runoff interval.
        for i in range(n_inner):
            val = interval_sum[i] * inv_substeps
            discharge_array[t, inner_idx[i]] = val if val > 0.0 else 0.0
