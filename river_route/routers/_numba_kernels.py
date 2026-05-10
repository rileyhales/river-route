import numba
import numpy as np


@numba.njit(cache=True, fastmath=True)
def linear_muskingum(
        q_t, discharge_array, downstream_indices,
        downstream_c1, downstream_c2, c3,
        n_rivers, n_steps, n_substeps,
):
    rhs = np.empty(n_rivers, dtype=np.float32)

    if n_substeps == 1:
        for t in range(n_steps):
            for i in range(n_rivers):
                rhs[i] = c3[i] * q_t[i]

            for i in range(n_rivers):
                q_old = q_t[i]
                q_new = rhs[i]
                q_t[i] = q_new
                discharge_array[t, i] = q_new
                downstream_idx = downstream_indices[i]
                if downstream_idx >= 0:
                    rhs[downstream_idx] += downstream_c2[i] * q_old + downstream_c1[i] * q_new
        return

    interval_sum = np.empty(n_rivers, dtype=np.float32)
    inv_substeps = np.float32(1.0 / n_substeps)

    for t in range(n_steps):
        for i in range(n_rivers):
            interval_sum[i] = 0.0

        for _ in range(n_substeps):
            for i in range(n_rivers):
                rhs[i] = c3[i] * q_t[i]

            for i in range(n_rivers):
                q_old = q_t[i]
                q_new = rhs[i]
                q_t[i] = q_new
                interval_sum[i] += q_new
                downstream_idx = downstream_indices[i]
                if downstream_idx >= 0:
                    rhs[downstream_idx] += downstream_c2[i] * q_old + downstream_c1[i] * q_new
        for i in range(n_rivers):
            discharge_array[t, i] = interval_sum[i] * inv_substeps


@numba.njit(cache=True, fastmath=True)
def linear_muskingum_qlateral(
        q_t, discharge_array, downstream_indices,
        downstream_c1, downstream_c2, c3,
        n_rivers, n_steps, n_substeps,
        qlateral, c4_dt,
):
    rhs = np.empty(n_rivers, dtype=np.float32)

    if n_substeps == 1:
        for t in range(n_steps):
            for i in range(n_rivers):
                rhs[i] = c3[i] * q_t[i] + c4_dt[i] * qlateral[t, i]

            for i in range(n_rivers):
                q_old = q_t[i]
                q_new = rhs[i]
                q_t[i] = q_new
                discharge_array[t, i] = q_new
                downstream_idx = downstream_indices[i]
                if downstream_idx >= 0:
                    rhs[downstream_idx] += downstream_c2[i] * q_old + downstream_c1[i] * q_new
        return

    interval_sum = np.empty(n_rivers, dtype=np.float32)
    q_ext_t = np.empty(n_rivers, dtype=np.float32)
    inv_substeps = np.float32(1.0 / n_substeps)

    for t in range(n_steps):
        for i in range(n_rivers):
            interval_sum[i] = 0.0
            q_ext_t[i] = c4_dt[i] * qlateral[t, i]

        for _ in range(n_substeps):
            for i in range(n_rivers):
                rhs[i] = c3[i] * q_t[i] + q_ext_t[i]

            for i in range(n_rivers):
                q_old = q_t[i]
                q_new = rhs[i]
                q_t[i] = q_new
                interval_sum[i] += q_new
                downstream_idx = downstream_indices[i]
                if downstream_idx >= 0:
                    rhs[downstream_idx] += downstream_c2[i] * q_old + downstream_c1[i] * q_new
        for i in range(n_rivers):
            discharge_array[t, i] = interval_sum[i] * inv_substeps


@numba.njit(cache=True, fastmath=True)
def linear_muskingum_qexternal(
        q_t, discharge_array, downstream_indices,
        downstream_c1, downstream_c2, c3,
        n_rivers, n_steps, n_substeps,
        qexternal,
):
    rhs = np.empty(n_rivers, dtype=np.float32)

    if n_substeps == 1:
        for t in range(n_steps):
            for i in range(n_rivers):
                rhs[i] = c3[i] * q_t[i]

            for i in range(n_rivers):
                q_old = q_t[i]
                q_new = rhs[i]
                q_t[i] = q_new
                qext_i = qexternal[t, i]
                discharge_array[t, i] = q_new + qext_i
                downstream_idx = downstream_indices[i]
                if downstream_idx >= 0:
                    rhs[downstream_idx] += (
                            downstream_c2[i] * q_old + downstream_c1[i] * q_new +
                            (downstream_c1[i] + downstream_c2[i]) * qext_i
                    )
        return

    interval_sum = np.empty(n_rivers, dtype=np.float32)
    inv_substeps = np.float32(1.0 / n_substeps)

    for t in range(n_steps):
        for i in range(n_rivers):
            interval_sum[i] = 0.0

        for _ in range(n_substeps):
            for i in range(n_rivers):
                rhs[i] = c3[i] * q_t[i]

            for i in range(n_rivers):
                q_old = q_t[i]
                q_new = rhs[i]
                q_t[i] = q_new
                qext_i = qexternal[t, i]
                interval_sum[i] += q_new + qext_i
                downstream_idx = downstream_indices[i]
                if downstream_idx >= 0:
                    rhs[downstream_idx] += (
                            downstream_c2[i] * q_old
                            + downstream_c1[i] * q_new
                            + (downstream_c1[i] + downstream_c2[i]) * qext_i
                    )
        for i in range(n_rivers):
            discharge_array[t, i] = interval_sum[i] * inv_substeps
