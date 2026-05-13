import numba
import numpy as np

__all__ = [
    'static_muskingum',
    'static_muskingum_vlateral',
    'static_muskingum_qexternal',
    'dynamic_muskingum_vlateral',
]


@numba.njit(cache=True, fastmath=True)
def static_muskingum(
        *,
        q_t,  # Array shape (n_rivers,) of discharge at current time step, updated in-place
        discharge_array,  # Array shape (n_steps, n_rivers) to write discharge time series into
        downstream_indices,  # Array shape (n_rivers,) of downstream river indices, -1 for no downstream
        downstream_c1,  # Array shape (n_rivers,) of c1 downstream of river at index i
        downstream_c2,  # Array shape (n_rivers,) of c2 downstream of river at index i
        c3,  # Array shape (n_rivers,) of c3 for river at index i, used in forward substitution sweep
        n_rivers,  # integer number of rivers in the network
        n_steps,  # integer number of time steps to route
        n_substeps,  # integer number of routing substeps per runoff value
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
    return


# todo it ought to be faster to solve qlateral * c4dt. if not, update all naming and verbiage to vlateral because its
#  not really q and that's confusing.
@numba.njit(cache=True, fastmath=True)
def static_muskingum_vlateral(
        *,
        q_t,  # Array shape (n_rivers,) of discharge at current time step, updated in-place
        discharge_array,  # Array shape (n_steps, n_rivers) to write discharge time series into
        downstream_indices,  # Array shape (n_rivers,) of downstream river indices, -1 for no downstream
        downstream_c1,  # Array shape (n_rivers,) of c1 downstream of river at index i
        downstream_c2,  # Array shape (n_rivers,) of c2 downstream of river at index i
        c3,  # Array shape (n_rivers,) of c3 for river at index i, the order of the solution pass
        n_rivers,  # integer number of rivers in the network
        n_steps,  # integer number of time steps to route
        n_substeps,  # integer number of routing substeps per runoff value
        qlateral,  # Array shape (n_steps, n_rivers) of lateral inflow time series for each river
        c4_dt,  # Array shape (n_rivers,) of c4 * dt_routing for river at index i
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
    return


@numba.njit(cache=True, fastmath=True)
def static_muskingum_qexternal(
        *,
        q_t,  # Array shape (n_rivers,) of discharge at current time step, updated in-place
        discharge_array,  # Array shape (n_steps, n_rivers) to write discharge time series into
        downstream_indices,  # Array shape (n_rivers,) of downstream river indices, -1 for no downstream
        downstream_c1,  # Array shape (n_rivers,) of c1 downstream of river at index i
        downstream_c2,  # Array shape (n_rivers,) of c2 downstream of river at index i
        c3,  # Array shape (n_rivers,) of c3 for river at index i, the order of the solution pass
        n_rivers,  # integer number of rivers in the network
        n_steps,  # integer number of time steps to route
        n_substeps,  # integer number of routing substeps per runoff value
        qexternal,  # Array shape (n_steps, n_rivers) of external discharge time series for each river
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
                            downstream_c2[i] * q_old
                            + downstream_c1[i] * q_new
                            + (downstream_c1[i] + downstream_c2[i]) * qext_i
                    )
        return

    interval_sum = np.empty(n_rivers, dtype=np.float32)
    q_ext_t = np.empty(n_rivers, dtype=np.float32)
    inv_substeps = np.float32(1.0 / n_substeps)

    for t in range(n_steps):
        for i in range(n_rivers):
            interval_sum[i] = 0.0
            q_ext_t[i] = qexternal[t, i]

        for _ in range(n_substeps):
            for i in range(n_rivers):
                rhs[i] = c3[i] * q_t[i]

            for i in range(n_rivers):
                q_old = q_t[i]
                q_new = rhs[i]
                q_t[i] = q_new
                qext_i = q_ext_t[i]
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
    return


@numba.njit(cache=True, fastmath=True)
def dynamic_muskingum_vlateral(
        *,
        q_t,  # Array shape (n_rivers,) of discharge at current time step, updated in-place
        discharge_array,  # Array shape (n_steps, n_rivers) to write discharge time series into
        downstream_indices,  # Array shape (n_rivers,) of downstream river indices, -1 for no downstream
        alpha,  # Array shape (n_rivers,) of alpha for river at index i, used to compute K_i
        beta,  # Array shape (n_rivers,) of beta for river at index i, used to compute K_i
        x,  # Array shape (n_rivers,) of x for river at index i, used to compute Muskingum coefficients
        dt_routing,  # integer routing timestep in seconds, used to compute Muskingum coefficients
        dt_runoff,  # integer runoff timestep in seconds, used to scale vlateral to Q per routing substep
        n_rivers,  # integer number of rivers in the network
        n_steps,  # integer number of time steps to route
        n_substeps,  # integer number of routing substeps per runoff value
        vlateral,  # Array shape (n_steps, n_rivers) of lateral volume time series for each river
):
    """
    Nonlinear Muskingum with lateral inflow. Each substep, K_i is recomputed
    as K_i = alpha_i * max(q_t[i], qmin)^beta_i and the Muskingum
    coefficients are rebuilt from that K before the forward-substitution sweep.
    """
    rhs = np.empty(n_rivers, dtype=np.float32)
    c1 = np.empty(n_rivers, dtype=np.float32)
    c2 = np.empty(n_rivers, dtype=np.float32)
    c3 = np.empty(n_rivers, dtype=np.float32)
    interval_sum = np.empty(n_rivers, dtype=np.float32)
    q_ext_t = np.empty(n_rivers, dtype=np.float32)

    inv_substeps = np.float32(1.0 / n_substeps)
    inv_dt_runoff = np.float32(1.0 / dt_runoff)
    qmin = np.float32(1e-6)

    for t in range(n_steps):
        for i in range(n_rivers):
            interval_sum[i] = 0.0
            q_ext_t[i] = vlateral[t, i] * inv_dt_runoff

        for _ in range(n_substeps):
            for i in range(n_rivers):
                q = q_t[i] if q_t[i] > qmin else qmin
                k_i = alpha[i] * q ** beta[i]
                dt_div_k = dt_routing / k_i
                two_x = 2.0 * x[i]
                two_one_minus_x = 2.0 * (1.0 - x[i])
                denom = dt_div_k + two_one_minus_x
                c1[i] = (dt_div_k - two_x) / denom
                c2[i] = (dt_div_k + two_x) / denom
                c3[i] = (two_one_minus_x - dt_div_k) / denom

            for i in range(n_rivers):
                rhs[i] = c3[i] * q_t[i] + (c1[i] + c2[i]) * q_ext_t[i]

            for i in range(n_rivers):
                q_old = q_t[i]
                q_new = rhs[i]
                q_t[i] = q_new
                interval_sum[i] += q_new
                downstream_idx = downstream_indices[i]
                if downstream_idx >= 0:
                    rhs[downstream_idx] += c2[downstream_idx] * q_old + c1[downstream_idx] * q_new

        for i in range(n_rivers):
            discharge_array[t, i] = interval_sum[i] * inv_substeps
    return
