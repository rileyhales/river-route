import numba
import numpy as np

__all__ = ['static_channel', 'static_vlateral', 'static_vlateral_expanded', 'dynamic_vlateral']

# Routing kernels cover a continuous range of indices, which could be the whole network, so that code does not
# change between the single-threaded and multithreaded cases.

# To facilitate multithreading, the kernels use a "pass" contract to define what each thread is allowed to do.
# Each pass is one call to the kernel, and the pass contract is:
#
#   unthreaded   -- one block covering the network, outlet = -1, cut_target empty. Identical to sweeping the
#                   whole array, so there is no separate serial code path to keep in sync.
#   region pass  -- one block, the contiguous range of a subtree in DFS computation order (see
#                   network.streams.assign_regions). Upstream-closed, so it needs nothing from outside its range and
#                   can run concurrently with the other regions. ``outlet`` is the single river whose
#                   downstream lies outside the block; rather than pushing into another thread's rhs, its
#                   contribution is buffered into boundary[region] for every routing step.
#   main stem    -- the blocks left between the regions, swept in increasing index order after every region
#                   pass has finished. cut_target says which river each buffered contribution drains into, and
#                   it is injected before the blocks are swept for that step.
#
# Coupling between rivers only ever runs downstream, so a region's whole time series can be computed before the
# main stem starts. That is why this needs one barrier for the simulation, not one per time step.


@numba.njit(cache=True, fastmath=True, nogil=True)
def static_channel(
    *,
    q_t,  # Array (n_rivers,) of discharge at current time step, updated in-place
    discharge_array,  # Array (n_steps, n_rivers) to write discharge time series into
    downstream_indices,  # Array (n_rivers,) of downstream river indices, -1 for no downstream
    downstream_c1,  # Array (n_rivers,) of c1 downstream of river at index i
    downstream_c2,  # Array (n_rivers,) of c2 downstream of river at index i
    c3,  # Array (n_rivers,) of c3 for river at index i, used in forward substitution sweep
    n_rivers,  # integer number of rivers in the network
    n_steps,  # integer number of time steps to route
    n_substeps,  # integer number of routing substeps per runoff value
    block_starts,  # Array (n_blocks,) of first river index of each block this call routes
    block_stops,  # Array (n_blocks,) of one past the last river index of each block
    outlet,  # index of the river whose push is buffered instead of applied, -1 when there is none
    region,  # row of boundary this call writes its outlet contribution to
    cut_target,  # Array (n_cuts,) of rivers to inject buffered contributions into; EMPTY for a region pass
    boundary,  # Array (n_regions, n_steps * n_substeps) of outlet contributions handed between passes
):
    """Route the rivers in the given blocks for the whole simulation; see the module header for the pass contract."""
    rhs = np.empty(n_rivers, dtype=np.float32)
    n_blocks = block_starts.shape[0]
    n_cuts = cut_target.shape[0]
    step = 0

    if n_substeps == 1:
        for t in range(n_steps):
            for b in range(n_blocks):
                for i in range(block_starts[b], block_stops[b]):
                    rhs[i] = c3[i] * q_t[i]

            for c in range(n_cuts):
                target = cut_target[c]
                if target >= 0:
                    rhs[target] += boundary[c, step]

            for b in range(n_blocks):
                for i in range(block_starts[b], block_stops[b]):
                    q_old = q_t[i]
                    q_new = rhs[i]
                    q_t[i] = q_new
                    # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
                    discharge_array[t, i] = q_new if q_new > 0.0 else np.float32(0.0)
                    downstream_idx = downstream_indices[i]
                    if downstream_idx >= 0:
                        push = downstream_c2[i] * q_old + downstream_c1[i] * q_new
                        if i == outlet:
                            boundary[region, step] = push
                        else:
                            rhs[downstream_idx] += push
            step += 1
        return

    interval_sum = np.empty(n_rivers, dtype=np.float32)
    inv_substeps = np.float32(1.0 / n_substeps)

    for t in range(n_steps):
        for b in range(n_blocks):
            for i in range(block_starts[b], block_stops[b]):
                interval_sum[i] = 0.0

        for _ in range(n_substeps):
            for b in range(n_blocks):
                for i in range(block_starts[b], block_stops[b]):
                    rhs[i] = c3[i] * q_t[i]

            for c in range(n_cuts):
                target = cut_target[c]
                if target >= 0:
                    rhs[target] += boundary[c, step]

            for b in range(n_blocks):
                for i in range(block_starts[b], block_stops[b]):
                    q_old = q_t[i]
                    q_new = rhs[i]
                    q_t[i] = q_new
                    interval_sum[i] += q_new
                    downstream_idx = downstream_indices[i]
                    if downstream_idx >= 0:
                        push = downstream_c2[i] * q_old + downstream_c1[i] * q_new
                        if i == outlet:
                            boundary[region, step] = push
                        else:
                            rhs[downstream_idx] += push
            step += 1

        for b in range(n_blocks):
            for i in range(block_starts[b], block_stops[b]):
                val = interval_sum[i] * inv_substeps
                # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
                discharge_array[t, i] = val if val > 0.0 else np.float32(0.0)
    return


@numba.njit(cache=True, fastmath=True, nogil=True)
def static_vlateral(
    *,
    q_t,  # Array (n_rivers,) of discharge at current time step, updated in-place
    discharge_array,  # Array (n_steps, n_rivers) to write discharge time series into
    downstream_indices,  # Array (n_rivers,) of downstream river indices, -1 for no downstream
    downstream_c1,  # Array (n_rivers,) of c1 downstream of river at index i
    downstream_c2,  # Array (n_rivers,) of c2 downstream of river at index i
    c3,  # Array (n_rivers,) of c3 for river at index i, the order of the solution pass
    n_rivers,  # integer number of rivers in the network
    n_steps,  # integer number of time steps to route
    n_substeps,  # integer number of routing substeps per runoff value
    vlateral,  # Array (n_steps, n_rivers) of lateral inflow time series for each river
    c4_dt,  # Array (n_rivers,) of c4 * dt_routing for river at index i
    block_starts,  # Array (n_blocks,) of first river index of each block this call routes
    block_stops,  # Array (n_blocks,) of one past the last river index of each block
    outlet,  # index of the river whose push is buffered instead of applied, -1 when there is none
    region,  # row of boundary this call writes its outlet contribution to
    cut_target,  # Array (n_cuts,) of rivers to inject buffered contributions into; EMPTY for a region pass
    boundary,  # Array (n_regions, n_steps * n_substeps) of outlet contributions handed between passes
):
    """Route the rivers in the given blocks for the whole simulation; see the module header for the pass contract."""
    rhs = np.empty(n_rivers, dtype=np.float32)
    n_blocks = block_starts.shape[0]
    n_cuts = cut_target.shape[0]
    step = 0

    if n_substeps == 1:
        for t in range(n_steps):
            for b in range(n_blocks):
                for i in range(block_starts[b], block_stops[b]):
                    rhs[i] = c3[i] * q_t[i] + c4_dt[i] * vlateral[t, i]

            for c in range(n_cuts):
                target = cut_target[c]
                if target >= 0:
                    rhs[target] += boundary[c, step]

            for b in range(n_blocks):
                for i in range(block_starts[b], block_stops[b]):
                    q_old = q_t[i]
                    q_new = rhs[i]
                    q_t[i] = q_new
                    # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
                    discharge_array[t, i] = q_new if q_new > 0.0 else np.float32(0.0)
                    downstream_idx = downstream_indices[i]
                    if downstream_idx >= 0:
                        push = downstream_c2[i] * q_old + downstream_c1[i] * q_new
                        if i == outlet:
                            boundary[region, step] = push
                        else:
                            rhs[downstream_idx] += push
            step += 1
        return

    interval_sum = np.empty(n_rivers, dtype=np.float32)
    q_ext_t = np.empty(n_rivers, dtype=np.float32)
    inv_substeps = np.float32(1.0 / n_substeps)

    for t in range(n_steps):
        for b in range(n_blocks):
            for i in range(block_starts[b], block_stops[b]):
                interval_sum[i] = 0.0
                q_ext_t[i] = c4_dt[i] * vlateral[t, i]

        for _ in range(n_substeps):
            for b in range(n_blocks):
                for i in range(block_starts[b], block_stops[b]):
                    rhs[i] = c3[i] * q_t[i] + q_ext_t[i]

            for c in range(n_cuts):
                target = cut_target[c]
                if target >= 0:
                    rhs[target] += boundary[c, step]

            for b in range(n_blocks):
                for i in range(block_starts[b], block_stops[b]):
                    q_old = q_t[i]
                    q_new = rhs[i]
                    q_t[i] = q_new
                    interval_sum[i] += q_new
                    downstream_idx = downstream_indices[i]
                    if downstream_idx >= 0:
                        push = downstream_c2[i] * q_old + downstream_c1[i] * q_new
                        if i == outlet:
                            boundary[region, step] = push
                        else:
                            rhs[downstream_idx] += push
            step += 1

        for b in range(n_blocks):
            for i in range(block_starts[b], block_stops[b]):
                val = interval_sum[i] * inv_substeps
                # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
                discharge_array[t, i] = val if val > 0.0 else np.float32(0.0)
    return


@numba.njit(cache=True, fastmath=True)
def static_vlateral_expanded(
    *,
    q,  # Array (n_reaches,) of per-reach discharge state, updated in-place (instantaneous end-of-step value)
    substeps_per_reach,  # Array (n_reaches,) of temporal substeps to route+average each reach (>= 1)
    discharge_array,  # Array (n_steps, n_rivers) to write discharge time series into (per original river)
    parent_index,  # Array (n_reaches,) of the original river index a reach belongs to (output + vlateral)
    downstream_index,  # Array (n_reaches,) of downstream reach indices, -1 for no downstream
    downstream_c1,  # Array (n_reaches,) of c1 of reach r's downstream reach, pre-gathered for sequential push
    downstream_c2,  # Array (n_reaches,) of c2 of reach r's downstream reach, pre-gathered for sequential push
    c3,  # Array (n_reaches,) of c3 for reach r, used in the forward substitution sweep
    c4_dt,  # Array (n_reaches,) of (c4 / dt_runoff) * lateral_scale for reach r: lateral VOLUME -> rate forcing
    vlateral,  # Array (n_steps, n_rivers) of lateral inflow time series for each original river
    n_reaches,  # integer number of expanded reaches in the network
    n_steps,  # integer number of time steps to route
):
    """
    Unified stabilized Muskingum with lateral inflow over an EXPANDED network (see network.streams.expand_network). Both
    stability levers are handled by one code path:

        - Subdivision (too-long reaches) is already materialized: a split river is a contiguous chain of reaches,
          each with c1/c2/c3 built from k/N. The chain routes through the forward-substitution sweep exactly like
          any other reaches, so subdivided-reach state needs no special handling. Such reaches usually have
          substeps_per_reach == 1, but a river needing BOTH levers carries its substep count on every sub-reach.
        - Substepping (too-short reaches) is done per reach: a reach with substeps_per_reach == S routes S internal
          Muskingum iterations at dt = period/S against a HELD period inflow forcing, and reports the average of the
          S substep outflows. substeps_per_reach == 1 is the degenerate case whose average is the single
          instantaneous value, reducing exactly to static_muskingum_vlateral.

    Reach state ``q`` is the instantaneous end-of-step outflow (this is what couples to downstream, identical to the
    original kernel). The REPORTED discharge is the per-reach substep average. discharge_array is written per
    original river via parent_index; because each river's reaches are contiguous and its outlet is processed last,
    the last write for a river is its outlet's value (instantaneous for subdivided rivers, averaged for substepped).

    Caller responsibilities (kernel does no validation to stay a pure hot path; see network.streams.expand_network):
        - c1, c2, c3, c4_dt are per reach and built for dt = period / substeps_per_reach[r].
        - c4_dt already folds in the reach's lateral_scale (1 / subdivisions), so lateral volume is conserved.
        - q has length n_reaches, seeded by broadcasting each river's initial discharge across its reaches.
        - substeps_per_reach[r] >= 1.

    Note: a substepped reach reports its substep MEAN but couples its instantaneous endpoint to downstream (held
    across the downstream's substeps). Total routed volume is conserved; only the per-period shape differs.
    """
    rhs = np.empty(n_reaches, dtype=np.float32)

    for t in range(n_steps):
        # seed each reach's inflow forcing with its lateral inflow (held constant across its substeps)
        for r in range(n_reaches):
            rhs[r] = c4_dt[r] * vlateral[t, parent_index[r]]

        for r in range(n_reaches):
            s = substeps_per_reach[r]
            c3r = c3[r]
            q_old = q[r]
            forcing = rhs[r]  # lateral + upstream inflow terms pushed in earlier this sweep; held for all substeps
            if s == 1:
                q_new = forcing + c3r * q_old
                reported = q_new
            else:
                acc = np.float32(0.0)
                q_s = q_old
                for _ in range(s):
                    q_s = forcing + c3r * q_s
                    acc += q_s
                q_new = q_s
                reported = acc / np.float32(s)
            q[r] = q_new
            # outlet is processed last, so its (clamped) value wins per river
            # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
            discharge_array[t, parent_index[r]] = reported if reported > 0.0 else np.float32(0.0)
            downstream_idx = downstream_index[r]
            if downstream_idx >= 0:
                rhs[downstream_idx] += downstream_c2[r] * q_old + downstream_c1[r] * q_new
    return


@numba.njit(cache=True, fastmath=True, nogil=True)
def dynamic_vlateral(
    *,
    q_t,  # Array (n_rivers,) of discharge at current time step, updated in-place
    discharge_array,  # Array (n_steps, n_rivers) to write discharge time series into
    downstream_indices,  # Array (n_rivers,) of downstream river indices, -1 for no downstream
    alpha,  # Array (n_rivers,) of alpha for river at index i, used to compute K_i
    beta,  # Array (n_rivers,) of beta for river at index i, used to compute K_i
    x,  # Array (n_rivers,) of x for river at index i, used to compute Muskingum coefficients
    dt_routing,  # integer routing timestep in seconds, used to compute Muskingum coefficients
    dt_runoff,  # integer runoff timestep in seconds, used to scale vlateral to Q per routing substep
    n_rivers,  # integer number of rivers in the network
    n_steps,  # integer number of time steps to route
    n_substeps,  # integer number of routing substeps per runoff value
    vlateral,  # Array (n_steps, n_rivers) of lateral volume time series for each river
    block_starts,  # Array (n_blocks,) of first river index of each block this call routes
    block_stops,  # Array (n_blocks,) of one past the last river index of each block
    outlet,  # index of the river whose state is buffered instead of pushed, -1 when there is none
    region,  # row of boundary_old/boundary_new this call writes its outlet state to
    cut_target,  # Array (n_cuts,) of rivers to inject buffered contributions into; EMPTY for a region pass
    boundary_old,  # Array (n_regions, n_steps * n_substeps) of each region outlet's pre-step discharge
    boundary_new,  # Array (n_regions, n_steps * n_substeps) of each region outlet's post-step discharge
):
    """
    Nonlinear Muskingum with lateral inflow. Each substep, K_i is recomputed
    as K_i = alpha_i * max(q_t[i], qmin)^beta_i and the Muskingum
    coefficients are rebuilt from that K before the forward-substitution sweep.

    Routes the rivers in the given blocks; see the module header for the pass contract. This kernel buffers its
    outlet's q_old and q_new rather than the finished contribution, because a push is weighted by the
    coefficients of the river being pushed INTO, and those are rebuilt from that river's own discharge every
    substep. A region cannot know them: they belong to the main stem, which has not run yet. Handing over the
    two discharges instead lets the main stem apply its own freshly computed c1/c2, which is exactly what an
    unpartitioned sweep does.
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
    n_blocks = block_starts.shape[0]
    n_cuts = cut_target.shape[0]
    step = 0

    for t in range(n_steps):
        for b in range(n_blocks):
            for i in range(block_starts[b], block_stops[b]):
                interval_sum[i] = 0.0
                q_ext_t[i] = vlateral[t, i] * inv_dt_runoff

        for _ in range(n_substeps):
            for b in range(n_blocks):
                for i in range(block_starts[b], block_stops[b]):
                    q = q_t[i] if q_t[i] > qmin else qmin
                    k_i = alpha[i] * q ** beta[i]
                    dt_div_k = dt_routing / k_i
                    two_x = 2.0 * x[i]
                    two_one_minus_x = 2.0 * (1.0 - x[i])
                    denom = dt_div_k + two_one_minus_x
                    c1[i] = (dt_div_k - two_x) / denom
                    c2[i] = (dt_div_k + two_x) / denom
                    c3[i] = (two_one_minus_x - dt_div_k) / denom

            for b in range(n_blocks):
                for i in range(block_starts[b], block_stops[b]):
                    rhs[i] = c3[i] * q_t[i] + (c1[i] + c2[i]) * q_ext_t[i]

            # inject what the region passes produced, weighted by this pass's own coefficients for the target
            for c in range(n_cuts):
                target = cut_target[c]
                if target >= 0:
                    rhs[target] += c2[target] * boundary_old[c, step] + c1[target] * boundary_new[c, step]

            for b in range(n_blocks):
                for i in range(block_starts[b], block_stops[b]):
                    q_old = q_t[i]
                    q_new = rhs[i]
                    q_t[i] = q_new
                    interval_sum[i] += q_new
                    downstream_idx = downstream_indices[i]
                    if downstream_idx >= 0:
                        if i == outlet:
                            boundary_old[region, step] = q_old
                            boundary_new[region, step] = q_new
                        else:
                            rhs[downstream_idx] += c2[downstream_idx] * q_old + c1[downstream_idx] * q_new
            step += 1

        for b in range(n_blocks):
            for i in range(block_starts[b], block_stops[b]):
                val = interval_sum[i] * inv_substeps
                # todo clamping negative discharge to zero is a stopgap; fix the root-cause instability
                discharge_array[t, i] = val if val > 0.0 else np.float32(0.0)
    return
