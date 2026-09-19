"""Tests for the Network class: parsing, the cached routing schedule, stability analysis, and subdivision."""

import numpy as np
import pandas as pd
import pytest
from conftest import build_network

import river_route as rr

DT = 3600


def configs_for(tmp_path, df: pd.DataFrame, **kwargs) -> rr.Configs:
    """Write a parameter table and return a Configs pointing at it, since every class takes only a Configs."""
    params_file = tmp_path / f'params_{abs(hash(df.to_csv()))}.parquet'
    df.to_parquet(params_file, index=False)
    return rr.Configs(params_file=str(params_file), log=False, **kwargs)


def chain(tmp_path, k, x=0.2) -> rr.Network:
    """A chain of len(k) rivers, each draining into the next, with the given per river k."""
    k = np.asarray(k, dtype=float)
    n = k.shape[0]
    df = pd.DataFrame(
        {'river_id': np.arange(1, n + 1), 'next_river_id': list(range(2, n + 1)) + [-1], 'k': k, 'x': np.full(n, x)}
    )
    return rr.Network.from_configs(configs_for(tmp_path, df))


def window(net: rr.Network, dt: float = DT) -> tuple[np.ndarray, np.ndarray]:
    """The travel time range a single sub-reach of each river must land in to be stable for dt."""
    x = net.x.astype(np.float64)
    return dt / (2 * (1 - x)), np.where(x > 0, dt / (2 * x), np.inf)


# ── parsing and topology ────────────────────────────────────────────────────


def test_network_reads_a_params_file(tmp_path):
    synthetic = build_network(tmp_path / 'net', n_rivers=5)
    net = rr.Network.from_configs(rr.Configs(params_file=str(synthetic.params_file), log=False))
    assert len(net) == 5
    np.testing.assert_array_equal(net.river_ids, np.arange(1, 6))
    np.testing.assert_array_equal(net.downstream_indices, [1, 2, 3, 4, -1])
    assert net.alpha is None and net.beta is None


def test_network_rejects_a_broken_table(tmp_path):
    base = {'river_id': [1, 2], 'next_river_id': [2, -1], 'k': [1.0, 1.0], 'x': [0.2, 0.2]}
    with pytest.raises(ValueError, match='missing required column'):
        rr.Network.from_configs(configs_for(tmp_path, pd.DataFrame(base).drop(columns='k')))
    with pytest.raises(ValueError, match='duplicate river IDs'):
        rr.Network.from_configs(configs_for(tmp_path, pd.DataFrame({**base, 'river_id': [1, 1]})))
    with pytest.raises(ValueError, match='not in the river_id column'):
        rr.Network.from_configs(configs_for(tmp_path, pd.DataFrame({**base, 'next_river_id': [9, -1]})))
    with pytest.raises(ValueError, match='topologically sorted'):  # river 1 drains into river 2, listed above it
        rr.Network.from_configs(
            configs_for(tmp_path, pd.DataFrame({**base, 'river_id': [2, 1], 'next_river_id': [-1, 2]}))
        )
    with pytest.raises(ValueError, match='alpha column'):
        alpha = pd.DataFrame({**base, 'alpha': [0.0, 1.0], 'beta': [1.0, 1.0]})
        rr.Network.from_configs(configs_for(tmp_path, alpha, coeff='dynamic'))


def test_network_builds_the_same_thing_directly_and_from_configs(tmp_path):
    """from_configs only reads the values off a Configs; the direct constructor takes those same values."""
    df = pd.DataFrame({'river_id': [1, 2], 'next_river_id': [2, -1], 'k': [3600.0, 3600.0], 'x': [0.2, 0.2]})
    configs = configs_for(tmp_path, df)
    direct = rr.Network(configs.params_file, var_river_id='river_id', coeff='static')
    from_configs = rr.Network.from_configs(configs)
    np.testing.assert_array_equal(direct.river_ids, from_configs.river_ids)
    np.testing.assert_array_equal(direct.k, from_configs.k)
    np.testing.assert_array_equal(direct.downstream_indices, from_configs.downstream_indices)
    with pytest.raises(TypeError, match='from_configs takes a Configs'):
        rr.Network.from_configs(str(configs.params_file))
    with pytest.raises(ValueError, match='params_file is required'):
        rr.Network(None)


def test_routing_schedule_is_cached_per_thread_count(tmp_path):
    net = rr.Network.from_configs(
        rr.Configs(params_file=str(build_network(tmp_path / 'net', n_rivers=8).params_file), log=False)
    )
    assert net.routing_schedule(threads=1) is net.routing_schedule(threads=1)
    # without a pool to run regions on, the schedule is the single whole-network job whatever threads says
    jobs, cut_target = net.routing_schedule(threads=4, concurrent=False)
    assert len(jobs) == 1 and cut_target.shape == (0,)
    assert jobs[0][0][0] == 0 and jobs[0][1][0] == len(net)


# ── stability ───────────────────────────────────────────────────────────────


def test_stability_report_counts_each_failure_direction(tmp_path):
    # k=3600 stable (2kx=1440 <= 3600 <= 2k(1-x)=5760); k=36000 too long; k=900 too short
    report = chain(tmp_path, [3600.0, 36000.0, 900.0]).stability_report(DT)
    assert (report.n_rivers, report.n_stable, report.n_too_long, report.n_too_short) == (3, 1, 1, 1)
    assert report.n_resolved_by_split == 1  # the too-long river splits; the too-short one cannot
    assert report.n_unresolvable == 1
    assert report.subreaches_to_create == report.total_subreaches - report.n_rivers
    assert 'too short for dt' in str(report)


def test_stability_reports_add_up_across_networks(tmp_path):
    """evaluate_network.py sweeps many parameter files; reports add so the sweep needs no manual tallying."""
    a = chain(tmp_path, [3600.0, 36000.0]).stability_report(DT)
    b = chain(tmp_path, [900.0]).stability_report(DT)
    total = sum([a, b])
    assert total.n_rivers == 3
    assert total.n_too_long == 1 and total.n_too_short == 1
    assert total.max_subreaches == max(a.max_subreaches, b.max_subreaches)
    with pytest.raises(ValueError, match='different dt'):
        a + chain(tmp_path, [3600.0]).stability_report(1800)


def test_check_stability_actions(tmp_path):
    net = chain(tmp_path, [3600.0, 36000.0])
    assert net.check_stability(DT, action='ignore') is None
    assert chain(tmp_path, [3600.0]).check_stability(DT) is None  # nothing unstable, nothing reported
    with pytest.raises(ValueError, match='not Muskingum-stable'):
        net.check_stability(DT, action='raise')
    assert net.check_stability(DT, action='warn').n_too_long == 1


def test_substeps_required_only_exceeds_one_for_too_short_rivers(tmp_path):
    net = chain(tmp_path, [3600.0, 36000.0, 900.0])
    # river 3: 2*k*(1-x) = 1440, so dt must shrink by ceil(3600/1440) = 3
    np.testing.assert_array_equal(net.substeps_required(DT), [1, 1, 3])


# ── subdivision ─────────────────────────────────────────────────────────────


@pytest.mark.parametrize('mode', ['uniform', 'nonuniform'])
def test_stabilize_conserves_travel_time_and_lateral_inflow(tmp_path, mode):
    net = chain(tmp_path, [3600.0, 36000.0, 25000.0, 900.0])
    sub = net.stabilize(DT, mode=mode)
    starts = sub.reach_indptr[:-1]
    # each river's pieces must sum back to its own k, and its lateral shares to exactly one river's worth
    np.testing.assert_allclose(np.add.reduceat(sub.k.astype(np.float64), starts), net.k, rtol=1e-5)
    np.testing.assert_allclose(np.add.reduceat(sub.lateral_scale.astype(np.float64), starts), 1.0, rtol=1e-5)


@pytest.mark.parametrize('mode', ['uniform', 'nonuniform'])
def test_stabilize_makes_every_piece_of_a_fixed_river_stable(tmp_path, mode):
    net = chain(tmp_path, [3600.0, 36000.0, 25000.0, 20000.0, 900.0])
    sub = net.stabilize(DT, mode=mode)
    lo, hi = window(net)
    inside = (sub.k >= np.repeat(lo, sub.subdivisions) * (1 - 1e-5)) & (
        sub.k <= np.repeat(hi, sub.subdivisions) * (1 + 1e-5)
    )
    per_river = np.logical_and.reduceat(inside, sub.reach_indptr[:-1])
    assert per_river[sub.resolvable].all()


def test_stabilize_does_not_add_a_reach_to_a_river_that_divides_the_bound_exactly(tmp_path):
    """k and x are float32, so a river whose k is an exact multiple of the upper bound must not round up."""
    net = chain(tmp_path, [36000.0])  # x=0.2 -> largest stable piece is 9000, and 36000 / 9000 is exactly 4
    sub = net.stabilize(DT, mode='uniform')
    assert sub.subdivisions.tolist() == [4]
    np.testing.assert_allclose(sub.k, 9000.0, rtol=1e-5)


def test_uniform_and_nonuniform_agree_on_count_but_not_on_boundaries(tmp_path):
    """Feasibility depends only on the piece count, so the modes fix the same rivers with the same reach count."""
    net = chain(tmp_path, [25000.0])
    uniform = net.stabilize(DT, mode='uniform')
    nonuniform = net.stabilize(DT, mode='nonuniform')
    np.testing.assert_array_equal(uniform.subdivisions, nonuniform.subdivisions)
    np.testing.assert_allclose(uniform.k, 25000 / 3, rtol=1e-5)
    # nonuniform packs the largest stable piece (9000) first and leaves the remainder last
    np.testing.assert_allclose(nonuniform.k, [9000.0, 9000.0, 7000.0], rtol=1e-5)


def test_nonuniform_evens_out_a_river_whose_remainder_would_be_too_short(tmp_path):
    """k=20000 packs as 9000 + 9000 + 2000, but 2000 is below the 2250 lower bound, so it falls back to uniform."""
    sub = chain(tmp_path, [20000.0]).stabilize(DT, mode='nonuniform')
    np.testing.assert_allclose(sub.k, 20000 / 3, rtol=1e-5)


def test_too_short_rivers_are_kept_whole_and_flagged(tmp_path):
    net = chain(tmp_path, [3600.0, 900.0])
    sub = net.stabilize(DT)
    assert sub.subdivisions.tolist() == [1, 1]
    assert sub.resolvable.tolist() == [True, False]
    assert sub.too_short.tolist() == [False, True]
    np.testing.assert_allclose(sub.k[1], 900.0, rtol=1e-5)  # kept at its own k, not split into unstable pieces
    assert sub.substeps_required[1] == 3


def test_stabilized_layout_is_a_topologically_sorted_dag(tmp_path):
    net = chain(tmp_path, [36000.0, 3600.0, 25000.0])
    sub = net.stabilize(DT)
    # every reach feeds a later one or nothing, which is what lets a kernel sweep in array order
    assert np.all((sub.downstream_index < 0) | (sub.downstream_index > np.arange(sub.n_reaches)))
    # the last reach of each river is its outlet, is numbered 0, and drains into the next river's first reach
    np.testing.assert_array_equal(sub.outlet_index, sub.reach_indptr[1:] - 1)
    np.testing.assert_array_equal(sub.subreach_number[sub.outlet_index], 0)
    np.testing.assert_array_equal(sub.downstream_index[sub.outlet_index[:-1]], sub.reach_indptr[1:-1])
    assert sub.downstream_index[sub.outlet_index[-1]] == -1
    np.testing.assert_array_equal(sub.reach_river_id[sub.outlet_index], net.river_ids)


def test_stabilized_state_round_trips_per_river(tmp_path):
    sub = chain(tmp_path, [36000.0, 3600.0]).stabilize(DT)
    state = np.array([1.5, 2.5], dtype=np.float32)
    np.testing.assert_array_equal(sub.collapse_state(sub.broadcast_state(state)), state)
    with pytest.raises(ValueError, match='for 2 rivers'):
        sub.broadcast_state(np.zeros(3, dtype=np.float32))


def test_stabilize_with_explicit_weights_apportions_k_by_length(tmp_path):
    """Explicit weights are how sub-reaches are matched to real segment geometry rather than to a count."""
    net = chain(tmp_path, [36000.0, 3600.0])
    sub = net.stabilize(DT, weights=[[3, 1], [1]])
    assert sub.mode == 'weights'
    # the geometry is built as asked even though the 27000 piece is far over the 9000 bound; that river is
    # reported as unstable rather than collapsed, which is how the weights path differs from the automatic modes
    assert sub.subdivisions.tolist() == [2, 1]
    np.testing.assert_allclose(sub.k, [27000.0, 9000.0, 3600.0], rtol=1e-5)
    assert sub.resolvable.tolist() == [False, True]


def test_stabilize_rejects_bad_arguments(tmp_path):
    net = chain(tmp_path, [3600.0])
    with pytest.raises(ValueError, match='dt must be positive'):
        net.stabilize(0)
    with pytest.raises(ValueError, match='mode must be'):
        net.stabilize(DT, mode='linear')
    with pytest.raises(ValueError, match='one sequence per river'):
        net.stabilize(DT, weights=[[1], [1]])
    with pytest.raises(ValueError, match='strictly positive'):
        net.stabilize(DT, weights=[[1, -1]])


# ── how Router uses it ──────────────────────────────────────────────────────


def test_router_builds_and_reuses_one_network(tmp_path):
    synthetic = build_network(tmp_path / 'net', n_rivers=5)
    configs = rr.Configs(params_file=str(synthetic.params_file), discharge_files=[synthetic.path('q.nc')], log=False)
    router = rr.Router(configs)
    assert router.network is router.network  # built once, then reused
    np.testing.assert_array_equal(router.network.river_ids, np.arange(1, 6))


def test_router_accepts_a_prebuilt_network(tmp_path):
    """One parsed and partitioned network can back many simulations, which is what a calibration loop needs."""
    synthetic = build_network(tmp_path / 'net', n_rivers=5)
    configs = rr.Configs(params_file=str(synthetic.params_file), discharge_files=[synthetic.path('q.nc')], log=False)
    network = rr.Network.from_configs(configs)
    router = rr.Router(configs)
    router.network = network
    assert router.network is network


def test_router_takes_a_network_at_construction(tmp_path):
    """The same network given up front, which is what a loop over many Routers does."""
    synthetic = build_network(tmp_path / 'net', n_rivers=5)
    configs = rr.Configs(params_file=str(synthetic.params_file), discharge_files=[synthetic.path('q.nc')], log=False)
    network = rr.Network.from_configs(configs)
    assert rr.Router(configs, network=network).network is network
    with pytest.raises(TypeError, match='must be a Network'):
        rr.Router(configs, network='not a network')
