"""
Routing tests built on small synthetic networks. These run anywhere, with no downloaded test data, and cover
the numerical behavior of the kernels plus every check that stands between a config file and an njit kernel.
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from conftest import SyntheticNetwork, build_network, write_params, write_vlateral

import river_route as rr


def route_vlateral(network: SyntheticNetwork, out_name: str = 'q.nc', **kwargs) -> np.ndarray:
    """Route the network's lateral inflow file and return the discharge array that was written."""
    out = network.path(out_name)
    kwargs.setdefault('params_file', str(network.params_file))
    kwargs.setdefault('vlateral_files', [str(network.vlateral_file)])
    kwargs.setdefault('channel_state_init_file', str(network.state_file))
    kwargs.setdefault('dt_routing', network.dt_runoff)
    rr.Router(forcing='vlateral', discharge_files=[out], log=False, progress_bar=False, **kwargs).route()
    with xr.open_dataset(out) as ds:
        return ds['Q'].values


# ── numerics ────────────────────────────────────────────────────────────────


def test_stable_parameters_conserve_mass(network: SyntheticNetwork):
    """All of the water put into the headwater must leave the outlet when every reach is stable for dt."""
    q = route_vlateral(network)
    outlet_volume = q[:, -1].sum() * network.dt_runoff
    assert outlet_volume == pytest.approx(network.inflow_volume, rel=1e-4)


def test_zero_forcing_gives_zero_discharge(network: SyntheticNetwork):
    zero_file = network.directory / 'vlateral_zero.nc'
    write_vlateral(zero_file, np.zeros((24, network.n_rivers)), np.arange(1, network.n_rivers + 1))
    q = route_vlateral(network, out_name='q_zero.nc', vlateral_files=[str(zero_file)])
    assert np.all(q == 0)


def test_routing_attenuates_downstream(network: SyntheticNetwork):
    """A pulse routed down a chain must arrive later and lower at each successive river."""
    q = route_vlateral(network)
    peaks = q.max(axis=0)
    peak_times = q.argmax(axis=0)
    assert np.all(np.diff(peaks) < 0), f'peak discharge should attenuate downstream, got {peaks}'
    assert np.all(np.diff(peak_times) >= 0), f'peak should arrive later downstream, got {peak_times}'


def test_channel_routing_decays_to_zero(tmp_path):
    """Channel only routing from a uniform initial state drains the network and never goes negative."""
    net = build_network(tmp_path / 'net')
    state_file = net.directory / 'state_10.parquet'
    pd.DataFrame({'Q': np.full(net.n_rivers, 10.0)}).to_parquet(state_file, index=False)
    out = net.path('q_channel.nc')

    rr.Router(
        forcing='channel',
        params_file=str(net.params_file),
        discharge_files=[out],
        channel_state_init_file=str(state_file),
        dt_routing=900,
        dt_total=3600 * 24,
        log=False,
        progress_bar=False,
    ).route()

    with xr.open_dataset(out) as ds:
        q = ds['Q'].values
    assert q.shape == (24 * 4, net.n_rivers)
    assert np.all(q >= 0)
    assert q[-1].sum() < q[0].sum()


def test_final_state_roundtrip(network: SyntheticNetwork):
    """Routing two files in sequence must match saving state after the first and restarting from it."""
    ids = np.arange(1, network.n_rivers + 1)
    volumes = np.zeros((48, network.n_rivers), dtype=np.float32)
    volumes[:3, 0] = 10.0 * network.dt_runoff
    first = network.directory / 'part1.nc'
    second = network.directory / 'part2.nc'
    write_vlateral(first, volumes, ids)
    write_vlateral(second, np.zeros((48, network.n_rivers)), ids)

    both = [network.path('q_both_1.nc'), network.path('q_both_2.nc')]
    rr.Router(
        forcing='vlateral',
        params_file=str(network.params_file),
        vlateral_files=[str(first), str(second)],
        discharge_files=both,
        channel_state_init_file=str(network.state_file),
        log=False,
        progress_bar=False,
    ).route()

    state_between = network.path('state_between.parquet')
    rr.Router(
        forcing='vlateral',
        params_file=str(network.params_file),
        vlateral_files=[str(first)],
        discharge_files=[network.path('q_split_1.nc')],
        channel_state_init_file=str(network.state_file),
        channel_state_final_file=state_between,
        log=False,
        progress_bar=False,
    ).route()
    rr.Router(
        forcing='vlateral',
        params_file=str(network.params_file),
        vlateral_files=[str(second)],
        discharge_files=[network.path('q_split_2.nc')],
        channel_state_init_file=state_between,
        log=False,
        progress_bar=False,
    ).route()

    with xr.open_dataset(both[1]) as ds_all, xr.open_dataset(network.path('q_split_2.nc')) as ds_split:
        np.testing.assert_allclose(ds_all['Q'].values, ds_split['Q'].values, rtol=1e-5, atol=1e-6)


def test_route_twice_is_repeatable(network: SyntheticNetwork):
    """A second route() call must start from the configured initial state, not the previous run's final state."""
    router = rr.Router(
        forcing='vlateral',
        params_file=str(network.params_file),
        vlateral_files=[str(network.vlateral_file)],
        discharge_files=[network.path('q_first.nc')],
        channel_state_init_file=str(network.state_file),
        log=False,
        progress_bar=False,
    )
    router.route()
    with xr.open_dataset(network.path('q_first.nc')) as ds:
        first = ds['Q'].values.copy()

    router.cfg.discharge_files = [network.path('q_second.nc')]
    router.route()
    with xr.open_dataset(network.path('q_second.nc')) as ds:
        second = ds['Q'].values

    np.testing.assert_array_equal(first, second)


# ── time options ────────────────────────────────────────────────────────────


def test_dt_discharge_averages_the_output(network: SyntheticNetwork):
    """Output at a coarser dt_discharge must be the mean of the runoff timestep values it spans."""
    fine = route_vlateral(network, out_name='q_fine.nc')
    coarse = route_vlateral(network, out_name='q_coarse.nc', dt_discharge=network.dt_runoff * 6)
    assert coarse.shape[0] == fine.shape[0] // 6
    expected = fine.reshape(coarse.shape[0], 6, network.n_rivers).mean(axis=1)
    np.testing.assert_allclose(coarse, expected, rtol=1e-5, atol=1e-6)


def test_dt_total_shorter_than_input_truncates(network: SyntheticNetwork):
    """A dt_total covering part of the input file routes that part and writes matching dates and rows."""
    steps = 12
    q = route_vlateral(network, out_name='q_short.nc', dt_total=network.dt_runoff * steps)
    assert q.shape == (steps, network.n_rivers)
    full = route_vlateral(network, out_name='q_full.nc')
    np.testing.assert_allclose(q, full[:steps], rtol=1e-5, atol=1e-6)


def test_dt_total_longer_than_input_raises(network: SyntheticNetwork):
    with pytest.raises(ValueError, match='steps of lateral inflow'):
        route_vlateral(network, out_name='q_long.nc', dt_total=network.dt_runoff * (network.n_steps + 10))


# ── input validation: the checks that keep bad arrays out of the njit kernels ──


def test_vlateral_with_wrong_river_count_raises(network: SyntheticNetwork):
    """A file covering fewer rivers than the network would be read past its end by the kernel."""
    partial = network.directory / 'vlateral_partial.nc'
    write_vlateral(partial, np.zeros((24, 3)), np.arange(1, 4))
    with pytest.raises(ValueError, match='lateral inflow for 3 rivers'):
        route_vlateral(network, out_name='q_partial.nc', vlateral_files=[str(partial)])


def test_vlateral_in_wrong_river_order_raises(network: SyntheticNetwork):
    """Column order is positional in the kernel, so a reordered file would route water down the wrong reach."""
    ids = np.arange(1, network.n_rivers + 1)
    reversed_file = network.directory / 'vlateral_reversed.nc'
    write_vlateral(reversed_file, np.zeros((24, network.n_rivers)), ids[::-1].copy())
    with pytest.raises(ValueError, match='different order'):
        route_vlateral(network, out_name='q_reversed.nc', vlateral_files=[str(reversed_file)])


def test_vlateral_with_different_rivers_raises(network: SyntheticNetwork):
    ids = np.arange(101, 101 + network.n_rivers)
    other = network.directory / 'vlateral_other.nc'
    write_vlateral(other, np.zeros((24, network.n_rivers)), ids)
    with pytest.raises(ValueError, match='not the same rivers'):
        route_vlateral(network, out_name='q_other.nc', vlateral_files=[str(other)])


def test_initial_state_wrong_length_raises(network: SyntheticNetwork):
    """A short state array would be written past its end by the kernel; deep validation reports it first."""
    short_state = network.directory / 'state_short.parquet'
    pd.DataFrame({'Q': np.zeros(network.n_rivers - 2)}).to_parquet(short_state, index=False)
    with pytest.raises(ValueError, match='same number of rows'):
        route_vlateral(network, out_name='q_short_state.nc', channel_state_init_file=str(short_state))


def test_initial_state_wrong_length_raises_without_deep_validation(network: SyntheticNetwork):
    """With the file content checks off, the read must still refuse a state that does not fit the network."""
    short_state = network.directory / 'state_short.parquet'
    pd.DataFrame({'Q': np.zeros(network.n_rivers - 2)}).to_parquet(short_state, index=False)
    with pytest.raises(ValueError, match='state file must have one row per river'):
        route_vlateral(
            network, out_name='q_short_state.nc', channel_state_init_file=str(short_state), deep_validation=False
        )


def test_dispatch_rejects_mismatched_arrays(network: SyntheticNetwork):
    """The dispatch guard is the last line of defense and is checked directly."""
    from river_route.routers._kernel_registry import dispatch

    router = rr.Router(
        forcing='vlateral',
        params_file=str(network.params_file),
        vlateral_files=[str(network.vlateral_file)],
        discharge_files=[network.path('q_unused.nc')],
        dt_routing=network.dt_runoff,
        log=False,
        progress_bar=False,
    )
    router.cfg.validate()
    router._set_vectors_from_params()
    router._set_connectivity_vectors()
    router.num_runoff_steps = 4
    q_t = np.zeros(network.n_rivers - 1, dtype=np.float32)
    with pytest.raises(ValueError, match='q_t has shape'):
        dispatch(
            router,
            q_t=q_t,
            discharge_array=np.zeros((4, network.n_rivers), dtype=np.float32),
            vlateral=np.zeros((4, network.n_rivers), dtype=np.float32),
        )


def test_params_missing_column_raises(tmp_path):
    net = build_network(tmp_path / 'net')
    bad = tmp_path / 'bad_params.parquet'
    pd.read_parquet(net.params_file).drop(columns=['k']).to_parquet(bad, index=False)
    with pytest.raises(ValueError, match='missing k column'):
        rr.Router(
            forcing='vlateral',
            params_file=str(bad),
            vlateral_files=[str(net.vlateral_file)],
            discharge_files=[net.path('q.nc')],
            log=False,
            progress_bar=False,
        ).route()


def test_params_not_topologically_sorted_raises(tmp_path):
    net = build_network(tmp_path / 'net')
    unsorted = tmp_path / 'unsorted.parquet'
    pd.read_parquet(net.params_file).iloc[::-1].to_parquet(unsorted, index=False)
    with pytest.raises(ValueError, match='topologically sorted'):
        rr.Router(
            forcing='vlateral',
            params_file=str(unsorted),
            vlateral_files=[str(net.vlateral_file)],
            discharge_files=[net.path('q.nc')],
            log=False,
            progress_bar=False,
        ).route()


def test_params_x_out_of_range_raises(tmp_path):
    net = build_network(tmp_path / 'net')
    bad = tmp_path / 'bad_x.parquet'
    write_params(bad, n_rivers=net.n_rivers, x=0.9)
    with pytest.raises(ValueError, match='x column must be in the range'):
        rr.Router(
            forcing='vlateral',
            params_file=str(bad),
            vlateral_files=[str(net.vlateral_file)],
            discharge_files=[net.path('q.nc')],
            log=False,
            progress_bar=False,
        ).route()


def test_deep_validation_can_be_skipped(tmp_path):
    """deep_validation=False must skip the file content checks, not the structural ones."""
    net = build_network(tmp_path / 'net')
    bad = tmp_path / 'bad_x.parquet'
    write_params(bad, n_rivers=net.n_rivers, x=0.9)
    rr.Router(
        forcing='vlateral',
        params_file=str(bad),
        vlateral_files=[str(net.vlateral_file)],
        discharge_files=[net.path('q.nc')],
        deep_validation=False,
        unstable_coefficients='ignore',
        log=False,
        progress_bar=False,
    ).route()


# ── stability of the muskingum coefficients ─────────────────────────────────


def route_capturing_logs(net, log_file, **kwargs) -> str:
    """Route with logging directed to a file and return what was written. The router's logger does not
    propagate to the root logger, so caplog cannot see it."""
    rr.Router(
        forcing='vlateral',
        params_file=str(net.params_file),
        vlateral_files=[str(net.vlateral_file)],
        discharge_files=[net.path('q.nc')],
        channel_state_init_file=str(net.state_file),
        log=True,
        log_level='WARNING',
        log_stream=str(log_file),
        progress_bar=False,
        **kwargs,
    ).route()
    return log_file.read_text()


def test_unstable_coefficients_warn_by_default(tmp_path):
    """dt outside 2kx <= dt <= 2k(1-x) oscillates; the default must say so rather than route silently."""
    net = build_network(tmp_path / 'net', k=600.0, x=0.1, n_steps=24)  # 2k(1-x) = 1080 < dt_routing 3600
    logs = route_capturing_logs(net, tmp_path / 'unstable.log')
    assert 'not Muskingum-stable' in logs
    assert '5 of 5 rivers' in logs


def test_unstable_coefficients_can_raise(tmp_path):
    net = build_network(tmp_path / 'net', k=600.0, x=0.1, n_steps=24)
    with pytest.raises(ValueError, match='not Muskingum-stable'):
        rr.Router(
            forcing='vlateral',
            params_file=str(net.params_file),
            vlateral_files=[str(net.vlateral_file)],
            discharge_files=[net.path('q.nc')],
            channel_state_init_file=str(net.state_file),
            unstable_coefficients='raise',
            log=False,
            progress_bar=False,
        ).route()


def test_stable_network_does_not_warn(tmp_path):
    net = build_network(tmp_path / 'net', n_steps=24)
    logs = route_capturing_logs(net, tmp_path / 'stable.log')
    assert 'not Muskingum-stable' not in logs


# ── config handling ─────────────────────────────────────────────────────────


@pytest.mark.parametrize('threads', [0, -1, 2.0, True])
def test_invalid_threads_raise(network: SyntheticNetwork, threads):
    cfg = rr.Configs(params_file=str(network.params_file), discharge_files=[network.path('q.nc')], threads=threads)
    with pytest.raises(ValueError, match='threads must be'):
        cfg.validate()


def test_threads_without_a_thread_pool_route_single_threaded(network: SyntheticNetwork):
    """threads > 1 without a thread_pool creates no pool and routes the same discharge as a single thread."""
    single = route_vlateral(network, out_name='q_single.nc')
    threaded = route_vlateral(network, out_name='q_threaded.nc', threads=2)
    np.testing.assert_array_equal(threaded, single)


def test_route_leaves_a_given_thread_pool_open(network: SyntheticNetwork):
    """A thread_pool passed to route() belongs to the caller: it is used, then left for the with block to close."""
    router = rr.Router(
        forcing='vlateral',
        params_file=str(network.params_file),
        vlateral_files=[str(network.vlateral_file)],
        channel_state_init_file=str(network.state_file),
        dt_routing=network.dt_runoff,
        discharge_files=[network.path('q_pool.nc')],
        threads=2,
        log=False,
        progress_bar=False,
    )
    with ThreadPoolExecutor(2) as pool:
        router.route(thread_pool=pool)
        assert pool.submit(lambda: 1).result() == 1


def test_unknown_config_key_names_the_key(network: SyntheticNetwork):
    with pytest.raises(ValueError, match='Unrecognized config key'):
        rr.Router(forcing='channel', params_file=str(network.params_file), var_vlateral='vlateral')


def test_unknown_config_key_suggests_a_close_match(network: SyntheticNetwork):
    with pytest.raises(ValueError, match="did you mean 'dt_routing'"):
        rr.Router(forcing='channel', params_file=str(network.params_file), dt_routeing=3600)


def test_example_config_template_loads(network: SyntheticNetwork, tmp_path):
    """Every key in the shipped template must be a real config option."""
    import yaml

    template_path = Path(__file__).resolve().parent.parent / 'examples' / 'config.yaml'
    with open(template_path) as f:
        template = yaml.safe_load(f)
    template = {k: v for k, v in template.items() if v not in ('', [''], None)}
    template.update(
        params_file=str(network.params_file),
        discharge_files=[str(tmp_path / 'q.nc')],
        vlateral_files=[str(network.vlateral_file)] if template.get('forcing') == 'vlateral' else [],
        channel_state_init_file=str(network.state_file),
        dt_routing=900,
        dt_total=3600,
    )
    template.pop('discharge_dir', None)
    rr.Configs.from_mapping(template).validate()


def test_duplicate_input_basenames_raise(network: SyntheticNetwork, tmp_path):
    """Two inputs with the same name would resolve to one output file and silently overwrite each other."""
    nested = tmp_path / 'nested'
    nested.mkdir()
    write_vlateral(nested / 'vlateral.nc', np.zeros((24, network.n_rivers)), np.arange(1, network.n_rivers + 1))
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    with pytest.raises(ValueError, match='duplicate names'):
        rr.Router(
            forcing='vlateral',
            params_file=str(network.params_file),
            vlateral_files=[str(network.vlateral_file), str(nested / 'vlateral.nc')],
            discharge_dir=str(out_dir),
            log=False,
            progress_bar=False,
        )


def test_each_router_has_one_log_handler(network: SyntheticNetwork):
    """Loggers were named by id(self), which CPython reuses, so handlers accumulated across instances."""
    for _ in range(8):
        router = rr.Router(
            forcing='channel',
            params_file=str(network.params_file),
            discharge_files=[network.path('q.nc')],
            channel_state_init_file=str(network.state_file),
            dt_routing=3600,
            dt_total=3600,
            log=True,
            log_level='WARNING',
            progress_bar=False,
        )
        assert len(router.logger.handlers) == 1
        del router


def test_router_logger_does_not_propagate(network: SyntheticNetwork):
    """The router owns its handler, so propagating to the root logger would print every message twice."""
    router = rr.Router(
        forcing='channel',
        params_file=str(network.params_file),
        discharge_files=[network.path('q.nc')],
        channel_state_init_file=str(network.state_file),
        dt_routing=3600,
        dt_total=3600,
        log=True,
        progress_bar=False,
    )
    assert router.logger.propagate is False
