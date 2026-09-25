"""
Properties routing must have whatever the reference says: steady state is a fixed point, doubling the runoff doubles
the discharge, a restart continues the run it was saved from, water is conserved, and different ways of routing the
same runoff agree. Answers that should be identical must agree within TOLERANCE (see reference_solution.assert_same).
Most tests route the Willamette River to take seconds; those comparing the package's own runs use the whole Columbia.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from reference_solution import TOLERANCE, Basin, assert_same, route, route_run

import river_route as rr
from river_route.router import writers

DT = 3600  # the runoff's time step, and the routing step unless a test says otherwise


def accumulate_downstream(values: np.ndarray, downstream_indices: np.ndarray) -> np.ndarray:
    """Each river's value plus the values of every river upstream of it."""
    total = np.array(values, dtype=np.float64)
    for river, downstream in enumerate(downstream_indices):
        if downstream >= 0:
            total[downstream] += total[river]
    return total


def outlet_of(basin: Basin) -> int:
    return int(np.flatnonzero(basin.network.downstream_indices < 0)[0])


def read_discharge(package: Path, files: list[str]) -> np.ndarray:
    """A run's discharge files from the package, joined in time."""
    series = []
    for file in files:
        with xr.open_zarr(package / file) as discharge:
            series.append(discharge['Q'].to_numpy())
    return np.concatenate(series, axis=1)


################################################
# Correctness beyond the reference
################################################


@pytest.mark.parametrize(('network_type', 'dt_routing'), [('standard', 3600), ('standard', 1800), ('stabilized', 3600)])
def test_steady_state_is_a_fixed_point(network_type: str, dt_routing: int, willamette: Basin, tmp_path: Path) -> None:
    """
    With runoff constant in time, a river's inflow plus its own runoff rate is a fixed point of its recurrence, because
    c1 + c2 + c3 = 1. Starting every river there, its discharge must stay there for every step.
    """
    dates, volumes = willamette.months[0]
    constant = np.repeat(np.maximum(volumes.mean(axis=1, keepdims=True), 0), volumes.shape[1], axis=1)
    rate = constant[:, 0].astype(np.float64) / DT
    steady = accumulate_downstream(rate, willamette.network.downstream_indices)
    state = steady
    if network_type == 'stabilized':
        # the j-th of a river's N equal pieces carries its inflow plus j/N of its own runoff
        pieces, _, _ = willamette.network.conditioning(dt_routing)
        inflow = steady - rate
        state = np.concatenate([inflow[r] + rate[r] * np.arange(1, n + 1) / n for r, n in enumerate(pieces)])
    pd.DataFrame({'Q': state.astype(np.float32)}).to_parquet(tmp_path / 'steady.parquet')
    routed = route(
        tmp_path,
        [(dates, constant)],
        params_file=willamette.params_file,
        network_type=network_type,
        dt_routing=dt_routing,
        channel_state_init_file=tmp_path / 'steady.parquet',
    )
    assert_same(routed.discharge[0], np.repeat(steady[:, np.newaxis], constant.shape[1], axis=1))


@pytest.mark.parametrize('network_type', ['standard', 'stabilized'])
def test_doubling_the_runoff_doubles_the_discharge(network_type: str, willamette: Basin, tmp_path: Path) -> None:
    dates, volumes = willamette.months[0]
    once = route(tmp_path, [(dates, volumes)], params_file=willamette.params_file, network_type=network_type)
    twice = route(tmp_path, [(dates, 2 * volumes)], params_file=willamette.params_file, network_type=network_type)
    assert_same(twice.discharge[0], 2 * once.discharge[0])
    assert_same(twice.final_state, 2 * once.final_state)


@pytest.mark.parametrize('name', ['static_standard', 'static_stabilized'])
def test_a_restart_continues_the_run_it_was_saved_from(
    name: str, package: Path, manifest: dict, tmp_path: Path
) -> None:
    """Route March and April, then May and June from the state saved after April, and compare with the reference run
    that routed all four months at once."""
    run = manifest['runs'][name]
    runoff_files = run['configs']['runoff_files']
    first = run | {'configs': run['configs'] | {'runoff_files': runoff_files[:2]}}
    route_run(package, first, [tmp_path / 'march', tmp_path / 'april'], tmp_path / 'april.parquet', writers.null_writer)
    resumed = run | {
        'configs': run['configs']
        | {'runoff_files': runoff_files[2:], 'channel_state_init_file': tmp_path / 'april.parquet'}
    }
    continued = {}

    def keep_discharge(router, dates, discharge_array, discharge_file, runoff_file=''):
        continued[Path(discharge_file).name] = discharge_array.copy()

    route_run(package, resumed, [tmp_path / 'may', tmp_path / 'june'], tmp_path / 'june.parquet', keep_discharge)
    assert_same(
        np.concatenate([continued['may'], continued['june']], axis=1), read_discharge(package, run['discharge'][2:])
    )
    assert_same(pd.read_parquet(tmp_path / 'june.parquet')['Q'], pd.read_parquet(package / run['final_state'])['Q'])


def test_months_in_sequence_match_one_combined_series(package: Path, manifest: dict, tmp_path: Path) -> None:
    """Route the four months of catchment runoff as one file, and compare with the reference run that routed them as
    four files in sequence, handing the channel state from each to the next."""
    run = manifest['runs']['catchment_file']
    months = list(rr.CatchmentRunoff().generator([package / path for path in run['configs']['runoff_files']]))
    dates = np.concatenate([month_dates for month_dates, _, _ in months])
    volumes = np.concatenate([month.runoff for _, month, _ in months], axis=1)
    routed = route(tmp_path, [(dates, volumes)], params_file=package / run['configs']['params_file'])
    assert_same(routed.discharge[0], read_discharge(package, run['discharge']))


def test_headwaters_follow_the_muskingum_formula(willamette: Basin, tmp_path: Path) -> None:
    """A stable headwater has no inflow, so its discharge is the recurrence q = c3 q + (c1 + c2) R, computed here in
    float64 from its k and x."""
    network = willamette.network
    dates, volumes = willamette.months[0]
    routed = route(tmp_path, [(dates, volumes)], params_file=willamette.params_file)
    has_upstream = np.zeros(network.size, dtype=bool)
    has_upstream[network.downstream_indices[network.downstream_indices >= 0]] = True
    too_long, too_short = network.unstable_mask(DT)
    rivers = np.flatnonzero(~has_upstream & ~too_long & ~too_short)
    assert rivers.size, 'the basin has no stable headwater to check'
    dt_div_k = DT / network.k[rivers].astype(np.float64)
    two_x = 2 * network.x[rivers].astype(np.float64)
    denominator = dt_div_k + 2 - two_x
    c1_plus_c2, c3 = 2 * dt_div_k / denominator, (2 - two_x - dt_div_k) / denominator
    q = np.zeros(rivers.size)
    expected = np.empty((rivers.size, volumes.shape[1]))
    for step in range(volumes.shape[1]):
        q = c3 * q + c1_plus_c2 * volumes[rivers, step] / DT
        expected[:, step] = np.maximum(q, 0)
    assert_same(routed.discharge[0][rivers], expected)


@pytest.mark.parametrize('network_type', ['stabilized', 'standard'])
def test_water_is_conserved(network_type: str, willamette: Basin, tmp_path: Path) -> None:
    """Route two months of runoff, then two months of none so the basin drains: all the runoff must leave through the
    outlet."""
    wet = [(dates, np.maximum(volumes, 0)) for dates, volumes in willamette.months[:2]]
    dry = [(dates, np.zeros_like(volumes)) for dates, volumes in willamette.months[2:]]
    routed = route(tmp_path, wet + dry, params_file=willamette.params_file, network_type=network_type)
    volume_in = sum(volumes.sum(dtype=np.float64) for _, volumes in wet)
    volume_out = sum(discharge[outlet_of(willamette)].sum(dtype=np.float64) for discharge in routed.discharge) * DT
    assert abs(volume_out - volume_in) <= TOLERANCE * volume_in, f'{volume_out:.6g} m³ left of {volume_in:.6g} m³'


def test_a_pulse_keeps_its_volume_and_arrives_after_the_travel_time(willamette: Basin, tmp_path: Path) -> None:
    """
    Put one pulse of runoff into the headwater farthest from the outlet. All of it must reach the outlet, and its center
    of mass must arrive after the Muskingum travel time.

    A discrete Muskingum reach delays what flows in at its top by exactly k, so each river below the headwater adds its
    k. Runoff enters along a reach instead, and one reach delays it by k (1 - x). The headwater is split into N equal
    pieces that each take 1/N of its runoff, so its own runoff is delayed by k ((1 - x) / N + (N - 1) / 2N) on average.
    """
    network = willamette.network
    downstream = network.downstream_indices
    below = np.zeros(network.size)  # seconds from a river's bottom to the outlet: the k of every river below it
    for river in range(network.size - 1, -1, -1):
        if downstream[river] >= 0:
            below[river] = network.k[downstream[river]] + below[downstream[river]]
    headwater = int(np.argmax(below + network.k))
    k, x = float(network.k[headwater]), float(network.x[headwater])
    pieces = int(network.conditioning(DT)[0][headwater])
    travel_time = below[headwater] + k * ((1 - x) / pieces + (pieces - 1) / (2 * pieces))
    months = [(dates, np.zeros_like(volumes)) for dates, volumes in willamette.months]
    pulse = 1e6  # m³
    months[0][1][headwater, 0] = pulse
    routed = route(tmp_path, months, params_file=willamette.params_file, network_type='stabilized')
    outflow = np.concatenate([discharge[outlet_of(willamette)] for discharge in routed.discharge]).astype(np.float64)
    assert abs(outflow.sum() * DT - pulse) <= TOLERANCE * pulse
    step_ends = np.arange(1, outflow.size + 1) * DT
    delay = (outflow * step_ends).sum() / outflow.sum() - DT / 2  # the pulse enters at a constant rate over one step
    assert abs(delay - travel_time) <= TOLERANCE * travel_time, (
        f'arrived after {delay / 3600:.2f} h against a travel time of {travel_time / 3600:.2f} h'
    )


################################################
# Consistency between routing paths
################################################


def test_gridded_and_catchment_file_runoff_route_the_same(package: Path, manifest: dict) -> None:
    grid, catchment = manifest['runs']['static_standard'], manifest['runs']['catchment_file']
    assert_same(read_discharge(package, catchment['discharge']), read_discharge(package, grid['discharge']))


def test_threaded_routing_is_repeatable(willamette: Basin, tmp_path: Path) -> None:
    first = route(tmp_path, willamette.months, params_file=willamette.params_file, threads=4)
    second = route(tmp_path, willamette.months, params_file=willamette.params_file, threads=4)
    assert len(first.router.routing_jobs) > 1, 'the basin was routed as one region, so no threads were used'
    for first_discharge, second_discharge in zip(first.discharge, second.discharge, strict=True):
        np.testing.assert_array_equal(first_discharge, second_discharge)
    np.testing.assert_array_equal(first.final_state, second.final_state)


@pytest.mark.parametrize('name', ['static_standard', 'static_stabilized'])
def test_threaded_and_single_threaded_runs_agree(name: str, package: Path, manifest: dict) -> None:
    single, threaded = manifest['runs'][name], manifest['runs'][f'{name}_4_threads']
    assert_same(read_discharge(package, threaded['discharge']), read_discharge(package, single['discharge']))


def test_split_networks_route_like_kernel_stabilized_networks(willamette: Basin, tmp_path: Path) -> None:
    """
    A network stabilized by writing its sub-reaches as rows must route like network_type='stabilized', which splits
    them inside the kernel. The written network does not sub-cycle rivers too short for the routing step, so only
    rivers with no such river upstream of them are compared.
    """
    split_network = rr.Network(willamette.params_file).stabilize(DT)
    split = route(tmp_path, network=split_network, **willamette.gridded(months=2))
    in_kernel = route(tmp_path, network_type='stabilized', **willamette.gridded(months=2))
    _, too_short = willamette.network.unstable_mask(DT)
    comparable = accumulate_downstream(too_short, willamette.network.downstream_indices) == 0
    assert comparable.sum() > willamette.network.size / 2
    for split_discharge, kernel_discharge in zip(split.discharge, in_kernel.discharge, strict=True):
        assert_same(split_discharge[comparable], kernel_discharge[comparable])


def test_daily_output_is_the_mean_of_each_day(willamette: Basin, tmp_path: Path) -> None:
    hourly = route(tmp_path, willamette.months[:1], params_file=willamette.params_file)
    daily = route(tmp_path, willamette.months[:1], params_file=willamette.params_file, dt_discharge=86400)
    n_rivers, n_hours = hourly.discharge[0].shape
    assert_same(daily.discharge[0], hourly.discharge[0].reshape(n_rivers, n_hours // 24, 24).mean(axis=2))
    np.testing.assert_array_equal(daily.dates[0], hourly.dates[0][::24])


def test_dt_total_routes_the_first_steps(willamette: Basin, tmp_path: Path) -> None:
    full = route(tmp_path, willamette.months[:1], params_file=willamette.params_file)
    first_days = route(tmp_path, willamette.months[:1], params_file=willamette.params_file, dt_total=10 * 86400)
    assert_same(first_days.discharge[0], full.discharge[0][:, : 10 * 24])


def test_ensemble_members_start_from_the_same_state(willamette: Basin, tmp_path: Path) -> None:
    ensemble = route(tmp_path, willamette.months, params_file=willamette.params_file, runoff_processing_mode='ensemble')
    members = [route(tmp_path, [month], params_file=willamette.params_file) for month in willamette.months]
    for ensemble_discharge, member in zip(ensemble.discharge, members, strict=True):
        assert_same(ensemble_discharge, member.discharge[0])
    assert_same(ensemble.final_state, np.mean([member.final_state for member in members], axis=0))
