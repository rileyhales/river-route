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


@pytest.mark.parametrize('dt_routing', [3600, 1800, 'auto'])
def test_steady_state_is_a_fixed_point(dt_routing: int | str, willamette: Basin, tmp_path: Path) -> None:
    """
    With runoff constant in time, a river's inflow plus its own runoff rate is a fixed point of its recurrence, because
    c1 + c2 + c3 = 1. Starting every river there, its discharge must stay there for every step.
    """
    dates, volumes = willamette.months[0]
    constant = np.repeat(np.maximum(volumes.mean(axis=1, keepdims=True), 0), volumes.shape[1], axis=1)
    rate = constant[:, 0].astype(np.float64) / DT
    steady = accumulate_downstream(rate, willamette.network.downstream_indices)
    state, river_ids = steady, willamette.network.river_ids
    if dt_routing == 'auto':
        # each sub-reach of a river carries its inflow plus the shares of its runoff of it and the sub-reaches above it
        sub_reaches, _, _ = willamette.network.conditioning(DT)
        k_reach, _ = willamette.network.sub_reach_k_and_x(sub_reaches)
        share = k_reach / np.repeat(willamette.network.k.astype(np.float64), sub_reaches)
        first_sub_reach = np.concatenate(([0], np.cumsum(sub_reaches)[:-1]))
        runoff_above = np.cumsum(share) - np.repeat(
            np.cumsum(share)[first_sub_reach] - share[first_sub_reach], sub_reaches
        )
        state = np.repeat(steady - rate, sub_reaches) + np.repeat(rate, sub_reaches) * runoff_above
        river_ids = np.repeat(river_ids, sub_reaches)
    pd.DataFrame({'riverId': river_ids, 'Q': state.astype(np.float32)}).to_parquet(tmp_path / 'steady.parquet')
    routed = route(
        tmp_path,
        [(dates, constant)],
        network_file=willamette.network_file,
        dt_routing=dt_routing,
        channel_state_init_file=tmp_path / 'steady.parquet',
    )
    assert_same(routed.discharge[0], np.repeat(steady[:, np.newaxis], constant.shape[1], axis=1))


@pytest.mark.parametrize('dt_routing', [DT, 'auto'])
def test_doubling_the_runoff_doubles_the_discharge(dt_routing: int | str, willamette: Basin, tmp_path: Path) -> None:
    dates, volumes = willamette.months[0]
    once = route(tmp_path, [(dates, volumes)], network_file=willamette.network_file, dt_routing=dt_routing)
    twice = route(tmp_path, [(dates, 2 * volumes)], network_file=willamette.network_file, dt_routing=dt_routing)
    assert_same(twice.discharge[0], 2 * once.discharge[0])
    assert_same(twice.final_state, 2 * once.final_state)


@pytest.mark.parametrize('name', ['static_standard', 'static_dt_routing'])
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

    def keep_discharge(router, dates, discharge_array, discharge_file, runoff_file='', *, thread_pool=None, threads=1):
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
    routed = route(tmp_path, [(dates, volumes)], network_file=package / run['configs']['network_file'])
    assert_same(routed.discharge[0], read_discharge(package, run['discharge']))


@pytest.mark.parametrize('dt_routing', [3600, 900])
def test_headwaters_follow_the_muskingum_formula(dt_routing: int, willamette: Basin, tmp_path: Path) -> None:
    """
    A headwater whose coefficients are all non-negative has no inflow, so its discharges follow the recurrence
    q = c3 q + (c1 + c2) R, computed here in float64 from its k and x. The discharge written for each runoff step of m
    routing steps is its discharge at the end of the runoff step, q_m.
    """
    network = willamette.network
    dates, volumes = willamette.months[0]
    routed = route(tmp_path, [(dates, volumes)], network_file=willamette.network_file, dt_routing=dt_routing)
    has_upstream = np.zeros(network.size, dtype=bool)
    has_upstream[network.downstream_indices[network.downstream_indices >= 0]] = True
    too_long, too_short = network.invalid_mask(dt_routing)
    rivers = np.flatnonzero(~has_upstream & ~too_long & ~too_short)
    assert rivers.size, 'the basin has no headwater with non-negative coefficients to check'
    dt_div_k = dt_routing / network.k[rivers].astype(np.float64)
    two_x = 2 * network.x[rivers].astype(np.float64)
    denominator = dt_div_k + 2 - two_x
    c1_plus_c2, c3 = 2 * dt_div_k / denominator, (2 - two_x - dt_div_k) / denominator
    steps = DT // dt_routing
    discharges = np.zeros((rivers.size, volumes.shape[1] * steps + 1))
    for step in range(volumes.shape[1] * steps):
        discharges[:, step + 1] = c3 * discharges[:, step] + c1_plus_c2 * volumes[rivers, step // steps] / DT
    assert_same(routed.discharge[0][rivers], discharges[:, steps::steps])


@pytest.mark.parametrize('dt_routing', ['auto', DT])
def test_water_is_conserved(dt_routing: int | str, willamette: Basin, tmp_path: Path) -> None:
    """Route two months of runoff, then two months of none so the basin drains: all the runoff must leave through the
    outlet."""
    wet = [(dates, np.maximum(volumes, 0)) for dates, volumes in willamette.months[:2]]
    dry = [(dates, np.zeros_like(volumes)) for dates, volumes in willamette.months[2:]]
    routed = route(tmp_path, wet + dry, network_file=willamette.network_file, dt_routing=dt_routing)
    volume_in = sum(volumes.sum(dtype=np.float64) for _, volumes in wet)
    volume_out = sum(discharge[outlet_of(willamette)].sum(dtype=np.float64) for discharge in routed.discharge) * DT
    assert abs(volume_out - volume_in) <= TOLERANCE * volume_in, f'{volume_out:.6g} m³ left of {volume_in:.6g} m³'


def test_a_pulse_keeps_its_volume_and_arrives_after_the_travel_time(willamette: Basin, tmp_path: Path) -> None:
    """
    Put one pulse of runoff into the headwater farthest from the outlet. All of it must reach the outlet, and its center
    of mass must arrive after the Muskingum travel time.

    A discrete Muskingum reach delays what flows in at its top by exactly k, so each river below the headwater adds its
    k. Runoff enters along a reach instead, and one reach delays it by k (1 - x). A headwater split into the
    sub-reaches of Network.sub_reach_k_and_x, each taking its share of the river's k and of its runoff, delays the
    runoff of each sub-reach by k_i (1 - x_i) through it and by the k of every sub-reach below it. The outlet writes
    its discharge at the end of each step.
    """
    network = willamette.network
    downstream = network.downstream_indices
    below = np.zeros(network.size)  # seconds from a river's bottom to the outlet: the k of every river below it
    for river in range(network.size - 1, -1, -1):
        if downstream[river] >= 0:
            below[river] = network.k[downstream[river]] + below[downstream[river]]
    headwater, outlet = int(np.argmax(below + network.k)), outlet_of(willamette)
    sub_reaches = network.conditioning(DT)[0]
    k_reach, x_reach = network.sub_reach_k_and_x(sub_reaches)
    first_sub_reach = int(sub_reaches[:headwater].sum())
    k, x = (values[first_sub_reach : first_sub_reach + sub_reaches[headwater]] for values in (k_reach, x_reach))
    below_in_river = k[::-1].cumsum()[::-1] - k  # the k of the headwater's sub-reaches below each one
    travel_time = below[headwater] + np.sum(k / k.sum() * (k * (1 - x) + below_in_river))
    months = [(dates, np.zeros_like(volumes)) for dates, volumes in willamette.months]
    pulse = 1e6  # m³
    months[0][1][headwater, 0] = pulse
    routed = route(tmp_path, months, network_file=willamette.network_file)
    outflow = np.concatenate([discharge[outlet] for discharge in routed.discharge]).astype(np.float64)
    assert abs(outflow.sum() * DT - pulse) <= TOLERANCE * pulse
    ends = (np.arange(outflow.size) + 1) * DT
    delay = (outflow * ends).sum() / outflow.sum() - DT / 2  # the pulse enters at a constant rate over one step
    assert abs(delay - travel_time) <= TOLERANCE * travel_time, (
        f'arrived after {delay / 3600:.2f} h against a travel time of {travel_time / 3600:.2f} h'
    )


@pytest.mark.parametrize(('dt_routing', 'routed_reaches'), [(DT, 2), ('auto', 3)])
def test_a_long_river_delays_and_spreads_a_pulse_by_its_own_moments(
    dt_routing: int | str, routed_reaches: int, tmp_path: Path
) -> None:
    """
    Whatever the step, the response of a Muskingum reach to what flows in at its top has mean k and variance
    k^2 (1 - 2x), and reaches in series add their moments. A river of k = 20 h and x = 0.2 is too long for a step of
    1 h, and dt_routing auto routes it in two sub-reaches. Equal ones would need x = 1/2 - 2 (1/2 - 0.2) = -0.1 to keep
    its variance, and Muskingum x is never negative, so they are linear reservoirs of 14.5 h and 5.5 h, whose squares
    sum to the river's 240 h^2. A pulse of runoff on a headwater flows into it, and the moments of the river's discharge
    exceed the headwater's by its own.
    """
    k, x = 20.0, 0.2  # hours
    table = pd.DataFrame(
        {
            'riverId': [1, 2],
            'nextRiverId': [2, -1],
            'muskingumK': np.float32([2 * DT, k * DT]),
            'muskingumX': np.float32([x, x]),
            'riverIndex': np.int32([0, 1]),
            'upstreamCount': np.int32([0, 1]),
        }
    )
    table.to_parquet(tmp_path / 'chain.parquet')
    dates = np.datetime64('2000-01-01T00:00:00', 's') + np.arange(1000) * np.timedelta64(DT, 's')
    volumes = np.zeros((2, dates.size), dtype=np.float32)
    volumes[0, 0] = 1e6
    routed = route(tmp_path, [(dates, volumes)], network_file=tmp_path / 'chain.parquet', dt_routing=dt_routing)
    assert routed.router.channel_state.shape[0] == routed_reaches

    def mean_and_variance(discharge: np.ndarray) -> tuple[float, float]:
        weights = discharge.astype(np.float64) / discharge.sum(dtype=np.float64)
        hours = np.arange(discharge.size)
        mean = float((weights * hours).sum())
        return mean, float((weights * (hours - mean) ** 2).sum())

    mean_in, variance_in = mean_and_variance(routed.discharge[0][0])
    mean_out, variance_out = mean_and_variance(routed.discharge[0][1])
    np.testing.assert_allclose(mean_out - mean_in, k, rtol=1e-4)
    np.testing.assert_allclose(variance_out - variance_in, k**2 * (1 - 2 * x), rtol=1e-3)


def test_gridded_and_catchment_file_runoff_route_the_same(package: Path, manifest: dict) -> None:
    grid, catchment = manifest['runs']['static_standard'], manifest['runs']['catchment_file']
    assert_same(read_discharge(package, catchment['discharge']), read_discharge(package, grid['discharge']))


def test_threaded_routing_is_repeatable(willamette: Basin, tmp_path: Path) -> None:
    first = route(tmp_path, willamette.months, network_file=willamette.network_file, threads=4)
    second = route(tmp_path, willamette.months, network_file=willamette.network_file, threads=4)
    routing_blocks, _ = first.router.network.routing_blocks(threads=4)
    assert len(routing_blocks) > 1, 'the basin was routed as one block, so no threads were used'
    for first_discharge, second_discharge in zip(first.discharge, second.discharge, strict=True):
        np.testing.assert_array_equal(first_discharge, second_discharge)
    np.testing.assert_array_equal(first.final_state, second.final_state)


@pytest.mark.parametrize('name', ['static_standard', 'static_dt_routing'])
def test_threaded_and_single_threaded_runs_agree(name: str, package: Path, manifest: dict) -> None:
    single, threaded = manifest['runs'][name], manifest['runs'][f'{name}_4_threads']
    assert_same(read_discharge(package, threaded['discharge']), read_discharge(package, single['discharge']))


def test_split_networks_route_like_dt_routing_auto(willamette: Basin, tmp_path: Path) -> None:
    """
    A network stabilized by writing its sub-reaches as rows must route like dt_routing auto, which routes them as
    sub-reaches inside the kernel. The written network does not route rivers too short for the routing step in
    subcycles, so only rivers with no such river upstream of them are compared.
    """
    split_file = rr.Network(willamette.network_file).write_stabilized(DT, tmp_path / 'split.parquet')
    split = route(tmp_path, **(willamette.gridded(months=2) | {'network_file': split_file}))
    in_kernel = route(tmp_path, **willamette.gridded(months=2))
    _, too_short = willamette.network.invalid_mask(DT)
    comparable = accumulate_downstream(too_short, willamette.network.downstream_indices) == 0
    assert comparable.sum() > willamette.network.size / 2
    for split_discharge, kernel_discharge in zip(split.discharge, in_kernel.discharge, strict=True):
        assert_same(split_discharge[comparable], kernel_discharge[comparable])


@pytest.mark.parametrize('dt_routing', [DT, 'auto'])
def test_dynamic_without_exponents_routes_like_static(dt_routing: int | str, willamette: Basin, tmp_path: Path) -> None:
    """
    With exponents of 0 a river's k = dynamicKCoefficient * Q ** dynamicKExponent and
    x = 1/2 - dynamicXCoefficient * Q ** dynamicXExponent are its coefficients whatever it carries, so with
    dynamicKCoefficient = k and dynamicXCoefficient = 1/2 - x the dynamic coefficients of every step are the static
    ones. With dt_routing auto the plan of each river is then the static one too: the same sub-reaches of the same
    shares, and at every step the static river's subcycles.
    """
    table = pd.read_parquet(willamette.network_file)
    dynamic = table.assign(
        dynamicKCoefficient=table['muskingumK'],
        dynamicKExponent=np.float32(0),
        dynamicXCoefficient=np.float32(0.5) - table['muskingumX'],
        dynamicXExponent=np.float32(0),
        dynamicMinQ=np.float32(1),
        dynamicMaxQ=np.float32(1),
    )
    dynamic.to_parquet(tmp_path / 'dynamic.parquet')
    options = {'dt_routing': dt_routing}
    static = route(tmp_path, willamette.months[:1], network_file=willamette.network_file, **options)
    dynamic = route(
        tmp_path, willamette.months[:1], network_file=tmp_path / 'dynamic.parquet', parameters='dynamic', **options
    )
    np.testing.assert_array_equal(
        dynamic.router.compute_plan.sub_reach_starts, static.router.compute_plan.sub_reach_starts
    )
    assert_same(dynamic.discharge[0], static.discharge[0])
    assert_same(dynamic.final_state, static.final_state)


def test_dt_routing_auto_keeps_dynamic_discharge_from_going_negative(willamette: Basin, tmp_path: Path) -> None:
    """
    Muskingum coefficients that are all non-negative route non-negative runoff into non-negative discharge. With
    dt_routing auto, rivers whose k and x change with discharge are split and subcycled so that theirs are at every
    step, for laws whose k is each river's own at 10 m³/s and falls with discharge, and whose x rises with it.
    """
    table = pd.read_parquet(willamette.network_file)
    reference_q = 10.0
    table.assign(
        dynamicKCoefficient=(table['muskingumK'] * reference_q**0.3).astype(np.float32),
        dynamicKExponent=np.float32(-0.3),
        dynamicXCoefficient=np.float32(0.2 / reference_q**0.2),  # x = 0.3 at 10 m³/s
        dynamicXExponent=np.float32(-0.2),
        dynamicMinQ=np.float32(0.1),
        dynamicMaxQ=np.float32(10_000),
    ).to_parquet(tmp_path / 'dynamic.parquet')
    assert all(np.all(volumes >= 0) for _, volumes in willamette.months[:1])
    fixed_step, auto = (
        route(
            tmp_path,
            willamette.months[:1],
            network_file=tmp_path / 'dynamic.parquet',
            parameters='dynamic',
            dt_routing=dt_routing,
        )
        for dt_routing in (DT, 'auto')
    )
    assert np.any(fixed_step.discharge[0] < 0), 'every river is valid routed whole, so this checks nothing'
    assert np.all(auto.discharge[0] >= 0)
    assert np.all(auto.final_state >= 0)


@pytest.mark.parametrize('restrict_x_range', [True, False])
def test_restrict_x_range_routes_no_x_above_its_limit(
    restrict_x_range: bool, willamette: Basin, tmp_path: Path
) -> None:
    """
    With restrict_x_range, rivers of x = 0.49 are routed with x = 0.45 by static and dynamic parameters alike, and
    without it with their own x.
    """
    table = pd.read_parquet(willamette.network_file).assign(muskingumX=np.float32(0.49))
    table.assign(
        dynamicKCoefficient=table['muskingumK'],
        dynamicKExponent=np.float32(0),
        dynamicXCoefficient=np.float32(0.01),  # 1/2 - 0.49
        dynamicXExponent=np.float32(0),
        dynamicMinQ=np.float32(1),
        dynamicMaxQ=np.float32(1),
    ).to_parquet(tmp_path / 'x_of_049.parquet')
    routed_x = 0.45 if restrict_x_range else 0.49
    table.assign(muskingumX=np.float32(routed_x)).to_parquet(tmp_path / 'routed_x.parquet')
    options = {'network_file': tmp_path / 'x_of_049.parquet', 'dt_routing': DT, 'restrict_x_range': restrict_x_range}
    expected = route(
        tmp_path,
        willamette.months[:1],
        network_file=tmp_path / 'routed_x.parquet',
        dt_routing=DT,
        restrict_x_range=False,
    )
    for parameters in ('static', 'dynamic'):
        routed = route(tmp_path, willamette.months[:1], parameters=parameters, **options)
        assert_same(routed.discharge[0], expected.discharge[0])
        np.testing.assert_array_equal(routed.router.network.x.max(), np.float32(routed_x))


def test_daily_output_is_the_discharge_at_the_end_of_each_day(willamette: Basin, tmp_path: Path) -> None:
    hourly = route(tmp_path, willamette.months[:1], network_file=willamette.network_file)
    daily = route(tmp_path, willamette.months[:1], network_file=willamette.network_file, dt_discharge=86400)
    assert_same(daily.discharge[0], hourly.discharge[0][:, 23::24])
    np.testing.assert_array_equal(daily.dates[0], hourly.dates[0][::24])


def test_dt_total_routes_the_first_steps(willamette: Basin, tmp_path: Path) -> None:
    full = route(tmp_path, willamette.months[:1], network_file=willamette.network_file)
    first_days = route(tmp_path, willamette.months[:1], network_file=willamette.network_file, dt_total=10 * 86400)
    assert_same(first_days.discharge[0], full.discharge[0][:, : 10 * 24])


def test_ensemble_members_start_from_the_same_state(willamette: Basin, tmp_path: Path) -> None:
    ensemble = route(
        tmp_path, willamette.months, network_file=willamette.network_file, runoff_processing_mode='ensemble'
    )
    members = [route(tmp_path, [month], network_file=willamette.network_file) for month in willamette.months]
    for ensemble_discharge, member in zip(ensemble.discharge, members, strict=True):
        assert_same(ensemble_discharge, member.discharge[0])
    assert_same(ensemble.final_state, np.mean([member.final_state for member in members], axis=0))
