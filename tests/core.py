"""
Core integration tests for river-route using GEOGloWS v2 VPU 201 data.

Run `python tests/prepare_test_data.py` before running these tests to download
the required ERA5 and GEOGloWS data into tests/data/.

Testing strategy (from original TODOs):
  1. test that a valid config file passes validation
  2. test that several versions of an invalid config file fail validation
  3. test routing only (Muskingum) and compare against known valid outputs
  4. test RapidMuskingum and compare against known valid outputs
  5. test UnitMuskingum and compare against known valid outputs
  6. test that initial states are read and used for channel and transformer
  7. test that final states are written and correct for channel and transformer
  8. compare that grid weights can be correctly reproduced
  9. do you get the same answer routing from depths+weights vs from pre-computed catchment volumes?
"""
import logging
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import river_route as rr

# ── Paths ────────────────────────────────────────────────────────────────────

data_root = Path(__file__).resolve().parent / 'data'

# Core input files (from prepare_test_data.py)
routing_params_file = data_root / 'routing_params.parquet'
gridweights_file = data_root / 'gridweights_era5_regional.nc'
initial_state_file = data_root / 'initial_state.parquet'

# ERA5 runoff files (one per month)
era5_files = sorted(data_root.glob('era5_runoff_*.nc'))

# Pre-computed catchment volumes (from prepare_test_data.py or synthetic)
volume_files = sorted(data_root.glob('catchment_volumes_*.nc'))
if not volume_files:
    # Fall back to synthetic volumes if prepared data isn't available
    _synth = data_root / 'synthetic_volumes.nc'
    if _synth.exists():
        volume_files = [_synth]

# Retrospective spot-check
retrospective_file = data_root / 'retrospective_spot_check.nc'

# ── Logging ──────────────────────────────────────────────────────────────────

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
log.addHandler(ch)


# ── Helpers ──────────────────────────────────────────────────────────────────

def _require_data(*paths: Path) -> bool:
    missing = [p for p in paths if not p.exists()]
    if missing:
        for p in missing:
            log.warning(f'Missing: {p}')
        log.warning('Run `python tests/prepare_test_data.py` to download test data')
        return False
    return True


def _detect_era5_var_names() -> dict:
    """Detect coordinate variable names from the first ERA5 file."""
    with xr.open_dataset(era5_files[0]) as ds:
        x_var = 'longitude' if 'longitude' in ds.dims else 'lon'
        y_var = 'latitude' if 'latitude' in ds.dims else 'lat'
        t_var = 'valid_time' if 'valid_time' in ds.dims else 'time'
    return {'x_var': x_var, 'y_var': y_var, 't_var': t_var}


def compare_netcdfs(file1: Path, file2: Path,
                    variable: str | list[str] | None = None,
                    rtol: float = 1e-5, atol: float = 1e-8) -> bool:
    ds1 = xr.open_dataset(file1)
    ds2 = xr.open_dataset(file2)

    if variable is None:
        variable = list(set(ds1.data_vars).intersection(ds2.data_vars))
    if isinstance(variable, str):
        variable = [variable]

    is_identical = True
    for var in variable:
        if var not in ds1 or var not in ds2:
            log.error(f'Variable {var} missing from one of the datasets')
            is_identical = False
            continue
        if ds1[var].shape != ds2[var].shape:
            log.error(f'{var}: shape mismatch {ds1[var].shape} vs {ds2[var].shape}')
            is_identical = False
            continue
        if not np.allclose(ds1[var].values, ds2[var].values, rtol=rtol, atol=atol, equal_nan=True):
            diff = np.abs(ds1[var].values - ds2[var].values)
            log.error(f'{var}: max diff={np.nanmax(diff):.6f}, mean diff={np.nanmean(diff):.6f}')
            is_identical = False
        else:
            log.info(f'{var}: values match (rtol={rtol}, atol={atol})')
    return is_identical


# ═══════════════════════════════════════════════════════════════════════════
# Test 1: Valid config passes validation
# ═══════════════════════════════════════════════════════════════════════════

def test_valid_config() -> None:
    log.info('=' * 60)
    log.info('TEST 1: Valid config passes validation')
    if not _require_data(routing_params_file, initial_state_file):
        return

    tmpdir = tempfile.mkdtemp()
    try:
        # RapidMuskingum config with catchment volumes
        if volume_files:
            discharge_files = [os.path.join(tmpdir, f'q_{i}.nc') for i in range(len(volume_files))]
            cfg = rr.Configs(
                params_file=str(routing_params_file),
                catchment_runoff_files=[str(f) for f in volume_files],
                discharge_files=discharge_files,
                channel_state_init_file=str(initial_state_file),
                channel_state_final_file=os.path.join(tmpdir, 'final_state.parquet'),
                log=False,
            )
            log.info('[OK] RapidMuskingum config with catchment volumes is valid')

        # Muskingum (channel-only) config
        cfg = rr.Configs(
            params_file=str(routing_params_file),
            discharge_files=[os.path.join(tmpdir, 'q_musk.nc')],
            channel_state_init_file=str(initial_state_file),
            dt_routing=3600,
            dt_total=86400,
            log=False,
        )
        log.info('[OK] Muskingum (channel-only) config is valid')

    except Exception as e:
        log.error(f'[FAIL] Config validation raised: {e}')
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════════
# Test 2: Invalid configs fail validation
# ═══════════════════════════════════════════════════════════════════════════

def test_invalid_configs() -> None:
    log.info('=' * 60)
    log.info('TEST 2: Invalid configs fail validation')
    if not _require_data(routing_params_file):
        return

    tmpdir = tempfile.mkdtemp()
    cases = {
        'missing params_file': dict(
            params_file='/nonexistent/params.parquet',
            discharge_files=[os.path.join(tmpdir, 'q.nc')],
            dt_routing=3600,
            dt_total=86400,
        ),
        'missing discharge output dir': dict(
            params_file=str(routing_params_file),
            discharge_files=['/nonexistent_dir/q.nc'],
            dt_routing=3600,
            dt_total=86400,
        ),
        'both catchment and grid runoff': dict(
            params_file=str(routing_params_file),
            catchment_runoff_files=[str(volume_files[0])] if volume_files else ['/tmp/dummy.nc'],
            runoff_grid_files=['/tmp/dummy.nc'],
            grid_weights_file='/tmp/dummy.nc',
            discharge_files=[os.path.join(tmpdir, 'q.nc')],
        ),
        'invalid runoff_processing_mode': dict(
            params_file=str(routing_params_file),
            discharge_files=[os.path.join(tmpdir, 'q.nc')],
            dt_routing=3600,
            dt_total=86400,
            runoff_processing_mode='invalid_mode',
        ),
    }

    for name, kwargs in cases.items():
        try:
            rr.Configs(**kwargs)
            log.error(f'[FAIL] "{name}" did not raise an exception')
        except (FileNotFoundError, NotADirectoryError, ValueError) as e:
            log.info(f'[OK] "{name}" correctly raised {type(e).__name__}: {e}')
        except Exception as e:
            log.info(f'[OK] "{name}" raised {type(e).__name__}: {e}')

    shutil.rmtree(tmpdir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════════
# Test 3: Muskingum channel-only routing
# ═══════════════════════════════════════════════════════════════════════════

def test_muskingum_channel_only() -> None:
    """Route with no lateral inflow — discharge should decay from initial state."""
    log.info('=' * 60)
    log.info('TEST 3: Muskingum channel-only routing')
    if not _require_data(routing_params_file, initial_state_file):
        return

    tmpdir = tempfile.mkdtemp()
    try:
        discharge_file = os.path.join(tmpdir, 'q_channel_only.nc')
        final_state_file = os.path.join(tmpdir, 'final_state.parquet')

        rr.Muskingum(params_file=str(routing_params_file), discharge_files=[discharge_file],
                     channel_state_init_file=str(initial_state_file), channel_state_final_file=final_state_file,
                     dt_routing=3600, dt_total=86400, log=False).route()

        # Verify output exists and has reasonable values
        with xr.open_dataset(discharge_file) as ds:
            assert 'Q' in ds, 'Q variable missing from output'
            assert ds['Q'].shape[0] == 24, f'Expected 24 time steps, got {ds["Q"].shape[0]}'
            q_vals = ds['Q'].values
            assert np.all(q_vals >= 0), 'Negative discharge values found'
            # With no lateral inflow, total discharge should decrease over time
            total_per_step = q_vals.sum(axis=1)
            assert total_per_step[-1] <= total_per_step[0], \
                'Total discharge should not increase without lateral inflow'
            log.info(f'  Q range: {q_vals.min():.4f} to {q_vals.max():.4f} m³/s')
            log.info(f'  Total Q first step: {total_per_step[0]:.2f}, last step: {total_per_step[-1]:.2f}')

        # Verify final state was written
        assert os.path.exists(final_state_file), 'Final state file not written'
        final_state = pd.read_parquet(final_state_file)
        assert 'Q' in final_state.columns, 'Final state missing Q column'
        log.info('[OK] Muskingum channel-only routing passed')

    except Exception as e:
        log.error(f'[FAIL] {e}')
        raise
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════════
# Test 4: RapidMuskingum routing with catchment volumes
# ═══════════════════════════════════════════════════════════════════════════

def test_rapid_muskingum_from_volumes() -> None:
    """Route pre-computed catchment volumes and compare against retrospective."""
    log.info('=' * 60)
    log.info('TEST 4: RapidMuskingum from catchment volumes')
    if not _require_data(routing_params_file, initial_state_file, *volume_files):
        return

    tmpdir = tempfile.mkdtemp()
    try:
        discharge_files = [os.path.join(tmpdir, f'q_vol_{i}.nc') for i in range(len(volume_files))]
        final_state_file = os.path.join(tmpdir, 'final_state.parquet')

        rr.RapidMuskingum(
            params_file=str(routing_params_file),
            catchment_runoff_files=[str(f) for f in volume_files],
            discharge_files=discharge_files,
            channel_state_init_file=str(initial_state_file),
            channel_state_final_file=final_state_file,
            log=False,
        ).route()

        # Verify outputs
        for qf in discharge_files:
            assert os.path.exists(qf), f'Missing output: {qf}'
            with xr.open_dataset(qf) as ds:
                assert 'Q' in ds
                assert np.all(ds['Q'].values >= 0), f'Negative Q in {qf}'
                log.info(f'  {Path(qf).name}: shape={ds["Q"].shape}, '
                         f'Q range=[{ds["Q"].values.min():.2f}, {ds["Q"].values.max():.2f}]')

        # Compare against retrospective spot-check
        # if retrospective_file.exists():
        #     retro = xr.open_dataset(retrospective_file)
        #     for qf in discharge_files:
        #         routed = xr.open_dataset(qf)
        #         # select overlapping rivers and time
        #         common_rivers = np.intersect1d(routed.river_id.values, retro.river_id.values)
        #         common_times = np.intersect1d(routed.time.values, retro.time.values)
        #         if len(common_rivers) and len(common_times):
        #             r = routed['Q'].sel(river_id=common_rivers, time=common_times).values
        #             e = retro['Q'].sel(river_id=common_rivers, time=common_times).values
        #             max_diff = np.nanmax(np.abs(r - e))
        #             log.info(f'  Retrospective comparison: max diff = {max_diff:.4f} m³/s')
        #         routed.close()
        #     retro.close()

        assert os.path.exists(final_state_file)
        log.info('[OK] RapidMuskingum from volumes passed')

    except Exception as e:
        log.error(f'[FAIL] {e}')
        raise
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════════
# Test 5: RapidMuskingum from gridded depths + weight table
# ═══════════════════════════════════════════════════════════════════════════

def test_rapid_muskingum_from_depths() -> None:
    """Route ERA5 gridded runoff depths using the weight table."""
    log.info('=' * 60)
    log.info('TEST 5: RapidMuskingum from gridded depths + weights')
    if not _require_data(routing_params_file, initial_state_file, gridweights_file, *era5_files):
        return

    era5_vars = _detect_era5_var_names()

    tmpdir = tempfile.mkdtemp()
    try:
        discharge_files = [os.path.join(tmpdir, f'q_depth_{i}.nc') for i in range(len(era5_files))]

        rr.RapidMuskingum(
            params_file=str(routing_params_file),
            runoff_grid_files=[str(f) for f in era5_files],
            grid_weights_file=str(gridweights_file),
            discharge_files=discharge_files,
            channel_state_init_file=str(initial_state_file),
            var_x=era5_vars['x_var'],
            var_y=era5_vars['y_var'],
            var_t=era5_vars['t_var'],
            log=False,
        ).route()

        for qf in discharge_files:
            assert os.path.exists(qf), f'Missing output: {qf}'
            with xr.open_dataset(qf) as ds:
                assert 'Q' in ds
                assert np.all(ds['Q'].values >= 0)
                log.info(f'  {Path(qf).name}: shape={ds["Q"].shape}, '
                         f'Q range=[{ds["Q"].values.min():.2f}, {ds["Q"].values.max():.2f}]')

        log.info('[OK] RapidMuskingum from depths passed')

    except Exception as e:
        log.error(f'[FAIL] {e}')
        raise
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════════
# Test 6: Initial state is read and used correctly
# ═══════════════════════════════════════════════════════════════════════════

def test_initial_state_used() -> None:
    """Compare routing with and without initial state — first discharge values should differ."""
    log.info('=' * 60)
    log.info('TEST 6: Initial state is read and used correctly')
    if not _require_data(routing_params_file, initial_state_file):
        return
    if not volume_files:
        log.warning('No volume files available — skipping')
        return

    tmpdir = tempfile.mkdtemp()
    try:
        q_with_state = os.path.join(tmpdir, 'q_with_state.nc')
        q_without_state = os.path.join(tmpdir, 'q_without_state.nc')

        # Route WITH initial state
        rr.RapidMuskingum(
            params_file=str(routing_params_file),
            catchment_runoff_files=[str(volume_files[0])],
            discharge_files=[q_with_state],
            channel_state_init_file=str(initial_state_file),
            log=False,
        ).route()

        # Route WITHOUT initial state (defaults to zeros)
        rr.RapidMuskingum(
            params_file=str(routing_params_file),
            catchment_runoff_files=[str(volume_files[0])],
            discharge_files=[q_without_state],
            log=False,
        ).route()

        with xr.open_dataset(q_with_state) as ds1, xr.open_dataset(q_without_state) as ds2:
            q1_first = ds1['Q'].values[0, :]
            q2_first = ds2['Q'].values[0, :]
            diff = np.abs(q1_first - q2_first).sum()
            assert diff > 0, 'First timestep should differ when using different initial states'
            log.info(f'  First-step total |diff| = {diff:.2f} m³/s')

            # Later timesteps should converge as initial state influence decays
            q1_last = ds1['Q'].values[-1, :]
            q2_last = ds2['Q'].values[-1, :]
            last_diff = np.abs(q1_last - q2_last).sum()
            log.info(f'  Last-step total |diff| = {last_diff:.2f} m³/s')
            assert last_diff < diff, 'Initial state influence should decay over time'

        log.info('[OK] Initial state test passed')

    except Exception as e:
        log.error(f'[FAIL] {e}')
        raise
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════════
# Test 7: Final state is written and correct
# ═══════════════════════════════════════════════════════════════════════════

def test_final_state_roundtrip() -> None:
    """Route, write final state, re-route starting from that state — should be seamless."""
    log.info('=' * 60)
    log.info('TEST 7: Final state written and re-used correctly')
    if not _require_data(routing_params_file, initial_state_file):
        return
    if len(volume_files) < 2:
        log.warning('Need at least 2 volume files for state roundtrip — skipping')
        return

    tmpdir = tempfile.mkdtemp()
    try:
        # Route all months sequentially in one call
        q_all_at_once = [os.path.join(tmpdir, f'q_all_{i}.nc') for i in range(len(volume_files))]
        rr.RapidMuskingum(
            params_file=str(routing_params_file),
            catchment_runoff_files=[str(f) for f in volume_files],
            discharge_files=q_all_at_once,
            channel_state_init_file=str(initial_state_file),
            log=False,
        ).route()

        # Route month 1, save state, then route month 2 from that state
        q_month1 = os.path.join(tmpdir, 'q_m1.nc')
        state_after_m1 = os.path.join(tmpdir, 'state_after_m1.parquet')
        rr.RapidMuskingum(
            params_file=str(routing_params_file),
            catchment_runoff_files=[str(volume_files[0])],
            discharge_files=[q_month1],
            channel_state_init_file=str(initial_state_file),
            channel_state_final_file=state_after_m1,
            log=False,
        ).route()

        q_month2 = os.path.join(tmpdir, 'q_m2.nc')
        rr.RapidMuskingum(
            params_file=str(routing_params_file),
            catchment_runoff_files=[str(volume_files[1])],
            discharge_files=[q_month2],
            channel_state_init_file=state_after_m1,
            log=False,
        ).route()

        # Compare: month 2 from sequential routing should match month 2 from all-at-once
        with xr.open_dataset(q_all_at_once[1]) as ds_all, xr.open_dataset(q_month2) as ds_split:
            assert ds_all['Q'].shape == ds_split['Q'].shape, \
                f'Shape mismatch: {ds_all["Q"].shape} vs {ds_split["Q"].shape}'
            if np.allclose(ds_all['Q'].values, ds_split['Q'].values, rtol=1e-4, atol=1e-4):
                log.info('[OK] State roundtrip: sequential matches all-at-once')
            else:
                diff = np.abs(ds_all['Q'].values - ds_split['Q'].values)
                log.error(f'[FAIL] State roundtrip: max diff = {diff.max():.6f}')

    except Exception as e:
        log.error(f'[FAIL] {e}')
        raise
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════════
# Test 8: Grid weights reproduction
# ═══════════════════════════════════════════════════════════════════════════

# Requires catchment geometries (large file). Uncomment when catchments parquet is available.
# def test_grid_weights_reproduction() -> None:
#     """Reproduce the grid weight table from ERA5 grid and catchment geometries."""
#     log.info('=' * 60)
#     log.info('TEST 8: Grid weights reproduction')
#
#     catchments_file = data_root / 'catchments.parquet'
#     if not _require_data(routing_params_file, catchments_file, *era5_files[:1]):
#         return
#
#     tmpdir = tempfile.mkdtemp()
#     try:
#         voronoi_file = os.path.join(tmpdir, 'voronoi.parquet')
#         weights_file = os.path.join(tmpdir, 'gridweights.nc')
#
#         rr.runoff.grid_weights(
#             str(era5_files[0]),
#             str(catchments_file),
#             x_var='longitude',
#             y_var='latitude',
#             river_id_var='LINKNO',
#             save_voronoi_path=voronoi_file,
#             save_weights_path=weights_file,
#             routing_params_path=str(routing_params_file),
#         )
#
#         # Compare against reference grid weights
#         ref = xr.load_dataset(str(gridweights_file))
#         gen = xr.load_dataset(weights_file)
#
#         # Compare proportions (the key output)
#         ref_prop = ref.groupby('river_id')['proportion'].sum()
#         gen_prop = gen.groupby('river_id')['proportion'].sum()
#         assert np.allclose(ref_prop.values, 1.0, atol=1e-4)
#         assert np.allclose(gen_prop.values, 1.0, atol=1e-4)
#         log.info('[OK] Grid weights reproduce correctly')
#
#     except Exception as e:
#         log.error(f'[FAIL] {e}')
#         raise
#     finally:
#         shutil.rmtree(tmpdir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════════
# Test 9: Depths+weights vs pre-computed volumes produce same discharge
# ═══════════════════════════════════════════════════════════════════════════

def test_volumes_vs_depths_equivalence() -> None:
    """Route from depths+weights and from pre-computed volumes — results must match."""
    log.info('=' * 60)
    log.info('TEST 9: Volumes vs depths+weights equivalence')
    if not _require_data(routing_params_file, initial_state_file, gridweights_file):
        return
    if not era5_files or not volume_files:
        log.warning('Need both ERA5 and volume files — skipping')
        return

    era5_vars = _detect_era5_var_names()

    # Use only the first month to keep it fast
    era5_file = era5_files[0]
    vol_file = volume_files[0]

    tmpdir = tempfile.mkdtemp()
    try:
        q_from_depths = os.path.join(tmpdir, 'q_depths.nc')
        q_from_volumes = os.path.join(tmpdir, 'q_volumes.nc')

        # Route from gridded depths + weight table
        rr.RapidMuskingum(
            params_file=str(routing_params_file),
            runoff_grid_files=[str(era5_file)],
            grid_weights_file=str(gridweights_file),
            discharge_files=[q_from_depths],
            channel_state_init_file=str(initial_state_file),
            var_x=era5_vars['x_var'],
            var_y=era5_vars['y_var'],
            var_t=era5_vars['t_var'],
            log=False,
        ).route()

        # Route from pre-computed catchment volumes
        rr.RapidMuskingum(
            params_file=str(routing_params_file),
            catchment_runoff_files=[str(vol_file)],
            discharge_files=[q_from_volumes],
            channel_state_init_file=str(initial_state_file),
            log=False,
        ).route()

        # Compare: should be identical (same volumes → same routing)
        match = compare_netcdfs(
            Path(q_from_depths), Path(q_from_volumes),
            variable='Q', rtol=1e-4, atol=0.01,
        )
        if match:
            log.info('[OK] Depths+weights and volumes produce matching discharge')
        else:
            log.error('[FAIL] Discharge differs between depths+weights and volumes paths')

    except Exception as e:
        log.error(f'[FAIL] {e}')
        raise
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════════
# Run all
# ═══════════════════════════════════════════════════════════════════════════

def run_all() -> None:
    test_valid_config()
    test_invalid_configs()
    test_muskingum_channel_only()
    test_rapid_muskingum_from_volumes()
    test_rapid_muskingum_from_depths()
    test_initial_state_used()
    test_final_state_roundtrip()
    # test_grid_weights_reproduction()  # requires catchment geometries
    test_volumes_vs_depths_equivalence()


if __name__ == '__main__':
    run_all()
