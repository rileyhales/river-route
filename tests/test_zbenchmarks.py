import os
import shutil
import tempfile
import time
import tracemalloc

import pytest
import river_route as rr
from conftest import ERA5_FILES, ERA5_KWARGS, RFSv2ConfigsData
from river_route.runoff import runoff_to_qlateral


def _fmt_mem(bytes_val):
    if bytes_val >= 1 << 30:
        return f'{bytes_val / (1 << 30):.2f} GB'
    if bytes_val >= 1 << 20:
        return f'{bytes_val / (1 << 20):.2f} MB'
    return f'{bytes_val / (1 << 10):.2f} KB'


def _report(label, elapsed, peak_mem):
    print(f'\n  {label}')
    print(f'    elapsed:    {elapsed:.3f} s')
    print(f'    peak memory: {_fmt_mem(peak_mem)}')


# ── Routing benchmarks ──────────────────────────────────────────────────────

def test_bench_route_from_qlateral(vpu: RFSv2ConfigsData):
    """Benchmark: route 1 month from pre-computed qlateral."""
    tmpdir = tempfile.mkdtemp()
    try:
        discharge_file = os.path.join(tmpdir, 'q.nc')

        tracemalloc.start()
        t0 = time.perf_counter()

        rr.RapidMuskingum(
            params_file=str(vpu.rr2_params_file),
            qlateral_files=[vpu.qlateral_files[0]],
            discharge_files=[discharge_file],
            log=False,
            progress_bar=False,
        ).route()

        elapsed = time.perf_counter() - t0
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        _report(f'route_from_qlateral vpu={vpu.number}', elapsed, peak)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_bench_route_from_depths(vpu: RFSv2ConfigsData):
    """Benchmark: route 1 month of ERA5 through grid weights."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    tmpdir = tempfile.mkdtemp()
    try:
        discharge_file = os.path.join(tmpdir, 'q.nc')

        tracemalloc.start()
        t0 = time.perf_counter()

        rr.RapidMuskingum(
            params_file=str(vpu.rr2_params_file),
            grid_weights_file=str(vpu.grid_weights_file),
            grid_runoff_files=[ERA5_FILES[0]],
            discharge_files=[discharge_file],
            log=False,
            progress_bar=False,
            **ERA5_KWARGS,
        ).route()

        elapsed = time.perf_counter() - t0
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        _report(f'route_from_depths vpu={vpu.number}', elapsed, peak)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_bench_runoff_to_qlateral(vpu: RFSv2ConfigsData):
    """Benchmark: convert 1 month ERA5 gridded runoff to qlateral."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    tracemalloc.start()
    t0 = time.perf_counter()

    runoff_to_qlateral(
        ERA5_FILES[0],
        grid_weights_file=str(vpu.grid_weights_file),
        var_runoff='ro',
        var_x='longitude',
        var_y='latitude',
        var_t='valid_time',
        as_volumes=True,
    )

    elapsed = time.perf_counter() - t0
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    _report(f'runoff_to_qlateral vpu={vpu.number}', elapsed, peak)


# ── End-to-end benchmark ────────────────────────────────────────────────────

def test_bench_end_to_end(vpu: RFSv2ConfigsData):
    """Benchmark: full pipeline - ERA5 grid runoff -> qlateral -> routed discharge."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')

    tmpdir = tempfile.mkdtemp()
    try:
        discharge_file = os.path.join(tmpdir, 'q.nc')

        tracemalloc.start()
        t0 = time.perf_counter()

        ds_ql = runoff_to_qlateral(
            ERA5_FILES[0],
            grid_weights_file=str(vpu.grid_weights_file),
            var_runoff='ro',
            var_x='longitude',
            var_y='latitude',
            var_t='valid_time',
            as_volumes=True,
        )

        qlateral_file = os.path.join(tmpdir, 'ql.nc')
        ds_ql.to_netcdf(qlateral_file)

        rr.RapidMuskingum(
            params_file=str(vpu.rr2_params_file),
            qlateral_files=[qlateral_file],
            discharge_files=[discharge_file],
            log=False,
            progress_bar=False,
        ).route()

        elapsed = time.perf_counter() - t0
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        _report(f'end_to_end vpu={vpu.number}', elapsed, peak)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
