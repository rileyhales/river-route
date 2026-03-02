"""Shared fixtures and helpers for the river-route test suite."""
from dataclasses import dataclass
from glob import glob
from pathlib import Path

import pytest

# ── Path constants ───────────────────────────────────────────────────────────
TESTS_DIR = Path(__file__).resolve().parent
DATA_DIR = TESTS_DIR / 'data'

# ERA5 is global gridded data — shared across VPUs
ERA5_DIR = DATA_DIR / 'era5'
ERA5_KWARGS = dict(var_y='latitude', var_x='longitude', var_t='valid_time')

# Exclude data helper scripts and legacy solutions from test collection
collect_ignore_glob = ['data/scripts/*', 'data/solutions/*']


# ── VPU dataclass ────────────────────────────────────────────────────────────

@dataclass
class VPUData:
    """Per-VPU test data paths. Fields are None when data is not available."""
    vpu_id: str
    params_file: Path | None = None
    grid_weights_file: Path | None = None
    discharge_dir: Path | None = None
    catchment_volumes_dir: Path | None = None
    hydrography_dir: Path | None = None

    def __str__(self) -> str:
        return f'vpu={self.vpu_id}'


def _resolve_or_none(path: Path) -> Path | None:
    """Return the path if it exists, otherwise None."""
    return path if path.exists() else None


def _first_glob(directory: Path, pattern: str) -> Path | None:
    """Return the first file matching *pattern* inside *directory*, or None."""
    matches = sorted(directory.glob(pattern))
    return matches[0] if matches else None


def discover_vpus() -> list[VPUData]:
    """Scan data directories for vpu=* subdirs and build a VPUData for each."""
    vpu_ids: set[str] = set()

    for parent in ('routing-configs', 'discharge', 'catchment-volumes', 'hydrography'):
        parent_dir = DATA_DIR / parent
        if parent_dir.is_dir():
            for child in parent_dir.iterdir():
                if child.is_dir() and child.name.startswith('vpu='):
                    vpu_ids.add(child.name.removeprefix('vpu='))

    vpus = []
    for vpu_id in sorted(vpu_ids):
        rc_dir = DATA_DIR / 'routing-configs' / f'vpu={vpu_id}'
        vpus.append(VPUData(
            vpu_id=vpu_id,
            params_file=_first_glob(rc_dir, '*params*.parquet') if rc_dir.is_dir() else None,
            grid_weights_file=_first_glob(rc_dir, 'gridweights*.nc') if rc_dir.is_dir() else None,
            discharge_dir=_resolve_or_none(DATA_DIR / 'discharge' / f'vpu={vpu_id}'),
            catchment_volumes_dir=_resolve_or_none(DATA_DIR / 'catchment-volumes' / f'vpu={vpu_id}'),
            hydrography_dir=_resolve_or_none(DATA_DIR / 'hydrography' / f'vpu={vpu_id}'),
        ))

    return vpus


# Discover once at import time
ALL_VPUS = discover_vpus()


# ── pytest hook: auto-parametrize tests with a `vpu` parameter ──────────────

def pytest_generate_tests(metafunc):
    if 'vpu' in metafunc.fixturenames:
        metafunc.parametrize('vpu', ALL_VPUS, ids=[str(v) for v in ALL_VPUS])


# ── Helpers ──────────────────────────────────────────────────────────────────

def era5_files() -> list[str]:
    return sorted(glob(str(ERA5_DIR / 'era5_1940*.nc')))


def known_discharge_files(vpu: VPUData) -> list[str]:
    if vpu.discharge_dir is None:
        return []
    return sorted(glob(str(vpu.discharge_dir / '*era5_1940*.nc')))


def known_catchment_volume_files(vpu: VPUData) -> list[str]:
    if vpu.catchment_volumes_dir is None:
        return []
    return sorted(glob(str(vpu.catchment_volumes_dir / 'volumes_1940*.nc')))


def skip_if_vpu_missing(vpu: VPUData, *fields: str) -> None:
    """Skip the test if any named field on the VPU dataclass is None or doesn't exist."""
    for field in fields:
        value = getattr(vpu, field)
        if value is None:
            pytest.skip(f'VPU {vpu.vpu_id}: {field} not available')
        if isinstance(value, Path) and not value.exists():
            pytest.skip(f'VPU {vpu.vpu_id}: {field} path does not exist: {value}')
