"""Shared fixtures and helpers for the river-route test suite."""
import os
from dataclasses import dataclass
from glob import glob
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
DATA_DIR = TESTS_DIR / 'data'

# provided by the zip downloaded from s3
ERA5_DIR = DATA_DIR / 'era5'
DISCHARGE_DIR = DATA_DIR / 'discharge'
CATCHMENT_RUNOFF_DIR = DATA_DIR / 'catchment-runoff'
# obtained from s3 based on which vpus are given as having solutions
PARAMS_DIR = DATA_DIR / 'routing-configs'
HYDROGRAPHY_DIR = DATA_DIR / 'hydrography'

ERA5_FILES = sorted(glob(str(ERA5_DIR / 'era5*.nc')))
ERA5_KWARGS = dict(var_y='latitude', var_x='longitude', var_t='valid_time')


# ── VPU dataclass ────────────────────────────────────────────────────────────

@dataclass
class VPUData:
    """Per-VPU test data paths. All fields are required and paths should exist"""
    number: int
    params_file: Path
    grid_weights_file: Path
    discharge_dir: Path
    catchment_volumes_dir: Path
    hydrography_dir: Path
    discharge_files: list[str] = None  # populated in __post_init__
    catchment_files: list[str] = None  # populated in __post_init__

    def __str__(self) -> str:
        return f'vpu={self.number}'

    def __post_init__(self):
        self.discharge_files = list(sorted(glob(str(self.discharge_dir / 'discharge*.nc'))))
        self.catchment_files = list(sorted(glob(str(self.catchment_volumes_dir / 'inflowing_*.nc'))))

    def valid(self) -> bool:
        """Check that all paths exist and file lists are non-empty."""
        path_fields = ('params_file', 'grid_weights_file', 'discharge_dir', 'catchment_volumes_dir', 'hydrography_dir')
        for field in path_fields:
            value = getattr(self, field)
            if not value.exists():
                print(f'VPU {self.number}: {field} path does not exist: {value}')
                return False

        if not self.discharge_files:
            print(f'VPU {self.number}: no discharge files found in {self.discharge_dir}')
            return False
        if not self.catchment_files:
            print(f'VPU {self.number}: no catchment files found in {self.catchment_volumes_dir}')
            return False

        return True


def find_test_units() -> list[VPUData]:
    """find the vpus of routing configs, then prepare them for river-route version 2 format tests"""
    # vpus = sorted([os.path.basename(glob(str(PARAMS_DIR / 'vpu=*'))))
    vpus = [os.path.basename(path) for path in glob(str(PARAMS_DIR / 'vpu=*'))]
    testable_sets = []
    for vpu in sorted(vpus):
        testable_sets.append(
            VPUData(
                number=int(vpu.split('=')[1]),
                params_file=PARAMS_DIR / vpu / 'params.parquet',
                grid_weights_file=PARAMS_DIR / vpu / 'grid_weights.parquet',
                discharge_dir=DISCHARGE_DIR / vpu,
                catchment_volumes_dir=CATCHMENT_RUNOFF_DIR / vpu,
                hydrography_dir=HYDROGRAPHY_DIR / vpu,
            )
        )
    testable_sets = [vpu for vpu in testable_sets if vpu.valid()]
    return testable_sets


# Discover once at import time
TEST_CASES = find_test_units()

# ── pytest hook: auto-parametrize tests with a `vpu` parameter ──────────────

def pytest_generate_tests(metafunc):
    if 'vpu' in metafunc.fixturenames:
        metafunc.parametrize('vpu', TEST_CASES, ids=[str(v) for v in TEST_CASES])
