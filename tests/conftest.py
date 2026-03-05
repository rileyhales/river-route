"""Shared fixtures and helpers for the river-route test suite."""
import os
from dataclasses import dataclass
from glob import glob
from pathlib import Path

import pandas as pd

TESTS_DIR = Path(__file__).resolve().parent
DATA_DIR = TESTS_DIR / 'data'

# provided by the zip downloaded from s3
ERA5_DIR = DATA_DIR / 'era5'
DISCHARGE_DIR = DATA_DIR / 'discharge'
QLATERAL_DIR = DATA_DIR / 'qlateral'
# obtained from s3 based on which vpus are given as having solutions
CONFIGS_DIR = DATA_DIR / 'routing-configs'
HYDROGRAPHY_DIR = DATA_DIR / 'hydrography'

ERA5_FILES = sorted(glob(str(ERA5_DIR / 'era5*.nc')))
ERA5_KWARGS = dict(var_y='latitude', var_x='longitude', var_t='valid_time')


@dataclass
class RFSv2ConfigsData:
    """Per-VPU test data paths. All fields are required and paths should exist"""
    number: int
    configs_dir: Path
    hydrography_dir: Path
    discharge_dir: Path
    qlateral_dir: Path
    discharge_files: list[str] = None  # populated in __post_init__
    qlateral_files: list[str] = None  # populated in __post_init__

    def __str__(self) -> str:
        return f'vpu={self.number}'

    def __post_init__(self) -> None:
        self.discharge_files = list(sorted(glob(str(self.discharge_dir / 'discharge*.nc'))))
        self.qlateral_files = list(sorted(glob(str(self.qlateral_dir / 'qlateral_*.nc'))))

        # convert the rfs v2 configs to river-route v2 formats
        pdf = pd.read_parquet(self.rr1_params_file)
        cdf = pd.read_parquet(self.rr1_connectivity_file)
        (
            pdf
            .merge(cdf, left_on='river_id', right_on='river_id', how='left')
            .rename(columns={'ds_river_id': 'downstream_river_id'})
            [['river_id', 'downstream_river_id', 'k', 'x']]
            .to_parquet(self.rr2_params_file, index=False)
        )
        return

    @property
    def rr2_params_file(self) -> Path:
        return self.configs_dir / 'rr2_params.parquet'

    @property
    def rr1_params_file(self) -> Path:
        return self.configs_dir / 'routing_parameters.parquet'

    @property
    def rr1_connectivity_file(self) -> Path:
        return self.configs_dir / 'connectivity.parquet'

    @property
    def rr1_grid_weights_file(self) -> Path:
        return self.configs_dir / f'gridweights_ERA5_vpu={self.number}.nc'

    @property
    def catchments(self) -> Path:
        return self.hydrography_dir / f'catchments_{self.number}.parquet'

    @property
    def streams(self) -> Path:
        return self.hydrography_dir / f'streams_{self.number}.parquet'

    def valid(self) -> bool:
        """Check that all paths exist and file lists are non-empty."""
        path_fields = ('configs_dir', 'discharge_dir', 'qlateral_dir', 'hydrography_dir')
        for field in path_fields:
            value = getattr(self, field)
            if not value.exists():
                print(f'VPU {self.number}: {field} path does not exist: {value}')
                return False

        if not self.discharge_files:
            print(f'VPU {self.number}: no discharge files found in {self.discharge_dir}')
            return False
        if not self.qlateral_files:
            print(f'VPU {self.number}: no qlateral files found in {self.qlateral_dir}')
            return False

        if not self.rr1_params_file.exists():
            print(f'VPU {self.number}: params file does not exist')
            return False
        if not self.rr1_connectivity_file.exists():
            print(f'VPU {self.number}: connectivity file does not exist')
            return False
        if not self.rr1_grid_weights_file.exists():
            print(f'VPU {self.number}: grid weights file does not exist')
            return False
        if not self.catchments.exists():
            print(f'VPU {self.number}: catchments file does not exist')
            return False
        if not self.streams.exists():
            print(f'VPU {self.number}: streams file does not exist')
            return False

        return True


def find_test_units() -> list[RFSv2ConfigsData]:
    """find the vpus of routing configs, then prepare them for river-route version 2 format tests"""
    # vpus = sorted([os.path.basename(glob(str(CONFIGS_DIR / 'vpu=*'))))
    vpus = [os.path.basename(path) for path in glob(str(CONFIGS_DIR / 'vpu=*'))]
    testable_sets = []
    for vpu in sorted(vpus):
        testable_sets.append(
            RFSv2ConfigsData(
                number=int(vpu.split('=')[1]),
                configs_dir=CONFIGS_DIR / vpu,
                discharge_dir=DISCHARGE_DIR / vpu,
                qlateral_dir=QLATERAL_DIR / vpu,
                hydrography_dir=HYDROGRAPHY_DIR / vpu,
            )
        )
    testable_sets = [vpu for vpu in testable_sets if vpu.valid()]
    return testable_sets


# ── pytest hook: auto-parametrize tests with a `vpu` parameter ──────────────
def pytest_generate_tests(metafunc):
    if 'vpu' in metafunc.fixturenames:
        metafunc.parametrize('vpu', TEST_CASES, ids=[str(v) for v in TEST_CASES])


TEST_CASES = find_test_units()
