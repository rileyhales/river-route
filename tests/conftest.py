"""Shared fixtures and helpers for the river-route test suite."""

import os
from dataclasses import dataclass
from glob import glob
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd
import pytest

TESTS_DIR = Path(__file__).resolve().parent
DATA_DIR = TESTS_DIR / 'data'

# provided by the zip downloaded from s3
ERA5_DIR = DATA_DIR / 'era5'
DISCHARGE_DIR = DATA_DIR / 'discharge'
VLATERAL_DIR = DATA_DIR / 'vlateral'
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
    vlateral_dir: Path
    discharge_files: list[str] = None  # populated in __post_init__
    vlateral_files: list[str] = None  # populated in __post_init__

    def __str__(self) -> str:
        return f'vpu={self.number}'

    def __post_init__(self) -> None:
        self.discharge_files = list(sorted(glob(str(self.discharge_dir / 'discharge*.nc'))))
        self.vlateral_files = list(sorted(glob(str(self.vlateral_dir / 'vlateral_*.nc'))))

    def prepare(self) -> None:
        """Convert rfs v2 configs to river-route v3 formats. Call after valid() confirms files exist."""
        pdf = pd.read_parquet(self.rr1_params_file)
        cdf = pd.read_parquet(self.rr1_connectivity_file)
        (
            pdf.merge(cdf, left_on='river_id', right_on='river_id', how='left')
            .rename(columns={'ds_river_id': 'next_river_id'})[['river_id', 'next_river_id', 'k', 'x']]
            .to_parquet(self.rr2_params_file, index=False)
        )

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
    def grid_weights_file(self) -> Path:
        return self.configs_dir / f'gridweights_ERA5_vpu={self.number}.nc'

    @property
    def catchments(self) -> Path:
        return self.hydrography_dir / f'catchments_{self.number}.parquet'

    @property
    def streams(self) -> Path:
        return self.hydrography_dir / f'streams_{self.number}.gpkg'

    def valid(self) -> bool:
        """Check that all paths exist and file lists are non-empty."""
        path_fields = ('configs_dir', 'discharge_dir', 'vlateral_dir', 'hydrography_dir')
        for field in path_fields:
            value = getattr(self, field)
            if not value.exists():
                print(f'VPU {self.number}: {field} path does not exist: {value}')
                return False

        if not self.discharge_files:
            print(f'VPU {self.number}: no discharge files found in {self.discharge_dir}')
            return False
        if not self.vlateral_files:
            print(f'VPU {self.number}: no vlateral files found in {self.vlateral_dir}')
            return False

        if not self.rr1_params_file.exists():
            print(f'VPU {self.number}: params file does not exist')
            return False
        if not self.rr1_connectivity_file.exists():
            print(f'VPU {self.number}: connectivity file does not exist')
            return False
        if not self.grid_weights_file.exists():
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
    vpus = [os.path.basename(path) for path in glob(str(DISCHARGE_DIR / 'vpu=*'))]
    testable_sets = []
    for vpu in sorted(vpus):
        testable_sets.append(
            RFSv2ConfigsData(
                number=int(vpu.split('=')[1]),
                configs_dir=CONFIGS_DIR / vpu,
                discharge_dir=DISCHARGE_DIR / vpu,
                vlateral_dir=VLATERAL_DIR / vpu,
                hydrography_dir=HYDROGRAPHY_DIR / vpu,
            )
        )
    testable_sets = [vpu for vpu in testable_sets if vpu.valid()]
    for vpu in testable_sets:
        vpu.prepare()
    return testable_sets


# ── synthetic networks: small, self-contained fixtures that need no downloaded data ──


@dataclass
class SyntheticNetwork:
    """A small river network written to disk, used to exercise routing without the downloaded test data."""

    directory: Path
    n_rivers: int
    k: float
    x: float
    n_steps: int
    dt_runoff: int
    params_file: Path
    state_file: Path
    vlateral_file: Path
    inflow_volume: float

    def path(self, name: str) -> str:
        return str(self.directory / name)


def write_params(path: Path, n_rivers: int = 5, k: float = 3600.0, x: float = 0.2, **columns) -> pd.DataFrame:
    """Write a params file for a single chain of rivers, each flowing into the next, sorted upstream to down."""
    df = pd.DataFrame(
        {
            'river_id': np.arange(1, n_rivers + 1, dtype=np.int64),
            'next_river_id': np.append(np.arange(2, n_rivers + 1), -1).astype(np.int64),
            'k': np.full(n_rivers, k, dtype=np.float64),
            'x': np.full(n_rivers, x, dtype=np.float64),
            **columns,
        }
    )
    df.to_parquet(path, index=False)
    return df


def write_vlateral(path: Path, volumes: np.ndarray, river_ids: np.ndarray, dt: int = 3600) -> None:
    """Write a lateral inflow netCDF with a CF encoded time axis, matching what Runoff produces."""
    volumes = np.asarray(volumes, dtype=np.float32)
    with nc.Dataset(str(path), mode='w', format='NETCDF4') as ds:
        ds.createDimension('time', volumes.shape[0])
        ds.createDimension('river_id', volumes.shape[1])
        time_var = ds.createVariable('time', 'f8', ('time',))
        time_var.units = 'seconds since 2000-01-01 00:00:00'
        time_var[:] = np.arange(volumes.shape[0]) * dt
        id_var = ds.createVariable('river_id', 'i8', ('river_id',))
        id_var[:] = river_ids
        vlateral = ds.createVariable('vlateral', 'f4', ('time', 'river_id'))
        vlateral[:] = volumes
        vlateral.units = 'm3'
    return


def build_network(
    directory: Path,
    n_rivers: int = 5,
    k: float = 3600.0,
    x: float = 0.2,
    n_steps: int = 240,
    dt_runoff: int = 3600,
    pulse_steps: int = 3,
    pulse_rate: float = 10.0,
) -> SyntheticNetwork:
    """
    Build a chain network plus a lateral inflow file holding a short pulse into the headwater river.

    The default k and x are stable for the default dt_runoff: 2*k*x = 1440 <= 3600 <= 2*k*(1-x) = 5760.
    The window is long enough for the pulse to drain completely so that mass balance can be checked.
    """
    directory.mkdir(parents=True, exist_ok=True)
    params_file = directory / 'params.parquet'
    write_params(params_file, n_rivers=n_rivers, k=k, x=x)

    state_file = directory / 'state_zero.parquet'
    pd.DataFrame({'Q': np.zeros(n_rivers)}).to_parquet(state_file, index=False)

    volumes = np.zeros((n_steps, n_rivers), dtype=np.float32)
    volumes[:pulse_steps, 0] = pulse_rate * dt_runoff  # a constant rate held for pulse_steps, as a volume per step
    vlateral_file = directory / 'vlateral.nc'
    write_vlateral(vlateral_file, volumes, np.arange(1, n_rivers + 1), dt=dt_runoff)

    return SyntheticNetwork(
        directory=directory,
        n_rivers=n_rivers,
        k=k,
        x=x,
        n_steps=n_steps,
        dt_runoff=dt_runoff,
        params_file=params_file,
        state_file=state_file,
        vlateral_file=vlateral_file,
        inflow_volume=float(volumes.sum()),
    )


@pytest.fixture
def network(tmp_path: Path) -> SyntheticNetwork:
    """A stable 5 river chain with a 3 hour pulse of lateral inflow and a 240 hour routing window."""
    return build_network(tmp_path / 'network')


# ── pytest hook: auto-parametrize tests with a `vpu` parameter ──────────────
def pytest_generate_tests(metafunc):
    if 'vpu' in metafunc.fixturenames:
        metafunc.parametrize('vpu', TEST_CASES, ids=[str(v) for v in TEST_CASES])


TEST_CASES = find_test_units()
