"""Tests for river_route.Configs validation."""
import os
import shutil
import tempfile

import pytest

import river_route as rr
from conftest import ERA5_FILES, ERA5_KWARGS, VPUData


def test_valid_muskingum_config(vpu: VPUData):
    """Construct a Muskingum config without error."""
    tmpdir = tempfile.mkdtemp()
    try:
        rr.Configs(
            params_file=str(vpu.params_file),
            discharge_files=[os.path.join(tmpdir, 'q.nc')],
            dt_routing=3600,
            dt_total=86400,
            log=False,
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_valid_rapid_muskingum_config(vpu: VPUData):
    """Construct a RapidMuskingum config without error."""
    if not ERA5_FILES:
        pytest.skip('Missing ERA5 files')
    tmpdir = tempfile.mkdtemp()
    try:
        rr.Configs(
            params_file=str(vpu.params_file),
            grid_runoff_files=[ERA5_FILES[0]],
            grid_weights_file=str(vpu.grid_weights_file),
            discharge_dir=tmpdir,
            log=False,
            **ERA5_KWARGS,
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_invalid_config_missing_params_file():
    tmpdir = tempfile.mkdtemp()
    try:
        with pytest.raises(FileNotFoundError):
            rr.Configs(
                params_file='/nonexistent/params.parquet',
                discharge_files=[os.path.join(tmpdir, 'q.nc')],
                dt_routing=3600,
                dt_total=86400,
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_invalid_config_bad_discharge_dir(vpu: VPUData):
    with pytest.raises(NotADirectoryError):
        rr.Configs(
            params_file=str(vpu.params_file),
            discharge_files=['/nonexistent_dir/q.nc'],
            dt_routing=3600,
            dt_total=86400,
        )


def test_invalid_config_bad_processing_mode(vpu: VPUData):
    tmpdir = tempfile.mkdtemp()
    try:
        with pytest.raises(ValueError):
            rr.Configs(
                params_file=str(vpu.params_file),
                discharge_files=[os.path.join(tmpdir, 'q.nc')],
                dt_routing=3600,
                dt_total=86400,
                runoff_processing_mode='invalid_mode',
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
