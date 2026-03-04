"""Tests for river_route.uhkernels — SCS kernel builders and UnitHydrograph convolution."""
import os
import shutil
import tempfile

import numpy as np
import pandas as pd
import pytest

from river_route.uhkernels import SCSTriangular, SCSCurvilinear, UnitHydrograph

# ═════════════════════════════════════════════════════════════════════════════
# SCS Kernel Builders — shared tests for both triangular and curvilinear
# ═════════════════════════════════════════════════════════════════════════════

KERNEL_CLASSES = [SCSTriangular, SCSCurvilinear]


@pytest.fixture(params=KERNEL_CLASSES, ids=lambda cls: cls.__name__)
def kernel_cls(request):
    return request.param


def test_volume_conservation(kernel_cls):
    """sum(kernel[:, j] * tr) should equal area[j] for each basin."""
    tr = 3600.0
    tc = np.array([7200.0, 14400.0, 3600.0])
    area = np.array([1e6, 5e6, 2e5])
    uh = kernel_cls(tr, tc, area)
    integrated = uh.kernel.sum(axis=0) * tr
    np.testing.assert_allclose(
        integrated, area, rtol=1e-6,
        err_msg=f'{kernel_cls.__name__}: kernel does not conserve volume',
    )
    return


def test_kernel_shape(kernel_cls):
    """Kernel should be 2D with shape (n_steps, n_basins)."""
    tr = 3600.0
    tc = np.array([7200.0, 14400.0])
    area = np.array([1e6, 5e6])
    uh = kernel_cls(tr, tc, area)
    assert uh.kernel.ndim == 2
    assert uh.kernel.shape[1] == 2
    assert uh.kernel.shape[0] >= 1
    return


# ═════════════════════════════════════════════════════════════════════════════
# UnitHydrograph — convolution correctness
# ═════════════════════════════════════════════════════════════════════════════

def _make_kernel_file(kernel: np.ndarray, tmpdir: str) -> str:
    """Save a kernel array as parquet and return the path."""
    path = os.path.join(tmpdir, 'kernel.parquet')
    pd.DataFrame(kernel.T).to_parquet(path)
    return path


def test_convolve_vs_convolve_incrementally():
    """Full-pass convolve and step-by-step convolve_incrementally should produce identical results."""
    n_basins = 4
    n_ks = 3
    n_steps = 10
    np.random.seed(123)
    kernel = np.random.rand(n_ks, n_basins)
    lateral = np.random.rand(n_steps, n_basins)

    tmpdir = tempfile.mkdtemp()
    try:
        path = _make_kernel_file(kernel, tmpdir)

        # Full pass
        uh_full = UnitHydrograph(path)
        result_full = uh_full.convolve(lateral)

        # Incremental
        uh_inc = UnitHydrograph(path)
        result_inc = np.zeros_like(result_full)
        for t in range(n_steps):
            result_inc[t] = uh_inc.convolve_incrementally(lateral[t])

        np.testing.assert_allclose(result_full, result_inc, rtol=1e-12)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return


def test_convolve_impulse_response():
    """A unit impulse through a known kernel should reproduce the kernel column."""
    n_basins = 2
    kernel = np.array([[1.0, 0.5], [0.5, 0.3], [0.0, 0.2]])  # (3, 2)
    lateral = np.zeros((5, n_basins))
    lateral[0, :] = 1.0  # impulse at t=0

    tmpdir = tempfile.mkdtemp()
    try:
        path = _make_kernel_file(kernel, tmpdir)
        uh = UnitHydrograph(path)
        result = uh.convolve(lateral)
        # First n_ks rows should equal the kernel values
        np.testing.assert_allclose(result[:3], kernel, rtol=1e-12)
        # Remaining rows should be zero
        np.testing.assert_allclose(result[3:], 0.0, atol=1e-15)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return
