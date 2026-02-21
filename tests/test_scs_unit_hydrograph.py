import numpy as np

from river_route.transformers import SCSUnitHydrograph


def test_kernel_shape_changes_with_dt() -> None:
    tc = np.array([7200.0], dtype=np.float64)
    area = np.array([1000.0], dtype=np.float64)
    coarse = SCSUnitHydrograph(tc=tc, area=area, dt=3600.0)
    fine = SCSUnitHydrograph(tc=tc, area=area, dt=1800.0)
    assert coarse.kernel.shape[1] == 1
    assert fine.kernel.shape[1] == 1
    # tb = 2.67 * (0.6*tc + dt/2)
    expected_coarse = int(np.ceil((2.67 * (0.6 * tc[0] + 3600.0 / 2.0)) / 3600.0))
    expected_fine = int(np.ceil((2.67 * (0.6 * tc[0] + 1800.0 / 2.0)) / 1800.0))
    assert coarse.kernel.shape[0] == expected_coarse
    assert fine.kernel.shape[0] == expected_fine


def test_kernel_is_interval_average_flow() -> None:
    tc = np.array([7200.0], dtype=np.float64)
    area = np.array([12345.0], dtype=np.float64)
    transformer = SCSUnitHydrograph(tc=tc, area=area, dt=1800.0)
    kernel = transformer.kernel[:, 0]

    # For 1 m runoff depth, sum(flow * dt) equals basin area (m^3).
    volume_from_kernel = np.sum(kernel * transformer.dt)
    assert np.isclose(volume_from_kernel, area[0], rtol=0.0, atol=1e-9)


def test_requires_area_vector() -> None:
    tc = np.array([100.0, 200.0], dtype=np.float64)
    area = np.array([1.0], dtype=np.float64)
    try:
        SCSUnitHydrograph(tc=tc, area=area, dt=60.0)
    except ValueError as e:
        assert 'same length' in str(e)
    else:
        raise AssertionError('Expected ValueError for tc/area length mismatch')
