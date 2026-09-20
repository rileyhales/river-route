"""
Tests for Configs.discharge_dtype, which narrows the discharge buffer to float16.

The routing math stays float32 and the channel state is never narrowed, so these check both halves: the kernels
really store float16 bit patterns, and each writer puts them in the file as the format can hold them.
"""

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
import xarray as xr
from conftest import SyntheticNetwork

import river_route as rr
from river_route.router import writers

# float16 keeps 11 significand bits, so round to nearest bounds the relative error at 2^-11.
TOLERANCE = 2**-11


def route(network: SyntheticNetwork, discharge_dtype: str, out_name: str, writer=None, **kwargs) -> str:
    out = network.path(out_name)
    router = rr.Router(
        rr.Configs(
            forcing='vlateral',
            params_file=str(network.params_file),
            vlateral_files=[str(network.vlateral_file)],
            channel_state_init_file=str(network.state_file),
            dt_routing=network.dt_runoff,
            discharge_files=[out],
            discharge_dtype=discharge_dtype,
            log=False,
            progress_bar=False,
            **kwargs,
        )
    )
    if writer is not None:
        router.set_discharge_writer(writer)
    router.route()
    return out


def route_netcdf(network: SyntheticNetwork, discharge_dtype: str, out_name: str) -> str:
    return route(network, discharge_dtype, out_name, writer=writers.netcdf_writer)


def read_netcdf(path: str) -> np.ndarray:
    with xr.open_dataset(path) as ds:
        return ds['Q'].transpose('time', 'river_id').values


def assert_matches(narrowed: np.ndarray, expected: np.ndarray) -> None:
    """Values agree to float16's rounding bound, measured against the scale of the series: float16 holds no
    precision below 6.1e-5, where its subnormals begin, so the bound is applied to the peak, not per value."""
    np.testing.assert_allclose(narrowed, expected, rtol=TOLERANCE, atol=expected.max() * TOLERANCE)


def test_float32_is_the_default():
    assert rr.Configs(params_file='x.parquet').discharge_dtype == 'float32'


def test_kernels_fill_a_narrowed_buffer(network: SyntheticNetwork):
    """The array the kernels write, and hand to the writer, is half width float16 bit patterns."""
    seen = {}

    def capture(router, dates, discharge_array, discharge_file, runoff_file=''):
        seen['dtype'] = discharge_array.dtype
        seen['itemsize'] = discharge_array.itemsize

    route(network, 'float16', 'q.nc', writer=capture)
    assert seen['dtype'] == np.uint16
    assert seen['itemsize'] == 2


def test_netcdf_widens_to_float32(network: SyntheticNetwork):
    """netCDF has no half type, so the file is float32 and holds the rounded values."""
    expected = read_netcdf(route_netcdf(network, 'float32', 'q32.nc'))
    narrowed_file = route_netcdf(network, 'float16', 'q16.nc')
    with xr.open_dataset(narrowed_file) as ds:
        assert ds['Q'].dtype == np.float32
        assert ds['Q'].dims == ('river_id', 'time')
    assert expected.any()
    assert_matches(read_netcdf(narrowed_file), expected)


def test_zarr_and_parquet_keep_float16(network: SyntheticNetwork):
    """Both formats hold half floats natively, so a narrowed run halves the file instead of widening."""
    expected = read_netcdf(route_netcdf(network, 'float32', 'q32.nc'))
    zarr_out = route(network, 'float16', 'q.zarr', writer=writers.zarr_writer)
    parquet_out = route(network, 'float16', 'q.parquet', writer=writers.parquet_writer)
    with xr.open_zarr(zarr_out) as ds:
        assert ds['Q'].dtype == np.float16
        assert_matches(ds['Q'].transpose('time', 'river_id').values.astype(np.float32), expected)
    frame = pd.read_parquet(parquet_out).set_index('river_id')
    assert pq.read_schema(parquet_out).field(1).type == 'halffloat'
    assert_matches(frame.to_numpy().T.astype(np.float32), expected)


def test_float32_runs_still_write_float32(network: SyntheticNetwork):
    zarr_out = route(network, 'float32', 'q.zarr', writer=writers.zarr_writer)
    with xr.open_zarr(zarr_out) as ds:
        assert ds['Q'].dtype == np.float32
        assert ds['Q'].dims == ('river_id', 'time')


def test_narrowed_dtype_cannot_be_resampled(network: SyntheticNetwork):
    with pytest.raises(ValueError, match='cannot be averaged'):
        route(network, 'float16', 'q.nc', dt_runoff=network.dt_runoff, dt_discharge=network.dt_runoff * 2)
