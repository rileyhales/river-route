# Test Data

All tests use data under `tests/data/`, with shared fixtures and parametrization defined in
`tests/conftest.py`.

## Directory Layout

```
tests/data/
├── era5/                                # shared gridded runoff inputs (NetCDF)
├── discharge/
│   └── vpu=<id>/discharge_*.nc          # reference routed discharge
├── qlateral/
│   └── vpu=<id>/qlateral_*.nc           # reference lateral inflow at river level
├── routing-configs/
│   └── vpu=<id>/
│       ├── routing_parameters.parquet   # RRv1 params
│       ├── connectivity.parquet         # RRv1 connectivity
│       └── gridweights_ERA5_vpu=<id>.nc
└── hydrography/
    └── vpu=<id>/
        ├── catchments_<id>.parquet
        └── streams_<id>.gpkg
```

## How `conftest.py` Works

`tests/conftest.py` discovers test VPUs from `tests/data/discharge/vpu=*` directories.
For each candidate VPU it builds an `RFSv2ConfigsData` object, then validates required files.
Only valid VPUs are used.

During `RFSv2ConfigsData.__post_init__`, fixture setup also:

1. Builds `rr2_params.parquet` from `routing_parameters.parquet + connectivity.parquet`
2. Builds a unit-hydrograph kernel file:
   `kernel.scs_triangular.kfactor=5.tr=3600.npz`
3. Collects sorted `discharge_files` and `qlateral_files`

Any test function with a `vpu` argument is automatically parametrized by
`pytest_generate_tests`, with IDs like `vpu=718`.

### `RFSv2ConfigsData`

| Field | Type | Description |
|---|---|---|
| `number` | `int` | VPU number |
| `configs_dir` | `Path` | `tests/data/routing-configs/vpu=<id>/` |
| `hydrography_dir` | `Path` | `tests/data/hydrography/vpu=<id>/` |
| `discharge_dir` | `Path` | `tests/data/discharge/vpu=<id>/` |
| `qlateral_dir` | `Path` | `tests/data/qlateral/vpu=<id>/` |
| `discharge_files` | `list[str]` | Populated from `discharge*.nc` |
| `qlateral_files` | `list[str]` | Populated from `qlateral_*.nc` |

Useful computed properties include `rr1_params_file`, `rr1_connectivity_file`,
`rr2_params_file`, `grid_weights_file`, `catchments`, and `streams`.

## Test Modules and Data Use

All tests using `vpu: RFSv2ConfigsData` require the base VPU fixture to validate.
That validation requires configs, hydrography, discharge, and qlateral data to exist.

| Test file | Uses `vpu` fixture | Uses ERA5 files | Uses qlateral files | Uses discharge reference |
|---|:---:|:---:|:---:|:---:|
| `test_muskingum.py` | x |  |  |  |
| `test_rapid_muskingum.py` | x | x | x | x |
| `test_runoff.py` | x | x | x |  |
| `test_tools.py` | x |  |  |  |
| `test_unit_muskingum.py` | x |  |  |  |
| `test_uhkernels.py` |  |  |  |  |

## CI Workflow

`.github/workflows/tests.yaml` runs on manual dispatch and executes:

1. `./tests/download_test_data.sh` (cached in GitHub Actions `tests/data`)
2. `pip install ".[test]"`
3. `pytest tests/ -v -s`

The test job matrix runs Python `3.12`, `3.13`, and `3.14`.

## Running Tests Locally

```bash
./tests/download_test_data.sh
pytest tests/ -v -s
```

`tests/download_test_data.sh` downloads `routing-test-data.zip`, extracts it to `tests/data/`,
then copies additional VPU-specific routing-config and hydrography files from S3.

## Adding a New VPU

To include a new VPU in the parametrized test set:

1. Add `tests/data/discharge/vpu=<id>/` with at least one `discharge*.nc`
2. Add `tests/data/qlateral/vpu=<id>/` with at least one `qlateral_*.nc`
3. Add `tests/data/routing-configs/vpu=<id>/` with:
   `routing_parameters.parquet`, `connectivity.parquet`, and `gridweights_ERA5_vpu=<id>.nc`
4. Add `tests/data/hydrography/vpu=<id>/` with:
   `catchments_<id>.parquet` and `streams_<id>.gpkg`

No code changes are required if file names and directory structure match these conventions.
