# Test Data

All test data lives in `tests/data/` at the project root. Tests resolve paths through
`tests/conftest.py`, which auto-discovers VPUs from Hive-style `vpu=` directory partitioning.

## Directory Layout

```
tests/data/
├── era5/                              # global ERA5 forcing — shared across VPUs
│   ├── era5_194001.nc
│   ├── era5_194002.nc
│   └── era5_194003.nc
├── routing-configs/
│   ├── params.parquet                 # routing params (river_id, downstream_river_id, k, x)
│   └── gridweights_era5.nc            # grid weight table
├── discharge/
│   ├── discharge_era5_194001.nc       # known-good routed discharge
│   ├── ...
│   └── discharge_era5_194012.nc
├── catchment-volumes/
│   ├── volumes_194001.nc              # pre-computed catchment volumes (from grid_to_catchment)
│   ├── volumes_194002.nc
│   └── volumes_194003.nc
├── hydrography/
│   ├── catchments_718.parquet         # catchment boundaries (linkno, geometry)
│   └── streams_718.gpkg
```

## VPU Auto-Discovery

At import time, `conftest.py` scans `routing-configs/`, `discharge/`, `catchment-volumes/`,
and `hydrography/` for `vpu=*` subdirectories, unions the VPU IDs, and builds a `RFSv2ConfigsData`
dataclass for each. File discovery uses globs: `*params*.parquet` for params, `gridweights*.nc`
for weights.

Any test function with a `vpu` parameter is automatically parametrized over all discovered VPUs
via the `pytest_generate_tests` hook. Test IDs become e.g. `test_foo[vpu=718]`.

### `RFSv2ConfigsData` Dataclass

| Field                   | Type   | Description                                      |
|-------------------------|--------|--------------------------------------------------|
| `vpu_id`                | `str`  | The VPU identifier (e.g. `"718"`)                |
| `params_file`           | `Path` | Routing parameters parquet                       |
| `grid_weights_file`     | `Path` | Grid weight table NetCDF                         |
| `discharge_dir`         | `Path` | Directory of known-good discharge files          |
| `catchment_volumes_dir` | `Path` | Directory of pre-computed catchment volume files |
| `hydrography_dir`       | `Path` | Directory of hydrography files                   |

Fields are `None` when the corresponding data is not available for that VPU.

### `conftest.py` Exports

| Export                              | Type                     | Description                                                          |
|-------------------------------------|--------------------------|----------------------------------------------------------------------|
| `RFSv2ConfigsData`                  | `dataclass`              | Per-VPU test data paths                                              |
| `ALL_VPUS`                          | `list[RFSv2ConfigsData]` | All discovered VPUs                                                  |
| `ERA5_DIR`                          | `Path`                   | Directory containing ERA5 forcing files                              |
| `ERA5_KWARGS`                       | `dict`                   | `{'var_y': 'latitude', 'var_x': 'longitude', 'var_t': 'valid_time'}` |
| `era5_files()`                      | `list[str]`              | Sorted list of ERA5 file paths matching `era5_1940*.nc`              |
| `known_discharge_files(vpu)`        | `list[str]`              | Sorted list of known discharge file paths for a VPU                  |
| `known_catchment_volume_files(vpu)` | `list[str]`              | Sorted list of catchment volume file paths for a VPU                 |
| `skip_if_vpu_missing(vpu, ...)`     | `None`                   | Skips the test if any named VPU field is unavailable                 |

## Data Availability Tiers

Tests are designed to run with whatever data is available. Missing data causes individual tests
to skip, not fail.

| Tier   | Files                                                                         | Size    | Source        |
|--------|-------------------------------------------------------------------------------|---------|---------------|
| Core   | `routing-configs/vpu=N/params.parquet`                                        | ~312 KB | Cloud storage |
| Medium | `routing-configs/vpu=N/gridweights*.nc`                                       | ~625 KB | Cloud storage |
| Large  | `era5/`, `discharge/vpu=N/`, `catchment-volumes/vpu=N/`, `hydrography/vpu=N/` | ~2.5 GB | Cloud storage |

All tiers are downloaded in CI via `tests/download_test_data.sh` from a cloud storage URL
configured as the `DATA_URL` repository variable.

## Test Dependency Matrix

| Test file                 | `params_file` | `grid_weights_file` | `era5/` | `discharge/` | `catchment-volumes/` | `hydrography/` |
|---------------------------|:-------------:|:-------------------:|:-------:|:------------:|:--------------------:|:--------------:|
| `test_config.py`          |       x       |          x          |    x    |              |                      |                |
| `test_metrics.py`         |               |                     |         |              |                      |                |
| `test_muskingum.py`       |       x       |                     |         |              |                      |                |
| `test_rapid_muskingum.py` |       x       |          x          |    x    |      x       |          x           |                |
| `test_runoff.py`          |       x       |          x          |    x    |              |                      |       x        |
| `test_tools.py`           |       x       |          x          |         |              |                      |                |
| `test_unit_muskingum.py`  |       x       |                     |         |              |                      |                |

## CI Integration

The GitHub Actions workflow (`.github/workflows/tests.yaml`) runs the test suite on Python
3.12, 3.13, and 3.14:

1. `tests/download_test_data.sh` downloads the test data tarball from the `DATA_URL` repository variable
2. The tarball extracts preserving `vpu=` directory structure into `tests/data/`
3. `pytest tests/ -v` runs with auto-discovered VPUs

The download script is idempotent — it skips if `routing-configs/vpu=*` directories already exist.

## Adding a New VPU

To add test coverage for a new VPU (e.g. VPU 203):

1. Create the data directories under `tests/data/`:
    ```
    routing-configs/vpu=203/
    discharge/vpu=203/
    catchment-volumes/vpu=203/    # optional
    hydrography/vpu=203/          # optional
    ```

2. Populate with data files:
    - `routing-configs/vpu=203/*params*.parquet` — routing parameters (must contain
      `river_id`, `downstream_river_id`, `k`, `x`)
    - `routing-configs/vpu=203/gridweights*.nc` — grid weight table
    - `discharge/vpu=203/discharge_era5_1940*.nc` — known-good discharge
    - `catchment-volumes/vpu=203/volumes_1940*.nc` — pre-computed catchment volumes

3. No code changes needed. `discover_vpus()` finds the new VPU automatically and all
   VPU-parametrized tests run for it. Tests skip gracefully for any missing data files.

4. Re-package test data for CI:
    ```bash
    ./tests/package_test_data.sh 718 203
    ```

## Shell Scripts

### `tests/package_test_data.sh`

Packages test data into a tarball for uploading as a GitHub release asset. Accepts VPU IDs
as arguments (defaults to `718`).

```bash
./tests/package_test_data.sh              # package VPU 718 only
./tests/package_test_data.sh 718 203      # package multiple VPUs
```

### `tests/download_test_data.sh`

Downloads test data from cloud storage into `tests/data/`. Used by CI and developers who
don't have local data. Uses the default S3 URL, or override with `DATA_URL`.

```bash
./tests/download_test_data.sh
```
