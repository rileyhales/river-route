## Finding Inputs and Config Files at Runtime

Instead of manually preparing config files in advance, you may want to generate them in your code which executes the
routing. This is useful when you have a large number of routing runs to perform or if you want to automate the process.
Depending on your preference, you may want to generate many config files in advance or store them for repeatability and
future use.

The following code snippet demonstrates how to identify the essential input arguments and pass them as keyword arguments
to the `TeleportMuskingum` class. You could alternatively write the inputs to a YAML or JSON file and use that config file
instead.

```python
import glob
import os

import river_route as rr

root_dir = '/path/to/root/directory'
vpu_name = 'sample-project'

configs = os.path.join(root_dir, 'configs', vpu_name)
params_file = os.path.join(configs, 'params.parquet')

volume_files = sorted(glob.glob(f'/path/to/runoffs/directory/*.nc'))

outputs = os.path.join(root_dir, 'outputs', vpu_name)
output_files = [os.path.join(outputs, f'Qout_{os.path.basename(f)}') for f in volume_files]
os.makedirs(outputs, exist_ok=True)

m = (
    rr
    .TeleportMuskingum(**{
        'routing_params_file': params_file,
        'catchment_volumes_files': volume_files,
        'discharge_files': output_files,
    })
    .route()
)
```

## Customizing Outputs

You can override the default function used by `river-route` when writing routed flows to disk.
By default, routed discharges are written to netCDF.

A single netCDF is not ideal for all use cases, so you can override it to store your data how you prefer. Some examples
of reasons you would want to do this include appending the outputs to an existing file, writing values to a
database, or to add metadata or attributes to the file.

You can override the `write_discharges` method directly in your code or use the `set_write_discharges` method. Your custom
function must accept exactly 4 arguments:

1. `dates`: datetime array for rows in the discharge array.
2. `discharge_array`: routed discharge array with shape `(time, river_id)`.
3. `discharge_file`: path to the output file.
4. `runoff_file`: path to the runoff input used to produce this output.

As an example, you might want to write output as Parquet instead.

```python title="Write Routed Flows to Parquet"
import pandas as pd
import xarray as xr

import river_route as rr


def custom_write_discharges(dates, discharge_array, discharge_file: str, runoff_file: str) -> None:
    with xr.open_dataset(runoff_file) as runoff_ds:
        river_ids = runoff_ds['river_id'].values
    df = pd.DataFrame(discharge_array, index=pd.to_datetime(dates), columns=river_ids)
    df.to_parquet(discharge_file)
    return


(
    rr
    .TeleportMuskingum('../../examples/config.yaml')
    .set_write_discharges(custom_write_discharges)
    .route()
)
```

```python title="Write Routed Flows to SQLite"
import pandas as pd
import sqlite3
import xarray as xr

import river_route as rr


def write_discharges_to_sqlite(dates, discharge_array, discharge_file: str, runoff_file: str) -> None:
    with xr.open_dataset(runoff_file) as runoff_ds:
        river_ids = runoff_ds['river_id'].values
    df = pd.DataFrame(discharge_array, index=pd.to_datetime(dates), columns=river_ids)
    conn = sqlite3.connect(discharge_file)
    df.to_sql('routed_flows', conn, if_exists='replace')
    conn.close()
    return


(
    rr
    .TeleportMuskingum('config.yaml')
    .set_write_discharges(write_discharges_to_sqlite)
    .route()
)
```

```python title="Append Routed Flows to Existing netCDF"
import os

import xarray as xr

import river_route as rr


def append_to_existing_file(dates, discharge_array, discharge_file: str, runoff_file: str) -> None:
    ensemble_number = os.path.basename(runoff_file).split('_')[1]
    ds = xr.load_dataset(discharge_file)
    ds['Q'].loc[dict(ensemble=ensemble_number)] = discharge_array
    ds.to_netcdf(discharge_file)
    return


(
    rr
    .TeleportMuskingum('config.yaml')
    .set_write_discharges(append_to_existing_file)
    .route()
)
```

```python title="Save a Subset of the Routed Flows"
import pandas as pd
import xarray as xr

import river_route as rr


def save_partial_results(dates, discharge_array, discharge_file: str, runoff_file: str) -> None:
    with xr.open_dataset(runoff_file) as runoff_ds:
        river_ids = runoff_ds['river_id'].values
    df = pd.DataFrame(discharge_array, index=pd.to_datetime(dates), columns=river_ids)
    river_ids_to_save = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    df = df[river_ids_to_save]
    df.to_parquet(discharge_file)
    return


(
    rr
    .TeleportMuskingum('config.yaml')
    .set_write_discharges(save_partial_results)
    .route()
)
```

## Clark-Muskingum Advanced Features

### Setting Time-Area Curves Programmatically

Instead of providing a `time_area_file`, you can set per-catchment time-area curves directly as a NumPy array.
The array must have shape `(n_points, n_catchments)` with cumulative area fractions at evenly spaced normalized
times over `[0, 1]`.

```python
import numpy as np
import river_route as rr

# Build or load your curves
# shape: (n_curve_points, n_catchments)
curves = np.load('my_time_area_curves.npy')

(
    rr
    .ClarkMuskingum('config_clark.yaml')
    .set_time_area_curves(curves)
    .route()
)
```

If neither `time_area_file` nor `set_time_area_curves()` is used, `ClarkMuskingum` falls back to the SCS
dimensionless time-area curve applied uniformly to all catchments.

### Pre-computing Lateral Inflows

For repeated routing runs over the same input volumes (e.g., sensitivity analysis, calibration), you can
pre-compute the Clark UH-transformed lateral inflows once and reuse them. This avoids re-applying the
convolution on every run.

```python
import river_route as rr

router = rr.ClarkMuskingum('config_clark.yaml')

# Step 1: pre-compute and write UH-transformed inflows to disk
precomputed_files = ['lateral_1.nc', 'lateral_2.nc']
router.compute_lateral_inflows(precomputed_files)

# Step 2: route using the precomputed inflows directly
(
    rr
    .ClarkMuskingum(
        routing_params_file='params.parquet',
        catchment_volumes_files=precomputed_files,
        discharge_files=['discharge_1.nc', 'discharge_2.nc'],
        precomputed_lateral_inflows=True,
    )
    .route()
)
```

The `compute_lateral_inflows()` method writes netCDF files in the same format as catchment volume files.
Setting `precomputed_lateral_inflows=True` tells the router to use the volumes directly as lateral inflow
rates rather than applying the UH transform again.

### Precomputing Unit Hydrographs

Unit hydrographs are computed automatically the first time `route()` is called. If you want to inspect them
or trigger precomputation early (e.g., before a batch of runs), call `precompute_unit_hydrographs()` directly.
Note that the router must have already read the routing parameters and time metadata.

```python
import river_route as rr

router = (
    rr
    .ClarkMuskingum('config_clark.yaml')
)
# precompute_unit_hydrographs is called automatically during route(),
# but can also be called explicitly after initialization
router.route()
print('UH matrix shape:', router._uh_matrix.shape)
```
