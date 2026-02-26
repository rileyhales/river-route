## Customizing Outputs

By default, `river-route` writes routed discharges to netCDF. You can override this behavior by
passing a custom function to `set_write_discharges()`. This is useful when you want to write to a
different format, append to an existing file, write to a database, or add custom metadata.

Your custom function must accept exactly four arguments:

1. `dates`: datetime array for rows in the discharge array.
2. `discharge_array`: routed discharge array with shape `(time, river_id)`.
3. `discharge_file`: path to the output file provided in your config.
4. `runoff_file`: path to the runoff input used to produce this output.

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
    .RapidMuskingum('config.yaml')
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
    .RapidMuskingum('config.yaml')
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
    .RapidMuskingum('config.yaml')
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
    .RapidMuskingum('config.yaml')
    .set_write_discharges(save_partial_results)
    .route()
)
```
