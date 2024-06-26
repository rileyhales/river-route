## Writing Routed Flows to Disc

You can override the default function used by `river-route` when writing routed flows to disc. The MuskingumCunge class
formats the discharge data into a Pandas DataFrame and then calls the `write_outflows` method. By default, this function
writes the dataframe to a netCDF file.

A single netCDF is not ideal for all use cases so you can override it to store your data how you prefer. Some examples
of reasons you would want to do this include appending the outputs to an existing file, writing the DataFrame to a
database, or to add metadata or attributes to the file.

You can override the `write_outflows` method directly in your code or use the `set_write_outflows` method. Your custom
function should accept exactly 3 keyword arguments:

1. `df`: a Pandas DataFrame with a datetime index, river id numbers as column labels, and float discharge values.
2. `output_file`: a string with the path to the output file.
3. `runoff_file`: a string with the path to the runoff file used to produce this output.

As an example, you might want to write the output DataFrame to a Parquet file instead.

```python title="Write Routed Flows to Parquet"
import pandas as pd

import river_route as rr


def custom_write_outflows(df: pd.DataFrame, output_file: str, runoff_file: str) -> None:
    df.to_parquet(output_file)
    pass


(
    rr
    .MuskingumCunge('config.yaml')
    .set_write_outflows(custom_write_outflows)
    .route()
)
```

```python title="Write Routed Flows to SQLite"
import pandas as pd
import sqlite3

import river_route as rr


def write_outflows_to_sqlite(df: pd.DataFrame, output_file: str, runoff_file: str) -> None:
    conn = sqlite3.connect(output_file)
    df.to_sql('routed_flows', conn, if_exists='replace')
    conn.close()
    return


(
    rr
    .MuskingumCunge('config.yaml')
    .set_write_outflows(write_outflows_to_sqlite)
    .route()
)
```

```python title="Append Routed Flows to Existing netCDF"
import os

import pandas as pd
import xarray as xr

import river_route as rr


def append_to_existing_file(df: pd.DataFrame, output_file: str, runoff_file: str) -> None:
    ensemble_number = os.path.basename(runoff_file).split('_')[1]
    with xr.open_dataset(output_file) as ds:
        ds.sel(ensemble=ensemble_number).Q = df.values
        ds.to_netcdf(output_file)
    return


(
    rr
    .MuskingumCunge('config.yaml')
    .set_write_outflows(append_to_existing_file)
    .route()
)
```

```python title="Save a Subset of the Routed Flows"
import pandas as pd

import river_route as rr


def save_partial_results(df: pd.DataFrame, output_file: str, runoff_file: str) -> None:
    river_ids_to_save = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    df = df[river_ids_to_save]
    df.to_parquet(output_file)
    return


(
    rr
    .MuskingumCunge('config.yaml')
    .set_write_outflows(save_partial_results)
    .route()
)
```
