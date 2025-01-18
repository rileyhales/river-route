## Finding Inputs and Config Files at Runtime

Instead of manually preparing config files in advance, you may want to generate them in your code which executes the 
routing. This is useful when you have a large number of routing runs to perform or if you want to automate the process.
Depending on your preference, you may want to generate many config files in advance or store them for repeatability and 
future use.

The following code snippet demonstrates how to identify the essential input arguments and pass them as keyword arguments 
to the `Muskingum` class. You could alternatively write the inputs to a YAML or JSON file and use that config file 
instead.

```python
import glob
import os

import river_route as rr

root_dir = '/path/to/root/directory'
vpu_name = 'sample-project'

configs = os.path.join(root_dir, 'configs', vpu_name)
params_file = os.path.join(configs, 'params.parquet')
connectivity_file = os.path.join(configs, 'connectivity.parquet')

volume_files = sorted(glob.glob(f'/path/to/runoffs/directory/*.nc'))

outputs = os.path.join(root_dir, 'outputs', vpu_name)
output_files = [os.path.join(outputs, f'Qout_{os.path.basename(f)}') for f in volume_files]
os.makedirs(outputs, exist_ok=True)

m = (
    rr
    .Muskingum(**{
        'routing_params_file': params_file,
        'connectivity_file': connectivity_file,
        'catchment_volumes_file': volume_files,
        'outflow_file': output_files,
    })
    .route()
)
```

## Customizing Output File Type and Structure

You can override the default function used by `river-route` when writing routed flows to disc. The Muskingum class
formats the discharge data into a Pandas DataFrame and then calls the `write_outflows` method. By default, this function
writes the dataframe to a netCDF file.

A single netCDF is not ideal for all use cases, so you can override it to store your data how you prefer. Some examples
of reasons you would want to do this include appending the outputs to an existing file, writing the DataFrame to a
database, or to add metadata or attributes to the file.

You can override the `write_outflows` method directly in your code or use the `set_write_outflows` method. Your custom
function should accept exactly 3 keyword arguments:

1. `df`: a Pandas DataFrame with a datetime index, river id numbers as column labels, and float discharge values.
2. `outflow_file`: a string with the path to the output file.
3. `runoff_file`: a string with the path to the runoff file used to produce this output.

As an example, you might want to write the output DataFrame to a Parquet file instead.

```python title="Write Routed Flows to Parquet"
import pandas as pd

import river_route as rr


def custom_write_outflows(df: pd.DataFrame, outflow_file: str, runoff_file: str) -> None:
    df.to_parquet(outflow_file)
    return


(
    rr
    .Muskingum('../../examples/config.yaml')
    .set_write_outflows(custom_write_outflows)
    .route()
)
```

```python title="Write Routed Flows to SQLite"
import pandas as pd
import sqlite3

import river_route as rr


def write_outflows_to_sqlite(df: pd.DataFrame, outflow_file: str, runoff_file: str) -> None:
    conn = sqlite3.connect(outflow_file)
    df.to_sql('routed_flows', conn, if_exists='replace')
    conn.close()
    return


(
    rr
    .Muskingum('config.yaml')
    .set_write_outflows(write_outflows_to_sqlite)
    .route()
)
```

```python title="Append Routed Flows to Existing netCDF"
import os

import pandas as pd
import xarray as xr

import river_route as rr


def append_to_existing_file(df: pd.DataFrame, outflow_file: str, runoff_file: str) -> None:
    ensemble_number = os.path.basename(runoff_file).split('_')[1]
    with xr.open_dataset(outflow_file) as ds:
        ds.sel(ensemble=ensemble_number).Q = df.values
        ds.to_netcdf(outflow_file)
    return


(
    rr
    .Muskingum('config.yaml')
    .set_write_outflows(append_to_existing_file)
    .route()
)
```

```python title="Save a Subset of the Routed Flows"
import pandas as pd

import river_route as rr


def save_partial_results(df: pd.DataFrame, outflow_file: str, runoff_file: str) -> None:
    river_ids_to_save = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    df = df[river_ids_to_save]
    df.to_parquet(outflow_file)
    return


(
    rr
    .Muskingum('config.yaml')
    .set_write_outflows(save_partial_results)
    .route()
)
```
