## Configuration File

`river-route` computations are controlled by config values passed as keyword arguments or from a YAML/JSON file.
The three routers (`Muskingum`, `RapidMuskingum`, `UnitMuskingum`) share a set of base config keys and each
adds router-specific keys.

1. Paths to input and output files
2. Routing and timestep options
3. Input/output variable names
4. Logging options

## Minimum Required Inputs

All routing classes require the following 2 configuration options:

- `params_file` - path to the [routing parameters file](io-file-schema.md#routing-parameters) (parquet)
- `discharge_dir` - directory where [routed discharge](io-file-schema.md#routed-discharge) output files are written (netCDF, named after input files)

`Muskingum` (channel routing only, no lateral inflows) requires:

- `channel_state_init_file` - parquet state file to initialize discharge
- `dt_routing` - routing timestep in seconds
- `dt_total` - total simulation duration in seconds

`RapidMuskingum` requires:

- one water input source:
    - `catchment_runoff_files`, or
    - `runoff_grid_files` plus `grid_weights_file`

`UnitMuskingum` requires:

- `transformer_kernel_file` - pre-computed parquet kernel
- one water input source:
    - `catchment_runoff_files`, or
    - `runoff_grid_files` plus `grid_weights_file`

## Example Configuration YAML

The general template below covers all routers with annotations. Minimal templates
are also available in the `examples/` directory:
`config_muskingum.yaml`, `config_rapid_muskingum.yaml`, `config_unit_muskingum.yaml`.

```yaml title="config.yaml"
{% include-markdown "../../examples/config.yaml" %}
```

## Required Config Keys

| Config key                     | Description                                            |  Muskingum   |              RapidMuskingum               |               UnitMuskingum               |
|--------------------------------|--------------------------------------------------------|:------------:|:-----------------------------------------:|:-----------------------------------------:|
| **core**                       |                                                        |              |                                           |                                           |
| `params_file`                  | Routing parameters parquet.                            | **Required** |               **Required**                |               **Required**                |
| **state**                      |                                                        |              |                                           |                                           |
| `channel_state_init_file`      | Parquet with initial channel state (column `Q`)        | **Required** |          optional - default to 0          |          optional - default to 0          |
| `channel_state_final_file`     | Path to save final channel state                       |   optional   |                 optional                  |                 optional                  |
| **output**                     |                                                        |              |                                           |                                           |
| `discharge_dir`                | Directory for output discharge files                   | one of `discharge_dir` / `discharge_files` required | one of `discharge_dir` / `discharge_files` required | one of `discharge_dir` / `discharge_files` required |
| `discharge_files`              | Explicit output paths (alternative to `discharge_dir`) | one of `discharge_dir` / `discharge_files` required | one of `discharge_dir` / `discharge_files` required | one of `discharge_dir` / `discharge_files` required |
| **input data**                 |                                                        |              |                                           |                                           |
| `catchment_runoff_files`       | Per-catchment runoff time series (m³ or m)             |              |                _Option 1_                 |                _Option 1_                 |
| `runoff_grid_files`            | Gridded runoff depths; need `grid_weights_file`        |              |                _Option 2_                 |                _Option 2_                 |
| `grid_weights_file`            | Table to convert depth grids to lateral inflow         |              |                _Option 2_                 |                _Option 2_                 |
| **unit hydrograph**            |                                                        |              |                                           |                                           |
| `transformer_kernel_file`      | Pre-computed convolution kernel parquet                |              |                                           |               **Required**                |
| `transformer_state_init_file`  | Parquet with initial transformer state                 |              |                                           |                 optional                  |
| `transformer_state_final_file` | Path to save final transformer state                   |              |                                           |                 optional                  |
| **time**                       |                                                        |              |                                           |                                           |
| `dt_total`                     | Total simulation duration in seconds                   | **Required** | optional - [time docs](time-variables.md) | optional - [time docs](time-variables.md) |
| `dt_discharge`                 | Output timestep (s)                                    |   optional   | optional - [time docs](time-variables.md) | optional - [time docs](time-variables.md) |
| `dt_runoff`                    | Runoff data timestep (s)                               |              | optional - [time docs](time-variables.md) | optional - [time docs](time-variables.md) |
| `dt_routing`                   | Routing computational timestep (s)                     | **Required** | optional - [time docs](time-variables.md) | optional - [time docs](time-variables.md) |
| `start_datetime`               | Simulation start date for output timestamps            |   optional   |                                           |                                           |

## Optional Configs with Defaults

| Config Key                      | Description                                                                | Default                                       |
|---------------------------------|----------------------------------------------------------------------------|-----------------------------------------------|
| `log`                           | Enable or disable logging                                                  | `True`                                        |
| `progress_bar`                  | Show tqdm progress bar                                                     | `True`                                        |
| `log_level`                     | Logger level (`'INFO'`, `'DEBUG'`, etc.)                                   | `'INFO'`                                      |
| `log_stream`                    | `'stdout'` or a file path                                                  | `'stdout'`                                    |
| `log_format`                    | Python logging format string                                               | `'%(levelname)s - %(asctime)s - %(message)s'` |
| `var_river_id`                  | River ID dimension name in files                                           | `'river_id'`                                  |
| `var_discharge`                 | Discharge variable name in output                                          | `'Q'`                                         |
| `var_catchment_runoff_variable` | Runoff variable name in `catchment_runoff_files`                           | `'runoff'`                                    |
| `var_runoff_depth`              | Depth variable name in `runoff_grid_files`                                 | `'ro'`                                        |
| `var_x`                         | X-dimension name in depth grids                                            | `'x'`                                         |
| `var_y`                         | Y-dimension name in depth grids                                            | `'y'`                                         |
| `var_t`                         | Time dimension name in depth grids                                         | `'time'`                                      |
| `runoff_accumulation_type`      | Whether lateral inflows are `'incremental'` or `'cumulative'`              | `'incremental'`                               |
| `runoff_processing_mode`        | `'sequential'` (carry state forward) or `'ensemble'` (average final state) | `'sequential'`                                |
