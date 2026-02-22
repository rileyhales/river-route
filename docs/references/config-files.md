## Configuration File

`river-route` computations are controlled by config values passed as keyword arguments or from a YAML/JSON file.
The three routers (`Muskingum`, `RapidMuskingum`, `UnitMuskingum`) share a set of base config keys and each
adds router-specific keys. `RapidMuskingum` and `UnitMuskingum` are parallel routers — they do not share
input data config keys, so a config file written for one cannot be used with the other without changes.

1. Paths to input and output files
2. Routing and timestep options
3. Input/output variable names
4. Logging options

## Minimum Required Inputs

`Muskingum` (channel routing only, no lateral inflows) requires:

- `routing_params_file` - parquet with columns: `river_id`, `downstream_river_id`, `k`, `x`
- `channel_state_file` - parquet state file to initialize discharge
- `discharge_files` - output netCDF path as a single-element list

`RapidMuskingum` requires:

- `routing_params_file` - path to the [routing parameters file](io-file-schema.md#routing-parameters) (parquet)
- one water input source:
    - `lateral_volume_files`, or
    - `runoff_depth_grids` plus `grid_weights_file`
- `discharge_files` - path(s) where [routed discharge](io-file-schema.md#routed-discharge) output files will be saved (netCDF)

`UnitMuskingum` requires:

- `routing_params_file` - parquet with columns: `river_id`, `downstream_river_id`, `k`, `x`
- one transformer source:
    - `transformer_kernel_file` - pre-computed parquet kernel, or
    - a transformer injected via `set_transformer()` before calling `route()`
- one water input source:
    - `lateral_depth_files`, or
    - `runoff_depth_grids` plus `grid_weights_file`
- `discharge_files`

## Example Configuration YAML

```yaml title="RapidMuskingum Config"
{ % include-markdown "../../examples/config.yaml" % }
```

## Required Config Keys

| Config key                     | Description                                     |  Muskingum   |             RapidMuskingum             |               UnitMuskingum               |
|--------------------------------|-------------------------------------------------|:------------:|:-----------------------------------------:|:-----------------------------------------:|
| **core**                       |                                                 |              |                                           |                                           |
| `routing_params_file`          | Routing parameters parquet.                     | **Required** |               **Required**                |               **Required**                |
| **state**                      |                                                 |              |                                           |                                           |
| `channel_state_file`           | Parquet with initial channel state (column `Q`) | **Required** |          optional - default to 0          |          optional - default to 0          |
| `final_channel_state_file`     | Path to save final channel state                |   optional   |                 optional                  |                 optional                  |
| **output**                     |                                                 |              |                                           |                                           |
| `discharge_files`              | Output discharge file path(s)                   | **Required** |               **Required**                |               **Required**                |
| **input data**                 |                                                 |              |                                           |                                           |
| `lateral_volume_files`         | Per-catchment runoff volume time series (m³)    |              |                _Option 1_                 |                                           |
| `lateral_depth_files`          | Per-catchment runoff depth time series (m)      |              |                                           |                _Option 1_                 |
| `runoff_depth_grids`           | Gridded runoff depths; need `grid_weights_file` |              |                _Option 2_                 |                _Option 2_                 |
| `grid_weights_file`            | Table to convert depth grids to lateral inflow  |              |                _Option 2_                 |                _Option 2_                 |
| **unit hydrograph**            |                                                 |              |                                           |                                           |
| `transformer_kernel_file`      | Pre-computed convolution kernel parquet         |              |                                           |               **Required**                |
| `transformer_state_file`       | Parquet with initial transformer state          |              |                                           |                 optional                  |
| `final_transformer_state_file` | Path to save final transformer state            |              |                                           |                 optional                  |
| **time**                       |                                                 |              |                                           |                                           |
| `dt_total`                     | Total simulation duration in seconds            | **Required** | optional - [time docs](time-variables.md) | optional - [time docs](time-variables.md) |
| `dt_discharge`                 | Output timestep (s)                             |   optional   | optional - [time docs](time-variables.md) | optional - [time docs](time-variables.md) |
| `dt_runoff`                    | Runoff data timestep (s)                        |              | optional - [time docs](time-variables.md) | optional - [time docs](time-variables.md) |
| `dt_routing`                   | Routing computational timestep (s)              | **Required** | optional - [time docs](time-variables.md) | optional - [time docs](time-variables.md) |
| `start_datetime`               | Simulation start date for output timestamps     |   optional   |                                           |                                           |

## Optional Configs with Defaults

| Config Key             | Description                                               | Default                                       |
|------------------------|-----------------------------------------------------------|-----------------------------------------------|
| `log`                  | Enable or disable logging                                 | True                                          |
| `progress_bar`         | Show tqdm progress bar                                    | True                                          |
| `log_level`            | Logger level (`'INFO'`, `'DEBUG'`, etc.)                  | `'INFO'`                                      |
| `log_stream`           | `'stdout'` or a file path                                 | `'stdout'`                                    |
| `log_format`           | Python logging format string                              | `'%(asctime)s - %(levelname)s - %(message)s'` |
| `var_river_id`         | River ID dimension name in files.                         | `'river_id'`                                  |
| `var_discharge`        | Discharge variable name in output.                        | `'Q'`                                         |
| `var_catchment_volume` | Volume variable name in `lateral_volume_files`.           | `'volume'`                                    |
| `var_runoff_depth`     | Depth variable name in depth grid files.                  | `'ro'`                                        |
| `var_x`                | X-dimension name in depth grids.                          | `'x'`                                         |
| `var_y`                | Y-dimension name in depth grids.                          | `'y'`                                         |
| `var_t`                | Time dimension name in depth grids.                       | `'time'`                                      |
| `runoff_type`          | Whether the lateral inflows are incremental or cumulative | `'incremental'`                               |                                                            |                             d                              |
| `input_type`           | [advanced tutorials](../tutorial/routing-ensembles.md)    | `'sequential'` or `'ensemble'`                |
