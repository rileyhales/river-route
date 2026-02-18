## Configuration File

`river-route` computations are controlled by config values passed as keyword arguments or from a YAML/JSON file.
Both routers (`Muskingum`, `ClarkMuskingum`) share the same base config keys, and Clark adds a few optional keys.

1. Paths to input and output files
2. Routing and timestep options
3. Input/output variable names
4. Logging options

## Minimum Required Inputs

Every run needs:

- `routing_params_file` - path to the [routing parameters file](io-file-schema.md#routing-parameters) (parquet)
- one water input source:
    - `catchment_volumes_files`, or
    - `runoff_depths_files` plus `weight_table_file`
- `discharge_files` - path(s) where [routed discharge](io-file-schema.md#routed-discharge) output files will be saved (netCDF)

## Example YAML File

Muskingum example:

```yaml title="Config File Example"
{ % include-markdown "../../examples/config.yaml" % }
```

Clark-Muskingum example:

```yaml title="Clark Config File Example"
{ % include-markdown "../../examples/config_clark.yaml" % }
```

## Config Options Table

The following table is the full set of recognized keys.

| Parameter Name                | Required        | Type                         | Description                                                               |
|-------------------------------|-----------------|------------------------------|---------------------------------------------------------------------------|
| `routing_params_file`         | Yes             | file path                    | Routing parameters parquet file.                                          |
| `discharge_files`             | Yes             | file path or list[file path] | Output netCDF discharge file(s). Count must match input file count.       |
| `catchment_volumes_files`     | Conditionally   | file path or list[file path] | Input netCDF volume file(s). Use this or `runoff_depths_files`.           |
| `runoff_depths_files`         | Conditionally   | file path or list[file path] | Input netCDF runoff depth file(s). Use this or `catchment_volumes_files`. |
| `weight_table_file`           | Conditionally   | file path                    | Required if `runoff_depths_files` is provided.                            |
| `input_type`                  | No              | string                       | `sequential` or `ensemble`. Default: `sequential`.                        |
| `runoff_type`                 | No              | string                       | `incremental` or `cumulative`. Default: `incremental`.                    |
| `dt_total`                    | No              | integer                      | Total simulation duration in seconds. Defaults to input duration.         |
| `dt_routing`                  | No              | integer                      | Routing sub-step in seconds. Defaults to `dt_runoff`.                     |
| `dt_discharge`                | No              | integer                      | Output discharge timestep in seconds. Defaults to `dt_runoff`.            |
| `initial_state_file`          | No              | file path                    | Parquet state file used to initialize routing state.                      |
| `final_state_file`            | No              | file path                    | Path where final state parquet is written.                                |
| `time_area_file`              | No (Clark only) | file path                    | Parquet time-area histograms for `ClarkMuskingum`.                        |
| `precomputed_lateral_inflows` | No (Clark only) | boolean                      | If `true`, treat inputs as UH-transformed lateral inflows.                |
| `log`                         | No              | boolean                      | Enable or disable logging. Default: `true`.                               |
| `progress_bar`                | No              | boolean                      | Show tqdm progress bar. Default: same as `log`.                           |
| `log_level`                   | No              | string                       | Logger level. Default: `INFO`.                                            |
| `log_stream`                  | No              | string                       | `stdout` or a file path. Default: `stdout`.                               |
| `log_format`                  | No              | string                       | Python logging format string.                                             |
| `var_x`                       | No              | string                       | X-dimension name in runoff depth grids. Default: `x`.                     |
| `var_y`                       | No              | string                       | Y-dimension name in runoff depth grids. Default: `y`.                     |
| `var_t`                       | No              | string                       | Time dimension name in runoff depth grids. Default: `time`.               |
| `var_runoff_depth`            | No              | string                       | Runoff depth variable name. Default: `ro`.                                |
| `var_catchment_volume`        | No              | string                       | Catchment volume variable name. Default: `volume`.                        |
| `var_river_id`                | No              | string                       | River ID dimension/variable name. Default: `river_id`.                    |
| `var_discharge`               | No              | string                       | Discharge variable name in outputs. Default: `Q`.                         |

| Config key                  | Muskingum | LumpedMuskingum | ClarkMuskingum |
|-----------------------------|-----------|-----------------|----------------|
| routing_params_file         | Required  | Required        | Required       |
| discharge_files             | Required  | Required        | Required       |
| catchment_volumes_files     |           | Required        |                |
| runoff_depths_files         |           | optional        | Required       |
| weight_table_file           |           | optional        | optional       |
| time_area_file              |           |                 | optional       |
| precomputed_lateral_inflows |           |                 | optional       |
| dt_routing                  | Required  | optional        | optional       |
| dt_discharge                | optional  | optional        | optional       |
| dt_runoff                   |           | optional        | optional       |
| initial_state_file          | Required  | optional        | optional       |
| final_state_file            | optional  | optional        | optional       |
| input_type                  |           | optional        | optional       |
| runoff_type                 |           | optional        | optional       |
| log                         | optional  | optional        | optional       |
| progress_bar                | optional  | optional        | optional       |
| log_level                   | optional  | optional        | optional       |
| log_stream                  | optional  | optional        | optional       |
| log_format                  | optional  | optional        | optional       |
| var_river_id                | optional  | optional        | optional       |
| var_discharge               | optional  | optional        | optional       |
| var_x                       | optional  | optional        | optional       |
| var_y                       | optional  | optional        | optional       |
| var_t                       | optional  | optional        | optional       |
| var_catchment_volume        | optional  | optional        | optional       |
| var_runoff_depth            | optional  | optional        | optional       |

Required = required, optional = optional, blank = not used.
