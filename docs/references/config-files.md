## Configuration File

`river-route` computations are controlled by config values passed as keyword arguments or from a YAML/JSON file.
The three routers (`Muskingum`, `RapidMuskingum`, `UnitMuskingum`) share a set of base config keys and each
adds router-specific keys.

## Minimum Required Inputs

All routing classes require the following 2 configuration options:

- `params_file` - path to the [routing parameters file](io-file-schema.md#routing-parameters) (parquet)
- `discharge_dir` - directory where [routed discharge](io-file-schema.md#routed-discharge) output files are written (netCDF, named after input files)

`Muskingum` (channel routing only, no lateral inflows) also requires:

- `channel_state_init_file` - parquet state file to initialize discharge
- `dt_routing` - routing timestep in seconds
- `dt_total` - total simulation duration in seconds

`RapidMuskingum` also requires:

- one water input source:
    - `qlateral_files`, or
    - `grid_runoff_files` plus `grid_weights_file`

`UnitMuskingum` also requires:

- `uh_kernel_file` - pre-computed scipy sparse npz kernel
- one water input source:
    - `qlateral_files`, or
    - `grid_runoff_files` plus `grid_weights_file`

## Required Config Keys

| Config key                 | Description                      |  Muskingum   |             RapidMuskingum              |              UnitMuskingum              |
|----------------------------|----------------------------------|:------------:|:---------------------------------------:|:---------------------------------------:|
| **core**                   |                                  |              |                                         |                                         |
| `params_file`              | Routing parameters parquet.      | **Required** |              **Required**               |              **Required**               |
| **state**                  |                                  |              |                                         |                                         |
| `channel_state_init_file`  | Path to initial channel state    | **Required** |         optional - default to 0         |         optional - default to 0         |
| `channel_state_final_file` | Path to save final channel state |   optional   |                optional                 |                optional                 |
| **output**                 |                                  |              |                                         |                                         |
| `discharge_dir`            | Directory for output  files      |  _Option 1_  |               _Option 1_                |               _Option 1_                |
| `discharge_files`          | Explicit output paths            |  _Option 2_  |               _Option 2_                |               _Option 2_                |
| **input data**             |                                  |              |                                         |                                         |
| `qlateral_files`           | Per-catchment runoff time series |              |               _Option 1_                |               _Option 1_                |
| `grid_runoff_files`        | Gridded runoff depths            |              |               _Option 2_                |               _Option 2_                |
| `grid_weights_file`        | Converts depth grids to qlateral |              |               _Option 2_                |               _Option 2_                |
| **unit hydrograph**        |                                  |              |                                         |                                         |
| `uh_kernel_file`           | Pre-computed convolution kernel  |              |                                         |              **Required**               |
| `uh_state_init_file`       | Parquet with initial UH state    |              |                                         |                optional                 |
| `uh_state_final_file`      | Path to save final UH state      |              |                                         |                optional                 |
| **time**                   |                                  |              |                                         |                                         |
| `start_datetime`           | Simulation start date            |   optional   |                                         |                                         |
| `dt_total`                 | Total simulation duration        | **Required** | optional - [time docs](time-options.md) | optional - [time docs](time-options.md) |
| `dt_discharge`             | Output timestep                  |   optional   | optional - [time docs](time-options.md) | optional - [time docs](time-options.md) |
| `dt_runoff`                | Runoff data timestep             |              | optional - [time docs](time-options.md) | optional - [time docs](time-options.md) |
| `dt_routing`               | Routing computational timestep   | **Required** | optional - [time docs](time-options.md) | optional - [time docs](time-options.md) |

## Optional Configs with Defaults

| Config Key               | Description                                        | Default                                       |
|--------------------------|----------------------------------------------------|-----------------------------------------------|
| `log`                    | Enable or disable logging                          | `True`                                        |
| `progress_bar`           | Show tqdm progress bar                             | `True`                                        |
| `log_level`              | Logger level, defaults to between INFO and WARNING | `'PROGRESS'`                                  |
| `log_stream`             | `'stdout'` or a file path                          | `'stdout'`                                    |
| `log_format`             | Python logging format string                       | `'%(levelname)s - %(asctime)s - %(message)s'` |
| `var_river_id`           | River ID dimension name in files                   | `'river_id'`                                  |
| `var_discharge`          | Discharge variable name in output                  | `'Q'`                                         |
| `var_grid_runoff`        | Runoff variable name in `grid_runoff_files`        | `'ro'`                                        |
| `var_x`                  | X-dimension name in depth grids                    | `'x'`                                         |
| `var_y`                  | Y-dimension name in depth grids                    | `'y'`                                         |
| `var_t`                  | Time dimension name in depth grids                 | `'time'`                                      |
| `grid_accumulation_type` | Is runoff grid `'incremental'` or `'cumulative'`   | `'incremental'`                               |
| `runoff_processing_mode` | Are runoff `'sequential'` or `'ensemble'` inputs   | `'sequential'`                                |

## Example Configuration YAMLs

The general template in YAML format lists all keys with comments. Template files for specific routers are available in the examples directory
`config_muskingum.yaml`, `config_rapid_muskingum.yaml`, `config_unit_muskingum.yaml`.

### General Config File

```yaml title="config.yaml"
{% include-markdown "../../examples/config.yaml" %}
```

### Muskingum Config File

```yaml title="config_muskingum.yaml"
{% include-markdown "../../examples/config_muskingum.yaml" %}
```

### Rapid Muskingum Config File

```yaml title="config_rapid_muskingum.yaml"
{% include-markdown "../../examples/config_rapid_muskingum.yaml" %}
```

### Unit Muskingum Config File

```yaml title="config_unit_muskingum.yaml"
{% include-markdown "../../examples/config_unit_muskingum.yaml" %}
```
