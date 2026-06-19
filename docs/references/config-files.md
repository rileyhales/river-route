## Configuration File

`river-route` computations are controlled by config values passed as keyword arguments or from a YAML/JSON file.
There is a single `Router`. The routing procedure it runs is set by three selector keys (`coeff`, `forcing`,
`network`), and the required config keys depend on which selections you make.

### Routing procedure selectors

- `coeff` - `'static'` (constant Muskingum K from columns `k`, `x`) or `'dynamic'` (nonlinear K = alpha\*Q^beta
  from columns `alpha`, `beta`, `x`). Default `'static'`.
- `forcing` - one of `'channel'` (channel routing only, no inflows), `'lateral'` (lateral runoff inflow), or
  `'external'` (planned, not yet available). A single value; specifying multiple forcings is planned for later. Default `'channel'`.
- `network` - `'standard'` (one reach per river). `'expanded'` (automatically subdivide and substep unstable
  reaches for numerical stability) is planned and not yet available. Default `'standard'`.

## Minimum Required Inputs

Every routing procedure requires the following 2 configuration options:

- `params_file` - path to the [routing parameters file](io-file-schema.md#routing-parameters) (parquet)
- One of two options for specifying where the [routed discharge](io-file-schema.md#routed-discharge) output is written
    - `discharge_dir` - a string path to a directory where outputs are saved based on the names of the inputs
    - `discharge_files` - list of explicit paths for each output file, one per input file required.

## Required config keys by selection

Beyond the always-required keys above, additional keys are required depending on your selector choices.

**`forcing: channel` (channel routing only, no inflows)** also requires:

- `channel_state_init_file` - parquet state file to initialize discharge
- `dt_routing` - routing timestep in seconds
- `dt_total` - total simulation duration in seconds

**`forcing: lateral`** also requires a water input source:

- `qlateral_files`, or
- `grid_runoff_files` plus `grid_weights_file`

  Time keys for forced procedures (`dt_total`, `dt_discharge`, `dt_runoff`, `dt_routing`, `start_datetime`)
  are resolved from the inputs where possible; see the [time options](time-options.md).

**`coeff` selection** determines the required `params_file` columns:

- `coeff: static` requires columns `k`, `x`
- `coeff: dynamic` requires columns `k`, `x`, `alpha`, `beta` (the K formula uses `alpha`, `beta`, `x`, but a `k` column must still be present)

The following table lists where each remaining key applies.

| Config key                 | Description                      | Required when                                       |
|----------------------------|----------------------------------|-----------------------------------------------------|
| **core**                   |                                  |                                                     |
| `params_file`              | Routing parameters parquet.      | always                                              |
| **state**                  |                                  |                                                     |
| `channel_state_init_file`  | Path to initial channel state    | `forcing: channel` (else optional, default 0)          |
| `channel_state_final_file` | Path to save final channel state | optional                                            |
| **output**                 |                                  |                                                     |
| `discharge_dir`            | Directory for output  files      | _Option 1_                                          |
| `discharge_files`          | Explicit output paths            | _Option 2_                                          |
| **input data**             |                                  |                                                     |
| `qlateral_files`           | Per-catchment runoff time series | `forcing: lateral`, _Option 1_                      |
| `grid_runoff_files`        | Gridded runoff depths            | `forcing: lateral`, _Option 2_                      |
| `grid_weights_file`        | Converts depth grids to qlateral | `forcing: lateral`, _Option 2_                      |
| **time**                   |                                  |                                                     |
| `start_datetime`           | Simulation start date            | optional                                            |
| `dt_total`                 | Total simulation duration        | `forcing: channel` (else [time docs](time-options.md)) |
| `dt_discharge`             | Output timestep                  | optional - [time docs](time-options.md)             |
| `dt_runoff`                | Runoff data timestep             | optional - [time docs](time-options.md)             |
| `dt_routing`               | Routing computational timestep   | `forcing: channel` (else [time docs](time-options.md)) |

## Optional configs with defaults

| Config Key               | Description                                        | Default                                       |
|--------------------------|----------------------------------------------------|-----------------------------------------------|
| `coeff`                  | Muskingum K source: `'static'` or `'dynamic'`      | `'static'`                                    |
| `forcing`                | Inflow forcing: `'channel'`, `'lateral'`, `'external'` | `'channel'`                                   |
| `network`                | Reach handling: `'standard'` (`'expanded'` planned) | `'standard'`                                  |
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

The general template in YAML format lists all keys with comments. Template files for specific
procedures are available in the examples directory `config_muskingum.yaml`, `config_rapid_muskingum.yaml`, `config_unit_muskingum.yaml`.

### General Config File

```yaml title="config.yaml"
{% include-markdown "../../examples/config.yaml" %}
```

### Channel-only

```yaml title="config_muskingum.yaml"
{% include-markdown "../../examples/config_muskingum.yaml" %}
```

### Lateral runoff (Rapid-style)

```yaml title="config_rapid_muskingum.yaml"
{% include-markdown "../../examples/config_rapid_muskingum.yaml" %}
```

### Unit hydrograph routing (planned)

!!! warning "Not yet available in v3"
    Unit hydrograph routing is planned for a later v3 release and is not yet available. The example
    below is retained for reference.

```yaml title="config_unit_muskingum.yaml"
{% include-markdown "../../examples/config_unit_muskingum.yaml" %}
```
