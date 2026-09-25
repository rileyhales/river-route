## Configuration File

`river-route` computations are controlled by a `Configs` object, built from keyword arguments or read from a
JSON file with `Configs.from_json`.
All routing runs through `Router`. The procedure it runs is set by the selector keys (`coeff`, `forcing`,
`transform`, `runoff_type`, `network_type`), and the required config keys depend on which selections you make.

`Router` takes a `Configs` and nothing else. `Network` and the Runoff classes take the options they need as ordinary
arguments and each has a `from_configs` classmethod that reads those same values off a `Configs`; `Router` builds
them that way. `examples/config.json` below lists every option.

### Routing procedure selectors

- `coeff` - `'static'` (constant Muskingum K from columns `k`, `x`) or `'dynamic'` (nonlinear K = alpha\*Q^beta
  from columns `alpha`, `beta`, `x`). Default `'static'`.
- `forcing` - one of `'channel'` (channel routing only, no inflows) or `'runoff'` (runoff enters the rivers in addition to routing).
  A single value. Default `'channel'`.
- `transform` - `'uniform'` or `'unit_hydrograph'`, the runoff transformation applied under lateral forcing.
  Only read when `forcing` is `'runoff'`. Default `'uniform'`.
- `runoff_type` - the form of `runoff_files`: `'catchment'` (already aggregated to catchments, read by
  `CatchmentRunoff`), `'gaussian_grid'` (a grid with x and y dimensions, read by `GaussianGridRunoff`), or
  `'reduced_gaussian_grid'` (a grid with one cell dimension, read by `ReducedGaussianGridRunoff`, not implemented yet).
  Required when `forcing` is `'runoff'`, with no default.
- `network_type` - `'standard'` (one reach per river) or `'stabilized'` (each river too long for
  `dt_routing` is routed as the fewest equal sub-reaches in series that are each Muskingum-stable, and each river
  too short for it is sub-cycled in the fewest equal steps of its own that are). Default `'standard'`.
  `'stabilized'` needs static coefficients, and multiplies the routing work by the average
  number of sub-reaches and substeps per river. A sub-cycled river interpolates its upstream inflow linearly within
  each routing step. Without `dt_routing` it routes at the largest divisor of `dt_runoff` at which every river can
  be made stable. The state files then hold one value per sub-reach, each
  river's sub-reaches upstream to downstream; a state file with one value per river seeds all of its sub-reaches.

The selectors together name one kernel; see [kernels](kernels.md). A combination with no kernel yet raises
`NotImplementedError` naming it and listing those that exist.

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

**`forcing: runoff`** also requires a water input source:

- `runoff_files` and `runoff_type`
- `grid_weights_file` when `runoff_type` is `gaussian_grid` or `reduced_gaussian_grid`. It must not be set for
  `catchment`.

  Time keys for forced procedures (`dt_total`, `dt_discharge`, `dt_runoff`, `dt_routing`, `start_datetime`)
  are resolved from the inputs where possible; see the [time options](time-options.md).

**`coeff` selection** determines the required `params_file` columns:

- `coeff: static` requires columns `k`, `x`
- `coeff: dynamic` requires columns `k`, `x`, `alpha`, `beta` (the K formula uses `alpha`, `beta`, `x`, but a `k` column must still be present)

The following table lists where each remaining key applies.

| Config key                 | Description                      | Required when                                          |
|----------------------------|----------------------------------|--------------------------------------------------------|
| **core**                   |                                  |                                                        |
| `params_file`              | Routing parameters parquet.      | always                                                 |
| **state**                  |                                  |                                                        |
| `channel_state_init_file`  | Path to initial channel state    | `forcing: channel` (else optional, default 0)          |
| `channel_state_final_file` | Path to save final channel state | optional                                               |
| **output**                 |                                  |                                                        |
| `discharge_dir`            | Directory for output  files      | _Option 1_                                             |
| `discharge_files`          | Explicit output paths            | _Option 2_                                             |
| **input data**             |                                  |                                                        |
| `runoff_files`             | Runoff read as the `runoff_type` | `forcing: runoff`                                      |
| `grid_weights_file`        | Aggregates grids to catchments   | `runoff_type` a grid                                   |
| **unit hydrograph**        |                                  |                                                        |
| `uh_kernel_file`           | Unit hydrograph kernel (npz)     | `transform: unit_hydrograph`                           |
| `uh_state_init_file`       | Initial unit hydrograph state    | optional                                               |
| `uh_state_final_file`      | Path to save final UH state      | optional                                               |
| **time**                   |                                  |                                                        |
| `start_datetime`           | Simulation start date            | optional                                               |
| `dt_total`                 | Total simulation duration        | `forcing: channel` (else [time docs](time-options.md)) |
| `dt_discharge`             | Output timestep                  | optional - [time docs](time-options.md)                |
| `dt_runoff`                | Runoff data timestep             | optional - [time docs](time-options.md)                |
| `dt_routing`               | Routing computational timestep   | `forcing: channel` (else [time docs](time-options.md)) |

## Optional configs with defaults

| Config Key               | Description                                            | Default                                       |
|--------------------------|--------------------------------------------------------|-----------------------------------------------|
| `coeff`                  | Muskingum K source: `'static'` or `'dynamic'`          | `'static'`                                    |
| `forcing`                | Inflow forcing: `'channel'`, `'runoff'`                | `'channel'`                                   |
| `transform`              | Runoff transform: `'uniform'`, `'unit_hydrograph'`     | `'uniform'`                                   |
| `network_type`           | Reach handling: `'standard'` or `'stabilized'`         | `'standard'`                                  |
| `log`                    | Enable or disable logging                              | `True`                                        |
| `progress_bar`           | Show tqdm progress bar                                 | `True`                                        |
| `log_level`              | Logger level, defaults to between INFO and WARNING     | `'PROGRESS'`                                  |
| `log_stream`             | `'stdout'` or a file path                              | `'stdout'`                                    |
| `log_format`             | Python logging format string                           | `'%(levelname)s - %(asctime)s - %(message)s'` |
| `var_river_id`           | River ID dimension name in files                       | `'river_id'`                                  |
| `var_discharge`          | Discharge variable name in output                      | `'Q'`                                         |
| `var_grid_runoff`        | Runoff variable name in grid `runoff_files`            | `'ro'`                                        |
| `var_x`                  | X-dimension name in gaussian grids                     | `'x'`                                         |
| `var_y`                  | Y-dimension name in gaussian grids                     | `'y'`                                         |
| `var_cell`               | Cell dimension name in reduced gaussian grids          | `'cell'`                                      |
| `var_t`                  | Time dimension name in depth grids                     | `'time'`                                      |
| `grid_accumulation_type` | Is runoff grid `'incremental'` or `'cumulative'`       | `'incremental'`                               |
| `runoff_processing_mode` | Are runoff `'sequential'` or `'ensemble'` inputs       | `'sequential'`                                |
| `runoff_depth_unit`      | Unit of grid runoff depths, else read from the file    | `None`                                        |
| `force_positive_runoff`  | Clip negative grid runoff depths to zero               | `False`                                       |
| `force_uniform_timesteps` | Resample irregular grid runoff to the first timestep  | `True`                                        |
| `as_volumes`             | `GaussianGridRunoff` prepares volumes instead of depths            | `False`                                       |
| `unstable_coefficients`  | `'warn'`, `'raise'`, or `'ignore'` unstable rivers     | `'warn'`                                      |

## Validation

`Router.route()` validates the configs with `Configs.validate_routing` before it computes anything, and
`GaussianGridRunoff.from_configs` validates them with `Configs.validate_runoff` before it reads the weight table. Configs are
frozen, so once they pass they are not checked again. A `GaussianGridRunoff` built directly, without a `Configs`, has nothing
to validate and so runs neither.

`Configs.deep_validate()` reads the params file, grid weights, and initial state and checks their columns,
types, and value ranges, that the network is topologically sorted, and that the weight table proportions sum to
1 per river. Nothing calls it for you, because it reads every input file, which is the same work routing is
about to do. Run it once on inputs you have not checked before:

```python
rr.Configs.from_json('config.json').deep_validate()
```

`unstable_coefficients` controls what happens when a river's parameters are not Muskingum-stable for the
routing timestep, which requires `2*k*x <= dt_routing <= 2*k*(1-x)`. Outside that window the solution for
that river oscillates and negative discharges are clamped to zero, which does not conserve mass. The
default `'warn'` logs how many rivers are affected; `'raise'` refuses to route; `'ignore'` is silent. `Router`
applies it by calling `Network.check_stability`. Use `Network.unstable_mask(dt)` to inspect a network before
routing it, and `Network.stabilize(dt)` to build the stabilized network, which adds sub-reaches until every
one routes stably at that dt. `network_conditioning='stabilized'` routes that network, and sub-cycles the rivers
too short for `dt_routing`, which splitting cannot fix, in `Network.substeps(dt)` steps each.

## Example Configuration

The general template lists every key.

```json title="config.json"
--8<-- "examples/config.json"
```
