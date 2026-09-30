## Configuration File

`river-route` computations are controlled by a `Configs` object, built from keyword arguments or read from a
JSON file with `Configs.from_json`.
All routing runs through `Router`. The procedure it runs is set by the selector keys (`coefficients`, `forcing`,
`transform`, `network_type`), and the required config keys depend on which selections you make.

`Router` takes a `Configs`, and optionally the `Network` and the Runoff it would otherwise build for itself:
`Router(configs, network=..., runoff=...)`. `Network` and the Runoff classes take the options they need as ordinary
arguments and each has a `from_configs` classmethod that reads those same values off a `Configs`; `Router` builds
them that way when it is not given them. `examples/config.json` below lists every option.

### Routing procedure selectors

- `coefficients` - `'static'` (constant Muskingum K from columns `muskingumK`, `muskingumX`) or `'dynamic'`
  (nonlinear K = dynamicAlpha\*Q^dynamicBeta from columns `dynamicAlpha`, `dynamicBeta`, `muskingumX`). Default
  `'static'`.
- `forcing` - `'channel'` (channel routing only, no inflows), or the form of the `runoff_files` whose runoff enters the
  rivers in addition to routing: `'catchment'` (already aggregated to catchments, read by `CatchmentRunoff`),
  `'grid'` (a grid with x and y dimensions, read by `GridRunoff`), or `'ecmwf_grib'` (ECMWF GRIB
  files on a reduced gaussian grid, read by `ECMWFGribReducedGrid`). A single value. Default `'channel'`.
- `transform` - `'uniform'`, the only option: each step's catchment runoff enters its river at a constant rate over
  the step. Default `'uniform'`.
- `network_type` - `'standard'` (one reach per river) or `'stabilized'` (each river too long for
  `dt_routing` is routed in substeps, the fewest equal sub-reaches in series that are each Muskingum-stable, and each
  river too short for it is routed in subcycles, the fewest equal steps of its own that are). Default `'standard'`.
  `'stabilized'` needs static coefficients, and multiplies the routing work by the average
  number of substeps and subcycles per river. A river routed in subcycles interpolates its upstream inflow linearly
  within each routing step. Without `dt_routing` it routes at the largest divisor of `dt_runoff` at which every river can
  be made stable. The state files then hold one value per sub-reach, each river's sub-reaches upstream to
  downstream, in the order a final state file is written. A state file with one value per river is refused, and so
  is a run whose sub-reaches change between its runoff files, which setting `dt_routing` prevents.

The selectors together choose the routing method and the form of runoff it reads; see [kernels](kernels.md). The one
combination no routing method routes yet, `'dynamic'` coefficients on a `'stabilized'` network, raises
`NotImplementedError` before any runoff is read.

## Minimum Required Inputs

Every routing procedure requires the following 2 configuration options:

- `network_file` - path to the [network file](io-file-schema.md#network-file) (parquet)
- One of two options for specifying where the [routed discharge](io-file-schema.md#routed-discharge) output is written
    - `discharge_dir` - a string path to a directory where outputs are saved based on the names of the inputs. Each
      output is named `discharge_<input name>.zarr`, or `discharge.zarr` for channel routing, the zarr store that the
      default writer writes. Give `discharge_files` to name outputs for another writer.
    - `discharge_files` - list of explicit paths for each output file, one per input file required.

## Required config keys by selection

Beyond the always-required keys above, additional keys are required depending on your selector choices.

**`forcing: channel` (channel routing only, no inflows)** also requires:

- `channel_state_init_file` - parquet state file to initialize discharge
- `dt_routing` - routing timestep in seconds
- `dt_total` - total simulation duration in seconds

**`forcing: catchment`, `grid`, or `ecmwf_grib`** also requires a water input source:

- `runoff_files`
- `grid_weights_file` when `forcing` is `grid` or `ecmwf_grib`. It must not be set for
  `catchment`.

  Time keys for forced procedures (`dt_total`, `dt_discharge`, `dt_runoff`, `dt_routing`) are resolved from the
  inputs where possible; see the [time options](time-options.md). `start_datetime` is only read by channel routing,
  which has no input dates to copy.

**`coefficients` selection** determines the required `network_file` columns:

- `coefficients: static` requires columns `muskingumK`, `muskingumX`
- `coefficients: dynamic` also requires columns `dynamicAlpha`, `dynamicBeta` (the K formula uses `dynamicAlpha`,
  `dynamicBeta`, `muskingumX`, but a `muskingumK` column must still be present)

Every network file also has the columns `riverId`, `nextRiverId`, `riverIndex`, and `upstreamCount`; see the
[File Schemas reference](io-file-schema.md#network-file).

The following table lists where each remaining key applies.

| Config key                 | Description                      | Required when                                          |
|----------------------------|----------------------------------|--------------------------------------------------------|
| **core**                   |                                  |                                                        |
| `network_file`              | Network file parquet.      | always                                                 |
| **state**                  |                                  |                                                        |
| `channel_state_init_file`  | Path to initial channel state    | `forcing: channel` (else optional, default 0)          |
| `channel_state_final_file` | Path to save final channel state | optional                                               |
| **output**                 |                                  |                                                        |
| `discharge_dir`            | Directory for output  files      | _Option 1_                                             |
| `discharge_files`          | Explicit output paths            | _Option 2_                                             |
| **input data**             |                                  |                                                        |
| `runoff_files`             | Runoff read as the `forcing`     | `forcing` not `channel`                                |
| `grid_weights_file`        | Aggregates grids to catchments   | `forcing` a grid                                       |
| **time**                   |                                  |                                                        |
| `start_datetime`           | Channel routing start date       | optional, `forcing: channel` only                      |
| `dt_total`                 | Total simulation duration        | `forcing: channel` (else [time docs](time-options.md)) |
| `dt_discharge`             | Output timestep                  | optional - [time docs](time-options.md)                |
| `dt_runoff`                | Runoff data timestep             | optional - [time docs](time-options.md)                |
| `dt_routing`               | Routing computational timestep   | `forcing: channel` (else [time docs](time-options.md)) |

## Optional configs with defaults

| Config Key               | Description                                            | Default                                       |
|--------------------------|--------------------------------------------------------|-----------------------------------------------|
| `coefficients`           | Muskingum K source: `'static'` or `'dynamic'`          | `'static'`                                    |
| `forcing`                | Channel routing, or the form of the runoff files       | `'channel'`                                   |
| `transform`              | Runoff transform: `'uniform'`                          | `'uniform'`                                   |
| `network_type`           | Reach handling: `'standard'` or `'stabilized'`         | `'standard'`                                  |
| `log`                    | Enable or disable logging                              | `True`                                        |
| `progress_bar`           | Show tqdm progress bar                                 | `True`                                        |
| `log_level`              | Logger level, defaults to between INFO and WARNING     | `'PROGRESS'`                                  |
| `log_stream`             | `'stdout'` or a file path                              | `'stdout'`                                    |
| `log_format`             | Python logging format string                           | `'%(levelname)s - %(asctime)s - %(message)s'` |
| `var_river_id`           | River ID dimension name in files                       | `'river_id'`                                  |
| `var_discharge`          | Discharge variable name in output                      | `'Q'`                                         |
| `var_grid_runoff`        | Runoff variable name in grid `runoff_files`            | `'ro'`                                        |
| `var_x`                  | X-dimension name in `grid` runoff files                | `'x'`                                         |
| `var_y`                  | Y-dimension name in `grid` runoff files                | `'y'`                                         |
| `var_cell`               | Cell dimension name in reduced gaussian grids          | `'cell'`                                      |
| `var_t`                  | Time dimension name in depth grids                     | `'time'`                                      |
| `grid_accumulation_type` | Is runoff grid `'incremental'` or `'cumulative'`       | `'incremental'`                               |
| `runoff_processing_mode` | Are runoff `'sequential'` or `'ensemble'` inputs       | `'sequential'`                                |
| `runoff_depth_unit`      | Unit of grid runoff depths, else read from the file    | `None`                                        |
| `force_positive_runoff`  | Clip negative grid runoff depths to zero               | `False`                                       |
| `force_uniform_timesteps` | Resample irregular grid runoff to the first timestep  | `True`                                        |
| `as_volumes`             | `GridRunoff` prepares volumes instead of depths        | `False`                                       |
| `unstable_coefficients`  | `'warn'`, `'raise'`, or `'ignore'` unstable rivers     | `'warn'`                                      |

## Validation

`Router.route()` validates the configs with `Configs.validate_routing` before it computes anything, and
`GridRunoff.from_configs` validates them with `Configs.validate_runoff` before it reads the weight table. Configs are
frozen, so once they pass they are not checked again. Passing one of the two does not pass the other. A `GridRunoff`
built directly, without a `Configs`, has nothing to validate and so runs neither.

`Configs.deep_validate()` reads the network file, grid weights, and initial state and checks their columns,
types, and value ranges, that the network is topologically sorted, and that the weight table proportions sum to
1 per river. Nothing calls it for you, because reading those files is the same work routing is about to do. Run it
once on inputs you have not checked before:

```python
rr.Configs.from_json('config.json').deep_validate()
```

`unstable_coefficients` controls what happens when a river's parameters are not Muskingum-stable for the
routing timestep, which requires `2*k*x <= dt_routing <= 2*k*(1-x)`. Outside that window the solution for
that river oscillates and negative discharges are clamped to zero, which does not conserve mass. The
default `'warn'` issues a warning saying how many rivers are affected; `'raise'` refuses to route; `'ignore'` is
silent. It applies to static coefficients, which the Router checks by calling `Network.check_stability`. Dynamic
coefficients change as each river routes, so there is nothing to check before routing starts.

Use `Network.unstable_mask(dt)` to inspect a network before routing it. `network_type='stabilized'` routes each river
too long for `dt_routing` in `Network.substeps(dt)` sub-reaches inside the kernel, and each river too short for it,
which substeps cannot fix, in `Network.subcycles(dt)` steps of its own. `Network.stabilize(dt)` instead rewrites the
network in place with those sub-reaches as rows of their own, for the rivers substeps can fix, and
`Network.write_stabilized(dt)` saves that network as a network table.

## Example Configuration

The general template lists every key.

```json title="config.json"
--8<-- "examples/config.json"
```
