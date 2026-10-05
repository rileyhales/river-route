## Configuration File

`river-route` computations are controlled by a `Configs` object, built from keyword arguments or read from a
JSON file with `Configs.from_json`.
All routing runs through `Router`. The procedure it runs is set by the selector keys (`parameterization`, `forcing`,
`transform`, `force_validity`), and the required config keys depend on which selections you make.

`Router` takes only a `Configs`. `Network` and the Runoff classes take the options they need as ordinary arguments
and each has a `from_configs` classmethod that reads those same values off a `Configs`, which is how `Router` builds
them. `examples/config.json` below lists every option.

### Routing procedure selectors

- `parameterization` - how each river's Muskingum k and x are set, which also chooses the routing method:
  `'static'` builds constant coefficients from columns `muskingumK`, `muskingumX`, and `'power_law'` recomputes k from
  the river's discharge every routing step (variable parameter Muskingum), k = powerLawAlpha\*Q^powerLawBeta from
  columns `powerLawAlpha`, `powerLawBeta`, with column `muskingumX`. Default `'static'`.
- `forcing` - `'channel'` (channel routing only, no inflows), or the form of the `runoff_files` whose runoff enters the
  rivers in addition to routing: `'catchment'` (already aggregated to catchments, read by `CatchmentRunoff`),
  `'grid'` (a grid with x and y dimensions, read by `GridRunoff`), or `'ecmwf_grib'` (ECMWF GRIB
  files on a reduced gaussian grid, read by `ECMWFGribReducedGrid`). A single value. Default `'channel'`.
- `transform` - `'uniform'`, the only option: each step's catchment runoff enters its river at a constant rate over
  the step. Default `'uniform'`.
- `force_validity` - `True` (the default) makes every Muskingum coefficient non-negative, the condition for valid
  routing: each river too long for `dt_routing` is routed in substeps, the fewest equal sub-reaches in series whose
  coefficients are all non-negative, each with the x that keeps the river's diffusion, and each river too short for
  it is routed in subcycles, the fewest equal steps of its own whose coefficients are. `False` routes each river whole
  as given. `True` needs the `'static'` parameterization, and multiplies the routing work by the average
  number of substeps and subcycles per river. A river routed in subcycles interpolates its upstream inflow linearly
  within each routing step. Without `dt_routing` it routes at the largest divisor of `dt_runoff` at which every
  coefficient of every river can be made non-negative, which for `x <= 1/3` is `dt_runoff` itself. The state files then hold one value per sub-reach, each river's sub-reaches upstream to
  downstream, in the order a final state file is written, with each row's `riverId` repeated for every sub-reach of
  its river. A state file with one value per river is refused, and so is a run whose sub-reaches change between its
  runoff files, which setting `dt_routing` prevents.

The selectors together choose the routing method and the form of runoff it reads. The one combination no routing
method routes yet, a parameterization other than `'static'` with `force_validity`, raises `NotImplementedError` before
any runoff is read, so route `'power_law'` with `force_validity=False`.

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

  Time keys for forced procedures (`dt_total`, `dt_discharge`, `dt_routing`) are resolved from the inputs where
  possible, and the runoff time step is always read from the runoff files; see the
  [time options](time-options.md). `start_datetime` is only read by channel routing, which has no input dates to copy.

**`parameterization` selection** determines the required `network_file` columns:

- `parameterization: static` requires columns `muskingumK`, `muskingumX`
- `parameterization: power_law` also requires columns `powerLawAlpha`, `powerLawBeta` (k is recomputed from them and
  `muskingumX` is read, but a `muskingumK` column must still be present)

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
| `dt_routing`               | Routing computational timestep   | `forcing: channel` (else [time docs](time-options.md)) |

## Optional configs with defaults

The `var_` options name the variables of the gridded runoff files that are read. Every other file, read or written,
uses fixed names, such as `riverId` for the river id.

| Config Key               | Description                                              | Default                                       |
|--------------------------|----------------------------------------------------------|-----------------------------------------------|
| `parameterization`       | Muskingum k and x: `'static'` or `'power_law'`           | `'static'`                                    |
| `forcing`                | Channel routing, or the form of the runoff files         | `'channel'`                                   |
| `transform`              | Runoff transform: `'uniform'`                            | `'uniform'`                                   |
| `force_validity`         | Route in substeps and subcycles so coefficients are >= 0 | `True`                                        |
| `log`                    | Enable or disable logging                                | `True`                                        |
| `progress_bar`           | Show tqdm progress bar                                   | `True`                                        |
| `log_level`              | Logger level, defaults to between INFO and WARNING       | `'PROGRESS'`                                  |
| `log_stream`             | `'stdout'` or a file path                                | `'stdout'`                                    |
| `log_format`             | Python logging format string                             | `'%(levelname)s - %(asctime)s - %(message)s'` |
| `var_grid_runoff`        | Runoff variable name in grid `runoff_files`              | `'ro'`                                        |
| `var_x`                  | X-dimension name in `grid` runoff files                  | `'x'`                                         |
| `var_y`                  | Y-dimension name in `grid` runoff files                  | `'y'`                                         |
| `var_t`                  | Time dimension name in depth grids                       | `'time'`                                      |
| `grid_accumulation_type` | Is runoff grid `'incremental'` or `'cumulative'`         | `'incremental'`                               |
| `runoff_processing_mode` | Are runoff `'sequential'` or `'ensemble'` inputs         | `'sequential'`                                |
| `runoff_depth_unit`      | Unit of grid runoff depths, else read from the file      | `None`                                        |
| `force_positive_runoff`  | Clip negative grid runoff depths to zero                 | `False`                                       |
| `as_volumes`             | `GridRunoff` prepares volumes instead of depths          | `False`                                       |
| `invalid_coefficients`   | `'warn'`, `'raise'`, or `'ignore'` negative coefficients | `'warn'`                                      |

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

`invalid_coefficients` controls what happens when a river has a negative Muskingum coefficient at the routing
timestep. Every coefficient is non-negative, and the routing valid, when `2*k*x <= dt_routing <= 2*k*(1-x)`. Outside
that window routing is still stable and still conserves mass, but the river's discharge can dip below zero ahead of a rise (`c1 < 0`, a river
too long for the step) or alternate from step to step and go negative (`c3 < 0`, a river too short for it). Negative
discharge is written as routed, never clamped, so the written discharge conserves mass and shows the river that
produced it. The default `'warn'` issues a warning saying how many rivers are affected; `'raise'` refuses to route;
`'ignore'` is silent. It applies to the `'static'` parameterization, which routing checks with
`Network.check_validity`. The coefficients of the other parameterizations change as each river routes, so there is
nothing to check before routing starts.

Use `Network.invalid_mask(dt)` to inspect a network before routing it. `force_validity=True` routes each river too
long for `dt_routing` in `Network.substeps(dt)` sub-reaches inside the kernel, each with the x of
`Network.substep_x`, `1/2 - N*(1/2 - x)` for N sub-reaches, which keeps the river's mean travel time and its travel time
variance `k**2 * (1 - 2x)`. At `x = 0.2` two sub-reaches make every river too long for any step non-negative. Each
river too short for `dt_routing`, which substeps cannot fix, is routed in `Network.subcycles(dt)` steps of its own. `Network.stabilize(dt)` instead rewrites the
network in place with those sub-reaches as rows of their own, for the rivers substeps can fix, and
`Network.write_stabilized(dt)` saves that network as a network table. Such a network routes with the `channel`,
`grid`, and `ecmwf_grib` forcings, whose runoff is shared among the sub-reaches of each river; `catchment` forcing
refuses it, since a catchment runoff file has one row per river.

## Example Configuration

The general template lists every key.

```json title="config.json"
--8<-- "examples/config.json"
```
