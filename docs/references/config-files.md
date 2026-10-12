## Configuration File

`river-route` computations are controlled by a `Configs` object, built from keyword arguments or read from a
JSON file with `Configs.from_json`.
All routing runs through `Router`, which routes Muskingum with static coefficients: each river's k and x are fixed for
the whole run. The procedure it runs is set by the `forcing` selector and the compute options (`dt_routing`,
`conditioning`), and the required config keys depend on which selections you make.

`Router` takes only a `Configs`. `Network` and the Runoff classes take the options they need as ordinary arguments
and each has a `from_configs` classmethod that reads those same values off a `Configs`, which is how `Router` builds
them. `examples/config.json` below lists every option.

### Routing procedure selector

- `forcing` - `'channel'` (channel routing only, no inflows), or the form of the `forcing_files` whose runoff enters
  the rivers in addition to routing: `'catchment'` (already aggregated to catchments, read by `CatchmentRunoff`),
  `'grid'` (a grid with x and y dimensions, read by `GridRunoff`), or `'ecmwf_grib'` (ECMWF GRIB
  files on a reduced gaussian grid, read by `ECMWFGribReducedGrid`). A single value. Default `'channel'`. Each step's
  catchment runoff enters its river at a constant rate over the step.

Each river's coefficients are built from the `muskingumK` and `muskingumX` columns of the network file.

### Compute options

`dt_routing` and `conditioning` settle how each river computes, in two jobs: the first sets each river's step, and the
second makes every coefficient of each river non-negative at it, the condition for valid routing.

- `dt_routing` - the step each river computes in. `'auto'` (the default) gives each river its own valid step: the
  runoff step over the fewest subcycles, equal steps of its own, that keep its coefficients non-negative. A step in
  seconds, a divisor of the runoff step, is instead the step of every river. Setting it to the runoff step routes each
  river in one step of each runoff step.
- `conditioning` - `True` (the default) splits each river too long for its step into the fewest sub-reaches in series
  whose coefficients are all non-negative, which keep the river's mean travel time and its diffusion, and routes each
  river too short for its step in subcycles of it, the fewest whose coefficients are. `False` routes every river whole
  at its step: no river is split or given a shorter step, and a river with a negative coefficient at it is reported.
  `dt_routing: 'auto'` requires it, since a river too long for the runoff step has no valid step of its own until it is
  split. It multiplies the routing work by the average number of sub-reaches per river, while the subcycles of each
  sub-reach are combined into one step. A river routed in subcycles interpolates its
  upstream inflow linearly within each runoff step; see the [time options](time-options.md#how-each-river-computes).
  The state files then hold one value per sub-reach, each river's sub-reaches upstream to downstream, in the order a
  final state file is written, with each row's `riverId` repeated for every sub-reach of its river. A state file with
  one value per river is refused, and so is a run whose sub-reaches change between its forcing files of different time
  steps.
- `restrict_x_range` - `True` (the default) routes every river with an x from 0 to 0.45, lowering any larger x to
  0.45. Every coefficient of a river is
  non-negative only when its step is within `(1 - 2x) * k` of its travel time, so an x near 0.5 needs many
  sub-reaches and subcycles to find one; 0.45 holds them to a few. `False` routes x up to 0.5 as given.

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
- `dt_discharge`, or `dt_routing` in seconds - the timestep at which rivers hand their discharge downstream
- `dt_total` - total simulation duration in seconds

**`forcing: catchment`, `grid`, or `ecmwf_grib`** also requires a water input source:

- `forcing_files`
- `grid_weights_file` when `forcing` is `grid` or `ecmwf_grib`. It must not be set for
  `catchment`.

  Time keys for forced procedures (`dt_total`, `dt_discharge`) are resolved from the inputs where
  possible, and the runoff time step is always read from the forcing files; see the
  [time options](time-options.md). `start_datetime` is only read by channel routing, which has no input dates to copy.

Every network file has the columns `riverId`, `nextRiverId`, `muskingumK`, `muskingumX`, `riverIndex`, and
`upstreamCount`; see the [File Schemas reference](io-file-schema.md#network-file).

The following table lists where each remaining key applies.

| Config key                 | Description                      | Required when                                          |
|----------------------------|----------------------------------|--------------------------------------------------------|
| **core**                   |                                  |                                                        |
| `network_file`             | Network file parquet.            | always                                                 |
| **state**                  |                                  |                                                        |
| `channel_state_init_file`  | Path to initial channel state    | `forcing: channel` (else optional, default 0)          |
| `channel_state_final_file` | Path to save final channel state | optional                                               |
| **output**                 |                                  |                                                        |
| `discharge_dir`            | Directory for output  files      | _Option 1_                                             |
| `discharge_files`          | Explicit output paths            | _Option 2_                                             |
| **input data**             |                                  |                                                        |
| `forcing_files`            | Runoff read as the `forcing`     | `forcing` not `channel`                                |
| `grid_weights_file`        | Aggregates grids to catchments   | `forcing` a grid                                       |
| **time**                   |                                  |                                                        |
| `start_datetime`           | Channel routing start date       | optional, `forcing: channel` only                      |
| `dt_total`                 | Total simulation duration        | `forcing: channel` (else [time docs](time-options.md)) |
| `dt_discharge`             | Output timestep                  | optional - [time docs](time-options.md)                |
| `dt_routing`               | `'auto'`, or every river's step  | optional - [time docs](time-options.md)                |
| `conditioning`             | Make every coefficient valid     | optional, required by `dt_routing: auto`               |

## Optional configs with defaults

The `var_` options name the variables of the gridded forcing files that are read. Every other file, read or written,
uses fixed names, such as `riverId` for the river id.

| Config Key               | Description                                              | Default                                       |
|--------------------------|----------------------------------------------------------|-----------------------------------------------|
| `forcing`                | Channel routing, or the form of the forcing files        | `'channel'`                                   |
| `dt_routing`             | `'auto'` each river's own valid step, or seconds         | `'auto'`                                      |
| `conditioning`           | Sub-reaches and subcycles for rivers invalid at the step | `True`                                        |
| `restrict_x_range`       | Route x from 0 to 0.45 rather than up to 0.5             | `True`                                        |
| `log`                    | Enable or disable logging                                | `True`                                        |
| `progress_bar`           | Show tqdm progress bar                                   | `True`                                        |
| `log_level`              | Logger level, defaults to between INFO and WARNING       | `'PROGRESS'`                                  |
| `log_stream`             | `'stdout'` or a file path                                | `'stdout'`                                    |
| `log_format`             | Python logging format string                             | `'%(levelname)s - %(asctime)s - %(message)s'` |
| `var_grid_runoff`        | Runoff variable name in grid `forcing_files`             | `'ro'`                                        |
| `var_x`                  | X-dimension name in `grid` forcing files                 | `'x'`                                         |
| `var_y`                  | Y-dimension name in `grid` forcing files                 | `'y'`                                         |
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

`invalid_coefficients` controls what happens when a river has a negative Muskingum coefficient at its own
step. Every coefficient is non-negative, and the routing valid, when `2*k*x <= step <= 2*k*(1-x)`. Outside
that window routing is still stable and still conserves mass, but the river's discharge can dip below zero ahead of a rise (`c1 < 0`, a river
too long for the step) or alternate from step to step and go negative (`c3 < 0`, a river too short for it). Negative
discharge is written as routed, never clamped, so the written discharge conserves mass and shows the river that
produced it. The default `'warn'` issues a warning saying how many rivers are affected; `'raise'` refuses to route;
`'ignore'` is silent. Routing checks the coefficients with `Network.check_validity` before it starts.

Use `Network.invalid_mask(dt)` to inspect a network before routing it. `dt_routing: auto` gives each river the
`Network.subcycles(dt_runoff)` steps of its own that it needs, and conditioning routes each river too long for its
step in `Network.sub_reaches(dt)` sub-reaches inside the kernel, laid out by
`Network.sub_reach_k_and_x`, which keep the river's mean travel time and its travel time variance `k**2 * (1 - 2x)`:
equal sub-reaches with `x = 1/2 - N*(1/2 - x)` wherever that is not negative, and otherwise, since Muskingum x is
never negative, sub-reaches with x = 0 whose unequal travel times keep the variance. At `x = 0.2` two sub-reaches of
0.724 and 0.276 of k make every river too long for any step non-negative. Each sub-reach takes its share of the
river's runoff. Each river too short for a step given in seconds, which sub-reaches cannot fix, is routed in
`Network.subcycles(dt)` steps of its own. `Network.conditioning(dt)` gives both. `Network.stabilize(dt)` instead rewrites the
network in place with those sub-reaches as rows of their own, for the rivers sub-reaches can fix, and
`Network.write_stabilized(dt)` saves that network as a network table. Such a network routes with the `channel`,
`grid`, and `ecmwf_grib` forcings, whose runoff is shared among the sub-reaches of each river; `catchment` forcing
refuses it, since a catchment runoff file has one row per river.

## Example Configuration

The general template lists every key.

```json title="config.json"
--8<-- "examples/config.json"
```
