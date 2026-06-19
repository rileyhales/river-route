## Migrating from v2 to v3

!!! note
    v3 changes **how you select a routing procedure**, not how the routing is computed. The Muskingum math,
    the input/output file schemas, the time options, and the runoff tools are all unchanged. If your v2 results
    are correct, the same configuration will produce the same results in v3.

---

## The change in one sentence

The separate router classes (`Muskingum`, `RapidMuskingum`, `UnitMuskingum`) are replaced by a **single
`Router`** whose behavior you describe with configuration options.

```python
import river_route as rr

# v2 — the class you pick determines the routing procedure
rr.Muskingum("config.yaml").route()
rr.RapidMuskingum("config.yaml").route()

# v3 — one class; the config describes the routing procedure
rr.Router("config.yaml").route()
```

---

## Why this changed

Routing options in `river-route` are **independent of one another**. Whether you route with constant or
nonlinear parameters, whether you add lateral runoff, and how you treat the river network for numerical
stability are separate choices that can be combined freely. Expressing every combination as its own class
means the number of classes multiplies as options are added — and most combinations would never be written.

v3 separates the two kinds of choices:

- **Choices that change the computation** — the coefficient model, the inflow terms, and the network
  representation — are now configuration *options*, resolved to the correct numerical kernel internally.
- **Choices that only change orchestration** — where the inflow comes from, sequential vs. ensemble
  processing, and how output is written — were already plain options and stay that way.

The result is one router that scales to new options by adding a config value, not a class.

---

## The new mental model

| v2                                                   | v3                                                              |
|------------------------------------------------------|-----------------------------------------------------------------|
| Choose the **class** that matches your procedure.    | Use `Router` and **describe** the procedure with config keys.   |
| Procedure is fixed by the class you imported.        | Procedure is resolved from the config when you call `.route()`. |
| New procedures require new classes.                  | New procedures are new config values.                           |

Three new configuration keys describe a routing procedure:

| Key       | Values                                       | Default  | Chooses                                                          |
|-----------|----------------------------------------------|----------|------------------------------------------------------------------|
| `coeff`   | `static`, `dynamic`                          | `static` | Constant Muskingum K, or nonlinear K recomputed from discharge.  |
| `forcing` | `channel`, `lateral`, `external`             | `channel`| The inflow term; `channel` is routing only. One value at a time.  |
| `network` | `standard`, `expanded`                       | `standard` | One reach per river, or a stability-expanded reach network.    |

!!! tip
    `forcing` takes a single value: `channel` (no inflow), `lateral`, or `external`. Use `forcing: channel`
    for channel-only routing. Specifying more than one forcing term is planned for a later release.

---

## Mapping your v2 router to v3

| v2 router        | v3 equivalent                          | Config keys                          |
|------------------|----------------------------------------|--------------------------------------|
| `Muskingum`      | `Router(cfg)` or `Router(cfg, forcing='channel')` | `coeff: static`, `forcing: channel`     |
| `RapidMuskingum` | `Router(cfg, forcing='lateral')`       | `coeff: static`, `forcing: lateral` |
| `UnitMuskingum`  | *Deferred — see [Unit hydrograph routing](#unit-hydrograph-routing) below* | *unchanged for now*                  |

You can set the selectors two equivalent ways. Pick whichever fits your workflow.

**1. In the config file** — add the selector keys alongside your existing config:

```yaml title="config.yaml"
params_file: /path/to/params.parquet
discharge_dir: /path/to/output/
qlateral_files:
  - /path/to/runoff.nc
coeff: static
forcing: lateral
network: standard
```

```python
import river_route as rr

rr.Router("config.yaml").route()
```

**2. As keyword arguments** — override or supply the selectors in code:

```python
import river_route as rr

(
    rr
    .Router(
        "config.yaml",
        coeff="static",
        forcing="lateral",
    )
    .route()
)
```

Keyword arguments override values from the config file, so the same config can drive different procedures.

---

## Required parameters by selection

The columns required in your `params_file` depend on the `coeff` you select. The base requirements
(`river_id`, `downstream_river_id`) and the [file schemas](../references/io-file-schema.md) are unchanged
from v2.

| Selection            | Required `params_file` columns | Inflow source                                                  |
|----------------------|--------------------------------|----------------------------------------------------------------|
| `coeff: static`      | `k`, `x`                       | —                                                              |
| `coeff: dynamic`     | `k`, `x`, `alpha`, `beta`      | —                                                              |
| `forcing: lateral` | (as per `coeff`)               | `qlateral_files`, or `grid_runoff_files` + `grid_weights_file` |

`Router` validates the config keys for the procedure you described up front and reports what is missing.
The `params_file` columns are checked when the file is read; if you select a combination that has no kernel
yet (see below), `dispatch` raises `NotImplementedError` and lists the implemented combinations.

---

## Newly available in v3

Because the procedure is now described rather than fixed by a class, v3 makes combinations available that did
not have a class in v2. These are opt-in; existing configurations keep their v2 behavior.

Available now:

- **Nonlinear coefficients** — `coeff: dynamic` recomputes K from discharge each step
  (K = α·Qᵝ) using `alpha`/`beta` columns. Currently wired only with `forcing: lateral`; pairing
  `coeff: dynamic` with the default `forcing: channel` raises `NotImplementedError`.

Planned for a later v3 release (selecting these today raises `NotImplementedError`):

- **Stability-expanded networks** — `network: expanded` to automatically subdivide or sub-step reaches that
  are otherwise numerically unstable for your time step.
- **External discharge forcing** — `forcing: external` to add an external discharge series to the channel
  and propagate it downstream (for example, inserting gauge or reservoir-release discharge).
- **Multiple forcings** — combining inflow terms (e.g. lateral runoff plus an inserted discharge series)
  in a single run.
- **Unit hydrograph routing** — see [below](#unit-hydrograph-routing).

---

## Backwards compatibility

!!! warning "v3 is a hard cut"
    The v2 router classes are **removed** in v3 — there are no deprecation shims. Importing or calling
    `rr.Muskingum`, `rr.RapidMuskingum`, or `rr.UnitMuskingum` will raise `AttributeError`. Update every call
    site to `rr.Router(...)` using the [mapping above](#mapping-your-v2-router-to-v3).

```python
import river_route as rr

# v2 — no longer exists
rr.RapidMuskingum("config.yaml").route()        # AttributeError

# v3
rr.Router("config.yaml", forcing="lateral").route()
```

---

## Command line interface

The per-router subcommands are replaced by a single `route` command that reads the selectors from your config
file.

```bash
# v2
rr Muskingum config.yaml
rr RapidMuskingum config.yaml

# v3
rr route config.yaml
```

---

## What is *not* changing

These are identical to v2 — no action needed:

- **Routing math** — the Muskingum equations and forward-substitution solver
  (see [Math Derivations](../references/math.md)).
- **File schemas** — `params_file`, grid weights, state files, and discharge output formats
  (see [File Schemas](../references/io-file-schema.md)).
- **Time options** — `dt_routing`, `dt_runoff`, `dt_discharge`, `dt_total`
  (see [Time Variables](../references/time-options.md)).
- **Runoff tools** — grid-weighting and runoff-to-qlateral conversion in `river_route.runoff`.
- **Custom output writers** — `set_write_discharges(...)` works on `Router` exactly as before and still
  returns the instance for chaining.
- **Processing modes** — `runoff_processing_mode: sequential | ensemble` is unchanged.

---

## Unit hydrograph routing

Unit hydrograph routing (`UnitMuskingum` in v2) is **being migrated in a later phase of the v3 line** and is
intentionally out of scope for this first transition. Its place in the `Router` model — and how its
convolution state is handled alongside channel state — will be defined and documented separately.

!!! warning
    If your workflow depends on unit hydrograph routing, hold on that part of your pipeline until the unit
    hydrograph migration guidance is published. The `Muskingum` and `RapidMuskingum` procedures described
    above are unaffected.
