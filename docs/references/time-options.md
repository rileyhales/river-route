## Time Options

Routing simulations have four different time steps, all given in seconds.

- `dt_routing`: The routing computation step. This is the dt of the Muskingum equation and the
  only one you need to specify.
- `dt_runoff`: Interval between runoff inputs. If you don't provide it, it will be identified
  when the runoff file is opened.
- `dt_discharge`: interval over which to average discharge to write to disc. Must be greater
  than or equal to `dt_runoff`.
- `dt_total`: total simulation duration.

The most important time step is `dt_routing`. All other time steps are derived from this and the runoff inputs.

The following rules apply:

1. You must route each runoff increment at least 1 time so `dt_routing` must be less than or equal to `dt_runoff`.
2. `dt_routing` must be an integer divisor of `dt_runoff` because runoff distributions won't be resampled.
3. `dt_discharge` must be an integer multiple of `dt_runoff` because discharge outputs are
   averaged over runoff intervals.
4. `dt_total` is equal to `dt_runoff` multiplied by the number of runoff time steps. Recession routing is not available.

## Router Defaults

- `Muskingum`: requires `dt_total` and `dt_routing`; `dt_discharge` defaults to `dt_routing`.
- `RapidMuskingum` and `UnitMuskingum`:
    - `dt_runoff` defaults to the runoff-file timestep.
    - `dt_discharge` defaults to `dt_runoff`.
    - `dt_total` defaults to `dt_runoff * number_of_timesteps`.
    - `dt_routing` defaults to `dt_runoff`.

## Required Relationships

For transform routers (`RapidMuskingum`, `UnitMuskingum`):

```
dt_total >= dt_discharge >= dt_runoff >= dt_routing
dt_total % dt_discharge == 0
dt_discharge % dt_runoff == 0
dt_runoff % dt_routing == 0
dt_total == dt_runoff * number_of_timesteps
```

For channel-only `Muskingum`:

```
dt_total >= dt_discharge >= dt_routing
dt_total % dt_discharge == 0
dt_discharge % dt_routing == 0
```

## Practical Notes

1. `dt_routing` should be chosen for numerical stability and performance, not just convenience.
2. If you need a different runoff timestep, resample runoff before routing.
3. To route longer than your runoff record, pad runoff with zeros upstream of routing.

!!! warning "Left- vs right-aligned timestamps"
    Runoff timestamps can represent interval starts (left-aligned) or interval ends (right-aligned).
    Example: an hourly value labeled `17:00` may represent either `17:00-18:00` or `16:00-17:00`.
    Keep this convention consistent to avoid off-by-one timing errors.
