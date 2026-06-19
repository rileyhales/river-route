## Time Options

Routing simulations have four different time steps, all given in seconds.

- `dt_routing`: The routing computation step, dt in the Muskingum equation. **The only one you need to specify**.
- `dt_runoff`: Interval between runoff inputs. If you don't provide it, it will be identified when the runoff file is opened.
- `dt_discharge`: Interval over which to average discharge to write to disc. Must be greater than or equal to `dt_runoff`.
- `dt_total`: Total simulation duration.

The most important time step is `dt_routing`. All other time steps are derived from this and the runoff inputs.

The following rules apply:

1. You must route each runoff increment at least 1 time so `dt_routing` must be less than or equal to `dt_runoff`.
2. `dt_routing` must be an integer divisor of `dt_runoff` because runoff distributions won't be resampled.
3. `dt_discharge` must be an integer multiple of `dt_runoff` because discharge outputs are averaged over runoff intervals.
4. By default `dt_total` is `dt_runoff` multiplied by the runoff record length. If you set it explicitly, it must be an integer multiple of `dt_runoff` (and of `dt_discharge`). Recession routing is not available.

## Router Defaults

- Channel-only routing (`forcing: channel`): requires `dt_total` and `dt_routing`; `dt_discharge` defaults to `dt_routing`.
- Routing with forcing (e.g. `forcing: lateral`):
    - `dt_runoff` defaults to the runoff-file timestep.
    - `dt_discharge` defaults to `dt_runoff`.
    - `dt_total` defaults to `dt_runoff * number_of_timesteps`.
    - `dt_routing` defaults to `dt_runoff`.

## Required Relationships

For routing with forcing (e.g. `forcing: lateral`):

```
dt_total >= dt_discharge >= dt_runoff >= dt_routing
dt_total % dt_discharge == 0
dt_discharge % dt_runoff == 0
dt_runoff % dt_routing == 0
```

If you do not supply `dt_total`, it defaults to `dt_runoff * number_of_timesteps`.

For channel-only routing (`forcing: channel`):

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
