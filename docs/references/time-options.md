## Time Options

Routing simulations have four different time steps, all given in seconds.

- `dt_routing`: The routing computation step, dt in the Muskingum equation. Channel routing requires it, and routing
  with forcing defaults it from `dt_runoff`.
- `dt_runoff`: Interval between runoff inputs. It is not a config: it is read from the time steps of each runoff file.
- `dt_discharge`: Interval over which to average discharge to write to disk. Must be greater than or equal to `dt_runoff`.
- `dt_total`: Total simulation duration.

The most important time step is `dt_routing`. `dt_discharge` and `dt_total` are derived from the runoff inputs unless
the configs give them.

The following rules apply:

1. You must route each runoff increment at least 1 time so `dt_routing` must be less than or equal to `dt_runoff`.
2. `dt_routing` must be an integer divisor of `dt_runoff` because runoff distributions won't be resampled.
3. `dt_discharge` must be an integer multiple of `dt_runoff` because discharge outputs are averaged over runoff intervals.
4. By default `dt_total` is `dt_runoff` multiplied by the runoff record length. If you set it explicitly, it must be an integer 
   multiple of `dt_runoff` (and of `dt_discharge`). Recession routing is not available.
5. Every step of a runoff file must be as long as its first, which is its `dt_runoff`. Runoff is never resampled, so a
   file with irregular time steps is refused, and so is a file with a single time step, which has no step to read,
   or with dates that do not increase.

## Router Defaults

- Channel-only routing (`forcing: channel`): requires `dt_total` and `dt_routing`; `dt_discharge` defaults to `dt_routing`.
- Routing with forcing (e.g. `forcing: catchment`):
    - `dt_runoff` is the time step of each runoff file.
    - `dt_discharge` defaults to `dt_runoff`.
    - `dt_total` defaults to `dt_runoff * number_of_timesteps`.
    - `dt_routing` defaults to `dt_runoff`, or with `network_type: stabilized` to the largest divisor of `dt_runoff`
      at which every river can be made stable.

## Required Relationships

For routing with forcing (e.g. `forcing: catchment`):

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
    Runoff timestamps can be given as interval starts (left-aligned) or interval ends (right-aligned).
    Example: an hourly value labeled `17:00` may represent either the hour `17:00-18:00` or `16:00-17:00`.
    Keep this convention consistent to avoid off-by-one timing errors.
