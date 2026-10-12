## Time Options

Routing simulations have four different time steps, all given in seconds.

- `dt_runoff`: Interval between runoff inputs. It is not a config: it is read from the time steps of each forcing file.
  Every river hands its discharge to the river downstream at the end of each runoff step, and records it there.
- `dt_routing`: `'auto'`, the default, or the step in seconds every river computes in, `dt` in the Muskingum
  equation, a divisor of `dt_runoff`. With `'auto'` each river computes in its own valid step, the fewest steps of
  each runoff step that keep its coefficients valid, which `use_subcycles` gives it, and `use_subreaches` splits
  the rivers too long for the runoff step into sub-reaches. `'auto'` requires both.
- `dt_discharge`: Channel routing only. Interval between the discharges written to disk, and, with no runoff, the step
  at which every river hands its discharge downstream. It defaults to `dt_routing` when that is given in seconds.
- `dt_total`: Channel routing only. Total simulation duration, an integer multiple of `dt_discharge`.

Routing with forcing routes every step of each forcing file and writes the discharge at the end of each, labeled with
the time the step begins, so `dt_discharge` and `dt_total` are refused with any `forcing` but `channel`.

## How Each River Computes

Rivers are routed one at a time, upstream to downstream, each river's whole series before the next, so every river
computes in a plan of its own. Each runoff step of a river is taken in subcycles, steps of its own, with its upstream
inflow interpolated linearly across them and its runoff spread evenly over them. Each river's plan is settled in two
jobs: `dt_routing` sets its step, and `use_subreaches` and `use_subcycles` make every coefficient non-negative at it.

- With `dt_routing: auto`, the default, each river's step is its own: the runoff step over the fewest subcycles that
  keep its coefficients non-negative, which `use_subcycles` gives it. `use_subreaches` then splits a river too long
  for the runoff step into the fewest sub-reaches in series that keep its coefficients non-negative. `'auto'`
  requires both, so that every river it routes is valid where it can be.
- With `dt_routing` given in seconds, every river's step is exactly `dt_routing`, `dt_runoff / dt_routing` subcycles
  of each runoff step. With `use_subreaches`, a river too long for it is split into sub-reaches, and with
  `use_subcycles`, a river too short for it takes the fewest subcycles of `dt_routing` that keep its coefficients
  non-negative. Without them, no river is split or given a shorter step, so a river too long or too short for it
  keeps a negative coefficient and is reported, as `invalid_coefficients` says. Setting `dt_routing` to `dt_runoff`
  with both routes as `'auto'` does, and without them routes each river whole in one step of each runoff step.

The subcycles of each sub-reach are combined into one step of coefficients that give the
discharge of taking them one at a time, so they cost one step, and the sub-reaches of a river hand their discharge to
each other at the end of each runoff step, as rivers do.

The following rules apply:

1. `dt_routing`, when given in seconds, must be less than or equal to `dt_runoff` and an integer divisor of it,
   because runoff is never resampled.
2. Every step of each forcing file is routed, and the discharge at the end of each is written. Recession routing is not
   available.
3. Every step of a forcing file must be as long as its first, which is its `dt_runoff`. Runoff is never resampled, so a
   file with irregular time steps is refused, and so is a file with a single time step, which has no step to read,
   or with dates that do not increase.

## Router Defaults

- Channel-only routing (`forcing: channel`): requires `dt_total` and `dt_discharge`, or `dt_routing` in seconds. With
  no runoff, rivers hand their discharge downstream every `dt_discharge`, which defaults to `dt_routing` when it is
  given in seconds.
- Routing with forcing (e.g. `forcing: catchment`):
    - `dt_runoff` is the time step of each forcing file.
    - Every step of each forcing file is routed and written; `dt_discharge` and `dt_total` are refused.
    - `dt_routing` is `'auto'`, and `use_subreaches` and `use_subcycles` are `True`, so each river takes the fewest
      steps of each runoff step and the sub-reaches that keep it valid.

## Required Relationships

For routing with forcing (e.g. `forcing: catchment`):

```
dt_runoff >= dt_routing (when given)
dt_runoff % dt_routing == 0 (when given)
```

For channel-only routing (`forcing: channel`):

```
dt_total >= dt_discharge >= dt_routing (when given)
dt_total % dt_discharge == 0
dt_discharge % dt_routing == 0 (when given)
```

## Practical Notes

1. Leave `dt_routing` at `'auto'` unless every river must take the same step. A step in seconds keeps every river
   valid only with `use_subreaches`, which splits the rivers too long for it, and `use_subcycles`, which subcycles the
   rivers too short for it.
   A step shorter than a river's travel time does not route it more accurately:
   it does not change how much the river delays or spreads a flood wave, and it puts more rivers below the lower bound
   `2*k*x`.
2. If you need a different runoff timestep, resample runoff before routing.
3. To route longer than your runoff record, pad runoff with zeros upstream of routing.
4. To route part of a forcing file, trim the file before routing. To write fewer time steps, save only the ones you
   need with a [custom writer](../tutorial/advanced.md#customizing-outputs).

!!! warning "Left- vs right-aligned timestamps"
    Runoff timestamps can be given as interval starts (left-aligned) or interval ends (right-aligned).
    Example: an hourly value labeled `17:00` may represent either the hour `17:00-18:00` or `16:00-17:00`.
    Keep this convention consistent to avoid off-by-one timing errors.
