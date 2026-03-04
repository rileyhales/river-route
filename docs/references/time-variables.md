## Time Variables

Time requirements depend on which router you use:

**`Muskingum`** (channel-only routing, no lateral inflow) requires:

- `dt_routing` - the routing computational timestep in seconds
- `dt_total` - the total simulation duration in seconds

**`RapidMuskingum` and `UnitMuskingum`** (lateral inflow routers) auto-detect time parameters
from the input runoff files. All time options are optional — override them only when you need
a different routing sub-step or output frequency:

- `dt_routing` - the time interval, in seconds, between routing calculation steps. Defaults to the runoff timestep.
  It must be constant across all rivers and for the full simulation.
- `dt_discharge` - the time interval, in seconds, between discharge values which get written to disc. Defaults to
  the runoff timestep. It must be constant between all discharge time steps.
- `dt_total` - the total simulation duration in seconds. Defaults to the full length of the input runoff data.
- `dt_runoff` - the time interval, in seconds, between runoff volume values. Auto-detected from the input file.
  It must be constant between all time steps of runoff.

The following relationships must be true for all 4 time variables

```
dt_total >= dt_discharge >= dt_runoff >= dt_routing
dt_total % dt_discharge == 0
dt_discharge % dt_runoff == 0
dt_runoff % dt_routing == 0
dt_total == dt_runoff * number_of_time_steps
```

- The total time of the runoff data must be greater than or equal to the time interval between discharge values.
- The time interval between discharge values must be greater than or equal to the time interval between runoff values.
- The time interval between runoff values must be greater than or equal to the time interval between routing
  calculations.
- Each time interval must be evenly divisible by the next smallest time interval so that the loops of calculations can
  be automatically constructed.
- The length of the routing computations must be equal to the length of the runoff data.

## Right and Left Aligned Time Series

!!! warning
    Runoff data may be given by right- or left-aligned dates and intervals. For example, hourly runoff data may say 5pm and the runoff number may
    represent the runoff that occurred from 4pm to 5pm (right-aligned) or from 5pm to 6pm (left-aligned). The routing calculations are performed and
    interpretations are most naturally computed with right-aligned dates but can be valid with left-aligned dates. Beware the alignment of your
    data, or you might inadvertently average or shift data off by 1 interval/step.

## Limitations

1. You must route each interval of runoff data at least 1 time before the next interval of runoff data. You should
   resample the runoff data to larger or smaller time intervals using your own methods before routing.
2. You cannot run the routing calculations for longer than runoff data are available. If the runoff data cover 10 days,
   the routing computations will only last 10 days. Routing longer can be accomplished by padding the runoff data with
   zero values for more time steps.
