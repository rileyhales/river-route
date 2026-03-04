## Time Variables

There are 4 time variables that determine how routing is performed. All times are given in seconds. 
Once set at the start of the simulation, times must stay constant for all computations

- `dt_total` - the total simulation duration in seconds. Defaults to the full length of the input runoff data.
- `dt_discharge` - the interval between discharge values which get written to disc. Defaults to the runoff timestep. 
- `dt_runoff` - the interval between runoff volume values. Auto-detected from the input file.
- `dt_routing` - the interval between routing calculation steps. Defaults to the runoff timestep.

The following relationships must be true for all 4 time variables

```
dt_total >= dt_discharge >= dt_runoff >= dt_routing
dt_total % dt_discharge == 0
dt_discharge % dt_runoff == 0
dt_runoff % dt_routing == 0
dt_total == dt_runoff * number_of_time_steps
```

!!! warning "Left and Right aligned time labels"
    Runoff data may be given by right- or left-aligned dates and intervals. For example, hourly runoff data may say 5pm and the runoff number may
    represent the runoff that occurred from 4pm to 5pm (right-aligned) or from 5pm to 6pm (left-aligned). The routing calculations can be 
    performed with date labels of either type. You need to know this in advance to avoid an off-by-one error. 


- The total time of the runoff data must be greater than or equal to the time interval between discharge values.
- The time interval between discharge values must be greater than or equal to the time interval between runoff values.
- The time interval between runoff values must be greater than or equal to the time interval between routing
  calculations.
- Each time interval must be evenly divisible by the next smallest time interval so that the loops of calculations can
  be automatically constructed.
- The length of the routing computations must be equal to the length of the runoff data.

1. You must route each interval of runoff data at least 1 time before the next interval of runoff data. You should
   resample the runoff data to larger or smaller time intervals using your own methods before routing.
2. You cannot run the routing calculations for longer than runoff data are available. If the runoff data cover 10 days,
   the routing computations will only last 10 days. Routing longer can be accomplished by padding the runoff data with
   zero values for more time steps.
