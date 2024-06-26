## Parameter Descriptions

Only 1 time option is a required input in the configuration file:

- `dt_routing` - the time interval, in seconds, between routing calculation steps. It must be constant across all rivers
  and for the full simulation.

By default, river-route reports routed flows as the average flow that occured over the same time frequency that the
runoff volumes are provided in. If you provide 3-hourly runoff volumes for each catchment. Your output will be the
average flow in each river segment over the same 3-hour interval. You can change this by specifying:

- `dt_outflow` - the time interval, in seconds, between outflow values which get written to disc. It must be constant
  between all time steps of outflow.

By way of information, 2 other time parameters are calculated and verified internally by river-route.

- `dt_runoff` - the time interval, in seconds, between runoff volume values. It must be constant between all time steps
  of runoff.
- `dt_total` - the total time, in seconds, of the runoff data. It is equal to the number of time steps multiplied
  by `dt_runoff`.

The following relationships must be true for all 4 time variables

```
dt_total >= dt_outflow >= dt_runoff >= dt_routing
```

- The total time of the runoff data must be greater than or equal to the time interval between outflow values.
- The time interval between outflow values must be greater than or equal to the time interval between runoff values.
- The time interval between runoff values must be greater than or equal to the time interval between routing
  calculations.

```
dt_total % dt_outflow == 0
dt_outflow % dt_runoff == 0
dt_runoff % dt_routing == 0
```

- Each time interval must be evenly divisible by the next smallest time interval so that the loops of calculations can
  be automatically constructed.

```
dt_total === dt_runoff * number_of_time_steps
```

- The length of the routing computations must be equal to the length of the runoff data.

## Limitations

1. You must route each interval of runoff data at least 1 time before the next interval of runoff data. You should
   resample the runoff data to larger or smaller time intervals using your own methods before routing.
2. You cannot run the routing calculations for longer than runoff data are available. If the runoff data cover 10 days,
   the routing computations will only last 10 days. Routing longer can be accomplished by padding the runoff data with
   zero values for more time steps.
