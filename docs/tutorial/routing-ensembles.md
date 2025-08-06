## Routing Runoff Ensembles

Ensemble simulations are usually 1) multiple simulations, 2) from the same starting state, 3) covering the same time period, and 4) using the same
model parameters. The runoff projections that are the members of the ensemble are described as iterations, realizations, formulations, or
perturbations. Routing an ensemble of runoffs has 1 major difference, compared to routing an individual runoff or a series of runoffs, when
determining the final state of the ensemble to use as the state to begin from for the next time step. There are 2 methods for routing an ensemble of 
runoffs in river-route.

1. In a loop or in concurrent/parallel jobs, route each member using same initial conditions and routing parameters for each job. Use the member number 
   to make output file names unique and more readily searchable and sortable. Provide all configurations as normal without special considerations. In 
   this situation, you create a new `rr.Muskingum` instance for each member which each read configs and parameter files. However, you can process each 
   member simultaneously which may be faster.
2. Provide a list of input and output files and specify `input_type=ensemble` in your configs. This will route each member in the order given by the
   input list, one at a time. Each member will be initialized from the same final state and use the same routing parameters. The final state is still 
   from the last routing step but is the average of the final step across all members, not just the first. In this situation, you create a single 
   instance of `rr.Muskingum`, read configs and parameter files only one time, but the members are processed sequentially which could be slow.

## Routing ensemble members simultaneously

Ensembl members are independent of each other and start from the same conditions so they can be computed simultaneously to possibly save time and 
more fully utilize computer capacity. There are many methods for concurrency or parallel processing in Python. This is only one example of how to 
implement it using `multiprocessing`. You will need to write some custom logic to:

1. Set the config variables and routing parameter files
2. Find all the volumes files, 1 for each member 
3. Determine unique output names for each so you don't overwrite or corrupt files
4. Submit each input/output file pair to a parallel processing framework

```python title="Parallel Routing Jobs for Ensemble Members"
from multiprocessing import Pool

import river_route as rr

connectivity_file = 'connectivity.parquet'
routing_params_file = 'routing_parameters.parquet'
volumes_files = ['volumes_member_1.nc',
                 'volumes_member_2.nc', ]
output_files = ['outflows_member_1.nc',
                'outflows_member_2.nc', ]


def route(input_file: str, output_file: str) -> None:
    rr.Muskingum(
        connectivity_file=connectivity_file,
        routing_params_file=routing_params_file,
        catchment_volumes_file=input_file,
        outflow_file=output_file,
    )


if __name__ == '__main__':
    with Pool() as pool:
        pool.starmap(route, zip(volumes_files, output_files))
```

## Customizing initial and final state files

Ensembles are often used for forecast simulations. Forecast lead times are typically several days or weeks long but new forecasts are generated one 
or more times a day. In this case, each day that your route data, you want to initialize at +24 hours from the previous forecast's initial time 
instead of the final time step. To get a final state at a specific time step instead of the final step, you should not provide a `final_state_file` 
value. Instead, write a custom output function which write the routed discharge to disc and also selects values at a specific time step to write to a 
final state file while the values are still in memory. This is the most efficient method since you avoid needing to load and filter the outputs in a 
separate process. Your custom function will need to know the datetime of the next simulation or the number of timesteps after initialization which 
corresponds to your next model run's start time. For more information about this, see the 
[advanced concepts section on custom outputs](./advanced-tutorial.md#customizing-outputs) and the 
[API documentation](../api.md#river_route.Muskingum.write_outflows).

After you have written these custom outputs to disc for each member, you can combine them into a single final state. You could calculate an average 
or median of the member's final states, assimilate gauge data with a Kalman filter, or use some other algorithm or machine learning to aggregate the 
individual members' final states.

```python title="Custom Outputs for Ensemble Member Init Files"
import glob
import os

import pandas as pd

import river_route as rr


def custom_output_writer(df, outflow_file, runoff_file):
    # df: a dataframe with datetime index and river_id numbers for columns with the routed discharge
    # outflow_file: the path to the output file provided by your config file
    # runoff_file: the path to the runoff file used to produce this output, if you need it

    # you probably want to include the member number in the output file name which could come from the outflow or runoff file
    member_number = os.path.basename(runoff_file)

    # option 1
    init_values = df.loc['2023-10-01 12:00:00'].T  # for if you know the exact time step to use
    init_values.to_parquet(f'member_init_from_{member_number}.parquet')  # write the next state to a file

    # option 2
    init_values = df.iloc[24, :].T  # for if you know the number of time steps after initialization to use
    init_values.to_parquet(f'member_init_from_{member_number}.parquet')  # write the next state to a file
    
    # continue with writing the full outputs
    df.to_netcdf(outflow_file)  # for recommendations of how to format this, refer to the API docs and/or source code
    return


m = (
    rr
    .Muskingum('your_config_file.yaml')
    .set_write_outflows(custom_output_writer)  # set the custom output writer function
    .route()
)

# Now find all the member init files and average or otherwise combine them to get the next state
member_init_files = sorted(glob.glob('member_init_from_*.parquet'))
combo_init = pd.concat([pd.read_parquet(f) for f in member_init_files], axis=1).mean(axis=1)
combo_init.columns = ['Q', ]  # rename the column to 'Q'
combo_init['R'] = 0  # add a column for runoff values at this step. You could get actual values from the runoff file sources.
combo_init.to_parquet('ensemble_init_state.parquet')  # write the combined state to a file
```
