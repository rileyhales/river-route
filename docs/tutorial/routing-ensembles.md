## Routing Runoff Ensembles

Ensemble routing means running multiple runoff members over the same period using the same
routing parameters and initial conditions. The main difference from sequential single-member
routing is state handling:

- In sequential routing, each run initializes from the previous run's final state.
- In ensemble routing, all members generally start from the same initial state, then member
  final states are combined into a new initialization state by your chosen method.

There are two common ways to run ensembles in `river-route`.

Ensemble routing is supported by `RapidMuskingum` and `UnitMuskingum` (not the base `Muskingum`).

1. Run each member in a separate job (loop, multiprocessing, cluster workers). This is easiest
   to parallelize and gives full control over member-specific output paths.
2. Provide a list of input and output files and set `runoff_processing_mode=ensemble`. Members
   are routed sequentially from the same initial channel state, and the router sets the next
   `channel_state` to the mean of member final states.

## Customizing initial and final state files

Ensembles are common in forecasting, where you may need the next initialization state at a fixed
lead time (for example, +24 h) instead of the final routed timestep. In that case, do not rely on
`channel_state_final_file` alone. Instead, write a custom output function that:

1. writes routed discharge outputs, and
2. extracts the specific timestep you want for the next initialization state while data are in memory.

This avoids a second read/filter pass over output files. Your custom function needs to know either
the target datetime or the timestep offset corresponding to the next forecast cycle. For more
information, see the [advanced concepts section](advanced.md#customizing-outputs).

After writing these member-specific init files to disk, combine them into a single final state.
You could calculate an average
or median of the member's final states, assimilate gauge data with a Kalman filter, or use some
other algorithm or machine learning to aggregate the individual members' final states.

```python title="Custom Outputs for Ensemble Member Init Files"
import glob
import os

import pandas as pd
import xarray as xr

import river_route as rr


def custom_output_writer(dates, discharge_array, discharge_file, runoff_file):
    # dates: datetime array for routed discharge rows
    # discharge_array: routed flows with shape (time, river_id)
    # discharge_file: the path to the output file provided by your config file
    # runoff_file: the path to the runoff file used to produce this output, if you need it

    with xr.open_dataset(runoff_file) as runoff_ds:
        river_ids = runoff_ds['river_id'].values
    df = pd.DataFrame(discharge_array, index=pd.to_datetime(dates), columns=river_ids)

    # you probably want to include the member number in the output file name which could come from the discharge or runoff file
    member_number = os.path.basename(runoff_file)

    # option 1
    init_values = df.loc['2023-10-01 12:00:00'].T  # for if you know the exact time step to use
    init_values.to_parquet(f'member_init_from_{member_number}.parquet')  # write the next state to a file

    # option 2
    init_values = df.iloc[24, :].T  # for if you know the number of time steps after initialization to use
    init_values.to_parquet(f'member_init_from_{member_number}.parquet')  # write the next state to a file

    # continue with writing the full outputs
    ds_out = xr.Dataset(
        data_vars={'Q': (('time', 'river_id'), discharge_array)},
        coords={'time': pd.to_datetime(dates), 'river_id': river_ids}
    )
    ds_out.to_netcdf(discharge_file)
    return


m = (
    rr
    .RapidMuskingum('your_config_file.yaml')
    .set_write_discharges(custom_output_writer)  # set the custom output writer function
    .route()
)

# Now find all the member init files and average or otherwise combine them to get the next state
member_init_files = sorted(glob.glob('member_init_from_*.parquet'))
combo_init = pd.concat([pd.read_parquet(f) for f in member_init_files], axis=1).mean(axis=1)
combo_init = combo_init.to_frame(name='Q')
combo_init.to_parquet('ensemble_init_state.parquet')  # write the combined state to a file
```
