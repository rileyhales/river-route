## Unit Hydrograph Routing with UnitMuskingum

`UnitMuskingum` applies a **unit hydrograph (UH) convolution** to each timestep of runoff before it
enters the Muskingum channel routing equations. Each runoff depth is convolved against a precomputed
kernel to distribute overland flow through time before it reaches the channel, rather than being
instantly placed at the inlet. The convolved lateral inflow is then routed between river segments
using the Muskingum method.

A kernel file must be provided via `transformer_kernel_file` in the config.

## Kernel File Format

The kernel is a scipy sparse npz file with shape `(n_basins, n_time_steps)`:

- each **row** is the discretized unit hydrograph for one basin (response to one unit of runoff)
- values represent the **average flow per unit input** (m²/s per m of runoff depth) over each `dt` interval
- each row should sum to `1.0` (unit volume conservation)

```python
import numpy as np
import scipy.sparse

# Example: build a triangular kernel manually and save it
n_basins = 100
n_steps = 12  # must cover the full unit hydrograph base time
kernel = np.zeros((n_basins, n_steps))

# ... populate kernel rows with your unit hydrograph values ...

scipy.sparse.save_npz('kernel.npz', scipy.sparse.csr_matrix(kernel))
```

## Routing with a Kernel File

```yaml
params_file: '/path/to/params.parquet'
qlateral_files: '/path/to/depths.nc'
discharge_dir: '/path/to/output/'
transformer_kernel_file: 'kernel.npz'
dt_routing: 3600
```

```python
import river_route as rr

rr.UnitMuskingum('config.yaml').route()
```

## Warm-Starting the Convolution State

The convolution state (the rolling accumulation buffer) can be persisted between runs to produce
a continuous simulation from multiple sequential input files.

```yaml
# First run: save the final convolution state
params_file: 'params.parquet'
qlateral_files: 'depths_period1.nc'
discharge_dir: 'output/'
transformer_kernel_file: 'kernel.npz'
transformer_state_final_file: 'state.parquet'
dt_routing: 3600
```

```yaml
# Second run: warm-start from the saved state
params_file: 'params.parquet'
qlateral_files: 'depths_period2.nc'
discharge_dir: 'output/'
transformer_kernel_file: 'kernel.npz'
transformer_state_init_file: 'state.parquet'
dt_routing: 3600
```

The state file has the same shape as the kernel `(n_basins, n_time_steps)` and records how much of
each basin's past runoff has not yet been discharged into the channel.
