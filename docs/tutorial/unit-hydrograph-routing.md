## Unit Hydrograph Routing with UnitMuskingum

`UnitMuskingum` convolves runoff depths with a precomputed unit hydrograph (UH) kernel before routing
through the Muskingum channel equations. Instead of placing runoff directly at the channel inlet each
timestep (as `RapidMuskingum` does), the UH spreads overland flow through time, producing a more
realistic lateral inflow hydrograph. The convolved lateral inflow is then routed between river segments
using the Muskingum method.

See the [Math Derivations](../references/math.md) page for the full equations.

## Building a Kernel

river-route provides two SCS (NRCS) dimensionless unit hydrograph builders that generate kernels from
basin properties:

- **`SCSTriangular`** — piecewise-linear triangular UH (NEH 630, Ch. 16)
- **`SCSCurvilinear`** — tabulated curvilinear UH with 33 ordinates (NEH 630, Ch. 16)

Both require three inputs:

| Parameter | Type             | Description                                      |
|-----------|------------------|--------------------------------------------------|
| `tc`      | 1D array (float) | Time of concentration per basin, in seconds      |
| `area`    | 1D array (float) | Basin (catchment) area per basin, in m²          |
| `tr`      | float            | Duration of the unit rainfall pulse, in seconds  |

```python
import numpy as np
import pandas as pd
import xarray as xr
import river_route as rr

# Load basin properties — tc and area must be in the same order as the routing params
params = pd.read_parquet('params.parquet')
weights = xr.open_dataset('grid_weights.nc').to_dataframe()

tc = params['k'].values * 5          # example: scale Muskingum k to estimate tc
area = (
    weights
    .groupby('river_id')['area_sqm']
    .sum()
    .reindex(params['river_id'])
    .values
)

# Build and save the kernel — tr must match the timestep of your runoff data
kernel = rr.uhkernels.SCSTriangular(tc=tc, area=area, tr=3600)
kernel.save('kernel.npz')
```

!!! warning "tr must match dt_runoff"
    The kernel's `tr` parameter **must equal the timestep of the runoff data** used during routing.
    Volume conservation requires `sum(kernel[:, j]) * tr == area[j]` for each basin. If `tr` does
    not match the actual runoff timestep, cumulative discharged volumes will be scaled by
    `dt_runoff / tr` — for example, a kernel built with `tr=3600` used with 3-hourly data will
    produce volumes 3× too large.

### Kernel file format

The kernel is saved as a scipy sparse npz file with shape `(n_kernel_steps, n_basins)`:

- Each **column** is the discretized unit hydrograph for one basin
- Each value is the average flow rate per unit depth of runoff (m²/s) over one `tr` interval
- Volume conservation: `sum(kernel[:, j]) * tr == area[j]`

You can also build a kernel manually and save it:

```python
import numpy as np
import scipy.sparse

n_steps = 5
n_basins = 100
kernel = np.zeros((n_steps, n_basins), dtype=np.float64)

# ... populate each column with your unit hydrograph ordinates ...

scipy.sparse.save_npz('kernel.npz', scipy.sparse.csr_matrix(kernel))
```

## Running the Router

```python
import river_route as rr

(
    rr
    .UnitMuskingum(
        params_file='params.parquet',
        uh_kernel_file='kernel.npz',
        qlateral_files=['depths.nc'],
        discharge_dir='output/',
    )
    .route()
)
```

The `qlateral_files` for `UnitMuskingum` contain **runoff depths in meters** (not volumes). The
`river_id` dimension must match the routing parameters file in both values and order. See the
[File Schemas](../references/io-file-schema.md) reference for the full specification.

Alternatively, you can route directly from gridded runoff using a weight table:

```python
(
    rr
    .UnitMuskingum(
        params_file='params.parquet',
        uh_kernel_file='kernel.npz',
        grid_runoff_files=['era5_runoff.nc'],
        grid_weights_file='grid_weights.nc',
        discharge_dir='output/',
    )
    .route()
)
```

## Warm-Starting the Convolution State

The UH convolution produces a tail that extends beyond the current simulation window. To get
continuous results across sequential input files, save and restore the convolution state:

```python
# First run — save the final UH and channel state
(
    rr
    .UnitMuskingum(
        params_file='params.parquet',
        uh_kernel_file='kernel.npz',
        qlateral_files=['depths_period1.nc'],
        discharge_dir='output/',
        channel_state_final_file='channel_state.parquet',
        uh_state_final_file='uh_state.parquet',
    )
    .route()
)

# Second run — warm-start from saved state
(
    rr
    .UnitMuskingum(
        params_file='params.parquet',
        uh_kernel_file='kernel.npz',
        qlateral_files=['depths_period2.nc'],
        discharge_dir='output/',
        channel_state_init_file='channel_state.parquet',
        uh_state_init_file='uh_state.parquet',
    )
    .route()
)
```

The UH state file is a parquet with shape `(n_basins, n_kernel_steps)` — one row per basin. It
records how much of each basin's recent runoff has not yet been discharged into the channel.

!!! note "Volume in the state buffer"
    At the end of a simulation, the UH state buffer may hold non-trivial volume that has not yet
    appeared in the discharge output. If you are checking cumulative volume conservation, account
    for the volume remaining in both the UH state and the channel state.

## Headwater Optimization

`UnitMuskingum` automatically identifies headwater segments (those with no upstream connections) and
excludes them from the Muskingum matrix solve. Headwater discharge is set directly to the UH
convolution output. This typically cuts the linear solve size roughly in half with no loss of accuracy.
