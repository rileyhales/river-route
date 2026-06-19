## Unit Hydrograph Routing

!!! warning "Planned — not yet available in v3"
    Unit hydrograph **routing** is being migrated in a later phase of the v3 line and is **not yet
    available**. Building UH kernels with `rr.uhkernels` (the first steps below) works today. The routing
    examples that follow use the v2 `UnitMuskingum` API for reference only; the v3 `rr.Router` interface for
    unit-hydrograph forcing will be defined and documented when the migration lands. See the
    [v2 → v3 migration guide](../migrating/v2-to-v3.md).

Unit hydrograph routing convolves runoff depths with a precomputed unit hydrograph (UH) kernel before routing
through the Muskingum channel equations. Instead of placing runoff directly at the channel inlet each
timestep (as `forcing: lateral` does), the UH spreads overland flow through time, producing a more
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
# v2 API shown for reference — does NOT run in v3 (UnitMuskingum was removed)
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
# v2 API shown for reference — does NOT run in v3 (UnitMuskingum was removed)
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
# v2 API shown for reference — does NOT run in v3 (UnitMuskingum was removed)
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

!!! warning "Planned — not yet available in v3"
    When UH routing lands, the UH state file will be a parquet with shape `(n_basins, n_kernel_steps)` —
    one row per basin — recording how much of each basin's recent runoff has not yet been discharged into
    the channel. The v3 `rr.Router` does not read or write `uh_state` files today.

!!! note "Volume in the state buffer"
    At the end of a simulation, the UH state buffer may hold non-trivial volume that has not yet
    appeared in the discharge output. If you are checking cumulative volume conservation, account
    for the volume remaining in both the UH state and the channel state.

## Headwater Optimization

!!! warning "Planned — not yet available in v3"
    When UH routing lands, the router will identify headwater segments (those with no upstream
    connections) and set their discharge directly to the UH convolution output, so they can be skipped
    in the routing step. This optimization is not implemented in v3 today.
