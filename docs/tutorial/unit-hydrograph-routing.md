## Unit Hydrograph Routing with UnitMuskingum

`UnitMuskingum` applies a **unit hydrograph (UH) transformation** to each timestep of runoff before it
enters the Muskingum channel routing equations. This adds catchment-scale travel time and attenuation so
that runoff is spread over time before reaching the channel, rather than being teleported instantly to the
inlet. Unlike `TeleportMuskingum`, it takes per-catchment runoff **depths** as input (not volumes) — the
transformer handles depth × area and the temporal convolution internally.

The router delegates all transformation logic to a pluggable **transformer** object. Two options are available:
use the built-in SCS triangular unit hydrograph, or inject your own transformer class.

## Additional Routing Parameters

`UnitMuskingum` requires two columns beyond the base routing parameters:

| Column     | Data Type | Description                                          |
|------------|-----------|------------------------------------------------------|
| `tc`       | float     | Time of concentration (seconds), must be ≥ 0        |
| `area_sqm` | float     | Catchment area in square metres, must be > 0         |

## Option 1: Built-in SCS Triangular Unit Hydrograph

The SCS triangular UH derives a per-basin kernel from `tc` and `area_sqm` using the standard SCS equations:

- Lag time: `tl = 0.6 * tc`
- Time to peak: `tp = tl + dt/2`
- Time to base: `tb = 2.67 * tp`
- Peak flow: `qp = 2 * area / tb`

The kernel is discretized as the average flow over each `dt` interval using the antiderivative of the
triangular hydrograph, so volume is conserved exactly.

Specify `uh_type: 'scs'` in your config file:

```yaml
routing_params_file: '/path/to/params.parquet'
lateral_depth_files: '/path/to/depths.nc'
discharge_files: '/path/to/discharge.nc'
uh_type: 'scs'
```

```python
import river_route as rr

(
    rr
    .UnitMuskingum('config.yaml')
    .route()
)
```

## Caching the Kernel

For large networks with fixed `tc` and `area_sqm`, computing the kernel matrix every run is wasteful.
Save it once and reload with `uh_kernel_file`:

```python
import river_route as rr
from river_route.transformers import SCSUnitHydrograph

# Compute and save the kernel once
transformer = SCSUnitHydrograph(dt=3600, tc=tc_array, area=area_array)
transformer.save_kernel('kernel.parquet')

# Subsequent runs load from disk — no uh_type needed
(
    rr
    .UnitMuskingum(
        routing_params_file='params.parquet',
        lateral_depth_files='volumes.nc',
        discharge_files='discharge.nc',
        uh_kernel_file='kernel.parquet',
    )
    .route()
)
```

Optionally persist and reload the transformer **state** (the rolling convolution buffer) across runs:

```python
# After routing, save the transformer state
router = rr.UnitMuskingum('config.yaml').route()
router._runoff_transformer.save_kernel('kernel.parquet')  # or skip if already saved

# Next run: warm-start both kernel and state
(
    rr
    .UnitMuskingum(
        routing_params_file='params.parquet',
        lateral_depth_files='volumes_next.nc',
        discharge_files='discharge_next.nc',
        uh_kernel_file='kernel.parquet',
        transformer_state_file='state.parquet',
    )
    .route()
)
```

## Option 2: Custom Transformer

You can supply any transformer that subclasses `AbstractBaseTransformer` and implements `_build_kernel()`.
The base class handles all state management and calls `transform()` at each routing timestep — you only need
to provide the kernel.

`_build_kernel()` must return a 2D array of shape `(n_time_steps, n_basins)` where:

- each column is the discretized unit hydrograph for one basin (response to one unit of runoff)
- values represent the **average flow per unit input** over each `dt` interval
- columns should sum to `1.0` (unit volume conservation)

```python
import numpy as np
import river_route as rr
from river_route.transformers import AbstractBaseTransformer


class LinearReservoirTransformer(AbstractBaseTransformer):
    """Simple linear reservoir: u(t) = (1/K) * exp(-t/K)"""

    def __init__(self, k_values: np.ndarray, dt: float) -> None:
        self.k_values = np.asarray(k_values, dtype=np.float64)  # shape (n_basins,)
        super().__init__(dt=dt)  # validates dt, calls _build_kernel(), resets state

    def _build_kernel(self) -> np.ndarray:
        max_steps = int(np.ceil(5 * self.k_values.max() / self.dt))
        t_edges = np.arange(max_steps + 1, dtype=np.float64)[:, np.newaxis] * self.dt
        # cumulative integral of (1/K)*exp(-t/K) = 1 - exp(-t/K)
        u_cum = 1.0 - np.exp(-t_edges / self.k_values)  # (max_steps+1, n_basins)
        return np.diff(u_cum, axis=0) / self.dt           # (max_steps, n_basins)


transformer = LinearReservoirTransformer(k_values=k_array, dt=3600)

(
    rr
    .UnitMuskingum('config.yaml')
    .set_transformer(transformer)
    .route()
)
```

The only contract `UnitMuskingum` requires of the transformer is the `transform(runoff_vector)` method,
which is fully implemented by `AbstractBaseTransformer`. Your subclass never needs to override `transform()`.

## Saving and Loading a Custom Kernel

Custom transformer kernels can be saved and reloaded the same way as built-in ones:

```python
# Save
transformer.save_kernel('my_kernel.parquet')

# Reload via Transformer (no uh_type needed; skips _build_kernel)
from river_route.transformers import Transformer

t = (
    Transformer
    .from_kernel(dt=3600, kernel='my_kernel.parquet')
    .set_state('my_state.parquet')   # optional warm-start
)

(
    rr
    .UnitMuskingum('config.yaml')
    .set_transformer(t)
    .route()
)
```
