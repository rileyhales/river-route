## Unit Hydrograph Routing with UnitMuskingum

`UnitMuskingum` applies a **unit hydrograph (UH) transformation** to each timestep of runoff before it
enters the Muskingum channel routing equations. This adds catchment-scale travel time and attenuation so
that runoff is spread over time before reaching the channel, rather than being instantly placed at the
inlet. The router delegates all transformation logic to a pluggable **transformer** object. A transformer must be
provided before calling `route()`, via one of two paths:

- **`transformer_kernel_file` in config** — load a pre-computed kernel parquet at runtime
- **`set_transformer()` in Python** — inject a built or custom transformer directly

## Option A: Pre-computed Kernel File

Build a kernel once in Python, save it to disk, then point the config at the file. This is the recommended
approach for production runs where the network and timestep are fixed.

```python
import numpy as np
import river_route as rr
from river_route.transformers import SCSUnitHydrograph

# Read tc and area from your routing params (or any other source)
import pandas as pd

df = pd.read_parquet('params.parquet')

# Build and save the kernel once
transformer = SCSUnitHydrograph(dt=3600, tc=df['tc'].values, area=df['area_sqm'].values)
transformer.save_kernel('kernel.parquet')
```

Then route using only config:

```yaml
routing_params_file: '/path/to/params.parquet'
lateral_depth_files: '/path/to/depths.nc'
discharge_files: '/path/to/discharge.nc'
transformer_kernel_file: 'kernel.parquet'
```

```python
import river_route as rr

rr.UnitMuskingum('config.yaml').route()
```

Optionally persist and reload the transformer **state** (the rolling convolution buffer) to warm-start
a subsequent run:

```python
# After routing, save the transformer state
router = rr.UnitMuskingum('config.yaml').route()
router._Transformer.save_state('state.parquet')

# Next run: warm-start both kernel and state
(
    rr
    .UnitMuskingum(
        routing_params_file='params.parquet',
        lateral_depth_files='volumes_next.nc',
        discharge_files='discharge_next.nc',
        transformer_kernel_file='kernel.parquet',
        transformer_state_file='state.parquet',
    )
    .route()
)
```

## Option B: Dependency Injection

Build a transformer in the same script and inject it with `set_transformer()`. Use this when building
the kernel and routing in a single script, or when using a custom transformer class.

```python
import river_route as rr
from river_route.transformers import SCSUnitHydrograph

transformer = SCSUnitHydrograph(dt=3600, tc=tc_array, area=area_array)

(
    rr
    .UnitMuskingum(
        routing_params_file='params.parquet',
        lateral_depth_files='depths.nc',
        discharge_files='discharge.nc',
    )
    .set_transformer(transformer)
    .route()
)
```

## Custom Transformers

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
        return np.diff(u_cum, axis=0) / self.dt  # (max_steps, n_basins)


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

# Reload via Transformer (skips _build_kernel)
from river_route.transformers import Transformer

t = (
    Transformer
    .from_kernel(dt=3600, kernel='my_kernel.parquet')
    .set_state('my_state.parquet')  # optional warm-start
)

(
    rr
    .UnitMuskingum('config.yaml')
    .set_transformer(t)
    .route()
)
```
