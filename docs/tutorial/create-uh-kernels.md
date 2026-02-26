## Generating Unit Hydrograph Kernels

`UnitMuskingum` requires a precomputed convolution kernel provided as a parquet file. This tutorial
explains the kernel format and shows how to build one using the SCS unit hydrograph methods and
custom approaches.

## Kernel File Format

The kernel file is a parquet with shape `(n_basins, n_time_steps)`:

- each **row** is the discretized unit hydrograph for one basin
- each value is the **average flow per unit depth of runoff** (m²/s) over one routing timestep `dt`
- each row must sum to `1.0` — i.e., one unit of runoff depth produces one unit of integrated flow:
  `sum(row) * dt == area` ... or equivalently `sum(row_normalized) == 1.0` when normalized by area
- column labels are integer time step indices (0, 1, 2, …)
- rows are ordered to match the `river_id` column in the routing params file

```python
import numpy as np
import pandas as pd

# kernel shape: (n_basins, n_time_steps)
# Each value: average flow (m²/s) per meter of runoff depth over one dt interval
pd.DataFrame(kernel).to_parquet('kernel.parquet')
```

## SCS Triangular Unit Hydrograph

The SCS triangular UH approximates the catchment response as a triangle in time. The required
catchment attributes are time-of-concentration `tc` (seconds) and contributing area `area` (m²).

```python
import numpy as np
import pandas as pd


def scs_triangular_kernel(dt: float, tc: np.ndarray, area: np.ndarray) -> np.ndarray:
    """
    Build an SCS triangular unit hydrograph kernel.

    Parameters
    ----------
    dt   : routing timestep in seconds
    tc   : 1D array of time-of-concentration values in seconds, one per basin
    area : 1D array of basin areas in m², one per basin

    Returns
    -------
    kernel : 2D array of shape (n_basins, n_time_steps)
             Each row sums to 1.0 (unit volume conservation when multiplied by area/dt).
    """
    tc = np.asarray(tc, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)

    tl = 0.6 * tc                  # lag time
    tp = tl + dt / 2.0             # time to peak
    tb = 2.67 * tp                 # base time
    qp = 2.0 * area / tb           # peak flow (m²/s)

    n_steps = np.ceil(tb / dt).astype(int)
    max_steps = int(n_steps.max())

    # Time edges for each interval: shape (max_steps+1, n_basins)
    i_row = np.arange(max_steps + 1, dtype=np.float64)[:, np.newaxis]
    t_edges = np.minimum(i_row * dt, tb)

    # Cumulative integral of the triangular UH
    rising  = qp * t_edges ** 2 / (2.0 * tp)
    falling = (qp * tp / 2.0
               + qp / (tb - tp) * (tb * (t_edges - tp) - (t_edges ** 2 - tp ** 2) / 2.0))
    u_cum = np.where(t_edges <= tp, rising, falling)

    # Average flow over each interval: shape (max_steps, n_basins)
    kernel = np.diff(u_cum, axis=0) / dt   # m²/s per m depth

    return kernel.T  # transpose to (n_basins, n_time_steps)


# Example usage
import pandas as pd

df = pd.read_parquet('params.parquet')
dt = 3600.0  # seconds

kernel = scs_triangular_kernel(dt=dt, tc=df['tc'].values, area=df['area_sqm'].values)
pd.DataFrame(kernel).to_parquet('kernel.parquet')
```

## SCS Curvilinear Unit Hydrograph

The SCS curvilinear UH uses the NRCS dimensionless unit hydrograph table (NEH Part 630, Ch. 16),
which gives a more realistic S-shaped rise and gradual recession. The base time extends to `5 * tp`
versus `2.67 * tp` for the triangular form.

```python
import numpy as np
import pandas as pd

# NRCS dimensionless UH table (NEH Part 630, Ch. 16, Table 16-1)
_DIM_T = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
    2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
    4.0, 4.5, 5.0,
])
_DIM_Q = np.array([
    0.000, 0.015, 0.075, 0.160, 0.280, 0.430, 0.600, 0.770, 0.890, 0.970,
    1.000, 0.980, 0.920, 0.840, 0.750, 0.660, 0.560, 0.460, 0.390, 0.330,
    0.280, 0.207, 0.147, 0.107, 0.077, 0.055, 0.040, 0.029, 0.021, 0.015,
    0.011, 0.005, 0.000,
])
# Precomputed cumulative integral of the dimensionless UH at each breakpoint
_DIM_I_CUM = np.concatenate([
    [0.0],
    np.cumsum((_DIM_Q[:-1] + _DIM_Q[1:]) / 2.0 * np.diff(_DIM_T)),
])


def scs_curvilinear_kernel(dt: float, tc: np.ndarray, area: np.ndarray) -> np.ndarray:
    """
    Build an SCS curvilinear unit hydrograph kernel using the NRCS dimensionless UH table.

    Parameters
    ----------
    dt   : routing timestep in seconds
    tc   : 1D array of time-of-concentration values in seconds, one per basin
    area : 1D array of basin areas in m², one per basin

    Returns
    -------
    kernel : 2D array of shape (n_basins, n_time_steps)
    """
    tc = np.asarray(tc, dtype=np.float64)
    area = np.asarray(area, dtype=np.float64)

    tl = 0.6 * tc
    tp = tl + dt / 2.0
    tb = _DIM_T[-1] * tp   # base time = 5 * tp

    I_dim = _DIM_I_CUM[-1]         # total dimensionless integral ≈ 1.333
    qp = area / (I_dim * tp)       # peak flow (m²/s)

    n_steps = np.ceil(tb / dt).astype(int)
    max_steps = int(n_steps.max())

    i_row = np.arange(max_steps + 1, dtype=np.float64)[:, np.newaxis]
    t_edges = np.minimum(i_row * dt, tb)             # (max_steps+1, n_basins)

    # Evaluate cumulative dimensionless integral at each edge via interpolation
    t_over_tp = t_edges / tp
    i_cum_dim = np.interp(
        t_over_tp.ravel(), _DIM_T, _DIM_I_CUM, right=I_dim
    ).reshape(t_over_tp.shape)

    u_cum = i_cum_dim * (qp * tp)                    # (max_steps+1, n_basins)
    kernel = np.diff(u_cum, axis=0) / dt             # (max_steps, n_basins)

    return kernel.T  # (n_basins, n_time_steps)


# Example usage
df = pd.read_parquet('params.parquet')
dt = 3600.0

kernel = scs_curvilinear_kernel(dt=dt, tc=df['tc'].values, area=df['area_sqm'].values)
pd.DataFrame(kernel).to_parquet('kernel.parquet')
```

## Custom Unit Hydrograph

Any physically based unit hydrograph can be discretized into a kernel. The only requirement is that
each row sums to `1.0`. Here is an example using a linear reservoir UH:

```python
import numpy as np
import pandas as pd


def linear_reservoir_kernel(dt: float, k_values: np.ndarray) -> np.ndarray:
    """
    Instantaneous unit hydrograph for a linear reservoir: u(t) = (1/K) * exp(-t/K).

    Parameters
    ----------
    dt       : routing timestep in seconds
    k_values : 1D array of reservoir storage constants in seconds, one per basin

    Returns
    -------
    kernel : 2D array of shape (n_basins, n_time_steps), rows sum to 1.0
    """
    k_values = np.asarray(k_values, dtype=np.float64)
    max_steps = int(np.ceil(5 * k_values.max() / dt))

    t_edges = np.arange(max_steps + 1, dtype=np.float64)[:, np.newaxis] * dt
    u_cum = 1.0 - np.exp(-t_edges / k_values)       # (max_steps+1, n_basins)
    kernel = np.diff(u_cum, axis=0)                  # (max_steps, n_basins), sums to ~1.0

    return kernel.T  # (n_basins, n_time_steps)


k_values = pd.read_parquet('params.parquet')['k_reservoir'].values
kernel = linear_reservoir_kernel(dt=3600.0, k_values=k_values)
pd.DataFrame(kernel).to_parquet('kernel.parquet')
```

## Verifying Volume Conservation

After building a kernel, verify that each row integrates to the correct volume:

```python
import numpy as np
import pandas as pd

kernel = pd.read_parquet('kernel.parquet').to_numpy()  # (n_basins, n_time_steps)
dt = 3600.0

row_sums = kernel.sum(axis=1) * dt          # should equal area[j] for each basin
print(f'Min row integral: {row_sums.min():.4f}')
print(f'Max row integral: {row_sums.max():.4f}')

# Or if you normalize to a unitless kernel (rows sum to 1.0):
unitless = kernel / kernel.sum(axis=1, keepdims=True)
assert np.allclose(unitless.sum(axis=1), 1.0), 'Rows do not sum to 1.0'
```
