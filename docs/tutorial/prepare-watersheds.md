## Generating Channel Routing Configuration Files

Configuration files (YAML or JSON) and the routing parameters parquet are the two static inputs that
describe your river network. This tutorial shows how to build them programmatically.

## Routing Parameters File

The routing parameters parquet has one row per river segment. Required columns:

| Column                | Type    | Description                                                           |
|-----------------------|---------|-----------------------------------------------------------------------|
| `river_id`            | int64   | Unique integer ID for each segment                                    |
| `downstream_river_id` | int64   | ID of the downstream segment; `-1` (or any negative value) at outlets |
| `k`                   | float64 | Muskingum K — wave travel time through the segment in seconds         |
| `x`                   | float64 | Muskingum X — attenuation factor, 0 ≤ x ≤ 0.5                         |

Rows **must be in topological order** — every upstream segment before its downstream neighbor.

If you are using `UnitMuskingum`, the parameters file should also include catchment attributes needed
to build the UH kernel (e.g., `tc` and `area_sqm`). See
[Generating Unit Hydrograph Kernels](create-uh-kernels.md).

### Building from a GeoDataFrame

```python
import geopandas as gpd
import pandas as pd

rivers = gpd.read_file('rivers.gpkg')

# Compute K from channel length and estimated wave speed
rivers['k'] = rivers['length_m'] / rivers['wave_speed_m_s']

# X is dimensionless; start with 0.3 as a first guess
rivers['x'] = 0.3

# Topological sort: upstream segments first
from river_route import tools

df = rivers[['river_id', 'downstream_river_id', 'k', 'x']].copy()

# Ensure topological ordering (river-route requires this)
# Use a graph library or iterative upstream-first sort
import networkx as nx

G = nx.DiGraph()
for _, row in df.iterrows():
    if row['downstream_river_id'] > 0:
        G.add_edge(row['river_id'], row['downstream_river_id'])
    else:
        G.add_node(row['river_id'])

ordered_ids = list(nx.topological_sort(G))
df = df.set_index('river_id').loc[ordered_ids].reset_index()

df.to_parquet('params.parquet', index=False)
```

### Muskingum K Estimation

K is the wave travel time in seconds. A first estimate comes from channel geometry:

```python
# Manning-based wave speed
# v_wave ≈ (5/3) * v_mean  for wide rectangular channel
# v_mean = (1/n) * R^(2/3) * S^(1/2)
import numpy as np

n_manning = 0.035  # Manning's roughness (-)
slope = rivers['slope']  # channel slope (m/m)
hydraulic_radius = 1.0  # rough estimate (m) for bankfull

v_mean = (1 / n_manning) * (hydraulic_radius ** (2 / 3)) * np.sqrt(slope)
v_wave = (5 / 3) * v_mean  # kinematic wave speed

rivers['k'] = rivers['length_m'] / v_wave
```

## Config File

Config values may be passed as a YAML/JSON file, as keyword arguments, or both. The most portable
approach for production pipelines is to generate them programmatically.

### Generating a Config Programmatically

```python
import glob
import os

import river_route as rr

root_dir = '/path/to/root/directory'
vpu_name = 'sample-vpu'

params_file = os.path.join(root_dir, 'configs', vpu_name, 'params.parquet')

volume_files = sorted(glob.glob(os.path.join(root_dir, 'volumes', vpu_name, '*.nc')))

output_dir = os.path.join(root_dir, 'outputs', vpu_name)
os.makedirs(output_dir, exist_ok=True)
output_files = [os.path.join(output_dir, f'Qout_{os.path.basename(f)}') for f in volume_files]

(
    rr
    .RapidMuskingum(
        routing_params_file=params_file,
        lateral_volume_files=volume_files,
        discharge_files=output_files,
        dt_routing=3600,
    )
    .route()
)
```

### Writing a YAML Config File

```python
import yaml

config = {
    'routing_params_file': '/path/to/params.parquet',
    'lateral_volume_files': ['/path/to/volumes_jan.nc', '/path/to/volumes_feb.nc'],
    'discharge_files': ['/path/to/discharge_jan.nc', '/path/to/discharge_feb.nc'],
    'dt_routing': 3600,
    'channel_state_init_file': '/path/to/initial_state.parquet',
    'channel_state_final_file': '/path/to/final_state.parquet',
}

with open('config.yaml', 'w') as f:
    yaml.dump(config, f, default_flow_style=False)
```

## Grid Weights for Gridded Runoff

If your runoff data is on a regular grid rather than pre-aggregated to catchments, `river-route` can
aggregate it using a weight table. Weights are computed once from the intersection of grid cells with
catchment boundaries and reused across all routing runs on the same network and grid.

```python
# Conceptual steps — use your preferred GIS/intersection tool
# 1. Intersect grid cell polygons with catchment polygons
# 2. For each (grid_cell, catchment) pair, record the intersection area
# 3. Save as a parquet with columns: river_id, row, col, weight (area fraction)

grid_weights_df.to_parquet('grid_weights.parquet')
```

Then reference it in the config:

```yaml
routing_params_file: 'params.parquet'
runoff_depth_grids: 'runoff.nc'
grid_weights_file: 'grid_weights.parquet'
discharge_files: 'discharge.nc'
dt_routing: 3600
```

Weights need to be recomputed if the grid resolution, grid extent, or catchment boundaries change.
