"""
Convert river-route v2 inputs to v3 inputs. v3 requires the rivers of the params file in depth-first search (DFS)
order: every river's upstream rivers are the rows immediately before it, so each river is the last row of the
contiguous range of rows holding its whole upstream watershed. v2 only required upstream before downstream, so a v2
params file that is not in DFS order is reordered, and every file with one row per river (grid weights, qlateral
files, channel state files) is converted in that same river order so they stay aligned with the params file.

    python examples/migrate_v2_to_v3.py \
        --params v2/params.parquet --out-params v3/params.parquet \
        --weights v2/weights.nc --out-weights v3/weights.nc \
        --qlateral v2/qlateral_2020.nc v2/qlateral_2021.nc --out-catchment-runoff-dir v3/ \
        --channel-state v2/state.parquet --out-channel-state v3/state.parquet \
        --config v2/config.yaml --out-config v3/config.json
"""

import argparse
import json
from dataclasses import fields
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import river_route as rr

V2_UNIT_HYDROGRAPH_KEYS = ('uh_kernel_file', 'uh_state_init_file', 'uh_state_final_file')


def dfs_order(river_ids: np.ndarray, next_river_ids: np.ndarray) -> np.ndarray:
    """
    The positions of the rivers in DFS order: each river after all of its upstream rivers, which are the rows
    immediately before it. Upstream rivers and basin outlets keep their original relative order.

    river_ids: (n,) id of each river
    next_river_ids: (n,) id of the river each drains into, negative at a basin outlet
    """
    n = river_ids.shape[0]
    if next_river_ids.shape != (n,):
        raise ValueError('river_ids and next_river_ids must have the same length')
    position = pd.Series(np.arange(n), index=river_ids)
    if not position.index.is_unique:
        raise ValueError('river_id values must be unique')
    is_outlet = next_river_ids < 0
    downstream = np.full(n, n, dtype=np.int64)  # n is a root that every basin outlet drains into
    downstream[~is_outlet] = position.reindex(next_river_ids[~is_outlet]).to_numpy(dtype=float, na_value=-1)
    if np.any(downstream < 0):
        raise ValueError('next_river_id values reference ids not present in river_id')

    upstream = np.argsort(downstream, kind='stable')  # every river, grouped by the river it drains into
    indptr = np.concatenate(([0], np.cumsum(np.bincount(downstream, minlength=n + 1))))
    next_upstream = indptr[:-1].copy()
    stack = np.empty(n + 1, dtype=np.int64)
    stack[0], depth = n, 1
    order = np.empty(n, dtype=np.int64)
    emitted = 0
    # each river is pushed onto the stack once and emitted once, and the root is pushed and popped once
    for _ in range(2 * n + 2):
        if depth == 0:
            break
        river = stack[depth - 1]
        if next_upstream[river] < indptr[river + 1]:
            stack[depth] = upstream[next_upstream[river]]
            next_upstream[river] += 1
            depth += 1
            continue
        depth -= 1
        if river != n:
            order[emitted] = river
            emitted += 1
    if emitted != n:
        raise ValueError(f'{n - emitted:,} rivers are not upstream of any basin outlet; the network has a cycle')
    return order


def convert_params(params_path: str, output_path: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Rename downstream_river_id to next_river_id and write the rivers in DFS order.
    Returns the river ids in their v2 order and in their v3 order.
    """
    params = pd.read_parquet(params_path).rename(columns={'downstream_river_id': 'next_river_id'})
    missing = {'river_id', 'next_river_id', 'k', 'x'} - set(params.columns)
    if missing:
        raise ValueError(f'{params_path} is missing required columns: {missing}')
    v2_river_ids = params['river_id'].to_numpy()
    if not rr.network.streams.is_dfs_ordered(rr.Network(params).downstream_indices).all():
        params = params.iloc[dfs_order(v2_river_ids, params['next_river_id'].to_numpy())].reset_index(drop=True)
    params.to_parquet(output_path, index=False)
    return v2_river_ids, params['river_id'].to_numpy()


def positions_of(river_ids: np.ndarray, file_river_ids: np.ndarray, path: str) -> np.ndarray:
    """The position in a file of each river of river_ids, which must be exactly the rivers of the file."""
    if np.setxor1d(river_ids, file_river_ids).size:
        raise ValueError(f'{path} does not have exactly the rivers of the params file')
    return pd.Series(np.arange(file_river_ids.shape[0]), index=file_river_ids).loc[river_ids].to_numpy()


def convert_grid_weights(weights_path: str, river_ids: np.ndarray, output_path: str) -> np.ndarray:
    """
    Write the weight table with its rows in the river order of the v3 params file, the order each river's runoff is
    read in. Returns the catchment area of each river in m².
    """
    with xr.open_dataset(weights_path) as ds:
        weights = ds.load()
    rows = weights['river_id'].to_numpy()
    first_row = pd.Series(np.arange(rows.shape[0])).groupby(rows).first()
    rank = pd.Series(np.arange(river_ids.shape[0]), index=river_ids).reindex(first_row.index)
    if rank.isna().any() or first_row.shape[0] != river_ids.shape[0]:
        raise ValueError(f'{weights_path} does not have exactly the rivers of the params file')
    weights.isel(index=np.argsort(rank.loc[rows].to_numpy(), kind='stable')).to_netcdf(output_path)
    catchment_area = weights[['river_id', 'area_sqm']].to_dataframe().groupby('river_id')['area_sqm'].sum()
    return catchment_area.loc[river_ids].to_numpy()


def convert_qlateral(qlateral_path: str, river_ids: np.ndarray, catchment_area: np.ndarray, output_path: str) -> None:
    """
    Write a v2 qlateral file, the runoff volume (m³) of each catchment per step with dimensions (time, river_id), as
    a v3 catchment runoff file with dimensions (river_id, time) in the river order of the v3 params file.
    catchment_area is only read when catchment runoff holds depths, so it may be NaN for these volumes.
    """
    with xr.open_dataset(qlateral_path) as ds:
        rows = positions_of(river_ids, ds['river_id'].to_numpy(), qlateral_path)
        rr.CatchmentRunoff().to_netcdf(
            output_path,
            dates=ds['time'].to_numpy(),
            catchment_runoff=ds['qlateral'].transpose('river_id', 'time').to_numpy()[rows],
            river_ids=river_ids,
            catchment_area=catchment_area,
        )


def convert_channel_state(state_path: str, v2_river_ids: np.ndarray, river_ids: np.ndarray, output_path: str) -> None:
    """Write a channel state file, one Q per river in v2 params file order, in the river order of the v3 file."""
    state = pd.read_parquet(state_path, columns=['Q'])
    if state.shape[0] != v2_river_ids.shape[0]:
        raise ValueError(f'{state_path} has {state.shape[0]} rows but the params file has {v2_river_ids.shape[0]}')
    state.iloc[positions_of(river_ids, v2_river_ids, state_path)].to_parquet(output_path, index=False)


def convert_config(config_path: str, output_path: str) -> None:
    """
    Write a v2 YAML or JSON config as a v3 JSON config. The forcing is the form of the runoff files the v2 config
    routed. Paths are copied unchanged, so point them at the converted files. Reading YAML needs pyyaml.
    """
    with open(config_path, encoding='utf-8') as f:
        if Path(config_path).suffix in ('.yml', '.yaml'):
            import yaml

            config = yaml.safe_load(f)
        else:
            config = json.load(f)
    unit_hydrograph = [key for key in V2_UNIT_HYDROGRAPH_KEYS if config.pop(key, None)]
    if unit_hydrograph:
        raise ValueError(f'{config_path} sets {unit_hydrograph}: the unit hydrograph transform is not in v3')
    qlateral_files = config.pop('qlateral_files', None) or []
    grid_runoff_files = config.pop('grid_runoff_files', None) or []
    if qlateral_files and grid_runoff_files:
        raise ValueError(f'{config_path} sets both qlateral_files and grid_runoff_files')
    config['forcing'] = 'catchment' if qlateral_files else 'grid' if grid_runoff_files else 'channel'
    config['runoff_files'] = qlateral_files or grid_runoff_files
    unknown = set(config) - {field.name for field in fields(rr.Configs) if field.init}
    if unknown:
        raise ValueError(f'{config_path} has keys that are not v2 configs: {sorted(unknown)}')
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(config, f, indent=2)


def main() -> None:
    parser = argparse.ArgumentParser(description='Convert river-route v2 inputs to v3 inputs')
    parser.add_argument('--params', help='v2 params parquet')
    parser.add_argument('--out-params', help='v3 params parquet to write, in DFS order')
    parser.add_argument('--weights', help='v2 grid weights netCDF')
    parser.add_argument('--out-weights', help='v3 grid weights netCDF to write')
    parser.add_argument('--qlateral', nargs='+', default=[], help='v2 qlateral netCDF files')
    parser.add_argument('--out-catchment-runoff-dir', help='directory to write the v3 catchment runoff files in')
    parser.add_argument('--channel-state', help='v2 channel state parquet')
    parser.add_argument('--out-channel-state', help='v3 channel state parquet to write')
    parser.add_argument('--config', help='v2 YAML or JSON config')
    parser.add_argument('--out-config', help='v3 JSON config to write')
    args = parser.parse_args()

    pairs = {
        '--params': (args.params, args.out_params),
        '--weights': (args.weights, args.out_weights),
        '--qlateral': (args.qlateral, args.out_catchment_runoff_dir),
        '--channel-state': (args.channel_state, args.out_channel_state),
        '--config': (args.config, args.out_config),
    }
    for flag, (source, output) in pairs.items():
        if bool(source) != bool(output):
            parser.error(f'{flag} and its output are given together: one is useless without the other')
    if (args.weights or args.qlateral or args.channel_state) and not args.params:
        parser.error('--params is required to put the files with one row per river in the v3 river order')

    if args.config:
        convert_config(args.config, args.out_config)
    if not args.params:
        return
    v2_river_ids, river_ids = convert_params(args.params, args.out_params)
    catchment_area = np.full(river_ids.shape[0], np.nan, dtype=np.float32)
    if args.weights:
        catchment_area = convert_grid_weights(args.weights, river_ids, args.out_weights)
    for qlateral_path in args.qlateral:
        output_path = Path(args.out_catchment_runoff_dir) / f'catchment_runoff_{Path(qlateral_path).stem}.nc'
        convert_qlateral(qlateral_path, river_ids, catchment_area, str(output_path))
    if args.channel_state:
        convert_channel_state(args.channel_state, v2_river_ids, river_ids, args.out_channel_state)


if __name__ == '__main__':
    main()
