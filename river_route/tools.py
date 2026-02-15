import logging
import os

import networkx as nx
import numpy as np
import pandas as pd
import scipy
import xarray as xr

logger = logging.getLogger(__name__)

__all__ = [
    'routing_files_from_RAPID',
    'connectivity_to_digraph',
    'get_adjacency_matrix',
]


def routing_files_from_RAPID(
        riv_bas_id: str,
        k: str,
        x: str,
        rapid_connect: str,
        out_params: str,
) -> None:
    """
    Generate river-route configuration files from input files for RAPID

    Args:
        riv_bas_id: Path to riv_bas_id CSV file
        k: Path to k CSV file
        x: Path to x CSV file
        rapid_connect: Path to rapid_connect CSV file
        out_params: Path to output routing parameters parquet file

    Returns:
        None
    """
    for f in [riv_bas_id, k, x, rapid_connect]:
        assert os.path.isfile(f), FileNotFoundError(f)

    params = pd.concat([
        pd.read_csv(riv_bas_id, header=None, names=['river_id']),
        pd.read_csv(x, header=None, names=['x']),
        pd.read_csv(k, header=None, names=['k']),
    ], axis=1)
    connectivity = (
        pd
        .read_csv(rapid_connect, header=None)
        .iloc[:, :2]
        .rename(columns={0: 'river_id', 1: 'downstream_river_id'})
    )
    params = params.merge(connectivity, on='river_id', how='left')
    if params['downstream_river_id'].isna().any():
        raise ValueError('Some river_id values in params were not found in rapid_connect')
    params['downstream_river_id'] = params['downstream_river_id'].astype(int)
    params.to_parquet(out_params)
    return


def inits_from_RAPID(
        init_file: str,
        out_init_file: str,
):
    """
    Generate river-route initial conditions file from RAPID final state file

    Args:
        init_file: Path to RAPID init file
        out_init_file: Path to output initial conditions file

    Returns:
        None
    """
    assert os.path.isfile(init_file), FileNotFoundError(init_file)
    (
        xr
        .open_dataset(init_file)
        .to_dataframe()
        .reset_index()
        [['rivid', 'Qout']]
        .rename(columns={'Qout': 'Q', 'rivid': 'river_id'})
        .set_index('river_id')
        .assign(R=0)
        .to_parquet(out_init_file)
    )
    return


def subset_configs_to_river(
        target_river: int,
        params: str,
        out_params: str,
        connectivity: str,
        out_connectivity: str,
        weights: str = None,
        out_weights: str = None,
) -> None:
    """
    Subset routing parameters, connectivity, and weight tables to only river upstream of a given id
    """
    pdf = pd.read_parquet(params)
    cdf = pd.read_parquet(connectivity)

    graph = connectivity_to_digraph(pdf['river_id'], pdf['downstream_river_id'])
    upstreams = list(nx.ancestors(graph, target_river))
    upstreams.append(target_river)

    pdf[pdf['river_id'].isin(upstreams)].to_parquet(out_params)
    cdf = cdf[cdf['river_id'].isin(upstreams)]
    cdf.loc[cdf['river_id'] == target_river, 'downstream_river_id'] = -1
    cdf.to_parquet(out_connectivity)

    if weights is not None and out_weights is not None:
        wdf = pd.read_csv(weights)
        wdf[wdf['river_id'].isin(upstreams)].to_csv(out_weights, index=False)
    return


def connectivity_to_digraph(river_ids: np.ndarray, downstream_ids: np.ndarray) -> nx.DiGraph:
    """
    Generate directed graph from the routing parameters file
    """
    graph = nx.DiGraph()
    graph.add_edges_from(zip(river_ids, downstream_ids))
    return graph


def get_adjacency_matrix(routing_params_file: str) -> scipy.sparse.csc_matrix:
    """
    Generate adjacency matrix from routing params file
    """
    params = pd.read_parquet(routing_params_file)
    river_id_col = params.columns[0]
    if 'downstream_river_id' not in params.columns:
        raise ValueError('routing_params_file must include a downstream river ID column')
    river_ids = params[river_id_col].to_numpy(dtype=np.int64, copy=False)
    downstream_ids = params['downstream_river_id'].to_numpy(dtype=np.int64, copy=False)
    river_index = {int(river_id): idx for idx, river_id in enumerate(river_ids.tolist())}

    row_indices: list[int] = []
    col_indices: list[int] = []
    for upstream_idx, downstream_river_id in enumerate(downstream_ids.tolist()):
        if downstream_river_id == -1:
            continue
        if downstream_river_id not in river_index:
            raise ValueError(f'Unknown downstream_river_id: {downstream_river_id}')
        downstream_idx = river_index[int(downstream_river_id)]
        if downstream_idx <= upstream_idx:
            raise ValueError('routing_params_file must be topologically sorted upstream to downstream')
        row_indices.append(downstream_idx)
        col_indices.append(upstream_idx)

    data = np.ones(len(row_indices), dtype=np.float64)
    return scipy.sparse.csc_matrix((data, (row_indices, col_indices)), shape=(river_ids.shape[0], river_ids.shape[0]))
