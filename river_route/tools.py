import logging

import networkx as nx
import numpy as np
import pandas as pd
import scipy
import xarray as xr

from .types import PathInput

logger = logging.getLogger(__name__)

__all__ = [
    'subset_configs_to_river',
    'connectivity_to_digraph',
    'adjacency_matrix',
]


def subset_configs_to_river(
        target_river: int,
        params: PathInput,
        out_params: PathInput,
        weights: PathInput | None = None,
        out_weights: PathInput | None = None,
) -> None:
    """
    Subset routing parameters and weight tables to the target river and all rivers upstream of it.

    The target river becomes the outlet (downstream_river_id set to -1) in the subset. If weight
    table paths are provided, the weight table is also filtered to only include matching rivers.

    Args:
        target_river: river_id of the river to subset to (becomes the outlet)
        params: path to the full routing parameters parquet file
        out_params: path to write the subsetted parameters parquet file
        weights: path to the full grid weights netCDF file (optional)
        out_weights: path to write the subsetted grid weights netCDF file (optional, required if weights is given)
    """
    pdf = pd.read_parquet(params)

    graph = connectivity_to_digraph(pdf['river_id'].values, pdf['downstream_river_id'].values)
    upstreams = list(nx.ancestors(graph, target_river))
    upstreams.append(target_river)

    subset = pdf[pdf['river_id'].isin(upstreams)].copy()
    subset.loc[subset['river_id'] == target_river, 'downstream_river_id'] = -1
    subset.to_parquet(out_params)

    if weights is not None and out_weights is not None:
        with xr.open_dataset(weights) as ds:
            mask = np.isin(ds['river_id'].values, list(upstreams))
            ds.isel(index=mask).to_netcdf(out_weights)
    return


def connectivity_to_digraph(river_ids: np.ndarray, downstream_ids: np.ndarray) -> nx.DiGraph:
    """
    Build a NetworkX DiGraph from river connectivity arrays.

    Each edge goes from a river to its downstream river (including the -1 sentinel for outlets).

    Args:
        river_ids: 1D array of river ID integers
        downstream_ids: 1D array of downstream river ID integers (-1 for outlets)

    Returns:
        Directed graph with edges from each river_id to its downstream_river_id
    """
    graph = nx.DiGraph()
    graph.add_edges_from(zip(river_ids, downstream_ids))
    return graph


def adjacency_matrix(river_ids: np.ndarray, downstream_ids: np.ndarray) -> scipy.sparse.csc_matrix:
    """
    Build a sparse adjacency matrix for the river network.

    Entry A[downstream_idx, upstream_idx] = 1 for each river that flows into a downstream river.
    Outlet rivers (downstream_id < 0) have no outgoing edges. The input arrays must be topologically
    sorted (upstream before downstream) — a ValueError is raised otherwise.

    Args:
        river_ids: 1D array of river ID integers, topologically sorted upstream to downstream
        downstream_ids: 1D array of downstream river ID integers (-1 for outlets)

    Returns:
        Sparse CSC matrix of shape (n, n) where n = len(river_ids)

    Raises:
        ValueError: if a downstream_id is not found in river_ids, or if the arrays are not
            topologically sorted
    """
    river_index = {int(river_id): idx for idx, river_id in enumerate(river_ids.tolist())}
    row_indices: list[int] = []
    col_indices: list[int] = []
    for upstream_idx, downstream_river_id in enumerate(downstream_ids.tolist()):
        if downstream_river_id < 0:
            continue
        if downstream_river_id not in river_index:
            raise ValueError(f'Unknown downstream_river_id: {downstream_river_id}')
        downstream_idx = river_index[int(downstream_river_id)]
        if downstream_idx <= upstream_idx:
            raise ValueError('params_file must be topologically sorted upstream to downstream')
        row_indices.append(downstream_idx)
        col_indices.append(upstream_idx)

    data = np.ones(len(row_indices), dtype=np.float64)
    return scipy.sparse.csc_matrix((data, (row_indices, col_indices)), shape=(river_ids.shape[0], river_ids.shape[0]))
