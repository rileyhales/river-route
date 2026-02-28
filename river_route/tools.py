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
    Subset routing parameters and weight tables to only rivers including and upstream of a given id
    """
    pdf = pd.read_parquet(params)

    graph = connectivity_to_digraph(pdf['river_id'].values, pdf['downstream_river_id'].values)
    upstreams = list(nx.ancestors(graph, target_river))
    upstreams.append(target_river)

    pdf[pdf['river_id'].isin(upstreams)].to_parquet(out_params)
    pdf.loc[pdf['river_id'] == target_river, 'downstream_river_id'] = -1
    pdf.to_parquet(out_params)

    if weights is not None and out_weights is not None:
        with xr.open_dataset(weights) as ds:
            mask = np.isin(ds['river_id'].values, list(upstreams))
            ds.isel(index=mask).to_netcdf(out_weights)
    return


def connectivity_to_digraph(river_ids: np.ndarray, downstream_ids: np.ndarray) -> nx.DiGraph:
    """
    Generate directed graph from the routing parameters file
    """
    graph = nx.DiGraph()
    graph.add_edges_from(zip(river_ids, downstream_ids))
    return graph


def adjacency_matrix(river_ids: np.ndarray, downstream_ids: np.ndarray) -> scipy.sparse.csc_matrix:
    """
    Generate adjacency matrix from routing params file
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
