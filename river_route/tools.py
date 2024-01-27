import logging
import os

import networkx as nx
import pandas as pd
import scipy

logger = logging.getLogger(__name__)

__all__ = [
    'configs_from_rapid',
    'connectivity_to_digraph',
    'connectivity_to_adjacency_matrix',
]


def configs_from_rapid(riv_bas_id: str,
                       k: str,
                       x: str,
                       rapid_connect: str,
                       out_params: str,
                       out_connectivity: str, ) -> None:
    """
    Generate river-route configuration files from input files for RAPID

    Args:
        riv_bas_id: Path to riv_bas_id CSV file
        k: Path to k CSV file
        x: Path to x CSV file
        rapid_connect: Path to rapid_connect CSV file
        out_params: Path to output routing parameters parquet file
        out_connectivity: Path to output connectivity parquet file

    Returns:
        None
    """
    for f in [riv_bas_id, k, x, rapid_connect]:
        assert os.path.isfile(f), FileNotFoundError(f)

    pd.concat([
        pd.read_csv(riv_bas_id, header=None, names=['rivid']),
        pd.read_csv(x, header=None, names=['x']),
        pd.read_csv(k, header=None, names=['k']),
    ], axis=1).to_parquet(out_params)
    (
        pd
        .read_csv(rapid_connect, header=None)
        .iloc[:, :2]
        .rename(columns={0: 'rivid', 1: 'downstream_rivid'})
        .to_parquet(out_connectivity)
    )
    return


def connectivity_to_digraph(connectivity_file: str) -> nx.DiGraph:
    """
    Generate directed graph from connectivity file
    """
    G = nx.DiGraph()
    df = pd.read_parquet(connectivity_file)
    if len(df.columns) == 3:  # ID, DownstreamID, Weight
        # check that the weights are all positive and sum to 1 for each unique ID
        id_col, downstream_col, weight_col = df.columns
        weight_check = df.groupby(id_col)[weight_col].sum()
        if not weight_check.ge(0).all():
            logger.error(f'Weights are not all positive')
            logger.debug('The following IDs have negative weights')
            logger.debug(weight_check[~weight_check.ge(0)].index.tolist())
            raise ValueError(f'Weights must be positive')
        if not weight_check.eq(1).all():
            logger.error(f'Weights do not sum to 1.0 for each unique ID')
            logger.debug('The following IDs have weights that do not sum to 1.0')
            logger.debug(weight_check[~weight_check.eq(1)].index.tolist())
            raise ValueError(f'Weights must sum to 1 for each unique ID')
        G.add_weighted_edges_from(df.values)
    elif len(df.columns) == 2:  # ID, DownstreamID
        G.add_edges_from(df.values)
    else:
        raise ValueError(f'Connectivity file should have 2 or 3 columns, not {len(df.columns)}')
    return G


def connectivity_to_adjacency_matrix(connectivity_file: str) -> scipy.sparse.csc_matrix:
    """
    Generate adjacency matrix from connectivity file
    """
    G = connectivity_to_digraph(connectivity_file)
    sorted_order = list(nx.topological_sort(G))
    if -1 in sorted_order:
        sorted_order.remove(-1)
    return scipy.sparse.csc_matrix(nx.to_scipy_sparse_array(G, nodelist=sorted_order).T)
