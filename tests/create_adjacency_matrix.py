from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse as sp

from river_route.tools import adjacency_matrix


def networkx_adjacency_matrix(params: pd.DataFrame) -> sp.csc_matrix:
    graph = nx.DiGraph()
    graph.add_edges_from(params[['river_id', 'downstream_river_id']].values)
    nx_adj_matrix = nx.adjacency_matrix(graph, nodelist=params['river_id'].values).astype(np.float64)
    return sp.csc_matrix(nx_adj_matrix).T


def run_adjacency_matrix_comparison() -> None:
    data_dir = Path(__file__).resolve().parent / 'data' / 'sample_watershed'
    parameters_file = data_dir / 'prepared_streams.parquet'
    params = pd.read_parquet(parameters_file)

    nx_adj_matrix = networkx_adjacency_matrix(params)
    adj_matrix = adjacency_matrix(params['river_id'].values, params['downstream_river_id'].values)

    diff = nx_adj_matrix - adj_matrix
    assert diff.nnz == 0, 'adjacency matrix does not match NetworkX reference'
    print('[OK] adjacency matrix matches NetworkX reference (scipy sparse)')

    # compute each 50 times and average the compute times
    speeds_nx = []
    speeds_rr = []
    for _ in range(50):
        start_time = pd.Timestamp.now()
        networkx_adjacency_matrix(params)
        speeds_nx.append((pd.Timestamp.now() - start_time).total_seconds())
        start_time = pd.Timestamp.now()
        adjacency_matrix(params['river_id'].values, params['downstream_river_id'].values)
        speeds_rr.append((pd.Timestamp.now() - start_time).total_seconds())

    print(f'[TIMING] NetworkX adjacency matrix computation time: {np.mean(speeds_nx):.4f} seconds')
    print(f'[TIMING] river_route adjacency matrix computation time: {np.mean(speeds_rr):.4f} seconds')


if __name__ == '__main__':
    run_adjacency_matrix_comparison()
