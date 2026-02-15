from __future__ import annotations

import networkx as nx
import pandas as pd

from river_route.tools import adjacency_matrix

parameters_file = './data/sample_watershed/prepared_streams.parquet'
params = pd.read_parquet(parameters_file)
graph = nx.DiGraph()
graph.add_edges_from(params[['river_id', 'downstream_river_id']].values)
nx_adj_matrix = nx.adjacency_matrix(graph, nodelist=params['river_id'].values)
# todo write to scipy sparse
adj_matrix = adjacency_matrix(params['river_id'].values, params['downstream_river_id'].values)
# todo compare the two matrices
