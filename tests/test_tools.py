"""Tests for river_route.tools — adjacency_matrix, connectivity_to_digraph, subset_configs_to_river."""
import os
import tempfile

import networkx as nx
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from conftest import RFSv2ConfigsData
from river_route.tools import adjacency_matrix, connectivity_to_digraph, subset_configs_to_river


# ═════════════════════════════════════════════════════════════════════════════
# adjacency_matrix
# ═════════════════════════════════════════════════════════════════════════════

def _networkx_adjacency_matrix(params: pd.DataFrame) -> sp.csc_matrix:
    graph = nx.DiGraph()
    graph.add_edges_from(params[['river_id', 'downstream_river_id']].values)
    nx_adj = nx.adjacency_matrix(graph, nodelist=params['river_id'].values).astype(np.float64)
    return sp.csc_matrix(nx_adj).T


def test_adjacency_matrix_matches_networkx(vpu: RFSv2ConfigsData):
    params = pd.read_parquet(vpu.rr2_params_file)
    nx_adj = _networkx_adjacency_matrix(params)
    adj = adjacency_matrix(params['river_id'].values, params['downstream_river_id'].values)
    diff = nx_adj - adj
    assert diff.nnz == 0, 'adjacency matrix does not match NetworkX reference'


def test_adjacency_matrix_shape(vpu: RFSv2ConfigsData):
    params = pd.read_parquet(vpu.rr2_params_file)
    n = len(params)
    adj = adjacency_matrix(params['river_id'].values, params['downstream_river_id'].values)
    assert adj.shape == (n, n)


def test_adjacency_matrix_outlets_have_no_outgoing_edges(vpu: RFSv2ConfigsData):
    params = pd.read_parquet(vpu.rr2_params_file)
    adj = adjacency_matrix(params['river_id'].values, params['downstream_river_id'].values)
    outlet_mask = params['downstream_river_id'].values == -1
    outlet_indices = np.where(outlet_mask)[0]
    # adj[downstream, upstream]=1, so an outlet's column (as upstream) should be all zeros
    for idx in outlet_indices:
        assert adj[:, idx].nnz == 0, f'Outlet at index {idx} has outgoing edges'


def test_adjacency_matrix_rejects_unsorted():
    # downstream index must be > upstream index (topological order)
    river_ids = np.array([10, 20, 30])
    downstream_ids = np.array([20, -1, 10])  # 30→10 violates order (index 2 → index 0)
    with pytest.raises(ValueError, match='topologically sorted'):
        adjacency_matrix(river_ids, downstream_ids)


def test_adjacency_matrix_rejects_unknown_downstream():
    river_ids = np.array([10, 20])
    downstream_ids = np.array([-1, 999])  # 999 not in river_ids
    with pytest.raises(ValueError, match='Unknown downstream_river_id'):
        adjacency_matrix(river_ids, downstream_ids)


# ═════════════════════════════════════════════════════════════════════════════
# connectivity_to_digraph
# ═════════════════════════════════════════════════════════════════════════════

def test_connectivity_to_digraph(vpu: RFSv2ConfigsData):
    params = pd.read_parquet(vpu.rr2_params_file)
    graph = connectivity_to_digraph(params['river_id'].values, params['downstream_river_id'].values)
    assert isinstance(graph, nx.DiGraph)
    assert graph.number_of_nodes() > 0
    assert graph.number_of_edges() > 0


def test_connectivity_to_digraph_simple():
    river_ids = np.array([1, 2, 3])
    downstream_ids = np.array([-1, 1, 1])
    graph = connectivity_to_digraph(river_ids, downstream_ids)
    assert graph.has_edge(2, 1)
    assert graph.has_edge(3, 1)
    assert graph.has_edge(1, -1)


# ═════════════════════════════════════════════════════════════════════════════
# subset_configs_to_river
# ═════════════════════════════════════════════════════════════════════════════

def test_subset_configs_to_river(vpu: RFSv2ConfigsData):
    """Subset to a known river; verify the target becomes the outlet and upstream rivers are included."""
    params = pd.read_parquet(vpu.rr2_params_file)
    # Pick a river that has upstream tributaries (not a headwater)
    target = params.loc[params['downstream_river_id'] == -1, 'river_id'].iloc[0]

    tmpdir = tempfile.mkdtemp()
    try:
        out_params = os.path.join(tmpdir, 'subset.parquet')
        subset_configs_to_river(target, str(vpu.rr2_params_file), out_params)

        sub = pd.read_parquet(out_params)
        assert target in sub['river_id'].values
        # Target river should be the outlet in the subset
        assert sub.loc[sub['river_id'] == target, 'downstream_river_id'].iloc[0] == -1
        # All subset rivers should exist in the original
        assert set(sub['river_id']).issubset(set(params['river_id']))
    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_subset_configs_to_river_with_weights(vpu: RFSv2ConfigsData):
    """Subset both params and grid weights; verify weight river_ids are a subset of params river_ids."""
    params = pd.read_parquet(vpu.rr2_params_file)
    target = params.loc[params['downstream_river_id'] == -1, 'river_id'].iloc[0]

    tmpdir = tempfile.mkdtemp()
    try:
        out_params = os.path.join(tmpdir, 'subset.parquet')
        out_weights = os.path.join(tmpdir, 'subset_weights.nc')
        subset_configs_to_river(target, str(vpu.rr2_params_file), out_params,
                                weights=str(vpu.grid_weights_file), out_weights=out_weights)

        import xarray as xr
        sub = pd.read_parquet(out_params)
        with xr.open_dataset(out_weights) as ds:
            weight_river_ids = set(ds['river_id'].values.tolist())
        assert weight_river_ids.issubset(set(sub['river_id'].values.tolist()))
    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)
