import glob
import logging
import os

import networkx as nx
import pandas as pd
import scipy
from natsort import natsorted
import json

logger = logging.getLogger(__name__)

__all__ = [
    'routing_files_from_RAPID',
    'connectivity_to_digraph',
    'connectivity_to_adjacency_matrix',
]


def routing_files_from_RAPID(riv_bas_id: str,
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


def job_file(
        routing_params_file: str = None,
        connectivity_file: str = None,
        runoff_file: str = None,
        outflow_file: str = None,
        adj_file: str = None,
        dt_routing: str = None,
        dt_outflows: str = None,
        min_q: str = None,
        initial_state_file: str = '',
        final_state_file: str = '',
        log: bool = True,
        log_stream: str = '',
        log_level: str = 'INFO',
        job_name: str = '',
        progress_bar: bool = True,
):
    """

    Args:
        routing_params_file: (string), Path to the routing parameters file.
        connectivity_file: (string), Path to the network connectivity file.
        runoff_file: (string), Path to the file with runoff values to be routed.
        outflow_file: (string), Path where the outflows file should be saved.
        adj_file: (string), Path where the adjacency matrix should be cached.
        dt_routing: (int), Time interval in seconds between routing computations.
        dt_outflows: (int), Time interval in seconds between writing flows to disc.
        min_q: (number), Minimum flow value allowed after each routing computation.
        initial_state_file: (string), Path to the file with initial state values.
        final_state_file: (string), Path to the file with final state values.
        log: (boolean), whether to display log messages defaulting to False
        log_stream: (string), the destination for logged messages: stdout stderr or a file path. default to stdout
        log_level: (string), Level of logging: either 'debug' 'info' 'warning' 'error' or 'critical'.
        job_name: (string), A name for this job to be printed in debug statements.
        progress_bar: (boolean), Indicates whether or not to show a progress bar in debug statements: true or false.

    Returns:

    """
    return


def job_files_from_directories(vpu_dir: str, inflows_dir: str, jobs_dir: str, states_dir: str, outputs_dir: str,
                               runoff_type: str = 'sequential', dt_routing: int = 900,
                               logs_dir: str = None, progress_bar: bool = True) -> None:
    """
    Generate river-route job files for all VPUs and for all runoff files in given directories
    """
    vpus = natsorted(glob.glob(os.path.join(vpu_dir, '*')))
    vpus = [os.path.basename(x) for x in vpus if os.path.isdir(x)]

    for vpu in vpus:
        runoffs = natsorted(glob.glob(os.path.join(inflows_dir, vpu, '*.nc')))
        # todo runoff_type == sequential, 1 start state for first runoff, 1 last state from last runoff end time

        if runoff_type == 'sequential':
            output_files = [os.path.join(outputs_dir, os.path.basename(r).replace('m3_', 'Qout_', 1)) for r in runoffs]
            start_state = os.path.join(states_dir, vpu, 'state_{date}.parquet')
            final_state = os.path.join(states_dir, vpu, 'state_{date}.parquet')

        for runoff in runoffs:
            # todo search for final state before start date
            # todo runoff_type == ensemble, 1 start state for each runoff, 1 last state for each runoff end time
            start_date = pd.to_datetime(os.path.basename(runoff).split('_')[2], format='%Y%m%d')
            start_state = os.path.join(states_dir, vpu, 'state_{date}.parquet')
            final_state = os.path.join(states_dir, vpu, 'state_{date}.parquet')

            job_file = os.path.join(jobs_dir, f'{vpu}_{runoff}.json')
            with open(job_file, 'w') as f:
                json.dumps({})
            job = {
                "routing_params_file": os.path.join(vpu_dir, vpu, 'params.parquet'),
                "connectivity_file": os.path.join(vpu_dir, vpu, 'connectivity.parquet'),
                "runoff_file": runoff,
                "outflow_file": os.path.join(outputs_dir, os.path.basename(runoff).replace('m3_', 'Qout_', 1)),
                "routing": "linear",
                "adj_file": os.path.join(vpu_dir, vpu, 'adj.npz'),
                "dt_routing": dt_routing,
                "min_q": 0,
                "start_state_file": start_state,
                "final_state_file": final_state,
                "log": True,
                "log_level": "DEBUG",
                "job_name": f'{vpu}_{start_date}_{runoff}',
                "progress_bar": True
            }
