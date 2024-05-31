from typing import List, Dict, Any
import json
from os.path import join, basename, isdir
import glob
from natsort import natsorted
import pandas as pd

__all__ = [
    'job_file',
]


def job_file(
        routing_params_file: str,
        connectivity_file: str,
        adj_file: str,
        runoff_file: str or List[str],
        outflow_file: str or List[str],
        dt_routing: int = 0,
        dt_outflows: int = 0,
        positive_flow: bool = True,
        routing: str = 'linear',
        nonlinear_routing_params_file: str = '',
        nonlinear_thresholds_file: str = '',
        initial_state_file: str = '',
        final_state_file: str = '',
        log: bool = True,
        log_stream: str = 'stdout',
        log_level: str = 'INFO',
        job_name: str = 'untitled_job',
        progress_bar: bool = True,
):
    """
    Constructs a dictionary with all necessary parameters for a river routing job

    Args:
        routing_params_file: (string), Path to the routing parameters file.
        connectivity_file: (string), Path to the network connectivity file.
        adj_file: (string), Path where the adjacency matrix should be cached.
        runoff_file: (string), Path to the file with runoff values to be routed.
        outflow_file: (string), Path where the outflows file should be saved.
        dt_routing: (int), Time interval in seconds between routing computations.
        dt_outflows: (int), Time interval in seconds between writing flows to disc.
        positive_flow: (boolean), Whether to enforce positive flows.
        routing: (string), Routing method to use: either 'linear' or 'nonlinear'.
        nonlinear_routing_params_file: (string), Path to the nonlinear routing parameters file.
        nonlinear_thresholds_file: (string), Path to the nonlinear routing thresholds file.
        initial_state_file: (string), Path to the file with initial state values.
        final_state_file: (string), Path to the file with final state values.
        log: (boolean), whether to display log messages defaulting to False
        log_stream: (string), the destination for logged messages: stdout stderr or a file path. default to stdout
        log_level: (string), Level of logging: either 'debug' 'info' 'warning' 'error' or 'critical'.
        job_name: (string), A name for this job to be printed in debug statements.
        progress_bar: (boolean), Indicates whether to show a progress bar in debug statements: true or false.

    Returns:
        (dict), A dictionary with parameters to be passed to MuskingumCunge class init.
    """
    if routing == 'nonlinear':
        assert nonlinear_routing_params_file and nonlinear_thresholds_file, \
            'Nonlinear routing requires a nonlinear routing parameters file and a nonlinear thresholds file'
    return {
        # Required Watershed Files
        'routing_params_file': routing_params_file,
        'connectivity_file': connectivity_file,
        'adj_file': adj_file,
        # Routing Input and Output
        'runoff_file': runoff_file,
        'outflow_file': outflow_file,
        # Compute Options - Optional
        'positive_flow': positive_flow,
        'dt_routing': dt_routing,
        'dt_outflows': dt_outflows,
        # Routing Method - Optional
        'routing': routing,
        'nonlinear_routing_params_file': nonlinear_routing_params_file,
        'nonlinear_thresholds_file': nonlinear_thresholds_file,
        # initial and final state files - Optional
        'initial_state_file': initial_state_file,
        'final_state_file': final_state_file,
        # simulation management and debugging - Optional
        'log': log,
        'progress_bar': progress_bar,
        'log_level': log_level,
        'log_stream': log_stream,
        'job_name': job_name,
    }


def job_files_from_directories(
        vpu_dir: str, inflows_dir: str, jobs_dir: str, outputs_dir: str, states_dir: str, logs_dir: str = None,
        initial_state_file: str = None, sim_type: str = 'sequential', dt_routing: int = 0, dt_outflows: int = 0,
        **kwargs
) -> List[Dict[str, Any]]:
    """
    Generate river-route job files for all VPUs and for all runoff files in given directories
    """
    vpus = natsorted(glob.glob(join(vpu_dir, '*')))
    vpus = [basename(x) for x in vpus if isdir(x)]

    job_files = []

    if sim_type == 'sequential':
        for vpu in vpus:
            routing_params_file = join(vpu_dir, vpu, 'params.parquet')
            connectivity_file = join(vpu_dir, vpu, 'connectivity.parquet')
            adj_file = join(vpu_dir, vpu, 'adj.npz')
            runvol_files = natsorted(glob.glob(join(inflows_dir, vpu, '*.nc')))
            output_files = [join(outputs_dir, basename(f).replace('m3_', 'Qout_')) for f in runvol_files]
            start_date = pd.to_datetime(basename(runvol_files[0]).split('_')[2], format='%Y%m%d')
            final_state = join(states_dir, vpu, f'finalstate_{vpu}_{start_date}.parquet')
            job_files.append(
                job_file(
                    routing_params_file=routing_params_file,
                    connectivity_file=connectivity_file,
                    adj_file=adj_file,
                    runoff_file=runvol_files,
                    outflow_file=output_files,
                    dt_routing=dt_routing,
                    dt_outflows=dt_outflows,
                    initial_state_file=initial_state_file,
                    final_state_file=final_state,
                    job_name=f'jobfile_{vpu}_{start_date}_{sim_type}',
                    **kwargs
                )
            )
    elif sim_type == 'ensemble':
        for vpu in vpus:
            routing_params_file = join(vpu_dir, vpu, 'params.parquet')
            connectivity_file = join(vpu_dir, vpu, 'connectivity.parquet')
            adj_file = join(vpu_dir, vpu, 'adj.npz')
            runvol_files = natsorted(glob.glob(join(inflows_dir, vpu, '*.nc')))
            output_files = [join(outputs_dir, basename(f).replace('m3_', 'Qout_')) for f in runvol_files]
            for runvol_file, output_file in zip(runvol_files, output_files):
                start_date = pd.to_datetime(basename(runvol_file).split('_')[2], format='%Y%m%d')
                final_state = join(states_dir, vpu, f'finalstate_{vpu}_{start_date}.parquet')
                job_files.append(
                    job_file(
                        routing_params_file=routing_params_file,
                        connectivity_file=connectivity_file,
                        adj_file=adj_file,
                        runoff_file=runvol_file,
                        outflow_file=output_file,
                        dt_routing=dt_routing,
                        dt_outflows=dt_outflows,
                        initial_state_file=initial_state_file,
                        final_state_file=final_state,
                        job_name=f'jobfile_{vpu}_{start_date}_{sim_type}',
                        **kwargs
                    )
                )

    for job in job_files:
        with open(join(jobs_dir, f'{job["job_name"]}.json'), 'w') as f:
            json.dump(job, f)

    return job_files
