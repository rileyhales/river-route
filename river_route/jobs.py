import glob
import json
from os.path import join, basename, isdir
from typing import List, Dict, Any

import pandas as pd
from natsort import natsorted

__all__ = [
    'job_configs',
    'job_files_from_directories',
]


def job_configs(
        routing_params_file: str,
        connectivity_file: str,
        runoff_volumes_file: str or List[str],
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
        progress_bar: bool = True,
        job_name: str = 'untitled_job',
        log_level: str = 'INFO',
        log_stream: str = 'stdout',
        var_runoff_volumes: str = 'ro_vol',
        var_river_id: str = 'river_id',
        var_discharge: str = 'Q',
):
    """
    Constructs a dictionary with all necessary parameters for a river routing job

    Args:
        routing_params_file: (string), Path to the routing parameters file.
        connectivity_file: (string), Path to the network connectivity file.

        runoff_volumes_file: (string), Path to the file with runoff values to be routed.
        outflow_file: (string), Path where the outflows file should be saved.

        dt_routing: (int), Time interval in seconds between routing computations.
        dt_outflows: (int), Time interval in seconds between writing flows to disc.
        positive_flow: (Bool), Whether to enforce positive flows.

        routing: (str), Routing method to use: either 'linear' or 'nonlinear'.
        nonlinear_routing_params_file: (str), Path to the nonlinear routing parameters file.
        nonlinear_thresholds_file: (str), Path to the nonlinear routing thresholds file.

        initial_state_file: (str), Path to the file with initial state values.
        final_state_file: (str), Path to the file with final state values.

        log: (Bool), whether to display log messages defaulting to False.
        progress_bar: (Bool), whether to show a progress bar in debug statements: true or false.
        job_name: (str), A name for this job to be printed in debug statements.
        log_level: (str), Level of logging: either 'debug' 'info' 'warning' 'error' or 'critical'.
        log_stream: (str), the destination for logged messages: stdout stderr or a file path. default to stdout.

        var_runoff_volumes: (str), the name of the runoff volumes variable in runoff_volumes_file
        var_river_id: (str), the name of the river id variable in all files
        var_discharge: (str), the name of the discharge variable in outflow_file
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
        # Routing Input and Output
        'runoff_file': runoff_volumes_file,
        'outflow_file': outflow_file,
        # Timestep Options - Optional
        'dt_routing': dt_routing,
        'dt_outflows': dt_outflows,
        # Routing Method - Optional
        'routing': routing,
        'positive_flow': positive_flow,
        'nonlinear_routing_params_file': nonlinear_routing_params_file,
        'nonlinear_thresholds_file': nonlinear_thresholds_file,
        # Initial and Final State files - Optional
        'initial_state_file': initial_state_file,
        'final_state_file': final_state_file,
        # Logging, Management, and Debugging - Optional
        'log': log,
        'progress_bar': progress_bar,
        'job_name': job_name,
        'log_level': log_level,
        'log_stream': log_stream,
        # Variable Names - Optional
        'var_runoff_volumes': var_runoff_volumes,
        'var_river_id': var_river_id,
        'var_discharge': var_discharge,
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

    jobs = []

    for vpu in vpus:
        routing_params_file = join(vpu_dir, vpu, 'params.parquet')
        connectivity_file = join(vpu_dir, vpu, 'connectivity.parquet')
        adj_file = join(vpu_dir, vpu, 'adj.npz')
        runvol_files = natsorted(glob.glob(join(inflows_dir, vpu, '*.nc')))
        output_files = [join(outputs_dir, basename(f).replace('m3_', 'Qout_')) for f in runvol_files]

        if sim_type == 'sequential':
            start_date = pd.to_datetime(basename(runvol_files[0]).split('_')[2], format='%Y%m%d')
            final_state = join(states_dir, vpu, f'finalstate_{vpu}_{start_date}.parquet')
            jobs.append(
                job_configs(routing_params_file=routing_params_file, connectivity_file=connectivity_file,
                            runoff_volumes_file=runvol_files, outflow_file=output_files, dt_routing=dt_routing,
                            dt_outflows=dt_outflows, initial_state_file=initial_state_file,
                            final_state_file=final_state, job_name=f'job_{vpu}_{start_date}_{sim_type}', **kwargs)
            )

        elif sim_type == 'ensemble':
            for runvol_file, output_file in zip(runvol_files, output_files):
                start_date = pd.to_datetime(basename(runvol_file).split('_')[2], format='%Y%m%d')
                final_state = join(states_dir, vpu, f'finalstate_{vpu}_{start_date}.parquet')
                jobs.append(
                    job_configs(routing_params_file=routing_params_file, connectivity_file=connectivity_file,
                                runoff_volumes_file=runvol_file, outflow_file=output_file, dt_routing=dt_routing,
                                dt_outflows=dt_outflows, initial_state_file=initial_state_file,
                                final_state_file=final_state, job_name=f'job_{vpu}_{start_date}_{sim_type}', **kwargs)
                )

    for job in jobs:
        with open(join(jobs_dir, f'{job["job_name"]}.json'), 'w') as f:
            json.dump(job, f)

    return jobs
