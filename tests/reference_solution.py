"""
The reference solution: one river basin routed with several months of ERA5 runoff, one file per month, packaged with
every input and output so the test suite can require river-route to regenerate it exactly. It is the correct answer by
definition, so rebuild it only when the routing is deliberately changed. Build it with

    python tests/reference_solution.py \
        --region-dir ~/data/rfsv3/routing/region=7020014250 \
        --outlet-river-id 720207783 \
        --runoff ~/data/era5/year=2000/era5_2000{03,04,05,06}.nc

which routes the Columbia River, river 720207783 at its mouth and the 29,404 rivers upstream of it in the level 2
basin covering Washington, Oregon, and southwestern Canada, through March to June 2000. It cuts the region's network
file, network.parquet, down to that basin with ``streams.subset_network_to_river``, and its grid weights to the same
rivers with the river id renamed from the region's river_id to riverId, copies the runoff files unmodified into
``inputs/``,
aggregates each runoff file into a catchment runoff file, routes every run of ``reference_runs`` into ``outputs/`` with
the months in sequence, and writes ``manifest.json`` describing how each output was made. Discharge is written by the
default zarr writer with every float32 mantissa bit kept, so nothing is rounded.

numba compiles for the machine it runs on, so a package built on one CPU architecture may differ in the last bit on
another, where fused multiply-adds are formed differently. The manifest records the machine it was built on, and the
tests compare to ``TOLERANCE`` rather than bit for bit, so they hold on any machine.

The module also holds what the tests share: ``assert_same`` for answers that should be identical, and ``route`` for
routing runoff read from files or shaped in memory and keeping what the router hands its discharge writer.
"""

import argparse
import contextlib
import json
import platform
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numba
import numpy as np
import pandas as pd
import xarray as xr
import zarr

import river_route as rr
from river_route.network import streams
from river_route.router import writers

PACKAGE = Path(__file__).parent / 'data' / 'reference_solution'
GRID_NAMES = {'var_x': 'longitude', 'var_y': 'latitude', 'var_t': 'valid_time'}
REGION_VAR_RIVER_ID = 'river_id'  # the river id of the region's grid weights, which the package names riverId
INPUT_PATH_OPTIONS = ('network_file', 'grid_weights_file', 'channel_state_init_file')
TOLERANCE = 10e-4  # answers that should be identical may differ by this share of each river's largest value
WILLAMETTE = 720184033  # the Willamette River where it joins the Columbia at Portland: 1,390 rivers


@dataclass
class Basin:
    """A basin cut from the package: its network file and weights, and its catchment runoff volumes for each month."""

    network_file: Path
    weights_file: Path
    runoff_files: list[Path]  # the package's gridded runoff, one file per month
    months: list[tuple[np.ndarray, np.ndarray]]  # (dates, (river, time) catchment runoff volumes) per runoff file
    network: rr.Network  # for reading the topology; build a fresh one for anything that changes it

    def gridded(self, months: int = 4) -> dict:
        """Configs options that route the basin's first ``months`` of gridded runoff."""
        return GRID_NAMES | {
            'network_file': self.network_file,
            'forcing': 'grid',
            'runoff_files': self.runoff_files[:months],
            'grid_weights_file': self.weights_file,
        }


def assert_same(actual: np.ndarray, desired: np.ndarray, tolerance: float = TOLERANCE) -> None:
    """
    Require two answers that should be identical to agree within ``tolerance`` of each river's largest value: each row
    of a (river, time) array is measured against its own largest value, and each value of a per river vector against
    itself. A floor of a millionth of the largest value anywhere keeps values near zero from demanding exact equality.
    """
    actual, desired = np.asarray(actual, dtype=np.float64), np.asarray(desired, dtype=np.float64)
    assert actual.shape == desired.shape, f'shapes differ: {actual.shape} and {desired.shape}'
    if not desired.size:
        return
    magnitude = np.abs(desired)
    scale = magnitude.max(axis=-1, keepdims=True) if desired.ndim > 1 else magnitude
    difference = np.abs(actual - desired) / np.maximum(scale, magnitude.max() * 1e-6 or 1.0)
    worst = np.unravel_index(np.argmax(difference), difference.shape)
    assert difference[worst] <= tolerance, (
        f'{np.count_nonzero(difference > tolerance)} of {difference.size} values differ by more than {tolerance} of '
        f"their river's largest value; the worst, at {worst}, is {actual[worst]} against {desired[worst]}"
    )


@dataclass
class Routed:
    """What the router handed its discharge writer for each runoff file, and the channel state it ended with."""

    dates: list[np.ndarray]
    discharge: list[np.ndarray]  # (river, time) per runoff file
    discharge_files: list[str]
    final_state: np.ndarray
    router: rr.Router


def route(
    directory: Path,
    files: list[tuple[np.ndarray, np.ndarray]] | None = None,
    *,
    threads: int = 1,
    writer=None,
    **options,
) -> Routed:
    """
    Route and keep what the router hands its discharge writer. ``files`` are catchment runoff volumes shaped in
    memory, one (dates, (river, time) volumes) pair per runoff file with rows in network file order, written to
    catchment runoff files and routed as ``forcing`` catchment; without them ``options`` name the runoff to read. Any
    file the run writes goes to ``directory``. ``writer`` also writes the discharge when given.
    """
    if files is not None:
        river_ids = pd.read_parquet(options['network_file'], columns=['riverId'])['riverId'].to_numpy()
        no_areas = np.full(river_ids.shape[0], np.nan)  # volumes are read without their catchment areas
        runoff_files = [directory / f'catchment_runoff_{i}.nc' for i in range(len(files))]
        for (dates, volumes), runoff_file in zip(files, runoff_files, strict=True):
            rr.CatchmentRunoff.to_netcdf(runoff_file, dates, volumes, river_ids, no_areas, as_volumes=True)
        options |= {'forcing': 'catchment', 'runoff_files': runoff_files}
    n_outputs = len(options.get('runoff_files', [])) or 1
    defaults = {
        'discharge_files': [directory / f'discharge_{i}.zarr' for i in range(n_outputs)],
        'unstable_coefficients': 'ignore',
        'log': False,
        'progress_bar': False,
    }
    configs = rr.Configs(**(defaults | options))
    routed = Routed([], [], [], np.zeros(0), None)

    def keep_discharge(router, dates, discharge_array, discharge_file, runoff_file='', *, thread_pool=None, threads=1):
        routed.dates.append(dates.copy())
        routed.discharge.append(discharge_array.copy())
        routed.discharge_files.append(Path(discharge_file).name)
        if writer is not None:
            writer(
                router, dates, discharge_array, discharge_file, runoff_file, thread_pool=thread_pool, threads=threads
            )

    routed.router = rr.Router(configs).set_discharge_writer(keep_discharge)
    with ThreadPoolExecutor(threads) if threads > 1 else contextlib.nullcontext() as pool:
        routed.router.route(thread_pool=pool, threads=threads)
    routed.final_state = np.asarray(routed.router.channel_state).copy()
    return routed


def reference_runs(network_file: str, weights: str, runoff: list[str], catchment_runoff: list[str]) -> dict[str, dict]:
    """
    Each run's thread count and configs, with its input paths relative to the package. They are built in this order,
    since later runs read the catchment runoff files and the standard run's final state.
    """
    common = {'network_file': network_file, 'dt_routing': 3600, 'unstable_coefficients': 'ignore'}
    grid = common | GRID_NAMES | {'forcing': 'grid', 'runoff_files': runoff, 'grid_weights_file': weights}
    stabilized = grid | {'network_type': 'stabilized'}
    catchment = common | {'forcing': 'catchment', 'runoff_files': catchment_runoff}
    channel = common | {
        'forcing': 'channel',
        'channel_state_init_file': 'outputs/static_standard/final_state.parquet',
        'dt_total': 30 * 86400,
        'start_datetime': '2000-07-01',
    }
    return {
        'static_standard': {'threads': 1, 'configs': grid},
        'static_stabilized': {'threads': 1, 'configs': stabilized},
        'static_standard_4_threads': {'threads': 4, 'configs': grid},
        'static_stabilized_4_threads': {'threads': 4, 'configs': stabilized},
        'catchment_file': {'threads': 1, 'configs': catchment},
        'channel': {'threads': 1, 'configs': channel},
    }


def route_run(package: Path, run: dict, discharge_files: list[Path], final_state_file: Path, writer) -> None:
    """Route one run from the package's files, writing its discharge with ``writer`` and its final state to a file."""
    options = {
        key: str(package / value) if key in INPUT_PATH_OPTIONS else value for key, value in run['configs'].items()
    }
    if 'runoff_files' in options:
        options['runoff_files'] = [str(package / path) for path in options['runoff_files']]
    configs = rr.Configs(
        **options,
        discharge_files=discharge_files,
        channel_state_final_file=final_state_file,
        log=False,
        progress_bar=False,
    )
    router = rr.Router(configs).set_discharge_writer(writer)
    threads = run['threads']
    with ThreadPoolExecutor(threads) if threads > 1 else contextlib.nullcontext() as pool:
        router.route(thread_pool=pool, threads=threads)


def write_catchment_runoff(package: Path, catchment_runoff: dict, runoff_file: str, path: Path) -> None:
    """Aggregate one of the package's gridded runoff files into a catchment runoff file with ``aggregate_to_file``."""
    runoff = rr.GridRunoff(
        package / catchment_runoff['grid_weights_file'], as_volumes=catchment_runoff['as_volumes'], **GRID_NAMES
    )
    runoff.aggregate_to_file(package / runoff_file, path)


def build(region_dir: Path, outlet_river_id: int | None, runoff: list[Path], package: Path) -> None:
    """Copy or cut the inputs into ``package``, route every reference run, and write the manifest."""
    if package.exists():
        raise FileExistsError(f'{package} exists; the reference solution is only rebuilt after it is deleted')
    network_source = region_dir / 'network.parquet'
    weights_source = next(region_dir.glob('gridweights_*.nc'))
    network_file, weights = f'inputs/{network_source.name}', f'inputs/{weights_source.name}'
    (package / 'inputs').mkdir(parents=True)
    if outlet_river_id is None:
        shutil.copy2(network_source, package / network_file)
    else:
        streams.subset_network_to_river(outlet_river_id, network_source, package / network_file)
    river_ids = pd.read_parquet(package / network_file, columns=['riverId'])['riverId'].to_numpy()
    with xr.open_dataset(weights_source) as region_weights:
        in_package = np.isin(region_weights[REGION_VAR_RIVER_ID].to_numpy(), river_ids)
        region_weights.isel(index=in_package).rename({REGION_VAR_RIVER_ID: 'riverId'}).to_netcdf(package / weights)
    for source in runoff:
        copy = shutil.copytree if source.is_dir() else shutil.copy2
        copy(source, package / 'inputs' / source.name)
    runoff_copies = [f'inputs/{source.name}' for source in runoff]

    (package / 'outputs').mkdir()
    catchment_runoff = {
        'grid_weights_file': weights,
        'as_volumes': True,
        'runoff_files': runoff_copies,
        'files': [f'outputs/catchment_runoff_{Path(path).stem}.nc' for path in runoff_copies],
    }
    for runoff_file, path in zip(catchment_runoff['runoff_files'], catchment_runoff['files'], strict=True):
        write_catchment_runoff(package, catchment_runoff, runoff_file, package / path)

    writers.ZARR_KEEPBITS = 23  # every float32 mantissa bit, so zarr_writer writes the discharge unrounded
    runs = reference_runs(network_file, weights, runoff_copies, catchment_runoff['files'])
    for name, run in runs.items():
        (package / 'outputs' / name).mkdir()
        months = [Path(path).stem for path in run['configs'].get('runoff_files', [])]
        run['discharge'] = [f'outputs/{name}/discharge_{month}.zarr' for month in months] or [
            f'outputs/{name}/discharge.zarr'
        ]
        run['final_state'] = f'outputs/{name}/final_state.parquet'
        discharge_files = [package / path for path in run['discharge']]
        route_run(package, run, discharge_files, package / run['final_state'], writers.zarr_writer)

    git = ['git', '-C', str(Path(__file__).parent)]
    manifest = {
        'description': f'region {region_dir.name.removeprefix("region=")}, '
        f'{"whole" if outlet_river_id is None else f"cut to river {outlet_river_id} and its upstreams"}, '
        f'routed with {", ".join(source.name for source in runoff)} in sequence',
        'sources': {
            'region_dir': str(region_dir),
            'outlet_river_id': outlet_river_id,
            'runoff': [str(source) for source in runoff],
        },
        'built': {
            'date': datetime.now().isoformat(timespec='seconds'),
            'river_route': rr.__version__,
            'git_commit': subprocess.run([*git, 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip(),
            'git_uncommitted_changes': bool(
                subprocess.run([*git, 'status', '--porcelain'], capture_output=True).stdout
            ),
            'python': platform.python_version(),
            'numpy': np.__version__,
            'numba': numba.__version__,
            'xarray': xr.__version__,
            'zarr': zarr.__version__,
            'machine': platform.machine(),
            'platform': platform.platform(),
        },
        'catchment_runoff': catchment_runoff,
        'runs': runs,
    }
    (package / 'manifest.json').write_text(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build the reference solution the test suite reproduces.')
    parser.add_argument('--region-dir', type=Path, required=True, help='folder with network.parquet and grid weights')
    parser.add_argument('--outlet-river-id', type=int, help='cut the region to this river and its upstreams')
    parser.add_argument('--runoff', type=Path, nargs='+', required=True, help='gridded runoff files, routed in order')
    parser.add_argument('--package', type=Path, default=PACKAGE, help='where to write the reference solution')
    args = parser.parse_args()
    runoff = [path.expanduser() for path in args.runoff]
    build(args.region_dir.expanduser(), args.outlet_river_id, runoff, args.package.expanduser())
