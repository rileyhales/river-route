"""The frozen Configs object that holds and validates every option of a routing run and of preparing its runoff."""

import json
import os
import types
from collections.abc import Mapping
from dataclasses import dataclass, field, fields
from typing import Any, ClassVar, Literal, Self, TypeAliasType, get_args, get_origin, get_type_hints

import numpy as np
import pandas as pd
import xarray as xr

from river_route._logging import build_logger
from river_route.types import PathInput, PathList

__all__ = ['Configs', 'is_dev_null']

_PATH_INPUT_TYPES: frozenset[type] = frozenset(get_args(PathInput))
_DEV_NULL: frozenset[str] = frozenset({os.devnull, '/dev/null'})


def is_dev_null(path: PathInput) -> bool:
    """True for the null device, which always counts as a path that exists and discards whatever is written."""
    return str(path) in _DEV_NULL


def _derive_valid_values(cls: type) -> dict[str, frozenset[str]]:
    """The allowed values of every field annotated with a Literal, by field name."""
    result = {}
    for name, hint in get_type_hints(cls).items():
        if not name.startswith('_') and get_origin(hint) is Literal:
            result[name] = frozenset(get_args(hint))
    return result


def _unalias(hint: Any) -> Any:
    """The type a ``type X = ...`` alias names, such as list[PathInput] for PathList, whose origin is otherwise None."""
    return hint.__value__ if isinstance(hint, TypeAliasType) else hint


def _derive_path_sets(cls: type) -> tuple[frozenset[str], frozenset[str]]:
    """The names of the fields annotated PathInput or PathInput | None, and of the fields annotated PathList."""
    single, lists = set(), set()
    for name, hint in get_type_hints(cls).items():
        if name.startswith('_'):
            continue
        hint = _unalias(hint)
        if get_origin(hint) is list:
            lists.add(name)  # PathList
        elif get_origin(hint) is types.UnionType and set(get_args(hint)) >= _PATH_INPUT_TYPES:
            single.add(name)  # PathInput or PathInput | None
    return frozenset(single), frozenset(lists)


@dataclass(kw_only=True, frozen=True)
class Configs:
    """
    Accepts and validates every option of a routing run and of preparing its runoff.
    """

    # annotate file path fields with PathInput or PathList
    # _derive_path_sets() will detect them by inspecting class annotations

    # Routing procedure selectors — the Router chooses the routing method and the runoff reader from these
    coefficients: Literal['static', 'dynamic'] = 'static'
    forcing: Literal['channel', 'catchment', 'grid', 'ecmwf_grib'] = 'channel'
    transform: Literal['uniform'] = 'uniform'
    network_type: Literal['standard', 'stabilized'] = 'standard'  # route as given, or add substeps/subcycles
    unstable_coefficients: Literal['warn', 'raise', 'ignore'] = 'warn'  # action when a river is not stable for dt

    # Network and routing descriptor
    network_file: PathInput | None = None

    # Core Routing Files
    discharge_dir: PathInput | None = None
    discharge_files: PathList = field(default_factory=list)  # optional override for explicit output paths
    channel_state_init_file: PathInput | None = None
    channel_state_final_file: PathInput | None = None

    # Types of runoff data handling
    grid_accumulation_type: Literal['incremental', 'cumulative'] = 'incremental'
    runoff_processing_mode: Literal['sequential', 'ensemble'] = 'sequential'
    runoff_depth_unit: str | None = None  # unit of the runoff depths; None reads the file attributes, else meters
    force_positive_runoff: bool = False  # clip negative runoff depths to zero
    force_uniform_timesteps: bool = True  # resample runoff with irregular timesteps to the first timestep
    as_volumes: bool = False  # prepare volumes (m³) instead of depths (m); routing always uses volumes

    # Runoff sources: catchment runoff files, or gridded runoff aggregated to catchments with a weight table
    runoff_files: PathList = field(default_factory=list)  # read as the forcing names
    grid_weights_file: PathInput | None = None  # required for the grid runoff types
    var_river_id: str = 'river_id'
    var_discharge: str = 'Q'
    var_grid_runoff: str = 'ro'
    var_x: str = 'x'  # grid x dimension
    var_y: str = 'y'  # grid y dimension
    var_cell: str = 'cell'  # ecmwf_grib cell dimension
    var_t: str = 'time'

    # Time options
    dt_routing: int = 0  # Interval in seconds between calculating discharges, <= dt_runoff
    dt_runoff: int = 0  # Interval in seconds between forcing values, >= dt_routing
    dt_discharge: int = 0  # Interval in seconds between discharge outputs, >= dt_runoff
    dt_total: int = 0  # Length in seconds of the simulation, >= dt_discharge
    start_datetime: str = '1970-01-01'

    # Misc behavior that users may want to override
    log: bool = True
    progress_bar: bool = True
    log_level: Literal['DEBUG', 'INFO', 'PROGRESS', 'WARNING', 'ERROR', 'CRITICAL'] = 'PROGRESS'
    log_stream: str = 'stdout'
    log_format: str = '%(levelname)s - %(asctime)s - %(message)s'

    # False until validate_routing, or validate_runoff, passes; neither satisfies the other
    _routing_validated: bool = field(default=False, init=False, repr=False, compare=False)
    _runoff_validated: bool = field(default=False, init=False, repr=False, compare=False)

    # the single path fields that name outputs, so their directory needs to exist rather than the file
    _OUTPUT_FILES: ClassVar[frozenset[str]] = frozenset({'channel_state_final_file'})
    # 2 options for specifying how the computed discharge files are saved
    _OUTPUT_DIRS: ClassVar[frozenset[str]] = frozenset({'discharge_dir'})
    _OUTPUT_FILE_LISTS: ClassVar[frozenset[str]] = frozenset({'discharge_files'})
    # extension of the outputs named in discharge_dir, the store written by zarr_writer, the default writer
    _DISCHARGE_SUFFIX: ClassVar[str] = '.zarr'
    # the network types the routing method of each coefficients option can route
    _NETWORK_TYPES_FOR_COEFFICIENTS: ClassVar[dict[str, frozenset[str]]] = {
        'static': frozenset({'standard', 'stabilized'}),
        'dynamic': frozenset({'standard'}),
    }

    # Populated at module level below
    _SINGLE_PATH_FIELDS: ClassVar[frozenset[str]]
    _LIST_PATH_FIELDS: ClassVar[frozenset[str]]
    _VALID_VALUES: ClassVar[dict[str, frozenset[str]]]

    def __post_init__(self) -> None:
        for name, allowed in self._VALID_VALUES.items():
            value = getattr(self, name)
            if value not in allowed:
                raise ValueError(f'{name} must be one of {sorted(allowed)}, got {value!r}')
        # turn off progress bar if logging was turned off but progress was left at default on
        object.__setattr__(self, 'progress_bar', bool(self.log) and bool(self.progress_bar))
        self._coerce_path_list_fields()
        self._absolutize_paths()
        self._resolve_discharge_dir()
        return

    @classmethod
    def from_json(cls, path: PathInput) -> Self:
        """
        Build Configs from a JSON file. Relative paths in it are made absolute against the working directory.

        Args:
            path: JSON file holding an object of config keys and values, as ``to_json`` writes

        Raises:
            ValueError: if the file has a key that is not a config option
        """
        with open(path, encoding='utf-8') as f:
            raw = json.load(f)
        cls._check_keys(raw)
        return cls(**raw)

    def to_dict(self) -> dict[str, Any]:
        """Every option by name, with discharge_files left empty when they were named from discharge_dir."""
        values = {f.name: getattr(self, f.name) for f in fields(self) if f.init}
        values = {key: list(value) if isinstance(value, list) else value for key, value in values.items()}
        if self.discharge_dir:
            values['discharge_files'] = []
        return values

    def to_json(self, path: PathInput) -> None:
        """
        Write every option to a JSON file that ``from_json`` reads back.

        Args:
            path: JSON file to write
        """
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def _check_keys(cls, raw: Mapping[str, Any]) -> None:
        known = {f.name for f in fields(cls) if f.init}
        unknown = sorted(set(raw) - known)
        if unknown:
            raise ValueError(f'Unrecognized config key(s): {", ".join(map(repr, unknown))}')
        return

    def _coerce_path_list_fields(self) -> None:
        """Normalize any list-of-paths field given as a single string to [str]."""
        for key in self._LIST_PATH_FIELDS:
            val = getattr(self, key)
            if isinstance(val, PathInput) and val:
                object.__setattr__(self, key, [str(val)])
        return

    def _absolutize_paths(self) -> None:
        """Convert all path fields to absolute path strings in-place. The null device is kept as it is named."""
        for key in self._SINGLE_PATH_FIELDS:
            val = getattr(self, key, None)
            if val:
                object.__setattr__(self, key, str(val) if is_dev_null(val) else os.path.abspath(val))
        for key in self._LIST_PATH_FIELDS:
            val = getattr(self, key, [])
            if val:
                object.__setattr__(self, key, [str(p) if is_dev_null(p) else os.path.abspath(p) for p in val])
        return

    def _verify_input_files_exist(self) -> None:
        """Raise FileNotFoundError for any set input path that does not exist."""
        input_single = self._SINGLE_PATH_FIELDS - self._OUTPUT_FILES - self._OUTPUT_DIRS
        input_list = self._LIST_PATH_FIELDS - self._OUTPUT_FILE_LISTS
        for key in input_single:
            val = getattr(self, key, None)
            if val and not os.path.exists(val):
                raise FileNotFoundError(f'{key} not found: {val}')
        for key in input_list:
            for path in getattr(self, key, []):
                if not os.path.exists(path):
                    raise FileNotFoundError(f'{key}: {path} not found')
        return

    def _resolve_discharge_dir(self) -> None:
        """Populate discharge_files from discharge_dir when explicit paths are not given."""
        if not self.discharge_dir:
            return
        if self.discharge_files:
            raise ValueError('Provide discharge_dir or discharge_files, not both')

        d = self.discharge_dir
        if is_dev_null(d):
            # every output is discarded, so one null device per input rather than names under a directory
            object.__setattr__(self, 'discharge_files', [d] * (len(self.runoff_files) or 1))
            return

        input_files = self.runoff_files
        if input_files:
            # each output is named for its input but takes the writer's extension, not the extension of the input
            stems = [os.path.splitext(os.path.basename(f))[0] for f in input_files]
            duplicates = sorted({stem for stem in stems if stems.count(stem) > 1})
            if duplicates:
                raise ValueError(
                    f'Input files with duplicate names would resolve to the same output file in discharge_dir: '
                    f'{", ".join(duplicates)}. Use discharge_files to give explicit output paths.'
                )
            discharge_files = [os.path.join(d, f'discharge_{stem}{self._DISCHARGE_SUFFIX}') for stem in stems]
        else:
            # channel routing has no runoff files to name its one output after
            discharge_files = [os.path.join(d, f'discharge{self._DISCHARGE_SUFFIX}')]
        object.__setattr__(self, 'discharge_files', discharge_files)
        return

    def _verify_output_directories_exist(self) -> None:
        """Raise NotADirectoryError for any output path whose parent directory does not exist."""
        paths: list[str] = []
        for key in self._OUTPUT_FILES:
            val = getattr(self, key, None)
            if val:
                paths.append(val)
        for key in self._OUTPUT_FILE_LISTS:
            paths.extend(getattr(self, key, []))
        for path in paths:
            if is_dev_null(path):  # devnull is a special case that is allowed to take outputs
                continue
            d = os.path.dirname(path)
            if not os.path.exists(d):
                raise NotADirectoryError(f'Directory not found for specified output: {path}')
        for key in self._OUTPUT_DIRS:
            val = getattr(self, key, None)
            if val and not is_dev_null(val) and not os.path.isdir(val):
                raise NotADirectoryError(f'Output directory not found: {val}')
        return

    def validate_routing(self) -> Self:
        """
        Validate the options for routing: the options the chosen procedure requires are set and consistent, input
        paths exist, and outputs are in existing directories. Called by Router.route. Returns immediately once the
        Configs has passed it. The contents of the input files are not read; call deep_validate for that.

        Raises:
            ValueError, FileNotFoundError, NotADirectoryError: if any option is missing or invalid
            NotImplementedError: if no routing method routes the chosen options yet
        """
        if not self.network_file:
            raise ValueError('network_file is required to route')
        if self._routing_validated:
            return self
        self._verify_input_files_exist()
        self._verify_output_directories_exist()
        # the Router picks the null writer for the whole run, so a mix of discarded and written outputs has no meaning
        null_outputs = [is_dev_null(f) for f in self.discharge_files]
        if any(null_outputs) and not all(null_outputs):
            raise ValueError('discharge_files mixes the null device with real paths; use one or the other')
        if self.forcing == 'channel':
            if not self.discharge_files:
                raise ValueError('Provide discharge_dir (or discharge_files for explicit output paths)')
            for key in ('channel_state_init_file', 'dt_routing', 'dt_total'):
                if not getattr(self, key, None):
                    raise ValueError(f'{key} is required for channel routing')
            runoff_source = [key for key in ('runoff_files', 'grid_weights_file') if getattr(self, key)]
            if runoff_source:
                build_logger(self, 'configs').warning(
                    f'forcing is channel, so {", ".join(runoff_source)} will be ignored and no runoff is routed'
                )
            if len(self.discharge_files) != 1:
                raise ValueError('Channel routing requires exactly one entry in discharge_files')
        else:
            if not self.discharge_files:
                raise ValueError('Provide discharge_dir (or discharge_files for explicit output paths)')
            if not self.runoff_files:
                raise ValueError('runoff_files is required for runoff forcing')
            if self.forcing == 'catchment' and self.grid_weights_file:
                raise ValueError('grid_weights_file is not used with forcing catchment')
            if self.forcing != 'catchment' and not self.grid_weights_file:
                raise ValueError(f'grid_weights_file is required with forcing {self.forcing}')
            if len(self.discharge_files) != len(self.runoff_files):
                raise ValueError('Number of resolved discharge output files must match number of input files')
            outputs = [f for f in self.discharge_files if not is_dev_null(f)]
            if len(set(outputs)) != len(outputs):
                raise ValueError('discharge_files contains duplicate paths; each input file needs a distinct output')
        # options that are consistent, but that no routing method routes yet
        if self.network_type not in self._NETWORK_TYPES_FOR_COEFFICIENTS[self.coefficients]:
            raise NotImplementedError(
                f'{self.coefficients} coefficients cannot route a {self.network_type} network yet'
            )
        object.__setattr__(self, '_routing_validated', True)
        return self

    def validate_runoff(self) -> Self:
        """
        Validate the options for preparing gridded runoff: grid_weights_file is set and every input path exists.
        Called by the grid runoff classes. Returns immediately once the Configs has passed it, which does not validate
        it for routing. The contents of the input files are not read; call deep_validate for that.

        Raises:
            ValueError, FileNotFoundError: if any option is missing or invalid
        """
        if not self.grid_weights_file:
            raise ValueError('grid_weights_file is required to prepare runoff')
        if self._runoff_validated:
            return self
        self._verify_input_files_exist()
        object.__setattr__(self, '_runoff_validated', True)
        return self

    def deep_validate(self) -> Self:
        """
        Validate the contents of the network file, the grid weight table, and the initial channel state, whichever are
        set, and their consistency with each other: the network file columns, types, and value ranges and that it is
        topologically sorted, that the grid weight table matches the network file and its proportions sum to 1 per
        river, and that the initial channel state has one row per river, or on a stabilized network one per
        sub-reach. The runoff files are not read.

        Nothing calls this for you. It reads every input file, which is the same work routing is about to do, so
        run it once on inputs you have not checked before rather than on every route. ``validate_routing`` and
        ``validate_runoff`` check the options and the paths only.

        Raises:
            ValueError: if any file or combination of files is invalid
        """
        network_df = self._deep_validate_network_file() if self.network_file else None
        if self.grid_weights_file:
            self._deep_validate_grid_weights_file(self.grid_weights_file, network_df)
        if self.channel_state_init_file:
            self._deep_validate_channel_state_init_file(network_df)
        return self

    def _deep_validate_network_file(self) -> pd.DataFrame:
        # network file parquet has columns: riverId, nextRiverId, muskingumK, muskingumX, riverIndex, upstreamCount
        # riverId should be non-null, positive, integer, and unique except for -1 (outlets/sinks)
        # nextRiverId should be non-null, integer, all -1 or positive, and exist in riverId (except for -1)
        # muskingumK should be positive float
        # muskingumX should be positive float less than or equal to 0.5
        try:
            network_df = pd.read_parquet(self.network_file)
        except Exception as e:
            raise ValueError('Error reading network file. Must be valid parquet file') from e
        for column in ('riverId', 'nextRiverId', 'muskingumK', 'muskingumX'):
            if column not in network_df.columns:
                raise ValueError(f'{self.network_file} missing {column} column')
        if network_df['riverId'].isna().any():
            raise ValueError(f'{self.network_file} riverId column contains null values')
        if not pd.api.types.is_integer_dtype(network_df['riverId']):
            raise ValueError(f'{self.network_file} riverId column must be integer type')
        if not network_df['riverId'].is_unique:
            raise ValueError(f'{self.network_file} riverId column must be unique')
        if np.any(network_df['riverId'] == -1):
            raise ValueError(f'{self.network_file} riverId must not be -1, which marks a basin outlet in nextRiverId')
        if network_df['nextRiverId'].isna().any():
            raise ValueError(f'{self.network_file} nextRiverId column contains null values')
        if not pd.api.types.is_integer_dtype(network_df['nextRiverId']):
            raise ValueError(f'{self.network_file} nextRiverId column must be integer type')
        downstream_ids = set(network_df['nextRiverId'].unique())
        river_ids = set(network_df['riverId'].unique())
        if not downstream_ids.issubset(river_ids.union({-1})):
            raise ValueError(f'{self.network_file} nextRiverId values must exist in riverId (except -1)')
        if network_df[['muskingumK', 'muskingumX']].isna().any(axis=None):
            raise ValueError(f'{self.network_file} muskingumK and muskingumX columns must not contain null values')
        if np.any(network_df['muskingumK'] <= 0):
            raise ValueError(f'{self.network_file} muskingumK column must be positive')
        if np.any(network_df['muskingumX'] < 0) or np.any(network_df['muskingumX'] > 0.5):
            raise ValueError(f'{self.network_file} muskingumX column must be in the range [0, 0.5]')

        # dynamic coefficients are rebuilt in the kernel from K = dynamicAlpha * Q ** dynamicBeta
        if self.coefficients == 'dynamic':
            for column in ('dynamicAlpha', 'dynamicBeta'):
                if column not in network_df.columns:
                    raise ValueError(
                        f'{self.network_file} missing {column} column required when coefficients is dynamic'
                    )
                if network_df[column].isna().any():
                    raise ValueError(f'{self.network_file} {column} column contains null values')
                if not pd.api.types.is_numeric_dtype(network_df[column]):
                    raise ValueError(f'{self.network_file} {column} column must be numeric type')
            if np.any(network_df['dynamicAlpha'] <= 0):
                raise ValueError(f'{self.network_file} dynamicAlpha column must be strictly positive')

        # check topological sort: every nextRiverId must appear later in the table than its upstream
        river_id_index = {int(river_id): i for i, river_id in enumerate(network_df['riverId'])}
        for upstream_idx, ds_id in enumerate(network_df['nextRiverId']):
            if int(ds_id) == -1:
                continue
            if river_id_index[int(ds_id)] <= upstream_idx:
                raise ValueError(f'{self.network_file} is not topologically sorted (upstream to downstream)')
        self._deep_validate_dfs_order_columns(network_df)
        return network_df

    def _deep_validate_dfs_order_columns(self, network_df: pd.DataFrame) -> None:
        # riverIndex numbers the rows one apart and upstreamCount counts the rivers upstream of each river, so in DFS
        # order a river's upstream watershed is exactly the rows from its riverIndex - upstreamCount to its riverIndex.
        # That holds when each count is the sum over the river's upstreams of their counts plus one, and each
        # upstream's own range of rows lies inside the range of the river it drains into.
        for column in ('riverIndex', 'upstreamCount'):
            if column not in network_df.columns:
                raise ValueError(f'{self.network_file} missing {column} column')
            if not pd.api.types.is_integer_dtype(network_df[column]):
                raise ValueError(f'{self.network_file} {column} column must be integer type')
        if np.any(np.diff(network_df['riverIndex'].to_numpy(dtype=np.int64)) != 1):
            raise ValueError(f'{self.network_file} riverIndex must increase by one from each row to the next')
        count = network_df['upstreamCount'].to_numpy(dtype=np.int64)
        downstream = pd.Index(network_df['riverId']).get_indexer(network_df['nextRiverId'])  # -1 at outlets
        upstream = np.flatnonzero(downstream >= 0)
        down = downstream[upstream]
        if not np.array_equal(count, np.bincount(down, weights=count[upstream] + 1, minlength=count.shape[0])):
            raise ValueError(f'{self.network_file} upstreamCount must count every river upstream of each river')
        if np.any(upstream - count[upstream] < down - count[down]):
            raise ValueError(
                f'{self.network_file} is not in DFS order: the rivers upstream of each river must be the rows '
                f'immediately before it'
            )
        return

    def _deep_validate_grid_weights_file(self, grid_weights_file: PathInput, network_df: pd.DataFrame | None) -> None:
        # weights should be netcdf with variables river_id, the cell index columns, x, y, area_sqm, proportion.
        # A reduced grid locates its cells with cell_index, a grid with x_index and y_index.
        rid = self.var_river_id
        try:
            ds = xr.load_dataset(grid_weights_file)
        except Exception as e:
            raise ValueError('Error reading grid weights file. Must be valid netCDF file') from e
        cell_columns = ('cell_index',) if self.forcing == 'ecmwf_grib' else ('x_index', 'y_index')
        expected_variables = (rid, *cell_columns, 'x', 'y', 'area_sqm', 'proportion')
        for variable in expected_variables:
            if variable not in ds:
                raise ValueError(f'Grid weights file missing {variable} variable')
        if np.any(ds[rid].isnull()):
            raise ValueError(f'Grid weights {rid} variable contains null values')
        if not pd.api.types.is_integer_dtype(ds[rid].dtype):
            raise ValueError(f'Grid weights {rid} variable must be integer type')
        if network_df is not None and not np.isin(ds[rid].values, network_df['riverId'].to_numpy()).all():
            raise ValueError(f'Grid weights {rid} values must exist in the riverId column of the network file')
        for variable in expected_variables[1:]:
            if np.any(ds[variable].isnull()):
                raise ValueError(f'Grid weights {variable} variable contains null values')
            if not pd.api.types.is_numeric_dtype(ds[variable].dtype):
                raise ValueError(f'Grid weights {variable} variable must be numeric type')
        if np.any(ds['area_sqm'] <= 0):
            raise ValueError('Grid weights area_sqm variable must be positive')
        if np.any(ds['proportion'] <= 0) or np.any(ds['proportion'] > 1):
            raise ValueError('Grid weights proportion variable must be in the range (0, 1]')
        # pandas groups in one hashed pass. xarray's groupby materializes a DataArray per group and concatenates
        # them, which on a weight table with a group per river costs more than the whole routing it validates.
        proportions_sum = pd.Series(ds['proportion'].values).groupby(ds[rid].values).sum()
        if not np.allclose(proportions_sum.to_numpy(), 1.0):
            raise ValueError('Grid weights proportion variable must sum to 1 for each river_id')
        return

    def _deep_validate_channel_state_init_file(self, network_df: pd.DataFrame | None) -> None:
        # initial channel state should be parquet with a column named Q, non-null, numeric, and non-negative.
        # it has one row per river of the network file in the same order, or on a stabilized network one row per
        # sub-reach, at least one per river, whose count depends on dt_routing and is checked when routing.
        try:
            state_df = pd.read_parquet(self.channel_state_init_file)
        except Exception as e:
            raise ValueError('Error reading initial state file. Must be valid parquet file') from e
        if 'Q' not in state_df.columns:
            raise ValueError('Initial state file missing Q column')
        if state_df['Q'].isna().any():
            raise ValueError('Initial state file Q column contains null values')
        if not pd.api.types.is_numeric_dtype(state_df['Q']):
            raise ValueError('Initial state file Q column must be numeric type')
        if np.any(state_df['Q'] < 0):
            raise ValueError('Initial state file Q column must be non-negative')
        if network_df is None:
            return
        if self.network_type == 'stabilized' and state_df.shape[0] < network_df.shape[0]:
            raise ValueError(
                f'Initial state file must have a row per sub-reach, at least one per river of {self.network_file}'
            )
        if self.network_type == 'standard' and state_df.shape[0] != network_df.shape[0]:
            raise ValueError(f'Initial state file must have the same number of rows as {self.network_file}')
        return


Configs._SINGLE_PATH_FIELDS, Configs._LIST_PATH_FIELDS = _derive_path_sets(Configs)
Configs._VALID_VALUES = _derive_valid_values(Configs)
