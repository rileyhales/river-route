import difflib
import json
import os
import types
from collections.abc import Mapping
from dataclasses import dataclass, field, fields
from typing import Any, ClassVar, Literal, Self, get_args, get_origin, get_type_hints

import numpy as np
import pandas as pd
import xarray as xr
import yaml

from river_route.types import PathInput, PathList

__all__ = ['Configs']

_PATH_INPUT_TYPES: frozenset[type] = frozenset(get_args(PathInput))


@dataclass(kw_only=True, frozen=True)
class Configs:
    """
    Frozen configuration options. Build one from keyword arguments with ``Configs(...)`` or from a file with
    ``Configs.from_file``, then pass it to ``Router`` or ``RunoffGaussianGrid``. Building a Configs normalizes it: file
    paths are made absolute, discharge_files are derived from discharge_dir, and selector values are checked. The
    options are validated when they are used: ``Router.route`` calls ``validate_routing`` and ``RunoffGaussianGrid``
    calls ``validate_runoff``.
    Either one checks that the required options are set, that input paths exist, and that output paths are in
    existing directories, and then marks the Configs as validated so it is not checked again.

    ``deep_validate`` reads the input files and checks their contents. Nothing calls it for you: it is the one
    validation that costs as much as the read it repeats, so it is a method to run once on inputs you have not
    checked before, not something a route pays for every time.

    A Configs is set once and never changed. There is no method to copy one with an option altered: build the
    Configs you want. ``to_json`` and ``to_yaml`` write the options to a file that ``from_file`` reads back.
    """

    # annotate file path fields with PathInput or PathList
    # _derive_path_sets() will detect them by inspecting class annotations

    # Routing procedure selectors — describe the procedure resolved to a kernel by router._kernel_registry
    coeff: Literal['static', 'dynamic'] = 'static'
    forcing: Literal['channel', 'vlateral'] = 'channel'
    transform: Literal['uniform', 'unit_hydrograph'] = 'uniform'
    network: Literal['standard', 'expanded'] = 'standard'
    routing_order: Literal['time', 'river'] = 'time'  # sweep every river per step, or each river's whole series

    # Core Routing Files
    params_file: PathInput | None = None
    discharge_dir: PathInput | None = None
    discharge_files: PathList = field(default_factory=list)  # optional override for explicit output paths
    channel_state_init_file: PathInput | None = None
    channel_state_final_file: PathInput | None = None

    # Time options
    dt_routing: int = 0
    dt_total: int = 0
    dt_discharge: int = 0
    dt_runoff: int = 0
    start_datetime: str = '1970-01-01'

    # For vlateral / runoff transformation - used by TransformMuskingum subclasses
    vlateral_files: PathList = field(default_factory=list)
    grid_runoff_files: PathList | None = field(default_factory=list)
    grid_weights_file: PathInput | None = None
    grid_accumulation_type: Literal['incremental', 'cumulative'] = 'incremental'
    runoff_processing_mode: Literal['sequential', 'ensemble'] = 'sequential'
    uh_kernel_file: PathInput | None = None
    uh_state_init_file: PathInput | None = None
    uh_state_final_file: PathInput | None = None

    # Gridded runoff preparation (RunoffGaussianGrid)
    runoff_depth_unit: str | None = None  # unit of the runoff depths; None reads the file attributes, else meters
    force_positive_runoff: bool = False  # clip negative runoff depths to zero
    force_uniform_timesteps: bool = True  # resample runoff with irregular timesteps to the first timestep
    as_volumes: bool = False  # prepare volumes (m³) instead of depths (m); routing always uses volumes

    # Validation behavior
    unstable_coefficients: Literal['warn', 'raise', 'ignore'] = 'warn'  # action when a river is not stable for dt

    # Misc behavior that users may want to override
    log: bool = True
    progress_bar: bool = True
    log_level: Literal['DEBUG', 'INFO', 'PROGRESS', 'WARNING', 'ERROR', 'CRITICAL'] = 'PROGRESS'
    log_stream: str = 'stdout'
    log_format: str = '%(levelname)s - %(asctime)s - %(message)s'
    var_river_id: str = 'river_id'
    var_discharge: str = 'Q'
    var_grid_runoff: str = 'ro'
    var_x: str = 'x'
    var_y: str = 'y'
    var_t: str = 'time'

    # False until validate_routing or validate_runoff passes
    _validated: bool = field(default=False, init=False, repr=False, compare=False)

    # special subset of auto-detected PathLists where the directory needs to exist, not the file
    _OUTPUT_FILES: ClassVar[frozenset[str]] = frozenset({'channel_state_final_file', 'uh_state_final_file'})
    # 2 options for specifying how the computed discharge files are saved
    _OUTPUT_DIRS: ClassVar[frozenset[str]] = frozenset({'discharge_dir'})
    _OUTPUT_FILE_LISTS: ClassVar[frozenset[str]] = frozenset({'discharge_files'})
    # where the runoff the router routes comes from, which only the runoff classes' readers take from these options

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
        # normalize paths given as strings to lists. make paths absolute.
        self._coerce_path_list_fields()
        self._absolutize_paths()
        # derive discharge_files from discharge_dir + input files when not explicitly provided
        self._resolve_discharge_dir()
        return

    # --- construction, copying, and serialization ---
    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> Self:
        """
        Build Configs from a mapping, reporting unrecognized keys by name instead of raising the
        dataclass TypeError. Suggests the closest valid key name for likely typos.

        Args:
            raw: mapping of config key to value, e.g. parsed from a YAML or JSON config file

        Raises:
            ValueError: if the mapping contains any key that is not a config option
        """
        cls._check_keys(raw)
        return cls(**raw)

    @classmethod
    def from_file(cls, path: PathInput) -> Self:
        """Build Configs from a .json, .yml, or .yaml file."""
        path = str(path)
        if path.endswith('.json'):
            return cls.from_json(path)
        if path.endswith(('.yml', '.yaml')):
            return cls.from_yaml(path)
        raise ValueError(f'Unrecognized config file type: {path}. Must be .json, .yml, or .yaml')

    @classmethod
    def from_json(cls, path: PathInput) -> Self:
        """Build Configs from a JSON file."""
        with open(path) as f:
            return cls.from_mapping(json.load(f))

    @classmethod
    def from_yaml(cls, path: PathInput) -> Self:
        """Build Configs from a YAML file."""
        with open(path) as f:
            return cls.from_mapping(yaml.safe_load(f) or {})

    def to_dict(self) -> dict[str, Any]:
        """
        Every option as a mapping that ``from_mapping`` builds the same Configs from. Discharge files derived from
        ``discharge_dir`` are left out so that the directory and the files are not both set.
        """
        values = {f.name: getattr(self, f.name) for f in fields(self) if f.init}
        values = {key: list(value) if isinstance(value, list) else value for key, value in values.items()}
        if self.discharge_dir:
            values['discharge_files'] = []
        return values

    def to_json(self, path: PathInput) -> None:
        """Write every option to a JSON file that ``from_file`` reads back."""
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    def to_yaml(self, path: PathInput) -> None:
        """Write every option to a YAML file that ``from_file`` reads back."""
        with open(path, 'w') as f:
            yaml.safe_dump(self.to_dict(), f, sort_keys=False)

    @classmethod
    def _check_keys(cls, raw: Mapping[str, Any]) -> None:
        known = {f.name for f in fields(cls) if f.init}
        unknown = sorted(set(raw) - known)
        if unknown:
            described = []
            for key in unknown:
                close = difflib.get_close_matches(key, sorted(known), n=1)
                described.append(f'{key!r}' + (f' (did you mean {close[0]!r}?)' if close else ''))
            raise ValueError(f'Unrecognized config key(s): {", ".join(described)}')
        return

    # --- path normalization and verification ---
    def _coerce_path_list_fields(self) -> None:
        """Normalize any list-of-paths field given as a single string to [str]."""
        for key in self._LIST_PATH_FIELDS:
            val = getattr(self, key)
            if isinstance(val, PathInput) and val:
                object.__setattr__(self, key, [str(val)])
        return

    def _absolutize_paths(self) -> None:
        """Convert all relative path fields to absolute paths in-place."""
        for key in self._SINGLE_PATH_FIELDS:
            val = getattr(self, key, None)
            if val:
                object.__setattr__(self, key, os.path.abspath(val))
        for key in self._LIST_PATH_FIELDS:
            val = getattr(self, key, [])
            if val:
                object.__setattr__(self, key, [os.path.abspath(p) for p in val])
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
        input_files = self.vlateral_files or self.grid_runoff_files or []
        if input_files:
            basenames = [os.path.basename(f) for f in input_files]
            duplicates = sorted({name for name in basenames if basenames.count(name) > 1})
            if duplicates:
                raise ValueError(
                    f'Input files with duplicate names would resolve to the same output file in discharge_dir: '
                    f'{", ".join(duplicates)}. Use discharge_files to give explicit output paths.'
                )
            discharge_files = [os.path.join(d, f'discharge_{name}') for name in basenames]
        else:
            # Muskingum (no lateral inflow files)
            discharge_files = [os.path.join(d, 'discharge.nc')]
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
            d = os.path.dirname(path)
            if not os.path.exists(d):
                raise NotADirectoryError(f'Directory not found for specified output: {path}')
        for key in self._OUTPUT_DIRS:
            val = getattr(self, key, None)
            if val and not os.path.isdir(val):
                raise NotADirectoryError(f'Output directory not found: {val}')
        return

    # --- validation ---
    def validate_routing(self) -> Self:
        """
        Validate the options for routing: the options the chosen procedure requires are set and consistent, input
        paths exist, and outputs are in existing directories. Called by Router.route. Returns immediately once the
        Configs has been validated. The contents of the input files are not read; call deep_validate for that.

        Raises:
            ValueError, FileNotFoundError, NotADirectoryError: if any option is missing or invalid
        """
        if not self.params_file:
            raise ValueError('params_file is required to route')
        if self._validated:
            return self
        self._verify_input_files_exist()
        self._verify_output_directories_exist()
        if self.forcing == 'channel':
            if not self.discharge_files:
                raise ValueError('Provide discharge_dir (or discharge_files for explicit output paths)')
            for key in ('channel_state_init_file', 'dt_routing', 'dt_total'):
                if not getattr(self, key, None):
                    raise ValueError(f'{key} is required for channel routing')
            if len(self.discharge_files) != 1:
                raise ValueError('Channel routing requires exactly one entry in discharge_files')
        else:
            if self.transform == 'unit_hydrograph' and not self.uh_kernel_file:
                raise ValueError('uh_kernel_file is required when transform is unit_hydrograph')
            if not self.discharge_files:
                raise ValueError('Provide discharge_dir (or discharge_files for explicit output paths)')
            vlateral = self.vlateral_files
            grids = self.grid_runoff_files and self.grid_weights_file
            if vlateral and grids:
                raise ValueError('Provide vlateral_files or grid_runoff_files with grid_weights_file, not both')
            if not vlateral and not grids:
                raise ValueError('Provide vlateral_files or grid_runoff_files with grid_weights_file')
            n_inputs = len(vlateral) + len(self.grid_runoff_files or [])
            if len(self.discharge_files) != n_inputs:
                raise ValueError('Number of resolved discharge output files must match number of input files')
            if len(set(self.discharge_files)) != len(self.discharge_files):
                raise ValueError('discharge_files contains duplicate paths; each input file needs a distinct output')
        object.__setattr__(self, '_validated', True)
        return self

    def validate_runoff(self) -> Self:
        """
        Validate the options for preparing gridded runoff: grid_weights_file is set and every input path exists.
        Called by RunoffGaussianGrid. Returns immediately once the Configs has been validated. The contents of the
        input files are not read; call deep_validate for that.

        Raises:
            ValueError, FileNotFoundError: if any option is missing or invalid
        """
        if not self.grid_weights_file:
            raise ValueError('grid_weights_file is required to prepare runoff')
        if self._validated:
            return self
        self._verify_input_files_exist()
        object.__setattr__(self, '_validated', True)
        return self

    def deep_validate(self) -> Self:
        """
        Validate the contents of every input file that is set and their consistency with each other: the params
        file columns, types, and value ranges and that it is topologically sorted, that the grid weight table
        matches the params file and its proportions sum to 1 per river, and that the initial channel state has one
        row per river.

        Nothing calls this for you. It reads every input file, which is the same work routing is about to do, so
        run it once on inputs you have not checked before rather than on every route. ``validate_routing`` and
        ``validate_runoff`` check the options and the paths only.

        Raises:
            ValueError: if any file or combination of files is invalid
        """
        params_df = self._deep_validate_params_file() if self.params_file else None
        if self.grid_weights_file:
            self._deep_validate_grid_weights_file(self.grid_weights_file, params_df)
        if self.channel_state_init_file:
            self._deep_validate_channel_state_init_file(params_df)
        return self

    def _deep_validate_params_file(self) -> pd.DataFrame:
        # params df should be parquet with columns river_id, next_river_id, k, x
        # river_id should be non-null, integer, and unique
        # next_river_id should be non-null, integer, all -1 or positive, and exist in river_id (except for -1)
        # k should be positive float
        # x should be positive float less than or equal to 0.5
        try:
            params_df = pd.read_parquet(self.params_file)
        except Exception as e:
            raise ValueError('Error reading params file. Must be valid parquet file') from e
        rid = self.var_river_id
        if rid not in params_df.columns:
            raise ValueError(f'{self.params_file} missing {rid} column')
        if 'next_river_id' not in params_df.columns:
            raise ValueError(f'{self.params_file} missing next_river_id column')
        if 'k' not in params_df.columns:
            raise ValueError(f'{self.params_file} missing k column')
        if 'x' not in params_df.columns:
            raise ValueError(f'{self.params_file} missing x column')
        if np.any(params_df[rid].isnull()):
            raise ValueError(f'{self.params_file} {rid} column contains null values')
        if not pd.api.types.is_integer_dtype(params_df[rid]):
            raise ValueError(f'{self.params_file} {rid} column must be integer type')
        if not params_df[rid].is_unique:
            raise ValueError(f'{self.params_file} {rid} column must be unique')
        if np.any(params_df['next_river_id'].isnull()):
            raise ValueError(f'{self.params_file} next_river_id column contains null values')
        if not pd.api.types.is_integer_dtype(params_df['next_river_id']):
            raise ValueError(f'{self.params_file} next_river_id column must be integer type')
        if np.any(params_df['next_river_id'] < -1):
            raise ValueError(f'{self.params_file} next_river_id column must be -1 or positive integers')
        downstream_ids = set(params_df['next_river_id'].unique())
        river_ids = set(params_df[rid].unique())
        if not downstream_ids.issubset(river_ids.union({-1})):
            raise ValueError(f'{self.params_file} next_river_id values must exist in {rid} (except -1)')
        if np.any(params_df['k'].isnull()) or np.any(params_df['x'].isnull()):
            raise ValueError(f'{self.params_file} k and x columns must not contain null values')
        if np.any(params_df['k'] <= 0):
            raise ValueError(f'{self.params_file} k column must be positive')
        if np.any(params_df['x'] < 0) or np.any(params_df['x'] > 0.5):
            raise ValueError(f'{self.params_file} x column must be in the range [0, 0.5]')

        # dynamic coefficients are rebuilt in the kernel from K = alpha * Q ** beta
        if self.coeff == 'dynamic':
            for column in ('alpha', 'beta'):
                if column not in params_df.columns:
                    raise ValueError(f'{self.params_file} missing {column} column required when coeff is dynamic')
                if np.any(params_df[column].isnull()):
                    raise ValueError(f'{self.params_file} {column} column contains null values')
                if not pd.api.types.is_numeric_dtype(params_df[column]):
                    raise ValueError(f'{self.params_file} {column} column must be numeric type')
            if np.any(params_df['alpha'] <= 0):
                raise ValueError(f'{self.params_file} alpha column must be strictly positive')

        # check topological sort: every next_river_id must appear later in the table than its upstream
        river_id_index = {int(river_id): i for i, river_id in enumerate(params_df[rid])}
        for upstream_idx, ds_id in enumerate(params_df['next_river_id']):
            if int(ds_id) < 0:
                continue
            if river_id_index[int(ds_id)] <= upstream_idx:
                raise ValueError(f'{self.params_file} is not topologically sorted (upstream to downstream)')
        return params_df

    def _deep_validate_grid_weights_file(self, grid_weights_file: PathInput, params_df: pd.DataFrame | None) -> None:
        # weights should be netcdf with variables river_id, x_index, y_index, x, y, area_sqm, proportion.
        rid = self.var_river_id
        try:
            ds = xr.load_dataset(grid_weights_file)
        except Exception as e:
            raise ValueError('Error reading grid weights file. Must be valid netCDF file') from e
        expected_variables = (rid, 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion')
        for variable in expected_variables:
            if variable not in ds:
                raise ValueError(f'Grid weights file missing {variable} variable')
        if np.any(ds[rid].isnull()):
            raise ValueError(f'Grid weights {rid} variable contains null values')
        if not pd.api.types.is_integer_dtype(ds[rid].dtype):
            raise ValueError(f'Grid weights {rid} variable must be integer type')
        if params_df is not None and not np.isin(ds[rid].values, params_df[rid].to_numpy()).all():
            raise ValueError(f'Grid weights {rid} values must exist in params {rid}')
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

    def _deep_validate_channel_state_init_file(self, params_df: pd.DataFrame | None) -> None:
        # initial channel state should be parquet with 1 column named Q.
        # Q should be non-null, numeric, and non-negative.
        # it should be exactly the same shape as the number of rows in the params file and in the same order.
        try:
            state_df = pd.read_parquet(self.channel_state_init_file)
        except Exception as e:
            raise ValueError('Error reading initial state file. Must be valid parquet file') from e
        if 'Q' not in state_df.columns:
            raise ValueError('Initial state file missing Q column')
        if np.any(state_df['Q'].isnull()):
            raise ValueError('Initial state file Q column contains null values')
        if not pd.api.types.is_numeric_dtype(state_df['Q']):
            raise ValueError('Initial state file Q column must be numeric type')
        if np.any(state_df['Q'] < 0):
            raise ValueError('Initial state file Q column must be non-negative')
        if params_df is not None and state_df.shape[0] != params_df.shape[0]:
            raise ValueError(f'Initial state file must have the same number of rows as {self.params_file}')
        return


def _derive_valid_values(cls: type) -> dict[str, frozenset[str]]:
    result = {}
    for name, hint in get_type_hints(cls).items():
        if name.startswith('_'):
            continue
        if get_origin(hint) is Literal:
            result[name] = frozenset(get_args(hint))
    return result


def _derive_path_sets(cls: type) -> tuple[frozenset[str], frozenset[str]]:
    single, lists = set(), set()
    for name, hint in get_type_hints(cls).items():
        if name.startswith('_'):
            continue
        origin = get_origin(hint)
        if origin is list:
            args = get_args(hint)
            if args and get_origin(args[0]) is Literal:
                continue  # selector list (e.g. forcing), not a list of file paths
            lists.add(name)
        elif origin is types.UnionType:
            non_none = [a for a in get_args(hint) if a is not type(None)]
            if len(non_none) == 1 and get_origin(non_none[0]) is list:
                lists.add(name)  # PathList | None
            elif set(get_args(hint)) >= _PATH_INPUT_TYPES:
                single.add(name)  # PathInput or PathInput | None
    return frozenset(single), frozenset(lists)


Configs._SINGLE_PATH_FIELDS, Configs._LIST_PATH_FIELDS = _derive_path_sets(Configs)
Configs._VALID_VALUES = _derive_valid_values(Configs)
