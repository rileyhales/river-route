import os
import types
from dataclasses import dataclass, field
from typing import Any, ClassVar, Literal, get_args, get_origin, get_type_hints, Self

import numpy as np
import pandas as pd
import xarray as xr

from river_route.types import PathInput, PathList, PathTypes

# Evaluated members of PathInput (str | Path) — used by _derive_path_sets for type inspection
_PATH_INPUT_TYPES: frozenset[type] = frozenset(get_args(PathInput))

__all__ = ['Configs', ]

_MISSING = object()  # sentinel for distinguishing between missing and None values in __contains__ and get() methods


@dataclass
class Configs:
    """
    Accepts and validates configuration options for routing and runoff transform simulations
    """
    # annotate file path fields with PathInput or PathList
    # _derive_path_sets() will detect them by inspecting class annotations

    # Core Routing Files
    params_file: PathInput | None = None
    discharge_files: PathList = field(default_factory=list)
    channel_state_init_file: PathInput | None = None
    channel_state_final_file: PathInput | None = None

    # Time options
    dt_routing: int = 0
    dt_total: int = 0
    dt_discharge: int = 0
    dt_runoff: int = 0
    start_datetime: str = '1970-01-01'

    # For runoff transformation - used by TransformMuskingum subclasses
    catchment_runoff_files: PathList = field(default_factory=list)
    runoff_grid_files: PathList | None = field(default_factory=list)
    grid_weights_file: PathInput | None = None
    transformer_kernel_file: PathInput | None = None
    transformer_state_init_file: PathInput | None = None
    transformer_state_final_file: PathInput | None = None

    # Misc behavior that users may want to override
    log: bool = True
    progress_bar: bool = True
    log_level: Literal['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'] = 'INFO'
    log_stream: str = 'stdout'
    log_format: str = '%(asctime)s - %(levelname)s - %(message)s'
    runoff_processing_mode: Literal['sequential', 'ensemble'] = 'sequential'
    runoff_accumulation_type: Literal['incremental', 'cumulative'] = 'incremental'
    var_river_id: str = 'river_id'
    var_discharge: str = 'Q'
    var_x: str = 'x'
    var_y: str = 'y'
    var_t: str = 'time'
    var_catchment_runoff_variable: str = 'runoff'
    var_runoff_depth: str = 'ro'

    # special subset of auto-detected PathLists where the directory needs to exist, not the file
    _OUTPUT_FILES: ClassVar[frozenset[str]] = frozenset({
        'channel_state_final_file',
        'transformer_state_final_file',
    })
    _OUTPUT_FILE_LISTS: ClassVar[frozenset[str]] = frozenset({
        'discharge_files',
    })

    _ALWAYS_REQUIRED: ClassVar[tuple[str, ...]] = ('params_file', 'discharge_files')

    # Populated at module level below
    _SINGLE_PATH_FIELDS: ClassVar[frozenset[str]]
    _LIST_PATH_FIELDS: ClassVar[frozenset[str]]
    _VALID_VALUES: ClassVar[dict[str, frozenset[str]]]

    def __setattr__(self, name: str, value: object) -> None:
        allowed = type(self).__dict__.get('_VALID_VALUES', {}).get(name)
        if allowed is not None and value not in allowed:
            raise ValueError(f'{name} must be one of {sorted(allowed)}, got {value!r}')
        object.__setattr__(self, name, value)

    def __post_init__(self) -> None:
        # turn off progress bar if logging was turned off but progress was left at default on
        self.progress_bar = bool(self.log) and bool(self.progress_bar)
        # normalize paths given as strings to lists. make paths absolute.
        self.coerce_path_list_fields()
        self.absolutize_paths()
        self.verify_input_files_exist()
        self.verify_output_directories_exist()
        # make sure required fields are set
        for key in self._ALWAYS_REQUIRED:
            if getattr(self, key) in (None, '', []):
                raise ValueError(f'Missing required config: {key}')

    # --- path normalization and verification ---
    def coerce_path_list_fields(self) -> None:
        """Normalize any list-of-paths field given as a single string to [str]."""
        for key in self._LIST_PATH_FIELDS:
            val = getattr(self, key)
            if isinstance(val, PathTypes) and val:
                setattr(self, key, [str(val)])

    def absolutize_paths(self) -> None:
        """Convert all relative path fields to absolute paths in-place."""
        for key in self._SINGLE_PATH_FIELDS:
            val = getattr(self, key, None)
            if val:
                setattr(self, key, os.path.abspath(val))
        for key in self._LIST_PATH_FIELDS:
            val = getattr(self, key, [])
            if val:
                setattr(self, key, [os.path.abspath(p) for p in val])

    def verify_input_files_exist(self) -> None:
        """Raise FileNotFoundError for any set input path that does not exist."""
        input_single = self._SINGLE_PATH_FIELDS - self._OUTPUT_FILES
        input_list = self._LIST_PATH_FIELDS - self._OUTPUT_FILE_LISTS
        for key in input_single:
            val = getattr(self, key, None)
            if val and not os.path.exists(val):
                raise FileNotFoundError(f'{key} not found: {val}')
        for key in input_list:
            for path in getattr(self, key, []):
                if not os.path.exists(path):
                    raise FileNotFoundError(f'{key}: {path} not found')

    def verify_output_directories_exist(self) -> None:
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
                raise NotADirectoryError(f'Output directory not found for specified output path: {path}')

    # allow some dict like accessors
    def get(self, key: str, default: Any = None) -> Any:
        val = getattr(self, key, None)
        return val if (val is not None and val != '' and val != []) else default

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def __contains__(self, key: str) -> bool:
        val = getattr(self, key, _MISSING)
        if val is _MISSING:
            return False
        return val is not None and val != '' and val != []

    def deep_validate(self) -> Self:
        """Perform deep validation of file contents and inter-file consistency. Raise ValueError if any issues found."""
        # params df should be parquet with columns river_id, downstream_river_id, k, x
        # river_id should be non-null, integer, and unique
        # downstream_river_id should be non-null, integer, all -1 or positive, and exist in river_id (except for -1)
        # k should be positive float
        # x should be positive float less than or equal to 0.5
        # todo check that the rivers are topologically sorted
        try:
            params_df = pd.read_parquet(self.params_file)
        except Exception as e:
            raise ValueError('Error reading params file. Must be valid parquet file') from e
        if 'river_id' not in params_df.columns:
            raise ValueError(f'{self.params_file} missing river_id column')
        if 'downstream_river_id' not in params_df.columns:
            raise ValueError(f'{self.params_file} missing downstream_river_id column')
        if 'k' not in params_df.columns:
            raise ValueError(f'{self.params_file} missing k column')
        if 'x' not in params_df.columns:
            raise ValueError(f'{self.params_file} missing x column')
        if np.any(params_df['river_id'].isnull()):
            raise ValueError(f'{self.params_file} river_id column contains null values')
        if not pd.api.types.is_integer_dtype(params_df['river_id']):
            raise ValueError(f'{self.params_file} river_id column must be integer type')
        if not params_df['river_id'].is_unique:
            raise ValueError(f'{self.params_file} river_id column must be unique')
        if np.any(params_df['downstream_river_id'].isnull()):
            raise ValueError(f'{self.params_file} downstream_river_id column contains null values')
        if not pd.api.types.is_integer_dtype(params_df['downstream_river_id']):
            raise ValueError(f'{self.params_file} downstream_river_id column must be integer type')
        if np.any(params_df['downstream_river_id'] < -1):
            raise ValueError(f'{self.params_file} downstream_river_id column must be -1 or positive integers')
        downstream_ids = set(params_df['downstream_river_id'].unique())
        river_ids = set(params_df['river_id'].unique())
        if not downstream_ids.issubset(river_ids.union({-1})):
            raise ValueError(f'{self.params_file} downstream_river_id values must exist in river_id (except -1)')
        if np.any(params_df['k'] <= 0):
            raise ValueError(f'{self.params_file} k column must be positive')
        if np.any(params_df['x'] < 0) or np.any(params_df['x'] > 0.5):
            raise ValueError(f'{self.params_file} x column must be in the range [0, 0.5]')

        # weights should be netcdf with variables river_id, x_index, y_index, x, y, area_sqm, proportion.
        if self.grid_weights_file:
            try:
                ds = xr.load_dataset(self.grid_weights_file)
            except Exception as e:
                raise ValueError('Error reading grid weights file. Must be valid netCDF file') from e
            expected_variables = ('river_id', 'x_index', 'y_index', 'x', 'y', 'area_sqm', 'proportion',)
            for variable in expected_variables:
                if variable not in ds:
                    raise ValueError(f'Grid weights file missing {variable} variable')
            if np.any(ds['river_id'].isnull()):
                raise ValueError('Grid weights river_id variable contains null values')
            if not pd.api.types.is_integer_dtype(ds['river_id'].dtype):
                raise ValueError('Grid weights river_id variable must be integer type')
            if not set(ds['river_id'].values).issubset(river_ids):
                raise ValueError('Grid weights river_id values must exist in params river_id')
            for variable in expected_variables[1:]:
                if np.any(ds[variable].isnull()):
                    raise ValueError(f'Grid weights {variable} variable contains null values')
                if not pd.api.types.is_numeric_dtype(ds[variable].dtype):
                    raise ValueError(f'Grid weights {variable} variable must be numeric type')
            if np.any(ds['area_sqm'] <= 0):
                raise ValueError('Grid weights area_sqm variable must be positive')
            if np.any(ds['proportion'] <= 0) or np.any(ds['proportion'] > 1):
                raise ValueError('Grid weights proportion variable must be in the range (0, 1]')
            proportions_sum = ds.groupby('river_id')['proportion'].sum()
            if not np.allclose(proportions_sum.values, 1.0):
                raise ValueError('Grid weights proportion variable must sum to 1 for each river_id')

        # initial channel state should be parquet with 1 column named Q.
        # Q should be non-null, numeric, and non-negative.
        # it should be exactly the same shape as the number of rows in the params file and in the same order.
        if self.channel_state_init_file:
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
            if state_df.shape[0] != params_df.shape[0]:
                raise ValueError(f'Initial state file must have the same number of rows as {self.params_file}')
        return self


def _derive_valid_values(cls: type) -> dict[str, frozenset[str]]:
    """
    Inspect evaluated type hints for any field annotated with Literal and return a mapping
    of field name -> frozenset of allowed values. Used to populate _VALID_VALUES at import time.
    """
    result = {}
    for name, hint in get_type_hints(cls).items():
        if name.startswith('_'):
            continue
        if get_origin(hint) is Literal:
            result[name] = frozenset(get_args(hint))
    return result


def _derive_path_sets(cls: type) -> tuple[frozenset[str], frozenset[str]]:
    """
    Inspect evaluated type hints and return (single_path_fields, list_path_fields) for any
    field typed with PathInput (str | Path) or PathList (list[str | Path]), including Optional
    variants. Fields beginning with '_' are ignored.
    """
    single, lists = set(), set()
    for name, hint in get_type_hints(cls).items():
        if name.startswith('_'):
            continue
        origin = get_origin(hint)
        if origin is list:
            lists.add(name)
        elif origin is types.UnionType:
            non_none = [a for a in get_args(hint) if a is not type(None)]
            if len(non_none) == 1 and get_origin(non_none[0]) is list:
                lists.add(name)  # PathList | None
            elif _PATH_INPUT_TYPES <= set(get_args(hint)):
                single.add(name)  # PathInput or PathInput | None
    return frozenset(single), frozenset(lists)


Configs._SINGLE_PATH_FIELDS, Configs._LIST_PATH_FIELDS = _derive_path_sets(Configs)
Configs._VALID_VALUES = _derive_valid_values(Configs)
