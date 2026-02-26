import dataclasses
import os
import types
from dataclasses import dataclass, field
from typing import Any, ClassVar, Iterator, Literal, get_args, get_origin, get_type_hints

from ..types import PathInput, PathList, PathTypes

# Evaluated members of PathInput (str | Path) — used by _derive_path_sets for type inspection
_PATH_INPUT_TYPES: frozenset[type] = frozenset(get_args(PathInput))

__all__ = ['Configs', ]

_MISSING = object()  # sentinel for distinguishing between missing and None values in __contains__ and get() methods


@dataclass
class Configs:
    # Rule: annotate ALL file path fields with PathInput (or PathList) so _derive_path_sets() can
    # detect them via type inspection. Non-path fields use concrete types.

    # Core Routing Files
    routing_params_file: PathInput | None = None
    discharge_files: PathList = field(default_factory=list)
    channel_state_init_file: PathInput | None = None
    channel_state_final_file: PathInput | None = None

    # Time options
    dt_routing: int = 0
    dt_total: int = 0
    dt_discharge: int = 0
    dt_runoff: int = 0
    start_datetime: str = '1970-01-01'

    # For runoff transformation - used by TransformRouter subclasses
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

    def get(self, key: str, default: Any = None) -> Any:
        val = getattr(self, key, None)
        return val if (val is not None and val != '' and val != []) else default

    def update(self, data: dict[str, Any]) -> None:
        for k, v in data.items():
            setattr(self, k, v)

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        setattr(self, key, value)

    def __delitem__(self, key: str) -> None:
        fields = {f.name: f for f in dataclasses.fields(self)}
        if key not in fields:
            raise KeyError(key)
        f = fields[key]
        if f.default is not dataclasses.MISSING:
            setattr(self, key, f.default)
        elif f.default_factory is not dataclasses.MISSING:  # type: ignore[misc]
            setattr(self, key, f.default_factory())
        else:
            raise KeyError(f'Cannot delete required field: {key}')

    def __contains__(self, key: str) -> bool:
        val = getattr(self, key, _MISSING)
        if val is _MISSING:
            return False
        return val is not None and val != '' and val != []

    def keys(self) -> Iterator[str]:
        for f in dataclasses.fields(self):
            val = getattr(self, f.name)
            if val not in (None, '', []):
                yield f.name

    def items(self) -> Iterator[tuple[str, Any]]:
        for k in self.keys():
            yield k, getattr(self, k)


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
