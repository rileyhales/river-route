"""The form's fields, described from the Configs annotations, and the Configs options and Python of their values."""

import glob
import os
from pathlib import Path
from typing import Any, Literal, get_args, get_origin, get_type_hints

from ..configs import Configs

__all__ = ['describe_fields', 'form_values', 'python_of', 'shown_path', 'threads_from', 'typed_path']

RUNOFF_FORCINGS = ('catchment', 'grid', 'ecmwf_grib')
GRID_FORCINGS = ('grid', 'ecmwf_grib')
FORCING_GROUP = 'forcing options'  # the group of the forcing's own options, whose Python comment names the forcing
FORM_GROUPS = (  # each group's comment and fields, as the Configs of examples/example_usage.py has them
    (None, ('network_file',)),
    ('primary modeling choices', ('forcing',)),
    ('Compute Options', ('dt_routing', 'use_subreaches', 'use_subcycles', 'invalid_coefficients')),
    ('model state files', ('channel_state_init_file', 'channel_state_final_file')),
    ('Outputs Location', ('discharge_dir', 'discharge_to_devnull')),
    (
        FORCING_GROUP,
        (
            'start_datetime',
            'dt_discharge',
            'dt_total',
            'forcing_files',
            'grid_weights_file',
            'var_runoff',
            'var_x',
            'var_y',
            'var_t',
            'accumulation_type',
            'runoff_processing_mode',
            'force_positive_runoff',
        ),
    ),
    ('logging and validation options', ('progress_bar', 'log_level')),
)
FORCINGS_OF_FIELD = {  # the fields only some forcings read; the form shows every other field for every forcing
    'start_datetime': ('channel',),
    'dt_discharge': ('channel',),
    'dt_total': ('channel',),
    'forcing_files': RUNOFF_FORCINGS,
    'grid_weights_file': GRID_FORCINGS,
    'var_runoff': ('grid',),
    'var_x': ('grid',),
    'var_y': ('grid',),
    'var_t': ('grid',),
    'accumulation_type': ('grid',),
    'runoff_processing_mode': RUNOFF_FORCINGS,
    'force_positive_runoff': GRID_FORCINGS,
}
# forcing_files is text in the form, read by forcing_paths, and blank time steps are left unset
FORM_DEFAULTS = {'forcing_files': '', 'discharge_to_devnull': False, 'dt_discharge': '', 'dt_total': ''}


def describe_fields() -> list[dict[str, Any]]:
    """
    Each field of the form in order: its name, the comment of its group (None for the first), the forcings it is shown
    for (every one when None), its kind, its choices, its default, and what the file dialog picks for it. The
    groups follow the Configs of examples/example_usage.py, and a field the example leaves out follows the fields of the
    group it belongs with. Choices and defaults are read from the Configs annotations, so they never go stale.
    """
    hints = get_type_hints(Configs)
    defaults = Configs().to_dict() | FORM_DEFAULTS
    described = []
    for group, names in FORM_GROUPS:
        for name in names:
            hint, default = hints.get(name), defaults[name]
            literal = get_origin(hint) is Literal
            described.append(
                {
                    'name': name,
                    'group': group,
                    'forcings': FORCINGS_OF_FIELD.get(name),
                    'kind': 'choice' if literal else 'bool' if isinstance(default, bool) else 'text',
                    'choices': list(get_args(hint)) if literal else None,
                    'default': '' if default is None else default,
                    'pick': pick_of(name),
                }
            )
    return described


def pick_of(name: str) -> str | None:
    """
    What the file dialog picks for a field, read from the path fields of Configs: an output directory is a folder, an
    output file is a name in a folder, a list of files is a folder of them, as a glob pattern, and any other path is an
    existing file. None for a field that is not a path.
    """
    if name in Configs._OUTPUT_DIRS:
        return 'folder'
    if name in Configs._OUTPUT_FILES:
        return 'save'
    if name in Configs._LIST_PATH_FIELDS:
        return 'pattern'
    return 'file' if name in Configs._SINGLE_PATH_FIELDS else None


def shown_path(path: str, home: str, many: bool) -> str:
    """
    A path as its field shows it: relative to ``home`` when it is inside it, and as it is otherwise, or when it is
    blank. With ``many``, ``path`` is paths separated by os.pathsep, as forcing_files holds them, each shown so.
    """
    items = path.split(os.pathsep) if many else [path]
    shown = (
        os.path.relpath(item, home) if os.path.isabs(item) and Path(item).is_relative_to(home) else item
        for item in items
    )
    return os.pathsep.join(shown)


def typed_path(text: str, home: str, many: bool) -> str:
    """
    The absolute path of a field's text: relative to ``home``, absolute, or starting with ~, and blank for blank text.
    With ``many``, ``text`` is paths separated by os.pathsep, each read so.
    """
    items = text.split(os.pathsep) if many else [text]
    typed = (
        os.path.normpath(os.path.join(home, os.path.expanduser(item.strip()))) if item.strip() else '' for item in items
    )
    return os.pathsep.join(typed)


def forcing_paths(text: str) -> list[str]:
    """
    The forcing files of the form's forcing_files: paths or glob patterns separated by os.pathsep, each pattern standing
    for the files it matches in sorted order.

    Raises:
        ValueError: if a pattern matches no file
    """
    paths: list[str] = []
    for item in text.split(os.pathsep):
        pattern = os.path.expanduser(item.strip())
        matches = sorted(glob.glob(pattern)) if any(character in pattern for character in '*?[') else [pattern]
        if not matches:
            raise ValueError(f'forcing_files pattern {pattern} matches no file')
        paths.extend(path for path in matches if path)
    return paths


def form_values(form: dict[str, Any]) -> dict[str, Any]:
    """
    The Configs options of the form's values. An empty field keeps its default, a leading ~ in a path is the home
    directory, forcing_files is read by forcing_paths, and discharge_to_devnull discards the discharge.

    Raises:
        ValueError: if a whole number of seconds is not one, or a forcing_files pattern matches no file
    """
    hints = get_type_hints(Configs)
    values = {
        name: os.path.expanduser(value) if name in Configs._SINGLE_PATH_FIELDS else value
        for name, value in form.items()
        if value != ''
    }
    if 'forcing_files' in values:
        values['forcing_files'] = forcing_paths(values['forcing_files'])
    if values.pop('discharge_to_devnull', False):
        values['discharge_dir'] = os.devnull
    for name, value in values.items():
        if hints.get(name) is int or (name == 'dt_routing' and value != 'auto'):
            values[name] = int(value)
    return values


def threads_from(text: str) -> int:
    """
    The threads ``text`` names, a whole number, 1 or more.

    Raises:
        ValueError: if it names anything else
    """
    if not (text.strip().isdigit() and int(text) >= 1):
        raise ValueError(f'threads must be a whole number, 1 or more, not {text!r}')
    return int(text)


def python_of(options: dict[str, Any], threads: int) -> str:
    """
    The Python that routes ``options`` on ``threads`` threads as examples/example_usage.py does: under an
    ``if __name__ == '__main__':`` block, a Configs of the options of each group of FORM_GROUPS under its comment, the
    comment of FORCING_GROUP followed by the forcing chosen, then any others, routed by a Router on a ThreadPoolExecutor
    when threads is more than 1, else on no pool. Routing requires discharge_dir unless discharge_files is given, so
    it is always named, blank until it is chosen.
    """
    if 'discharge_files' not in options:
        options = {'discharge_dir': ''} | options
    grouped = [(comment, [name for name in names if name in options]) for comment, names in FORM_GROUPS]
    named = {name for _, names in grouped for name in names}
    grouped.append((None, [name for name in options if name not in named]))
    arguments = []
    for comment, names in grouped:
        if names and comment == FORCING_GROUP:
            arguments.append(f'        # {comment} - {options["forcing"]}')
        elif names and comment is not None:
            arguments.append(f'        # {comment}')
        arguments.extend(line for name in names for line in argument_lines(name, options[name]))
    lines = [
        'import contextlib',
        'from concurrent.futures import ThreadPoolExecutor',
        '',
        'import river_route as rr',
        '',
        "if __name__ == '__main__':",
        f'    threads = {threads}',
        '',
        '    conf = rr.Configs(',
        *arguments,
        '    )',
        '',
        '    pool_context = ThreadPoolExecutor(threads) if threads > 1 else contextlib.nullcontext()',
        '    with pool_context as pool:',
        '        rr.Router(conf).route(thread_pool=pool, threads=threads)',
    ]
    return '\n'.join(lines) + '\n'


def argument_lines(name: str, value: Any) -> list[str]:
    """The lines of a keyword argument of Configs: one, or one for each item of a list of more than one."""
    if isinstance(value, list) and len(value) > 1:
        return [f'        {name}=[', *(f'            {item!r},' for item in value), '        ],']
    return [f'        {name}={value!r},']
