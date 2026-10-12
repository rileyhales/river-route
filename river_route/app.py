"""A window that sets up Configs, queues them, routes each with rr route in a child process, and shows its output."""

import ctypes
import dataclasses
import functools
import glob
import itertools
import json
import os
import queue
import re
import subprocess
import sys
import threading
import time
import tkinter as tk
from collections.abc import Callable
from pathlib import Path
from tkinter import filedialog, font, ttk
from typing import Any, Literal, get_args, get_origin, get_type_hints

from .configs import Configs, is_dev_null

__all__ = ['open_window']

RUNOFF_FORCINGS = ('catchment', 'grid', 'ecmwf_grib')
GRID_FORCINGS = ('grid', 'ecmwf_grib')
FORCING_GROUP = 'forcing options'  # the group of the forcing's own options, whose Python comment names the forcing
FORM_GROUPS = (  # each group's comment and fields, as the Configs of examples/example_usage.py has them
    (None, ('network_file',)),
    ('primary modeling choices', ('forcing',)),
    (
        'Compute Options',
        ('dt_routing', 'conditioning', 'restrict_x_range', 'invalid_coefficients', 'dt_discharge', 'dt_total'),
    ),
    (
        FORCING_GROUP,
        (
            'start_datetime',
            'forcing_files',
            'grid_weights_file',
            'var_grid_runoff',
            'var_x',
            'var_y',
            'var_t',
            'grid_accumulation_type',
            'runoff_processing_mode',
            'force_positive_runoff',
        ),
    ),
    ('model state files', ('channel_state_init_file', 'channel_state_final_file')),
    ('where to place the results', ('discharge_dir', 'discharge_to_devnull')),
    ('logging and validation options', ('progress_bar', 'log_level')),
)
FORCINGS_OF_FIELD = {  # the fields only some forcings read; the form shows every other field for every forcing
    'start_datetime': ('channel',),
    'forcing_files': RUNOFF_FORCINGS,
    'grid_weights_file': GRID_FORCINGS,
    'var_grid_runoff': ('grid',),
    'var_x': ('grid',),
    'var_y': ('grid',),
    'var_t': ('grid',),
    'grid_accumulation_type': GRID_FORCINGS,
    'runoff_processing_mode': RUNOFF_FORCINGS,
    'force_positive_runoff': GRID_FORCINGS,
}
# forcing_files is text in the form, read by forcing_paths, and blank time steps are calculated from the inputs
FORM_DEFAULTS = {'forcing_files': '', 'discharge_to_devnull': False, 'dt_discharge': '', 'dt_total': ''}
PICKERS: dict[str, Callable[..., str]] = {  # the native dialog that picks each kind of path field
    'file': filedialog.askopenfilename,
    'save': filedialog.asksaveasfilename,
    'folder': filedialog.askdirectory,
    'pattern': filedialog.askdirectory,
}
CONFIG_FILE_TYPES = [('JSON', '*.json')]
PALETTE = {  # river-route-app's dark theme, from its src/style/global.css, with its color-mix() blends worked out
    'deepest': '#080d1d',  # the window
    'surface': '#101730',  # fields, and the bar across the top
    'elevated': '#151f3d',  # buttons, and the job chosen
    'hover': '#1d2a52',  # a button under the pointer
    'card': '#333b52',  # the form: surface with 18% white
    'section': '#242e4b',  # a group's heading: elevated with 8% white
    'border': '#2a3b6d',
    'border_light': '#3a4f8c',  # scrollbars
    'text': '#e4e9ff',
    'text_secondary': '#b6c2ea',
    'muted': '#8291bf',
    'white': '#ffffff',  # the text of the colored buttons
    'accent': '#7f8cff',
    'accent_active': '#8390ff',  # the accent under the pointer, 3% brighter
    'accent_border': '#444c8f',  # the accent with 35% black
    'accent_disabled': '#464f92',  # the accent at 52% opacity
    'success': '#4ac08e',
    'error': '#f27487',
    'error_active': '#fe7a8e',  # the error color under the pointer, 5% brighter
    'warning': '#f2bb63',
    'code': '#070b17',  # the preview and the log
    'code_text': '#dbe5ff',
    'code_border': '#334155',  # the preview's header: its border, the format chosen, and Copy
    'code_format': '#64748b',  # a format not chosen
    'code_format_chosen': '#e2e8f0',
    'copy_text': '#cbd5e1',
    'line_number': '#475569',
    'log_text': '#cbd5e1',
    'log_level': '#94a3b8',
    'log_warning': '#fbbf24',
    'log_error': '#fb7185',
}
FONT_FAMILIES = {  # river-route-app's open source fonts, each before the standard family Tk falls back on everywhere
    'main': ('Manrope', 'Helvetica'),
    'mono': ('JetBrains Mono', 'Courier'),
}
FONTS: dict[str, tuple[str, int, Literal['normal', 'bold']]] = {  # each font's family, size, and weight:
    'TkDefaultFont': ('main', 14, 'normal'),  # river-route-app's sizes in pixels, as Tk's points
    'TkTextFont': ('main', 14, 'normal'),
    'TkMenuFont': ('main', 14, 'normal'),
    'TkHeadingFont': ('main', 12, 'bold'),
    'TkFixedFont': ('mono', 13, 'normal'),
    'title': ('main', 18, 'bold'),
    'section': ('main', 12, 'bold'),
    'label': ('main', 13, 'bold'),
    'button': ('main', 14, 'bold'),
    'line_number': ('mono', 12, 'normal'),
}
LOG_COLORS = {  # the color of each log line that starts with the key
    'PROGRESS': 'log_level',
    'WARNING': 'log_warning',
    'ERROR': 'log_error',
    'CRITICAL': 'log_error',
    'Traceback': 'log_error',
}
CODE_COLORS = {  # the colors of the parts of the preview's Python or JSON
    'comment': '#60709f',
    'key': '#93c5fd',
    'string': '#bca6ff',
    'keyword': '#8ec3ff',
    'number': '#6fd8ff',
    'method': '#93c5fd',
    'type': '#fde68a',
    'punctuation': '#94a3b8',
}
CODE_TOKENS = re.compile(  # the parts of the preview's Python or JSON that CODE_COLORS colors
    r'(?P<comment>#.*)|(?P<key>"[^"\n]*"(?=\s*:))|(?P<string>\'[^\'\n]*\'|"[^"\n]*")'
    r'|(?P<keyword>\b(?:import|from|as|with|True|False|None|true|false|null)\b)|(?P<number>\b\d+(?:\.\d+)?\b)'
    r'|(?P<method>\.\w+)|(?P<type>\b[A-Z]\w*)|(?P<punctuation>[\[\](){},=:])'
)
PAGES = ('Config', 'Run Queue')
PREVIEW_FORMATS = ('Python', 'JSON')
QUEUE_COLUMNS = (  # each column of the queue: its name, heading, width, and whether it takes the width left over
    ('number', '#', 70, False),
    ('name', 'config', 280, True),
    ('status', 'status', 120, False),
    ('elapsed', 'elapsed', 110, False),
)
STATUS_COLORS = {'pending': 'muted', 'running': 'warning', 'complete': 'success', 'failed': 'error', 'stopped': 'muted'}
FINISHED = ('complete', 'failed', 'stopped')
LOGO = Path(__file__).with_name('logo.svg')
NETWORK_FILE_NAME = 'network.parquet'  # the network file the window selects when the home directory holds one
APPLICATION_NAME = 'River Route'  # the name macOS shows for the app, as in the Dock when the pointer is over it
CORE_SERVICES = '/System/Library/Frameworks/CoreServices.framework/CoreServices'
CORE_FOUNDATION = '/System/Library/Frameworks/CoreFoundation.framework/CoreFoundation'
FOUNDATION = '/System/Library/Frameworks/Foundation.framework/Foundation'
OBJECTIVE_C_RUNTIME = '/usr/lib/libobjc.A.dylib'
CF_STRING_ENCODING_UTF8 = 0x08000100
LS_DEFAULT_SESSION = -2  # LaunchServices' kLSDefaultSessionID, the session of this process
LINE_ENDINGS = re.compile(rb'(\r|\n)')
POLL_MILLISECONDS = 100  # how often the window shows what the run has printed
WHEEL_PIXELS = 40  # how far the form scrolls for each notch of a mouse wheel, as far as ttk's scrollable widgets do


def describe_fields() -> list[dict[str, Any]]:
    """
    Each field of the form in order: its name, the comment of its group (None for the first), the forcings it is shown
    for (every one when None), its kind, its choices, its default, and what the file dialog picks for it. The
    groups follow the Configs of examples/example_usage.py, and a field the example leaves out follows the fields of the
    group it belongs with. Choices and defaults are read from the Configs annotations, so they never go stale.
    """
    hints = get_type_hints(Configs)
    defaults = {f.name: f.default for f in dataclasses.fields(Configs) if f.init} | FORM_DEFAULTS
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
    when threads is more than 1, else on no pool.
    """
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


def read_output(process: subprocess.Popen[bytes], output: queue.Queue[tuple[str, str]]) -> None:
    """
    Put what the process prints on ``output`` as line and replace events until it exits, then its exit code as an end
    event. A line ended by a carriage return, as tqdm redraws its bar, is replaced by the line after it.
    """
    stdout = process.stdout
    assert stdout is not None
    replaces = False  # the line before ended with a carriage return, so the next takes its place
    unfinished = b''
    for chunk in iter(lambda: stdout.read(65536), b''):
        *finished, unfinished = LINE_ENDINGS.split(unfinished + chunk)
        for line, ending in zip(finished[::2], finished[1::2], strict=True):
            output.put(('replace' if replaces else 'line', line.decode(errors='replace')))
            replaces = ending == b'\r'
    if unfinished:
        output.put(('replace' if replaces else 'line', unfinished.decode(errors='replace')))
    output.put(('end', str(process.wait())))


def apply_theme(root: tk.Tk) -> dict[str, font.Font]:
    """
    Style every widget as river-route-app styles its page, on ttk's clam theme, whose colors can be set on every
    platform: in PALETTE's colors, and in the fonts of FONTS, each in the first family of FONT_FAMILIES this machine
    has, or else in the standard family Tk falls back on. Returns the fonts made, which the window keeps, since a font
    is deleted once nothing refers to it.
    """
    installed = set(font.families(root))
    families = {
        kind: next((name for name in names if name in installed), names[-1]) for kind, names in FONT_FAMILIES.items()
    }
    fonts: dict[str, font.Font] = {}
    for name, (kind, size, weight) in FONTS.items():
        if name.startswith('Tk'):
            font.nametofont(name).configure(family=families[kind], size=size, weight=weight)
        else:
            fonts[name] = font.Font(root, family=families[kind], size=size, weight=weight)
    p = PALETTE
    style = ttk.Style(root)
    style.theme_use('clam')
    style.configure(
        '.',
        background=p['deepest'],
        foreground=p['text'],
        bordercolor=p['border'],
        troughcolor=p['deepest'],
        fieldbackground=p['surface'],
        lightcolor=p['surface'],
        darkcolor=p['surface'],
        insertcolor=p['text'],
        selectbackground=p['accent'],
        selectforeground=p['white'],
        arrowcolor=p['muted'],
        focuscolor=p['accent'],
    )
    style.configure('TLabel', foreground=p['text_secondary'])
    style.configure('TSeparator', background=p['border'])
    style_fields(style, root)
    style_buttons(style, fonts)
    style_panes(style, fonts)
    root.configure(background=p['deepest'])
    return fonts


def style_fields(style: ttk.Style, root: tk.Tk) -> None:
    """Color the fields, text, menus, and checkboxes, with the accent around the one typed in, and the scrollbars."""
    p = PALETTE
    focused = [('focus', p['accent'])]
    style.configure('TEntry', padding=(10, 7), foreground=p['text'])
    style.map(
        'TEntry',
        fieldbackground=[('disabled', p['deepest'])],
        foreground=[('disabled', p['muted'])],
        bordercolor=focused,
        lightcolor=focused,
    )
    style.configure('TCombobox', padding=(10, 7), foreground=p['text'], background=p['surface'], arrowsize=14)
    style.map(
        'TCombobox',
        fieldbackground=[('readonly', p['surface'])],
        foreground=[('readonly', p['text'])],
        selectbackground=[('readonly', p['surface'])],
        selectforeground=[('readonly', p['text'])],
        background=[('active', p['surface'])],
        arrowcolor=[('active', p['accent'])],
        bordercolor=focused,
        lightcolor=focused,
    )
    for option, color in (('background', 'surface'), ('foreground', 'text'), ('selectBackground', 'accent')):
        root.option_add(f'*TCombobox*Listbox.{option}', p[color])  # the menu a combobox drops down
    root.option_add('*TCombobox*Listbox.selectForeground', p['white'])
    root.option_add('*TCombobox*Listbox.font', 'TkTextFont')
    style.configure('TSpinbox', padding=(10, 7), foreground=p['text'], background=p['surface'], arrowsize=12)
    style.map('TSpinbox', arrowcolor=[('active', p['accent'])], bordercolor=focused, lightcolor=focused)
    for name, background in (('TCheckbutton', 'deepest'), ('Card.TCheckbutton', 'card')):
        style.configure(
            name,
            focusthickness=0,
            background=p[background],
            indicatorbackground=p['surface'],
            indicatorforeground=p['white'],
            upperbordercolor=p['border'],
            lowerbordercolor=p['border'],
        )
        style.map(name, background=[('active', p[background])], indicatorbackground=[('selected', p['accent'])])
    for orient, sticky in (('Vertical', 'ns'), ('Horizontal', 'ew')):  # thin and without arrows, as a browser draws
        name, thumb = f'{orient}.TScrollbar', f'{orient}.Scrollbar.thumb'
        style.layout(
            name, [(f'{orient}.Scrollbar.trough', {'sticky': sticky, 'children': [(thumb, {'sticky': 'nswe'})]})]
        )
        style.configure(name, arrowsize=8, gripcount=0, background=p['border_light'], bordercolor=p['deepest'])
        style.configure(name, troughcolor=p['deepest'], lightcolor=p['border_light'], darkcolor=p['border_light'])
        style.map(name, background=[('active', p['muted'])])


def style_buttons(style: ttk.Style, fonts: dict[str, font.Font]) -> None:
    """
    Color the buttons: a plain one raised off the window, Accent.TButton for a page's action, Danger.TButton for Stop,
    and the buttons that choose a page or the preview's format, the one chosen marked, and Copy.
    """
    p = PALETTE
    for name, background, foreground, border, active, disabled in (
        ('TButton', 'elevated', 'text', 'border', 'hover', 'elevated'),
        ('Accent.TButton', 'accent', 'white', 'accent_border', 'accent_active', 'accent_disabled'),
        ('Danger.TButton', 'error', 'white', 'error', 'error_active', 'elevated'),
    ):
        style.configure(name, padding=(14, 7), focusthickness=0, font=fonts['button'], bordercolor=p[border])
        style.configure(name, background=p[background], foreground=p[foreground], lightcolor=p[background])
        style.configure(name, darkcolor=p[background])
        states = [('disabled', p[disabled]), ('active', p[active])]
        style.map(name, background=states, lightcolor=states, darkcolor=states)
        style.map(name, foreground=[('disabled', p['muted'])], bordercolor=[('disabled', p['border'])])
    for name, background, foreground, typeface in (
        ('Page.TButton', 'surface', 'text_secondary', fonts['label']),
        ('CurrentPage.TButton', 'accent', 'white', fonts['label']),
        ('Format.TButton', 'code', 'code_format', fonts['button']),
        ('CurrentFormat.TButton', 'code_border', 'code_format_chosen', fonts['button']),
        ('Copy.TButton', 'code_border', 'copy_text', 'TkDefaultFont'),
    ):
        color = p[background]  # its background and its borders, which keep their color under the pointer
        style.configure(name, padding=(14, 6), focusthickness=0, font=typeface, foreground=p[foreground])
        style.configure(name, background=color, bordercolor=color, lightcolor=color, darkcolor=color)
        unchanged = [('active', color)]
        style.map(name, background=unchanged, bordercolor=unchanged, lightcolor=unchanged, darkcolor=unchanged)
        style.map(name, foreground=[('active', p['text'])])


def style_panes(style: ttk.Style, fonts: dict[str, font.Font]) -> None:
    """Color the bar across the top, the form's card, the preview's header, and the queue."""
    p = PALETTE
    style.configure('Bar.TFrame', background=p['surface'])
    style.configure('Bar.TLabel', background=p['surface'], foreground=p['text_secondary'])
    style.configure('Title.TLabel', background=p['surface'], foreground=p['text'], font=fonts['title'])
    style.configure('Card.TFrame', background=p['card'])
    style.configure('Card.TLabel', background=p['card'], foreground=p['text'], font=fonts['label'])
    style.configure('Field.TLabel', foreground=p['text'], font=fonts['label'])  # a field's name off the card
    style.configure('Code.TFrame', background=p['code'])
    style.configure('Code.TSeparator', background=p['code_border'])
    style.configure(
        'Treeview',
        rowheight=40,
        background=p['deepest'],
        fieldbackground=p['deepest'],
        foreground=p['text'],
        bordercolor=p['deepest'],
        lightcolor=p['deepest'],
        darkcolor=p['deepest'],
    )
    style.map('Treeview', background=[('selected', p['elevated'])])
    style.configure(
        'Treeview.Heading',
        padding=(12, 8),
        font=fonts['section'],
        background=p['deepest'],
        foreground=p['muted'],
        bordercolor=p['border'],
        lightcolor=p['deepest'],
        darkcolor=p['deepest'],
    )
    style.map('Treeview.Heading', background=[('active', p['deepest'])])


@dataclasses.dataclass
class Job:
    """
    A run of the queue: ``config``, a config file's path or the options of a Configs, routed by rr route. Its status is
    one of STATUS_COLORS, and its lines are what it printed, each with the level that colors it.
    """

    number: int
    name: str
    config: str | dict[str, Any]
    status: str = 'pending'
    lines: list[tuple[str, str | None]] = dataclasses.field(default_factory=list)
    started: float | None = None  # the time.monotonic() it started at
    ended: float | None = None  # and ended at

    def elapsed(self) -> str:
        """How long it has routed, or routed, as m:ss or h:mm:ss, and blank before it starts."""
        if self.started is None:
            return ''
        ended = time.monotonic() if self.ended is None else self.ended
        minutes, seconds = divmod(int(ended - self.started), 60)
        hours, minutes = divmod(minutes, 60)
        return f'{hours}:{minutes:02}:{seconds:02}' if hours else f'{minutes}:{seconds:02}'


class RouteWindow:
    """
    The window, a page at a time beneath a bar of the logo and a button for each page. The Config page holds a home
    directory that every file dialog opens in and the form of describe_fields, a section for each group headed by its
    comment, with the fields the chosen forcing does not read hidden, above Load config, Save config, and Add to queue,
    and beside it the form's Configs as Python or JSON. The Run Queue page lists the jobs queued, which Run routes one
    after another with rr route in a child process so that Stop can end it, and beside them the log of the job chosen.
    """

    def __init__(self, root: tk.Tk, home: str, threads: int) -> None:
        self.root = root
        self.fonts = apply_theme(root)
        self.fields = describe_fields()
        self.variables: dict[str, tk.StringVar | tk.BooleanVar] = {}
        self.rows: dict[str, list[ttk.Widget]] = {}  # the widgets of each field's row of the form
        self.headings: dict[str | None, tk.Frame] = {}  # the band above each group, holding its comment
        self.kept: dict[str, Any] = {}  # the options of the config file loaded that the form has no field for
        self.form_name = 'form'  # the name of the config file loaded or saved last, which names the form's jobs
        self.home = tk.StringVar(value=home)  # the folder every file dialog opens in, and paths are shown relative to
        self.shown: dict[str, tk.StringVar] = {}  # the text of each path field, its absolute path relative to home
        self.converting = False  # whether a path field's absolute path and its text are being matched to each other
        self.threads = tk.StringVar(value=str(threads))  # a resource of each run, so never saved in a config file
        self.preview_format = tk.StringVar(value=PREVIEW_FORMATS[0])
        self.preview_requested = False  # whether the preview is shown again once Tk is idle
        self.preview_text = ''  # what the preview shows, without its line numbers, which Copy copies
        self.jobs: list[Job] = []
        self.job_numbers = itertools.count(1)
        self.running_job: Job | None = None
        self.queue_running = False  # whether the next pending job starts when the running one ends
        self.stop_requested = False  # whether Stop ended the running job
        self.process: subprocess.Popen[bytes] | None = None
        self.output: queue.Queue[tuple[str, str]] = queue.Queue()
        self.logo, self.icon = None, None
        if tk.TkVersion >= 9.0:  # Tk reads SVG from 9.0 on
            self.logo = tk.PhotoImage(file=LOGO, format='svg -scaletoheight 40')
            self.icon = tk.PhotoImage(file=LOGO, format='svg -scaletoheight 256')
            root.iconphoto(True, self.icon)
        root.title(APPLICATION_NAME)
        width, height = root.winfo_screenwidth(), root.winfo_screenheight()
        root.geometry(f'{width}x{int(height * 0.9)}+0+0')
        self._lay_out_pages()
        self.variables['forcing'].trace_add('write', lambda *_: self.show_fields_of_forcing())
        self.variables['discharge_to_devnull'].trace_add('write', lambda *_: self.show_discharge_dir())
        for variable in (*self.variables.values(), self.threads, self.preview_format):
            variable.trace_add('write', lambda *_: self.request_preview())
        self.show_fields_of_forcing()
        self._select_network_file_in(home)
        self.home.trace_add('write', lambda *_: self._show_paths())
        self.show_preview()
        root.bind_all('<MouseWheel>', self._scroll_by_wheel)
        if tk.TkVersion >= 9.0:  # trackpads scroll by TouchpadScroll events from Tk 9.0 on
            root.bind_all('<TouchpadScroll>', self._scroll_by_touchpad)
        root.protocol('WM_DELETE_WINDOW', self.close)
        root.after(POLL_MILLISECONDS, self.show_output)
        root.after_idle(self._name_application)  # once Tk has registered the window with the operating system

    def _name_application(self) -> None:
        if not name_application():
            print(
                f'WARNING - macOS did not take the name {APPLICATION_NAME}, so it shows python instead', file=sys.stderr
            )

    def _select_network_file_in(self, home: str) -> None:
        """
        Select the file of the home directory named exactly NETWORK_FILE_NAME as network_file, if it holds one. Its
        listing is searched rather than the path tested, since a case-insensitive file system finds Network.parquet too.
        """
        path = os.path.join(home, NETWORK_FILE_NAME)
        network_file = self.variables['network_file']
        assert isinstance(network_file, tk.StringVar), 'network_file is a path, so its value is text'
        if NETWORK_FILE_NAME in os.listdir(home) and os.path.isfile(path):
            network_file.set(path)

    def _lay_out_pages(self) -> None:
        """The bar across the top, and beneath it the pages of PAGES, each in the same place, the Config page shown."""
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(2, weight=1)
        self.page_buttons: dict[str, ttk.Button] = {}
        self._lay_out_bar()
        ttk.Separator(self.root).grid(row=1, column=0, sticky='ew')
        self.pages: dict[str, ttk.Frame] = {}
        for name in PAGES:
            page = ttk.Frame(self.root)
            page.grid(row=2, column=0, sticky='nsew')
            page.columnconfigure((0, 1), weight=1, uniform='halves')  # each half of a page takes half the window
            page.rowconfigure(0, weight=1)
            self.pages[name] = page
        self._lay_out_form(self.pages['Config'])
        self._lay_out_preview(self.pages['Config'])
        self._lay_out_queue(self.pages['Run Queue'])
        self._lay_out_log(self.pages['Run Queue'])
        self.show_page('Config')

    def _lay_out_bar(self) -> None:
        """The logo and title, and a button for each page."""
        bar = ttk.Frame(self.root, padding=(24, 12), style='Bar.TFrame')
        bar.grid(row=0, column=0, sticky='ew')
        if self.logo is not None:
            ttk.Label(bar, image=self.logo, style='Bar.TLabel').pack(side='left', padx=(0, 12))
        ttk.Label(bar, text=APPLICATION_NAME, style='Title.TLabel').pack(side='left', padx=(0, 28))
        for name in PAGES:
            self.page_buttons[name] = ttk.Button(bar, text=name, command=functools.partial(self.show_page, name))
            self.page_buttons[name].pack(side='left', padx=(0, 4))

    def show_page(self, name: str) -> None:
        """Show the page of PAGES named ``name``, its button in the accent."""
        self.pages[name].tkraise()
        for page, button in self.page_buttons.items():
            button.configure(style='CurrentPage.TButton' if page == name else 'Page.TButton')

    def _lay_out_form(self, page: ttk.Frame) -> None:
        """
        The home directory, and under it the form, scrolling on a card, above Load config, Save config, and Add to
        queue, on the Config page's left half.
        """
        side = ttk.Frame(page)
        side.grid(row=0, column=0, sticky='nsew')
        actions = ttk.Frame(side, padding=(24, 16))
        actions.pack(side='bottom', fill='x')  # packed first, so the form scrolls above it rather than covering it
        ttk.Button(actions, text='Load config', command=self.load_config).pack(side='left')
        ttk.Button(actions, text='Save config', command=self.save_config).pack(side='left', padx=(8, 0))
        add = ttk.Button(actions, text='Add to queue', command=self.add_form_to_queue, style='Accent.TButton')
        add.pack(side='right')
        home = ttk.Frame(side, padding=(24, 24, 20, 0))  # its right edge lines up with the card's, beside the scrollbar
        home.pack(side='top', fill='x')
        home.columnconfigure(1, weight=1)
        ttk.Label(home, text='home directory', style='Field.TLabel').grid(row=0, column=0, padx=(0, 16))
        ttk.Entry(home, textvariable=self.home).grid(row=0, column=1, sticky='ew')
        ttk.Button(home, text='Browse', command=self.browse_home).grid(row=0, column=2, padx=(6, 0))
        p = PALETTE
        self.form_canvas = tk.Canvas(side, background=p['card'], yscrollincrement=1, highlightthickness=1)
        self.form_canvas.configure(highlightbackground=p['border'], highlightcolor=p['border'])
        scrollbar = ttk.Scrollbar(side, orient='vertical', command=self.form_canvas.yview)
        self.form_canvas.configure(yscrollcommand=scrollbar.set)
        form = ttk.Frame(self.form_canvas, padding=20, style='Card.TFrame')
        window = self.form_canvas.create_window((0, 0), window=form, anchor='nw')
        form.bind('<Configure>', lambda _: self.form_canvas.configure(scrollregion=self.form_canvas.bbox('all')))
        self.form_canvas.bind('<Configure>', lambda event: self.form_canvas.itemconfigure(window, width=event.width))
        scrollbar.pack(side='right', fill='y', padx=(4, 8), pady=(20, 0))
        self.form_canvas.pack(side='left', fill='both', expand=True, padx=(24, 0), pady=(20, 0))
        self._build_form(form)

    def _lay_out_preview(self, page: ttk.Frame) -> None:
        """
        The form's Configs in the format of PREVIEW_FORMATS chosen above it, which Copy copies, on the right half of the
        Config page.
        """
        frame = ttk.Frame(page, style='Code.TFrame')
        frame.grid(row=0, column=1, sticky='nsew')
        ttk.Separator(frame, orient='vertical').pack(side='left', fill='y')
        choices = ttk.Frame(frame, padding=(16, 8), style='Code.TFrame')
        choices.pack(side='top', fill='x')
        ttk.Separator(frame, style='Code.TSeparator').pack(side='top', fill='x')
        self.format_buttons: dict[str, ttk.Button] = {}
        for name in PREVIEW_FORMATS:
            self.format_buttons[name] = ttk.Button(
                choices, text=name, style='Format.TButton', command=functools.partial(self.preview_format.set, name)
            )
            self.format_buttons[name].pack(side='left', padx=(0, 2))
        ttk.Button(choices, text='Copy', command=self.copy_preview, style='Copy.TButton').pack(side='right')
        self.preview = self._add_text(frame, 'none', PALETTE['code_text'])
        for kind, color in CODE_COLORS.items():
            self.preview.tag_configure(kind, foreground=color)
        self.preview.tag_configure('line_number', foreground=PALETTE['line_number'], font=self.fonts['line_number'])
        self.preview.tag_configure('ERROR', foreground=PALETTE['log_error'])

    def _lay_out_queue(self, page: ttk.Frame) -> None:
        """
        The jobs of the queue, under Load configs, Remove, Requeue, and Clear finished, and above the threads each job
        routes with, Run, and Stop, on the left half of the Run Queue page.
        """
        side = ttk.Frame(page)
        side.grid(row=0, column=0, sticky='nsew')
        edits = ttk.Frame(side, padding=(24, 16))
        edits.pack(side='top', fill='x')
        for text, command in (
            ('Load configs', self.load_configs),
            ('Remove', self.remove_job),
            ('Requeue', self.requeue_job),
            ('Clear finished', self.clear_finished),
        ):
            ttk.Button(edits, text=text, command=command).pack(side='left', padx=(0, 8))
        actions = ttk.Frame(side, padding=(24, 16))
        actions.pack(side='bottom', fill='x')  # packed first, so the queue fills what is left above it
        self.stop_button = ttk.Button(actions, text='Stop', command=self.stop_queue, style='Danger.TButton')
        self.stop_button.state(['disabled'])
        self.stop_button.pack(side='right')
        self.run_button = ttk.Button(actions, text='Run', command=self.run_queue, style='Accent.TButton')
        self.run_button.pack(side='right', padx=(0, 8))
        threads = ttk.Spinbox(actions, from_=1, to=os.cpu_count() or 1, textvariable=self.threads, width=4)
        threads.pack(side='right', padx=(0, 16))
        ttk.Label(actions, text='threads').pack(side='right', padx=(0, 8))
        self.tree = ttk.Treeview(side, columns=[column[0] for column in QUEUE_COLUMNS], show='headings')
        self.tree.configure(selectmode='browse')
        for column, heading, width, stretch in QUEUE_COLUMNS:
            self.tree.heading(column, text=heading.upper(), anchor='w')
            self.tree.column(column, width=width, stretch=stretch, anchor='w')
        for status, color in STATUS_COLORS.items():
            self.tree.tag_configure(status, foreground=PALETTE[color])
        self.tree.bind('<<TreeviewSelect>>', lambda _: self.show_job_log())
        scrollbar = ttk.Scrollbar(side, orient='vertical', command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side='right', fill='y', padx=(4, 8))
        self.tree.pack(side='left', fill='both', expand=True, padx=(24, 0))

    def _lay_out_log(self, page: ttk.Frame) -> None:
        """The log of the job chosen in the queue, its lines colored by level, on the Run Queue page's right half."""
        frame = ttk.Frame(page, style='Code.TFrame')
        frame.grid(row=0, column=1, sticky='nsew')
        ttk.Separator(frame, orient='vertical').pack(side='left', fill='y')
        self.log = self._add_text(frame, 'word', PALETTE['log_text'])
        for level, color in LOG_COLORS.items():
            self.log.tag_configure(level, foreground=PALETTE[color])

    @staticmethod
    def _add_text(frame: ttk.Frame, wrap: Literal['none', 'word'], color: str) -> tk.Text:
        """
        A text in ``color`` on the code background that can be read but not typed in, filling ``frame`` beside its
        scrollbar, and above another when its lines are not wrapped.
        """
        p = PALETTE
        text = tk.Text(
            frame,
            wrap=wrap,
            font='TkFixedFont',
            state='disabled',
            padx=16,
            pady=12,
            spacing1=2,
            spacing3=2,
            borderwidth=0,
            highlightthickness=0,
            background=p['code'],
            foreground=color,
            insertbackground=p['text'],
            selectbackground=p['accent'],
            selectforeground=p['white'],
        )
        scrollbar = ttk.Scrollbar(frame, orient='vertical', command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side='right', fill='y')
        if wrap == 'none':
            across = ttk.Scrollbar(frame, orient='horizontal', command=text.xview)
            text.configure(xscrollcommand=across.set)
            across.pack(side='bottom', fill='x')
        text.pack(side='left', fill='both', expand=True)
        return text

    def _build_form(self, form: ttk.Frame) -> None:
        """Each field's row, with a heading above each group but the first."""
        form.columnconfigure(1, weight=1)
        grid_row, group = 0, None
        for field in self.fields:
            if field['group'] != group:
                group = field['group']
                self.headings[group] = self._add_heading(form, grid_row, group)
                grid_row += 1
            browse = None if field['pick'] is None else (lambda field=field: self.browse(field))
            self.rows[field['name']] = self._add_row(form, grid_row, field['name'], self._control(form, field), browse)
            grid_row += 1

    def _add_heading(self, form: ttk.Frame, grid_row: int, comment: str) -> tk.Frame:
        """A band across the form with a bar of the accent at its left, holding its comment in capitals."""
        p = PALETTE
        band = tk.Frame(form, background=p['section'])
        band.grid(row=grid_row, column=0, columnspan=3, sticky='ew', pady=(20, 10))
        tk.Frame(band, background=p['accent'], width=4).pack(side='left', fill='y')
        heading = tk.Label(
            band, text=comment.upper(), font=self.fonts['section'], background=p['section'], foreground=p['text']
        )
        heading.pack(side='left', padx=10, pady=7)
        return band

    def _control(self, form: ttk.Frame, field: dict[str, Any]) -> ttk.Widget:
        """The widget that holds a field's value, a checkbox, a menu of its choices, or text, and its variable."""
        if field['kind'] == 'bool':
            boolean = tk.BooleanVar(value=field['default'])
            self.variables[field['name']] = boolean
            return ttk.Checkbutton(form, variable=boolean, style='Card.TCheckbutton')
        text = tk.StringVar(value=str(field['default']))
        self.variables[field['name']] = text
        if field['kind'] != 'choice' and field['pick'] is None:
            return ttk.Entry(form, textvariable=text)
        if field['kind'] != 'choice':  # a path: its absolute path is kept in text, and shown relative to home
            name, shown = field['name'], tk.StringVar(value='')
            self.shown[name] = shown
            text.trace_add('write', lambda *_: self._show_path(name))
            shown.trace_add('write', lambda *_: self._take_typed_path(name))
            return ttk.Entry(form, textvariable=shown)
        menu = ttk.Combobox(form, textvariable=text, values=field['choices'], state='readonly')
        # the form scrolls under the pointer, where a combobox would otherwise take the scroll to change its choice
        menu.bind('<MouseWheel>', self._scroll_by_wheel)
        if tk.TkVersion >= 9.0:
            menu.bind('<TouchpadScroll>', self._scroll_by_touchpad)
        return menu

    def _show_paths(self) -> None:
        for name in self.shown:
            self._show_path(name)

    def _show_path(self, name: str) -> None:
        """Show a path field's absolute path relative to home, unless it was just read from what was typed there."""
        if self.converting:
            return
        home = os.path.abspath(os.path.expanduser(self.home.get()))
        self.converting = True
        self.shown[name].set(shown_path(str(self.variables[name].get()), home, name in Configs._LIST_PATH_FIELDS))
        self.converting = False

    def _take_typed_path(self, name: str) -> None:
        """Keep the absolute path of what was typed in a path field, leaving the text as it was typed."""
        if self.converting:
            return
        home = os.path.abspath(os.path.expanduser(self.home.get()))
        path = self.variables[name]
        assert isinstance(path, tk.StringVar), f'{name} is a path, so its value is text'
        self.converting = True
        path.set(typed_path(self.shown[name].get(), home, name in Configs._LIST_PATH_FIELDS))
        self.converting = False

    @staticmethod
    def _add_row(
        form: ttk.Frame, grid_row: int, text: str, control: ttk.Widget, browse: Callable[[], None] | None
    ) -> list[ttk.Widget]:
        """Lay out a row of the form, its label, its control, and a Browse button when it is a path, and return them."""
        label = ttk.Label(form, text=text, style='Card.TLabel')
        label.grid(row=grid_row, column=0, sticky='w', padx=(0, 16), pady=6)
        sticky = 'w' if isinstance(control, ttk.Checkbutton) else 'ew'
        control.grid(row=grid_row, column=1, sticky=sticky, pady=6)
        if browse is None:
            return [label, control]
        button = ttk.Button(form, text='Browse', command=browse)
        button.grid(row=grid_row, column=2, padx=(6, 0), pady=6)
        return [label, control, button]

    def _scroll_by_wheel(self, event: tk.Event) -> str | None:
        """Scroll the form WHEEL_PIXELS for each notch of a mouse wheel, of 120, while the pointer is over it."""
        return self._scroll_form(event, round(-event.delta / 120 * WHEEL_PIXELS))

    def _scroll_by_touchpad(self, event: tk.Event) -> str | None:
        """
        Scroll the form as far as a trackpad moves while the pointer is over it: Tk packs the vertical pixels into the
        low 16 bits of the delta, signed.
        """
        low = event.delta & 0xFFFF
        return self._scroll_form(event, -(low - 0x10000 if low >= 0x8000 else low))

    def _scroll_form(self, event: tk.Event, pixels: int) -> str | None:
        """Scroll the form ``pixels`` when the pointer is over it, and stop the event there; leave any other alone."""
        under_pointer = self.root.winfo_containing(event.x_root, event.y_root)
        if under_pointer is None or not str(under_pointer).startswith(str(self.form_canvas)):
            return None
        if self.form_canvas.yview() != (0.0, 1.0):  # scrolled only when it is taller than the window
            self.form_canvas.yview_scroll(pixels, 'units')
        return 'break'

    def show_fields_of_forcing(self) -> None:
        """Show the fields the chosen forcing reads, and the heading of each group that shows any."""
        forcing = self.variables['forcing'].get()
        groups_shown = set()
        for field in self.fields:
            shown = field['forcings'] is None or forcing in field['forcings']
            if shown:
                groups_shown.add(field['group'])
            for widget in self.rows[field['name']]:
                widget.grid() if shown else widget.grid_remove()
        for group, band in self.headings.items():
            band.grid() if group in groups_shown else band.grid_remove()

    def show_discharge_dir(self) -> None:
        """Disable discharge_dir and its Browse button while the discharge is discarded."""
        state = 'disabled' if self.variables['discharge_to_devnull'].get() else '!disabled'
        for widget in self.rows['discharge_dir'][1:]:
            widget.state([state])

    def browse_home(self) -> None:
        picked = filedialog.askdirectory(parent=self.root, initialdir=self.home.get(), title='home directory')
        if picked:
            self.home.set(picked)

    def browse(self, field: dict[str, Any]) -> None:
        """Pick a field's path with the native dialog of its kind, opened in the home directory."""
        picked = PICKERS[field['pick']](parent=self.root, initialdir=self.home.get(), title=field['name'])
        if not picked:
            return
        if field['pick'] == 'pattern' and not picked.endswith('.zarr'):
            picked = os.path.join(picked, '*')  # every forcing file in the folder, where a zarr store is one file
        path = self.variables[field['name']]
        assert isinstance(path, tk.StringVar), f'{field["name"]} is a path, so its value is text'
        path.set(picked)

    def values_shown(self) -> dict[str, Any]:
        """The value of each field the form shows, as form_values takes them."""
        forcing = self.variables['forcing'].get()
        discarding = self.variables['discharge_to_devnull'].get()
        values = {}
        for field in self.fields:
            name = field['name']
            hidden = field['forcings'] is not None and forcing not in field['forcings']
            if hidden or (name == 'discharge_dir' and discarding):
                continue
            value = self.variables[name].get()
            values[name] = value.strip() if isinstance(value, str) else value
        return values

    def options(self) -> dict[str, Any]:
        """The options of the form's Configs, over those of the config file loaded that the form has no field for."""
        return self.kept | form_values(self.values_shown())

    def configs(self) -> Configs:
        """The Configs of the form, with the options of the config file loaded that the form has no field for."""
        return Configs(**self.options())

    def load_config(self) -> None:
        """Fill the form from a config file, keeping the options the form has no field for to route and save with it."""
        path = filedialog.askopenfilename(
            parent=self.root, initialdir=self.home.get(), title='Load config', filetypes=CONFIG_FILE_TYPES
        )
        if not path:
            return
        try:
            options = Configs.from_json(path).to_dict()
        except (ValueError, TypeError, OSError) as error:
            self._show_preview_text(f'{path}: {error}', 'ERROR')
            return
        self.form_name = Path(path).name
        self.kept = {name: value for name, value in options.items() if name not in self.variables}
        discarding = options['discharge_dir'] is not None and is_dev_null(options['discharge_dir'])
        options |= {'discharge_to_devnull': discarding, 'forcing_files': os.pathsep.join(options['forcing_files'])}
        options['discharge_dir'] = None if discarding else options['discharge_dir']
        for name, variable in self.variables.items():
            value = None if name in Configs._TIME_STEPS and options[name] == 0 else options[name]  # 0 is calculated
            if isinstance(variable, tk.BooleanVar):
                variable.set(bool(value))
            else:
                variable.set('' if value is None else str(value))

    def save_config(self) -> None:
        """Write the form's Configs to a config file, which rr route and Load config read, or show why it cannot."""
        path = filedialog.asksaveasfilename(
            parent=self.root,
            initialdir=self.home.get(),
            title='Save config',
            defaultextension='.json',
            filetypes=CONFIG_FILE_TYPES,
        )
        if not path:
            return
        try:
            self.configs().to_json(path)
        except (ValueError, TypeError, OSError) as error:
            self._show_preview_text(str(error), 'ERROR')
            return
        self.form_name = Path(path).name

    def request_preview(self) -> None:
        """Show the preview again once Tk is idle, once for all the changes made until then."""
        if not self.preview_requested:
            self.preview_requested = True
            self.root.after_idle(self.show_preview)

    def show_preview(self) -> None:
        """
        Show the form's Configs as the Python that routes it on the threads chosen, leaving out the options of the
        config file loaded that are their defaults, or as the JSON Save config writes, or show why the form makes none.
        """
        self.preview_requested = False
        chosen = self.preview_format.get()
        for name, button in self.format_buttons.items():
            button.configure(style='CurrentFormat.TButton' if name == chosen else 'Format.TButton')
        try:
            options = self.options()
            configs = Configs(**options)
            threads = threads_from(self.threads.get())
        except (ValueError, TypeError, OSError) as error:
            self._show_preview_text(str(error), 'ERROR')
            return
        if chosen == 'JSON':
            self._show_preview_text(json.dumps(configs.to_dict(), indent=2))
            return
        defaults = Configs().to_dict()
        shown = {name: value for name, value in options.items() if name not in self.kept or value != defaults[name]}
        self._show_preview_text(python_of(shown, threads))

    def _show_preview_text(self, text: str, level: str | None = None) -> None:
        """
        Show ``text`` in the preview, colored by ``level``, or when there is none as code, each line numbered and its
        parts colored by CODE_COLORS.
        """
        self.preview_text = text
        self.preview.configure(state='normal')
        self.preview.delete('1.0', 'end')
        if level is not None:
            self.preview.insert('1.0', text, (level,))
        for row, line in enumerate(text.splitlines() if level is None else [], start=1):
            number = f'{row:>4}  '
            self.preview.insert('end', number, ('line_number',))
            self.preview.insert('end', f'{line}\n')
            for match in CODE_TOKENS.finditer(line):
                kind = match.lastgroup
                assert kind is not None, 'every match is one of the named kinds of CODE_TOKENS'
                self.preview.tag_add(kind, f'{row}.{len(number) + match.start()}', f'{row}.{len(number) + match.end()}')
        self.preview.configure(state='disabled')

    def copy_preview(self) -> None:
        """Put the preview's text, without its line numbers, on the clipboard."""
        self.root.clipboard_clear()
        self.root.clipboard_append(self.preview_text)

    def add_form_to_queue(self) -> None:
        """Queue the form's Configs as a job named after the config file loaded or saved last, and show the queue."""
        try:
            options = self.configs().to_dict()
        except (ValueError, TypeError, OSError) as error:
            self._show_preview_text(str(error), 'ERROR')
            return
        self.add_job(self.form_name, options)
        self.show_page('Run Queue')

    def load_configs(self) -> None:
        """Queue config files, each a job named after its file, which rr route reads when the job runs."""
        paths = filedialog.askopenfilenames(
            parent=self.root, initialdir=self.home.get(), title='Load configs', filetypes=CONFIG_FILE_TYPES
        )
        for path in paths:
            self.add_job(Path(path).name, path)

    def add_job(self, name: str, config: str | dict[str, Any]) -> None:
        """Add a pending job of ``config``, a config file's path or a Configs' options, to the end of the queue."""
        job = Job(next(self.job_numbers), name, config)
        self.jobs.append(job)
        self.tree.insert('', 'end', iid=str(job.number))
        self._show_job_row(job)

    def chosen_job(self) -> Job | None:
        """The job chosen in the queue, or None."""
        chosen = self.tree.selection()
        return next((job for job in self.jobs if str(job.number) in chosen), None)

    def remove_job(self) -> None:
        """Remove the job chosen from the queue, unless it is running."""
        job = self.chosen_job()
        if job is None or job is self.running_job:
            return
        self.jobs.remove(job)
        self.tree.delete(str(job.number))
        self.show_job_log()

    def requeue_job(self) -> None:
        """Put the job chosen back at the end of the queue, pending, once it has finished."""
        job = self.chosen_job()
        if job is None or job.status not in FINISHED:
            return
        self.jobs.remove(job)
        self.jobs.append(job)
        job.status, job.lines, job.started, job.ended = 'pending', [], None, None
        self.tree.move(str(job.number), '', 'end')
        self._show_job_row(job)
        self.show_job_log()

    def clear_finished(self) -> None:
        """Remove every finished job from the queue."""
        finished = [job for job in self.jobs if job.status in FINISHED]
        self.jobs = [job for job in self.jobs if job.status not in FINISHED]
        for job in finished:
            self.tree.delete(str(job.number))
        self.show_job_log()

    def run_queue(self) -> None:
        """Route the queue's pending jobs one after another, from the first, until none is left or Stop is pressed."""
        if self.running_job is None:
            self.start_next_job(following=False)

    def start_next_job(self, following: bool) -> None:
        """
        Route the first pending job with rr route in a child process on the threads chosen, choosing it in the queue
        when ``following`` or when no job is chosen, or end the queue's run when no job is pending.
        """
        job = next((job for job in self.jobs if job.status == 'pending'), None)
        if job is None:
            self._set_queue_running(False)
            return
        if following or not self.tree.selection():
            self.tree.selection_set(str(job.number))
            self.tree.see(str(job.number))
        job.status, job.started, self.running_job = 'running', time.monotonic(), job
        self._show_job_row(job)
        self.show_job_log()
        try:
            threads = threads_from(self.threads.get())
        except ValueError as error:
            self._add_line(job, str(error), 'ERROR')
            self._finish(job, 'failed')
            self._set_queue_running(False)
            return
        config = job.config
        source = config if isinstance(config, str) else '-'  # a config file, or - for the Configs given on stdin
        process = subprocess.Popen(
            [sys.executable, '-m', 'river_route._cli', 'route', source, '--threads', str(threads)],
            stdin=subprocess.PIPE if source == '-' else subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=0,
        )
        if source == '-':
            assert process.stdin is not None, 'a job of the form is given its Configs on stdin'
            process.stdin.write(json.dumps(config).encode())
            process.stdin.close()
        self.process = process
        self.output = queue.Queue()  # each run has its own, so a run's last events never reach the next
        threading.Thread(target=read_output, args=(process, self.output), daemon=True).start()
        self._set_queue_running(True)

    def stop_queue(self) -> None:
        """End the running job, and the queue's run with it, leaving the jobs after it pending."""
        self.queue_running = False
        if self.process is not None:
            self.stop_requested = True
            self.process.terminate()

    def show_output(self) -> None:
        """
        Show what the running job has printed since the last call, end it once it has exited, show how long it has
        routed, and call again in POLL_MILLISECONDS.
        """
        output = self.output  # the running job's, which a job started when it ends replaces
        for _ in range(output.qsize()):
            kind, text = output.get_nowait()
            job = self.running_job
            assert job is not None, 'only a running job prints'
            if kind == 'end':
                self._end_job(job, int(text))
                break
            if kind == 'replace':
                self._replace_line(job, text)
            else:
                self._add_line(job, text)
        if self.running_job is not None:
            self._show_job_row(self.running_job)
        self.root.after(POLL_MILLISECONDS, self.show_output)

    def _end_job(self, job: Job, code: int) -> None:
        """
        Finish the running job by its exit code, or as stopped when Stop ended it, and start the next pending job while
        the queue runs, following it in the queue when the job that ended was the one chosen.
        """
        following = job is self.chosen_job()
        self._finish(job, 'complete' if code == 0 else 'stopped' if self.stop_requested else 'failed')
        if self.queue_running:
            self.start_next_job(following)
        else:
            self._set_queue_running(False)

    def _finish(self, job: Job, status: str) -> None:
        job.status, job.ended = status, time.monotonic()
        self._show_job_row(job)
        self.running_job, self.process, self.stop_requested = None, None, False

    def _set_queue_running(self, running: bool) -> None:
        self.queue_running = running
        self.run_button.state(['disabled' if running else '!disabled'])
        self.stop_button.state(['!disabled' if running else 'disabled'])

    def _show_job_row(self, job: Job) -> None:
        """Show a job's number, name, status, and how long it has routed in its row of the queue, colored by status."""
        values = (job.number, job.name, job.status.capitalize(), job.elapsed())
        self.tree.item(str(job.number), values=values, tags=(job.status,))

    def show_job_log(self) -> None:
        """Show the log of the job chosen in the queue, or nothing when none is chosen."""
        job = self.chosen_job()
        self.log.configure(state='normal')
        self.log.delete('1.0', 'end')
        for line, level in job.lines if job is not None else []:
            self.log.insert('end', f'{line}\n', (level,) if level else ())
        self.log.configure(state='disabled')
        self.log.see('end')

    def _add_line(self, job: Job, line: str, level: str | None = None) -> None:
        """Add a line to a job's log, colored by ``level`` or by the level it starts with, shown if it is chosen."""
        level = level or next((key for key in LOG_COLORS if line.startswith(key)), None)
        job.lines.append((line, level))
        if job is self.chosen_job():
            self.log.configure(state='normal')
            self.log.insert('end', f'{line}\n', (level,) if level else ())
            self.log.configure(state='disabled')
            self.log.see('end')

    def _replace_line(self, job: Job, line: str) -> None:
        """Replace the last line of a job's log with ``line``, as tqdm redraws its bar."""
        if job.lines:
            job.lines.pop()
            if job is self.chosen_job():
                self.log.configure(state='normal')
                self.log.delete('end-2l linestart', 'end-1c')
                self.log.configure(state='disabled')
        self._add_line(job, line)

    def close(self) -> None:
        self.stop_queue()
        self.root.destroy()


def name_bundle() -> bool:
    """
    Set CFBundleName, in this process's bundle information, to APPLICATION_NAME, which names the application menu of the
    macOS menu bar and its About, Hide, and Quit items, which Tk builds from it as it starts, so it is set before Tk
    starts. A process outside an app bundle holds that information in a mutable dictionary, set through the Objective-C
    runtime by ctypes; a python inside one, whose dictionary is fixed, keeps its own name. Returns whether it was set,
    and True anywhere but macOS, where there is no such name to set.
    """
    if sys.platform != 'darwin':
        return True
    try:
        ctypes.cdll.LoadLibrary(FOUNDATION)  # defines NSBundle and NSString
        runtime = ctypes.cdll.LoadLibrary(OBJECTIVE_C_RUNTIME)
        runtime.objc_getClass.restype, runtime.objc_getClass.argtypes = ctypes.c_void_p, [ctypes.c_char_p]
        runtime.sel_registerName.restype, runtime.sel_registerName.argtypes = ctypes.c_void_p, [ctypes.c_char_p]
        message = ctypes.c_void_p, ctypes.c_void_p  # each message's receiver and selector, then its arguments
        send = ctypes.CFUNCTYPE(ctypes.c_void_p, *message)(('objc_msgSend', runtime))
        send_text = ctypes.CFUNCTYPE(ctypes.c_void_p, *message, ctypes.c_char_p)(('objc_msgSend', runtime))
        ask = ctypes.CFUNCTYPE(ctypes.c_bool, *message, ctypes.c_void_p)(('objc_msgSend', runtime))
        send_pair = ctypes.CFUNCTYPE(None, *message, ctypes.c_void_p, ctypes.c_void_p)(('objc_msgSend', runtime))
        selector = runtime.sel_registerName

        def string(text: str) -> int:  # a retained NSString, so that no autorelease pool is needed
            allocated = send(runtime.objc_getClass(b'NSString'), selector(b'alloc'))
            return send_text(allocated, selector(b'initWithUTF8String:'), text.encode())

        information = send(
            send(runtime.objc_getClass(b'NSBundle'), selector(b'mainBundle')), selector(b'infoDictionary')
        )
        set_object = selector(b'setObject:forKey:')
        if not ask(information, selector(b'respondsToSelector:'), set_object):  # fixed, inside an app bundle
            return False
        send_pair(information, set_object, string(APPLICATION_NAME), string('CFBundleName'))
        return True
    except OSError, AttributeError:  # a macOS without the runtime calls
        return False


def name_application() -> bool:
    """
    Register this process with macOS as APPLICATION_NAME, the name the Dock shows when the pointer is over it, rather
    than python, the name of the executable that macOS gives a process outside an app bundle. Only LaunchServices'
    private _LSSetApplicationInformationItem renames a process, the call Chromium names its processes with, so it is
    made through ctypes, once Tk has registered the window. Returns whether macOS took the name, and True anywhere else,
    where there is no such name to set.
    """
    if sys.platform != 'darwin':
        return True
    try:
        services = ctypes.cdll.LoadLibrary(CORE_SERVICES)
        create_string = ctypes.cdll.LoadLibrary(CORE_FOUNDATION).CFStringCreateWithCString
        create_string.restype = ctypes.c_void_p
        create_string.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_uint32]
        current_application = services._LSGetCurrentApplicationASN
        current_application.restype = ctypes.c_void_p
        set_information = services._LSSetApplicationInformationItem
        set_information.restype = ctypes.c_int32
        set_information.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
        name = create_string(None, APPLICATION_NAME.encode(), CF_STRING_ENCODING_UTF8)
        display_name_key = ctypes.c_void_p.in_dll(services, '_kLSDisplayNameKey')
        return set_information(LS_DEFAULT_SESSION, current_application(), display_name_key, name, None) == 0
    except OSError, AttributeError, ValueError:  # a macOS without the private call, which the window reports
        return False


def open_window(home: str | None = None, threads: int = 1) -> None:
    """
    Open the window, its file dialogs opening in ``home``, or in the working directory when it is None, its runs
    routing on ``threads`` threads until it is changed, and return once it is closed.

    Raises:
        NotADirectoryError: if ``home`` is not a directory
        ValueError: if threads is not 1 or more
    """
    home = os.path.abspath(os.path.expanduser(home)) if home is not None else os.getcwd()
    if not os.path.isdir(home):
        raise NotADirectoryError(f'home directory not found: {home}')
    threads_from(str(threads))
    name_bundle()  # before Tk starts, which names the menu bar's application menu from it
    root = tk.Tk()
    RouteWindow(root, home, threads)
    root.mainloop()
