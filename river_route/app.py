"""A window that sets up a Configs, routes it with rr route in a child process, and shows everything it prints."""

import ctypes
import glob
import json
import os
import queue
import re
import subprocess
import sys
import threading
import tkinter as tk
from collections.abc import Callable
from dataclasses import fields
from pathlib import Path
from tkinter import filedialog, font, ttk
from typing import Any, Literal, get_args, get_origin, get_type_hints

from .configs import Configs, is_dev_null

__all__ = ['open_window']

RUNOFF_FORCINGS = ('catchment', 'grid', 'ecmwf_grib')
GRID_FORCINGS = ('grid', 'ecmwf_grib')
FORM_GROUPS = (  # each group's comment, color, and fields, as the Configs of examples/example_usage.py has them
    (None, None, ('network_file',)),
    ('primary modeling choices', 'sky', ('parameters', 'forcing', 'transform', 'dt_routing')),
    ('channel routing options', 'teal', ('start_datetime',)),  # the forcing's own options: for channel, then the rest
    (
        'the runoff forcing data and how to read it',
        'teal',
        (
            'runoff_files',
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
    ('model state files', 'violet', ('channel_state_init_file', 'channel_state_final_file')),
    ('where to place the results', 'emerald', ('discharge_dir', 'discharge_to_devnull')),
    ('Compute Options', 'amber', ('restrict_x_range', 'invalid_coefficients', 'dt_discharge', 'dt_total')),
    ('logging and validation options', 'rose', ('progress_bar', 'log_level')),
)
FORCINGS_OF_FIELD = {  # the fields only some forcings read; the form shows every other field for every forcing
    'start_datetime': ('channel',),
    'runoff_files': RUNOFF_FORCINGS,
    'grid_weights_file': GRID_FORCINGS,
    'var_grid_runoff': ('grid',),
    'var_x': ('grid',),
    'var_y': ('grid',),
    'var_t': ('grid',),
    'grid_accumulation_type': GRID_FORCINGS,
    'runoff_processing_mode': RUNOFF_FORCINGS,
    'force_positive_runoff': GRID_FORCINGS,
}
# runoff_files is text in the form, read by runoff_paths, and blank time steps are calculated from the inputs
FORM_DEFAULTS = {'runoff_files': '', 'discharge_to_devnull': False, 'dt_discharge': '', 'dt_total': ''}
PICKERS: dict[str, Callable[..., str]] = {  # the native dialog that picks each kind of path field
    'file': filedialog.askopenfilename,
    'save': filedialog.asksaveasfilename,
    'folder': filedialog.askdirectory,
    'pattern': filedialog.askdirectory,
}
CONFIG_FILE_TYPES = [('JSON', '*.json')]
PALETTE = {  # jsrr's colors, from its src/styles/tokens.css, so the two apps look alike
    'ground': '#0a0e17',
    'surface': '#111827',
    'raised': '#1e293b',
    'border': '#334155',
    'muted': '#64748b',
    'label': '#94a3b8',
    'faint': '#b6c2d3',
    'dim': '#dde5ee',
    'text': '#f8fafc',
    'accent': '#38bdf8',
    'accent_dim': '#0ea5e9',
    'on_accent': '#04121b',
    'warning': '#f5a623',
    'danger': '#fca5a5',
}
HUES = {  # the colors of each group's heading: the bar at its left, the band behind it, and its comment
    'sky': ('#38bdf8', '#0c2b40', '#e0f2fe'),
    'amber': ('#fbbf24', '#33260b', '#fef3c7'),
    'violet': ('#a78bfa', '#251b4d', '#ede9fe'),
    'emerald': ('#34d399', '#0b3026', '#d1fae5'),
    'teal': ('#2dd4bf', '#0b302d', '#ccfbf1'),
    'rose': ('#fb7185', '#3b1520', '#ffe4e6'),
}
FONT_SIZES = {'TkDefaultFont': 13, 'TkTextFont': 13, 'TkMenuFont': 13, 'TkFixedFont': 12, 'title': 24, 'heading': 16}
LOG_COLORS = {  # the color of each log line that starts with the key
    'PROGRESS': 'accent',
    'WARNING': 'warning',
    'ERROR': 'danger',
    'CRITICAL': 'danger',
    'Traceback': 'danger',
}
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
    Each field of the form in order: its name, the comment and color of its group (None for the first), the forcings it
    is shown for (every one when None), its kind, its choices, its default, and what the file dialog picks for it. The
    groups follow the Configs of examples/example_usage.py, and a field the example leaves out follows the fields of the
    group it belongs with. Choices and defaults are read from the Configs annotations, so they never go stale.
    """
    hints = get_type_hints(Configs)
    defaults = {f.name: f.default for f in fields(Configs) if f.init} | FORM_DEFAULTS
    described = []
    for group, hue, names in FORM_GROUPS:
        for name in names:
            hint, default = hints.get(name), defaults[name]
            literal = get_origin(hint) is Literal
            described.append(
                {
                    'name': name,
                    'group': group,
                    'hue': hue,
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
    blank. With ``many``, ``path`` is paths separated by os.pathsep, as runoff_files holds them, each shown so.
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


def runoff_paths(text: str) -> list[str]:
    """
    The runoff files of the form's runoff_files: paths or glob patterns separated by os.pathsep, each pattern standing
    for the files it matches in sorted order.

    Raises:
        ValueError: if a pattern matches no file
    """
    paths: list[str] = []
    for item in text.split(os.pathsep):
        pattern = os.path.expanduser(item.strip())
        matches = sorted(glob.glob(pattern)) if any(character in pattern for character in '*?[') else [pattern]
        if not matches:
            raise ValueError(f'runoff_files pattern {pattern} matches no file')
        paths.extend(path for path in matches if path)
    return paths


def form_values(form: dict[str, Any]) -> dict[str, Any]:
    """
    The Configs options of the form's values. An empty field keeps its default, a leading ~ in a path is the home
    directory, runoff_files is read by runoff_paths, and discharge_to_devnull discards the discharge.

    Raises:
        ValueError: if a whole number of seconds is not one, or a runoff_files pattern matches no file
    """
    hints = get_type_hints(Configs)
    values = {
        name: os.path.expanduser(value) if name in Configs._SINGLE_PATH_FIELDS else value
        for name, value in form.items()
        if value != ''
    }
    if 'runoff_files' in values:
        values['runoff_files'] = runoff_paths(values['runoff_files'])
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


def apply_theme(root: tk.Tk) -> tuple[font.Font, font.Font]:
    """
    Color every widget in PALETTE on ttk's clam theme, whose colors can be set on every platform, unlike the native
    themes, and size its text by FONT_SIZES: the options in muted labels, the headings larger and bold above them, and
    the title largest. Returns the fonts of the title and of the headings, which the window keeps, since a font is
    deleted once nothing refers to it.
    """
    p = PALETTE
    style = ttk.Style(root)
    style.theme_use('clam')
    for name in ('TkDefaultFont', 'TkTextFont', 'TkMenuFont', 'TkFixedFont'):
        font.nametofont(name).configure(size=FONT_SIZES[name])
    title, heading = font.nametofont('TkDefaultFont').copy(), font.nametofont('TkDefaultFont').copy()
    title.configure(size=FONT_SIZES['title'], weight='bold')
    heading.configure(size=FONT_SIZES['heading'], weight='bold')
    style.configure(
        '.',
        background=p['surface'],
        foreground=p['dim'],
        bordercolor=p['border'],
        troughcolor=p['surface'],
        fieldbackground=p['raised'],
        lightcolor=p['raised'],
        darkcolor=p['raised'],
        insertcolor=p['text'],
        selectbackground=p['accent'],
        selectforeground=p['on_accent'],
        arrowcolor=p['faint'],
        focuscolor=p['accent'],
    )
    style.configure('TLabel', foreground=p['label'])
    style.configure('Title.TLabel', foreground=p['text'], font=title)
    style.configure('TSeparator', background=p['border'])
    style.configure('Vertical.TScrollbar', background=p['raised'], bordercolor=p['surface'], gripcount=0, arrowsize=12)
    style.map('Vertical.TScrollbar', background=[('active', p['border'])])
    style_controls(style, root)
    root.configure(background=p['surface'])
    return title, heading


def style_controls(style: ttk.Style, root: tk.Tk) -> None:
    """Color the controls of the form in PALETTE: its text fields, menus, checkboxes, and buttons."""
    p = PALETTE
    style.configure('TEntry', padding=6, foreground=p['text'])
    style.map(
        'TEntry',
        fieldbackground=[('disabled', p['surface'])],
        foreground=[('disabled', p['muted'])],
        bordercolor=[('focus', p['accent'])],
        lightcolor=[('focus', p['accent'])],
    )
    style.configure('TCombobox', padding=6, foreground=p['text'], arrowsize=14)
    style.map(
        'TCombobox',
        fieldbackground=[('readonly', p['raised'])],
        foreground=[('readonly', p['text'])],
        selectbackground=[('readonly', p['raised'])],
        selectforeground=[('readonly', p['text'])],
        background=[('active', p['raised'])],
        bordercolor=[('focus', p['accent'])],
        arrowcolor=[('active', p['accent'])],
    )
    for option, color in (('background', 'raised'), ('foreground', 'text'), ('selectBackground', 'accent')):
        root.option_add(f'*TCombobox*Listbox.{option}', p[color])  # the menu a combobox drops down
    root.option_add('*TCombobox*Listbox.selectForeground', p['on_accent'])
    root.option_add('*TCombobox*Listbox.font', 'TkTextFont')
    style.configure('TSpinbox', padding=6, foreground=p['text'], arrowsize=12)
    style.map('TSpinbox', bordercolor=[('focus', p['accent'])], arrowcolor=[('active', p['accent'])])
    style.configure(
        'TCheckbutton',
        focusthickness=0,
        indicatorbackground=p['raised'],
        indicatorforeground=p['on_accent'],
        upperbordercolor=p['border'],
        lowerbordercolor=p['border'],
    )
    style.map('TCheckbutton', background=[('active', p['surface'])], indicatorbackground=[('selected', p['accent'])])
    style.configure('TButton', padding=(14, 8), focusthickness=0, background=p['raised'], foreground=p['text'])
    style.map(
        'TButton',
        background=[('disabled', p['surface']), ('active', p['border'])],
        foreground=[('disabled', p['muted'])],
    )
    style.configure('Accent.TButton', background=p['accent'], foreground=p['on_accent'], bordercolor=p['accent'])
    style.map(
        'Accent.TButton',
        background=[('disabled', p['raised']), ('active', p['accent_dim'])],
        foreground=[('disabled', p['muted'])],
        bordercolor=[('disabled', p['border'])],
    )


class RouteWindow:
    """
    The window: on its left half, the logo, a home directory that every file dialog opens in, and the form of
    describe_fields, a section for each group headed by its comment, with the fields the chosen forcing does not read
    hidden; under them Load config, Save config, the threads to route with, Run, and Stop; and on its right half, a log
    of everything the run prints, which rr route prints in a child process so that Stop can end it.
    """

    def __init__(self, root: tk.Tk, home: str, threads: int) -> None:
        self.root = root
        self.fonts = apply_theme(root)
        self.fields = describe_fields()
        self.variables: dict[str, tk.StringVar | tk.BooleanVar] = {}
        self.rows: dict[str, list[ttk.Widget]] = {}  # the widgets of each field's row of the form
        self.headings: dict[str | None, tk.Frame] = {}  # the band above each group, holding its comment
        self.kept: dict[str, Any] = {}  # the options of the config file loaded that the form has no field for
        self.home = tk.StringVar(value=home)  # the folder every file dialog opens in, and paths are shown relative to
        self.shown: dict[str, tk.StringVar] = {}  # the text of each path field, its absolute path relative to home
        self.converting = False  # whether a path field's absolute path and its text are being matched to each other
        self.threads = tk.StringVar(value=str(threads))  # a resource of each run, so never saved in a config file
        self.running = False
        self.process: subprocess.Popen[bytes] | None = None
        self.output: queue.Queue[tuple[str, str]] = queue.Queue()
        self.logo, self.icon = None, None
        if tk.TkVersion >= 9.0:  # Tk reads SVG from 9.0 on
            self.logo = tk.PhotoImage(file=LOGO, format='svg -scaletoheight 48')
            self.icon = tk.PhotoImage(file=LOGO, format='svg -scaletoheight 256')
            root.iconphoto(True, self.icon)
        root.title(APPLICATION_NAME)
        width, height = root.winfo_screenwidth(), root.winfo_screenheight()
        root.geometry(f'{width}x{int(height * 0.9)}+0+0')
        root.columnconfigure((0, 1), weight=1, uniform='halves')  # the form takes half the window and the log the other
        root.rowconfigure(0, weight=1)
        self._lay_out_side()
        self._lay_out_log()
        self.variables['forcing'].trace_add('write', lambda *_: self.show_fields_of_forcing())
        self.variables['discharge_to_devnull'].trace_add('write', lambda *_: self.show_discharge_dir())
        self.show_fields_of_forcing()
        self._select_network_file_in(home)
        self.home.trace_add('write', lambda *_: self._show_paths())
        root.bind_all('<MouseWheel>', self._scroll_by_wheel)
        if tk.TkVersion >= 9.0:  # trackpads scroll by TouchpadScroll events from Tk 9.0 on
            root.bind_all('<TouchpadScroll>', self._scroll_by_touchpad)
        root.protocol('WM_DELETE_WINDOW', self.close)
        root.after(POLL_MILLISECONDS, self.show_output)
        root.after_idle(self._name_application)  # once Tk has registered the window with the operating system

    def _name_application(self) -> None:
        if not name_application():
            self._append(f'WARNING - macOS did not take the name {APPLICATION_NAME}, so it shows python instead')

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

    def _lay_out_side(self) -> None:
        """The form, scrolling, above the buttons, on the left half of the window."""
        side = ttk.Frame(self.root)
        side.grid(row=0, column=0, sticky='nsew')
        actions = ttk.Frame(side, padding=(28, 14))
        actions.pack(side='bottom', fill='x')  # packed first, so the form scrolls above it rather than covering it
        ttk.Separator(side).pack(side='bottom', fill='x')
        ttk.Button(actions, text='Load config', command=self.load_config).pack(side='left')
        ttk.Button(actions, text='Save config', command=self.save_config).pack(side='left', padx=(8, 0))
        self.stop_button = ttk.Button(actions, text='Stop', command=self.stop_run, state='disabled')
        self.stop_button.pack(side='right')
        self.run_button = ttk.Button(actions, text='Run', command=self.start_run, style='Accent.TButton')
        self.run_button.pack(side='right', padx=(0, 8))
        threads = ttk.Spinbox(actions, from_=1, to=os.cpu_count() or 1, textvariable=self.threads, width=4)
        threads.pack(side='right', padx=(0, 18))
        ttk.Label(actions, text='threads').pack(side='right', padx=(0, 8))
        self.form_canvas = tk.Canvas(side, highlightthickness=0, background=PALETTE['surface'], yscrollincrement=1)
        scrollbar = ttk.Scrollbar(side, orient='vertical', command=self.form_canvas.yview)
        self.form_canvas.configure(yscrollcommand=scrollbar.set)
        form = ttk.Frame(self.form_canvas, padding=(28, 24))
        window = self.form_canvas.create_window((0, 0), window=form, anchor='nw')
        form.bind('<Configure>', lambda _: self.form_canvas.configure(scrollregion=self.form_canvas.bbox('all')))
        self.form_canvas.bind('<Configure>', lambda event: self.form_canvas.itemconfigure(window, width=event.width))
        scrollbar.pack(side='right', fill='y')
        self.form_canvas.pack(side='left', fill='both', expand=True)
        self._build_form(form)

    def _lay_out_log(self) -> None:
        """The log, on the right half of the window, its lines colored by level."""
        p = PALETTE
        frame = ttk.Frame(self.root)
        frame.grid(row=0, column=1, sticky='nsew')
        self.log = tk.Text(
            frame,
            wrap='word',
            font='TkFixedFont',
            state='disabled',
            padx=20,
            pady=16,
            borderwidth=0,
            highlightthickness=0,
            background=p['ground'],
            foreground=p['dim'],
            insertbackground=p['text'],
            selectbackground=p['accent'],
            selectforeground=p['on_accent'],
        )
        for level, color in LOG_COLORS.items():
            self.log.tag_configure(level, foreground=p[color])
        scrollbar = ttk.Scrollbar(frame, orient='vertical', command=self.log.yview)
        self.log.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side='right', fill='y')
        self.log.pack(side='left', fill='both', expand=True)

    def _build_form(self, form: ttk.Frame) -> None:
        """
        The logo and title, the home directory, then each field's row, with a heading above each group but the first.
        """
        form.columnconfigure(1, weight=1)
        header = ttk.Frame(form)
        header.grid(row=0, column=0, columnspan=3, sticky='w', pady=(0, 18))
        if self.logo is not None:
            ttk.Label(header, image=self.logo).pack(side='left', padx=(0, 14))
        ttk.Label(header, text=APPLICATION_NAME, style='Title.TLabel').pack(side='left')
        self._add_row(form, 1, 'home directory', ttk.Entry(form, textvariable=self.home), self.browse_home)
        grid_row, group = 2, None
        for field in self.fields:
            if field['group'] != group:
                group = field['group']
                self.headings[group] = self._add_heading(form, grid_row, group, field['hue'])
                grid_row += 1
            browse = None if field['pick'] is None else (lambda field=field: self.browse(field))
            self.rows[field['name']] = self._add_row(form, grid_row, field['name'], self._control(form, field), browse)
            grid_row += 1

    def _add_heading(self, form: ttk.Frame, grid_row: int, comment: str, hue: str) -> tk.Frame:
        """A band across the form, tinted in the group's hue with a bar of it at its left, holding its comment."""
        bar, band_color, text = HUES[hue]
        band = tk.Frame(form, background=band_color)
        band.grid(row=grid_row, column=0, columnspan=3, sticky='ew', pady=(26, 10))
        tk.Frame(band, background=bar, width=6).pack(side='left', fill='y')
        heading = tk.Label(band, text=comment, font=self.fonts[1], background=band_color, foreground=text, anchor='w')
        heading.pack(side='left', fill='x', expand=True, padx=14, pady=8)
        return band

    def _control(self, form: ttk.Frame, field: dict[str, Any]) -> ttk.Widget:
        """The widget that holds a field's value, a checkbox, a menu of its choices, or text, and its variable."""
        if field['kind'] == 'bool':
            boolean = tk.BooleanVar(value=field['default'])
            self.variables[field['name']] = boolean
            return ttk.Checkbutton(form, variable=boolean)
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
        label = ttk.Label(form, text=text)
        label.grid(row=grid_row, column=0, sticky='w', padx=(0, 16), pady=4)
        sticky = 'w' if isinstance(control, ttk.Checkbutton) else 'ew'
        control.grid(row=grid_row, column=1, sticky=sticky, pady=4)
        if browse is None:
            return [label, control]
        button = ttk.Button(form, text='Browse', command=browse)
        button.grid(row=grid_row, column=2, padx=(8, 0), pady=4)
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
            picked = os.path.join(picked, '*')  # every runoff file in the folder, where a zarr store is one file
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

    def configs(self) -> Configs:
        """The Configs of the form, with the options of the config file loaded that the form has no field for."""
        return Configs(**(self.kept | form_values(self.values_shown())))

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
            self._clear_log()
            self._append(f'{path}: {error}', 'ERROR')
            return
        self.kept = {name: value for name, value in options.items() if name not in self.variables}
        discarding = options['discharge_dir'] is not None and is_dev_null(options['discharge_dir'])
        options |= {'discharge_to_devnull': discarding, 'runoff_files': os.pathsep.join(options['runoff_files'])}
        options['discharge_dir'] = None if discarding else options['discharge_dir']
        for name, variable in self.variables.items():
            value = None if name in Configs._TIME_STEPS and options[name] == 0 else options[name]  # 0 is calculated
            if isinstance(variable, tk.BooleanVar):
                variable.set(bool(value))
            else:
                variable.set('' if value is None else str(value))

    def save_config(self) -> None:
        """Write the form's Configs to a config file, which rr route and Load config read, or log why it cannot."""
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
            self._clear_log()
            self._append(str(error), 'ERROR')

    def start_run(self) -> None:
        """
        Route the form's Configs with rr route in a child process, on a thread pool when threads is more than 1, or log
        why Configs refuses them.
        """
        if self.running:
            return
        self._clear_log()
        try:
            configs = self.configs().validate_routing()
            threads = threads_from(self.threads.get())
        except (ValueError, TypeError, OSError) as error:
            self._append(str(error), 'ERROR')
            return
        process = subprocess.Popen(
            [sys.executable, '-m', 'river_route._cli', 'route', '-', '--threads', str(threads)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=0,
        )
        assert process.stdin is not None
        process.stdin.write(json.dumps(configs.to_dict()).encode())
        process.stdin.close()
        self.process = process
        self.output = queue.Queue()  # each run has its own, so a run's last events never reach the next
        threading.Thread(target=read_output, args=(process, self.output), daemon=True).start()
        self._set_running(True)

    def stop_run(self) -> None:
        if self.running and self.process is not None:
            self.process.terminate()

    def show_output(self) -> None:
        """Show what the run has printed since the last call, and call again in POLL_MILLISECONDS."""
        for _ in range(self.output.qsize()):
            kind, text = self.output.get_nowait()
            if kind == 'end':
                self._set_running(False)
            elif kind == 'replace':
                self._replace_last_line(text)
            else:
                self._append(text)
        self.root.after(POLL_MILLISECONDS, self.show_output)

    def close(self) -> None:
        self.stop_run()
        self.root.destroy()

    def _set_running(self, running: bool) -> None:
        self.running = running
        self.run_button.state(['disabled' if running else '!disabled'])
        self.stop_button.state(['!disabled' if running else 'disabled'])

    def _append(self, line: str, level: str | None = None) -> None:
        """Add a line to the log, colored by ``level``, or by the level it starts with."""
        level = level or next((key for key in LOG_COLORS if line.startswith(key)), None)
        self.log.configure(state='normal')
        self.log.insert('end', f'{line}\n', (level,) if level else ())
        self.log.configure(state='disabled')
        self.log.see('end')

    def _replace_last_line(self, line: str) -> None:
        self.log.configure(state='normal')
        self.log.delete('end-2l linestart', 'end-1c')
        self.log.configure(state='disabled')
        self._append(line)

    def _clear_log(self) -> None:
        self.log.configure(state='normal')
        self.log.delete('1.0', 'end')
        self.log.configure(state='disabled')


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
