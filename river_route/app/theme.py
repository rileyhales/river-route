import tkinter as tk
from tkinter import font, ttk
from typing import Literal

__all__ = ['CODE_COLORS', 'PALETTE', 'apply_theme']

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
