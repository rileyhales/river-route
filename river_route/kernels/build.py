"""
Compile the C kernel, SOURCES, with the clang of LLVM on this machine: into the shared library the kernel is loaded
from, for this machine's processor, or into the WebAssembly modules of the npm package in wasm/. ``rr compile`` runs
it. With ``checks``, the C must also compile without one warning of clang's strictest set and pass clang-tidy and
clang's static analyzer, the checks it is developed under.
"""

import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

__all__ = ['LIBRARY', 'build', 'build_wasm']

HERE = Path(__file__).parent
SOURCES = (HERE / 'c' / 'routing.c', HERE / 'c' / 'checks.c')  # which include routing.h, checks.h, and simd.h
LIBRARY = HERE / 'build' / ('routing.dylib' if sys.platform == 'darwin' else 'routing.so')
EXPORTS = ('rr_route_rivers', 'rr_describe_problem')
# the processor the library is compiled for: this machine's, since it is compiled on the machine that routes with it
NATIVE = {'arm64': '-mcpu=native', 'aarch64': '-mcpu=native'}.get(platform.machine().lower(), '-march=native')
# each module's features beyond WASM_TARGET: relaxed SIMD's fused multiply-adds, or none for runtimes without them
WASM_MODULES = {'routing-relaxed.wasm': ('-mrelaxed-simd',), 'routing.wasm': ()}
WASM_TARGET = (
    '--target=wasm32',  # 32 bit WebAssembly, which every browser runs
    '-msimd128',  # 128 bit vectors, which every browser has run since 2023
    '-mbulk-memory',  # lets clang fill memory with one instruction rather than call a memset there is no library for
    '-ffreestanding',  # no C library: the kernel calls none
    '-nostdlib',
)
WASM_LINK = (
    '-Wl,--no-entry',  # a library of functions, with no main
    *(f'-Wl,--export={name}' for name in EXPORTS),
    '-Wl,--export=__heap_base',  # where free memory starts, after the stack, for a caller to lay arrays out from
    '-Wl,--max-memory=4294967296',  # all 4 GB a 32 bit module can address
)
STRICT_WARNINGS = (
    '-Weverything',  # every warning clang has, each an error
    '-Werror',
    '-Wno-declaration-after-statement',  # C89's rule; NASA's rule 6 declares each variable in its smallest scope
    '-Wno-unsafe-buffer-usage',  # safe buffers ask for C++'s std::span in place of every pointer, which C cannot have
    '-Wno-padded',  # the interface's structs are laid out for clarity, not size
)
TIDY_CHECKS = ','.join(
    (
        '-*',
        'bugprone-*',
        'cert-*',
        'clang-analyzer-*',
        'misc-*',
        'performance-*',
        'portability-*',
        'readability-*',
        '-misc-use-internal-linkage',  # the exported functions are external so that Python and JavaScript can call them
        '-bugprone-easily-swappable-parameters',  # routines take steps and counts side by side, each checked first
        '-readability-identifier-length',  # the Muskingum coefficients keep their names, c1 to c4
        '-readability-function-cognitive-complexity',  # the checks of a call read in one place
        '-readability-magic-numbers',  # the lane numbers of a shuffle and the eight steps at a time
        '-misc-include-cleaner',
        '-performance-enum-size',  # the problems of routing.h are a C enum, an int in every language that calls it
        '-readability-trailing-comma',  # initializers end without a trailing comma, as the Python's calls do
    )
)


def build(llvm_bin: str | Path | None = None, checks: bool = False) -> Path:
    """
    Compile the C kernel into LIBRARY and return its path.

    Args:
        llvm_bin: the folder of LLVM's clang (and clang-tidy, for the checks), or None for LLVM_BIN in the environment,
            or else the clang on the PATH
        checks: also compile with every warning of clang as an error, and run clang-tidy and clang's static analyzer

    Raises:
        FileNotFoundError: if no clang is found
        subprocess.CalledProcessError: if a compile or a check fails
    """
    clang = _find_clang(llvm_bin)
    if checks:
        _check(clang)
    warnings = STRICT_WARNINGS if checks else ('-Wall', '-Wextra')
    LIBRARY.parent.mkdir(exist_ok=True)
    _run(clang, '-std=c17', '-O3', NATIVE, '-shared', '-fPIC', *warnings, '-o', LIBRARY, *SOURCES)
    return LIBRARY


def build_wasm(folder: str | Path, llvm_bin: str | Path | None = None, checks: bool = False) -> list[Path]:
    """
    Compile the C kernel into the WebAssembly modules of WASM_MODULES in ``folder`` and return their paths. clang links
    them with LLVM's wasm-ld, which it looks for beside itself and then on the PATH.

    Args:
        folder: the folder the modules are written into, such as the npm package's wasm/
        llvm_bin: the folder of LLVM's clang (and clang-tidy, for the checks), or None for LLVM_BIN in the environment,
            or else the clang on the PATH
        checks: also compile with every warning of clang as an error, and run clang-tidy and clang's static analyzer

    Raises:
        NotADirectoryError: if folder is not a folder
        FileNotFoundError: if no clang is found
        subprocess.CalledProcessError: if a compile or a check fails
    """
    folder = Path(folder)
    if not folder.is_dir():
        raise NotADirectoryError(f'the WebAssembly modules are written into a folder, and {folder} is not one')
    clang = _find_clang(llvm_bin)
    if checks:
        _check(clang)
    warnings = STRICT_WARNINGS if checks else ('-Wall', '-Wextra')
    for name, features in WASM_MODULES.items():
        _run(clang, '-std=c17', '-O3', *WASM_TARGET, *features, *warnings, *WASM_LINK, '-o', folder / name, *SOURCES)
    return [folder / name for name in WASM_MODULES]


def _find_clang(llvm_bin: str | Path | None) -> Path:
    folder = llvm_bin or os.environ.get('LLVM_BIN')
    clang = Path(folder) / 'clang' if folder else shutil.which('clang')
    if clang is None or not Path(clang).exists():
        raise FileNotFoundError(
            "the C kernel is compiled with LLVM's clang: install LLVM, or give its folder as LLVM_BIN"
        )
    return Path(clang)


def _check(clang: Path) -> None:
    tidy = clang.with_name('clang-tidy')
    tidy_options = ('--quiet', '--header-filter=.*', '--warnings-as-errors=*', f'--checks={TIDY_CHECKS}')  # headers too
    _run(tidy, *SOURCES, *tidy_options, '--', '-std=c17')
    for source in SOURCES:  # one at a time: the analyzer writes one output, here none, for each
        _run(clang, '--analyze', '-std=c17', '-Xclang', '-analyzer-werror', '-o', os.devnull, source)


def _run(*command: str | Path) -> None:
    print(' '.join(str(part) for part in command), flush=True)
    subprocess.run([str(part) for part in command], check=True)
