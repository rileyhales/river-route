import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

__all__ = ['LIBRARY', 'build']

HERE = Path(__file__).parent
SOURCES = (HERE / 'c' / 'routing.c', HERE / 'c' / 'checks.c')  # which include routing.h, checks.h, and simd.h
LIBRARY = HERE / 'build' / ('routing.dylib' if sys.platform == 'darwin' else 'routing.so')
EXPORTS = ('rr_route_rivers', 'rr_describe_problem')
NATIVE = {'arm64': '-mcpu=native', 'aarch64': '-mcpu=native'}.get(platform.machine().lower(), '-march=native')
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


def build(llvm_bin: str | None = None, checks: bool = False, wasm: str | None = None) -> None:
    """
    Compile the C kernel into LIBRARY, or into the WebAssembly modules of WASM_MODULES in the folder wasm, printing each
    command it runs.

    Args:
        llvm_bin: the folder of LLVM's clang (and clang-tidy, for the checks), or None for LLVM_BIN in the environment,
            or else the clang on the PATH
        checks: also compile with every warning of clang as an error, and run clang-tidy and clang's static analyzer
        wasm: the folder the WebAssembly modules are written into, such as the npm package's wasm/, or None to compile
            LIBRARY

    Raises:
        FileNotFoundError: if no clang is found, or no clang-tidy beside it for the checks
        NotADirectoryError: if wasm is not a folder
        subprocess.CalledProcessError: if a compile or a check fails
    """
    folder = llvm_bin or os.environ.get('LLVM_BIN')
    clang = Path(folder) / 'clang' if folder else shutil.which('clang')
    if clang is None or not Path(clang).exists():
        raise FileNotFoundError(
            "the C kernel is compiled with LLVM's clang: install LLVM, or give its folder as LLVM_BIN"
        )
    clang = Path(clang)
    tidy = clang.with_name('clang-tidy')
    if checks and not tidy.exists():
        raise FileNotFoundError(f'the checks run clang-tidy, and there is none at {tidy}')
    if wasm is not None and not Path(wasm).is_dir():
        raise NotADirectoryError(f'the WebAssembly modules are written into a folder, and {wasm} is not one')

    warnings = STRICT_WARNINGS if checks else ('-Wall', '-Wextra')
    commands: list[tuple[str | Path, ...]] = []
    if checks:
        # --header-filter=.* checks the headers too
        tidy_options = ('--quiet', '--header-filter=.*', '--warnings-as-errors=*', f'--checks={TIDY_CHECKS}')
        commands.append((tidy, *SOURCES, *tidy_options, '--', '-std=c17'))
        for source in SOURCES:  # one at a time: the analyzer writes one output, here none, for each
            commands.append((clang, '--analyze', '-std=c17', '-Xclang', '-analyzer-werror', '-o', os.devnull, source))
    if wasm is None:
        LIBRARY.parent.mkdir(exist_ok=True)
        commands.append((clang, '-std=c17', '-O3', NATIVE, '-shared', '-fPIC', *warnings, '-o', LIBRARY, *SOURCES))
    else:
        for name, features in WASM_MODULES.items():
            wasm_flags = (*WASM_TARGET, *features, *warnings, *WASM_LINK)
            commands.append((clang, '-std=c17', '-O3', *wasm_flags, '-o', Path(wasm) / name, *SOURCES))

    for command in commands:
        print(' '.join(str(part) for part in command), flush=True)
        subprocess.run([str(part) for part in command], check=True)
    print(f'Compiled the C routing kernel into {LIBRARY if wasm is None else wasm}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='build.py', description='Compile the C routing kernel with LLVM')
    parser.add_argument('--llvm-bin', type=str, default=None, help='Folder of LLVM clang (default: LLVM_BIN, or PATH)')
    parser.add_argument('--checks', action='store_true', help='Also require clang-tidy, and its static analyzer')
    parser.add_argument('--wasm', type=str, default=None, metavar='FOLDER', help='Compile the WASM into FOLDER instead')
    args = parser.parse_args()
    try:
        build(args.llvm_bin, args.checks, args.wasm)
    except (FileNotFoundError, NotADirectoryError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
