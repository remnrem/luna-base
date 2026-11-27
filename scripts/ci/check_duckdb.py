#!/usr/bin/env python3
"""Check native binaries actually embed DuckDB, then run the integration suite."""
import argparse
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import tarfile
import tempfile
import zipfile


def check_macos_package(binary_dir):
    binary_dir = Path(binary_dir).resolve()
    from bundle_macos import dependencies, system
    # Check the entire package, including transitive libraries, before
    # launching. A fresh HOME alone cannot hide runner-installed dylibs.
    for binary in [binary_dir / n for n in ("luna", "destrat", "destrat-sqlite", "behead")] + list(binary_dir.rglob("*.dylib")):
        for dependency in dependencies(binary):
            if binary.suffix == ".dylib" and dependency == "@rpath/" + binary.name:
                continue  # LC_ID_DYLIB
            if system(dependency):
                continue
            if not dependency.startswith("@loader_path/"):
                raise RuntimeError(f"Non-relocatable dependency in {binary}: {dependency}")
            resolved = (binary.parent / dependency[len("@loader_path/"):]).resolve()
            if not resolved.is_relative_to(binary_dir) or not resolved.is_file():
                raise RuntimeError(f"Missing bundled dependency in {binary}: {dependency}")


def check_static_lightgbm(binary_dir):
    system = platform.system()
    if system == "Windows":
        raise RuntimeError("Static LightGBM policy is for macOS/Linux only")
    command = ["otool", "-L"] if system == "Darwin" else ["readelf", "-d"]
    for name in ("luna", "destrat", "destrat-sqlite", "behead"):
        binary = Path(binary_dir) / name
        dependencies = subprocess.check_output([*command, str(binary)], text=True)
        # OpenMP-disabled LightGBM and FFTW must not resolve from the runner.
        if re.search(r"(?:lib_lightgbm|libomp|libgomp|libiomp|libfftw)[^\s/]*\.(?:dylib|so)", dependencies, re.IGNORECASE):
            raise RuntimeError(f"Unexpected shared native dependency in {binary}: {dependencies}")


def check(binary_dir, library, home=None, static_lightgbm=False, linux_baseline=False):
    binary_dir = Path(binary_dir).resolve()
    if linux_baseline:
        from check_linux_abi import check as check_linux_abi
        check_linux_abi(binary_dir)
    if static_lightgbm:
        check_static_lightgbm(binary_dir)
    suffix = ".exe" if os.name == "nt" else ""
    env = dict(os.environ, LUNA_TEST_BIN_DIR=str(binary_dir))
    if home:
        env.update(HOME=home, USERPROFILE=home)
        if platform.system() == "Darwin":
            check_macos_package(binary_dir)
    version = subprocess.check_output([str(binary_dir / ("luna" + suffix)), "-v"], stderr=subprocess.STDOUT, text=True, env=env)
    if not re.search(r"DuckDB:\s+v?2\.0\.1(?:\s|$)", version):
        raise RuntimeError("Expected embedded DuckDB 2.0.1: " + version)
    for name in ("luna", "destrat"):
        binary = binary_dir / (name + suffix)
        command = {"Darwin": ["otool", "-L"], "Linux": ["readelf", "-d"],
                   "Windows": ["objdump", "-p"]}[platform.system()]
        deps = subprocess.check_output([*command, str(binary)], text=True)
        if re.search(r"(?:lib)?duckdb[^\s/]*\.(?:dylib|dll|so)", deps, re.IGNORECASE):
            raise RuntimeError(f"{binary} depends on a DuckDB shared library: {deps}")
    subprocess.run([sys.executable, "duckdb/test_parity.py", str(Path(library).resolve())], env=env, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-dir", default=".")
    parser.add_argument("--library", required=True)
    parser.add_argument("--static-lightgbm", action="store_true", help="Require static LightGBM/FFTW and no OpenMP runtime")
    parser.add_argument("--linux-baseline", action="store_true", help="Require glibc <= 2.28 and static GCC/C++ runtimes")
    parser.add_argument("--archive", help="Extract and test a freshly packaged release")
    args = parser.parse_args()
    if args.archive:
        with tempfile.TemporaryDirectory() as tmp:
            extracted = Path(tmp) / "package"
            extracted.mkdir()
            if args.archive.endswith(".zip"):
                with zipfile.ZipFile(args.archive) as archive:
                    archive.extractall(extracted)
            else:
                with tarfile.open(args.archive) as archive:
                    archive.extractall(extracted, filter="data")
            name = "luna.exe" if os.name == "nt" else "luna"
            binaries = list(extracted.rglob(name))
            if len(binaries) != 1:
                raise RuntimeError("Expected exactly one Luna binary in release")
            home = Path(tmp) / "home"
            home.mkdir()
            check(binaries[0].parent, args.library, str(home), args.static_lightgbm, args.linux_baseline)
    else:
        check(args.bin_dir, args.library, static_lightgbm=args.static_lightgbm, linux_baseline=args.linux_baseline)
