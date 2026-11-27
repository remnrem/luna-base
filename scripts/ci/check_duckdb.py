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


def check(binary_dir, library, home=None):
    binary_dir = Path(binary_dir).resolve()
    suffix = ".exe" if os.name == "nt" else ""
    version = subprocess.check_output([str(binary_dir / ("luna" + suffix)), "-v"], stderr=subprocess.STDOUT, text=True)
    if not re.search(r"DuckDB:\s+v?2\.0\.1(?:\s|$)", version):
        raise RuntimeError("Expected embedded DuckDB 2.0.1: " + version)
    for name in ("luna", "destrat"):
        binary = binary_dir / (name + suffix)
        command = {"Darwin": ["otool", "-L"], "Linux": ["readelf", "-d"],
                   "Windows": ["objdump", "-p"]}[platform.system()]
        deps = subprocess.check_output([*command, str(binary)], text=True)
        if re.search(r"(?:lib)?duckdb[^\s/]*\.(?:dylib|dll|so)", deps, re.IGNORECASE):
            raise RuntimeError(f"{binary} depends on a DuckDB shared library: {deps}")
    env = dict(os.environ, LUNA_TEST_BIN_DIR=str(binary_dir))
    if home:
        env.update(HOME=home, USERPROFILE=home)
    subprocess.run([sys.executable, "duckdb/test_parity.py", str(Path(library).resolve())], env=env, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-dir", default=".")
    parser.add_argument("--library", required=True)
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
            check(binaries[0].parent, args.library, str(home))
    else:
        check(args.bin_dir, args.library)
