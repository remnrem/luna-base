#!/usr/bin/env python3
"""Reject Linux release binaries that exceed the glibc baseline or need GCC runtimes."""
import argparse
from pathlib import Path
import re
import subprocess

BINARIES = ("luna", "destrat", "destrat-sqlite", "behead")
SYSTEM_LIBRARIES = {
    "libc.so.6", "libm.so.6", "libpthread.so.0", "libdl.so.2", "librt.so.1",
    "libresolv.so.2", "libutil.so.1", "libanl.so.1", "ld-linux-x86-64.so.2",
}


def validate(dynamic, versions, baseline=(2, 28)):
    needed = re.findall(r"\(NEEDED\).*?\[([^\]]+)\]", dynamic)
    if not needed:
        raise RuntimeError("Expected a dynamically linked glibc executable")
    unexpected = sorted(set(needed) - SYSTEM_LIBRARIES)
    if unexpected:
        raise RuntimeError("Unexpected shared dependencies: " + ", ".join(unexpected))
    if re.search(r"\((?:RPATH|RUNPATH)\)", dynamic):
        raise RuntimeError("Release executable retains a runtime library search path")
    if "GLIBC_PRIVATE" in versions:
        raise RuntimeError("Executable depends on private glibc interfaces")
    if re.search(r"(?:GLIBCXX|CXXABI)_\d", versions):
        raise RuntimeError("Executable still requires a versioned external C++ runtime")
    tags = set(re.findall(r"\bGLIBC_[A-Za-z0-9_.]+", versions))
    if any(not re.fullmatch(r"GLIBC_\d+(?:\.\d+)+", tag) for tag in tags):
        raise RuntimeError("Unsupported named glibc ABI requirement")
    required = {tuple(map(int, item.split("."))) for item in re.findall(r"\bGLIBC_(\d+(?:\.\d+)+)\b", versions)}
    if not required:
        raise RuntimeError("No glibc version requirements found")
    too_new = sorted(v for v in required if v > baseline)
    if too_new:
        raise RuntimeError("glibc requirements exceed baseline: " + ", ".join(".".join(map(str, v)) for v in too_new))
    return needed, max(required)


def check(directory, baseline=(2, 28)):
    for name in BINARIES:
        binary = Path(directory) / name
        dynamic = subprocess.check_output(["readelf", "--wide", "--dynamic", str(binary)], text=True)
        versions = subprocess.check_output(["readelf", "--wide", "--version-info", str(binary)], text=True)
        needed, highest = validate(dynamic, versions, baseline)
        print(f"{name}: maximum GLIBC_{'.'.join(map(str, highest))}; system dependencies: {', '.join(needed)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", default=".")
    args = parser.parse_args()
    check(args.directory)
