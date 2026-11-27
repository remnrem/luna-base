#!/usr/bin/env python3
"""Build a pinned, PIC-enabled DuckDB SDK; CMake is needed only here, not by Luna."""
import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile

VERSION = "2.0.1"
COMMIT = "8d78114e90a53823803660bf821f6073ac87da8c"
EXTENSIONS = ("core_functions", "parquet", "json")


def run(*args, **kwargs):
    subprocess.run(list(map(str, args)), check=True, **kwargs)


def build(args):
    source, output = args.source.resolve(), args.prefix.resolve()
    build_dir = args.build_dir.resolve()
    if not source.exists():
        source.parent.mkdir(parents=True, exist_ok=True)
        run("git", "clone", "--depth", "1", "--branch", "v" + VERSION,
            "https://github.com/duckdb/duckdb.git", source)
    commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], encoding="utf-8").strip()
    if commit != COMMIT:
        raise SystemExit(f"Expected DuckDB {VERSION} at {COMMIT}, found {commit}")
    configure = [args.cmake, "-S", source, "-B", build_dir,
                 "-DDUCKDB_EXPLICIT_VERSION=v" + VERSION,
                 "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
                 "-DBUILD_UNITTESTS=OFF", "-DBUILD_SHELL=OFF", "-DENABLE_JEMALLOC=OFF",
                 "-DSKIP_EXTENSIONS=icu", "-DSTATICALLY_LINK_EXTENSIONS=" + ";".join(EXTENSIONS),
                 "-DDUCKDB_CAPABILITIES=none", "-DENABLE_EXTENSION_LOAD=OFF",
                 "-DENABLE_BUILTIN_HTTPLIB=OFF", "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF"]
    if args.generator:
        configure += ["-G", args.generator]
    run(*configure)
    targets = ["duckdb_static", *[e + "_extension" for e in EXTENSIONS], "duckdb"]
    run(args.cmake, "--build", build_dir, "--config", "Release", "--parallel", args.jobs,
        "--target", *targets)
    (output / "include").mkdir(parents=True, exist_ok=True)
    (output / "lib").mkdir(exist_ok=True)
    for header in ("duckdb.h", "duckdb_static_extension.h"):
        shutil.copy2(source / "src/include" / header, output / "include" / header)
    shutil.copy2(source / "LICENSE", output / "DUCKDB-LICENSE")
    notices = []
    for path in sorted((source / "third_party").rglob("*")):
        if path.is_file() and path.name.upper().startswith(("LICENSE", "COPYING", "NOTICE")):
            notices.append(str(path.relative_to(source)) + "\n" + path.read_text(encoding="utf-8", errors="replace"))
    (output / "DUCKDB-THIRD-PARTY-NOTICES.txt").write_text("\n\n".join(notices), encoding="utf-8")
    # Explicit registration roots all extension archives, without relying on
    # platform-specific whole-archive flags or global constructors.
    with tempfile.TemporaryDirectory(prefix="duckdb-bundle-", dir=build_dir) as temp:
        work = Path(temp)
        loader = work / "static_extension_loader.c"
        run(os.sys.executable, source / "scripts/generate_static_extension_loader.py",
            "--output", loader, *EXTENSIONS)
        run(os.environ.get("CC", "cc"), "-O2", "-fPIC", "-DDUCKDB_STATIC_BUILD",
            "-I" + str(output / "include"), "-c", loader, "-o", work / "loader.o")
        archives = [build_dir / "src/libduckdb_static.a"]
        archives += [build_dir / "extension" / e / ("lib" + e + "_extension.a") for e in EXTENSIONS]
        names = []
        for i, archive in enumerate(archives):
            name = f"input{i}.a"
            shutil.copy2(archive, work / name)
            names.append(name)
        if platform.system() == "Darwin":
            run("libtool", "-static", "-o", work / "bundle.a",
                work / "loader.o", *[work / n for n in names])
        else:
            script = "CREATE bundle.a\nADDMOD loader.o\n"
            script += "".join("ADDLIB " + n + "\n" for n in names)
            script += "SAVE\nEND\n"
            run(os.environ.get("AR", "ar"), "-M", input=script, text=True, cwd=work)
            run(os.environ.get("RANLIB", "ranlib"), work / "bundle.a")
        shutil.copy2(work / "bundle.a", output / "lib/libduckdb_static.a")
    (output / "include/luna_duckdb_static.h").write_text(
        "// This SDK includes an explicit static extension registry.\n", encoding="utf-8")
    # A matching shared library is for Python ctypes parity checks only.
    shared_name = "duckdb.dll" if os.name == "nt" else (
        "libduckdb.dylib" if platform.system() == "Darwin" else "libduckdb.so")
    shared = build_dir / "src" / shared_name
    if os.name == "nt" and not shared.exists():
        shared = build_dir / "src/libduckdb.dll"  # MinGW's default DLL prefix
    shutil.copy2(shared, output / "lib" / shared_name)
    manifest = {"version": VERSION, "commit": COMMIT, "extensions": EXTENSIONS,
                "pic": True, "platform": platform.platform(),
                "compiler": subprocess.check_output([os.environ.get("CXX", "c++"), "--version"], encoding="utf-8", errors="replace")}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Staged static DuckDB {VERSION}: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--generator")
    parser.add_argument("--jobs", type=int, default=2)
    build(parser.parse_args())
