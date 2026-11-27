#!/usr/bin/env python3
"""Bundle non-system dylibs and make a macOS CLI package relocatable."""
import argparse
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess


def dependencies(binary):
    output = subprocess.check_output(["otool", "-L", str(binary)], text=True)
    return [line.strip().split(" (compatibility", 1)[0]
            for line in output.splitlines()[1:]]


def rpaths(binary):
    output = subprocess.check_output(["otool", "-l", str(binary)], text=True)
    return re.findall(r"cmd LC_RPATH\s+cmdsize \d+\s+path (.*?) \(offset", output)


def system(path):
    return path.startswith(("/usr/lib/", "/System/Library/"))


def bundle(directory, search_dirs=()):
    directory = Path(directory).resolve()
    binaries = [directory / name for name in ("luna", "destrat", "destrat-sqlite", "behead")]
    lib = directory / "lib"
    sources = {}
    queue = [(binary, binary) for binary in binaries]
    processed = []

    def resolve(name, source):
        def expand(value):
            return value.replace("@loader_path", str(source.parent)).replace("@executable_path", str(directory))
        if name.startswith("@rpath/"):
            candidates = [Path(expand(p)) / name[len("@rpath/"):] for p in rpaths(source)]
            candidates += [Path(p) / Path(name).name for p in search_dirs]
        else:
            candidates = [Path(expand(name))]
        for candidate in candidates:
            if candidate.is_file():
                return candidate.resolve()
        raise RuntimeError(f"Cannot resolve {name} required by {source}; tried {candidates}")

    while queue:
        destination, source = queue.pop(0)
        if destination in processed:
            continue
        processed.append(destination)
        for name in dependencies(source):
            # otool includes LC_ID_DYLIB as the first entry for dylibs.
            if source.suffix == ".dylib" and Path(name).name == source.name:
                continue
            if system(name):
                continue
            original = resolve(name, source)
            if original.suffix != ".dylib":
                raise RuntimeError(f"Unsupported non-system library: {original}")
            target = lib / original.name
            previous = sources.get(target)
            if previous and previous != original:
                if hashlib.sha256(previous.read_bytes()).digest() != hashlib.sha256(original.read_bytes()).digest():
                    raise RuntimeError(f"Conflicting library names: {previous}, {original}")
            if not previous:
                lib.mkdir(exist_ok=True)
                sources[target] = original
                shutil.copy2(original, target)
                target.chmod(target.stat().st_mode | 0o200)
                queue.append((target, original))
            relative = os.path.relpath(target, destination.parent)
            subprocess.run(["install_name_tool", "-change", name, "@loader_path/" + relative,
                            str(destination)], check=True)
        if destination.parent == lib:
            subprocess.run(["install_name_tool", "-id", "@rpath/" + destination.name,
                            str(destination)], check=True)
        # Remove build-machine search paths so tests cannot fall back to them.
        for path in rpaths(destination):
            if path.startswith("/"):
                subprocess.run(["install_name_tool", "-delete_rpath", path, str(destination)], check=True)
    # Mach-O edits invalidate signatures, including signatures on copied dylibs.
    for binary in reversed(processed):
        subprocess.run(["codesign", "--force", "--sign", "-", str(binary)], check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory")
    parser.add_argument("--search-dir", action="append", default=[])
    args = parser.parse_args()
    bundle(args.directory, args.search_dir)
