#!/usr/bin/env python3
"""Download and assemble Luna's static DuckDB SDK in ./duckdb-sdk (no DuckDB build)."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile

VERSION = "2.0.1"
COMMIT = "8d78114e90a53823803660bf821f6073ac87da8c"
EXTENSIONS = ("core_functions", "parquet", "json")
HEADERS = ("duckdb.h", "duckdb_v2.h", "duckdb_extension.h", "duckdb_extension_v2.h", "duckdb_static_extension.h")
RELEASE = f"https://github.com/duckdb/duckdb/releases/download/v{VERSION}"
# Official release digests, pinned rather than trusting an unchecked download.
STATIC = {
    "osx-amd64": "c4458122a44182f8e790d2b1e731139d7f8437ce0b14fe5a0e461612979db76b",
    "osx-arm64": "adaada35202db67d010b06277320aa700f1313d9997b549a948349c9fd749a41",
    "linux-amd64": "2d5655fd5948dfb0cd16df5fc8f7cdf0b37948374cf0629af995cf4a5e3a4c82",
    "linux-arm64": "3d581ff935f1836363cc459e640737cdc1892ab3ed69e74853e93f01d9d071f0",
    "linux-amd64-musl": "16b18b670cc6e4c8b9cc6d5716dbb7eb76ec4c1043c022fe01d4e15cae0fd601",
    "linux-arm64-musl": "d83ef54db14dd00c30aefa4d462134ae1cf09783bf8879299c34b5e38746bcb3",
    "windows-mingw": "99a432db9e0edca8ddde279747f3dbd56605492f5a5d7eca879525a4c68f772c",
}
SHARED = {
    "osx-amd64": "e12f174ef01575590ec6e27a5d475308680a9e15292691623d9848c0fc5a9ba8",
    "osx-arm64": "d71e4800350f1135bc3f80a2e5cbcc01e199ca03759edf5a376db308b514b341",
    "linux-amd64": "e25502612fed533a243b61f23a2a574daad030e7e88906222a5f0ab30d06769e",
    "linux-arm64": "811e8b4d0226f5e531f583aadc798070ffb33b350e3823be563a2ebeb2c68593",
    "linux-amd64-musl": "8c6309a332e9839e00eb766ea90e671b1e46b186a94f4aae69d3593ea16709a6",
    "linux-arm64-musl": "229e78f13ba9ca71c75a537dc950b5908b492b18d0fdc1b894d308b28ea3cc75",
}
LICENSE_SHA = "075c33400ffcb0c586dd106a029d3e733e9da3693f8fcb9cebfc651c8a2f14e3"


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def native_platform():
    system, machine = platform.system(), platform.machine().lower()
    arch = {"x86_64": "amd64", "amd64": "amd64", "arm64": "arm64", "aarch64": "arm64"}.get(machine)
    if not arch:
        raise RuntimeError(f"Unsupported architecture: {machine}")
    if system == "Darwin":
        return "osx-" + arch
    if system == "Linux":
        musl = platform.libc_ver()[0] == "musl" or bool(list(Path("/lib").glob("ld-musl-*.so.1")))
        return "linux-" + arch + ("-musl" if musl else "")
    if system == "Windows" or system.startswith(("MINGW", "MSYS")):
        if arch != "amd64" or os.environ.get("MSYSTEM") != "UCRT64":
            raise RuntimeError("Windows requires the x86-64 MSYS2 UCRT64 shell; this SDK is not for MSVC")
        return "windows-mingw"
    raise RuntimeError(f"Unsupported operating system: {system}")


def run(command, **kwargs):
    subprocess.run([str(item) for item in command], check=True, **kwargs)


def fetch(url, expected, cache, filename):
    cache.mkdir(parents=True, exist_ok=True)
    target = cache / filename
    if target.exists() and digest(target) == expected:
        print(f"Using verified download: {target}", flush=True)
        return target
    print(f"Downloading {url}", flush=True)
    # Keep partial downloads separate; never reuse an incomplete archive.
    with tempfile.TemporaryDirectory(prefix="download-", dir=cache) as temp:
        partial = Path(temp) / "download"
        run(["curl", "--fail", "--location", "--proto", "=https", "--retry", "3", "--connect-timeout", "30",
             "--max-time", "600", "--output", partial, url])
        if digest(partial) != expected:
            raise RuntimeError(f"SHA-256 mismatch for {url}; download was not installed")
        os.replace(partial, target)
    return target


def extract_files(archive, names, directory):
    # Copy only known regular files, without extracting paths or archive symlinks.
    with tarfile.open(archive) as source:
        for name in names:
            candidates = [m for m in source.getmembers() if Path(m.name).name == name and m.isfile()]
            if len(candidates) != 1:
                raise RuntimeError(f"Expected exactly one {name} in {archive}")
            with source.extractfile(candidates[0]) as incoming, open(directory / name, "wb") as outgoing:
                shutil.copyfileobj(incoming, outgoing)


def validate_sdk(sdk, work, cxx, target):
    source = work / "probe.cpp"
    source.write_text(r'''
#include "duckdb.h"
#include <cstdio>
#include <cstring>
extern "C" int duckdb_register_static_extensions(void);
int main() {
  if (duckdb_register_static_extensions() || std::strcmp(duckdb_library_version(), "v2.0.1")) return 1;
  duckdb_database db = nullptr;
  duckdb_connection con = nullptr;
  if (duckdb_open(nullptr, &db) || duckdb_connect(db, &con)) return 2;
  const char *queries[] = {
    "SELECT sum(i), json_extract('{\"v\":42}', '$.v') FROM range(3) t(i)",
    "COPY (SELECT 42 AS v) TO 'probe.parquet' (FORMAT PARQUET)",
    "SELECT v FROM read_parquet('probe.parquet')"
  };
  int failed = 0;
  for (int i=0; i<3; ++i) {
    duckdb_result result{};
    if (duckdb_query(con, queries[i], &result)) {
      std::fprintf(stderr, "%s\n", duckdb_result_error(&result)); failed = 1;
    } else if (i==2 && (duckdb_row_count(&result)!=1 || duckdb_value_int64(&result,0,0)!=42)) failed = 1;
    duckdb_destroy_result(&result);
    if (failed) break;
  }
  duckdb_disconnect(&con); duckdb_close(&db);
  return failed;
}
''', encoding="utf-8")
    binary = work / ("probe.exe" if target == "windows-mingw" else "probe")
    libs = ["-pthread"]
    if target.startswith("linux-"):
        libs.append("-ldl")
    if target == "windows-mingw":
        libs += ["-lws2_32", "-lrstrtmgr"]
    run([*cxx, "-std=c++17", "-DDUCKDB_STATIC_BUILD", "-I" + str(sdk / "include"),
         source, sdk / "lib/libduckdb_static.a", *libs, "-o", binary])
    run([binary], cwd=work)


def install(args):
    target = native_platform()
    if args.with_shared and target not in SHARED:
        raise RuntimeError("The matching shared-library option is available for macOS/Linux only")
    prefix = args.prefix.absolute()
    cache = args.cache_dir.absolute() if args.cache_dir else prefix.parent / ".duckdb-downloads" / VERSION
    if prefix.exists() or prefix.is_symlink():
        manifest_path = prefix / "manifest.json"
        if manifest_path.is_file():
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if not isinstance(manifest, dict):
                raise RuntimeError(f"Unrecognized SDK manifest in {prefix}")
            files = manifest.get("sha256", {})
            required = {"include/" + name for name in HEADERS} | {"include/luna_duckdb_static.h", "lib/libduckdb_static.a", "DUCKDB-LICENSE"}
            if (manifest.get("version") == VERSION and manifest.get("platform") == target
                    and isinstance(files, dict) and required.issubset(files) and all((prefix / name).is_file() and digest(prefix / name) == value for name, value in files.items())
                    and (target != "windows-mingw" or manifest.get("msystem") == os.environ.get("MSYSTEM"))
                    and (not args.with_shared or manifest.get("with_shared"))):
                print(f"Verified SDK already installed: {prefix}")
                return prefix
        raise RuntimeError(f"Refusing to overwrite {prefix}; choose a different --prefix or move it aside")
    cc, cxx = shlex.split(os.environ.get("CC", "cc")), shlex.split(os.environ.get("CXX", "c++"))
    tools = ["curl", cc[0], cxx[0]] + (["libtool"] if target.startswith("osx-") else ["ar", "ranlib"])
    for tool in tools:
        if not shutil.which(tool):
            raise RuntimeError(f"Required tool not found: {tool}")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    filename = f"duckdb-static-libs-{target}.tar.gz"
    downloaded = fetch(f"{RELEASE}/{filename}", STATIC[target], cache, filename)
    license_file = fetch(f"https://raw.githubusercontent.com/duckdb/duckdb/{COMMIT}/LICENSE",
                         LICENSE_SHA, cache, "DUCKDB-LICENSE")
    with tempfile.TemporaryDirectory(prefix=".duckdb-sdk-", dir=prefix.parent) as temp:
        work = Path(temp)
        sdk = work / "sdk"
        (sdk / "include").mkdir(parents=True)
        (sdk / "lib").mkdir()
        names = ["libduckdb_static.a"] + [f"lib{ext}_extension.a" for ext in EXTENSIONS]
        extract_files(downloaded, names, work)
        extract_files(downloaded, HEADERS, sdk / "include")
        loader = work / "loader.c"
        code = '#include "duckdb_static_extension.h"\n'
        code += "\n".join(f"int32_t duckdb_extension_{ext}_describe(duckdb_extension_descriptor *);" for ext in EXTENSIONS)
        code += "\n#if defined(__GNUC__) || defined(__clang__)\n__attribute__((visibility(\"hidden\")))\n#endif\n"
        code += "int32_t duckdb_register_static_extensions(void) {\n"
        code += "\n".join(f"  if (duckdb_register_static_extension(duckdb_extension_{ext}_describe)) return 1;" for ext in EXTENSIONS)
        code += "\n  return 0;\n}\n"
        loader.write_text(code, encoding="utf-8")
        run([*cc, "-O2", "-fPIC", "-DDUCKDB_STATIC_BUILD", "-I" + str(sdk / "include"),
             "-c", loader, "-o", work / "loader.o"])
        output = sdk / "lib/libduckdb_static.a"
        if target.startswith("osx-"):
            run(["/usr/bin/libtool", "-static", "-o", output, work / "loader.o", *[work / n for n in names]])
        else:
            script = "CREATE bundle.a\nADDMOD loader.o\n" + "".join(f"ADDLIB {name}\n" for name in names) + "SAVE\nEND\n"
            run(["ar", "-M"], input=script, text=True, cwd=work)
            run(["ranlib", work / "bundle.a"])
            shutil.move(work / "bundle.a", output)
        (sdk / "include/luna_duckdb_static.h").write_text("// Explicit core_functions/Parquet/JSON registration.\n", encoding="utf-8")
        shutil.copy2(loader, sdk / "static_extension_loader.c")
        shutil.copy2(license_file, sdk / "DUCKDB-LICENSE")
        if args.with_shared:
            name = f"duckdb-shared-libs-{target}.tar.gz"
            shared = fetch(f"{RELEASE}/{name}", SHARED[target], cache, name)
            extract_files(shared, ["libduckdb.dylib" if target.startswith("osx-") else "libduckdb.so"], sdk / "lib")
        print("Checking native link, embedded version, JSON and Parquet round trip...", flush=True)
        validate_sdk(sdk, work, cxx, target)
        files = {p.relative_to(sdk).as_posix(): digest(p) for p in sdk.rglob("*") if p.is_file()}
        (sdk / "manifest.json").write_text(json.dumps({
            "version": VERSION, "commit": COMMIT, "platform": target, "extensions": EXTENSIONS,
            "source": f"{RELEASE}/{filename}", "source_sha256": STATIC[target],
            "with_shared": args.with_shared, "sha256": files, "msystem": os.environ.get("MSYSTEM"),
            "third_party_licenses": f"https://github.com/duckdb/duckdb/tree/{COMMIT}/third_party",
        }, indent=2) + "\n", encoding="utf-8")
        sdk.rename(prefix)
    return prefix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, default=Path("duckdb-sdk"), help="SDK directory relative to your current directory (default: ./duckdb-sdk)")
    parser.add_argument("--cache-dir", type=Path, help="Verified archive cache (default: .duckdb-downloads/2.0.1 beside the SDK)")
    parser.add_argument("--with-shared", action="store_true", help="Also download the macOS/Linux shared library for ctypes parity tests")
    args = parser.parse_args()
    try:
        installed = install(args)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError, tarfile.TarError) as error:
        parser.exit(1, f"DuckDB SDK installation failed: {error}\n")
    print(f"SDK ready: {installed}\nBuild with DUCKDB=1 DUCKDB_PATH={shlex.quote(str(installed))}")


if __name__ == "__main__":
    main()
