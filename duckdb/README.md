# DuckDB output

Optional typed DuckDB output for Luna, and `destrat`, a reader for both
SQLite (`-o`) and DuckDB (`-d`) result files.

`destrat` is the default reader and requires the DuckDB SDK to build. The
legacy SQLite-only reader is built as `destrat-sqlite` from `utils/reader2.cpp`
and does not require DuckDB. `make destrat-duckdb` remains available as a
compatibility target for the new reader. `DUCKDB=1` enables DuckDB output in
Luna; it is not required to build the new reader.

```sh
make -j8 DUCKDB=1 DUCKDB_PATH=/tmp/duckdb-2.0.1 luna destrat
./luna sample.lst -d results.duckdb -s 'STATS epoch'
./destrat results.duckdb +STATS -r CH E > stats.tsv
./destrat results.duckdb +STATS -r CH E --parquet stats.parquet
```

`DUCKDB_LINK=static` is the default. Set `DUCKDB_PATH` to a native SDK with
`lib/libduckdb_static.a` and `include/duckdb.h` (a flat SDK layout also works).
The archive must include the extensions you need; the build script below
bundles core functions, Parquet and JSON, with an explicit extension registry.
There is no fallback to shared linking if the static archive is missing.
For Windows builds use MSYS2 UCRT64, `ARCH=WINDOWS`, and a matching GCC/MinGW-w64 SDK, not MSVC libraries. `MINGW_PREFIX` follows the active MSYS2 environment (default `/ucrt64`).

For an existing shared-library SDK, opt in explicitly:

```sh
make DUCKDB=1 DUCKDB_LINK=shared DUCKDB_PATH=/path/to/sdk luna destrat
```

Shared Windows builds require `duckdb.dll` alongside the executables. Static
builds do not. `make luna` without `DUCKDB=1` builds Luna without DuckDB output;
the new `destrat` still needs its SDK. `make destrat-sqlite` needs no DuckDB SDK.

## Download an SDK for local builds

From the Luna checkout, run:

```sh
python3 scripts/download_duckdb.py
make -j8 ARCH=MAC DUCKDB=1 DUCKDB_PATH="$PWD/duckdb-sdk" FFTW=/opt/homebrew LGBM_PATH=/opt/homebrew
```

The installer puts `duckdb-sdk/` in **your current working directory**, regardless
of where the script itself lives. `--prefix` selects another location. It caches
verified downloads in `.duckdb-downloads/2.0.1/` beside the SDK and does not install
anything under your home directory or system library directories.

It downloads official, SHA-256-pinned DuckDB 2.0.1 static archives for the native
architecture, combines the engine/core_functions/Parquet/JSON archives, and
compiles a small C registration shim. It does not compile DuckDB itself and
needs no CMake or Git. Python 3, curl and the normal native C/C++ build tools are
required. The installer tests an actual native link and JSON/Parquet round trip
before installing. A matching verified SDK is reused; unrelated or incomplete
existing directories are never overwritten.

Supported targets are macOS Intel/ARM64, Linux x86-64/ARM64 (glibc or musl), and
Windows x86-64 in MSYS2 UCRT64. DuckDB's upstream archive is still named
`windows-mingw`; that is an asset name, not a request to use the deprecated
MINGW64 environment. The Windows archive is not compatible with MSVC. Use the matching `ARCH` and Make tool for your platform. Official prebuilt
Linux libraries still require a compatible system C/C++ runtime.

`--with-shared` additionally downloads the macOS/Linux shared library for the
ctypes parity tests. Production builds still link the static archive. Request
this on the initial install, or use a new prefix to add it later.

The downloader records upstream provenance and third-party license source links.
CI continues using the source builder below, including its aggregated license
notices and controlled build flags; the downloader does not change CI caches.

## Native dependency build and CI

`scripts/ci/build_duckdb.py` builds DuckDB 2.0.1 at the pinned commit
`8d78114e90a53823803660bf821f6073ac87da8c`. It uses CMake for the dependency
only; the normal Luna build continues to use Make. Git, Python, CMake and a
native C/C++ toolchain are needed for this dependency-build step:

```sh
python3 scripts/ci/build_duckdb.py --prefix .ci-deps/duckdb --source .ci-deps/duckdb-source --build-dir .ci-deps/duckdb-build
make DUCKDB=1 DUCKDB_PATH=.ci-deps/duckdb luna destrat destrat-sqlite
```

On macOS add `ARCH=MAC` to Make. With MSYS2 UCRT64 use
`--generator "MinGW Makefiles"` for the dependency, then `mingw32-make
ARCH=WINDOWS ...` for Luna. The dependency build is PIC-enabled for future
Python-extension reuse; it does not require Luna's unrelated `STATIC=1` option.

The SDK includes a combined static archive, headers, license, provenance
manifest and a matching shared library used only by the Python parity tests.
The engine has no external-extension loader or HTTP client: core functions,
Parquet and JSON work offline, without extension downloads.

GitHub CI builds and tests on Linux x86-64, macOS Intel/ARM64 and Windows
UCRT64. SDK caches are keyed by platform, architecture, compiler identity and
dependency-script revision. Development pushes to `luna-2.0-dev` run CI too.
Every job verifies the embedded version, checks the executables have no
DuckDB shared-library dependency, runs parity tests, packages both readers,
and repeats checks on the extracted archive with a fresh home directory.
macOS and Linux CI explicitly link static LightGBM built with OpenMP disabled
and PIC enabled. They also build static FFTW. This removes the LightGBM,
OpenMP and FFTW runtime-library requirements from those downloads; the normal
OS C/C++ runtime remains required. Native cache keys distinguish this setup
from older shared-library builds. Windows retains its existing DLL packaging,
but now builds with MSYS2 UCRT64 and `mingw-w64-ucrt-x86_64-*` packages.
LightGBM/FFTW caches live in `.ci-deps/windows-native` and include UCRT64 and
compiler identity in their keys; DuckDB's cache also records the MSYS2 environment.
The first migrated CI run rebuilds the dependencies rather than reusing MINGW64
artifacts. The CMake generator is still named `MinGW Makefiles`, as required by
the GCC toolchain; this does not select the deprecated MINGW64 environment.
This Windows GCC migration does not change luna-api's separate MSVC build.

For local macOS/Linux builds, `LGBM_LINK=shared` remains the default so existing
Homebrew installations work. To use a static LightGBM SDK, pass
`LGBM_LINK=static LGBM_PATH=/path/to/LightGBM`; the archive must be at
`lib/lib_lightgbm.a` or `lib_lightgbm.a`. The archive is linked explicitly even
if a shared library is present alongside it; a missing archive is an error.
Link-mode changes force relinking. CI builds the archive with
`BUILD_STATIC_LIB=ON`, `USE_OPENMP=OFF`, and
`CMAKE_POSITION_INDEPENDENT_CODE=ON`, and tests the native LightGBM wrapper.
`check_duckdb.py --static-lightgbm` rejects unexpected shared LightGBM/OpenMP
or FFTW dependencies before testing binaries or extracted archives.

macOS archives bundle non-system dylibs (including shared LightGBM and OpenMP
when needed) in `lib/`, rewrite references relative to each binary, and ad-hoc
sign the edited Mach-O files. Archive checks reject external runtime dependencies
before launching, so libraries installed on the runner cannot hide omissions.
A relocation regression test exercises a transitive dylib chain after its build
directory has been moved out of reach.

Archives are uploaded on development builds as well as releases. The
`check_duckdb.py` helper can also be run locally:

```sh
python3 scripts/ci/check_duckdb.py --library .ci-deps/duckdb/lib/libduckdb.so
```

Use `libduckdb.dylib` on macOS or `duckdb.dll` on Windows. `DUCKDB=1` also enables the native `libluna` result backends and `lunapi_t`
import/export methods; see [native API documentation](../lunapi/RESULTS-IO.md).
Python bindings and wheel integration remain separate.

## Output options

| option | file absent | existing SQLite | existing DuckDB |
|---|---|---|---|
| `-o FILE` | create SQLite | replace | replace with SQLite |
| `-d FILE` | create DuckDB | replace with DuckDB | replace |
| `-a FILE` | create SQLite | append | append (needs `DUCKDB=1`) |
| `-d FILE --append` | create DuckDB | error | append |

Formats are detected from the file header, not the extension. Conflicting
output options are rejected before any file is touched.

## Schema (version 1)

One table per command and factor set, as in `-t` text tables: factors and
variables are columns, factor levels are row values. Table names are `CMD` for
baseline output and `CMD__F1_F2` otherwise (factors sorted; tags included).

* bookkeeping columns: `_seq` (global emission order, continued on append),
  `_indiv`, `_exec` (refs into the catalog)
* factor columns: typed from Luna's factor declarations (`writer_t::set_types`):
  `E` and other integer factors BIGINT, real factors DOUBLE, labels VARCHAR;
  `T` is stored as `T_START`/`T_STOP` UBIGINT (time-points)
* variable columns: types from cmddefs (integer BIGINT, real DOUBLE, text
  VARCHAR, boolean BOOLEAN); a missing value or a variable not emitted is NULL

The catalog makes the file self-describing:

| table | contents |
|---|---|
| `luna_schema` | format `luna-typed-results`, version, Luna version, creation time |
| `luna_individuals` | `indiv_ref`, ID, source file; repeated IDs on append get new refs |
| `luna_executions` | `exec_ref`, individual, command, number, parameters, time |
| `luna_factors` | factor names in declaration order, with types |
| `luna_tables` | SQL name, command, factor set, cmddefs table identity |
| `luna_columns` | each column's SQL name, logical name, role and type |

There are no uniqueness constraints or indexes. A cell written twice for the
same row keeps the last value. New cmddefs variables are added to existing
tables on append; a changed type is an error. Undefined variables, unparseable
numeric levels, conversion failures and SQL name clashes (case-insensitive) halt
Luna, since they need fixes in cmddefs or the emitting command.

The earlier experimental DuckDB schema (SQLite-shaped `datapoints` tables) is
recognised and rejected with a request to regenerate it with `-d`.

## Reader

`destrat` takes the same options as the legacy `destrat-sqlite` and accepts any mix of
SQLite and DuckDB inputs. It shares one front end (strata selection, wide
merging, `-l`, `-t`, summaries) between two adapters:

* SQLite: the stored strata and datapoints, read as `destrat` reads them
* DuckDB: strata synthesised from each table's level combinations, numbered in
  order of first emission; only the requested columns, individuals, commands
  and level combinations are read

Values print as `destrat-sqlite` prints them from SQLite: whole-number reals as
integers, other reals with 15 significant digits. Wide and `-t` output match
`destrat-sqlite`. With `-l`, DuckDB inputs give the same values but list variables in
column order within a row, since per-cell emission order is not stored.
Summary timestamps use DuckDB's timestamp format.

`-p N` (or `-pN`) formats only floating-point values to N decimal places; integer,
boolean and text values (including dates and times) are unchanged. DuckDB
uses the declared column types; SQLite uses each value's stored type, so
whole-number reals stored as integers remain integers there.

`--parquet FILE` writes wide output with native types: `ID` VARCHAR, row
factors typed as declared (`T` as its `start_stop` label), and variable columns
typed from their values (BIGINT, DOUBLE, BOOLEAN or VARCHAR). Values are not
rounded by `-p`. SQLite files do not record integer vs real, so numeric
factors and variables whose values are all whole numbers become BIGINT there.
`--parquet` cannot be combined with `-l`, `-t`, summaries or `-z`.

## Tests

```sh
make -j8 DUCKDB=1 DUCKDB_PATH=.ci-deps/duckdb luna destrat destrat-sqlite
python3 duckdb/test_parity.py .ci-deps/duckdb/lib/libduckdb.so
```

This checks the catalog against the physical tables, cell-by-cell parity with
the SQLite writer, reader output against `destrat-sqlite` for both formats, Parquet
types and values, append, and output option conflicts.
