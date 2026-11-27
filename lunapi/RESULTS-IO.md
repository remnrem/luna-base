# Native result I/O for lunapi

The native library supports SQLite and, when built with `DUCKDB=1`, typed
DuckDB result databases and Parquet result-table export. EDF/slice storage and
`segsrv_t` are unchanged. Python bindings are maintained in `luna-api` and have
not been changed by this integration.

## Build and link

```sh
make -j8 ARCH=MAC FFTW=/opt/homebrew DUCKDB=1 DUCKDB_PATH=/path/to/sdk libluna
make ARCH=MAC FFTW=/opt/homebrew DUCKDB=1 DUCKDB_PATH=/path/to/sdk test-lunapi-io
./luna __LUNA_TESTS__ lunapi verbose
```

Select the appropriate platform, FFTW and LightGBM paths for your environment.
`DUCKDB_LINK=static` is the default. The SDK described in `duckdb/README.md`
includes the static extension registry and Parquet support. CMake is used for
the dependency build, not the normal Luna build.

`libluna.a` now contains the backend and reusable reader objects. It does not
embed the DuckDB archive: final consumers must link that archive and platform
libraries after libluna, as the Makefile does. Shared libluna links the dependency
itself. A `DUCKDB=0` library build needs no DuckDB SDK and retains SQLite I/O.
Feature-mode changes recreate the archive, removing stale backend members.

All consumers must rebuild against matching headers and archives: the native
result-cell variant now also supports `int64_t` and `bool`. Existing small
integers retain the `int` alternative. This is a native ABI change, not a
binary-compatible replacement for an already-built Python extension.

## API

The original signatures remain unambiguous for existing bindings:

```cpp
api->output_attach("results.db");                // SQLite, replace
api->write_db(path, id, source_file, tables);      // SQLite, replace
api->import_db(path);                             // Detect input from header
api->import_db(path, selected_ids);               // Empty set means all IDs
```

New entry points:

```cpp
lunapi_t::has_duckdb();                           // Compile-time capability
lunapi_t::duckdb_version();                       // Empty string if disabled
api->output_attach_format(path, "duckdb");        // "sqlite" or "duckdb"
api->write_db_format(path, id, source_file, tables, "duckdb");
api->export_parquet(database, output, command, strata, variables, ids);
```

`output_attach_format` and `write_db_format` are deliberately distinct names:
adding overloads would make the existing pybind method-pointer expressions
ambiguous before the separate Python binding update. Both replace their output;
append is not exposed by these native API methods. Unknown formats and requests
for a disabled backend fail before replacing the destination. Output paths must
be actual file paths, rather than empty strings or `:memory:`.

`write_db_format` consumes the existing column-major `rtables_return_t` format.
Every table needs a text ID column, factor columns named by the strata key, and
equal-length columns. Column positions can vary. `BL` or `.` denotes baseline;
other keys use the existing underscore-separated factor convention. Epoch E is
replayed through `writer.epoch()`. Interval T1/T2 columns contain native
time-points, not seconds, and are replayed through `writer.interval()`.
Imported endpoints retain signed 64-bit precision; legacy in-memory tables may
also supply exact integral doubles. The native return-table format rejects
endpoints above INT64_MAX; direct DuckDB-to-Parquet export retains UBIGINT.
Integer factor levels also retain 64-bit precision. Imports reject mixed
integer/real columns that would round a large integer during promotion.

For DuckDB replay, known cmddefs and standard factor types remain authoritative.
Unknown user variables get a per-output inferred schema without modifying global cmddefs.
Integer/real mixtures promote to real only if large integers remain exact;
incompatible text/bool/numeric mixtures fail. Unknown all-missing columns default
to VARCHAR. NULL/monostate means not written/missing; a table or variable with no
emitted values may therefore be absent from imported native tables. SQLite does
not retain a boolean schema; boolean values read back as integers there.

`import_db` builds a new result cache before replacing the old one. As before,
a nonexistent path returns an empty ID vector without replacing the cache.
An explicit ID selection with no matches clears the result cache. Library SQLite
reads are read-only and do not create the CLI reader's optional index. Numeric
SQLite values retain SQLite's stored integer/real distinction. Variable type
inference is scoped to each command/factor table, so an unrelated text variable
with the same name cannot corrupt integer results.

## Parquet selection and performance

```cpp
api->output_close();
api->export_parquet("results.duckdb", "psd.parquet",
                    "PSD", "CH_E_F", {"PSD"}, {});
```

Close active file output before opening that database for import/export.
Export selects exactly one command/factor table. Factors can be specified in
any order and are normalized to sorted order. Empty variable/ID sets mean all.
Output includes ID, sorted factor columns and sorted variable columns. T intervals
use the native API's T1/T2 fields. No column pivots, level filters, multiple-input
merging or whole-database directory exports are provided by this initial API.
Unsupported selections fail rather than silently choosing another table.

For typed DuckDB inputs, export uses `COPY (SELECT ...)` directly: it projects
selected columns, filters subjects, retains declared types and all-null columns,
and avoids materializing native return tables or a Python DataFrame. Rows where
all selected values are NULL are omitted (not written). Unknown IDs yield a typed
empty extract. This is an exact-factor extraction, not destrat's general pivot
and cross-stratum merge machinery; duplicate writes remain undefined usage.

For SQLite inputs, export uses the shared source adapter and native return tables,
then stages the selected table in an in-memory DuckDB connection before COPY.
This path is materialized, not bounded-memory streaming. Types are inferred from
emitted values; columns with no emitted values are absent. Empty selections
without an emitted table fail. Parquet values are not rounded for display.

Parquet output must differ from the input database, including hard-link aliases.
Parquet is an extracted result table, not a self-contained Luna database catalog.
It is not a waveform input or a replacement for EDF storage.

## Exceptions and lifecycle

File evaluation wraps subject metadata and values in a transaction. A failure
rolls back the current subject, closes file output, clears the file-output mode
and resets the native writer to safe in-memory operation. Earlier committed
subjects/results remain intact. Reattach output explicitly before further file
evaluation. Ordinary in-memory evaluation can be used immediately after failure.
This closes the failed output session rather than attempting to reuse cached
metadata from a rolled-back transaction.

In-memory and project evaluations also detach local result accumulators during
exception unwinding. Native failures propagate as C++ exceptions through the
API's existing throwing bail handler; they do not terminate the process.

## Verification

`make test-lunapi-io` builds a standalone consumer of libluna.a and checks 64-bit
round trips, Parquet readback/types/missingness, booleans, text, all-null typed
columns, selection and paths containing quotes/spaces. It also runs in a
SQLite-only build. The integrated `lunapi` test group additionally covers live
output, failed ACF rollback/recovery, replay, E/T factors and validation before
replacement. Existing `duckdb/test_parity.py` covers the CLI reader/writer.

Python bindings, wheel dependency staging, MSVC/NMake adaptation and luna-api CI
remain separate work. This native implementation does not claim those wheel
platforms have already been tested.
