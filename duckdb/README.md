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

Set `DUCKDB_PATH` to the unpacked SDK directory or installation prefix (the
platform SDK for Linux, macOS or Windows; for MinGW set `ARCH=WINDOWS` and ship
`duckdb.dll` alongside both executables). `make luna` selects the
default build again, without DuckDB.

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
make -j8 DUCKDB=1 DUCKDB_PATH=/tmp/duckdb-2.0.1 luna destrat destrat-sqlite
python3 duckdb/test_parity.py /tmp/duckdb-2.0.1/libduckdb.dylib
```

This checks the catalog against the physical tables, cell-by-cell parity with
the SQLite writer, reader output against `destrat-sqlite` for both formats, Parquet
types and values, append, and output option conflicts.
