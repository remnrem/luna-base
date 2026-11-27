"""Typed DuckDB output and reader checks. Run against a DUCKDB=1 build with
destrat and destrat-sqlite also built; pass the path to the matching DuckDB shared library:

    python3 duckdb/test_parity.py /tmp/duckdb-2.0.1/libduckdb.dylib

Checks:
  * the catalog (luna_schema/tables/columns/factors) describes the physical tables
  * every non-missing SQLite datapoint appears exactly once as a typed cell
  * destrat gives destrat-sqlite's wide output for SQLite and DuckDB inputs,
    and the same long-output values (variable order within a row may differ)
  * --parquet is typed and matches the TSV
  * -d --append continues refs and _seq; -a detects DuckDB by its header
  * conflicting output options leave existing files untouched
"""
import csv
import ctypes as c
import json
import math
import os
from pathlib import Path
import shutil
import sqlite3
import struct
import subprocess
import sys
import tempfile
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]
COMMANDS = "HEADERS & STATS epoch & ANNOTS & HYPNO"


def run(cwd, exe, *args, ok=True):
    proc = subprocess.run([ROOT / exe, *map(str, args)], cwd=cwd,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ok and proc.returncode:
        raise RuntimeError(f"{exe} {args}: {proc.stderr.decode()[-400:]}")
    return proc


def make_edf(path, epochs=6):
    # One channel; fixed dates and samples keep both runs comparable.
    def field(value, width):
        return str(value).encode("ascii").ljust(width)

    n = 100 * 30 * epochs
    header = [(0, 8), ("test", 80), ("test", 80), ("01.01.23", 8), ("00.00.00", 8),
              (512, 8), ("", 44), (30 * epochs, 8), (1, 8), (1, 4)]
    channel = [("EEG", 16), ("", 80), ("uV", 8), (-100, 8), (100, 8),
               (-32768, 8), (32767, 8), ("", 80), (100, 8), ("", 32)]
    samples = [int(12000 * math.sin(i / 13.0) + 3000 * math.sin(i / 3.1)) for i in range(n)]
    path.write_bytes(b"".join(field(*f) for f in header + channel)
                     + struct.pack(f"<{n}h", *samples))


class Duck:
    """The same embedded SDK as Luna, with no Python package or CLI dependency."""
    def __init__(self, lib, path):
        self.lib = lib
        self.db, self.con = c.c_void_p(), c.c_void_p()
        assert lib.duckdb_open(os.fsencode(str(path)), c.byref(self.db)) == 0
        assert lib.duckdb_connect(self.db, c.byref(self.con)) == 0

    def rows(self, query, tmp):
        path = tmp / "query.json"
        dest = str(path).replace("'", "''")
        sql = f"COPY ({query}) TO '{dest}' (FORMAT JSON, ARRAY TRUE)"
        if self.lib.duckdb_query(self.con, sql.encode(), None):
            raise RuntimeError(query)
        return json.loads(path.read_text())

    def close(self):
        self.lib.duckdb_disconnect(c.byref(self.con))
        self.lib.duckdb_close(c.byref(self.db))


def number(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def same_value(a, b):
    x, y = number(a), number(b)
    if x is not None and y is not None:
        return math.isclose(x, y, rel_tol=1e-12, abs_tol=1e-12)
    return str(a) == str(b)


def sqlite_cells(path):
    """(id, cmd, var, levels) -> Counter of values for non-missing datapoints."""
    db = sqlite3.connect(path)
    levels = {}
    for sid, fac, lvl in db.execute(
            "SELECT s.strata_id, f.factor_name, l.level_name FROM strata s "
            "JOIN levels l ON l.level_id = s.level_id JOIN factors f ON f.factor_id = l.factor_id"):
        if not fac.startswith("_") and fac not in ("E", "T"):
            levels.setdefault(sid, {})[fac] = lvl
    out = Counter()
    for ind, cmd, var, sid, epoch, start, stop, value in db.execute(
            "SELECT i.indiv_name, c.cmd_name, v.variable_name, d.strata_id, t.epoch, t.start, t.stop, d.value "
            "FROM datapoints d JOIN individuals i USING (indiv_id) JOIN commands c USING (cmd_id) "
            "JOIN variables v ON v.variable_id = d.variable_id LEFT JOIN timepoints t USING (timepoint_id) "
            "WHERE d.value IS NOT NULL"):
        key = dict(levels.get(sid, {}))
        if epoch is not None:
            key["E"] = str(epoch)
        if start is not None:
            key["T"] = f"{start}_{stop}"
        out[(ind, cmd, var, tuple(sorted(key.items())), str(value))] += 1
    db.close()
    return out


def typed_cells(duck, tmp):
    out = Counter()
    tables = duck.rows("SELECT table_ref, sql_name, command FROM luna_tables", tmp)
    for t in tables:
        cols = duck.rows(f"SELECT sql_name, name, role FROM luna_columns WHERE table_ref = {t['table_ref']} ORDER BY ordinal", tmp)
        rows = duck.rows(f"SELECT i.indiv_id AS _id, t.* FROM \"{t['sql_name']}\" t "
                         "JOIN luna_individuals i ON i.indiv_ref = t._indiv", tmp)
        for r in rows:
            key = {}
            for col in cols:
                v = r[col["sql_name"]]
                if col["role"] == "factor" and v is not None:
                    key[col["name"]] = str(v)
                elif col["role"] == "interval_start" and v is not None:
                    key[col["name"]] = f"{v}_{r[col['name'] + '_STOP']}"
            key = tuple(sorted(key.items()))
            for col in cols:
                if col["role"] == "variable" and r[col["sql_name"]] is not None:
                    out[(r["_id"], t["command"], col["name"], key, r[col["sql_name"]])] += 1
    return out


def match_cells(expected, actual):
    """Compare value multisets, numbers by value (SQLite stores whole reals as integers)."""
    def group(cells):
        g = {}
        for (ind, cmd, var, key, value), n in cells.items():
            g.setdefault((ind, cmd, var, key), []).extend([value] * n)
        return g
    e, a = group(expected), group(actual)
    assert e.keys() == a.keys(), sorted(set(e) ^ set(a))[:5]
    for k in e:
        ev = sorted(e[k], key=str)
        av = sorted(a[k], key=str)
        assert len(ev) == len(av), (k, ev, av)
        remaining = av.copy()
        for x in ev:
            for i, y in enumerate(remaining):
                if same_value(x, y):
                    remaining.pop(i)
                    break
            else:
                raise AssertionError((k, ev, av))


def test(tmp, lib):
    make_edf(tmp / "sample.edf")
    (tmp / "events.annot").write_text("Event\tx\t.\t1.5\t2.75\t.\nArousal\tx\t.\t40\t43.5\t.\n")
    (tmp / "stages.annot").write_text("".join(f"{s}\t.\t.\t{30*i}\t{30*(i+1)}\t.\n"
                                              for i, s in enumerate(["W", "N1", "N2", "N2", "R", "W"])))
    (tmp / "inputs.lst").write_text("".join(
        f"{subject}\t{tmp/'sample.edf'}\t{tmp/'events.annot'},{tmp/'stages.annot'}\n" for subject in ("s1", "s2")))
    run(tmp, "luna", tmp / "inputs.lst", "-o", "out.db", "-s", COMMANDS)
    run(tmp, "luna", tmp / "inputs.lst", "-d", "out.duckdb", "-s", COMMANDS)

    # catalog describes the physical tables
    duck = Duck(lib, tmp / "out.duckdb")
    schema = {r["key"]: r["value"] for r in duck.rows("SELECT * FROM luna_schema", tmp)}
    assert schema["format"] == "luna-typed-results" and schema["version"] == "1", schema
    described = duck.rows("SELECT t.sql_name AS t, c.sql_name AS c, c.sql_type AS ty FROM luna_columns c "
                          "JOIN luna_tables t USING (table_ref) ORDER BY t.sql_name, c.ordinal", tmp)
    physical = duck.rows("SELECT table_name AS t, column_name AS c, data_type AS ty FROM duckdb_columns() "
                         "WHERE schema_name = 'main' AND NOT internal AND table_name NOT LIKE 'luna_%' "
                         "ORDER BY table_name, column_index", tmp)
    assert described == physical
    assert {r["name"] for r in duck.rows("SELECT name FROM luna_factors", tmp)} >= {"E", "CH", "T"}

    # every SQLite value is a typed cell, and nothing else is
    match_cells(sqlite_cells(tmp / "out.db"), typed_cells(duck, tmp))
    duck.close()

    # Precision only formats floating-point values; dates/times and integer
    # counts retain their values. SQLite uses its stored per-value types.
    def header_row(filename, *options):
        rows = list(csv.DictReader(run(tmp, "destrat", filename, "+HEADERS", *options).stdout.decode().splitlines(),
                                   delimiter="\t", quoting=csv.QUOTE_NONE))
        return rows[0]
    for filename in ("out.db", "out.duckdb"):
        original = header_row(filename)
        formatted = header_row(filename, "-p", "4")
        assert header_row(filename, "-p4") == formatted
        for name in ("NR", "NS", "NS_ALL", "REC_DUR_HMS", "START_DATE", "START_TIME",
                     "STOP_DATE", "STOP_TIME", "TOT_DUR_HMS"):
            assert formatted[name] == original[name], (filename, name, formatted[name], original[name])
        if filename == "out.duckdb":
            for name in ("REC_DUR", "REC_DUR_SEC", "TOT_DUR_SEC"):
                assert formatted[name] == f"{float(original[name]):.4f}", (name, formatted[name])

    # reader: destrat-sqlite on SQLite == destrat on SQLite and DuckDB
    wide = [["+HEADERS"], ["+STATS", "-r", "E", "CH"], ["+STATS", "-r", "E", "-c", "CH"],
            ["+ANNOTS", "-r", "ANNOT"], ["+ANNOTS", "-r", "T", "INST", "ANNOT"], ["+HYPNO", "-r", "SS"],
            ["+STATS", "-r", "CH", "E", "-i", "s2"]]
    for args in wide:
        legacy = run(tmp, "destrat-sqlite", "out.db", *args).stdout
        assert legacy, args
        assert run(tmp, "destrat", "out.db", *args).stdout == legacy, ("sqlite", args)
        assert run(tmp, "destrat", "out.duckdb", *args).stdout == legacy, ("duckdb", args)
    for args in [["+STATS", "-r", "E", "CH", "-l"], ["+ANNOTS", "-r", "T", "INST", "ANNOT", "-l"]]:
        legacy = run(tmp, "destrat-sqlite", "out.db", *args).stdout
        assert run(tmp, "destrat", "out.db", *args).stdout == legacy, args
        def body(out):
            return sorted(line.split(b"\t", 1)[1] for line in out.splitlines())
        assert body(run(tmp, "destrat", "out.duckdb", *args).stdout) == body(legacy), args
    legacy = run(tmp, "destrat-sqlite", "out.db", "out.db", "+STATS", "-r", "E", "CH").stdout
    assert run(tmp, "destrat", "out.db", "out.duckdb", "+STATS", "-r", "E", "CH").stdout == legacy

    # typed parquet == TSV
    reader = Duck(lib, ":memory:")
    for args, types in [(["+STATS", "-r", "E", "CH"], {"ID": "VARCHAR", "E": "BIGINT", "CH": "VARCHAR", "MEAN": "DOUBLE"}),
                        (["+ANNOTS", "-r", "T", "INST", "ANNOT"], {"T": "VARCHAR", "START": "DOUBLE"})]:
        tsv = list(csv.reader(run(tmp, "destrat", "out.duckdb", *args).stdout.decode().splitlines(),
                              delimiter="\t", quoting=csv.QUOTE_NONE))
        pq = tmp / "extract.parquet"
        run(tmp, "destrat", "out.duckdb", *args, "--parquet", pq)
        where = str(pq).replace("'", "''")
        got = {r["column_name"]: r["column_type"] for r in reader.rows(f"DESCRIBE SELECT * FROM '{where}'", tmp)}
        assert all(got[k] == v for k, v in types.items()), got
        rows = reader.rows(f"SELECT * FROM '{where}'", tmp)
        assert tsv.pop(0) == list(rows[0])
        assert len(rows) == len(tsv)
        for row, line in zip(rows, tsv):
            assert all(same_value("NA" if v is None else v, t) for v, t in zip(row.values(), line)), (row, line)
    reader.close()
    assert run(tmp, "destrat", "out.duckdb", "+STATS", "-l", "--parquet", tmp / "x.parquet", ok=False).returncode

    # append: -d --append continues refs and _seq; -a detects DuckDB by its header
    shutil.copy(tmp / "out.duckdb", tmp / "app.db")
    run(tmp, "luna", tmp / "inputs.lst", "-d", "app.db", "--append", "-s", COMMANDS)
    run(tmp, "luna", tmp / "inputs.lst", "-a", "app.db", "-s", COMMANDS)
    duck = Duck(lib, tmp / "app.db")
    assert duck.rows("SELECT count(*) AS n, count(DISTINCT indiv_ref) AS r FROM luna_individuals", tmp)[0] == {"n": 6, "r": 6}
    seq = duck.rows("SELECT count(*) AS n, count(DISTINCT _seq) AS d FROM STATS__CH_E", tmp)[0]
    assert seq["n"] == seq["d"]
    duck.close()
    assert run(tmp, "luna", tmp / "inputs.lst", "-d", "out.db", "--append", "-s", "HEADERS", ok=False).returncode

    # Command errors close DuckDB before process shutdown, rolling back the
    # current individual even when earlier commands have buffered results.
    for args in [("-d", "failed.duckdb"), ("-a", "app.db")]:
        proc = run(tmp, "luna", tmp / "inputs.lst", *args, "-s", "STATS & ACF", ok=False)
        assert proc.returncode == 1, (proc.returncode, proc.stderr)
        assert proc.stderr.endswith(b"error : command requires parameter lag\n"), proc.stderr
        assert b"finishing" not in proc.stderr
        duck = Duck(lib, tmp / args[1])
        expected = 0 if args[0] == "-d" else 6
        assert duck.rows("SELECT count(*) AS n FROM luna_individuals", tmp)[0]["n"] == expected
        duck.close()

    # conflicting output options never touch existing files
    for flag in ["-o", "-a", "-t"]:
        for pair in [("-d", "sentinel", flag, "other"), (flag, "other", "-d", "sentinel")]:
            protected = tmp / "sentinel"
            protected.write_text("keep")
            proc = run(tmp, "luna", tmp / "inputs.lst", *pair, "-s", "HEADERS", ok=False)
            assert proc.returncode and protected.read_text() == "keep"
    print("catalog, cell parity, reader parity (SQLite and DuckDB), typed Parquet, append and option conflicts passed.")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as directory:
        test(Path(directory), c.CDLL(sys.argv[1]))
