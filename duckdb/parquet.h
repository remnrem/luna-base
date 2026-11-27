#pragma once
#include "connection.h"
#include "reader_source.h"
#include <vector>

namespace luna_duck {

inline std::string quote_sql(const std::string &s, char quote) {
  std::string out(1,quote);
  for (char c : s) { out += c; if (c == quote) out += c; }
  return out + quote;
}

// Wide extraction exported with native column types. Column types are
// "BIGINT", "DOUBLE", "BOOLEAN" or "VARCHAR"; each cell is converted to its
// column's type, and NULL cells stay NULL.
class parquet_table {
  connection db;
  duckdb_appender app = nullptr;
  std::vector<std::string> types;
public:
  std::string path;
  ~parquet_table() { close(); }
  void close() {
    if (app) duckdb_appender_destroy(&app);
    db.close();
  }
  void columns(const std::vector<std::string> &names, const std::vector<std::string> &column_types) {
    types = column_types;
    db.open(nullptr);
    std::string sql = "CREATE TABLE extracted(";
    for (size_t i=0; i<names.size(); ++i)
      sql += (i ? ", " : "") + quote_sql(names[i],'"') + " " + types[i];
    db.exec(sql+")");
    check(duckdb_appender_create(db.handle,nullptr,"extracted",&app),"creating Parquet export table");
  }
  void append(const std::vector<luna_reader::cell_t> &row) {
    using luna_reader::cell_t;
    for (size_t i=0; i<row.size(); ++i) {
      const cell_t &c = row[i];
      const std::string &t = types[i];
      duckdb_state status;
      if (c.kind == cell_t::null) status = duckdb_append_null(app);
      else if (t == "BIGINT") status = duckdb_append_int64(app, c.i);
      else if (t == "BOOLEAN") status = duckdb_append_bool(app, c.i != 0);
      else if (t == "DOUBLE") status = duckdb_append_double(app, c.kind == cell_t::real ? c.d : double(c.i));
      else {
        const std::string s = c.kind == cell_t::text ? c.s : luna_reader::cell_text(c);
        status = duckdb_append_varchar_length(app, s.data(), s.size());
      }
      check(status,"appending Parquet field");
    }
    check(duckdb_appender_end_row(app),"appending Parquet row");
  }
  void finish() {
    if (!app) throw std::runtime_error("no columns selected for Parquet output");
    check(duckdb_appender_close(app),"closing Parquet export table");
    db.exec("COPY extracted TO "+quote_sql(path,'\'')+" (FORMAT PARQUET)");
    close();
  }
};
}
