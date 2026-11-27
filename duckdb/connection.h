#pragma once

#include <duckdb.h>
#include <stdexcept>
#include <string>

namespace luna_duck {

inline void check(duckdb_state status, const char *message) {
  if (status != DuckDBSuccess)
    throw std::runtime_error(message ? message : "DuckDB operation failed");
}

struct result {
  duckdb_result data{};
  result() = default;
  result(const result &) = delete;
  ~result() { duckdb_destroy_result(&data); }
};

class connection {
public:
  duckdb_database db = nullptr;
  duckdb_connection handle = nullptr;
  connection() = default;
  connection(const connection &) = delete;
  ~connection() { close(); }

  void open(const char *path, bool readonly = false) {
    close();
    duckdb_config cfg = nullptr;
    check(duckdb_create_config(&cfg), "creating DuckDB configuration");
    char *message = nullptr;
    auto status = readonly ? duckdb_set_config(cfg, "access_mode", "READ_ONLY") : DuckDBSuccess;
    if (status == DuckDBSuccess) status = duckdb_open_ext(path, &db, cfg, &message);
    duckdb_destroy_config(&cfg);
    std::string error = message ? message : "opening DuckDB database";
    duckdb_free(message);
    check(status, error.c_str());
    check(duckdb_connect(db, &handle), "connecting to DuckDB");
  }
  void close() {
    if (handle) duckdb_disconnect(&handle);
    if (db) duckdb_close(&db);
  }
  void query(const std::string &sql, result &out) {
    auto status = duckdb_query(handle, sql.c_str(), &out.data);
    check(status, duckdb_result_error(&out.data));
  }
  void exec(const std::string &sql) { result r; query(sql, r); }
};

class statement {
  duckdb_prepared_statement stmt = nullptr;
public:
  statement() = default;
  statement(const statement &) = delete;
  ~statement() { duckdb_destroy_prepare(&stmt); }
  void prepare(connection &db, const std::string &sql) {
    auto status = duckdb_prepare(db.handle, sql.c_str(), &stmt);
    check(status, duckdb_prepare_error(stmt));
  }
  void bind(int i, int64_t v) { check(duckdb_bind_int64(stmt,i,v), "binding integer"); }
  void bind(int i, int v) { bind(i, int64_t(v)); }
  void bind(int i, uint64_t v) { check(duckdb_bind_uint64(stmt,i,v), "binding timepoint"); }
  void bind(int i, double v) { check(duckdb_bind_double(stmt,i,v), "binding double"); }
  void bind(int i, const std::string &v) {
    check(duckdb_bind_varchar_length(stmt,i,v.data(),v.size()), "binding text");
  }
  void bind(int i, std::nullptr_t) { check(duckdb_bind_null(stmt,i), "binding NULL"); }
  template<class... T> void insert(const T &...values) {
    int i = 0;
    (bind(++i, values), ...);
    execute();
  }
  void execute() {
    result r;
    auto status = duckdb_execute_prepared(stmt, &r.data);
    check(status, duckdb_result_error(&r.data));
  }
};

} // namespace luna_duck
