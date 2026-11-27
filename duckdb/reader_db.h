#pragma once
#include "connection.h"
#include <memory>

namespace luna_duck {

// The reader consumes strings for formatting and integer IDs for joins.
// Projecting columns to VARCHAR gives every result vector the same C API layout.
class cursor {
  result rows;
  duckdb_data_chunk chunk = nullptr;
  idx_t row = 0;
public:
  cursor(connection &db, std::string sql) {
    if (!sql.empty() && sql.back() == ';') sql.pop_back();
    db.query("SELECT COLUMNS(*)::VARCHAR FROM (" + sql + ") AS result_rows",rows);
  }
  ~cursor() { duckdb_destroy_data_chunk(&chunk); }
  bool next() {
    if (chunk && ++row < duckdb_data_chunk_get_size(chunk)) return true;
    duckdb_destroy_data_chunk(&chunk);
    chunk = duckdb_fetch_chunk(rows.data);
    if (auto err = duckdb_result_error(&rows.data)) throw std::runtime_error(err);
    row = 0;
    return chunk != nullptr;
  }
  bool is_null(int col) const {
    auto valid = duckdb_vector_get_validity(duckdb_data_chunk_get_vector(chunk,col));
    return valid && !duckdb_validity_row_is_valid(valid,row);
  }
  std::string text(int col) const {
    if (is_null(col)) return {};
    auto vec = duckdb_data_chunk_get_vector(chunk,col);
    auto &v = static_cast<duckdb_string_t *>(duckdb_vector_get_data(vec))[row];
    return std::string(duckdb_string_t_data(&v),duckdb_string_t_length(v));
  }
};

class reader_db {
  connection db;
public:
  void open(const std::string &file) { db.open(file.c_str(),true); }
  void close() { db.close(); }
  std::unique_ptr<cursor> prepare(const std::string &sql) { return std::make_unique<cursor>(db,sql); }
};
}
