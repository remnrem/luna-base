#pragma once

// Typed DuckDB output schema: names, version and type mappings shared by
// the writer and the reader adapter.
//
// Layout (schema version 1)
//
//   luna_schema(key, value)                       format, version, luna_version, created
//   luna_individuals(indiv_ref, indiv_id, source_file, added)
//   luna_executions(exec_ref, indiv_ref, command, cmd_number, parameters, executed)
//   luna_factors(ordinal, name, type)              declaration order; command factors as _CMD
//   luna_tables(table_ref, sql_name, command, factors, core_factors, description)
//   luna_columns(table_ref, ordinal, sql_name, name, role, type, sql_type, description)
//
// One result table per command and runtime factor set (tags included), e.g.
// STATS__CH_E. Rows hold factor levels and variables as columns, plus the
// bookkeeping columns _seq (global emission order), _indiv and _exec.
// Readers must take names, roles and types from luna_columns.

#include <string>

namespace luna_duck {

constexpr const char *schema_format = "luna-typed-results";
constexpr int schema_version = 1;

// luna_columns.role
constexpr const char *role_seq = "seq";
constexpr const char *role_indiv = "indiv";
constexpr const char *role_exec = "exec";
constexpr const char *role_factor = "factor";
constexpr const char *role_interval_start = "interval_start";
constexpr const char *role_interval_stop = "interval_stop";
constexpr const char *role_variable = "variable";

// Logical column types (luna_columns.type) and their DuckDB storage.
enum class ltype { integer, real, text, boolean, timepoint };

inline const char *ltype_name(ltype t)
{
  switch (t) {
  case ltype::integer: return "integer";
  case ltype::real: return "real";
  case ltype::text: return "text";
  case ltype::boolean: return "boolean";
  default: return "timepoint";
  }
}

inline const char *ltype_sql(ltype t)
{
  switch (t) {
  case ltype::integer: return "BIGINT";
  case ltype::real: return "DOUBLE";
  case ltype::text: return "VARCHAR";
  case ltype::boolean: return "BOOLEAN";
  default: return "UBIGINT";
  }
}

inline bool ltype_parse(const std::string &s, ltype *t)
{
  if (s == "integer") *t = ltype::integer;
  else if (s == "real") *t = ltype::real;
  else if (s == "text") *t = ltype::text;
  else if (s == "boolean") *t = ltype::boolean;
  else if (s == "timepoint") *t = ltype::timepoint;
  else return false;
  return true;
}

// Factor types are declared once, in writer_t::set_types() (db/db.cpp), and
// reach the writer through output_store_t::insert_factor(); the stored types
// are recorded in luna_columns. T (intervals) is stored as two timepoint
// columns.
inline std::string interval_start_name(const std::string &f) { return f + "_START"; }
inline std::string interval_stop_name(const std::string &f) { return f + "_STOP"; }

inline std::string quote_ident(const std::string &s)
{
  std::string out(1, '"');
  for (char c : s) { out += c; if (c == '"') out += c; }
  return out + '"';
}

inline std::string quote_literal(const std::string &s)
{
  std::string out(1, '\'');
  for (char c : s) { out += c; if (c == '\'') out += c; }
  return out + '\'';
}

} // namespace luna_duck
