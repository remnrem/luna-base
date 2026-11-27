#pragma once

// Catalog of the typed DuckDB schema: creation, recognition and loading.
// Errors are reported as std::runtime_error.

#include "connection.h"
#include "schema.h"
#include <cstdint>
#include <string>
#include <vector>

namespace luna_duck {

struct column_def {
  int ordinal = 0;
  std::string sql_name;     // physical column name
  std::string name;         // logical name (factor or variable; bookkeeping role name)
  std::string role;         // see schema.h role_*
  ltype type = ltype::text;
  std::string description;
};

struct table_def {
  int ref = 0;
  std::string sql_name;
  std::string command;
  std::vector<std::string> factors;   // runtime factor set, sorted, tags included
  std::string core_factors;           // cmddefs table identity (comma-separated)
  std::string description;
  std::vector<column_def> columns;    // in ordinal (= physical) order
};

struct factor_def {
  std::string name;                   // includes command factors (_CMD) and tags
  ltype type = ltype::text;
};

struct catalog_t {
  std::vector<factor_def> factors;    // in Luna's declaration order
  std::vector<table_def> tables;
  int64_t max_indiv = 0, max_exec = 0, max_seq = 0;
  int max_table = 0;
};

enum class schema_kind { empty, typed, experimental, unknown };

// Classify an open database. For a typed database, verifies the format and
// version, throwing if this Luna cannot read or extend it.
schema_kind inspect_schema(connection &db);

void create_catalog(connection &db, const std::string &luna_version, const std::string &created);

// Loads table and column definitions and verifies them against the physical
// tables; also returns the highest identifiers in use.
catalog_t load_catalog(connection &db);

std::string join(const std::vector<std::string> &v, const std::string &delim);
std::vector<std::string> split(const std::string &s, char delim);

} // namespace luna_duck
