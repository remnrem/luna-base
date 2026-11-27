#pragma once
#include "db/db.h"
#include "schema.h"

// Typed DuckDB output store (see duckdb/schema.h). With append = true the
// file must already hold a Luna typed results database (or be an empty
// DuckDB database); otherwise a new catalog is created.
std::unique_ptr<output_store_t> make_duckdb_output(bool append);

// Per-output schema additions for API-injected tables, never global cmddefs.
struct duckdb_output_schema {
  using key_t = std::pair<std::string,std::set<std::string>>;
  std::map<key_t,std::map<std::string,luna_duck::ltype>> variables;
  std::map<std::string,luna_duck::ltype> factors;
};
std::unique_ptr<output_store_t> make_duckdb_output(bool append, const duckdb_output_schema &schema);
