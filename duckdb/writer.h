#pragma once
#include "db/db.h"

// Typed DuckDB output store (see duckdb/schema.h). With append = true the
// file must already hold a Luna typed results database (or be an empty
// DuckDB database); otherwise a new catalog is created.
std::unique_ptr<output_store_t> make_duckdb_output(bool append);
