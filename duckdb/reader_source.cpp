#include "reader_source.h"
#include "format.h"
#include "db/sqlite3.h"

#include <cmath>
#include <stdexcept>

namespace luna_reader {

std::string cell_text(const cell_t &c)
{
  switch (c.kind) {
  case cell_t::null: return "";
  case cell_t::integer: return std::to_string(c.i);
  case cell_t::boolean: return c.i ? "1" : "0";
  case cell_t::text: return c.s;
  case cell_t::real: {
    // SQLite's NUMERIC affinity stores whole-number reals as integers.
    if (std::isfinite(c.d) && c.d == std::floor(c.d) && c.d > -9.2233720368547e18 && c.d < 9.2233720368547e18)
      return std::to_string(int64_t(c.d));
    char *formatted = sqlite3_mprintf("%!.15g", c.d);
    if (!formatted) throw std::bad_alloc();
    std::string s(formatted);
    sqlite3_free(formatted);
    return s;
  }
  }
  return "";
}

std::unique_ptr<source_t> open_source(const std::string &path, bool sqlite_index)
{
  switch (luna_duck::detect_format(path)) {
  case luna_duck::file_format::sqlite: return make_sqlite_source(path, sqlite_index);
  case luna_duck::file_format::duckdb: return make_typed_source(path);
  case luna_duck::file_format::absent: throw std::runtime_error("could not find stout file " + path);
  default: throw std::runtime_error(path + " is neither an SQLite nor a DuckDB Luna results database");
  }
}

} // namespace luna_reader
