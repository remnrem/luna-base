#pragma once

// File-format detection for output databases. Header-only and independent of
// the DuckDB library, so that builds without DUCKDB=1 can still recognise (and
// refuse) DuckDB files passed to -a.

#include <cstring>
#include <fstream>
#include <string>

namespace luna_duck {

enum class file_format { absent, empty, sqlite, duckdb, other };

// Inspect the leading bytes of the file, irrespective of its extension.
inline file_format detect_format(const std::string &path)
{
  std::ifstream in(path, std::ios::binary);
  if (!in) return file_format::absent;
  char header[16] = {0};
  in.read(header, sizeof header);
  const std::streamsize n = in.gcount();
  if (n == 0) return file_format::empty;
  if (n == 16 && std::memcmp(header, "SQLite format 3\0", 16) == 0) return file_format::sqlite;
  // DuckDB main header: 8-byte checksum, then the "DUCK" magic bytes.
  if (n >= 12 && std::memcmp(header + 8, "DUCK", 4) == 0) return file_format::duckdb;
  return file_format::other;
}

inline const char *format_name(file_format f)
{
  switch (f) {
  case file_format::absent: return "absent";
  case file_format::empty: return "empty";
  case file_format::sqlite: return "SQLite";
  case file_format::duckdb: return "DuckDB";
  default: return "unrecognised";
  }
}

} // namespace luna_duck
