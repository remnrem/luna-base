#pragma once

// Reader adapters for destrat (SQLite and DuckDB).
//
// The shared front end (duckdb/reader.cpp) selects commands, factors, levels
// and variables against Luna's strata model, held in writer_t's metadata maps.
// Each adapter fills that model from its own schema and streams the selected
// datapoints back in legacy extraction order:
//
//   ascending stratum ID, then emission order within the stratum
//
// so that "last emitted wins" when cells are merged, as in destrat.
//
//   SQLite (original schema)   strata/timepoints as stored; ORDER BY strata_id, rowid
//   typed DuckDB (schema v1)   strata synthesised per result table from its
//                              factor-level combinations, numbered by first
//                              emission (min _seq); E/T become timepoints;
//                              rows read in (stratum, _seq) order with only the
//                              requested columns and filters pushed into SQL.
//
// The experimental SQLite-shaped DuckDB schema is recognised and rejected.

#include "db/db.h"

#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace luna_reader {

// A typed value. NULL is an emitted missing value (SQLite); typed DuckDB
// sources do not produce NULL cells (an absent value there is not emitted).
struct cell_t {
  enum kind_t : uint8_t { null, integer, real, boolean, text } kind = null;
  int64_t i = 0;
  double d = 0;
  std::string s;
};

// Text as the SQLite reader renders it: integers and booleans as integers,
// whole-number reals as integers (SQLite NUMERIC affinity), other reals as
// %!.15g, text verbatim, NULL as "".
std::string cell_text(const cell_t &c);

struct datapoint_t {
  int indiv = 0, cmd = 0, var = 0, strata = -1;   // writer_t metadata IDs; strata -1 = root (NULL)
  int tp_key = -1;                                // distinct per timepoint within a source; -1 = none
  timepoint_t tp;
  cell_t value;
};

struct extract_request_t {
  const std::set<int> *strata = nullptr;   // matched strata
  bool baseline = false;                   // root (NULL) stratum only
  bool dump = false;                       // -t: any timepoint; include_root adds NULL strata
  bool include_root = false;
  bool timepoints = false;                 // rows with (true) / without (false) a timepoint
  const std::set<int> *inds = nullptr;     // empty: all
  const std::set<int> *vars = nullptr;     // empty: all
  const std::set<int> *cmds = nullptr;     // empty: all
};

class source_t {
public:
  virtual ~source_t() = default;
  virtual const char *format() const = 0;

  // Fill writer.{individuals,commands,factors,levels,strata,variables} and
  // their name maps. Individuals are unique by name.
  virtual void load(writer_t &w) = 0;

  virtual std::set<std::string> variable_names() = 0;
  virtual int64_t count_values() = 0;
  virtual packets_t enumerate_stratum(int strata_id) = 0;           // (indiv, cmd, var) -> count
  virtual std::map<int, std::set<int>> vars_by_stratum() = 0;      // root reported as 1
  // distinct (command ID, stratum ID or -1 for root) holding selected values
  virtual std::vector<std::pair<int, int>> dump_sources(const std::set<int> &vars, const std::set<int> &inds) = 0;
  virtual void extract(const extract_request_t &req, const std::function<void(const datapoint_t &)> &emit) = 0;

  // Declared factor kinds where the source knows them: 'i' integer, 'r' real,
  // 'n' numeric (integer/real unknown), 't' text.
  virtual std::map<std::string, char> factor_kinds() const = 0;
};

// Detects the format from the file header (not the extension).
std::unique_ptr<source_t> open_source(const std::string &path, bool sqlite_index = true);

std::unique_ptr<source_t> make_sqlite_source(const std::string &path, bool create_index = true);
std::unique_ptr<source_t> make_typed_source(const std::string &path);

} // namespace luna_reader
