// Typed DuckDB adapter (Luna typed results schema, see duckdb/schema.h).
//
// The catalog gives tables, columns, factors and types; individuals and
// command executions come from luna_individuals / luna_executions. Strata are
// synthesised from each table's distinct factor-level combinations (E and T
// excluded: they are timepoints, as in the SQLite model) and numbered in order
// of first emission (min _seq), which is the order SQLite assigned strata IDs.
// Extraction reads one table at a time, projecting only the requested
// variable columns, joining the selected level combinations (a temp table),
// filtering individuals/executions in SQL, and streaming rows in
// (stratum, _seq) order.

#include "reader_source.h"
#include "catalog.h"
#include "connection.h"
#include "reader_db.h"
#include "schema.h"

#include "helper/helper.h"

#include <algorithm>
#include <charconv>
#include <stdexcept>
#include <unordered_map>

namespace luna_reader {
namespace {
using namespace luna_duck;

// Streaming, typed result reader.
class stream_t {
  duckdb_prepared_statement stmt = nullptr;
  duckdb_result res{};
  bool have_result = false;
  duckdb_data_chunk chunk = nullptr;
  idx_t row = 0, size = 0;
  std::vector<duckdb_type> types;
  std::vector<void *> data;
  std::vector<uint64_t *> valid;

public:
  stream_t(connection &db, const std::string &sql)
  {
    if (duckdb_prepare(db.handle, sql.c_str(), &stmt) != DuckDBSuccess) {
      const std::string e = duckdb_prepare_error(stmt);
      throw std::runtime_error(e);
    }
    duckdb_pending_result pending;
    if (duckdb_pending_prepared_streaming(stmt, &pending) != DuckDBSuccess) {
      const std::string e = duckdb_pending_error(pending);
      duckdb_destroy_pending(&pending);
      throw std::runtime_error(e);
    }
    const auto state = duckdb_execute_pending(pending, &res);
    duckdb_destroy_pending(&pending);
    have_result = true;
    if (state != DuckDBSuccess) throw std::runtime_error(duckdb_result_error(&res));
    const idx_t n = duckdb_column_count(&res);
    for (idx_t i = 0; i < n; ++i) types.push_back(duckdb_column_type(&res, i));
    data.resize(n);
    valid.resize(n);
  }
  stream_t(const stream_t &) = delete;
  ~stream_t()
  {
    if (chunk) duckdb_destroy_data_chunk(&chunk);
    if (have_result) duckdb_destroy_result(&res);
    duckdb_destroy_prepare(&stmt);
  }

  bool next()
  {
    if (chunk && ++row < size) return true;
    while (true) {
      if (chunk) duckdb_destroy_data_chunk(&chunk);
      chunk = duckdb_fetch_chunk(res);
      if (!chunk) {
        const char *e = duckdb_result_error(&res);
        if (e) throw std::runtime_error(e);
        return false;
      }
      size = duckdb_data_chunk_get_size(chunk);
      row = 0;
      if (size == 0) continue;
      for (size_t i = 0; i < types.size(); ++i) {
        duckdb_vector v = duckdb_data_chunk_get_vector(chunk, i);
        data[i] = duckdb_vector_get_data(v);
        valid[i] = duckdb_vector_get_validity(v);
      }
      return true;
    }
  }

  bool null(size_t c) const { return valid[c] && !duckdb_validity_row_is_valid(valid[c], row); }

  int64_t integer(size_t c) const
  {
    switch (types[c]) {
    case DUCKDB_TYPE_BIGINT: return static_cast<int64_t *>(data[c])[row];
    case DUCKDB_TYPE_INTEGER: return static_cast<int32_t *>(data[c])[row];
    case DUCKDB_TYPE_UBIGINT: return int64_t(static_cast<uint64_t *>(data[c])[row]);
    case DUCKDB_TYPE_BOOLEAN: return static_cast<bool *>(data[c])[row];
    default: throw std::runtime_error("unexpected column type for an integer");
    }
  }
  uint64_t unsigned_integer(size_t c) const
  {
    if (types[c] == DUCKDB_TYPE_UBIGINT) return static_cast<uint64_t *>(data[c])[row];
    return uint64_t(integer(c));
  }
  double real(size_t c) const
  {
    if (types[c] == DUCKDB_TYPE_DOUBLE) return static_cast<double *>(data[c])[row];
    return double(integer(c));
  }
  bool boolean(size_t c) const { return integer(c) != 0; }
  std::string text(size_t c) const
  {
    if (types[c] != DUCKDB_TYPE_VARCHAR) throw std::runtime_error("unexpected column type for text");
    auto &v = static_cast<duckdb_string_t *>(data[c])[row];
    return std::string(duckdb_string_t_data(&v), duckdb_string_t_length(v));
  }

  cell_t cell(size_t c, ltype t) const
  {
    cell_t x;
    if (null(c)) return x;
    switch (t) {
    case ltype::integer: x.kind = cell_t::integer; x.i = integer(c); break;
    case ltype::real: x.kind = cell_t::real; x.d = real(c); break;
    case ltype::boolean: x.kind = cell_t::boolean; x.i = boolean(c); break;
    case ltype::timepoint: x.kind = cell_t::integer; x.i = int64_t(unsigned_integer(c)); break;
    default: x.kind = cell_t::text; x.s = text(c);
    }
    return x;
  }
};

std::string real_repr(double d)
{
  char buf[64];
  const auto r = std::to_chars(buf, buf + sizeof buf, d);
  return std::string(buf, r.ptr);
}

// Exact identity of a level combination (used as a map key).
std::string combo_key(const std::vector<cell_t> &v)
{
  std::string k;
  for (const auto &c : v) {
    k += char('0' + c.kind);
    switch (c.kind) {
    case cell_t::integer: k += std::to_string(c.i); break;
    case cell_t::real: k += real_repr(c.d); break;
    case cell_t::text: k += std::to_string(c.s.size()) + ':' + c.s; break;
    default: break;
    }
    k += '\x1f';
  }
  return k;
}

// A factor level as Luna writes it (writer_t::level(int/double) formatting).
std::string level_text(const cell_t &c)
{
  switch (c.kind) {
  case cell_t::integer: return std::to_string(c.i);
  case cell_t::real: return Helper::dbl2str(c.d);
  case cell_t::text: return c.s;
  default: return ".";
  }
}

std::string sql_literal(const cell_t &c)
{
  switch (c.kind) {
  case cell_t::integer: return std::to_string(c.i);
  case cell_t::real: return real_repr(c.d) + "::DOUBLE";
  case cell_t::text: return quote_literal(c.s);
  default: return "NULL";
  }
}

std::string id_list(const std::vector<int64_t> &ids)
{
  std::string r;
  for (auto id : ids) { if (!r.empty()) r += ','; r += std::to_string(id); }
  return r;
}

struct combo_t {
  std::vector<cell_t> values;   // non-time factor levels, in factor-column order
  int sid = 0;
};

struct ttable_t {
  const table_def *def = nullptr;
  std::vector<int> fcols;                       // non-time factor columns
  std::vector<ltype> ftypes;
  int e_col = -1, t_col = -1;                   // epoch / interval-start columns
  std::vector<std::pair<int, int>> vars;        // (column, writer variable ID)
  std::vector<combo_t> combos;
  std::unordered_map<std::string, int> sid_of;  // combo_key -> stratum
};

class typed_source final : public source_t {
  connection db;
  std::string path;
  catalog_t cat;
  std::vector<ttable_t> tables;
  std::unordered_map<int, std::vector<std::pair<int, int>>> sid_loc;   // stratum -> (table, combo)
  std::unordered_map<int64_t, int> indiv_of_ref, cmd_of_exec;
  std::map<int, std::vector<int64_t>> refs_of_indiv, execs_of_cmd;
  std::map<std::string, char> kinds;
  std::map<int64_t, int> epoch_keys;
  std::map<std::pair<uint64_t, uint64_t>, int> interval_keys;

  static std::string col(const ttable_t &t, int c) { return "t." + quote_ident(t.def->columns[c].sql_name); }

  std::string refs_filter(const std::set<int> &inds) const
  {
    std::vector<int64_t> r;
    for (int i : inds) { auto x = refs_of_indiv.find(i); if (x != refs_of_indiv.end()) r.insert(r.end(), x->second.begin(), x->second.end()); }
    return r.empty() ? "FALSE" : "t._indiv IN (" + id_list(r) + ")";
  }
  std::string execs_filter(const std::set<int> &cmds) const
  {
    std::vector<int64_t> r;
    for (int i : cmds) { auto x = execs_of_cmd.find(i); if (x != execs_of_cmd.end()) r.insert(r.end(), x->second.begin(), x->second.end()); }
    return r.empty() ? "FALSE" : "t._exec IN (" + id_list(r) + ")";
  }

  std::vector<std::pair<int, int>> selected_vars(const ttable_t &t, const std::set<int> &vars) const
  {
    std::vector<std::pair<int, int>> r;
    for (const auto &v : t.vars)
      if (vars.empty() || vars.count(v.second)) r.push_back(v);
    return r;
  }

  static std::string any_value(const ttable_t &t, const std::vector<std::pair<int, int>> &vars)
  {
    std::string r;
    for (const auto &v : vars) r += (r.empty() ? "(" : " OR ") + col(t, v.first) + " IS NOT NULL";
    return r + ")";
  }

  std::string factor_select(const ttable_t &t) const
  {
    std::string r;
    for (int c : t.fcols) r += col(t, c) + ", ";
    return r;
  }

  int sid_of_row(const ttable_t &t, const stream_t &s, size_t first) const
  {
    std::vector<cell_t> v;
    for (size_t i = 0; i < t.fcols.size(); ++i) v.push_back(s.cell(first + i, t.ftypes[i]));
    const auto f = t.sid_of.find(combo_key(v));
    if (f == t.sid_of.end()) throw std::runtime_error("level combination not found in " + t.def->sql_name);
    return f->second;
  }

  // Temp table of (sid, factor levels) for the selected combinations.
  void load_selection(const ttable_t &t, const std::vector<int> &combos)
  {
    std::string ddl = "CREATE OR REPLACE TEMP TABLE luna_sel(sid INTEGER";
    for (size_t i = 0; i < t.fcols.size(); ++i) ddl += ", c" + std::to_string(i) + " " + ltype_sql(t.ftypes[i]);
    db.exec(ddl + ")");
    duckdb_appender app = nullptr;
    if (duckdb_appender_create_ext(db.handle, "temp", "main", "luna_sel", &app) != DuckDBSuccess) {
      duckdb_appender_destroy(&app);
      throw std::runtime_error("creating stratum selection");
    }
    auto ok = [&](duckdb_state s) {
      if (s != DuckDBSuccess) {
        const std::string e = duckdb_appender_error(app) ? duckdb_appender_error(app) : "appending stratum selection";
        duckdb_appender_destroy(&app);
        throw std::runtime_error(e);
      }
    };
    for (int ci : combos) {
      const combo_t &c = t.combos[ci];
      ok(duckdb_append_int32(app, c.sid));
      for (const auto &v : c.values) {
        switch (v.kind) {
        case cell_t::integer: ok(duckdb_append_int64(app, v.i)); break;
        case cell_t::real: ok(duckdb_append_double(app, v.d)); break;
        case cell_t::text: ok(duckdb_append_varchar_length(app, v.s.data(), v.s.size())); break;
        default: ok(duckdb_append_null(app));
        }
      }
      ok(duckdb_appender_end_row(app));
    }
    ok(duckdb_appender_close(app));
    duckdb_appender_destroy(&app);
  }

  std::string selection_join(const ttable_t &t) const
  {
    std::string on;
    for (size_t i = 0; i < t.fcols.size(); ++i)
      on += (on.empty() ? "" : " AND ") + col(t, t.fcols[i]) + " IS NOT DISTINCT FROM s.c" + std::to_string(i);
    return " JOIN temp.luna_sel AS s ON " + (on.empty() ? std::string("TRUE") : on);
  }

  int tp_key(const timepoint_t &tp)
  {
    if (tp.is_epoch()) return epoch_keys.emplace(tp.epoch, int(epoch_keys.size() + interval_keys.size())).first->second;
    if (tp.is_interval())
      return interval_keys.emplace(std::make_pair(tp.start, tp.stop), int(epoch_keys.size() + interval_keys.size())).first->second;
    return -1;
  }

public:
  explicit typed_source(const std::string &p) : path(p)
  {
    db.open(p.c_str(), true);
    const auto kind = inspect_schema(db);
    if (kind == schema_kind::experimental)
      throw std::runtime_error(p + " uses the experimental Luna DuckDB schema, which is no longer supported; regenerate it with luna -d");
    if (kind != schema_kind::typed) throw std::runtime_error(p + " is not a Luna results database");
    cat = load_catalog(db);
  }
  ~typed_source() override { db.close(); }

  const char *format() const override { return "DuckDB"; }

  void load(writer_t &w) override
  {
    // individuals: unique by name, in order of first appearance
    {
      cursor c(db, "SELECT indiv_ref, indiv_id FROM luna_individuals ORDER BY indiv_ref");
      while (c.next()) {
        const int64_t ref = std::stoll(c.text(0));
        const std::string name = c.text(1);
        auto f = w.individuals_idmap.find(name);
        int id;
        if (f != w.individuals_idmap.end()) id = f->second;
        else {
          id = int(w.individuals.size()) + 1;
          indiv_t v;
          v.clear();
          v.indiv_id = id;
          v.indiv_name = name;
          w.individuals_idmap[name] = id;
          w.individuals.emplace(id, v);
        }
        indiv_of_ref[ref] = id;
        refs_of_indiv[id].push_back(ref);
      }
    }
    // commands: unique by name and number, as in the SQLite schema
    {
      std::map<std::string, int> by_key;
      cursor c(db, "SELECT exec_ref, command, cmd_number, parameters, executed FROM luna_executions ORDER BY exec_ref");
      while (c.next()) {
        const int64_t ref = std::stoll(c.text(0));
        const std::string name = c.text(1);
        const int number = c.is_null(2) ? 0 : std::stoi(c.text(2));
        const std::string key = name + "." + std::to_string(number);
        auto f = by_key.find(key);
        int id;
        if (f != by_key.end()) id = f->second;
        else {
          id = int(w.commands.size()) + 1;
          by_key[key] = id;
          command_t v;
          v.cmd_id = id;
          v.cmd_name = name;
          v.cmd_number = number;
          v.cmd_parameters = c.text(3);
          v.timestamp = c.text(4);
          w.commands.emplace(id, v);
        }
        w.commands_idmap[name] = id;
        cmd_of_exec[ref] = id;
        execs_of_cmd[id].push_back(ref);
      }
    }
    // factors, in Luna's declaration order
    auto add_factor = [&](const std::string &name, ltype type) {
      if (w.factors_idmap.count(name)) return;
      factor_t f(name);
      f.factor_id = int(w.factors.size()) + 1;
      f.is_numeric = type != ltype::text;
      w.factors_idmap[name] = f.factor_id;
      w.factors.emplace(f.factor_id, f);
      kinds[name] = type == ltype::integer ? 'i' : type == ltype::real ? 'r' : 't';
    };
    for (const auto &f : cat.factors) add_factor(f.name, f.type);
    for (const auto &t : cat.tables) {
      add_factor("_" + t.command, ltype::text);
      for (const auto &c : t.columns)
        if (c.role == role_factor || c.role == role_interval_start) add_factor(c.name, c.type);
    }

    // variables: one per command and name
    std::map<std::string, int> var_ids;
    tables.resize(cat.tables.size());
    for (size_t ti = 0; ti < cat.tables.size(); ++ti) {
      const table_def &d = cat.tables[ti];
      ttable_t &t = tables[ti];
      t.def = &d;
      for (const auto &c : d.columns) {
        if (c.role == role_factor && c.name == "E") t.e_col = c.ordinal;
        else if (c.role == role_factor) { t.fcols.push_back(c.ordinal); t.ftypes.push_back(c.type); }
        else if (c.role == role_interval_start) t.t_col = c.ordinal;
        else if (c.role == role_variable) {
          const std::string key = d.command + ":" + c.name;
          auto f = var_ids.find(key);
          int id;
          if (f != var_ids.end()) id = f->second;
          else {
            id = int(w.variables.size()) + 1;
            var_ids[key] = id;
            const auto cmd = w.commands_idmap.find(d.command);
            var_t v{id, cmd == w.commands_idmap.end() ? -1 : cmd->second, c.name, c.description.empty() ? "." : c.description};
            w.variables.emplace(id, v);
          }
          t.vars.emplace_back(c.ordinal, id);
        }
      }
    }

    // strata: distinct level combinations per table, numbered by first emission
    struct entry_t { int64_t first; int table; std::vector<cell_t> values; };
    std::vector<entry_t> entries;
    for (size_t ti = 0; ti < tables.size(); ++ti) {
      const ttable_t &t = tables[ti];
      std::string q = "SELECT " + factor_select(t) + "min(t._seq) FROM " + quote_ident(t.def->sql_name) + " AS t";
      if (!t.fcols.empty()) q += " GROUP BY ALL";
      stream_t s(db, q);
      while (s.next()) {
        if (s.null(t.fcols.size())) continue;   // empty table
        entry_t e;
        e.table = int(ti);
        for (size_t i = 0; i < t.fcols.size(); ++i) e.values.push_back(s.cell(i, t.ftypes[i]));
        e.first = s.integer(t.fcols.size());
        entries.push_back(std::move(e));
      }
    }
    std::sort(entries.begin(), entries.end(), [](const entry_t &a, const entry_t &b) { return a.first < b.first; });

    strata_t root;
    root.strata_id = 1;
    w.strata[1] = root;
    std::map<std::pair<int, std::string>, int> level_ids;
    auto level_of = [&](int fid, const std::string &text) {
      auto f = level_ids.find({fid, text});
      if (f != level_ids.end()) return w.levels.at(f->second);
      level_t l;
      l.level_id = int(w.levels.size()) + 1;
      l.factor_id = fid;
      l.level_name = text;
      level_ids[{fid, text}] = l.level_id;
      w.levels.emplace(l.level_id, l);
      return l;
    };
    std::map<std::vector<std::pair<int, int>>, int> sid_of_levels;
    int next_sid = 2;
    for (auto &e : entries) {
      ttable_t &t = tables[e.table];
      strata_t s;
      auto add = [&](const std::string &factor, const std::string &level) {
        const factor_t &f = w.factors.at(w.factors_idmap.at(factor));
        s.insert(level_of(f.factor_id, level), f);
      };
      add("_" + t.def->command, t.def->command);
      for (size_t i = 0; i < t.fcols.size(); ++i) add(t.def->columns[t.fcols[i]].name, level_text(e.values[i]));
      if (t.e_col >= 0) add("E", ".");
      if (t.t_col >= 0) add(t.def->columns[t.t_col].name, ".");
      std::vector<std::pair<int, int>> identity;
      for (const auto &l : s.levels) identity.emplace_back(l.first.factor_id, l.second.level_id);
      auto f = sid_of_levels.find(identity);
      int sid;
      if (f != sid_of_levels.end()) sid = f->second;
      else {
        sid = next_sid++;
        sid_of_levels[identity] = sid;
        s.strata_id = sid;
        w.strata[sid] = s;
      }
      t.sid_of[combo_key(e.values)] = sid;
      sid_loc[sid].emplace_back(e.table, int(t.combos.size()));
      t.combos.push_back({std::move(e.values), sid});
    }
  }

  std::set<std::string> variable_names() override
  {
    std::set<std::string> r;
    for (const auto &t : cat.tables)
      for (const auto &c : t.columns)
        if (c.role == role_variable) r.insert(c.name);
    return r;
  }

  int64_t count_values() override
  {
    std::string q;
    for (const auto &t : tables) {
      if (t.vars.empty()) continue;
      std::string sum;
      for (const auto &v : t.vars) sum += (sum.empty() ? "" : " + ") + std::string("count(") + col(t, v.first) + ")";
      q += (q.empty() ? "" : " UNION ALL ") + std::string("SELECT ") + sum + " AS n FROM " + quote_ident(t.def->sql_name) + " AS t";
    }
    if (q.empty()) return 0;
    cursor c(db, "SELECT coalesce(sum(n), 0) FROM (" + q + ")");
    return c.next() ? std::stoll(c.text(0)) : 0;
  }

  packets_t enumerate_stratum(int sid) override
  {
    packets_t packets;
    const auto loc = sid_loc.find(sid);
    if (loc == sid_loc.end()) return packets;
    std::map<std::tuple<int, int, int>, int64_t> counts;
    for (const auto &tc : loc->second) {
      const ttable_t &t = tables[tc.first];
      if (t.vars.empty()) continue;
      const combo_t &combo = t.combos[tc.second];
      std::string where;
      for (size_t i = 0; i < t.fcols.size(); ++i)
        where += (where.empty() ? "" : " AND ") + col(t, t.fcols[i]) + " IS NOT DISTINCT FROM " + sql_literal(combo.values[i]);
      std::string q = "SELECT t._indiv, t._exec";
      for (const auto &v : t.vars) q += ", count(" + col(t, v.first) + ")";
      q += " FROM " + quote_ident(t.def->sql_name) + " AS t" + (where.empty() ? "" : " WHERE " + where) + " GROUP BY 1, 2";
      stream_t s(db, q);
      while (s.next()) {
        const int indiv = indiv_of_ref.at(s.integer(0)), cmd = cmd_of_exec.at(s.integer(1));
        for (size_t i = 0; i < t.vars.size(); ++i) {
          const int64_t n = s.integer(2 + i);
          if (n) counts[std::make_tuple(indiv, cmd, t.vars[i].second)] += n;
        }
      }
    }
    for (const auto &c : counts) {
      packet_t p{};
      p.indiv_id = std::get<0>(c.first);
      p.cmd_id = std::get<1>(c.first);
      p.var_id = std::get<2>(c.first);
      p.value = value_t(int(c.second));
      packets.push_back(p);
    }
    return packets;
  }

  std::map<int, std::set<int>> vars_by_stratum() override
  {
    std::map<int, std::set<int>> r;
    for (const auto &t : tables) {
      if (t.vars.empty()) continue;
      std::string q = "SELECT " + factor_select(t);
      for (size_t i = 0; i < t.vars.size(); ++i) q += (i ? ", " : "") + std::string("count(") + col(t, t.vars[i].first) + ")";
      q += " FROM " + quote_ident(t.def->sql_name) + " AS t";
      if (!t.fcols.empty()) q += " GROUP BY ALL";
      stream_t s(db, q);
      while (s.next()) {
        const int sid = sid_of_row(t, s, 0);
        for (size_t i = 0; i < t.vars.size(); ++i)
          if (s.integer(t.fcols.size() + i)) r[sid].insert(t.vars[i].second);
      }
    }
    return r;
  }

  std::vector<std::pair<int, int>> dump_sources(const std::set<int> &vars, const std::set<int> &inds) override
  {
    std::set<std::pair<int, int>> r;
    for (const auto &t : tables) {
      const auto sel = selected_vars(t, vars);
      if (sel.empty() || t.combos.empty()) continue;
      std::string q = "SELECT DISTINCT " + factor_select(t) + "t._exec FROM " + quote_ident(t.def->sql_name) + " AS t WHERE " + any_value(t, sel);
      if (!inds.empty()) q += " AND " + refs_filter(inds);
      stream_t s(db, q);
      while (s.next()) r.emplace(cmd_of_exec.at(s.integer(t.fcols.size())), sid_of_row(t, s, 0));
    }
    return std::vector<std::pair<int, int>>(r.begin(), r.end());
  }

  void extract(const extract_request_t &req, const std::function<void(const datapoint_t &)> &emit) override
  {
    if (req.baseline || req.strata->empty()) return;   // no root stratum in this schema

    // selected combinations, grouped by table, tables in order of their first stratum
    std::map<int, std::vector<int>> by_table;
    std::map<int, int> first_sid;
    for (int sid : *req.strata) {
      const auto loc = sid_loc.find(sid);
      if (loc == sid_loc.end()) continue;
      for (const auto &tc : loc->second) {
        by_table[tc.first].push_back(tc.second);
        if (!first_sid.count(tc.first)) first_sid[tc.first] = sid;
      }
    }
    std::vector<int> order;
    for (const auto &x : by_table) order.push_back(x.first);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return first_sid[a] < first_sid[b]; });

    datapoint_t dp;
    for (int ti : order) {
      const ttable_t &t = tables[ti];
      const auto sel = selected_vars(t, *req.vars);
      if (sel.empty()) continue;
      const bool has_tp = t.e_col >= 0 || t.t_col >= 0;
      if (!req.dump && req.timepoints && !has_tp) continue;

      load_selection(t, by_table[ti]);
      std::string q = "SELECT s.sid, t._indiv, t._exec";
      q += t.e_col >= 0 ? ", " + col(t, t.e_col) : std::string(", NULL::BIGINT");
      q += t.t_col >= 0 ? ", " + col(t, t.t_col) + ", " + col(t, t.t_col + 1) : std::string(", NULL::UBIGINT, NULL::UBIGINT");
      for (const auto &v : sel) q += ", " + col(t, v.first);
      q += " FROM " + quote_ident(t.def->sql_name) + " AS t" + selection_join(t) + " WHERE " + any_value(t, sel);
      if (!req.dump && has_tp) {
        const std::string tp = t.e_col >= 0 ? col(t, t.e_col) : col(t, t.t_col);
        q += " AND " + tp + (req.timepoints ? " IS NOT NULL" : " IS NULL");
      }
      if (!req.inds->empty()) q += " AND " + refs_filter(*req.inds);
      if (!req.cmds->empty()) q += " AND " + execs_filter(*req.cmds);
      q += " ORDER BY s.sid, t._seq";

      stream_t s(db, q);
      while (s.next()) {
        dp.strata = int(s.integer(0));
        dp.indiv = indiv_of_ref.at(s.integer(1));
        dp.cmd = cmd_of_exec.at(s.integer(2));
        dp.tp = timepoint_t();
        if (!s.null(3)) dp.tp.epoch = int(s.integer(3));
        if (!s.null(4)) { dp.tp.start = s.unsigned_integer(4); dp.tp.stop = s.unsigned_integer(5); }
        dp.tp_key = tp_key(dp.tp);
        for (size_t i = 0; i < sel.size(); ++i) {
          if (s.null(6 + i)) continue;   // not emitted
          dp.var = sel[i].second;
          dp.value = s.cell(6 + i, t.def->columns[sel[i].first].type);
          emit(dp);
        }
      }
    }
    db.exec("DROP TABLE IF EXISTS temp.luna_sel");
  }

  std::map<std::string, char> factor_kinds() const override
  {
    auto k = kinds;
    k["E"] = 'i';
    for (const auto &t : cat.tables)
      for (const auto &c : t.columns)
        if (c.role == role_interval_start) k[c.name] = 't';   // T levels are start_stop labels
    return k;
  }
};

} // namespace

std::unique_ptr<source_t> make_typed_source(const std::string &path) { return std::make_unique<typed_source>(path); }

} // namespace luna_reader
