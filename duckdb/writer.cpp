// Typed DuckDB writer.
//
// writer_t forwards factors, levels, strata, timepoints, variables and values
// to this store by ID. The store rebuilds the -t organisation from those
// calls: one result table per command and runtime factor set, factor levels
// as row values and variables as typed columns.
//
// Rows are assembled as in -t mode: a row is open while consecutive values
// share the individual, command execution, stratum and timepoint; a change
// (or a commit) finalises it. Repeated writes to a cell within an open row
// overwrite it. Rows are appended in schema order, unset columns as NULL,
// through one buffered DuckDB appender per table.

#include "writer.h"
#include "catalog.h"
#include "connection.h"
#include "schema.h"

#include "cmddefs.h"
#include "defs/defs.h"
#include "helper/helper.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <memory>
#include <set>
#include <unordered_map>

namespace {
using namespace luna_duck;

std::string lower(std::string s)
{
  for (auto &c : s) c = char(std::tolower((unsigned char)c));
  return s;
}

std::string now_stamp()
{
  std::time_t t = std::time(nullptr);
  char buf[32];
  std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S", std::localtime(&t));
  return buf;
}

std::string trim(const std::string &s)
{
  const char *ws = " \t\n\r\v\f";
  const auto a = s.find_first_not_of(ws);
  if (a == std::string::npos) return "";
  return s.substr(a, s.find_last_not_of(ws) - a + 1);
}

// Text tokens that denote a missing value in a numeric or boolean column.
bool missing_token(const std::string &t)
{
  return t.empty() || t == "." || t == "NA" || t == "na" || t == "nan" || t == "NaN" || t == "-nan";
}

bool parse_int(const std::string &s, int64_t *out)
{
  const char *a = s.data(), *b = s.data() + s.size();
  if (a != b && *a == '+') ++a;
  const auto r = std::from_chars(a, b, *out);
  return r.ec == std::errc() && r.ptr == b && a != b;
}

bool parse_real(const std::string &s, double *out)
{
  if (s.empty()) return false;
  char *end = nullptr;
  *out = std::strtod(s.c_str(), &end);
  return end == s.c_str() + s.size();
}

std::string real_text(double d)
{
  char buf[64];
  const auto r = std::to_chars(buf, buf + sizeof buf, d);
  return std::string(buf, r.ptr);
}

struct cell_t {
  enum kind_t : uint8_t { null, integer, unsigned_integer, real, text, boolean } kind = null;
  int64_t i = 0;
  uint64_t u = 0;
  double d = 0;
  bool b = false;
  std::string s;
};

struct column_t {
  column_def def;
  vtype_t vtype = VT_UNDEFINED;   // variables only
};

struct table_t {
  table_def def;
  std::vector<column_t> cols;
  std::unordered_map<std::string, int> var_col;     // variable name -> column
  std::unordered_map<std::string, int> factor_col;  // factor name -> column (T: start column)
  std::set<std::string> used;                       // lower-cased SQL column names
  std::unordered_map<int, int> var_cache;           // writer variable ID -> column
  tfac_t core = tfac_t(std::string());
  bool validated = false;                           // checked against cmddefs this session
  duckdb_appender app = nullptr;
};

// Factor cells shared by every row in one (command, stratum).
struct binding_t {
  table_t *table = nullptr;
  std::vector<std::pair<int, cell_t>> levels;
  int epoch_col = -1, interval_col = -1;
};

struct row_t {
  table_t *table = nullptr;
  int indiv = -1, cmd = -1, strata = -1, tp = -1;
  int64_t exec = 0;
  std::vector<cell_t> cells;
};

class duck_writer final : public output_store_t {
  const bool append_mode;
  connection con;
  std::string filename;
  bool in_tx = false;

  statement ins_individual, ins_execution, ins_table, ins_column, ins_factor;
  std::set<std::string> recorded_factors;                   // already in luna_factors

  // records of writer_t's IDs (dense, 1-based)
  std::vector<std::string> var_names;
  std::vector<command_t> commands;
  std::vector<std::string> factor_names;
  std::unordered_map<std::string, ltype> factor_types;      // from writer_t::set_types()
  std::vector<std::pair<int, std::string>> level_records;   // factor ID, level
  std::vector<strata_t> strata;
  std::vector<timepoint_t> timepoints;

  std::unordered_map<int64_t, int64_t> exec_refs;           // (indiv, cmd) -> exec_ref
  std::vector<std::unique_ptr<table_t>> tables;
  std::unordered_map<std::string, table_t *> table_by_key;  // command + factor set
  std::set<std::string> used_table_names;                   // lower-cased
  std::unordered_map<int64_t, binding_t> bindings;          // (cmd, strata) -> binding

  int64_t next_indiv = 0, next_exec = 0, next_seq = 0;
  int next_table = 0;
  row_t row;

  //
  // errors: discard the open individual, then halt
  //

  [[noreturn]] void fail(const std::string &msg)
  {
    abort_session();
    Helper::halt("DuckDB output " + filename + ": " + msg);
    throw std::runtime_error(msg);
  }

  void ok(duckdb_state s, duckdb_appender app, const std::string &what)
  {
    if (s == DuckDBSuccess) return;
    const char *e = app ? duckdb_appender_error(app) : nullptr;
    fail(what + (e ? std::string(": ") + e : std::string()));
  }

  void sql(const std::string &q)
  {
    try { con.exec(q); } catch (const std::exception &e) { fail(e.what()); }
  }

  template <class... T> void insert(statement &st, const T &...v)
  {
    try { st.insert(v...); } catch (const std::exception &e) { fail(e.what()); }
  }

  void destroy_appenders()
  {
    for (auto &t : tables)
      if (t->app) duckdb_appender_destroy(&t->app);
  }

  // Leaves earlier committed individuals in place.
  void abort_session()
  {
    if (!con.handle) return;
    row = row_t();
    destroy_appenders();   // any flush lands in the open transaction
    if (in_tx) { duckdb_result r; duckdb_query(con.handle, "ROLLBACK", &r); duckdb_destroy_result(&r); }
    in_tx = false;
    con.close();
  }

  //
  // tables and columns
  //

  // Names derive from the command and factor set, and from factor and
  // variable names. DuckDB identifiers ignore case, so two different
  // identities can collide; that is a naming problem to fix in Luna.
  std::string claim_table_name(const std::string &name, const std::string &cmd, const std::vector<std::string> &factors)
  {
    if (!used_table_names.insert(lower(name)).second)
      fail("table name " + name + " for " + cmd + " [" + join(factors, ",")
           + "] clashes with an existing table (names compare case-insensitively); rename the factor or tag");
    return name;
  }

  std::string claim_column_name(table_t &t, const std::string &name, const std::string &role)
  {
    if (!t.used.insert(lower(name)).second)
      fail("column " + name + " (" + role + ") of " + t.def.command + " [" + (t.def.core_factors.empty() ? "baseline" : t.def.core_factors)
           + "] clashes with another column of that table (names compare case-insensitively); rename the factor, tag or variable");
    return name;
  }

  // Declared in writer_t::set_types(); undeclared factors (and tags) are text.
  ltype factor_kind(const std::string &f) const
  {
    const auto k = factor_types.find(f);
    return k == factor_types.end() ? ltype::text : k->second;
  }

  static ltype var_ltype(vtype_t v)
  {
    switch (v) {
    case VT_INTEGER: return ltype::integer;
    case VT_REAL: return ltype::real;
    case VT_BOOLEAN: return ltype::boolean;
    default: return ltype::text;
    }
  }

  static std::string column_sql(const column_t &c)
  {
    const bool bookkeeping = c.def.role == role_seq || c.def.role == role_indiv || c.def.role == role_exec;
    return quote_ident(c.def.sql_name) + " " + (bookkeeping ? "BIGINT" : ltype_sql(c.def.type));
  }

  void record_column(table_t &t, const column_t &c)
  {
    insert(ins_column, int64_t(t.def.ref), int64_t(c.def.ordinal), c.def.sql_name, c.def.name, c.def.role,
           std::string(ltype_name(c.def.type)),
           std::string(c.def.role == role_seq || c.def.role == role_indiv || c.def.role == role_exec ? "BIGINT" : ltype_sql(c.def.type)),
           c.def.description);
  }

  column_t &push_column(table_t &t, const std::string &name, const std::string &role, ltype type, const std::string &desc,
                        vtype_t vtype = VT_UNDEFINED)
  {
    column_t c;
    c.def.ordinal = int(t.cols.size());
    c.def.sql_name = claim_column_name(t, role == role_interval_start ? interval_start_name(name)
                                          : role == role_interval_stop ? interval_stop_name(name) : name, role);
    c.def.name = name;
    c.def.role = role;
    c.def.type = type;
    c.def.description = desc;
    c.vtype = vtype;
    t.cols.push_back(c);
    const int idx = c.def.ordinal;
    if (role == role_variable) t.var_col[name] = idx;
    else if (role == role_factor || role == role_interval_start) t.factor_col[name] = idx;
    return t.cols.back();
  }

  static std::string core_string(const tfac_t &core) { return join(std::vector<std::string>(core.fac.begin(), core.fac.end()), ","); }

  std::vector<std::string> registered_vars(const table_t &t) const
  {
    std::vector<std::string> out;
    for (const auto &v : globals::cmddefs().fetch_vars(t.def.command, t.def.core_factors))
      if (!globals::cmddefs().is_hidden_var(t.def.command, t.core, v)) out.push_back(v);
    return out;
  }

  table_t *create_table(const std::string &key, const std::string &cmd, const std::vector<std::string> &factors, const strata_t &s)
  {
    auto t = std::make_unique<table_t>();
    t->core = tfac_t(s);
    t->def.ref = ++next_table;
    t->def.command = cmd;
    t->def.factors = factors;
    t->def.core_factors = core_string(t->core);
    t->def.sql_name = claim_table_name(factors.empty() ? cmd : cmd + "__" + join(factors, "_"), cmd, factors);
    t->def.description = globals::cmddefs().fetch_desc_tbl(cmd, t->def.core_factors);
    if (t->def.description == ".") t->def.description.clear();

    push_column(*t, "_seq", role_seq, ltype::integer, "emission order");
    push_column(*t, "_indiv", role_indiv, ltype::integer, "luna_individuals.indiv_ref");
    push_column(*t, "_exec", role_exec, ltype::integer, "luna_executions.exec_ref");
    for (const auto &f : factors) {
      const ltype ft = factor_kind(f);
      if (ft == ltype::timepoint) {
        push_column(*t, f, role_interval_start, ltype::timepoint, "interval start (time-points)");
        push_column(*t, f, role_interval_stop, ltype::timepoint, "interval stop (time-points)");
      } else
        push_column(*t, f, role_factor, ft, "");
    }
    for (const auto &v : registered_vars(*t)) {
      const vtype_t vt = globals::cmddefs().fetch_type_var(cmd, t->core, v);
      if (vt == VT_UNDEFINED)
        fail("variable " + v + " of " + cmd + " [" + t->def.core_factors + "] has no type in cmddefs");
      push_column(*t, v, role_variable, var_ltype(vt), globals::cmddefs().fetch_desc_var(cmd, t->def.core_factors, v), vt);
    }

    std::string ddl = "CREATE TABLE " + quote_ident(t->def.sql_name) + "(";
    for (size_t i = 0; i < t->cols.size(); ++i) ddl += (i ? ", " : "") + column_sql(t->cols[i]);
    sql(ddl + ")");
    insert(ins_table, int64_t(t->def.ref), t->def.sql_name, cmd, join(factors, ","), t->def.core_factors, t->def.description);
    for (const auto &c : t->cols) record_column(*t, c);
    t->validated = true;

    table_t *p = t.get();
    table_by_key[key] = p;
    tables.push_back(std::move(t));
    return p;
  }

  // Add a typed variable column: flush and close the appender first.
  int add_variable_column(table_t &t, const std::string &var, vtype_t vt)
  {
    if (t.app) {
      ok(duckdb_appender_close(t.app), t.app, "flushing " + t.def.sql_name);
      duckdb_appender_destroy(&t.app);
    }
    column_t &c = push_column(t, var, role_variable, var_ltype(vt), globals::cmddefs().fetch_desc_var(t.def.command, t.def.core_factors, var), vt);
    sql("ALTER TABLE " + quote_ident(t.def.sql_name) + " ADD COLUMN " + column_sql(c));
    record_column(t, c);
    if (row.table == &t) row.cells.resize(t.cols.size());
    return c.def.ordinal;
  }

  // First use of a stored table in this session (append): reconcile with cmddefs.
  void validate(table_t &t)
  {
    if (t.validated) return;
    t.validated = true;
    for (const auto &c : t.cols)
      if ((c.def.role == role_factor || c.def.role == role_interval_start) && c.def.type != factor_kind(c.def.name))
        fail("factor " + c.def.name + " of " + t.def.sql_name + " is stored as " + ltype_name(c.def.type)
             + " but this Luna declares it " + ltype_name(factor_kind(c.def.name)) + " (writer_t::set_types)");
    for (const auto &v : registered_vars(t)) {
      const vtype_t vt = globals::cmddefs().fetch_type_var(t.def.command, t.core, v);
      if (vt == VT_UNDEFINED) continue;
      const auto c = t.var_col.find(v);
      if (c == t.var_col.end()) add_variable_column(t, v, vt);
      else if (t.cols[c->second].def.type != var_ltype(vt))
        fail("incompatible type for " + t.def.command + " [" + t.def.core_factors + "] " + v + ": stored as "
             + ltype_name(t.cols[c->second].def.type) + ", cmddefs now declares " + vtype_str(vt));
      else t.cols[c->second].vtype = vt;
    }
  }

  void load_existing()
  {
    catalog_t cat;
    try { cat = load_catalog(con); } catch (const std::exception &e) { fail(e.what()); }
    next_indiv = cat.max_indiv;
    next_exec = cat.max_exec;
    next_seq = cat.max_seq;
    next_table = cat.max_table;
    for (const auto &f : cat.factors) recorded_factors.insert(f.name);
    for (auto &d : cat.tables) {
      auto t = std::make_unique<table_t>();
      t->core = tfac_t(d.core_factors);
      for (const auto &cd : d.columns) {
        column_t c;
        c.def = cd;
        t->used.insert(lower(cd.sql_name));
        if (cd.role == role_variable) t->var_col[cd.name] = cd.ordinal;
        else if (cd.role == role_factor || cd.role == role_interval_start) {
          t->factor_col[cd.name] = cd.ordinal;
        }
        t->cols.push_back(c);
      }
      used_table_names.insert(lower(d.sql_name));
      const std::string key = d.command + '\x1f' + join(d.factors, ",");
      t->def = std::move(d);
      table_by_key[key] = t.get();
      tables.push_back(std::move(t));
    }
  }

  //
  // rows
  //

  binding_t &bind(int cmd_id, int strata_id)
  {
    const int64_t key = (int64_t(cmd_id) << 32) | uint32_t(strata_id);
    auto b = bindings.find(key);
    if (b != bindings.end()) return b->second;

    if (cmd_id < 1 || cmd_id >= int(commands.size())) fail("unknown command ID");
    if (strata_id < 1 || strata_id >= int(strata.size())) fail("unknown stratum ID");
    const std::string &cmd = commands[cmd_id].cmd_name;
    const strata_t &s = strata[strata_id];
    std::vector<std::string> factors;
    for (const auto &l : s.levels)
      if (!l.first.factor_name.empty() && l.first.factor_name[0] != '_') factors.push_back(l.first.factor_name);
    std::sort(factors.begin(), factors.end());

    const std::string tkey = cmd + '\x1f' + join(factors, ",");
    auto t = table_by_key.find(tkey);
    table_t *table = t == table_by_key.end() ? create_table(tkey, cmd, factors, s) : t->second;
    validate(*table);

    binding_t nb;
    nb.table = table;
    for (const auto &l : s.levels) {
      const std::string &f = l.first.factor_name;
      if (f.empty() || f[0] == '_') continue;
      const int col = table->factor_col.at(f);
      const ltype ft = table->cols[col].def.type;
      if (f == globals::epoch_strat && ft == ltype::integer) { nb.epoch_col = col; continue; }
      if (ft == ltype::timepoint) { nb.interval_col = col; continue; }
      cell_t c;
      const std::string &lv = l.second.level_name;
      if (ft == ltype::text) { c.kind = cell_t::text; c.s = lv; }
      else if (missing_token(trim(lv))) c.kind = cell_t::null;
      else if (ft == ltype::integer) {
        if (!parse_int(trim(lv), &c.i)) fail(cmd + ": level '" + lv + "' of factor " + f + " is not an integer");
        c.kind = cell_t::integer;
      } else {
        if (!parse_real(trim(lv), &c.d)) fail(cmd + ": level '" + lv + "' of factor " + f + " is not a number");
        c.kind = cell_t::real;
      }
      nb.levels.emplace_back(col, std::move(c));
    }
    return bindings.emplace(key, std::move(nb)).first->second;
  }

  int64_t exec_ref(int indiv, int cmd_id)
  {
    const int64_t key = (int64_t(indiv) << 32) | uint32_t(cmd_id);
    auto e = exec_refs.find(key);
    if (e != exec_refs.end()) return e->second;
    const command_t &c = commands.at(cmd_id);
    const int64_t ref = ++next_exec;
    insert(ins_execution, ref, int64_t(indiv), c.cmd_name, int64_t(c.cmd_number), c.cmd_parameters, now_stamp());
    exec_refs[key] = ref;
    return ref;
  }

  void set(int col, cell_t c) { row.cells[col] = std::move(c); }

  void start_row(binding_t &b, int indiv, int64_t exec, int cmd, int strata_id, int tp)
  {
    row.table = b.table;
    row.indiv = indiv; row.exec = exec; row.cmd = cmd; row.strata = strata_id; row.tp = tp;
    row.cells.assign(b.table->cols.size(), cell_t());
    for (const auto &l : b.levels)
      if (l.second.kind != cell_t::null) set(l.first, l.second);
    if (tp > 0) {
      const timepoint_t &t = timepoints.at(tp);
      if (b.epoch_col >= 0 && t.is_epoch()) { cell_t c; c.kind = cell_t::integer; c.i = t.epoch; set(b.epoch_col, c); }
      if (b.interval_col >= 0 && t.is_interval()) {
        cell_t a; a.kind = cell_t::unsigned_integer; a.u = t.start; set(b.interval_col, a);
        cell_t z; z.kind = cell_t::unsigned_integer; z.u = t.stop; set(b.interval_col + 1, z);
      }
    }
  }

  void finish_row()
  {
    table_t *t = row.table;
    if (!t) return;
    if (!t->app)
      ok(duckdb_appender_create(con.handle, nullptr, t->def.sql_name.c_str(), &t->app), t->app, "creating appender for " + t->def.sql_name);
    duckdb_appender a = t->app;
    ok(duckdb_append_int64(a, ++next_seq), a, "appending");
    ok(duckdb_append_int64(a, row.indiv), a, "appending");
    ok(duckdb_append_int64(a, row.exec), a, "appending");
    for (size_t i = 3; i < row.cells.size(); ++i) {
      const cell_t &c = row.cells[i];
      duckdb_state s;
      switch (c.kind) {
      case cell_t::integer: s = duckdb_append_int64(a, c.i); break;
      case cell_t::unsigned_integer: s = duckdb_append_uint64(a, c.u); break;
      case cell_t::real: s = duckdb_append_double(a, c.d); break;
      case cell_t::boolean: s = duckdb_append_bool(a, c.b); break;
      case cell_t::text: s = duckdb_append_varchar_length(a, c.s.data(), c.s.size()); break;
      default: s = duckdb_append_null(a);
      }
      ok(s, a, "appending to " + t->def.sql_name);
    }
    ok(duckdb_appender_end_row(a), a, "appending to " + t->def.sql_name);
    row.table = nullptr;
  }

  void flush_all()
  {
    finish_row();
    for (auto &t : tables)
      if (t->app) ok(duckdb_appender_flush(t->app), t->app, "flushing " + t->def.sql_name);
  }

  // Convert a value to its column's declared type; anything that cannot be
  // represented without loss is an error to be fixed in Luna.
  cell_t convert(const value_t &x, const column_t &col, const table_t &t)
  {
    cell_t c;
    if (x.missing) return c;
    auto bad = [&](const std::string &what) -> cell_t {
      fail(t.def.command + " [" + (t.def.core_factors.empty() ? "baseline" : t.def.core_factors) + "] " + col.def.name
           + " is declared " + ltype_name(col.def.type) + " but received " + what);
    };
    switch (col.def.type) {
    case ltype::real:
      if (x.numeric) { if (!std::isnan(x.d)) { c.kind = cell_t::real; c.d = x.d; } return c; }
      if (x.integer) { c.kind = cell_t::real; c.d = x.i; return c; }
      {
        const std::string s = trim(x.s);
        if (missing_token(s)) return c;
        if (!parse_real(s, &c.d)) return bad("text '" + x.s + "'");
        if (!std::isnan(c.d)) c.kind = cell_t::real;
        return c;
      }
    case ltype::integer:
      if (x.integer) { c.kind = cell_t::integer; c.i = x.i; return c; }
      if (x.numeric) {
        if (std::isnan(x.d)) return c;
        if (x.d != std::trunc(x.d) || !(x.d > -0x1p63 && x.d < 0x1p63)) return bad("non-integer value " + real_text(x.d));
        c.kind = cell_t::integer; c.i = int64_t(x.d); return c;
      }
      {
        const std::string s = trim(x.s);
        if (missing_token(s)) return c;
        if (parse_int(s, &c.i)) { c.kind = cell_t::integer; return c; }
        double d;
        if (parse_real(s, &d) && d == std::trunc(d) && d > -0x1p63 && d < 0x1p63) { c.kind = cell_t::integer; c.i = int64_t(d); return c; }
        return bad("text '" + x.s + "'");
      }
    case ltype::boolean:
      {
        if (x.integer || x.numeric) {
          const double d = x.integer ? double(x.i) : x.d;
          if (x.numeric && std::isnan(d)) return c;
          if (d != 0 && d != 1) return bad("value " + (x.integer ? std::to_string(x.i) : real_text(x.d)));
          c.kind = cell_t::boolean; c.b = d == 1; return c;
        }
        const std::string s = lower(trim(x.s));
        if (missing_token(s)) return c;
        if (s == "1" || s == "true") { c.kind = cell_t::boolean; c.b = true; return c; }
        if (s == "0" || s == "false") { c.kind = cell_t::boolean; c.b = false; return c; }
        return bad("text '" + x.s + "'");
      }
    default:   // text: numbers written to text columns keep their exact value
      c.kind = cell_t::text;
      if (x.integer) c.s = std::to_string(x.i);
      else if (x.numeric) { if (std::isnan(x.d)) { c.kind = cell_t::null; return c; } c.s = real_text(x.d); }
      else c.s = x.s;
      return c;
    }
  }

public:
  explicit duck_writer(bool append) : append_mode(append) {}

  ~duck_writer() override { abort_session(); }

  bool attach(const std::string &path, bool readonly, writer_t *caller) override
  {
    if (readonly) throw std::logic_error("DuckDB output requires write access");
    filename = Helper::expand(path);
    try { con.open(filename.c_str()); } catch (const std::exception &e) { fail(e.what()); }

    schema_kind kind = schema_kind::empty;
    try { kind = inspect_schema(con); } catch (const std::exception &e) { fail(e.what()); }
    if (kind == schema_kind::experimental)
      fail("this file uses the experimental Luna DuckDB schema, which is no longer supported; regenerate it with -d");
    if (kind == schema_kind::unknown) fail("not a Luna results database");
    if (kind == schema_kind::typed && !append_mode) fail("database already exists (internal error: expected a new file)");

    for (const char *n : {"luna_schema", "luna_individuals", "luna_executions", "luna_tables", "luna_columns"})
      used_table_names.insert(n);

    if (kind == schema_kind::empty) {
      sql("BEGIN TRANSACTION");
      try { create_catalog(con, trim(globals::version), now_stamp()); } catch (const std::exception &e) { fail(e.what()); }
      sql("COMMIT");
    } else
      load_existing();

    try {
      ins_individual.prepare(con, "INSERT INTO luna_individuals VALUES (?, ?, ?, ?)");
      ins_execution.prepare(con, "INSERT INTO luna_executions VALUES (?, ?, ?, ?, ?, ?)");
      ins_table.prepare(con, "INSERT INTO luna_tables VALUES (?, ?, ?, ?, ?, ?)");
      ins_column.prepare(con, "INSERT INTO luna_columns VALUES (?, ?, ?, ?, ?, ?, ?, ?)");
      ins_factor.prepare(con, "INSERT INTO luna_factors VALUES (?, ?, ?)");
    } catch (const std::exception &e) { fail(e.what()); }

    // slot 0 unused: writer_t IDs start at 1
    var_names.assign(1, "");
    commands.assign(1, command_t());
    factor_names.assign(1, "");
    level_records.assign(1, {0, ""});
    strata.assign(1, strata_t());
    timepoints.assign(1, timepoint_t());

    caller->set_types();
    return true;
  }

  bool attached() override { return con.handle != nullptr; }
  std::string name() const override { return filename; }

  bool dettach() override
  {
    if (!con.handle) return true;
    if (in_tx) { abort_session(); return true; }   // an individual was not committed: discard it
    flush_all();
    for (auto &t : tables)
      if (t->app) {
        ok(duckdb_appender_close(t->app), t->app, "closing " + t->def.sql_name);
        duckdb_appender_destroy(&t->app);
      }
    con.close();
    return true;
  }

  void begin() override
  {
    if (in_tx) return;
    flush_all();
    sql("BEGIN TRANSACTION");
    in_tx = true;
  }

  void commit() override
  {
    flush_all();
    if (!in_tx) return;
    sql("COMMIT");
    in_tx = false;
  }

  indiv_t insert_individual(const std::string &name, const std::string &file) override
  {
    const int64_t ref = ++next_indiv;
    insert(ins_individual, ref, name, file, now_stamp());
    indiv_t i;
    i.indiv_id = int(ref);
    i.indiv_name = name;
    i.file_name = file;
    return i;
  }

  var_t insert_variable(const std::string &name, const std::string &cmd, const std::string &label) override
  {
    var_names.push_back(name);
    return {int(var_names.size() - 1), -1, name, label};
  }

  command_t insert_command(const std::string &name, int number, const std::string &stamp, const std::string &params) override
  {
    command_t c;
    c.cmd_id = int(commands.size());
    c.cmd_name = name;
    c.cmd_number = number;
    c.cmd_parameters = params;
    c.timestamp = stamp;
    commands.push_back(c);
    return c;
  }

  factor_t insert_factor(const std::string &name, bool numeric, bool integer) override
  {
    const ltype type = name == globals::time_strat ? ltype::timepoint
                       : !numeric ? ltype::text : integer ? ltype::integer : ltype::real;
    factor_types[name] = type;
    // Declaration order defines the factor order readers use for labels.
    if (recorded_factors.insert(name).second)
      insert(ins_factor, int64_t(recorded_factors.size()), name, std::string(ltype_name(type)));
    factor_t f(name);
    f.is_numeric = numeric;
    f.factor_id = int(factor_names.size());
    factor_names.push_back(name);
    return f;
  }

  level_t insert_level(const std::string &name, int fid) override
  {
    level_t l;
    l.level_id = int(level_records.size());
    l.factor_id = fid;
    l.level_name = name;
    level_records.emplace_back(fid, name);
    return l;
  }

  strata_t insert_strata(const strata_t &s) override
  {
    strata_t e = s;
    e.strata_id = int(strata.size());
    strata.push_back(e);
    return e;
  }

  timepoint_t insert_epoch_timepoint(int e) override
  {
    timepoint_t t;
    t.timepoint_id = int(timepoints.size());
    t.epoch = e;
    timepoints.push_back(t);
    return t;
  }

  timepoint_t insert_interval_timepoint(const interval_t &span) override
  {
    timepoint_t t;
    t.timepoint_id = int(timepoints.size());
    t.start = span.start;
    t.stop = span.stop;
    timepoints.push_back(t);
    return t;
  }

  bool insert_value(int indiv, int cmd, int var, int sid, int tid, const value_t &x) override
  {
    const int strata_id = sid < 0 ? 1 : sid;   // -1: baseline (strata 1 is the empty stratum)
    const int tp = tid < 0 ? -1 : tid;
    const int64_t exec = exec_ref(indiv, cmd);

    if (row.table && (row.indiv != indiv || row.exec != exec || row.strata != strata_id || row.tp != tp)) finish_row();
    binding_t &b = bind(cmd, strata_id);
    if (!row.table) start_row(b, indiv, exec, cmd, strata_id, tp);
    table_t &t = *row.table;

    auto vc = t.var_cache.find(var);
    int col;
    if (vc != t.var_cache.end()) col = vc->second;
    else {
      const std::string &vname = var_names.at(var);
      auto c = t.var_col.find(vname);
      if (c != t.var_col.end()) col = c->second;
      else {
        const vtype_t vt = globals::cmddefs().fetch_type_var(t.def.command, t.core, vname);
        if (vt == VT_UNDEFINED)
          fail("variable " + vname + " of " + t.def.command + " [" + (t.def.core_factors.empty() ? "baseline" : t.def.core_factors)
               + "] is not defined in cmddefs");
        col = add_variable_column(t, vname, vt);
      }
      t.var_cache[var] = col;
    }
    set(col, convert(x, t.cols[col], t));
    return true;
  }
};
} // namespace

std::unique_ptr<output_store_t> make_duckdb_output(bool append) { return std::make_unique<duck_writer>(append); }
