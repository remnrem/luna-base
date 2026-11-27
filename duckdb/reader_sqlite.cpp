// SQLite adapter (original Luna schema): the stored strata model as is.

#include "reader_source.h"
#include "db/sqlwrap.h"
#include "helper/helper.h"

#include <stdexcept>

namespace luna_reader {
namespace {

std::string sql_ids(const std::set<int> &ids)
{
  std::string r;
  for (int id : ids) { if (!r.empty()) r += ','; r += std::to_string(id); }
  return r;
}

class sqlite_source final : public source_t {
  SQL db;
  std::string path;
  std::map<std::string,char> kinds;

  template <class F> void each(const std::string &q, F f)
  {
    sqlite3_stmt *s = db.prepare(q);
    if (!s) throw std::runtime_error("could not read from " + path + ": " + db.error_message() + " [" + q + "]");
    try { while (db.step(s)) f(s); }
    catch (...) { db.finalise(s); throw; }
    db.finalise(s);
  }

public:
  explicit sqlite_source(const std::string &p, bool create_index) : path(p)
  {
    if (create_index) {
      db.open(p);
      db.synchronous(false);
      // CLI optimization; library reads must not modify their input database.
      db.query("CREATE INDEX IF NOT EXISTS vIndex ON datapoints(strata_id);");
    } else db.open_readonly(p);
  }
  ~sqlite_source() override { db.close(); }

  const char *format() const override { return "SQLite"; }

  void load(writer_t &w) override
  {
    each("SELECT indiv_id,indiv_name FROM individuals;", [&](sqlite3_stmt *s) {
      indiv_t v;
      v.clear();
      v.indiv_id = db.get_int(s, 0);
      v.indiv_name = db.get_text(s, 1);
      w.individuals_idmap[v.indiv_name] = v.indiv_id;
      w.individuals.emplace(v.indiv_id, std::move(v));
    });
    each("SELECT cmd_id,cmd_name,cmd_number,cmd_timestamp,cmd_parameters FROM commands;", [&](sqlite3_stmt *s) {
      command_t v;
      v.cmd_id = db.get_int(s, 0);
      v.cmd_name = db.get_text(s, 1);
      v.cmd_number = db.get_int(s, 2);
      v.timestamp = db.get_text(s, 3);
      v.cmd_parameters = db.get_text(s, 4);
      w.commands_idmap[v.cmd_name] = v.cmd_id;
      w.commands.emplace(v.cmd_id, std::move(v));
    });
    each("SELECT factor_id,factor_name,is_numeric FROM factors;", [&](sqlite3_stmt *s) {
      factor_t v;
      v.factor_id = db.get_int(s, 0);
      v.factor_name = db.get_text(s, 1);
      v.is_numeric = db.get_int(s, 2) == 1;
      kinds[v.factor_name] = v.is_numeric ? 'n' : 't';
      w.factors_idmap[v.factor_name] = v.factor_id;
      w.factors.emplace(v.factor_id, std::move(v));
    });
    each("SELECT level_id,factor_id,level_name FROM levels;", [&](sqlite3_stmt *s) {
      level_t v;
      v.level_id = db.get_int(s, 0);
      v.factor_id = db.get_int(s, 1);
      v.level_name = db.get_text(s, 2);
      w.levels.emplace(v.level_id, std::move(v));
    });
    each("SELECT strata_id,level_id FROM strata;", [&](sqlite3_stmt *s) {
      const int id = db.get_int(s, 0);
      const int level = db.get_int(s, 1);
      auto &stratum = w.strata[id];
      stratum.strata_id = id;
      if (level) {
        const auto ll = w.levels.find(level);
        if (ll == w.levels.end()) Helper::halt("undefined stratum level in " + path);
        const auto ff = w.factors.find(ll->second.factor_id);
        if (ff == w.factors.end()) Helper::halt("undefined factor in " + path);
        stratum.insert(ll->second, ff->second);
      }
    });
    each("SELECT variable_id,variable_name,command_name,variable_label FROM variables;", [&](sqlite3_stmt *s) {
      var_t v;
      v.var_id = db.get_int(s, 0);
      v.var_name = db.get_text(s, 1);
      const auto cmd = w.commands_idmap.find(db.get_text(s, 2));
      v.cmd_id = cmd == w.commands_idmap.end() ? -1 : cmd->second;
      v.var_label = db.get_text(s, 3);
      w.variables.emplace(v.var_id, std::move(v));
    });
  }

  std::set<std::string> variable_names() override
  {
    std::set<std::string> names;
    each("SELECT DISTINCT variable_name FROM variables;", [&](sqlite3_stmt *s) { names.insert(db.get_text(s, 0)); });
    return names;
  }

  int64_t count_values() override
  {
    int64_t n = 0;
    each("SELECT count(*) FROM datapoints;", [&](sqlite3_stmt *s) { n = sqlite3_column_int64(s, 0); });
    return n;
  }

  packets_t enumerate_stratum(int sid) override
  {
    packets_t packets;
    each("SELECT indiv_id,cmd_id,variable_id,count(*) FROM datapoints WHERE strata_id "
           + (sid <= 1 ? std::string("IS NULL") : "= " + std::to_string(sid)) + " GROUP BY indiv_id,cmd_id,variable_id;",
         [&](sqlite3_stmt *s) {
           packet_t p{};
           p.indiv_id = db.get_int(s, 0);
           p.cmd_id = db.get_int(s, 1);
           p.var_id = db.get_int(s, 2);
           p.value = value_t(db.get_int(s, 3));
           packets.push_back(p);
         });
    return packets;
  }

  std::map<int, std::set<int>> vars_by_stratum() override
  {
    std::map<int, std::set<int>> r;
    each("SELECT DISTINCT coalesce(strata_id,1),variable_id FROM datapoints;",
         [&](sqlite3_stmt *s) { r[db.get_int(s, 0)].insert(db.get_int(s, 1)); });
    return r;
  }

  std::vector<std::pair<int, int>> dump_sources(const std::set<int> &vars, const std::set<int> &inds) override
  {
    std::vector<std::pair<int, int>> r;
    std::string q = "SELECT DISTINCT cmd_id,strata_id FROM datapoints WHERE variable_id IN (" + sql_ids(vars) + ")";
    if (!inds.empty()) q += " AND indiv_id IN (" + sql_ids(inds) + ")";
    each(q + " ORDER BY cmd_id,strata_id;", [&](sqlite3_stmt *s) {
      r.emplace_back(db.get_int(s, 0), db.is_null(s, 1) ? -1 : db.get_int(s, 1));
    });
    return r;
  }

  void extract(const extract_request_t &req, const std::function<void(const datapoint_t &)> &emit) override
  {
    const bool tps = req.dump || req.timepoints;
    std::string q = "SELECT d.indiv_id,d.cmd_id,d.variable_id,d.strata_id,d.timepoint_id,d.value";
    if (tps) q += ",t.epoch,t.start,t.stop";
    q += " FROM datapoints d";
    if (tps) q += " LEFT JOIN timepoints t ON t.timepoint_id=d.timepoint_id";
    if (req.dump) {
      q += " WHERE (";
      if (req.include_root) q += "d.strata_id IS NULL";
      if (!req.strata->empty()) {
        if (req.include_root) q += " OR ";
        q += "d.strata_id IN (" + sql_ids(*req.strata) + ")";
      }
      q += ")";
    } else {
      q += req.baseline ? " WHERE d.strata_id IS NULL" : " WHERE d.strata_id IN (" + sql_ids(*req.strata) + ")";
      q += req.timepoints ? " AND d.timepoint_id IS NOT NULL" : " AND d.timepoint_id IS NULL";
    }
    if (!req.inds->empty()) q += " AND d.indiv_id IN (" + sql_ids(*req.inds) + ")";
    if (!req.vars->empty()) q += " AND d.variable_id IN (" + sql_ids(*req.vars) + ")";
    if (!req.cmds->empty()) q += " AND d.cmd_id IN (" + sql_ids(*req.cmds) + ")";
    q += req.baseline ? " ORDER BY d.rowid;" : " ORDER BY d.strata_id,d.rowid;";

    datapoint_t dp;
    each(q, [&](sqlite3_stmt *s) {
      dp.indiv = db.get_int(s, 0);
      dp.cmd = db.get_int(s, 1);
      dp.var = db.get_int(s, 2);
      dp.strata = db.is_null(s, 3) ? -1 : db.get_int(s, 3);
      dp.tp_key = tps && !db.is_null(s, 4) ? db.get_int(s, 4) : -1;
      dp.tp = timepoint_t();
      if (tps && dp.tp_key >= 0) {
        dp.tp.epoch = db.is_null(s, 6) ? -1 : db.get_int(s, 6);
        if (!db.is_null(s, 7)) { dp.tp.start = db.get_uint64(s, 7); dp.tp.stop = db.get_uint64(s, 8); }
      }
      cell_t &c = dp.value;
      switch (sqlite3_column_type(s, 5)) {
      case SQLITE_INTEGER: c.kind = cell_t::integer; c.i = sqlite3_column_int64(s, 5); break;
      case SQLITE_FLOAT: c.kind = cell_t::real; c.d = sqlite3_column_double(s, 5); break;
      case SQLITE_NULL: c.kind = cell_t::null; break;
      default: c.kind = cell_t::text; c.s = db.get_text(s, 5);
      }
      emit(dp);
    });
  }

  std::map<std::string, char> factor_kinds() const override
  {
    return kinds;
  }
};

} // namespace

std::unique_ptr<source_t> make_sqlite_source(const std::string &path, bool create_index) { return std::make_unique<sqlite_source>(path, create_index); }

} // namespace luna_reader
