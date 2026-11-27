// Native result I/O. Python bindings deliberately live in the separate luna-api repository.
#include "lunapi/lunapi.h"
#include "duckdb/format.h"
#include "db/factor_types.h"
#include "duckdb/reader_source.h"
#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <climits>
#include <filesystem>
#include <limits>
#include <stdexcept>
#ifdef HAS_DUCKDB
#include "duckdb/writer.h"
#include "duckdb/parquet.h"
#include "duckdb/catalog.h"
#include "cmddefs.h"
#endif

namespace {
void validate_format(const std::string &format)
{
  if (format != "sqlite" && format != "duckdb")
    throw std::invalid_argument("result format must be sqlite or duckdb");
  if (format == "duckdb" && !lunapi_t::has_duckdb())
    throw std::runtime_error("DuckDB support not compiled in");
}

void validate_path(const std::string &path)
{
  if (path.empty() || path == ":memory:" || path == "-" || std::filesystem::is_directory(path))
    throw std::invalid_argument("result output requires a file path");
}

std::vector<std::string> table_factors(const std::string &strata)
{
  std::vector<std::string> out;
  if (strata == "." || strata == "BL") return out;
  if (!strata.empty() && strata.back()=='_') throw std::invalid_argument("invalid result factor set: "+strata);
  std::istringstream in(strata);
  std::string f;
  while (std::getline(in,f,'_')) {
    if (f.empty() || std::find(out.begin(),out.end(),f) != out.end())
      throw std::invalid_argument("invalid result factor set: " + strata);
    out.push_back(f);
  }
  if (out.empty()) throw std::invalid_argument("empty result factor set; use BL");
  return out;
}

std::string text(const rtable_elem_t &x)
{
  if (const auto p = std::get_if<std::string>(&x)) return *p;
  if (const auto p = std::get_if<int>(&x)) return std::to_string(*p);
  if (const auto p = std::get_if<int64_t>(&x)) return std::to_string(*p);
  if (const auto p = std::get_if<bool>(&x)) return *p ? "1" : "0";
  if (const auto p = std::get_if<double>(&x)) {
    std::ostringstream s; s.precision(std::numeric_limits<double>::max_digits10); s << *p; return s.str();
  }
  throw std::invalid_argument("missing factor or ID in result table");
}

uint64_t timepoint(const rtable_elem_t &x)
{
  if (const auto p = std::get_if<int>(&x)) {
    if (*p >= 0) return *p;
  } else if (const auto p = std::get_if<int64_t>(&x)) {
    if (*p >= 0) return *p;
  } else if (const auto p = std::get_if<double>(&x)) {
    // Legacy rtables represent interval endpoints as doubles. Reject imprecise inputs.
    if (std::isfinite(*p) && *p >= 0 && *p <= 0x1p53 && std::floor(*p) == *p) return uint64_t(*p);
  }
  throw std::invalid_argument("invalid integer epoch/interval endpoint in result table");
}

struct replay_table {
  std::string command;
  const rtable_return_t *table;
  size_t id_col;
  std::map<std::string,size_t> factors;
  std::vector<size_t> variables;
};

std::vector<replay_table> validate_tables(const rtables_return_t &tables)
{
  std::vector<replay_table> out;
  for (const auto &cc : tables) for (const auto &ss : cc.second) {
    const auto &cols = std::get<0>(ss.second);
    const auto &data = std::get<1>(ss.second);
    if (cc.first.empty() || cols.empty() || cols.size() != data.size())
      throw std::invalid_argument("invalid result table shape");
    std::map<std::string,size_t> positions;
    for (size_t c=0; c<cols.size(); ++c) {
      if (cols[c].empty() || !positions.emplace(cols[c],c).second || data[c].size() != data[0].size())
        throw std::invalid_argument("duplicate/empty columns or unequal result column lengths");
    }
    if (!positions.count("ID")) throw std::invalid_argument("result table requires an ID column");
    replay_table t{cc.first,&ss.second,positions.at("ID"),{}, {}};
    for (const auto &f : table_factors(ss.first)) {
      if (!positions.count(f) || f == "ID") throw std::invalid_argument("missing factor column: " + f);
      t.factors[f] = positions.at(f);
    }
    const bool t1 = t.factors.count(globals::time_strat+"1"), t2 = t.factors.count(globals::time_strat+"2");
    if (t.factors.count(globals::time_strat))
      throw std::invalid_argument("interval replay requires T1/T2 endpoint columns, not T");
    if (t1 != t2 || (t1 && t.factors.count(globals::epoch_strat)))
      throw std::invalid_argument("interval results require both T1/T2 and cannot also have E");
    for (size_t c=0; c<cols.size(); ++c)
      if (c != t.id_col && !t.factors.count(cols[c])) t.variables.push_back(c);
    for (size_t r=0; r<data[0].size(); ++r) {
      if (!std::holds_alternative<std::string>(data[t.id_col][r]))
        throw std::invalid_argument("result ID must be text");
      for (const auto &f : t.factors) {
        const auto &x = data[f.second][r];
        if (f.first == globals::epoch_strat) {
          auto e=timepoint(x); if (e > INT_MAX) throw std::invalid_argument("epoch exceeds native integer range");
        } else if (f.first == globals::time_strat+"1" || f.first == globals::time_strat+"2") timepoint(x);
        else text(x);
      }
      if (t1 && timepoint(data[t.factors.at(globals::time_strat+"2")][r]) < timepoint(data[t.factors.at(globals::time_strat+"1")][r]))
        throw std::invalid_argument("interval stop precedes start");
    }
    out.push_back(std::move(t));
  }
  return out;
}

#ifdef HAS_DUCKDB
using luna_duck::ltype;
ltype infer(const std::vector<rtable_elem_t> &values)
{
  bool seen=false; ltype type=ltype::text;
  for (const auto &v : values) {
    if (std::holds_alternative<std::monostate>(v)) continue;
    const ltype next=std::holds_alternative<std::string>(v) ? ltype::text
      : std::holds_alternative<double>(v) ? ltype::real
      : std::holds_alternative<bool>(v) ? ltype::boolean : ltype::integer;
    if (!seen) { type=next; seen=true; }
    else if (type != next) {
      if ((type==ltype::integer && next==ltype::real) || (type==ltype::real && next==ltype::integer)) type=ltype::real;
      else throw std::invalid_argument("incompatible mixed types in result column");
    }
  }
  if (type==ltype::real) for (const auto &v : values)
    if (const auto p=std::get_if<int64_t>(&v))
      if (*p < -9007199254740992LL || *p > 9007199254740992LL)
        throw std::invalid_argument("integer cannot be represented exactly in mixed real column");
  return type; // all-missing user columns default to text
}

duckdb_output_schema replay_schema(const std::vector<replay_table> &tables)
{
  duckdb_output_schema schema;
  std::set<std::string> standard;
  for (const auto &f : luna_standard_factors()) standard.insert(f.first);
  for (const auto &t : tables) {
    const auto &cols=std::get<0>(*t.table); const auto &data=std::get<1>(*t.table);
    std::set<std::string> factors;
    for (const auto &f : t.factors) {
      if (f.first==globals::time_strat+"1" || f.first==globals::time_strat+"2") { factors.insert(globals::time_strat); continue; }
      factors.insert(f.first);
      // Standard factor declarations remain authoritative. Custom factors are text
      // unless their actual values provide an unambiguous numeric type.
      if (standard.count(f.first) || globals::cmddefs().is_tag(f.first)) continue;
      auto type=infer(data[f.second]);
      if (type==ltype::boolean) type=ltype::integer;
      auto old=schema.factors.find(f.first);
      if (old!=schema.factors.end() && old->second!=type)
        throw std::invalid_argument("inconsistent factor types across result tables: "+f.first);
      schema.factors[f.first]=type;
    }
    auto &vars=schema.variables[{t.command,factors}];
    for (auto c : t.variables) vars[cols[c]]=infer(data[c]);
  }
  return schema;
}

// Exact-factor DuckDB extracts can go straight through COPY, retaining schema
// types (including all-null columns and unsigned interval endpoints).
void duckdb_parquet(const std::string &database, const std::string &path,
                    const std::string &command, const std::string &strata,
                    const std::set<std::string> &variables, const std::set<std::string> &ids)
{
  using namespace luna_duck;
  connection db; db.open(database.c_str(),true);
  if (inspect_schema(db)!=schema_kind::typed)
    throw std::invalid_argument("Parquet export requires a Luna typed results database");
  const auto cat=load_catalog(db);
  const table_def *table=nullptr;
  for (const auto &t : cat.tables) if (t.command==command) {
    std::set<std::string> factors;
    for (const auto &f : t.factors)
      if (f==globals::time_strat) { factors.insert(f+"1"); factors.insert(f+"2"); }
      else factors.insert(f);
    if ((factors.empty()?"BL":Helper::stringize(factors,"_"))==strata) { table=&t; break; }
  }
  if (!table) throw std::invalid_argument("no matching result table for Parquet export");
  std::map<std::string,std::string> factors, vars;
  for (const auto &c : table->columns) {
    if (c.role==role_factor) factors[c.name]=c.sql_name;
    else if (c.role==role_interval_start) factors[c.name+"1"]=c.sql_name;
    else if (c.role==role_interval_stop) factors[c.name+"2"]=c.sql_name;
    else if (c.role==role_variable && (variables.empty() || variables.count(c.name))) vars[c.name]=c.sql_name;
  }
  for (const auto &v : variables) if (!vars.count(v)) throw std::invalid_argument("unknown result variable: "+v);
  if (vars.empty()) throw std::invalid_argument("no result variables selected");
  std::set<std::string> output_names{"id"};
  auto claim=[&](std::string name) {
    std::transform(name.begin(),name.end(),name.begin(),[](unsigned char c){ return char(std::tolower(c)); });
    if (!output_names.insert(name).second) throw std::invalid_argument("duplicate Parquet output column: "+name);
  };
  std::string query="SELECT i.indiv_id AS \"ID\"", order="i.indiv_id", present;
  for (const auto &f : factors) {
    claim(f.first); const auto col="t."+quote_ident(f.second);
    query+=", "+col+" AS "+quote_ident(f.first); order+=", "+col;
  }
  for (const auto &v : vars) {
    claim(v.first); const auto col="t."+quote_ident(v.second);
    query+=", "+col+" AS "+quote_ident(v.first);
    present+=(present.empty()?"":" OR ")+col+" IS NOT NULL";
  }
  query+=" FROM "+quote_ident(table->sql_name)+" AS t JOIN luna_individuals AS i ON t._indiv=i.indiv_ref WHERE ("+present+")";
  if (!ids.empty()) {
    std::string selected;
    for (const auto &id : ids) selected+=(selected.empty()?"":",")+quote_sql(id,'\'');
    query+=" AND i.indiv_id IN ("+selected+")";
  }
  query+=" ORDER BY "+order+", t._seq";
  db.exec("COPY ("+query+") TO "+quote_sql(path,'\'')+" (FORMAT PARQUET)");
}

#endif

// Reuse the reader adapters without any CLI globals or process-level entry point.
rtables_t read_tables(const std::string &path, const std::set<std::string> &ids,
                     std::vector<std::string> *observed, const std::string &command="",
                     const std::set<std::string> &variables={}, const std::string &selected_strata="")
{
#ifdef HAS_DUCKDB
  auto source=luna_reader::open_source(path, false);
#else
  const auto format=luna_duck::detect_format(path);
  if (format==luna_duck::file_format::duckdb) throw std::runtime_error("DuckDB support not compiled in");
  if (format!=luna_duck::file_format::sqlite) throw std::invalid_argument("not an SQLite Luna results database: "+path);
  auto source=luna_reader::make_sqlite_source(path, false);
#endif
  writer_t metadata;
  source->load(metadata);
  const auto factor_kinds=source->factor_kinds();
  std::set<int> individuals, commands, vars, strata;
  for (const auto &i : metadata.individuals_idmap) if (ids.empty() || ids.count(i.first)) {
    individuals.insert(i.second); if (observed) observed->push_back(i.first);
  }
  if (!ids.empty() && individuals.empty()) return {};
  for (const auto &c : metadata.commands) if (command.empty() || c.second.cmd_name==command) commands.insert(c.first);
  if (!command.empty() && commands.empty()) throw std::invalid_argument("unknown result command: "+command);
  std::set<std::string> matched_variables;
  for (const auto &v : metadata.variables)
    if ((command.empty() || commands.count(v.second.cmd_id)) && (variables.empty() || variables.count(v.second.var_name))) {
      vars.insert(v.first); matched_variables.insert(v.second.var_name);
    }
  for (const auto &v : variables) if (!matched_variables.count(v))
    throw std::invalid_argument("unknown result variable for selected command: "+v);
  for (const auto &s : metadata.strata) {
    std::set<std::string> names;
    for (const auto &f : s.second.levels) {
      const auto &name=f.first.factor_name;
      if (name.empty() || name[0]=='_') continue;
      if (name==globals::time_strat) { names.insert(name+"1"); names.insert(name+"2"); }
      else names.insert(name);
    }
    const auto label=names.empty()?"BL":Helper::stringize(names,"_");
    if (selected_strata.empty() || selected_strata==label) strata.insert(s.first);
  }
  luna_reader::extract_request_t request;
  request.dump=true; request.include_root=selected_strata.empty() || selected_strata=="BL";
  request.inds=&individuals; request.cmds=&commands; request.vars=&vars; request.strata=&strata;
  retval_t ret;
  std::set<std::tuple<std::string,std::string,std::string>> boolean_columns;
  source->extract(request,[&](const luna_reader::datapoint_t &dp) {
    using luna_reader::cell_t;
    if (dp.value.kind==cell_t::null) return; // not written == missing
    const auto &cmd=metadata.commands.at(dp.cmd).cmd_name;
    const auto &var=metadata.variables.at(dp.var).var_name;
    const strata_t empty;
    const auto &s=dp.strata<0 ? empty : metadata.strata.at(dp.strata);
    const retval_factor_t fac(s,dp.tp);
    retval_strata_t level;
    for (const auto &f : s.levels) {
      const auto &name=f.first.factor_name;
      if (name.empty() || name[0]=='_' || name==globals::epoch_strat || name==globals::time_strat) continue;
      const auto kind=factor_kinds.find(name);
      const auto &value=f.second.level_name;
      int64_t integer=0;
      const char *begin=value.data();
      if (!value.empty() && *begin=='+') ++begin;
      const auto parsed=std::from_chars(begin,value.data()+value.size(),integer);
      if (f.first.is_numeric && parsed.ec==std::errc() && parsed.ptr==value.data()+value.size() &&
          kind!=factor_kinds.end() && (kind->second=='i' || kind->second=='n'))
        level.add(retval_factor_level_t(name,integer));
      else if (f.first.is_numeric) {
        double number=0;
        if (!Helper::str2dbl(value,&number)) throw std::invalid_argument("invalid numeric factor: "+name);
        level.add(retval_factor_level_t(name,number));
      } else level.add(retval_factor_level_t(name,value));
    }
    if (dp.tp.is_epoch()) level.add(retval_factor_level_t(globals::epoch_strat,dp.tp.epoch));
    else if (dp.tp.is_interval()) {
      if (dp.tp.start>uint64_t(INT64_MAX) || dp.tp.stop>uint64_t(INT64_MAX))
        throw std::overflow_error("interval exceeds the signed 64-bit native result-table range");
      level.add(retval_factor_level_t(globals::time_strat+"1",int64_t(dp.tp.start)));
      level.add(retval_factor_level_t(globals::time_strat+"2",int64_t(dp.tp.stop)));
    }
    const auto &id=metadata.individuals.at(dp.indiv).indiv_name;
    if (dp.value.kind==cell_t::text) ret.add(id,cmd,fac,var,level,dp.value.s);
    else if (dp.value.kind==cell_t::real) ret.add(id,cmd,fac,var,level,dp.value.d);
    else {
      ret.add(id,cmd,fac,var,level,dp.value.i);
      if (dp.value.kind==cell_t::boolean)
        boolean_columns.emplace(Helper::sanitize(cmd),fac.factors.empty()?"BL":Helper::sanitize(Helper::stringize(fac.factors,"_")),var);
    }
  });
  rtables_t result(ret);
  for (const auto &b : boolean_columns) {
    auto &t=result.tables.at(std::get<0>(b)).at(std::get<1>(b));
    auto c=std::find(t.cols.begin(),t.cols.end(),std::get<2>(b));
    if (c==t.cols.end()) continue;
    for (auto &v : t.data[c-t.cols.begin()]) {
      if (const auto p=std::get_if<int>(&v)) v=(*p!=0);
      else if (const auto p=std::get_if<int64_t>(&v)) v=(*p!=0);
    }
  }
  return result;
}
} // namespace

bool lunapi_t::has_duckdb()
{
#ifdef HAS_DUCKDB
  return true;
#else
  return false;
#endif
}
std::string lunapi_t::duckdb_version()
{
#ifdef HAS_DUCKDB
  return duckdb_library_version();
#else
  return "";
#endif
}

void lunapi_t::output_attach(const std::string &path) { output_attach_format(path,"sqlite"); }
void lunapi_t::output_attach_format(const std::string &path, const std::string &format)
{
  validate_format(format); validate_path(path);
  output_close();
  Helper::deleteFile(path);
  writer.clear();
  try {
#ifdef HAS_DUCKDB
    if (format=="duckdb") writer.attach(path,false,make_duckdb_output(false));
    else
#endif
      writer.attach(path);
    writer.set_types();
    _file_output_path=path; _file_output_plaintext=false;
  } catch (...) { output_close(); throw; }
}

void lunapi_t::write_db(const std::string &path, const std::string &id,
                        const std::string &file, const rtables_return_t &tables)
{ write_db_format(path,id,file,tables,"sqlite"); }

void lunapi_t::write_db_format(const std::string &path, const std::string &id,
                               const std::string &file, const rtables_return_t &tables,
                               const std::string &format)
{
  validate_format(format); validate_path(path);
  const auto replay=validate_tables(tables); // validate before touching the destination
#ifdef HAS_DUCKDB
  const auto schema=format=="duckdb" ? replay_schema(replay) : duckdb_output_schema{};
#endif
  output_close();
  writer.clear();
  Helper::deleteFile(path);
  try {
#ifdef HAS_DUCKDB
    if (format=="duckdb") writer.attach(path,false,make_duckdb_output(false,schema));
    else
#endif
      writer.attach(path);
    writer.set_types();
#ifdef HAS_DUCKDB
    if (format=="duckdb") for (const auto &f : schema.factors)
      if (f.second==ltype::integer || f.second==ltype::real) writer.numeric_factor(f.first,f.second==ltype::integer);
#endif
    writer.begin();
    int number=0;
    for (const auto &t : replay) {
      const auto &cols=std::get<0>(*t.table); const auto &data=std::get<1>(*t.table);
      writer.cmd(t.command,++number,""); writer.level(t.command,"_"+t.command);
      for (size_t r=0; r<data[0].size(); ++r) {
        const auto row_id=text(data[t.id_col][r]);
        writer.id(row_id.empty()?id:row_id,file);
        for (const auto &f : t.factors) {
          const auto &name=f.first; const auto &v=data[f.second][r];
          if (name==globals::epoch_strat) writer.epoch(int(timepoint(v)));
          else if (name==globals::time_strat+"1") writer.interval(interval_t(timepoint(v),timepoint(data[t.factors.at(globals::time_strat+"2")][r])));
          else if (name!=globals::time_strat+"2") writer.level(text(v),name);
        }
        for (const auto c : t.variables) {
          const auto &v=data[c][r];
          if (const auto p=std::get_if<double>(&v)) writer.value(cols[c],*p);
          else if (const auto p=std::get_if<int>(&v)) writer.value(cols[c],*p);
          // Numeric SQLite affinity and the DuckDB declared integer converter
          // both parse this exact decimal representation without narrowing.
          else if (const auto p=std::get_if<int64_t>(&v)) writer.value(cols[c],std::to_string(*p));
          else if (const auto p=std::get_if<bool>(&v)) writer.value(cols[c],int(*p));
          else if (const auto p=std::get_if<std::string>(&v)) writer.value(cols[c],*p);
        }
        writer.timeless();
        for (const auto &f : t.factors) writer.unlevel(f.first);
      }
      writer.unlevel("_"+t.command);
    }
    writer.commit();
    output_close(); writer.use_retval(nullptr); writer.set_types();
  } catch (...) {
    output_close(); writer.use_retval(nullptr); writer.set_types(); throw;
  }
}

std::vector<std::string> lunapi_t::import_db(const std::string &path)
{ return import_db(path,{}); }
std::vector<std::string> lunapi_t::import_db(const std::string &path, const std::set<std::string> &ids)
{
  std::vector<std::string> observed;
  if (!Helper::fileExists(path)) return observed; // preserve existing missing-path behavior
  auto result=read_tables(path,ids,&observed);
  rtables=std::move(result);
  logger << "  read data on " << observed.size() << " individuals from " << path << "\n";
  return observed;
}

void lunapi_t::export_parquet(const std::string &database, const std::string &path,
                              const std::string &command, const std::string &strata,
                              const std::set<std::string> &variables, const std::set<std::string> &ids)
{
#ifndef HAS_DUCKDB
  throw std::runtime_error("Parquet export requires DuckDB support compiled in");
#else
  validate_path(path);
  if (command.empty()) throw std::invalid_argument("Parquet export requires one command");
  const auto factors=table_factors(strata);
  const std::set<std::string> factor_set(factors.begin(),factors.end());
  const std::string key=factor_set.empty()?"BL":Helper::stringize(factor_set,"_");
  if (std::filesystem::absolute(database).lexically_normal()==std::filesystem::absolute(path).lexically_normal() ||
      (std::filesystem::exists(path) && std::filesystem::equivalent(database,path)))
    throw std::invalid_argument("Parquet output must differ from its input database");
  if (luna_duck::detect_format(database)==luna_duck::file_format::duckdb) {
    duckdb_parquet(database,path,command,key,variables,ids);
    return;
  }
  auto tables=read_tables(database,ids,nullptr,command,variables,key);
  const auto c=tables.tables.find(Helper::sanitize(command));
  if (c==tables.tables.end() || !c->second.count(key)) throw std::invalid_argument("no matching result table for Parquet export");
  const auto &t=c->second.at(key);
  std::vector<std::string> types;
  for (const auto &col : t.data) types.push_back(luna_duck::ltype_sql(infer(col)));
  luna_duck::parquet_table output;
  output.path=path; output.columns(t.cols,types);
  for (int r=0; r<t.nrows; ++r) {
    std::vector<luna_reader::cell_t> row;
    for (const auto &col : t.data) {
      const auto &v=col[r]; luna_reader::cell_t x;
      if (const auto p=std::get_if<std::string>(&v)) { x.kind=luna_reader::cell_t::text; x.s=*p; }
      else if (const auto p=std::get_if<double>(&v)) { x.kind=luna_reader::cell_t::real; x.d=*p; }
      else if (const auto p=std::get_if<int>(&v)) { x.kind=luna_reader::cell_t::integer; x.i=*p; }
      else if (const auto p=std::get_if<int64_t>(&v)) { x.kind=luna_reader::cell_t::integer; x.i=*p; }
      else if (const auto p=std::get_if<bool>(&v)) { x.kind=luna_reader::cell_t::boolean; x.i=*p; }
      row.push_back(std::move(x));
    }
    output.append(row);
  }
  output.finish();
#endif
}
