#include "catalog.h"
#include "reader_db.h"
#include <map>
#include <set>
#include <stdexcept>

namespace luna_duck {

std::string join(const std::vector<std::string> &v, const std::string &delim)
{
  std::string out;
  for (size_t i = 0; i < v.size(); ++i) { if (i) out += delim; out += v[i]; }
  return out;
}

std::vector<std::string> split(const std::string &s, char delim)
{
  std::vector<std::string> out;
  if (s.empty()) return out;
  size_t p = 0;
  while (true) {
    const size_t q = s.find(delim, p);
    out.push_back(s.substr(p, q == std::string::npos ? std::string::npos : q - p));
    if (q == std::string::npos) break;
    p = q + 1;
  }
  return out;
}

static std::set<std::string> main_tables(connection &db)
{
  std::set<std::string> names;
  cursor c(db, "SELECT table_name FROM duckdb_tables() WHERE schema_name = 'main'");
  while (c.next()) names.insert(c.text(0));
  return names;
}

schema_kind inspect_schema(connection &db)
{
  const auto names = main_tables(db);
  if (names.empty()) return schema_kind::empty;
  if (names.count("luna_schema")) {
    std::map<std::string, std::string> kv;
    cursor c(db, "SELECT key, value FROM luna_schema");
    while (c.next()) kv[c.text(0)] = c.text(1);
    if (kv["format"] != schema_format)
      throw std::runtime_error("luna_schema does not describe a Luna typed results database");
    int version = 0;
    try { version = std::stoi(kv["version"]); } catch (...) {}
    if (version != schema_version)
      throw std::runtime_error("Luna DuckDB schema version " + kv["version"] + " is not supported (this Luna writes and reads version "
                               + std::to_string(schema_version) + ")");
    return schema_kind::typed;
  }
  if (names.count("datapoints") && names.count("variables") && names.count("strata"))
    return schema_kind::experimental;
  return schema_kind::unknown;
}

void create_catalog(connection &db, const std::string &luna_version, const std::string &created)
{
  db.exec(
    "CREATE TABLE luna_schema(key VARCHAR NOT NULL, value VARCHAR);"
    "CREATE TABLE luna_individuals(indiv_ref BIGINT NOT NULL, indiv_id VARCHAR NOT NULL, source_file VARCHAR, added TIMESTAMP);"
    "CREATE TABLE luna_executions(exec_ref BIGINT NOT NULL, indiv_ref BIGINT NOT NULL, command VARCHAR NOT NULL,"
    " cmd_number INTEGER, parameters VARCHAR, executed TIMESTAMP);"
    "CREATE TABLE luna_factors(ordinal INTEGER NOT NULL, name VARCHAR NOT NULL, type VARCHAR NOT NULL);"
    "CREATE TABLE luna_tables(table_ref INTEGER NOT NULL, sql_name VARCHAR NOT NULL, command VARCHAR NOT NULL,"
    " factors VARCHAR NOT NULL, core_factors VARCHAR NOT NULL, description VARCHAR);"
    "CREATE TABLE luna_columns(table_ref INTEGER NOT NULL, ordinal INTEGER NOT NULL, sql_name VARCHAR NOT NULL,"
    " name VARCHAR NOT NULL, role VARCHAR NOT NULL, type VARCHAR NOT NULL, sql_type VARCHAR NOT NULL, description VARCHAR);");
  db.exec("INSERT INTO luna_schema VALUES ('format', " + quote_literal(schema_format) + "), ('version', '"
          + std::to_string(schema_version) + "'), ('luna_version', " + quote_literal(luna_version) + "), ('created', "
          + quote_literal(created) + ")");
}

catalog_t load_catalog(connection &db)
{
  catalog_t cat;
  std::map<int, size_t> by_ref;
  {
    cursor c(db, "SELECT table_ref, sql_name, command, factors, core_factors, description FROM luna_tables ORDER BY table_ref");
    while (c.next()) {
      table_def t;
      t.ref = std::stoi(c.text(0));
      t.sql_name = c.text(1);
      t.command = c.text(2);
      t.factors = split(c.text(3), ',');
      t.core_factors = c.text(4);
      t.description = c.text(5);
      by_ref[t.ref] = cat.tables.size();
      cat.max_table = std::max(cat.max_table, t.ref);
      cat.tables.push_back(std::move(t));
    }
  }
  {
    cursor c(db, "SELECT table_ref, ordinal, sql_name, name, role, type, description FROM luna_columns ORDER BY table_ref, ordinal");
    while (c.next()) {
      const auto t = by_ref.find(std::stoi(c.text(0)));
      if (t == by_ref.end()) throw std::runtime_error("luna_columns refers to an undefined table");
      column_def col;
      col.ordinal = std::stoi(c.text(1));
      col.sql_name = c.text(2);
      col.name = c.text(3);
      col.role = c.text(4);
      if (!ltype_parse(c.text(5), &col.type)) throw std::runtime_error("unknown column type '" + c.text(5) + "' in luna_columns");
      col.description = c.text(6);
      cat.tables[t->second].columns.push_back(std::move(col));
    }
  }

  // The catalog must describe the physical tables exactly.
  std::map<std::string, std::vector<std::pair<std::string, std::string>>> physical;
  {
    cursor c(db, "SELECT table_name, column_name, data_type FROM duckdb_columns() WHERE schema_name = 'main' ORDER BY table_name, column_index");
    while (c.next()) physical[c.text(0)].emplace_back(c.text(1), c.text(2));
  }
  for (const auto &t : cat.tables) {
    const auto p = physical.find(t.sql_name);
    if (p == physical.end()) throw std::runtime_error("catalog table " + t.sql_name + " is missing from the database");
    if (p->second.size() != t.columns.size())
      throw std::runtime_error("catalog and database disagree on the columns of " + t.sql_name);
    for (size_t i = 0; i < t.columns.size(); ++i) {
      const auto &col = t.columns[i];
      const std::string expected = (col.role == role_seq || col.role == role_indiv || col.role == role_exec) ? "BIGINT" : ltype_sql(col.type);
      if (col.ordinal != int(i) || p->second[i].first != col.sql_name || p->second[i].second != expected)
        throw std::runtime_error("catalog and database disagree on column " + col.sql_name + " of " + t.sql_name);
    }
  }

  auto max_of = [&](const std::string &sql) -> int64_t {
    cursor c(db, sql);
    if (!c.next() || c.is_null(0)) return 0;
    return std::stoll(c.text(0));
  };
  {
    cursor c(db, "SELECT name, type FROM luna_factors ORDER BY ordinal");
    while (c.next()) {
      factor_def f;
      f.name = c.text(0);
      if (!ltype_parse(c.text(1), &f.type)) throw std::runtime_error("unknown factor type '" + c.text(1) + "' in luna_factors");
      cat.factors.push_back(std::move(f));
    }
  }
  cat.max_indiv = max_of("SELECT max(indiv_ref) FROM luna_individuals");
  cat.max_exec = max_of("SELECT max(exec_ref) FROM luna_executions");
  std::string seq;
  for (const auto &t : cat.tables)
    for (const auto &col : t.columns)
      if (col.role == role_seq)
        seq += (seq.empty() ? "" : " UNION ALL ") + std::string("SELECT max(") + quote_ident(col.sql_name) + ") AS m FROM " + quote_ident(t.sql_name);
  if (!seq.empty()) cat.max_seq = max_of("SELECT max(m) FROM (" + seq + ")");
  return cat;
}

} // namespace luna_duck
