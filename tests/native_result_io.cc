// Standalone libluna consumer: verifies archive linkage and typed Parquet readback.
#include "lunapi/lunapi.h"
#include <chrono>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#ifdef HAS_DUCKDB
#include "duckdb/connection.h"
#endif

int main()
{
  auto *api=lunapi_t::inaugurate();
  const auto dir=std::filesystem::temp_directory_path()/
    ("luna-native-io-"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
  std::filesystem::create_directory(dir);
  try {
    const int64_t large=9007199254740993LL;
    rtables_return_t data;
    data["USER"]["CH_E"]={
      {"ID","CH","E","BIG","FLAG","LABEL","VALUE","EMPTY"},
      {{std::string("one"),std::string("two")},{std::string("EEG"),std::string("EEG")},
       {1,2},{large,large+1},{true,false},{std::string("21.58.17"),std::string("NA")},
       {std::monostate{},2.5},{std::monostate{},std::monostate{}}}};
    data["FACTOR"]["X"]={{"ID","X","BIG"},{{std::string("one")},{large},{large}}};
    data["INTERVAL"]["T1_T2"]={{"ID","T1","T2","BIG"},{{std::string("one")},{large},{large+50},{large}}};
    for (const std::string format : {"sqlite","duckdb"}) {
      if (format=="duckdb" && !api->has_duckdb()) continue;
      // Exercise filesystem paths containing quotes and spaces.
      const auto db=(dir/(format+" results' test.db")).string();
      const auto parquet=(dir/(format+" results' test.parquet")).string();
      api->write_db_format(db,"","",data,format);
      api->import_db(db);
      if (format=="duckdb") {
        const auto factor=api->table("FACTOR","X");
        const auto index=std::find(factor.cols.begin(),factor.cols.end(),"X")-factor.cols.begin();
        if (std::get<int64_t>(factor.data.at(index).at(0))!=large)
          throw std::runtime_error("integer factor lost 64-bit precision");
      }
      const auto interval=api->table("INTERVAL","T1_T2");
      const auto endpoint=std::find(interval.cols.begin(),interval.cols.end(),"T1")-interval.cols.begin();
      if (std::get<int64_t>(interval.data.at(endpoint).at(0))!=large)
        throw std::runtime_error("interval endpoint lost 64-bit precision");
      const auto t=api->table("USER","CH_E");
      const auto c=std::find(t.cols.begin(),t.cols.end(),"BIG")-t.cols.begin();
      if (t.nrows!=2 || std::get<int64_t>(t.data.at(c).at(0))!=large)
        throw std::runtime_error("64-bit result import changed a value");
#ifdef HAS_DUCKDB
      const std::set<std::string> variables=format=="duckdb"
        ? std::set<std::string>{"BIG","FLAG","LABEL","VALUE","EMPTY"}
        : std::set<std::string>{"BIG","FLAG","LABEL","VALUE"};
      api->export_parquet(db,parquet,"USER","CH_E",variables);
      luna_duck::connection check; check.open(nullptr);
      auto sql_string=[](const std::string &value) {
        std::string quoted="'";
        for (char ch : value) { quoted+=ch; if (ch=='\'') quoted+=ch; }
        return quoted+"'";
      };
      const auto quoted=sql_string(parquet);
      luna_duck::result result;
      check.query("SELECT BIG,FLAG,LABEL,VALUE FROM read_parquet("+quoted+") ORDER BY ID",result);
      if (duckdb_row_count(&result.data)!=2 || duckdb_value_int64(&result.data,0,0)!=large ||
          duckdb_value_int64(&result.data,0,1)!=large+1 || !duckdb_value_is_null(&result.data,3,0) ||
          duckdb_value_double(&result.data,3,1)!=2.5 || duckdb_column_type(&result.data,0)!=DUCKDB_TYPE_BIGINT ||
          duckdb_column_type(&result.data,2)!=DUCKDB_TYPE_VARCHAR || duckdb_column_type(&result.data,3)!=DUCKDB_TYPE_DOUBLE)
        throw std::runtime_error("Parquet types/values/missingness changed");
      if (format=="duckdb" && (duckdb_column_type(&result.data,1)!=DUCKDB_TYPE_BOOLEAN ||
          !duckdb_value_boolean(&result.data,1,0) || duckdb_value_boolean(&result.data,1,1)))
        throw std::runtime_error("Parquet boolean type changed");
      auto label=duckdb_value_varchar(&result.data,2,0);
      const bool exact=label && std::string(label)=="21.58.17";
      duckdb_free(label);
      if (!exact) throw std::runtime_error("Parquet text changed");
      if (format=="duckdb") {
        luna_duck::result empty;
        check.query("SELECT EMPTY FROM read_parquet("+quoted+")",empty);
        if (duckdb_column_type(&empty.data,0)!=DUCKDB_TYPE_VARCHAR || !duckdb_value_is_null(&empty.data,0,0))
          throw std::runtime_error("all-missing declared column lost");
        const auto filtered=(dir/"filtered.parquet").string();
        api->export_parquet(db,filtered,"USER","CH_E",{"BIG"},{"two"});
        luna_duck::result selected;
        check.query("SELECT BIG FROM read_parquet("+sql_string(filtered)+")",selected);
        if (duckdb_row_count(&selected.data)!=1 || duckdb_column_count(&selected.data)!=1 ||
            duckdb_value_int64(&selected.data,0,0)!=large+1)
          throw std::runtime_error("Parquet variable/subject selection changed");
      }
#endif
    }
    {
      const auto previous=api->results();
      const auto mixed=(dir/"mixed.db").string();
      rtables_return_t input;
      input["MIXED"]["E"]={{"ID","E","V"},{{std::string("one"),std::string("one")},{1,2},{large,2.5}}};
      api->write_db(mixed,"","",input);
      bool rejected=false;
      try { api->import_db(mixed); } catch (const std::overflow_error &) { rejected=true; }
      if (!rejected || api->results()!=previous)
        throw std::runtime_error("mixed real promotion lost precision or replaced the previous result cache");
    }
    api->output_close();
    lunapi_t::retire();
    std::filesystem::remove_all(dir);
    std::cout << "Native libluna result I/O passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    api->output_close(); lunapi_t::retire();
    std::filesystem::remove_all(dir);
    return 1;
  }
}
