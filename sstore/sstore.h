
//    --------------------------------------------------------------------
//
//    This file is part of Luna.
//
//    LUNA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Luna is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Luna. If not, see <http://www.gnu.org/licenses/>.
//
//    Please see LICENSE.txt for more details.
//
//    --------------------------------------------------------------------

#ifndef __SSTORE_H__
#define __SSTORE_H__

#include "db/sqlwrap.h"
#include <string>
#include <set>
#include "intervals/intervals.h"

struct sstore_key_t {
  std::string id;
  std::string lvl;
  std::string ch;
  sstore_key_t() { } 
  sstore_key_t( const std::string & id , 
		const std::string & lvl , 
		const std::string & ch ) : id(id), lvl(lvl), ch(ch) { }  

  sstore_key_t( const std::string & id ) : id(id) , lvl("") , ch("") { } 

  bool operator<( const sstore_key_t & rhs ) const { 
    if ( id < rhs.id ) return true;
    if ( id > rhs.id ) return false;
  
    if ( lvl < rhs.lvl ) return true;
    if ( lvl > rhs.lvl ) return false;
    
    return ch < rhs.ch;
  }
};

struct sstore_value_t { 

  sstore_value_t() { is_text = is_double = is_vector = false; } 

  bool is_text;
  bool is_double;
  bool is_vector;  
  std::string str_value;
  double      dbl_value;
  std::vector<double>  vec_value;  
};


struct sstore_data_t { 
  std::map<sstore_key_t,sstore_value_t> data;
};


struct sstore_t { 

  sstore_t( const std::string & );
  
  bool attached() { return sql.is_open(); }
  
  bool init();

  bool dettach();

  bool release();
  
  bool index();
  
  bool drop_index();

  // sets  
  void insert_base( const std::string & id , const std::string & value , const std::string * ch = NULL , const std::string * lvl = NULL );
  void insert_base( const std::string & id , const double      & value , const std::string * ch = NULL , const std::string * lvl = NULL);
  void insert_base( const std::string & id , const std::vector<double> & value , const std::string * ch = NULL , const std::string * lvl = NULL);
  
  // 1-based epoch codes
  void insert_epoch( const int e , const std::string & id , const std::string & value , const std::string * ch = NULL , const std::string * lvl = NULL);
  void insert_epoch( const int e , const std::string & id , const double      & value , const std::string * ch = NULL , const std::string * lvl = NULL);
  void insert_epoch( const int e , const std::string & id , const std::vector<double> & value , const std::string * ch = NULL , const std::string * lvl = NULL);

  // tp-based interval codes (i.e. 1 tp = 1e-12 of a second)
  void insert_interval( const uint64_t a , const uint64_t b , const std::string & id , const std::string & value , const std::string * ch = NULL , const std::string * lvl = NULL);  
  void insert_interval( const uint64_t a , const uint64_t b , const std::string & id , const double      & value , const std::string * ch = NULL , const std::string * lvl = NULL);  
  void insert_interval( const uint64_t a , const uint64_t b , const std::string & id , const std::vector<double> & value , const std::string * ch = NULL , const std::string * lvl = NULL);  

  // gets
  sstore_data_t fetch_base(); 
  
  sstore_data_t fetch_epoch( const int e ); 
  std::map<int,sstore_data_t> fetch_epochs(); 

  sstore_data_t fetch_interval( const interval_t & ); 
  std::map<interval_t,sstore_data_t> fetch_intervals(); 
  
  // summaries
  std::map<sstore_key_t,int> keys();
  std::map<sstore_key_t,int> keys_epoch();
  std::map<sstore_key_t,int> keys_interval();
      
private:

  
  SQL sql;

  std::string filename;
  
  
  //
  // Prepared queries
  //

  // sets
  
  sqlite3_stmt * stmt_insert_base;
  sqlite3_stmt * stmt_insert_epoch;
  sqlite3_stmt * stmt_insert_interval;    

  // gets

  sqlite3_stmt * stmt_fetch_base;

  sqlite3_stmt * stmt_fetch_epoch;
  sqlite3_stmt * stmt_fetch_all_epochs;

  sqlite3_stmt * stmt_fetch_interval;
  sqlite3_stmt * stmt_fetch_all_intervals;


  sqlite3_stmt * stmt_fetch_keys;
  sqlite3_stmt * stmt_fetch_keys_epochs;
  sqlite3_stmt * stmt_fetch_keys_intervals;

  // helper functiosn to 

  
};

#endif
