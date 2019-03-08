
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

#ifndef __LOCDB_H__
#define __LOCDB_H__

#include "sqlwrap.h"
#include "retval.h"
#include <string>
#include "intervals/intervals.h"


class writer_t;
extern writer_t writer;

struct value_t;
struct strata_t;
struct indiv_t;
struct var_t;
struct timepoint_t;
struct command_t;
struct factor_t;
struct level_t;



  //
  // Helper classes
  //
  

  // struct for: factors, levels, strata, variables, individuals, commands, timepoints, values

  struct factor_t {
      
    factor_t() { }   

    factor_t( const std::string & factor_name ) : factor_name( factor_name ) 
    {
      is_numeric = false;
      factor_id  = -1;
    }

    int factor_id;
    std::string factor_name;
    bool is_numeric;
    bool operator< ( const factor_t & rhs ) const { return factor_id < rhs.factor_id; }
  };
  
  struct level_t { 
    level_t() { level_id = -1; factor_id = -1; level_name = "."; } 
    int level_id;
    int factor_id;
    std::string level_name;  
    bool operator< ( const level_t & rhs ) const 
    {
      if ( factor_id == rhs.factor_id ) return level_id < rhs.level_id;
      return factor_id < rhs.factor_id;
    }  
    
  };
  
  struct indiv_t {
    void clear() { indiv_id = -1; indiv_name = ""; file_name = ""; }
    int indiv_id;
    std::string indiv_name;
    std::string file_name;
    bool operator< ( const indiv_t & rhs ) const { return indiv_id < rhs.indiv_id; } 
  };
  
  struct command_t {
    void clear() { cmd_id = -1; cmd_number = -1; cmd_parameters = ""; timestamp = ""; }
    int cmd_id;
    int cmd_number;
    std::string cmd_name;
    std::string cmd_parameters;
    std::string timestamp;
    bool operator< ( const command_t & rhs ) const { return cmd_id < rhs.cmd_id; } 
  };
  
  struct timepoint_t {
    timepoint_t() : timepoint_id(-1) , epoch(-1) , start(0), stop(0) { } 
    void timeless() { timepoint_id = -1; epoch = -1; start = 0; stop = 0; } 
    bool none() const { return epoch == -1 && start == 0 && stop == 0; } 
    int timepoint_id; 
    int epoch;     
    uint64_t start;
    uint64_t stop;
    bool is_epoch() const { return epoch != -1; } 
    bool is_interval() const { return ! ( start == 0 && stop == 0 ) ; } 

    std::string print() const 
    { 
      std::stringstream ss;
      if ( epoch != -1 ) ss << epoch;
      else if ( is_interval() ) ss << start << "-" << stop;
      else ss << ".";
      return ss.str();
    }

    bool operator< ( const timepoint_t & rhs ) const 
    { 
      // 'global timepoints' come first
      if ( epoch == -1 )
	{
	  if ( rhs.epoch == -1 ) return false;
	  return true;
	}
      // sample-point based comparison
      if ( epoch == rhs.epoch )
	{
	  if ( epoch == -1 ) return false; // both 'global' timepoints
	  if ( start == rhs.start ) return stop < rhs.stop;
	  return start < rhs.start;	
	}    
      // epoch-based comparison
      return epoch < rhs.epoch; 
    }
  };
  
  
  struct strata_t
  {
    strata_t() { strata_id = -1; } 
    strata_t (const int strata_id) : strata_id( strata_id ) { } 
    int strata_id;

    // factor_id --> level_id // i.e. can only have 1 level of each factor
    std::map<factor_t,level_t> levels;
    
    void clear() { levels.clear(); }
    void insert( const level_t & l , const factor_t & f ) { levels[f] = l; }
    bool empty() const { return levels.size() == 0; } 
    bool operator<( const strata_t & rhs ) const
    {
      if ( levels.size() != rhs.levels.size() ) return levels.size() < rhs.levels.size();
      std::map<factor_t,level_t>::const_iterator ii = levels.begin();
      std::map<factor_t,level_t>::const_iterator jj = rhs.levels.begin();
      while ( ii != levels.end() ) 
	{
	  if ( ii->first < jj->first ) return true;
	  if ( jj->first < ii->first ) return false;

	  if ( ii->second < jj->second ) return true;
	  if ( jj->second < ii->second ) return false;

	  ++ii; ++jj;
	}
      return false;
    }
    
    std::string factor_string() const;
    std::string level_string() const;    
    std::string factor_level_string() const;
    int matches( const std::set<int> & cvars , const std::set<int> & rvars );
    std::string print() const;
    std::string print_nocmd() const; // skip if factor starts w/ _
  
    bool drop( int factor_id )
    {            
      std::map<factor_t,level_t> levels_copy = levels;
      levels.clear();      
      std::map<factor_t,level_t>::const_iterator ii = levels_copy.begin();
      while ( ii != levels_copy.end() )
	{
	  if ( ii->first.factor_id != factor_id ) levels[ ii->first ] = ii->second;
	  ++ii;
	}      
      return true;
    }
    
    std::string num_print() const
    {
      std::stringstream ss;
      ss << "[" << strata_id << "]";
      std::map<factor_t,level_t>::const_iterator ii = levels.begin();
      while ( ii != levels.end() )
	{
	  ss << "~" << ii->second.level_id << "/" << ii->first.factor_id ;
	  ++ii;
	}
      return ss.str();
    }
    
  };


// Old version that included timepoint information as part of the strata
  /* struct strata_t */
  /* { */
  /*   strata_t() { timepoint = -1; }  */
  /*   strata_t (const int strata_id) : strata_id( strata_id ) , timepoint( -1 ) { }  */
  /*   int strata_id; */

  /*   // factor_id --> level_id // i.e. can only have 1 level of each factor */
  /*   std::map<factor_t,level_t> levels; */
    
  /*   // a further stratification by time-point?, -1 means none, otherwise */
  /*   // this points of a timepoint_id in timepoints[] */

  /*   int timepoint;  */

  /*   void clear() { levels.clear(); timepoint = -1; } */
  /*   void insert( const level_t & l , const factor_t & f ) { levels[f] = l; } */
  /*   bool empty() const { return levels.size() == 0 && timepoint == -1; }  */
  /*   bool operator<( const strata_t & rhs ) const */
  /*   { */
  /*     if ( timepoint < rhs.timepoint ) return true; */
  /*     if ( rhs.timepoint < timepoint ) return false; */
  /*     if ( levels.size() != rhs.levels.size() ) return levels.size() < rhs.levels.size(); */
  /*     std::map<factor_t,level_t>::const_iterator ii = levels.begin(); */
  /*     std::map<factor_t,level_t>::const_iterator jj = rhs.levels.begin(); */
  /*     while ( ii != levels.end() )  */
  /* 	{ */
  /* 	  if ( ii->first < jj->first ) return true; */
  /* 	  if ( jj->first < ii->first ) return false; */

  /* 	  if ( ii->second < jj->second ) return true; */
  /* 	  if ( jj->second < ii->second ) return false; */

  /* 	  ++ii; ++jj; */
  /* 	} */
  /*     return false; */
  /*   } */
    
  /*   std::string factor_string() const; */
  /*   std::string level_string() const;     */
  /*   std::string factor_level_string() const; */
  /*   int matches( const std::set<int> & cvars , const std::set<int> & rvars ); */
  /*   std::string print() const; */
  /*   void drop_time() { timepoint = -1; }  */
  /*   bool drop( int factor_id ) */
  /*   {             */
  /*     std::map<factor_t,level_t> levels_copy = levels; */
  /*     levels.clear();       */
  /*     std::map<factor_t,level_t>::const_iterator ii = levels_copy.begin(); */
  /*     while ( ii != levels_copy.end() ) */
  /* 	{ */
  /* 	  if ( ii->first.factor_id != factor_id ) levels[ ii->first ] = ii->second; */
  /* 	  ++ii; */
  /* 	}       */
  /*     return true; */
  /*   } */
    
  /*   std::string num_print() const */
  /*   { */
  /*     std::stringstream ss; */
  /*     ss << "[" << strata_id << "]~" << timepoint; */
  /*     std::map<factor_t,level_t>::const_iterator ii = levels.begin(); */
  /*     while ( ii != levels.end() ) */
  /* 	{ */
  /* 	  ss << "~" << ii->second.level_id << "/" << ii->first.factor_id ; */
  /* 	  ++ii; */
  /* 	} */
  /*     return ss.str(); */
  /*   } */
    
  /* }; */



  
struct var_t
{
  int var_id;
  int cmd_id;
  std::string var_name;
  std::string var_label;
  bool operator<( const var_t & rhs ) const 
  { 
    if ( cmd_id == rhs.cmd_id ) return var_id < rhs.var_id; 
    return cmd_id < rhs.cmd_id;
  }
};

struct value_t
{ 
  value_t( const std::string & s ) : numeric(false) , integer(false), missing( false ) , s(s) { } 
  value_t( double d ) : numeric(true) , integer(false) , missing( false ), d(d) { } 
  value_t( int i ) : numeric(false) , integer(true) , missing(false) , i(i) { } 
  value_t() : missing(true) { } 

  // numeric   integer   missing
  // T         F         F             double
  // F         T         F             int
  // F         F         F             string
  // F         F         T             missing
  
  bool is_string() const { return ! ( numeric || integer || missing ); } 
  bool is_numeric() const { return numeric; } 
  bool any_number() const { return numeric || integer; } 
  bool is_integer() const { return integer; } 
  bool is_missing() const { return missing; } 
  
  bool numeric;
  bool integer;
  bool missing;

  double d;
  std::string s;
  int i;

  void set( const double d_ )        { numeric = true;  integer = false; d = d_; missing = false; } 
  void set( const int i_ )           { numeric = true;  integer = true;  i = i_; missing = false; } 
  void set( const std::string & s_ ) { numeric = false; integer = false; s = s_; missing = false; } 
  void set_missing() { missing = true; } 

  std::string str() const 
  {
    std::stringstream ss;
    if ( missing ) ss << "NA";    
    else if ( numeric ) ss << d;
    else if ( integer ) ss << i;
    else ss << s;
    return ss.str();
  }

  // not even sure we need this in any context....
  bool operator<( const value_t & rhs ) const
  {
    if ( missing     && ! rhs.missing ) return false;
    if ( rhs.missing && ! missing ) return true;
    if ( missing     && rhs.missing ) return false;

    if ( any_number() && ! rhs.any_number() ) return true;
    if ( (! any_number() ) && rhs.any_number() ) return false;    
    
    if ( is_string() && rhs.is_string() ) return s < rhs.s;
    
    if ( numeric && rhs.numeric ) return d < rhs.d;
    if ( integer && rhs.integer ) return i < rhs.i;
    if ( numeric && rhs.integer ) return d < rhs.i;
    return i < rhs.d;

  }
};


struct packet_t 
{ 
  int indiv_id; 
  int cmd_id; 
  int var_id; 
  int strata_id; 
  int timepoint_id; 
  value_t value; 
}; 

typedef std::vector<packet_t> packets_t;

//
// Database with internal cache
//

class StratOutDBase {  
  
 public:
  
  StratOutDBase()
    {
      
    }
  
  ~StratOutDBase()
    {      
      dettach();
    }

  bool attach( const std::string & name , bool readonly , writer_t * caller );

  bool dettach();
  
  // on first attaching a database, we read all factor/level/encodings, etc.
  // and save back to the calling writer_t
  void read_all( writer_t * );  
  
  bool init();
  bool release();
  bool attached() { return sql.is_open(); }
  void check_version();

  std::string name() const { return filename; } 

  //
  // Writers
  //

  
  // core value-adders
  bool add_value( int indiv_id , int cmd_id , int var_id , int strata_id , int timepoint_id , const value_t & value );
  
  int add_command( const std::string & cmd_name , const std::string & cmd_parameters );
  int add_individual( const std::string & indiv_name );
  int add_factor( const std::string & factor_name );
  int add_level( const std::string & level_name , const std::string & factor_name );
  int add_variable( const std::string & variable_name , const std::string & variable_label );
  int add_epoch_timepoint( int e );
  int add_interval_timepoint( const interval_t & interval );
  int add_stratum( const strata_t & stratum );

  
  //
  // Query data
  //
  
  std::set<std::string> variable_names();
  std::set<std::string> indiv_names();
  
  void fetch( int strata_id , int time_mode , 
	      packets_t * , 
	      std::set<int> * indiv_id = NULL , std::set<int> * cmd_id = NULL , std::set<int> * var_id = NULL );

  packets_t enumerate( int strata_id );

  packets_t dump_all();

  packets_t dump_indiv( const int indiv_id );
  
  std::map<int,std::set<int> > dump_vars_by_strata();

  std::map<int,int> count_strata();

  // match on var or cmd, irrespective of cmd:var pairing
  std::set<int> all_matching_vars( const std::set<std::string> & vars );
  std::set<int> all_matching_cmds( const std::set<std::string> & cmds );


  int num_values();

  
  //
  // Other
  //

  bool index();
  bool drop_index();
  void begin() { sql.begin_exclusive(); }
  void commit() { sql.commit(); }
  

  //
  // setters
  //
  
  indiv_t   insert_individual( const std::string & indiv_name , const std::string & filename );
  var_t     insert_variable( const std::string & var_name, const std::string & cmd_name , const std::string & var_label );
  timepoint_t   insert_epoch_timepoint( const int epoch );
  timepoint_t   insert_interval_timepoint( const interval_t & interval );
  factor_t  insert_factor( const std::string & fac_name, const bool is_numeric = false );
  level_t   insert_level( const std::string & level_name, const std::string & fac_name );
  level_t   insert_level( const std::string & level_name, const int factor_id );
  strata_t  insert_strata( const strata_t & );
  command_t insert_command( const std::string & cmd_name , int , const std::string & timedate , const std::string & cmd_param );
  bool      insert_value( const int indiv_id , const int cmd_id , const int variable_id , const int strata_id , const int tp_id , const value_t & x );

  // fetchers

  
 private:
  
  //
  // Database connection
  //
  
  SQL sql;

  std::string filename;
  
  
  //
  // Prepared queries
  //
  
  sqlite3_stmt * stmt_insert_indiv;
  sqlite3_stmt * stmt_insert_factor;   
  sqlite3_stmt * stmt_insert_level;    
  sqlite3_stmt * stmt_insert_stratum;  
  sqlite3_stmt * stmt_insert_command;  
  sqlite3_stmt * stmt_insert_variable; 
  sqlite3_stmt * stmt_insert_timepoint;    
  sqlite3_stmt * stmt_insert_value;    

  sqlite3_stmt * stmt_dump_factors;	      
  sqlite3_stmt * stmt_dump_levels;	      
  sqlite3_stmt * stmt_dump_strata;	      
  sqlite3_stmt * stmt_dump_variables;	      
  sqlite3_stmt * stmt_dump_individuals;	      
  sqlite3_stmt * stmt_dump_timepoints;	      
  sqlite3_stmt * stmt_dump_commands;	      

  // bool are 0/1 integer
  sqlite3_stmt * stmt_dump_int_datapoints;
  sqlite3_stmt * stmt_dump_dbl_datapoints;
  sqlite3_stmt * stmt_dump_txt_datapoints;

  sqlite3_stmt * stmt_count_values;
  sqlite3_stmt * stmt_lookup_value_by_null_strata;
  sqlite3_stmt * stmt_lookup_value_by_strata;
  sqlite3_stmt * stmt_lookup_value_by_strata_and_timepoint;
  
  sqlite3_stmt * stmt_enumerate;
  sqlite3_stmt * stmt_enumerate_null_strata;
  sqlite3_stmt * stmt_dump_vars_by_strata;
  sqlite3_stmt * stmt_count_strata;
  sqlite3_stmt * stmt_match_vars;
  sqlite3_stmt * stmt_match_cmds;


};


struct factor_t;
struct indiv_t;
struct level_t;
struct command_t;
struct timepoint_t;
struct value_t;
struct var_t;

class writer_t 
{
  
  friend struct level_t;
  friend struct strata_t;
  friend class StratOutDBase;
  

    
 public:

  
  //
  // database
  //

  writer_t() { dbless = true; retval = NULL; } 

  bool attach( const std::string & filename , bool readonly = false )
  {

    dbless = false; retval = NULL;

    db.attach( filename , readonly , this );

    //
    // Ensure that default strata is set as '1' baseline
    //
    
    if ( ! readonly ) 
      {
	strata_t baseline;
	int dummy = get_strata_id( baseline );
	if ( dummy != 1 ) Helper::halt( "internal problem with root strata_id != 1" );
      }

    return db.attached();
  }

  void nodb() 
  { 
    // in an in-memory DB to store factor information, etc, but then set to 'dbless' 
    close();
    attach( ":memory:" );
    dbless = true; 
    retval = NULL;
  } 
  
  void use_retval( retval_t * r ) 
  { 
    // in an in-memory DB to store factor information, etc, but then set to write to a retval 
    close();
    attach( ":memory:" );
    dbless = false;  
    retval = r;
  } 
  
  
  
  std::string name() const { return dbless ? "." : db.name(); } 

  void index() { db.index(); } 
  void drop_index() { db.drop_index(); } 
  void begin() { return db.begin(); } 
  void commit() { return db.commit(); }
  void read_all() { db.read_all(this); }
  bool attached() { return db.attached(); }
  void fetch( int strata_id, int time_mode, packets_t * packets, std::set<int> * i = NULL, std::set<int> * c = NULL, std::set<int> * v = NULL )
  { return db.fetch( strata_id , time_mode, packets, i, c, v) ; }  

  packets_t enumerate( int strata_id ) { return db.enumerate( strata_id ); }

  std::map<int,std::set<int> > dump_vars_by_strata() { return db.dump_vars_by_strata(); }

  std::map<int,int> count_strata() { return db.count_strata(); }

  std::set<int> all_matching_vars( const std::set<std::string> & vars ) { return db.all_matching_vars( vars ); }
  std::set<int> all_matching_cmds( const std::set<std::string> & cmds ) { return db.all_matching_cmds( cmds ); }

  // open db and send to a retval
  static retval_t dump_to_retval( const std::string & dbname , const std::set<std::string> * = NULL , std::vector<std::string> * ids = NULL );

  bool close() 
  { 
    if ( ! attached() ) return false;
    clear(); 
    db.dettach();
    return true; 
  } 

  
  ~writer_t() { close(); } 
  

  //
  // Definitions
  //

    // Define factor types (not set current level)
  
  bool numeric_factor( const std::string & fac_name )
  {
    if ( factors_idmap.find( fac_name ) == factors_idmap.end() )
      {
	factor_t factor = db.insert_factor( fac_name , 1 );  // 1 -> numeric factor
	factors_idmap[ fac_name ] = factor.factor_id;
	factors[ factor.factor_id ] = factor;
      }
    return true;
  }
  
  bool string_factor( const std::string & fac_name )
  {
    if ( factors_idmap.find( fac_name ) == factors_idmap.end() )
      {
	factor_t factor = db.insert_factor( fac_name , 0 ); // 0 -> string factor  
	factors_idmap[ fac_name ] = factor.factor_id;
	factors[ factor.factor_id ] = factor;
      }
    return true;
  }

  bool var( const std::string & var_name , const std::string & var_label )
  {
    // use 'command.var' as the unique identifier
    std::string var_key = curr_command.cmd_name + ":" + var_name;
    if ( variables_idmap.find( var_key ) == variables_idmap.end() )
      {
	var_t var = db.insert_variable( var_name , curr_command.cmd_name , var_label );
	variables_idmap[ var_key ] = var.var_id;
	variables[ var.var_id ] = var;
      }
    return true;
  }
  

  //
  // writers
  //
  
  // Current command
  
  bool cmd( const std::string & cmd_name , const int cmd_number , const std::string & param )
  {
    
    // use 'command.number' as unqiue identifier
    
    std::string command_key = cmd_name + "." + Helper::int2str( cmd_number );

    // cached?
    if ( commands_idmap.find( command_key ) != commands_idmap.end() )
      curr_command = commands[ commands_idmap[ command_key ] ];
    else
      {
	// add to DB
	curr_command = db.insert_command( cmd_name , cmd_number , timestamp() , param );		
	// add to cache
	commands_idmap[ command_key ] = curr_command.cmd_id;
	commands[ curr_command.cmd_id ] = curr_command;
      }

    return true;

  }

  // Current individual

  bool id( const std::string & indiv_name , const std::string & file_name )
  {
    if ( individuals_idmap.find( indiv_name ) != individuals_idmap.end() )
      curr_indiv = individuals[ individuals_idmap[ indiv_name ] ];
    else
      {
	curr_indiv = db.insert_individual( indiv_name , file_name );
	individuals_idmap[ indiv_name ] = curr_indiv.indiv_id;
	individuals[ curr_indiv.indiv_id ] = curr_indiv;
      }
    return true;
  }
  
  // Current tag:  this is a LEVEL/FACTOR
  
  bool tag( const std::string & lvl_name , const std::string & fac_name )
  {

    if      ( fac_name == "." ) unlevel();
    else if ( lvl_name == "." ) unlevel( fac_name );
    else 
      {
	// ensure we have this set
	string_factor( fac_name );
	level( lvl_name , fac_name );
      }
    return true;
  }

  // Current factor
  bool level( const int level_name , const std::string & factor_name )
  {
    return level( Helper::int2str( level_name ) , factor_name );
  }

  bool level( const double level_name , const std::string & factor_name )
  {
    return level( Helper::dbl2str( level_name ) , factor_name );
  }
  
  bool level( const std::string & level_name , const std::string & factor_name )
  {

    // add factor (as string by default) if it doesn't already exist
    if ( factors_idmap.find( factor_name ) == factors_idmap.end() ) 
      string_factor( factor_name );
    
    factor_t factor = factors[ factors_idmap[ factor_name ] ];
    
    // for level, use level.factor as the lookup key
    std::string level_key = level_name + "." + factor_name ;
    
    // cached?
    if ( levels_idmap.find( level_key ) == levels_idmap.end() )
      {
	// if not, add to DB,
	level_t level = db.insert_level( level_name , factor.factor_id );
	// and then cache
	levels_idmap[ level_key ] = level.level_id;
	levels[ level.level_id ] = level;
      }

    // fetch from cache
    level_t level = levels[ levels_idmap[ level_key ] ];

    // swap/add to current strata
    curr_strata.insert( level , factor );
        
    return true;
  }
  
  bool unlevel( const std::string & factor_name )
  {
    // drop 'factor_name' from current stratification (curr_strata)

    // never added / no need to drop
    if ( factors_idmap.find( factor_name ) == factors_idmap.end() ) return false;

    // drop this level/factor from current strata
    curr_strata.drop( factors_idmap[ factor_name ] );
        
    return true;
  }
  
  bool unlevel() 
  {
    // set curr_strata to 'empty' 
    curr_strata.clear();
    return true;
  }
  
  int get_strata_id( const strata_t & s )
  {
    // when being set up, we always enter a first 'default' baseline
    // stratum this has level code of '0'
    
    // if this is here, it will have an ID
    if ( strata_idmap.find( s ) != strata_idmap.end() ) 
      {
	//std::cout << "found existing " << s.strata_id << " " << strata_idmap[ s ] << "\n";	
	return strata_idmap[ s ];
      }
    
    // if not, add to DB and track ID
    strata_t new_strata = db.insert_strata( s );
    //std::cout << "new id = " << new_strata.strata_id;
    strata_idmap[ new_strata ] = new_strata.strata_id;
    strata[ new_strata.strata_id ] = new_strata;
    return new_strata.strata_id;
  }

  // Time-points: epochs or intervals
  
  bool epoch( const int e )
  {

    if ( e == -1 ) 
      {
	curr_timepoint.timeless();
	return true;
      }
    
    std::string tp_key = Helper::int2str(e) + ":"; 
    if ( timepoints_idmap.find( tp_key ) != timepoints_idmap.end() )
      {
	curr_timepoint = timepoints[ timepoints_idmap[ tp_key ] ];	
      }
    else
      {
	// add new TP before attaching to current strata
	curr_timepoint = db.insert_epoch_timepoint( e );
       	timepoints_idmap[ tp_key ] = curr_timepoint.timepoint_id;
	timepoints[ curr_timepoint.timepoint_id ] = curr_timepoint;	
      }
    
    // add as a factor ('E') to the current strata -- although we only ever add a single 
    // dummy layer to this factor.  The lookup is done via the datapoints and timepoints 
    // tables
    
    level( "." , globals::epoch_strat );

    return true;
  }
  
  bool interval( const interval_t & interval )
  {
    
    if ( interval.start == 0 && interval.stop == 0 )
      {
	curr_timepoint.timeless();
	return true;
      }

    std::string tp_key = ":" + Helper::int2str(interval.start) + "-" + Helper::int2str( interval.stop ) ; 
    if ( timepoints_idmap.find( tp_key ) != timepoints_idmap.end() )
      curr_timepoint = timepoints[ timepoints_idmap[ tp_key ] ];
    else
      {
	curr_timepoint = db.insert_interval_timepoint( interval );
	timepoints_idmap[ tp_key ] = curr_timepoint.timepoint_id;
	timepoints[ curr_timepoint.timepoint_id ] = curr_timepoint;
      }
    
    // add as a factor ('T') to the current strata -- although we only ever add a single 
    // dummy layer to this factor.  The lookup is done via the datapoints and timepoints 
    // tables
    
    level( "." , globals::time_strat );

    return true;
  }
  
  bool uninterval() { unlevel( globals::epoch_strat ); return timeless(); } 

  bool unepoch() { unlevel( globals::time_strat ); return timeless(); } 
  
  bool timeless() 
  {
    unlevel( globals::epoch_strat );
    unlevel( globals::time_strat );    
    // set to timeless
    curr_timepoint.timeless();
    return true;
  }
  
  //
  // Value (to DB or stdout)
  //

  bool value( const std::string & var_name , double d , const std::string & desc = "" )
  {    
    if ( retval != NULL ) return to_retval( var_name , d );
    else if ( dbless ) return to_stdout( var_name , value_t( d ) ) ;
    if ( desc != "" ) var( var_name , desc );
    return value( var_name , value_t( d ) ) ;
  }

  bool value( const std::string & var_name , int i , const std::string & desc = "" ) 
  { 
    if ( retval != NULL ) return to_retval( var_name , i ); 
    else if ( dbless ) return to_stdout( var_name , value_t( i ) ) ; 
    if ( desc != "" ) var( var_name , desc ); 
    return value( var_name , value_t( i ) ) ; 
  } 
  
  bool value( const std::string & var_name , const std::string & s , const std::string & desc = "" )
  {
    if ( retval != NULL ) return to_retval( var_name , s );
    if ( dbless ) return to_stdout( var_name , value_t( s ) ); 
    if ( desc != "" ) var( var_name , desc );
    return value( var_name , value_t( s ) ) ;
  }
  
  bool missing_value( const std::string & var_name , const std::string & desc = "" )
  {
    if ( retval != NULL ) return to_retval( var_name ); // missing value 
    if ( dbless ) return to_stdout( var_name , value_t() ); 
    if ( desc != "" ) var( var_name , desc );
    return value( var_name , value_t() );
  }

  bool value( const std::string & var_name , const value_t & x )
  {

    // this should never be called in retval mode, but just in case... 
    if ( retval != NULL ) Helper::halt( "internal error in value(), should not get here" );

    if ( dbless ) return to_stdout( var_name , x );

    // use 'command.var' as the unique identifier                                                                                                     
    std::string var_key = curr_command.cmd_name + ":" + var_name;
    
    // should be already here, but incase it is not
    if ( variables_idmap.find( var_key ) == variables_idmap.end() )
      {
	var_t var = db.insert_variable( var_name , curr_command.cmd_name , "." );
        variables_idmap[ var_key ] = var.var_id;
        variables[ var.var_id ] = var;
      }      

    // check curr_strata is registered; add to DB if not
    curr_strata.strata_id = get_strata_id( curr_strata );    

    // store value    
    db.insert_value( curr_indiv.indiv_id , 
		     curr_command.cmd_id , 	
		     variables_idmap[ var_key ] , 		     
		     curr_strata.empty() ? -1 : curr_strata.strata_id , 
		     curr_timepoint.none() ? -1 : curr_timepoint.timepoint_id , 
		     x );
    
    return true;
  }


  bool to_stdout( const std::string & var_name , const value_t & x )  
  {
    std::cout << curr_indiv.indiv_name << "\t"
	      << curr_command.cmd_name ;
    
    if ( curr_strata.empty() ) std::cout << "\t.";
    else std::cout << "\t" << curr_strata.print_nocmd() ;
    
    if ( curr_timepoint.none() ) std::cout << "\t.";
    else std::cout << "\t" << curr_timepoint.print();
    
    std::cout << "\t" << var_name 
	      << "\t" << x.str() 
	      << "\n";
      
    return true;
  }

  
  
  bool to_retval( const std::string & var_name , double d )
  {

    retval->add( curr_indiv.indiv_name,
		 retval_cmd_t( curr_command.cmd_name ) , 
		 retval_factor_t( curr_strata , curr_timepoint ) ,
		 retval_var_t( var_name ) , 
		 retval_strata_t( curr_strata , curr_timepoint ) ,
		 d );
    
    return true;
  }


  bool to_retval( const std::string & var_name , int i )
  {

    retval->add( curr_indiv.indiv_name,
		 retval_cmd_t( curr_command.cmd_name ) , 
		 retval_factor_t( curr_strata , curr_timepoint ) ,
		 retval_var_t( var_name ) , 
		 retval_strata_t( curr_strata , curr_timepoint ) ,
		 i);

    return true;
  }



  bool to_retval( const std::string & var_name , const std::string & s  )
  {
    
    retval->add( curr_indiv.indiv_name,
		 retval_cmd_t( curr_command.cmd_name ) ,
		 retval_factor_t( curr_strata , curr_timepoint ) ,
		 retval_var_t( var_name ) , 
		 retval_strata_t( curr_strata , curr_timepoint ) ,
		 s );
    
    return true;
  }
  

  bool to_retval( const std::string & var_name )
  {
    // use special string code for missing data 'NA' 
    retval->add( curr_indiv.indiv_name,
		 retval_cmd_t( curr_command.cmd_name ) ,
		 retval_factor_t( curr_strata , curr_timepoint ) ,
		 retval_var_t( var_name ) ,
		 retval_strata_t( curr_strata , curr_timepoint ) ,
		 "NA" );
    
    return true;
  }






  
  //
  // readers
  //
  
  int num_factors() const { return factors.size(); } 
  int num_levels() const { return levels.size(); } 
  int num_variables() const { return variables.size(); }
  int num_strata() const { return strata.size(); }
  int num_commands() const { return commands.size(); }
  int num_individuals() const { return individuals.size(); } 
  int num_timepoints() const { return timepoints.size(); }
  int num_values() { return db.num_values(); } 

  // list all variable names
  std::set<std::string> variable_names() { return db.variable_names(); std::set<std::string> dummy; }
  std::set<std::string> indiv_names() { return db.indiv_names(); std::set<std::string> dummy; }
  
  // caches
  
  std::map<int,factor_t>  factors;
  std::map<int,level_t>   levels;
  std::map<int,var_t>     variables;
  std::map<int,strata_t>  strata;
  std::map<int,indiv_t>   individuals;
  std::map<int,command_t> commands;
  std::map<int,timepoint_t>   timepoints;
  
  // lookup of ID based on textual identifier
  std::map<std::string,int> factors_idmap;
  std::map<std::string,int> levels_idmap;
  std::map<std::string,int> variables_idmap;
  std::map<std::string,int> individuals_idmap;
  std::map<std::string,int> timepoints_idmap;
  std::map<strata_t,int>    strata_idmap;
  std::map<std::string,int> commands_idmap;

  void clear() 
  {
    factors.clear();     factors_idmap.clear();
    levels.clear();      levels_idmap.clear();
    variables.clear();   variables_idmap.clear();
    individuals.clear(); individuals_idmap.clear();
    commands.clear();    commands_idmap.clear();
    timepoints.clear();  timepoints_idmap.clear();
    strata.clear();      strata_idmap.clear();
    
    curr_indiv.clear();
    curr_strata.clear();
    curr_timepoint.timeless();
    curr_command.clear();
  }

  
 private:

  // primary data-store
  
  StratOutDBase db; 
  
  // write to std::cout, instead of to a DB

  bool dbless;

  // write to a retval_t, instead of a DB

  retval_t * retval;
  
  // Used when reading (do we need this??)
  
  // primary data-store:  indiv -> strat -> timepoint -> cmd -> variable -> value 
  //  std::map<indiv_t, std::map<strata_t, std::map<timepoint_t, std::map<command_t, std::map<var_t, value_t> > > > > data;
  
  // 'current' state when writing
  
  indiv_t       curr_indiv;
  command_t     curr_command;
  strata_t      curr_strata;
  timepoint_t   curr_timepoint;  
  
 
  
  // helper functions
  std::string timestamp()
    {
      time_t curr=time(0);
      std::string tdstamp = (std::string)ctime(&curr);
      // this ends with /n, but check just in case
      if ( tdstamp[ tdstamp.size() - 1 ] != '\n' ) tdstamp += '\n';
      return tdstamp.substr( 0 , tdstamp.size() - 1 ) ;
    }
    
};


#endif
