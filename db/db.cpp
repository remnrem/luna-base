
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

#include "db.h"

#include <iostream>
#include <set>
#include <utility>

extern writer_t writer;

std::string strata_t::factor_string() const
{
  if ( levels.size() == 0 ) return ".";
  std::string r;
  std::map<factor_t,level_t>::const_iterator ii = levels.begin();
  while ( ii != levels.end() )
    {
      if ( r != "" ) r += "/";	
      r += ii->first.factor_name;
      ++ii;
    }    
  return r;
}

std::string strata_t::level_string() const
{
  if ( levels.size() == 0 ) return ".";
  std::string r;
  std::map<factor_t,level_t>::const_iterator ii = levels.begin();
  while ( ii != levels.end() )
    {
      if ( r != "" ) r += "/";
      r += ii->second.level_name;	
      ++ii;
    }
  return r;
}

// for header variable names
std::string strata_t::factor_level_string() const
{
  if ( levels.size() == 0 ) return "";
  std::string r;
  std::map<factor_t,level_t>::const_iterator ii = levels.begin();
  while ( ii != levels.end() )
    {
      if ( r != "" ) r += ".";
      //std::string factor_name = writer.factors[ ii->factor_id ].factor_name ; 
      r += ii->first.factor_name + "_" + ii->second.level_name;
      ++ii;
    }    
  return r;
}

int strata_t::matches( const std::set<int> & cvars , 
		       const std::set<int> & rvars )            
{
  // +1 exact match (contains cvar and rvar and nothing else)
  //  0 doesn't contain all cvar and rvar
  // -1 contains cvar/rvar but also others 
  
  bool additional = false;
  int match = 0 ;
  
  std::map<factor_t,level_t>::const_iterator ll = levels.begin();
  while ( ll != levels.end() )
    { 
      if ( cvars.find( ll->first.factor_id ) == cvars.end() &&
	   rvars.find( ll->first.factor_id ) == rvars.end() )
	additional = true;
      else 
	++match;      
      ++ll;
    }

  if ( match < rvars.size() + cvars.size() ) return 0;
  return additional ? -1 : +1 ;
}

std::string strata_t::print() const
{
  if ( levels.size() == 0 ) return ".";

  std::stringstream ss;
  std::map<factor_t,level_t>::const_iterator aa = levels.begin();
  bool printed = false;
  while ( aa != levels.end() )
    {
      // skip if epoch/time-point
      if ( aa->first.factor_name == globals::epoch_strat || 
	   aa->first.factor_name == globals::time_strat ) { ++aa; continue; }  

      if ( printed ) ss << ";";
      ss << aa->first.factor_name << "/" << aa->second.level_name ; 
      printed = true;
      ++aa;
    }
  std::string rstr = ss.str();
  if ( rstr == "" ) return ".";
  return ss.str();
}    


bool StratOutDBase::attach( const std::string & n , bool readonly )
{
  
  if ( attached() ) dettach();
  
  if ( n == "-" || n == "." ) { dettach(); return false; } 

  sql.open(n); 
  
  sql.synchronous(false);
  
  filename = n;
  
  //
  // Tables
  //


  // Factors

  sql.query(" CREATE TABLE IF NOT EXISTS factors("
            "   factor_id   INTEGER PRIMARY KEY , "
	    "   factor_name VARCHAR(20) NOT NULL , "
	    "   is_numeric  INTEGER ) ; " ); 

  // Levels

  sql.query(" CREATE TABLE IF NOT EXISTS levels("
            "   level_id   INTEGER PRIMARY KEY , "
	    "   factor_id  INTEGER NOT NULL , "
	    "   level_name VARCHAR(20) ) ; " ); 


  // Strata (a specific combintion of factor levels, and optionally a time-point)
  // a -ve level_id means this is a timepoint_id (i.e. special case of a factor)
  // otherwise, lookup in the level/factor tables

  sql.query(" CREATE TABLE IF NOT EXISTS strata("
            "   strata_id    INTEGER NOT NULL , "
	    "   level_id     INTEGER NOT NULL ); " );
  
  
  // Variables
  
  sql.query(" CREATE TABLE IF NOT EXISTS variables("
            "   variable_id    INTEGER PRIMARY KEY , "
	    "   variable_name  VARCHAR(20) NOT NULL , "
	    "   command_name   VARCHAR(20) , "
	    "   variable_label VARCHAR(20) ); " ); 
  
  // Individuals
  
  sql.query(" CREATE TABLE IF NOT EXISTS individuals("
	    "   indiv_id    INTEGER PRIMARY KEY , "
	    "   indiv_name  VARCHAR(20) NOT NULL , "
	    "   file_name   VARCHAR(20) ); " ); 
  
  // Commands

  sql.query(" CREATE TABLE IF NOT EXISTS commands("
	    "   cmd_id          INTEGER PRIMARY KEY , "
	    "   cmd_name        VARCHAR(20) NOT NULL , "
	    "   cmd_number      INTEGER NOT NULL , "
	    "   cmd_timestamp   VARCHAR(20) NOT NULL , "
	    "   cmd_parameters  VARCHAR(20)  ); " ); 

  // Timepoints

  sql.query(" CREATE TABLE IF NOT EXISTS timepoints("
	    "   timepoint_id      INTEGER PRIMARY KEY , "
	    "   epoch         INTEGER , "
	    "   start         UNSIGNED BIG INT , "
	    "   stop          UNSIGNED BIG INT ); " );

  // Values

  sql.query(" CREATE TABLE IF NOT EXISTS datapoints("
	    "   indiv_id      INTEGER NOT NULL , "
	    "   cmd_id        INTEGER NOT NULL , "
	    "   variable_id   INTEGER NOT NULL , "	    
	    "   strata_id     INTEGER , "
	    "   timepoint_id  INTEGER , "
	    "   value         NUMERIC ); " );


  //
  // ensure index is dropped when writing
  //

  if ( ! readonly ) 
    drop_index();  
  
  //
  // Prepare some key queries
  //

  init();


  //
  // get any existing encodings
  //

  read_all();

       
  return true;

}



void StratOutDBase::read_all()
{


  //
  // Individuals
  //

  while ( sql.step( stmt_dump_individuals ) )
    {
      indiv_t indiv;
      indiv.indiv_id = sql.get_int( stmt_dump_individuals , 0 );
      indiv.indiv_name = sql.get_text( stmt_dump_individuals , 1 );
      writer.individuals[ indiv.indiv_id ] = indiv;
      writer.individuals_idmap[ indiv.indiv_name ] = indiv.indiv_id;
    }
  sql.reset( stmt_dump_individuals );
  
  
  //
  // Commands
  //

  while ( sql.step( stmt_dump_commands ) )
    {
      command_t cmd;
      cmd.cmd_id = sql.get_int( stmt_dump_commands , 0 );
      cmd.cmd_name = sql.get_text( stmt_dump_commands , 1 );
      cmd.cmd_number = sql.get_int( stmt_dump_commands , 2 );
      cmd.timestamp = sql.get_text( stmt_dump_commands , 3 );
      cmd.cmd_parameters = sql.get_text( stmt_dump_commands , 4 );
      writer.commands[ cmd.cmd_id ] = cmd;
      writer.commands_idmap[ cmd.cmd_name ] = cmd.cmd_id;
    }
  sql.reset( stmt_dump_commands );


  //
  // Factors
  //

  while ( sql.step( stmt_dump_factors ) )
    {
      factor_t factor;
      factor.factor_id = sql.get_int( stmt_dump_factors , 0 );
      factor.factor_name = sql.get_text( stmt_dump_factors , 1 );
      factor.is_numeric = sql.get_int( stmt_dump_factors , 2 ) == 1;
      writer.factors[ factor.factor_id ] = factor;
      writer.factors_idmap[ factor.factor_name ] = factor.factor_id;
    }
  sql.reset( stmt_dump_factors );


  //
  // Levels
  //

  while ( sql.step( stmt_dump_levels ) )
    {
      level_t level;
      level.level_id = sql.get_int( stmt_dump_levels , 0 );
      level.factor_id = sql.get_int( stmt_dump_levels , 1 );
      level.level_name = sql.get_text( stmt_dump_levels , 2 );
      writer.levels[ level.level_id ] = level;
      if ( writer.factors.find( level.factor_id ) == writer.factors.end() )
	Helper::halt( "internal error, undefined factor" );
      std::string level_key = level.level_name 
	+ "." + writer.factors[ level.factor_id ].factor_name ;
      writer.levels_idmap[ level_key ] = level.level_id;
    }
  sql.reset( stmt_dump_levels );


  //
  // Strata
  //

  while ( sql.step( stmt_dump_strata ) )
    {
            
      int strata_id = sql.get_int( stmt_dump_strata , 0 );
      int level_id  = sql.get_int( stmt_dump_strata , 1 );
      
      // note, dummy code of level_id == 0 means this is
      // the root strata, so do not add any levels in that 
      // case

      if ( level_id == 0 ) 
	{
	  // add root if needed
	  if ( writer.strata.find( strata_id ) == writer.strata.end() )
	    {
	      strata_t s;
	      s.strata_id = strata_id;
	      writer.strata[ strata_id ] = s;	
	    }
	  
	}
      
      else
	{
	  level_t level; 
	  factor_t factor; 
	  
	  level = writer.levels[ level_id ]; 
	  factor = writer.factors[ level.factor_id ];	    
	  
	  // existing strata?
	  if ( writer.strata.find( strata_id ) != writer.strata.end() )
	    {
	      writer.strata[ strata_id ].insert( level , factor );
	    }
	  else
	    {
	      strata_t s;
	      s.strata_id = strata_id;
	      
	      s.insert( level , factor );
	      
	      writer.strata[ strata_id ] = s;	
	    }
	}

    }
  sql.reset( stmt_dump_strata );


  //
  // Update strata idmap
  //
  
  writer.strata_idmap.clear();
  std::map<int,strata_t>::const_iterator ss = writer.strata.begin();
  while ( ss != writer.strata.end() )
    {
      writer.strata_idmap[ ss->second ] = ss->first;
      ++ss;
    }


  //
  // Timepoints
  //

  while ( sql.step( stmt_dump_timepoints ) )
    {
      timepoint_t timepoint;
      timepoint.timepoint_id = sql.get_int( stmt_dump_timepoints , 0 );

      bool has_epoch = ! sql.is_null( stmt_dump_timepoints , 1 );
      timepoint.epoch = has_epoch ? sql.get_int( stmt_dump_timepoints , 1 ) : -1 ;

      bool has_interval = ! sql.is_null( stmt_dump_timepoints , 2 );      
      timepoint.start = has_interval ? sql.get_int64( stmt_dump_timepoints , 2 ) : 0 ; 
      timepoint.stop  = has_interval ? sql.get_int64( stmt_dump_timepoints , 3 ) : 0 ; 
      
      std::string tp_key = ( has_epoch ? Helper::int2str( timepoint.epoch ) : "" ) + ":" 
	+ ( has_interval ? Helper::int2str( timepoint.start ) + "-" + Helper::int2str( timepoint.stop ) : "" ); 
      
      writer.timepoints[ timepoint.timepoint_id ] = timepoint;
      writer.timepoints_idmap[ tp_key ] = timepoint.timepoint_id;

    }

  sql.reset( stmt_dump_timepoints );


  //
  // Variables
  //

  while ( sql.step( stmt_dump_variables ) )
    {
      var_t var;
      var.var_id = sql.get_int( stmt_dump_variables , 0 );
      var.var_name = sql.get_text( stmt_dump_variables , 1 );
      var.var_label = sql.get_text( stmt_dump_variables , 3 );

      std::string command_name = sql.get_text( stmt_dump_variables , 2 );
      if ( writer.commands_idmap.find( command_name ) != writer.commands_idmap.end() ) 
	var.cmd_id = writer.commands_idmap[ command_name ];

      writer.variables[ var.var_id ] = var;
      writer.variables_idmap[ command_name + ":" + var.var_name ] = var.var_id;
    }
  sql.reset( stmt_dump_variables );
  
}

bool StratOutDBase::dettach()
{
  release();
  sql.close();
  return true;
}

std::set<std::string> StratOutDBase::variable_names()
{
  std::set<std::string> s;
  while ( sql.step( stmt_dump_variables ) )
    s.insert( sql.get_text( stmt_dump_variables , 1 ) );
  sql.reset( stmt_dump_variables );
  return s;
}

std::set<std::string> StratOutDBase::indiv_names()
{
  std::set<std::string> s;
  while ( sql.step( stmt_dump_individuals ) )
    s.insert( sql.get_text( stmt_dump_individuals , 1 ) );
  sql.reset( stmt_dump_individuals );
  return s;
}

bool StratOutDBase::init()
{

  // dumpers: from all except /values/

  stmt_dump_factors = sql.prepare( "SELECT * FROM factors;" );
  stmt_dump_levels = sql.prepare( "SELECT * FROM levels;" );
  stmt_dump_strata = sql.prepare( "SELECT * FROM strata;" );
  stmt_dump_variables = sql.prepare( "SELECT * FROM variables;" );
  stmt_dump_individuals = sql.prepare( "SELECT * FROM individuals;" );
  stmt_dump_timepoints = sql.prepare( "SELECT * FROM timepoints;" );
  stmt_dump_commands = sql.prepare( "SELECT * FROM commands;" );

  // queries
  stmt_count_values = sql.prepare( "SELECT count(1) FROM datapoints;" );
  stmt_lookup_value_by_null_strata = sql.prepare( "SELECT * FROM datapoints WHERE timepoint_id IS NULL AND strata_id IS NULL ; " );
  stmt_lookup_value_by_strata = sql.prepare( "SELECT * FROM datapoints WHERE timepoint_id IS NULL AND strata_id == :strata_id; " );
  stmt_lookup_value_by_strata_and_timepoint = sql.prepare( "SELECT * FROM datapoints WHERE timepoint_id IS NOT NULL AND strata_id == :strata_id; " );

  stmt_enumerate = sql.prepare( "SELECT indiv_id,cmd_id,variable_id,count(*) FROM datapoints WHERE strata_id == :strata_id GROUP BY indiv_id,cmd_id, variable_id;");
  stmt_enumerate_null_strata = sql.prepare( "SELECT indiv_id,cmd_id,variable_id,count(*) FROM datapoints WHERE strata_id IS NULL GROUP BY indiv_id,cmd_id, variable_id;");


  stmt_dump_vars_by_strata = sql.prepare( "SELECT DISTINCT strata_id , variable_id FROM datapoints;" );

  stmt_count_strata = sql.prepare( "SELECT strata_id,count(*) FROM datapoints GROUP BY strata_id ;" );

  stmt_match_vars = sql.prepare( "SELECT variable_id,variable_name FROM variables;" );
  stmt_match_cmds = sql.prepare( "SELECT cmd_id,cmd_name FROM commands;" );

  // inserters

  stmt_insert_indiv     = sql.prepare(" INSERT OR REPLACE INTO individuals ( indiv_name , file_name ) values( :indiv_name , :file_name ) ; ");
  stmt_insert_variable  = sql.prepare(" INSERT OR REPLACE INTO variables ( variable_name , command_name , variable_label ) values( :var_name, :cmd_name , :var_label ) ; ");
  stmt_insert_command   = sql.prepare(" INSERT OR REPLACE INTO commands ( cmd_name , cmd_number, cmd_timestamp, cmd_parameters ) "
				      " values( :cmd_name , :cmd_number, :cmd_timestamp, :cmd_parameters ) ; ");
  stmt_insert_factor    = sql.prepare(" INSERT OR REPLACE INTO factors ( factor_name , is_numeric ) values( :fac_name, :is_num ) ; ");
  stmt_insert_level     = sql.prepare(" INSERT OR REPLACE INTO levels ( level_name , factor_id ) values( :level_name, :fac_id ) ; ");
  stmt_insert_stratum   = sql.prepare(" INSERT OR REPLACE INTO strata ( strata_id , level_id ) values( :strata_id, :level_id ) ; ");
  stmt_insert_timepoint   = sql.prepare(" INSERT OR REPLACE INTO timepoints ( epoch , start , stop ) values( :epoch , :start , :stop ) ; ");  
  stmt_insert_value   = sql.prepare(" INSERT OR REPLACE INTO datapoints ( indiv_id, cmd_id, variable_id, strata_id, timepoint_id, value ) "
				    " values( :indiv_id, :cmd_id, :variable_id, :strata_id, :timepoint_id, :value ) ; ");  

  return true;
}

bool StratOutDBase::release()
{

  sql.finalise( stmt_insert_indiv );
  sql.finalise( stmt_insert_factor );   
  sql.finalise( stmt_insert_level);    
  sql.finalise( stmt_insert_stratum);  
  sql.finalise( stmt_insert_command);  
  sql.finalise( stmt_insert_variable); 
  sql.finalise( stmt_insert_timepoint);    
  sql.finalise( stmt_insert_value);    

  sql.finalise( stmt_dump_factors);	      
  sql.finalise( stmt_dump_levels);	      
  sql.finalise( stmt_dump_strata);	      
  sql.finalise( stmt_dump_variables);	      
  sql.finalise( stmt_dump_individuals);	      
  sql.finalise( stmt_dump_timepoints);	      
  sql.finalise( stmt_dump_commands);	      

  sql.finalise( stmt_lookup_value_by_strata);
  sql.finalise( stmt_lookup_value_by_strata_and_timepoint);
  sql.finalise( stmt_count_values );
  return true;
}


bool StratOutDBase::index()
{
  if ( ! attached() ) return false;
  sql.query( "CREATE INDEX IF NOT EXISTS vIndex ON datapoints(strata_id); " );
  // schema changed, so update prepared queries
  release();
  init();  
  return true;
}


bool StratOutDBase::drop_index()
{
  if ( ! attached() ) return false;
  sql.query( "DROP INDEX IF EXISTS vIndex;" );
  // schema changed, so update prepared queries
  release();
  init();
  
  return true;
}



indiv_t   StratOutDBase::insert_individual( const std::string & indiv_name , const std::string & file_name )
{
  
  sql.bind_text( stmt_insert_indiv , ":indiv_name" , indiv_name );
  sql.bind_text( stmt_insert_indiv , ":file_name" , file_name );
  sql.step( stmt_insert_indiv );
  sql.reset( stmt_insert_indiv );

  indiv_t indiv;
  indiv.indiv_name = indiv_name;
  indiv.file_name = file_name;
  indiv.indiv_id = sql.last_insert_rowid();
  return indiv;
}

var_t     StratOutDBase::insert_variable( const std::string & var_name, const std::string & cmd_name , const std::string & var_label )
{
  sql.bind_text( stmt_insert_variable , ":var_name" , var_name );
  sql.bind_text( stmt_insert_variable , ":cmd_name" , cmd_name );
  sql.bind_text( stmt_insert_variable , ":var_label" , var_label );
  sql.step( stmt_insert_variable );
  sql.reset( stmt_insert_variable );

  var_t var;
  var.var_id = sql.last_insert_rowid();
  var.var_name = var_name;
  var.var_label = var_label;
  return var;
}

timepoint_t   StratOutDBase::insert_epoch_timepoint( const int epoch )
{
  sql.bind_int( stmt_insert_timepoint , ":epoch" , epoch );
  sql.bind_null( stmt_insert_timepoint , ":start" );
  sql.bind_null( stmt_insert_timepoint , ":stop" );
  sql.step( stmt_insert_timepoint );
  sql.reset( stmt_insert_timepoint );

  timepoint_t timepoint;
  timepoint.timepoint_id = sql.last_insert_rowid();
  timepoint.epoch = epoch;
  return timepoint;
}


timepoint_t   StratOutDBase::insert_interval_timepoint( const interval_t & interval )
{
  sql.bind_null( stmt_insert_timepoint , ":epoch"  );
  sql.bind_int64( stmt_insert_timepoint , ":start" , interval.start );
  sql.bind_int64( stmt_insert_timepoint , ":stop" , interval.stop );
  sql.step( stmt_insert_timepoint );
  sql.reset( stmt_insert_timepoint );

  timepoint_t timepoint;
  timepoint.timepoint_id = sql.last_insert_rowid();
  timepoint.epoch = -1;
  timepoint.start = interval.start;
  timepoint.stop  = interval.stop;
  return timepoint;

}


factor_t  StratOutDBase::insert_factor( const std::string & fac_name, const bool is_numeric )
{
  sql.bind_text( stmt_insert_factor , ":fac_name" , fac_name );
  sql.bind_int( stmt_insert_factor , ":is_num" , is_numeric );
  sql.step( stmt_insert_factor );
  sql.reset( stmt_insert_factor );

  factor_t factor;
  factor.factor_id = sql.last_insert_rowid();
  factor.factor_name = fac_name;
  factor.is_numeric = is_numeric;
  return factor;
}

level_t   StratOutDBase::insert_level( const std::string & level_name, const std::string & fac_name )
{
  // find factor ID
  if ( writer.factors_idmap.find( fac_name ) == writer.factors_idmap.end() )
    Helper::halt( "need to enter factor before level" );
  const factor_t & factor = writer.factors[ writer.factors_idmap[ fac_name ] ];

  sql.bind_text( stmt_insert_level , ":level_name" , level_name );
  sql.bind_int( stmt_insert_level , ":fac_id" , factor.factor_id );
  sql.step( stmt_insert_level );
  sql.reset( stmt_insert_level );

  level_t level;  
  level.level_id = sql.last_insert_rowid();
  level.level_name = level_name;
  level.factor_id = factor.factor_id;
  return level;
}

level_t   StratOutDBase::insert_level( const std::string & level_name, const int factor_id )
{
  sql.bind_text( stmt_insert_level , ":level_name" , level_name );
  sql.bind_int( stmt_insert_level , ":fac_id" , factor_id );
  sql.step( stmt_insert_level );
  sql.reset( stmt_insert_level );

  level_t level;  
  level.level_id = sql.last_insert_rowid();
  level.level_name = level_name;
  level.factor_id = factor_id;
  return level;
}

strata_t  StratOutDBase::insert_strata( const strata_t & s )
{
  
  // we should always have all existing strata in cache, so next ID is:
  
  strata_t strata;
  strata.strata_id = writer.strata.size() + 1;
  strata.levels = s.levels;

  std::map<factor_t,level_t>::const_iterator ll = s.levels.begin();
  while ( ll != s.levels.end() )
    {        
      sql.bind_int( stmt_insert_stratum , ":strata_id" , strata.strata_id );
      sql.bind_int( stmt_insert_stratum , ":level_id" , ll->second.level_id );
      sql.step( stmt_insert_stratum );
      sql.reset( stmt_insert_stratum );
      ++ll;
    }
  
  // special case for root strata (i.e. no stratifying variables)
  // use level_code of 0

  if ( s.levels.size() == 0 )
    {
      sql.bind_int( stmt_insert_stratum , ":strata_id" , strata.strata_id );
      sql.bind_int( stmt_insert_stratum , ":level_id" , 0 );
      sql.step( stmt_insert_stratum );
      sql.reset( stmt_insert_stratum );
    }
  
  return strata;
}

command_t StratOutDBase::insert_command( const std::string & cmd_name , int cmd_number, const std::string & timedate , const std::string & cmd_param )
{
  
  sql.bind_text( stmt_insert_command , ":cmd_name" , cmd_name );
  sql.bind_int( stmt_insert_command , ":cmd_number" , cmd_number );
  sql.bind_text( stmt_insert_command , ":cmd_timestamp" , timedate );
  sql.bind_text( stmt_insert_command , ":cmd_parameteres" , cmd_param );
  sql.step( stmt_insert_command );
  sql.reset( stmt_insert_command );
  
  command_t command;
  command.cmd_id = sql.last_insert_rowid();
  command.cmd_name = cmd_name;
  command.cmd_number = cmd_number;
  command.timestamp = timedate;
  command.cmd_parameters = cmd_param;
  return command;
}

bool      StratOutDBase::insert_value( const int indiv_id , const int cmd_id , const int variable_id , 
				       const int strata_id , const int timepoint_id, 
				       const value_t & x )
{

  //  std::cout << " v s value " << variable_id << "\t" << strata_id << "\t" << x.d << "\n";
  //  std::cout << writer.strata[ strata_id ].num_print() << "\n";

  sql.bind_int( stmt_insert_value , ":indiv_id" , indiv_id );
  sql.bind_int( stmt_insert_value , ":cmd_id" , cmd_id );
  sql.bind_int( stmt_insert_value , ":variable_id" , variable_id );

  if ( strata_id == -1 ) 
    sql.bind_null( stmt_insert_value , ":strata_id" );
  else
    sql.bind_int( stmt_insert_value , ":strata_id" , strata_id );
  
  if ( timepoint_id == -1 )
    sql.bind_null( stmt_insert_value , ":timepoint_id" );
  else
    sql.bind_int( stmt_insert_value , ":timepoint_id" , timepoint_id );

  if      ( x.missing ) sql.bind_null( stmt_insert_value , ":value" );
  else if ( x.numeric ) sql.bind_double( stmt_insert_value , ":value" , x.d );
  else                  sql.bind_text( stmt_insert_value , ":value" , x.s );
  
  sql.step( stmt_insert_value );
  sql.reset( stmt_insert_value );
  
  return true;
}


int StratOutDBase::num_values() 
{
  sql.step( stmt_count_values );
  int n = sql.get_int( stmt_count_values , 0 );
  sql.reset( stmt_count_values );
  return n;
}


std::map<int,int> StratOutDBase::count_strata()
{
  std::map<int,int> ret;
  while ( sql.step( stmt_count_strata ) )
    ret[ sql.get_int( stmt_count_strata , 0 ) ] = sql.get_int(stmt_count_strata , 1 ) ;
  sql.reset( stmt_count_strata );
  return ret;
}

std::map<int,std::set<int> > StratOutDBase::dump_vars_by_strata()
{
  std::map<int,std::set<int> > r;
  while ( sql.step( stmt_dump_vars_by_strata ) )
    {
      int s = sql.get_int( stmt_dump_vars_by_strata , 0 );
      int v = sql.get_int( stmt_dump_vars_by_strata , 1 );
      if ( s == 0 ) s = 1;  // translate baseline strata to '1' (default)
      r[s].insert(v);
    }
  sql.reset( stmt_dump_vars_by_strata );
  return r;
}

packets_t StratOutDBase::enumerate( int strata_id )
{

  packets_t packets;

  if ( strata_id <= 1 ) 
    {
      
      while ( sql.step( stmt_enumerate_null_strata ) )
	{
	  packet_t packet;
	  packet.indiv_id = sql.get_int( stmt_enumerate_null_strata , 0);
	  packet.cmd_id   = sql.get_int( stmt_enumerate_null_strata , 1); 
	  packet.var_id   = sql.get_int( stmt_enumerate_null_strata , 2);
	  packet.strata_id = -1;
	  packet.timepoint_id = -1;
	  packet.value = value_t( sql.get_int( stmt_enumerate_null_strata , 3) );
	  packets.push_back( packet );
	}
      sql.reset( stmt_enumerate_null_strata );

    }
  else
    {
      
      sql.bind_int( stmt_enumerate , ":strata_id" , strata_id );

      while ( sql.step( stmt_enumerate ) )
	{
	  packet_t packet;
	  packet.indiv_id = sql.get_int( stmt_enumerate , 0);
	  packet.cmd_id   = sql.get_int( stmt_enumerate , 1); 
	  packet.var_id   = sql.get_int( stmt_enumerate , 2);
	  packet.strata_id = -1;
	  packet.timepoint_id = -1;
	  packet.value = value_t( sql.get_int( stmt_enumerate , 3) );
	  packets.push_back( packet );
	}
      sql.reset( stmt_enumerate );
    }
  
  return packets;
}


std::set<int> StratOutDBase::all_matching_vars( const std::set<std::string> & vars ) 
{ 
  std::set<int> ret;  
  while ( sql.step( stmt_match_vars ) )
    {
      int var_id = sql.get_int( stmt_match_vars , 0 );
      std::string var_name = sql.get_text( stmt_match_vars , 1 );
      if ( vars.find( var_name ) != vars.end() ) ret.insert( var_id );
    }
  sql.reset( stmt_match_vars );
  return ret;
}

std::set<int> StratOutDBase::all_matching_cmds( const std::set<std::string> & cmds ) 
{ 
  std::set<int> ret;
  while ( sql.step( stmt_match_cmds ) )
    {
      int cmd_id = sql.get_int( stmt_match_cmds , 0 );
      std::string cmd_name = sql.get_text( stmt_match_cmds , 1 );
      if ( cmds.find( cmd_name ) != cmds.end() ) ret.insert( cmd_id );
    }
  sql.reset( stmt_match_cmds );
  return ret;
}



void StratOutDBase::fetch( int strata_id , int time_mode, packets_t * packets, std::set<int> * indivs_id , std::set<int> * cmds_id , std::set<int> * vars_id )
{

  if ( packets == NULL ) return;
  
  // timemode
  //  0 only get null timepoints (i.e. no E/T stratification)
  //  1 only get non-null timepoints
  
  //0 indiv_id       
  //1 cmd_id         
  //2 variable_id    
  //4 strata_id     
  //5 timepoint_id
  //6 value         

  // baseline/root values
  if ( strata_id == -1 ) 
    {
          
      sqlite3_stmt * s = stmt_lookup_value_by_null_strata ;
      
      while ( sql.step( s ) )
	{
	  
	  packet_t packet;

	  packet.indiv_id = sql.get_int( s , 0);
	  if ( indivs_id != NULL  &&  indivs_id->find( packet.indiv_id ) == indivs_id->end() ) continue;

	  packet.cmd_id   = sql.get_int( s , 1);
	  if ( cmds_id   != NULL  &&  cmds_id->find( packet.cmd_id ) == cmds_id->end() ) continue;

	  packet.var_id   = sql.get_int( s , 2);
	  if ( vars_id   != NULL  &&  vars_id->find( packet.var_id ) == vars_id->end() ) continue;

	  // no strata/timepoints here
	  packet.strata_id = -1;
	  packet.timepoint_id = -1;
	  
	  // get as a string
	  packet.value = value_t( sql.get_text( s , 5) );
	  
	  // add to outout
	  packets->push_back( packet );
	  
	}
      sql.reset( s );
 
    }
  else
    {
      
      sqlite3_stmt * s = time_mode == 1 ? stmt_lookup_value_by_strata_and_timepoint : stmt_lookup_value_by_strata ;
      
      sql.bind_int( s , ":strata_id" , strata_id );
      
      while ( sql.step( s ) )
	{
	  packet_t packet;

	  packet.indiv_id = sql.get_int( s , 0);
	  if ( indivs_id != NULL && indivs_id->find( packet.indiv_id ) == indivs_id->end() ) continue;
	  
	  packet.cmd_id = sql.get_int( s , 1);
	  if ( cmds_id != NULL && cmds_id->find( packet.cmd_id ) == cmds_id->end() ) continue;

	  packet.var_id = sql.get_int( s , 2);	  
	  if ( vars_id != NULL && vars_id->find( packet.var_id ) == vars_id->end() ) continue;
	  
	  bool has_strata = ! sql.is_null( s , 3);
	  packet.strata_id = has_strata ? sql.get_int( s , 3) : -1;      
	  
	  packet.timepoint_id = time_mode == 1 ? sql.get_int( s , 4) : -1;      
	  
	  // get as a string always
	  packet.value = value_t( sql.get_text( s , 5) );

	  // add to output
	  packets->push_back( packet );
	  
	}
      sql.reset( s );
      
    }
}


