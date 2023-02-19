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


tfac_t strata_t::tfac() const
{

  tfac_t tfac("");

  std::map<factor_t,level_t>::const_iterator aa = levels.begin();
  while ( aa != levels.end() )
    {

      // skip commands
      if ( aa->first.factor_name[0] == '_' ) 
	{ ++aa; continue; } 
      
      // skip tags
      if ( globals::cmddefs().is_tag( aa->first.factor_name ) )
	{ ++aa; continue; }
      
      // otherwise, add (to ID which zfile_t to write to)
      tfac.fac.insert( aa->first.factor_name  );

      ++aa;
    }
  return tfac;
}

std::map<std::string,std::string> writer_t::faclvl() const
{
  // no commands, but includes TAGs and factors, as key-value map
  std::map<std::string,std::string> r;

  std::map<factor_t,level_t>::const_iterator aa = curr_strata.levels.begin();
  while ( aa != curr_strata.levels.end() )
    {
      // skip commands
      if ( aa->first.factor_name[0] == '_' ) { ++aa; continue; } 

      // epoch/time-point? levels stored separately
      if ( aa->first.factor_name == globals::epoch_strat || 
	   aa->first.factor_name == globals::time_strat )
	{
	  if ( curr_timepoint.none() ) r[ aa->first.factor_name ] = "." ;
	  else r[ aa->first.factor_name ] = curr_timepoint.print();
	}
      else      
	r[ aa->first.factor_name ] = aa->second.level_name ; 
      ++aa;
    }
  return r;
}

std::string strata_t::print_zfile_tag() const
{

  // factors only, (no commands, no TAGs) in a underscore-delim list

  if ( levels.size() == 0 ) return ""; // baseline

  std::stringstream ss;
  std::map<factor_t,level_t>::const_iterator aa = levels.begin();
  bool printed = false;
  while ( aa != levels.end() )
    {
      // skip commands
      if ( aa->first.factor_name[0] == '_' ) { ++aa; continue; } 
      
      if ( printed ) ss << "_";
      ss << aa->first.factor_name ;
      printed = true;
      ++aa;
    }
  std::string rstr = ss.str();
  return ss.str();
}
 
std::string strata_t::print_nocmd() const
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
      
      // skip commands
      if ( aa->first.factor_name[0] == '_' ) { ++aa; continue; } 

      if ( printed ) ss << ";";
      ss << aa->first.factor_name << "/" << aa->second.level_name ; 
      printed = true;
      ++aa;
    }
  std::string rstr = ss.str();
  if ( rstr == "" ) return ".";
  return ss.str();
}


bool StratOutDBase::attach( const std::string & n , bool readonly , writer_t * caller )
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

  read_all( caller );

       
  //
  // specify types for common stratifiers
  //

  caller->set_types();
  
  
  //
  // all done
  //

  return true;

}



void StratOutDBase::read_all( writer_t * w )
{

  //
  // dump everything into writer_t *except* datapoints
  //


  //
  // Individuals
  //

  while ( sql.step( stmt_dump_individuals ) )
    {
      indiv_t indiv;
      indiv.indiv_id = sql.get_int( stmt_dump_individuals , 0 );
      indiv.indiv_name = sql.get_text( stmt_dump_individuals , 1 );
      w->individuals[ indiv.indiv_id ] = indiv;
      w->individuals_idmap[ indiv.indiv_name ] = indiv.indiv_id;
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
      w->commands[ cmd.cmd_id ] = cmd;
      w->commands_idmap[ cmd.cmd_name ] = cmd.cmd_id;
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
      w->factors[ factor.factor_id ] = factor;
      w->factors_idmap[ factor.factor_name ] = factor.factor_id;
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
      w->levels[ level.level_id ] = level;
      if ( w->factors.find( level.factor_id ) == w->factors.end() )
	Helper::halt( "internal error, undefined factor" );
      std::string level_key = level.level_name 
	+ "." + w->factors[ level.factor_id ].factor_name ;
      w->levels_idmap[ level_key ] = level.level_id;
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
	  if ( w->strata.find( strata_id ) == w->strata.end() )
	    {
	      strata_t s;
	      s.strata_id = strata_id;
	      w->strata[ strata_id ] = s;	
	    }
	  
	}
      
      else
	{
	  level_t level; 
	  factor_t factor; 
	  
	  level = w->levels[ level_id ]; 
	  factor = w->factors[ level.factor_id ];	    
	  
	  // existing strata?
	  if ( w->strata.find( strata_id ) != w->strata.end() )
	    {
	      w->strata[ strata_id ].insert( level , factor );
	    }
	  else
	    {
	      strata_t s;
	      s.strata_id = strata_id;
	      
	      s.insert( level , factor );
	      
	      w->strata[ strata_id ] = s;	
	    }
	}

    }
  sql.reset( stmt_dump_strata );


  //
  // Update strata idmap
  //

  w->strata_idmap.clear();
  std::map<int,strata_t>::const_iterator ss = w->strata.begin();
  while ( ss != w->strata.end() )
    {
      w->strata_idmap[ ss->second ] = ss->first;
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
      timepoint.start = has_interval ? sql.get_uint64( stmt_dump_timepoints , 2 ) : 0 ; 
      timepoint.stop  = has_interval ? sql.get_uint64( stmt_dump_timepoints , 3 ) : 0 ; 
      
      std::string tp_key = ( has_epoch ? Helper::int2str( timepoint.epoch ) : "" ) + ":" 
	+ ( has_interval ? Helper::int2str( timepoint.start ) + "-" + Helper::int2str( timepoint.stop ) : "" ); 
      
      w->timepoints[ timepoint.timepoint_id ] = timepoint;
      w->timepoints_idmap[ tp_key ] = timepoint.timepoint_id;

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
      if ( w->commands_idmap.find( command_name ) != w->commands_idmap.end() ) 
	var.cmd_id = w->commands_idmap[ command_name ];

      w->variables[ var.var_id ] = var;
      w->variables_idmap[ command_name + ":" + var.var_name ] = var.var_id;
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

  // dumpers

  stmt_dump_factors = sql.prepare( "SELECT * FROM factors;" );
  stmt_dump_levels = sql.prepare( "SELECT * FROM levels;" );
  stmt_dump_strata = sql.prepare( "SELECT * FROM strata;" );
  stmt_dump_variables = sql.prepare( "SELECT * FROM variables;" );
  stmt_dump_individuals = sql.prepare( "SELECT * FROM individuals;" );
  stmt_dump_timepoints = sql.prepare( "SELECT * FROM timepoints;" );
  stmt_dump_commands = sql.prepare( "SELECT * FROM commands;" );

  // for datapoints, these are only pulled for a given individual (used only when making a retval_t)
  stmt_dump_int_datapoints = sql.prepare( "SELECT * FROM datapoints where indiv_id == :indiv_id AND typeof(value) == \"integer\" ;" );
  stmt_dump_dbl_datapoints = sql.prepare( "SELECT * FROM datapoints where indiv_id == :indiv_id AND typeof(value) == \"real\" ;" );
  stmt_dump_txt_datapoints = sql.prepare( "SELECT * FROM datapoints where indiv_id == :indiv_id AND typeof(value) == \"text\" ;" );
  
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

  sql.finalise( stmt_dump_int_datapoints);	      
  sql.finalise( stmt_dump_dbl_datapoints);	      
  sql.finalise( stmt_dump_txt_datapoints);	      

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
  sql.bind_uint64( stmt_insert_timepoint , ":start" , interval.start );
  sql.bind_uint64( stmt_insert_timepoint , ":stop" , interval.stop );
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
  sql.bind_text( stmt_insert_command , ":cmd_parameters" , cmd_param );
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

  if      ( x.missing ) sql.bind_null( stmt_insert_value ,   ":value" );
  else if ( x.numeric ) sql.bind_double( stmt_insert_value , ":value" , x.d );
  else if ( x.integer ) sql.bind_int( stmt_insert_value ,    ":value" , x.i );
  else                  sql.bind_text( stmt_insert_value ,   ":value" , x.s );
  
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



packets_t StratOutDBase::dump_all() 
{
  
  packets_t packets;

  while ( sql.step( stmt_dump_int_datapoints ) )
    {      
      packet_t packet;
      packet.indiv_id = sql.get_int( stmt_dump_int_datapoints , 0);
      packet.cmd_id   = sql.get_int( stmt_dump_int_datapoints , 1);
      packet.var_id   = sql.get_int( stmt_dump_int_datapoints , 2);
      bool has_strata = ! sql.is_null( stmt_dump_int_datapoints , 3);
      packet.strata_id = has_strata ? sql.get_int( stmt_dump_int_datapoints , 3) : -1;            
      bool has_tp = ! sql.is_null( stmt_dump_int_datapoints , 4); 
      packet.timepoint_id = has_tp ? sql.get_int( stmt_dump_int_datapoints , 4) : -1;
      // Integers
      packet.value = value_t( sql.get_int( stmt_dump_int_datapoints , 5) );
      packets.push_back( packet );
      
    }
  sql.reset( stmt_dump_int_datapoints );


  while ( sql.step( stmt_dump_dbl_datapoints ) )
    {      
      packet_t packet;
      packet.indiv_id = sql.get_int( stmt_dump_dbl_datapoints , 0);
      packet.cmd_id   = sql.get_int( stmt_dump_dbl_datapoints , 1);
      packet.var_id   = sql.get_int( stmt_dump_dbl_datapoints , 2);
      bool has_strata = ! sql.is_null( stmt_dump_dbl_datapoints , 3);
      packet.strata_id = has_strata ? sql.get_int( stmt_dump_dbl_datapoints , 3) : -1;      
      bool has_tp = ! sql.is_null( stmt_dump_dbl_datapoints , 4); 
      packet.timepoint_id = has_tp ? sql.get_int( stmt_dump_dbl_datapoints , 4) : -1;
      // doubles
      packet.value = value_t( sql.get_double( stmt_dump_dbl_datapoints , 5) );      
      packets.push_back( packet );      
    }
  sql.reset( stmt_dump_dbl_datapoints );


  while ( sql.step( stmt_dump_txt_datapoints ) )
    {     
      packet_t packet;
      packet.indiv_id = sql.get_int( stmt_dump_txt_datapoints , 0);
      packet.cmd_id   = sql.get_int( stmt_dump_txt_datapoints , 1);
      packet.var_id   = sql.get_int( stmt_dump_txt_datapoints , 2);
      bool has_strata = ! sql.is_null( stmt_dump_txt_datapoints , 3);
      packet.strata_id = has_strata ? sql.get_int( stmt_dump_txt_datapoints , 3) : -1;      
      bool has_tp = ! sql.is_null( stmt_dump_txt_datapoints , 4); 
      packet.timepoint_id = has_tp ? sql.get_int( stmt_dump_txt_datapoints , 4) : -1;
      // Strings
      packet.value = value_t( sql.get_text( stmt_dump_txt_datapoints , 5) );
      packets.push_back( packet );      
    }
  sql.reset( stmt_dump_txt_datapoints );

  return packets;

}


packets_t StratOutDBase::dump_indiv( const int indiv_id ) 
{
  
  packets_t packets;
  
  sql.bind_int( stmt_dump_int_datapoints , ":indiv_id" , indiv_id );
  while ( sql.step( stmt_dump_int_datapoints ) )
    {      
      packet_t packet;
      packet.indiv_id = sql.get_int( stmt_dump_int_datapoints , 0);
      packet.cmd_id   = sql.get_int( stmt_dump_int_datapoints , 1);
      packet.var_id   = sql.get_int( stmt_dump_int_datapoints , 2);
      bool has_strata = ! sql.is_null( stmt_dump_int_datapoints , 3);
      packet.strata_id = has_strata ? sql.get_int( stmt_dump_int_datapoints , 3) : -1;            
      bool has_tp = ! sql.is_null( stmt_dump_int_datapoints , 4); 
      packet.timepoint_id = has_tp ? sql.get_int( stmt_dump_int_datapoints , 4) : -1;
      // Integers
      packet.value = value_t( sql.get_int( stmt_dump_int_datapoints , 5) );
      packets.push_back( packet );
      
    }
  sql.reset( stmt_dump_int_datapoints );

  sql.bind_int( stmt_dump_dbl_datapoints , ":indiv_id" , indiv_id );
  while ( sql.step( stmt_dump_dbl_datapoints ) )
    {      
      packet_t packet;
      packet.indiv_id = sql.get_int( stmt_dump_dbl_datapoints , 0);
      packet.cmd_id   = sql.get_int( stmt_dump_dbl_datapoints , 1);
      packet.var_id   = sql.get_int( stmt_dump_dbl_datapoints , 2);
      bool has_strata = ! sql.is_null( stmt_dump_dbl_datapoints , 3);
      packet.strata_id = has_strata ? sql.get_int( stmt_dump_dbl_datapoints , 3) : -1;      
      bool has_tp = ! sql.is_null( stmt_dump_dbl_datapoints , 4); 
      packet.timepoint_id = has_tp ? sql.get_int( stmt_dump_dbl_datapoints , 4) : -1;
      // doubles
      packet.value = value_t( sql.get_double( stmt_dump_dbl_datapoints , 5) );      
      packets.push_back( packet );      
    }
  sql.reset( stmt_dump_dbl_datapoints );

  sql.bind_int( stmt_dump_txt_datapoints , ":indiv_id" , indiv_id );
  while ( sql.step( stmt_dump_txt_datapoints ) )
    {     
      packet_t packet;
      packet.indiv_id = sql.get_int( stmt_dump_txt_datapoints , 0);
      packet.cmd_id   = sql.get_int( stmt_dump_txt_datapoints , 1);
      packet.var_id   = sql.get_int( stmt_dump_txt_datapoints , 2);
      bool has_strata = ! sql.is_null( stmt_dump_txt_datapoints , 3);
      packet.strata_id = has_strata ? sql.get_int( stmt_dump_txt_datapoints , 3) : -1;      
      bool has_tp = ! sql.is_null( stmt_dump_txt_datapoints , 4); 
      packet.timepoint_id = has_tp ? sql.get_int( stmt_dump_txt_datapoints , 4) : -1;
      // Strings
      packet.value = value_t( sql.get_text( stmt_dump_txt_datapoints , 5) );
      packets.push_back( packet );      
    }
  sql.reset( stmt_dump_txt_datapoints );

  return packets;

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


bool writer_t::close() 
{ 

  // if in plaintext mode, close out if needed and clean up
  if ( plaintext ) 
    { 	
      if ( zfiles != NULL ) 
	{
	  
	  update_plaintext_curr_strata();

	  zfiles->close();

	  delete zfiles;

	  zfiles = NULL;

	}
      
    }



  // otherwise, handle any DB-related stuff
  if ( ! attached() ) return false;
  clear(); 
  db.dettach();
  
  return true; 
} 


retval_t writer_t::dump_to_retval( const std::string & dbname , const std::set<std::string> * persons , std::vector<std::string> * ids )
{

  //
  // Pulls out all information for one or more individuals/EDFs and
  // creates a single retval_t
  //

  retval_t retval;
  
  //
  // attach named database as read-only
  //
  
  const bool IS_READONLY = true;
  
  writer_t w;
  
  w.attach( dbname , IS_READONLY );


  //
  // Read all individuals, or a subset?
  //

  bool read_all_individuals = persons == NULL || persons->size() == 0 ; 
  
  
  //
  // Loop over each individual in the DB
  //

  std::map<std::string,int>::const_iterator ii = w.individuals_idmap.begin();
  while ( ii != w.individuals_idmap.end() )
    {

      const std::string & indiv_name = ii->first;

      if ( ! read_all_individuals )
	{
	  // not on the list, so skip
	  if ( persons->find( indiv_name ) == persons->end() ) 
	    {
	      ++ii;
	      continue;
	    }
	}
	

      //
      // track who is being read
      //

      if ( ids != NULL ) ids->push_back( indiv_name );

      //
      // get numeric indiv ID code
      //

      int indiv_id = w.individuals_idmap[ indiv_name ];
  
      //
      // separately dump all int, double and text values, so that appropriate retval_t types can be set
      //
      
      packets_t packets = w.db.dump_indiv( indiv_id );

      //
      // Convert packets_t to retval_t
      //  
      
      packets_t::const_iterator pp = packets.begin();
      while ( pp != packets.end() )
	{
	  
	  if ( pp->value.numeric )
	    retval.add( indiv_name , 
			retval_cmd_t(  w.commands[ pp->cmd_id ].cmd_name ) , 
			retval_factor_t( w.strata[ pp->strata_id] , w.timepoints[ pp->timepoint_id ] ) , 
			retval_var_t( w.variables[ pp->var_id ].var_name ) , 
			retval_strata_t( w.strata[ pp->strata_id] , w.timepoints[ pp->timepoint_id ] ) ,
			pp->value.d );
	  
	  else if ( pp->value.integer )
	    retval.add( indiv_name , 
			retval_cmd_t(  w.commands[ pp->cmd_id ].cmd_name ) , 
			retval_factor_t( w.strata[ pp->strata_id] , w.timepoints[ pp->timepoint_id ] ) , 
			retval_var_t( w.variables[ pp->var_id ].var_name ) , 
			retval_strata_t( w.strata[ pp->strata_id] , w.timepoints[ pp->timepoint_id ] ) ,
			pp->value.i );
	  else
	    retval.add( indiv_name , 
			retval_cmd_t(  w.commands[ pp->cmd_id ].cmd_name ) , 
			retval_factor_t( w.strata[ pp->strata_id] , w.timepoints[ pp->timepoint_id ] ) , 
			retval_var_t( w.variables[ pp->var_id ].var_name ) , 
			retval_strata_t( w.strata[ pp->strata_id] , w.timepoints[ pp->timepoint_id ] ) ,
			pp->value.s );

         
	  ++pp;
	}

      
      //
      // next individual
      //

      ++ii;

    }
  
  //
  // all done, return info
  //

  return retval;

}

bool writer_t::to_plaintext( const std::string & var_name , const value_t & x ) 
{

  // if trying to write to an ill-formed table, complain
  
  if ( curr_zfile == NULL ) 
    {
            
      if ( zfiles != NULL )
	{
	  zfiles->close();
	  delete zfiles;
	  zfiles = NULL;
	}

      Helper::halt( "internal error: null curr_zfile in writer_t: " + var_name + 
		    "\n -- output tables for this command have not yet been hooked up for '-t' mode output" 
		    "\n -- please re-run without -t (i.e. -o/-a or raw output to the console) " );
    }
  
  // write variable/value to buffer
  
  curr_zfile->set_value( var_name , x.str() );    

  return true;
}


void writer_t::update_plaintext_curr_strata()
{

  if ( zfiles == NULL ) return ;

  // get zfile-set for this individual, creating if it does not exist

  // note: need to get param into here... (i.e. when writer_t::cmd() function)

  // figure out which table (command/strata)  
  curr_zfile = zfiles->file( curr_command.cmd_name , NULL , curr_strata.print_zfile_tag() ) ;  

  // might not be a valid table (i.e. this could be the case if setting 
  // levels, e.g. A+B,  then when only level(A) is set, it will not be valid
  // this is fine, so we won't give an error yet;  but if somebody tries writing 
  // to a NULL,  writer_t::to_plaintext(), give an error then.

  if ( curr_zfile == NULL ) return;
  
  // set (all) levels 
  curr_zfile->set_stratum( faclvl() );

}


bool writer_t::level( const std::string & level_name , const std::string & factor_name )
{
  
  //std::cout << "lvl: " << factor_name << " " << level_name << "\n";
  
  // add factor (as string by default) if it doesn't already exist
  if ( factors_idmap.find( factor_name ) == factors_idmap.end() ) 
    {
      //std::cout << "setting str factor\n";
      string_factor( factor_name );
    }

  // std::cout << " factors_idmap size = " << factors_idmap.size() << "\n";
  
  // std::cout << "fidmap = " << factors_idmap[ factor_name ] << "\n";
  
  factor_t factor = factors[ factors_idmap[ factor_name ] ];

  //  std::cout << "fname " << factor.factor_name << "\n";
  
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
  //  std::cout << "lvl-set: " << level.level_name << " " << factor.factor_name << "\n";
  // swap/add to current strata
  curr_strata.insert( level , factor );
  
  // if needed, update which table to point to
  if ( plaintext ) update_plaintext_curr_strata();
  
  return true;
}

bool writer_t::numeric_factor( const std::string & fac_name )
  {
    if ( factors_idmap.find( fac_name ) == factors_idmap.end() )
      {
	factor_t factor = db.insert_factor( fac_name , 1 );  // 1 -> numeric factor
	//	std::cout << "ASSIGN: num " << fac_name << " " << factor.factor_id << "\n";
	factors_idmap[ fac_name ] = factor.factor_id;
	factors[ factor.factor_id ] = factor;
      }
    return true;
  }

bool writer_t::string_factor( const std::string & fac_name )
  {
    if ( factors_idmap.find( fac_name ) == factors_idmap.end() )
      {
	factor_t factor = db.insert_factor( fac_name , 0 ); // 0 -> string factor  
	//std::cout << "ASSIGN: str " << fac_name << " " << factor.factor_id << "\n";
	factors_idmap[ fac_name ] = factor.factor_id;
	factors[ factor.factor_id ] = factor;
      }
    return true;
  }


void writer_t::set_types()
{

  //  std::cout << "SETTING_TYPES()\n";
  
  numeric_factor( globals::epoch_strat );
  numeric_factor( globals::sample_strat );
  numeric_factor( globals::freq_strat );
  numeric_factor( globals::sec_strat );
  numeric_factor( globals::cycle_strat );
  string_factor( globals::band_strat );
  string_factor( globals::annot_strat );
  string_factor( globals::annot_instance_strat );
  string_factor( globals::annot_meta_strat );
  string_factor( globals::signal_strat );
  string_factor( globals::stage_strat );
  numeric_factor( globals::count_strat );
  numeric_factor( globals::time_strat );
  numeric_factor( globals::value_strat );

  numeric_factor( "EID" );
  numeric_factor( "IC" );
  numeric_factor( "TAP" );
  numeric_factor( "TH" );
  numeric_factor( "SPINDLE" );
  numeric_factor( "MSEC" );
  numeric_factor( "PHASE" );
  numeric_factor( "PSC" );
  numeric_factor( "SEG" );
}
