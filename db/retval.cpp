
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

#include "retval.h"

#include "db.h"

#include "sstore/sstore.h" 

#include <iostream>


retval_factor_t::retval_factor_t( const strata_t & s , const timepoint_t & tp )
{    

  std::map<factor_t,level_t>::const_iterator aa = s.levels.begin();
  while ( aa != s.levels.end() )
    {
      const factor_t & factor = aa->first;

      // skip E/T factors here (they will be added by the timepoint_t below) 
      // also skip _COMMANDS as they are represented separately
      if ( factor.factor_name == globals::epoch_strat || 
	   factor.factor_name == globals::time_strat  || 
	   factor.factor_name[0] == '_' ) { ++aa; continue; }  
  
     factors.insert( factor.factor_name );
     
     ++aa;
    }

  // any time-points;  split T into T1 and T2
  if ( tp.epoch != -1 ) add( globals::epoch_strat );
  else if ( tp.is_interval() ) { add( globals::time_strat + "1" ); add( globals::time_strat + "2" ); }

}

retval_strata_t::retval_strata_t( strata_t & strata , timepoint_t & tp )
{

  std::map<factor_t,level_t>::const_iterator aa = strata.levels.begin();

  while ( aa != strata.levels.end() )
    {
      const factor_t & factor = aa->first;
      const level_t & level = aa->second;

      // skip E/T factors here (they will be added by the timepoint_t below) 
      // also skip _COMMANDS as they are represented separately
      if ( factor.factor_name == globals::epoch_strat || 
	   factor.factor_name == globals::time_strat || 
	   factor.factor_name[0] == '_' ) { ++aa; continue; }  
      
      // try to maintain numeric encoding of numeric factors

      if ( factor.is_numeric )	
	{
	  double lvln = 0;
	  if ( ! Helper::str2dbl( level.level_name , &lvln ) ) 
	    Helper::halt( "problem converting level to numeric:" + factor.factor_name + " " + level.level_name ); 

	  add( retval_factor_level_t( aa->first.factor_name , lvln ) );
	}
      else
	add( retval_factor_level_t( factor.factor_name , level.level_name ) );

      ++aa;
    }

  // any time-points
  if ( tp.epoch != -1 )
    add( retval_factor_level_t( globals::epoch_strat , tp.epoch ) );
  else if ( tp.is_interval() )
    {
      add( retval_factor_level_t( globals::time_strat + "1" , (double)tp.start ) );
      add( retval_factor_level_t( globals::time_strat + "2" , (double)tp.stop ) );
    }
  
}



// to a sstore_t on disk
// void retval_t::write_sstore( const std::string & f )
// {

//   //
//   // this function only used by LW for a single cmd_t and a single
//   // indiv_t mapped to a single sstore_t
//   //
  
//   if ( data.size() != 1 ) 
//     Helper::halt( "internal error, expecting a single cmd_t here" );

//   //
//   // Open/create sstore_t
//   //

//   sstore_t ssdb( f );
  
//   ssdb.drop_index();

  
//   retval_data_t::iterator cc = data.begin();
//   while ( cc != data.end() )
//     {
      
//       const retval_cmd_t & cmd = cc->first;
      
//       // factors/tables for this command
      
//       std::map<retval_factor_t,
// 	std::map<retval_var_t, 
// 	std::map<retval_strata_t, 
// 	std::map<retval_indiv_t,retval_value_t > > > >::iterator ff = cc->second.begin();
//       while ( ff != cc->second.end() )
// 	{

// 	  const retval_factor_t & fac = ff->first; 
	  
	  
	  

	  
	  
// 	  // variables
	  
// 	  std::map<retval_var_t, 
// 	    std::map<retval_strata_t, 
// 	    std::map<retval_indiv_t,retval_value_t > > >::iterator vv = ff->second.begin();

// 	  while ( vv != ff->second.end() )
// 	    {
	      
// 	      const retval_var_t & var = vv->first;
	      
// 	      // strata
	      
// 	      std::map<retval_strata_t, std::map<retval_indiv_t,retval_value_t > >::iterator ss = vv->second.begin();
// 	      while ( ss != vv->second.end() )
// 		{
		  
// 		  const retval_strata_t & strata = ss->first;

// 		  // baseline, epoch or interval-level data?
		  
// 		  // epoch-level data?
// 		  retval_factor_level_t epoch_lvl = strata.find( "E" );		  
// 		  bool epoch_level = epoch_lvl.is_int ; 
		  
		  
// 		  // interval-level or event-level data?   T1/T2, N
// 		  retval_factor_level_t t1_lvl = strata.find( "T1" );		  
// 		  retval_factor_level_t t2_lvl = strata.find( "T2" );		  
// 		  bool interval1_level = t1_lvl.is_dbl && t2_lvl.is_dbl; 

// 		  retval_factor_level_t alt1_lvl = strata.find( "START" );		  
// 		  retval_factor_level_t alt2_lvl = strata.find( "STOP" );		  
// 		  bool interval2_level = alt1_lvl.is_dbl && alt2_lvl.is_dbl;

// 		  // channel?
//                   retval_factor_level_t ch_lvl = strata.find( "CH" );
// 		  bool has_channel = ch_lvl.is_str;
// 		  std::string ch_label = has_channel ? ch_lvl.str_level : "" ;
		  
// 		  // all non-channel, non-epoch, non-interval factors
		  
// 		  std::stringstream sss;
// 		  bool first = true;
// 		  std::set<retval_factor_level_t>::const_iterator ff = strata.factors.begin();
// 		  while ( ff != strata.factors.end() )
// 		    {
// 		      if ( ff->factor == "CH" ) { ++ff; continue; }
// 		      if ( ff->factor == "E" ) { ++ff; continue; }
// 		      if ( ff->factor == "T1" ) { ++ff; continue; }
// 		      if ( ff->factor == "T2" ) { ++ff; continue; }
// 		      if ( ff->factor == "N" ) { ++ff; continue; }
// 		      if ( ! first ) sss << ";";
// 		      sss << ff->print();
// 		      first = false;
// 		      ++ff;
// 		    }
		  
// 		  std::string lvl_label = sss.str();
// 		  bool has_lvl = lvl_label != "";

// 		  // individual
	
// 		  // expecting only a single individual
// 		  if ( ss->second.size() > 1 ) Helper::halt( "only expecting a single indiv_t here" );
		  
// 		  std::map<retval_indiv_t,retval_value_t>::iterator ii = ss->second.begin();
// 		  while ( ii != ss->second.end() )
// 		    {
		      
// 		      const retval_value_t & value = ii->second;
		      
// 		      // output
		      
// 		      if ( epoch_level )
// 			ssdb.insert_epoch( epoch_lvl.int_level , var.name  , value.print() , 
// 					   has_channel ? &ch_label : NULL , 
// 					   has_lvl ? &lvl_label : NULL ) ;		      
// 		      else if ( interval1_level ) // T1/T2
// 			ssdb.insert_interval( t1_lvl.dbl_level , t2_lvl.dbl_level , 
// 					      var.name  , value.print() , 
// 					      has_channel ? &ch_label : NULL , 
// 					      has_lvl ? &lvl_label : NULL ) ;		      
// 		      else if ( interval2_level ) // similar, but START/STOP 
// 			ssdb.insert_interval( alt1_lvl.dbl_level , alt2_lvl.dbl_level , 
// 					      var.name  , value.print() , 
// 					      has_channel ? &ch_label : NULL , 
// 					      has_lvl ? &lvl_label : NULL ) ;		      
// 		      else

// 			ssdb.insert_base( var.name  , value.print() , 
// 					  has_channel ? &ch_label : NULL , 
// 					  has_lvl ? &lvl_label : NULL ) ;		      
		      



// // 			ss.insert_interval( start , stop , annot , inst_lvl.str_level , NULL , NULL );
// // 		      else
// // 			ss.insert_base( annot , inst_lvl.str_level , NULL , NULL );

// // 		      std::cout << ii->first.name << "\t"
// // 				<< cmd.name << "\t"
// // 				<< fac.print() << "\t"
// // 				<< var.name << "\t"
// // 				<< strata.print() << "\t"
// // 				<< value.print() << "\n";
		  
// 		      ++ii; // next individual
// 		    }
// 		  ++ss; // next strata
// 		}
// 	      ++vv; // next variable
// 	    }
// 	 ++ff; // next factor/table
//          }
//       ++cc; // next command
//     }      
//       // all done

  
  
//   ssdb.index();
  
//   ssdb.dettach();

// }


void retval_t::dump()
{

  // std::map<retval_cmd_t, std::map<retval_factor_t,std::map<retval_var_t, std::map<retval_strata_t, std::map<retval_indiv_t,retval_value_t > > > > > retval_data_t;

  retval_data_t::iterator cc = data.begin();
  while ( cc != data.end() )
    {
      
      const retval_cmd_t & cmd = cc->first;
      
      //std::cout << "considering command " << cmd.name << "\n";

      // factors/tables for this command
      
      std::map<retval_factor_t,std::map<retval_var_t, std::map<retval_strata_t, std::map<retval_indiv_t,retval_value_t > > > >::iterator ff = cc->second.begin();
      while ( ff != cc->second.end() )
	{

	  const retval_factor_t & fac = ff->first; 
	  
	  // variables
	  
	  std::map<retval_var_t, std::map<retval_strata_t, std::map<retval_indiv_t,retval_value_t > > >::iterator vv = ff->second.begin();
	  while ( vv != ff->second.end() )
	    {

	      const retval_var_t & var = vv->first;

	      // strata
	      
	      std::map<retval_strata_t, std::map<retval_indiv_t,retval_value_t > >::iterator ss = vv->second.begin();
	      while ( ss != vv->second.end() )
		{
		  
		  const retval_strata_t & strata = ss->first;
		  
		  // individual
		  
		  std::map<retval_indiv_t,retval_value_t>::iterator ii = ss->second.begin();
		  while ( ii != ss->second.end() )
		    {
		      
		      const retval_value_t & value = ii->second;
		      
		      // output
		      std::cout << ii->first.name << "\t"
				<< cmd.name << "\t"
				<< fac.print() << "\t"
				<< var.name << "\t"
				<< strata.print() << "\t"
				<< value.print() << "\n";
		  
		      ++ii; // next individual
		    }
		  ++ss; // next strata
		}
	      ++vv; // next variable
	    }
	 ++ff; // next factor/table
         }
      ++cc; // next command
    }      
      // all done
}

