
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
      std::cout << "fac = " << factor.factor_name << "\n";
      if ( factor.is_numeric )	
	{
	  double lvln = 0;
	  if ( ! Helper::str2dbl( level.level_name , &lvln ) ) 
	    Helper::halt( "problem converting level to numeric:" + factor.factor_name + " " + level.level_name ); 

	  std::cout << "num " << lvln << "\n";

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



void retval_t::dump()
{
  // std::map<retval_cmd_t, std::map<retval_factor_t,std::map<retval_var_t, std::map<retval_strata_t, retval_value_t > > > > retval_data_t;

  retval_data_t::iterator cc = data.begin();
  while ( cc != data.end() )
    {
      
      const retval_cmd_t & cmd = cc->first;
      
      //std::cout << "considering command " << cmd.name << "\n";

      // factors/tables for this command
      
      std::map<retval_factor_t,std::map<retval_var_t, std::map<retval_strata_t, retval_value_t > > >::iterator ff = cc->second.begin();
      while ( ff != cc->second.end() )
	{

	  const retval_factor_t & fac = ff->first; 
	  
	  // variables
	  
	  std::map<retval_var_t, std::map<retval_strata_t, retval_value_t > >::iterator vv = ff->second.begin();
	  while ( vv != ff->second.end() )
	    {

	      const retval_var_t & var = vv->first;

	      // strata
	      
	      std::map<retval_strata_t, retval_value_t >::iterator ss = vv->second.begin();
	      while ( ss != vv->second.end() )
		{
		  
		  const retval_strata_t & strata = ss->first;
		  
		  const retval_value_t & value = ss->second;
		  
		  // output
		  std::cout << cmd.name << "\t"
			    << fac.print() << "\t"
			    << var.name << "\t"
			    << strata.print() << "\t"
			    << value.print() << "\n";
		  
		  
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

