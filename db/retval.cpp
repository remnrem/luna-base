
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

void retval_t::clear()
{
  data.clear();
  var_has_strings.clear();
  var_has_doubles.clear();
}


// add a double
void retval_t::add( const retval_indiv_t & id, 
		    const retval_cmd_t & cmd ,
		    const retval_factor_t & fac ,
		    const retval_var_t & var ,
		    const retval_strata_t & stratum ,
		    const double x )
{

  // std::cout << "add: "
  // 	    << id.name << " "
  // 	    << cmd.name << " "
  // 	    << fac.print() << " "
  // 	    << var.name << " "
  // 	    << stratum.print() << " "
  // 	    << x << "\n";
  
  var_has_doubles.insert( var.name );
  data[cmd][fac][var][stratum][id] = retval_value_t( x );
}

void retval_t::add( const retval_indiv_t & id , 
		    const retval_cmd_t & cmd ,
		    const retval_factor_t & fac ,
		    const retval_var_t & var ,
		    const retval_strata_t & stratum ,
		    int  x )
{
  data[cmd][fac][var][stratum][id] = retval_value_t( (int64_t)x );
}

void retval_t::add( const retval_indiv_t & id, 
		    const retval_cmd_t & cmd ,
		    const retval_factor_t & fac ,
		    const retval_var_t & var ,
		    const retval_strata_t & stratum ,
		    int64_t  x )
{
  data[cmd][fac][var][stratum][id] = retval_value_t( x );
}

void retval_t::add( const retval_indiv_t & id, 
		    const retval_cmd_t & cmd ,
		    const retval_factor_t & fac ,
		    const retval_var_t & var ,
		    const retval_strata_t & stratum ,
		    uint64_t  x )
{
  data[cmd][fac][var][stratum][id] = retval_value_t( (int64_t)x );
}

void retval_t::add( const retval_indiv_t & id, 
		    const retval_cmd_t & cmd ,
		    const retval_factor_t & fac ,
		    const retval_var_t & var ,
		    const retval_strata_t & stratum ,
		    const std::string & x )
{  
  var_has_strings.insert( var.name );
  data[cmd][fac][var][stratum][id] = retval_value_t( x );
}


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




//
// return in tabular form
//

std::map<std::string,std::map<std::string,rtable_t> > retval_t::make_tables() const 
{

  // list [ cmd ]
  //    list [ strata ]  e.g..  F_CH
  //            data.frame   :  cols = facs + vars  ; rows = lvls + values 
  
  // e.g.    F    CH   DENS  AMP
  //         11   C3   1.23  0.23
  //         16   C4   1.23  0.23
  
  // # of commands
  const int nc = data.size();

  // output
  std::map<std::string,std::map<std::string,rtable_t> > tables;
  
  // iterate over each command
  
  retval_data_t::const_iterator cc = data.begin();
  while ( cc != data.end() )
    {      
      
      const retval_cmd_t & cmd = cc->first;
      
      // number of virtual tables/factors for this command 
      
      const int nt = cc->second.size();
      
      // factors/tables
      
      std::map<retval_factor_t,
	       std::map<retval_var_t,
	std::map<retval_strata_t,
		 std::map<retval_indiv_t,
	retval_value_t > > > >::const_iterator tt = cc->second.begin();

      //
      // iterate over each table
      //

      while ( tt != cc->second.end() )
	{
	  
	  const retval_factor_t & table = tt->first; 
	  
	  //
	  // for this particular cmd/fac combination, make a data-frame
	  //
	  
	  //    cols = ID + factors + variables
	  //    rows = levels
	  //    items = values
	  
	  // in this table ( tt ):
	  // how many factors (i.e. F CH == 2 )
	  // how many variables

	  // but we need to split these by type
	  
	  const int nf = table.factors.size();
	  
	  const int nv = tt->second.size();
	  
	  const int ncols = 1 + nf + nv;  // 1 for ID 

	  //
	  // quickly scan for the a) number of rows
	  // and b) whether the factors are string, double or int
	  //

	  // str > dbl > int
	  std::set<std::string> int_factor, str_factor, dbl_factor;
	  
	  //
	  // to track rows, indiv/strata pairing
	  //
	  	  
	  std::set<retval_indiv_strata_t> rows;
	  
	  std::map<retval_var_t,
	    std::map<retval_strata_t,
	    std::map<retval_indiv_t,
	    retval_value_t > > >::const_iterator vv = tt->second.begin();

	  while ( vv != tt->second.end() )
	    {

	      std::map<retval_strata_t, 
		std::map<retval_indiv_t,retval_value_t > >::const_iterator ss = vv->second.begin();
	   
	      while ( ss != vv->second.end() )
		{
		  
		  const retval_strata_t & s = ss->first;
		  
		  // get rows (+ indivs)
		  std::map<retval_indiv_t,retval_value_t>::const_iterator ii = ss->second.begin();
		  while ( ii != ss->second.end() )
		    {
		      rows.insert( retval_indiv_strata_t( ii->first , s )  );
		      ++ii;
		    }
		      
		  // get factor types 
		  std::set<retval_factor_level_t>::const_iterator ll = s.factors.begin();
		  while ( ll != s.factors.end() )
		    {
		      if      ( ll->is_str ) str_factor.insert( ll->factor );
		      else if ( ll->is_dbl ) dbl_factor.insert( ll->factor );
		      else if ( ll->is_int ) int_factor.insert( ll->factor );
		      ++ll;
		    }
		      
		  ++ss;
		}
	      ++vv; 
	    }
	  
	  
	  //
	  // Now, we should know the number of rows, and whether a
	  // given factor is string, double or int
	  //
	  
	  const int nrows = rows.size();

	  //
	  // we now need to build a matrix of 'nrows' rows and 'ncols' colums (fac + vars)
	  //
	  
	  rtable_t df; 
	  
	  // // set class attribute for df
	  // PROTECT(cls = Rf_allocVector(STRSXP, 1)); // class attribute
	  // protect();
	  
	  // SET_STRING_ELT(cls, 0, Rf_mkChar( "data.frame" ));
	  // Rf_classgets(df, cls);
	  
	  // // col-names
	  // PROTECT(nam = Rf_allocVector(STRSXP, ncols)); // names attribute (column names)
	  // protect();
	  


	  //
	  // Add ID as column 1 
	  //

	  std::vector<std::string> id_col( nrows , "." );
	  
	  // populate w/ IDs
	  // consider all indiv/factor/level rows
	  int r_cnt = 0;
	  std::set<retval_indiv_strata_t>::const_iterator rr =  rows.begin();
	  while ( rr != rows.end() )
	    {	      	      
	      id_col[ r_cnt ] = rr->indiv.name ;
	      ++r_cnt;
	      ++rr;
	    }

	  df.add( "ID" , id_col );

	  //
	  // Add factors
	  //
	  
	  std::set<std::string>::const_iterator ff = table.factors.begin();
	  while ( ff != table.factors.end() )
	    {

	      bool is_str_factor = false , is_dbl_factor = false , is_int_factor = false ;
	      
	      if ( str_factor.find( *ff ) != str_factor.end() )
		is_str_factor = true; 
	      else if ( dbl_factor.find( *ff ) != dbl_factor.end() )
		is_dbl_factor = true;
	      else
		is_int_factor = true;
		
	      std::vector<std::string> strcol;
	      std::vector<int> intcol;
	      std::vector<double> dblcol;
	      std::vector<bool> missing( nrows , false );
		
	      if ( is_str_factor ) 
		strcol.resize( nrows );
	      else if ( is_dbl_factor )
		dblcol.resize( nrows );
	      else
		intcol.resize( nrows );
	      
	      // consider all indiv/factor/level rows
	      int r_cnt = 0;
	      std::set<retval_indiv_strata_t>::const_iterator rr =  rows.begin();
	      while ( rr != rows.end() )
		{
		  
		  //retval_indiv_t indiv = rr->indiv;
		  const retval_strata_t & strata = rr->strata;
		  const retval_factor_level_t & lvl = strata.find( *ff );
		  
		  // get value from 'fac', but bear in mind, it may
		  // be of different type (would be v. strange, but
		  // handle here just in case, w/ a cast)
		  
		  if ( is_str_factor )
		    {
		      if ( lvl.is_str )
			strcol[ r_cnt ] = lvl.str_level;
		      else if ( lvl.is_int )
			strcol[ r_cnt ]	= Helper::int2str(lvl.int_level);
		      else if ( lvl.is_dbl )
			strcol[ r_cnt ] = Helper::dbl2str(lvl.dbl_level);
		      else
			missing[ r_cnt ] = true;
		    }
		  else if ( is_dbl_factor )
		    {
		      if ( lvl.is_dbl )
			dblcol[ r_cnt ] = lvl.dbl_level;		      
		      else if ( lvl.is_int )
			dblcol[ r_cnt ] = lvl.int_level;
		      else
			missing[ r_cnt ] = true;
		    }
		  else if ( is_int_factor )
		    {
		      int i = lvl.int_level;
		      if ( lvl.is_int )
			intcol[ r_cnt ] = lvl.int_level;
		      else if ( lvl.is_dbl )
			intcol[ r_cnt ] = (int)lvl.dbl_level;
		      else
			missing[ r_cnt ] = true;
		    }
		  
		  ++r_cnt;
		  ++rr;
		}


	      if ( is_int_factor )
		df.add( *ff , intcol );
	      else if ( is_dbl_factor )
		df.add( *ff , dblcol );
	      else
		df.add( *ff , strcol );

	      // note -
	      // actually factor/levels should not have any missing data...
	      // but keep above as we use missing[] below
	      
	      ++ff;
	    }
	

	  //
	  // Repeat, as for factors, but now adding actual variables 
	  //
	  
	  vv = tt->second.begin();
	  
	  while ( vv != tt->second.end() )
	    {

	      const retval_var_t & var = vv->first;
		
	      // what type of variable is this?
	      // vv->is_string(), vv->is_double(), vv->is_int()
	      
	      bool var_is_string = var_has_strings.find( var.name ) != var_has_strings.end();
	      bool var_is_double = var_has_doubles.find( var.name ) != var_has_doubles.end();	      
	      
	      std::vector<std::string> strcol;
	      std::vector<double> dblcol;
	      std::vector<int> intcol;
	      std::vector<bool> missing( nrows );
		
	      // note - cases were we might have long long ints...
	      //   make int -> int64_t ? 
	      
	      if      ( var_is_string )
		strcol.resize( nrows );
	      else if ( var_is_double )
		dblcol.resize( nrows );
	      else
		intcol.resize( nrows );
	      
	      // consider all factor/level rows as before (based on same rows file)

	      int r_cnt = 0;

	      std::set<retval_indiv_strata_t>::iterator rr =  rows.begin();
	      while ( rr != rows.end() )
		{
		  
		  // i.e. we are ensuring that we are iterating in the
		  // same order as we previously did for each
		  // variable, so
		  
		  //
		  // does this variable have a non-missing value for
		  // this row/level, for this individual?
		  //
		  
		  std::map<retval_strata_t, std::map<retval_indiv_t,retval_value_t> >::const_iterator
		    yy = vv->second.find( rr->strata );
		    

	          // not present...
	          if ( yy == vv->second.end() )
		    {
		      missing[ r_cnt ] = true;		      		      
		    }
		  else // ...is present as a strata... check for this individual
		    {
		      
		      std::map<retval_indiv_t,retval_value_t>::const_iterator zz = yy->second.find( rr->indiv );
		      
		      // not present...
		      if ( zz == yy->second.end() ) 
			{
			  missing[ r_cnt ] = true;			  
			}
		      else
			{
			  
			  // because of how sqlite stores numeric values, a double may be cast as an int;
			  // therefore, some values for a double variable may in fact be stored as value.i (i.e. if 1.0, etc)
			  // therefore, we need to check for this special case, for data coming from db2retval at least
			  // (this will all be fine if coming from a luna eval() 
			  
			  if      ( var_is_string ) 
			    strcol[ r_cnt ] = zz->second.s ; 			  
			  else if ( var_is_double ) 
			    {
			      // special case
			      if ( zz->second.is_int )
				dblcol[ r_cnt ] = zz->second.i ;
			      else
				dblcol[ r_cnt ] = zz->second.d ;			      
			    }
			  else
			    intcol[ r_cnt ] = zz->second.i ;
			}

		    }
	    
		  // next row/lvl
		  ++r_cnt;
		  ++rr;
   	       }

	      // add this column to the df
	      if ( var_is_string )
		df.add( var.name , strcol , missing );
	      else if ( var_is_double )
		df.add( var.name , dblcol , missing );
	      else
		df.add( var.name , intcol , missing );
	      
	      
	      // next variable
	      ++vv;
   	    }

	  // add this data-frame to the t_list (i.e. all tables for this command)

	  // command (key1)
	  const std::string cmd_name = Helper::sanitize( cmd.name );

	  // label (factors, with _ delim) (key2)
	  std::string table_name = Helper::sanitize( Helper::stringize( table.factors , "_" ));
	  if ( table_name == "" ) table_name = "BL";
	  
	  // add to the return set
	  tables[ cmd_name ][ table_name ] = df;
	  
	  // Next virtual table
	  ++tt;
	}

      // next command
      ++cc;
      
    }

  // all done
  return tables;
  
}



