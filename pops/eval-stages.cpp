
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

#ifdef HAS_LGBM

#include "pops/pops.h"
#include "pops/indiv.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "stats/eigen_ops.h"
#include "db/db.h"

#include "edf/edf.h"

extern logger_t logger;
extern writer_t writer;


pops_indiv_t::pops_indiv_t( param_t & param , 
			    const std::string & file1, 
			    const std::string & file2 )
{
  
  // we'll call this, but should not have any impact
  pops_opt_t::set_options( param );
  
  logger << "  evaluating external stagings in " << file1 << "\n"
	 << "  against " << file2 << "\n";
  
  
  // treat as a 'trainer' (i.e. requires staging)
  trainer = true;

  // i.e. as 'S' will be populated
  has_staging = true;
  
  // build 'S' and 'E' from file1 (and ne)
  // build 'PS' from file2

  if ( ! Helper::fileExists( file1 ) )
    Helper::halt( "cannot open " + file1 );

  if ( ! Helper::fileExists( file2 ) )
    Helper::halt( "cannot open " + file2 );
  
  std::vector<std::string> exss;
  std::ifstream IN1( file1.c_str() , std::ios::in );
  while ( 1 ) 
    {
      std::string ss;
      IN1 >> ss;
      if ( IN1.eof() || IN1.bad() ) break;      
      exss.push_back( ss );      
    }
  IN1.close();
  
  logger << "  read " << exss.size() << " stages (as 'observed') from " << file1 << "\n";
  
  ne = exss.size();

  //
  // "predicted" 
  //
  
  std::vector<std::string> exss2;
  std::ifstream IN2( file2.c_str() , std::ios::in );
  while ( 1 ) 
    {
      std::string ss;
      IN2 >> ss;
      if ( IN2.eof() || IN2.bad() ) break;      
      exss2.push_back( ss );      
    }
  IN2.close();
  
  logger << "  read " << exss2.size() << " stages (as 'predicted') from " << file2 << "\n";
  
  
  const int ne2 = exss2.size();
  
  const int ne_both = ne < ne2 ? ne : ne2 ; 
  
  if ( ne != ne2 ) logger << "  *** warning -- found a different number of epochs across files\n";
  
  if ( ne_both != ne ) logger << "  *** only analysing the first " << ne_both << " epochs (assuming similar starts)\n";
  
  E.clear();
  PS.resize( ne , UNKNOWN );
  S.resize( ne , UNKNOWN );
  
  for (int e=0; e<ne2; e++)
    {

      E.push_back( e );

      if      ( exss[e] == "W" ) S[e] = 0;
      else if ( exss[e] == "R" ) S[e] = 1;
      else if ( exss[e] == "N1" ) S[e] = 2;
      else if ( exss[e] == "N2" ) S[e] = 3;
      else if ( exss[e] == "N3" ) S[e] = 4;
      else S[e] = POPS_UNKNOWN;
      
      if      ( exss2[e] == "W" ) PS[e] = 0;
      else if ( exss2[e] == "R" ) PS[e] = 1;
      else if ( exss2[e] == "N1" ) PS[e] = 2;
      else if ( exss2[e] == "N2" ) PS[e] = 3;
      else if ( exss2[e] == "N3" ) PS[e] = 4;
      else    PS[e] = POPS_UNKNOWN;
	    
    }
  
  eval_stages();
  
  
}


pops_indiv_t::pops_indiv_t( edf_t & edf ,
			    param_t & param , 
			    const std::string & file1 )
  
{

  // track the EDF
  pedf = &edf;

  pops_opt_t::set_options( param );

  logger << "  evaluating external staging in " << file1 << "\n";
  
  // treat as a 'trainer' (i.e. requires staging)
  trainer = true;
  
  // get any staging
  bool has_valid_staging = staging( edf , param );
  
  if ( ! has_valid_staging ) 
    Helper::halt( "no valid staging data found" );
    
  
  // build 'PS' from file1
  if ( ! Helper::fileExists( file1 ) )
    Helper::halt( "cannot open " + file1 );

  std::vector<std::string> exss;
  std::ifstream IN1( file1.c_str() , std::ios::in );
  while ( 1 ) 
    {
      std::string ss;
      IN1 >> ss;
      if ( IN1.eof() || IN1.bad() ) break;      
      exss.push_back( ss );      
    }
  IN1.close();
  
  logger << "  read " << exss.size() << " stages from " << file1 << "\n";
  
  // handle is external is different from internal... ??? 
  // always align from the start in any case.
  
  const int ne2 = exss.size();
  
  const int ne_both = ne < ne2 ? ne : ne2 ; 
  
  if ( ne != ne2 ) logger << "  *** warning -- found a different number of epochs in " << file1 << "\n";

  if ( ne_both != ne ) logger << "  *** only analysing the first " << ne_both << " epochs (assuming similar starts)\n";
  
  PS.resize( ne , UNKNOWN );
  for (int e=0; e<ne2; e++)
    {
      if      ( exss[e] == "W" ) PS[e] = 0;
      else if ( exss[e] == "R" ) PS[e] = 1;
      else if ( exss[e] == "N1" ) PS[e] = 2;
      else if ( exss[e] == "N2" ) PS[e] = 3;
      else if ( exss[e] == "N3" ) PS[e] = 4;
    }
  
  eval_stages();
  
}


void pops_indiv_t::eval_stages() 
{
  
  // track that we have no EDF
  pedf = NULL;

  // so that summarize() doesn't look for P[] to be populated

  pops_opt_t::eval_mode = true;
  
  summarize();
  
}


#endif



