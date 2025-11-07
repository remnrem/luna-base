
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

#include "stats/gpa.h"

#include "param.h"
#include "helper/logger.h"
#include "db/db.h"
#include "stats/eigen_ops.h"
#include "miscmath/crandom.h"
#include "helper/json.h"
#include "miscmath/miscmath.h"

using json = nlohmann::json;

extern logger_t logger;
extern writer_t writer;

gpa_t::gpa_t( param_t & param , const bool prep_mode )
{

  // invoked in one of two modes:
  //   1) prepare a binary data file  [ --gpa-prep ] --> prep_mode == T 
  //   2) apply associaiton models    [ --gpa ]

  // --gpa-prep can either make inputs or take (or make) a JSON specs file
  //    inputs=x,y  dat=z
  //    specs=x     day=z
  //    inputs=x,y  make-specs > specs.json
  
  if ( ! param.has( "make-specs" ) )
    bfile = param.requires( "dat" );
  
  dump_file = param.has( "dump" );
  
  infiles.clear();
  incvars.clear();
  excvars.clear();
  incfacs.clear();
  excfacs.clear();
  incfaclvls.clear();
  excfaclvls.clear();
  incnums.clear();
  excnums.clear();
  incgrps.clear();
  excgrps.clear();
  file2group.clear();
  file2fixed.clear();

  xz_incvars.clear();
  xz_incgrps.clear();

  y_incvars.clear();
  y_incgrps.clear();

  if ( prep_mode )
    {

      //
      // spec from JSON file? (allow both spec= and specs= forms)
      //
      
      if ( param.has( "spec" ) )
	{
	  std::vector<std::string> specs = param.strvector( "spec" );
	  for (int i=0; i<specs.size(); i++) parse( specs[i] );
	}
      
      if ( param.has( "specs" ) )
	{
	  std::vector<std::string> specs = param.strvector( "specs" );
	  for (int i=0; i<specs.size(); i++) parse( specs[i] );
	}


      //
      // additional, from command line:
      //

      std::set<std::string> files;
      if ( param.has( "inputs" ) )
	{
	  files = param.strset( "inputs" );
	  if ( files.size() == 0 ) return;
	}

      std::set<std::string>::const_iterator ff = files.begin();
      while ( ff != files.end() )
	{
	  std::vector<std::string> tok = Helper::parse( *ff , "|" );
	  if ( tok.size() == 1 )
	    Helper::halt( "expecting inputs=file|grp or file|grp|fac1|fac2|..." ); 
	  else if ( tok.size() == 2 )
	    {
	      const std::set<std::string> empty;
	      infiles[ tok[0] ] = empty;
	      file2group[ tok[0] ] = tok[1];
	    }
	  else if ( tok.size() > 2 )
	    {
	      file2group[ tok[0] ] = tok[1];

	      // to ensure that infiles gets populated, even if empty (i.e. if only fixed fac)
	      // so it gets loaded downstream
	      
	      const std::set<std::string> empty;
	      infiles[ tok[0] ] = empty;

	      for (int j=2;j<tok.size(); j++)
		{
		  
		  // add as fixed factor (if starting with inputs=res/spec.n2|spec|F|CH|SS/N2
		  // else, means to look for this col in the input
		  const bool is_fixed = tok[j].find( "/" ) != std::string::npos ;
		  
		  if ( is_fixed ) {
		    std::vector<std::string> tok2 = Helper::parse( tok[j] , "/" );
		    if ( tok2.size() == 2 ) 
		      file2fixed[ tok[0] ][ tok2[0] ] = tok2[1];
		    else
		      Helper::halt( "bad format for factor specification: " + *ff );
		  } else { 
		    infiles[ tok[0] ].insert( tok[j] );
		  }
		  
		}
	    }
	  
	  ++ff;
	}
    }
  

  //
  // make JSON? (and then quit) 
  //

  if ( param.has( "make-specs" ) )
    {
      if ( file2group.size() == 0 )
	Helper::halt( "no inputs specified" );
      if ( param.has( "spec" ) || param.has( "specs" ) )
	Helper::halt( "cannot combine make-specs and specs" );

      logger << "  writing .json specification (from inputs) to standard output\n";
      
      // start JSON
      std::cout << "{\n";

      // inputs
      std::cout << "  \"inputs\": [\n";
      
      // iterate over files
      int idx = 0;
      std::map<std::string,std::string>::const_iterator ff = file2group.begin(); 
      while ( ff != file2group.end() )
	{

	  //
	  // get vars from header row
	  //

	  if ( ! Helper::fileExists( ff->first ) )
	    Helper::halt( "could not open " + ff->first );
	  
	  std::ifstream IN1( ff->first.c_str() , std::ios::in );
	  std::string hdr;
	  Helper::safe_getline( IN1 , hdr );

	  // bad/empty file?
	  if ( IN1.eof() || hdr == "" )
	    {
	      IN1.close();
	      Helper::halt( ff->first + " is a bad/empty file" );
	    }

	  IN1.close();
	  
	  // assume tab-delimited
	  std::vector<std::string> tok = Helper::parse( hdr , "\t" );
	  if ( tok.size() == 0 )
	    Helper::halt( ff->first + " is a bad/empty file" );

	  if ( tok[0] != "ID" )
	    Helper::halt( "column 1 of " + ff->first + " is not ID" );
	  
	  int nfac = 0;
	  std::vector<std::string> labs;
	  
	  const std::set<std::string> & facs = infiles.find( ff->first )->second;

	  for (int i=1; i<tok.size(); i++)
	    {
	      if ( facs.find( tok[i] ) != facs.end() )
		++nfac;
	      else
		labs.push_back( tok[i] );	      
	    }

	  if ( labs.size() == 0 )
	    Helper::halt( "no non-factor variables in " + ff->first );

	  if ( nfac != facs.size() )
	    Helper::halt( "found " + Helper::int2str( nfac )
			  + " but expecting " + Helper::int2str( (int)facs.size() )
			  + " factors in " + ff->first );

	  
	  
	  logger << "  found " << labs.size()
		 << " variables and " << nfac
		 << " factors from " << ff->first << "\n";
	    
	  //
	  // write JSON
	  //
	  
	  std::cout << "    {\n";
	  
	  // group
	  std::cout << "      \"group\": \"" << ff->second << "\",\n";

	  // file
	  std::cout << "      \"file\": \"" << ff->first << "\",\n";
	  
	  // (opt) facs
	  if ( facs.size() != 0 )
	    {
	      std::cout << "      \"facs\": [";
	      std::set<std::string>::const_iterator ss = facs.begin();
	      while ( ss != facs.end() )
		{
		  int idx2 = 0;
		  while ( ss != facs.end() )
		    {
		      std::cout << " \"" << *ss << "\"";
		      if ( ++idx2 != facs.size() ) std::cout << ",";
		      ++ss;
		    }
		  std::cout << " ],\n";
		}
	    }
	  
	  // (opt) fixed-factors?
	  const std::map<std::string,std::string> & fixed = file2fixed.find( ff->first )->second;
	  if ( fixed.size() != 0 )
	    {
              std::cout << "      \"fixed\": [";
              std::map<std::string,std::string>::const_iterator ss = fixed.begin();
              while ( ss != fixed.end() )
                {
                  int idx2 = 0;
                  while ( ss != fixed.end() )
                    {
                      std::cout << " {\"" << ss->first << "\""
				<< ": \"" << ss->second << "\" }";
                      if ( ++idx2 != fixed.size() ) std::cout << ",";
                      ++ss;
                    }
                  std::cout << " ],\n";
                }
            }
	  
	  // vars
	  std::cout << "      \"vars\": [";   
	  
	  std::vector<std::string>::const_iterator vv = labs.begin();
	  int idx2 = 0;
	  while ( vv != labs.end() )
	    {
	      std::cout << " \"" << *vv << "\"";
	      if ( ++idx2 != labs.size() ) std::cout << ",";
	      ++vv;
	    }
	  std::cout << " ]\n";
	  
	  // close out
	  std::cout << "    }";
	  if ( ++idx < file2group.size() ) std::cout << ",";
	  
	  // next file
	  std::cout << "\n";
	  
	  ++ff;
	}
      
      
      // end-inputs
      std::cout << "  ]\n";

      // close out
      std::cout << "}\n";

      // all done
      return;
      
    }
  
  //
  // include/exclude vars (allow for partial population during spec JSON stage)
  //
  
  if ( param.has( "vars" ) )
    incvars = Helper::combine( incvars , param.strset( "vars" ) );
      
  if ( param.has( "xvars" ) )
    excvars = Helper::combine( excvars , param.strset( "xvars" ) );

  //
  // include vars based on presence/absence of a factor
  //

  if ( param.has( "facs" ) )
    incfacs = Helper::combine( incfacs, param.strset( "facs" ) );
  if ( param.has( "xfacs" ) )
    excfacs = Helper::combine( excfacs, param.strset( "xfacs" ) );


  //
  // include/exclude based on (file) group
  //
  
  if ( param.has( "grps" ) )
    incgrps = Helper::combine( incgrps, param.strset( "grps" ) );
  
  if ( param.has( "xgrps" ) )
    excgrps = Helper::combine( excgrps, param.strset( "xgrps" ) );


  //
  // other options not supported in prep-mode (nvars & faclvls)
  //
  
  //
  // fac/lvl pairs (for now, can only do on read)
  //

  if ( prep_mode )
    {
      if (  param.has( "faclvls" ) ||  param.has( "xfaclvls" ) )
        Helper::halt( "cannot specify faclvls/xfaclvls with --gpa-prep" );
    }

  
  if ( param.has( "faclvls" ) )
    {
      // faclvl=CH/CZ|FZ|OZ,B/SIGMA|GAMMA
      std::vector<std::string> tok = param.strvector( "faclvls" );
      for (int i=0; i<tok.size(); i++)
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[i] , "/" );
	  if ( tok2.size() != 2 ) Helper::halt( "bad fac/lvl1|lvl2 specification" );
	  incfaclvls[ tok2[0] ] = Helper::vec2set( Helper::parse( tok2[1] , "|" ) );
	}
    }

  if ( param.has( "xfaclvls" ) )
    {
      // xfaclvl=CH/CZ|FZ|OZ,B/SIGMA|GAMMA
      std::vector<std::string> tok = param.strvector( "xfaclvls" );
      for (int i=0; i<tok.size(); i++)
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[i] , "/" );
	  if ( tok2.size() != 2 ) Helper::halt( "bad fac/lvl1|lvl2 specification" );
	  excfaclvls[ tok2[0] ] = Helper::vec2set( Helper::parse( tok2[1] , "|" ) );
	}
    }


  // lvars/xlvars only allowed in read/run mode
  // nvars/xnvars only allowed in read/run mode
  if ( prep_mode )
    {
      if (  param.has( "nvars" ) ||  param.has( "xnvars" ) )
	Helper::halt( "cannot specify nvars/xnvars with --gpa-prep" );

      if (  param.has( "lvars" ) ||  param.has( "lnvars" ) )
	Helper::halt( "cannot specify lvars/xlvars with --gpa-prep" );	    
    }


  
  //
  // include/exclude vars (allow for partial population during spec JSON stage
  //
  
  if ( param.has( "lvars" ) )
    inclvars = param.strset( "lvars" );
  
  if ( param.has( "xlvars" ) )
    exclvars = param.strset( "xvars" );  
  
  if ( param.has( "nvars" ) )
    {
      std::vector<std::string> tok = param.strvector( "nvars" );
      for (int i=0; i<tok.size(); i++)
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[i] , "-" );
	  if ( tok2.size() == 1 )
	    {
	      int n1; // --> convert to 0-base
	      if ( Helper::str2int( tok2[0] , &n1 ) )
		{
		  if ( n1 >= 1 ) incnums.push_back( std::make_pair( n1-1 , n1-1 ) );		    
		}
	      else Helper::halt( "bad values for nvars" );
	    }
	  else if ( tok2.size() == 2 )
	    {
	      int n1, n2;
              if ( Helper::str2int( tok2[0] , &n1 ) && Helper::str2int( tok2[1] , &n2 ) )
		{
		  if ( n1 >= 1 && n2 >= 1 )
		    incnums.push_back( std::make_pair( n1-1 , n2-1 ) );
		}
              else Helper::halt( "bad values for nvars" );
	    }
	  else Helper::halt( "bad values for nvars" );

	}
    }
  
  if ( param.has( "xnvars" ) )
    {
      std::vector<std::string> tok = param.strvector( "xnvars" );
      for (int i=0; i<tok.size(); i++)
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[i] , "-" );
	  if ( tok2.size() == 1 )
	    {
	      int n1;
	      if ( Helper::str2int( tok2[0] , &n1 ) )
		{
		  if ( n1 >= 1 ) 
		    excnums.push_back( std::make_pair( n1-1 , n1-1 ) );
		}
	      else Helper::halt( "bad values for xnvars" );
	    }
	  else if ( tok2.size() == 2 )
	    {
	      int n1, n2;
              if ( Helper::str2int( tok2[0] , &n1 ) && Helper::str2int( tok2[1] , &n2 ) )
		{
		  if ( n1 >= 1 && n2 >= 1 )
		    excnums.push_back( std::make_pair( n1-1 , n2-1 ) );
		}
              else Helper::halt( "bad values for xnvars" );
	    }
	  else Helper::halt( "bad values for xnvars" );

	}
    }


  //
  // any other X, Y, Z vars?  track separately: if requested, these will
  // be included, unless explicitly excluded;   i.e. they will not feature
  // as part of the AND logic for requiring vars/grps etc
  //

  if ( param.has( "X" ) ) xz_incvars = Helper::combine( xz_incvars , param.strset( "X" ) );
  if ( param.has( "Z" ) ) xz_incvars = Helper::combine( xz_incvars , param.strset( "Z" ) );
  if ( param.has( "Y" ) ) y_incvars = Helper::combine( y_incvars , param.strset( "Y" ) );

  if ( param.has( "Xg" ) ) xz_incgrps = Helper::combine( xz_incgrps , param.strset( "Xg" ) );
  if ( param.has( "Zg" ) ) xz_incgrps = Helper::combine( xz_incgrps , param.strset( "Zg" ) );
  if ( param.has( "Yg" ) ) y_incgrps = Helper::combine( y_incgrps , param.strset( "Yg" ) );

  //
  // if strata flag (for stats) then ensure those variables are pulled to
  //

  // if ( param.has( "S" ) ) s_incvars = Helper::combine( s_incvars , param.strset( "S" ) );
  // if ( param.has( "Sg" ) ) s_incgrps = Helper::combine( s_incgrps , param.strset( "Sg" ) );
  //  n.b. means s_incgrps() is not used/needed now
  //       we *only* have stratifying vars for strata (stats-only) for subset
  
  if ( param.has( "strata" ) )
    {
      if ( ! param.has( "stats" ) )
	Helper::halt( "can only specify strata with stats" );
      s_incvars = Helper::combine( s_incvars , param.strset( "strata" ) );
    }

  if ( param.has( "subset" ) ) { // need to strip any '-' from starts
    std::set<std::string> s1 = param.strset( "subset" );
    std::set<std::string> s2;
    std::set<std::string>::const_iterator ss = s1.begin();
    while ( ss != s1.end() )
      {
	if ( ss->substr(0,1) == "-" || ss->substr(0,1) == "+" ) s2.insert( ss->substr(1) );
	else s2.insert( *ss );
	++ss;
      }
    s_incvars = Helper::combine( s_incvars , s2 ) ;
  }
  
  //
  // criteria to drop bad/empty cols (default, at least 5 non-missing,
  //   at least 5% of obs w/ values), although in --gpa-prep mode, we
  //   take everything by default
  //

  n_req = param.has( "n-req" ) ? param.requires_int( "n-req" ) : ( prep_mode ? 0 : 5 ) ;
  n_prop = param.has( "n-prop" ) ? param.requires_dbl( "n-prop" ) : ( prep_mode ? 0 : 0.05 ) ; 
  
  retain_cols = param.has( "retain-cols" ); 
  retain_rows = param.has( "retain-rows" );
  if ( prep_mode && retain_rows ) Helper::halt( "cannot use retain-rows in --gpa-prep" );

 
  
  //
  // misc options
  //

  verbose = param.has( "verbose" ) ? param.yesno( "verbose" ) : false;

  show_xfacs = param.has( "X-factors" ) ? param.yesno( "X-factors" ) : false; 
    
  //
  // read data / make binary file? 
  //
  
  
  if ( prep_mode )
    prep();
  else
    {
      
      // intention to perform association is denoted by the presence
      // of X or request for perms
      
      const bool request_assoc = param.has( "nreps" )
	|| param.has( "X" ) || param.has( "all-by-all" );
      
      if ( request_assoc
	   && ( param.has( "retain-cols" ) || param.has( "retain-rows" ) ) )
	Helper::halt( "can only use retain-cols or retain-rows when not running association (no X)" );

      // read data in  ( and this does filtering of columns)
      
      logger << "  reading binary data from " << bfile << "\n";

      read();
      
      //
      // store original N 
      //

      const int ni_orig1 = X.rows();


      // secondarily, subset to a smaller # of rows (e.g. case-only analysis)

      if ( param.has( "subset" ) || param.has( "inc-ids" ) || param.has( "ex-ids")  )
	{
	  std::set<std::string> sub_ids, exc_ids, sub_cols;
	  
	  if ( param.has( "inc-ids" ) ) sub_ids = param.strset( "inc-ids" );
	  if ( param.has( "ex-ids" ) ) exc_ids = param.strset( "ex-ids" );
	  if ( param.has( "subset" ) ) sub_cols = param.strset( "subset" );
	  std::set<int> rows;
	  std::map<int,bool> cols;
	  const int ni = X.rows();
	  const int nv = X.cols();

	  const bool has_inc = sub_ids.size(); // only include these
	  const bool has_exc = exc_ids.size(); // do not include these
	  	  
	  // build ID rows
	  for (int i=0; i<ni; i++)
	    {
	      if ( ( ! has_inc) || sub_ids.find( ids[i] ) != sub_ids.end() ) 
		{
		  if ( ( ! has_exc ) || exc_ids.find( ids[i] ) == exc_ids.end() )
		    rows.insert( i );
		}	      
	    }

	  
	  // build cols (i.e. searching for non-null value to include, not NaN or missing)
	  // allows for each term to be a pos or neg match subset=-MALE implies MALE == 0 (-->F)
	  //  where subset=MALE or subset=+MALE implies --> M
	  std::map<std::string,bool> sub_cols_found;
	  std::set<std::string>::const_iterator ss = sub_cols.begin();
	  while ( ss != sub_cols.end() ) { sub_cols_found[*ss] = false; ++ss; }
	  for (int j=0; j<nv; j++)
	    {	      
	      if ( sub_cols.find( vars[j] ) != sub_cols.end() )
		{
		  cols[j] = true;
		  sub_cols_found[ vars[j] ] = true;
		}
	      else if ( sub_cols.find( "+" + vars[j] ) != sub_cols.end() )
		{
		  cols[j] = true;
		  sub_cols_found[ "+" + vars[j] ] = true;
		}
	      else if ( sub_cols.find( "-" + vars[j] ) != sub_cols.end() )
		{
		  cols[j] = false;
		  sub_cols_found[ "-" + vars[j] ] = true;
		}	      
	    }

	  // did we find all requsted sub-sets?
	  std::map<std::string,bool>::const_iterator ff = sub_cols_found.begin();
	  while ( ff != sub_cols_found.end() )
	    {
	      if ( ! ff->second )
		Helper::halt( "could not find subset variable " + ff->first );
	      ++ff;
	    }
	  
	  // do the subsetting
	  subset( rows, cols );
	}

       
      //
      // store original N (post any selection 
      //

      const int ni_orig2 = X.rows();

     
      //
      // select predictors (X) and covariates (Z) - assume everything else is
      // a DV (i.e. sleep metric) unless explicitly told so; allow these to be 
      // specified either as variables (X,Z) or groups (Xg, Zg)
      //

	
      if ( param.has( "X" ) || param.has( "Xg" ) )
	{	  
	  ivs.clear();
	  std::set<std::string> v     = param.strset( "X" );
	  std::set<std::string> vgrps = param.strset( "Xg" );
	  std::set<std::string> found;
	  std::set<std::string> gfound;
	  const int nv = vars.size();
	  for (int j=0; j<nv; j++) {
	    if ( v.find( vars[j] ) != v.end() ) {
	      ivs.push_back( j );
	      found.insert( vars[j] );
	    }	      
	    if ( vgrps.find( var2group[ vars[j] ] ) != vgrps.end() ) {
	      ivs.push_back( j );
	      gfound.insert( vars[j] );
	    }	      
	  }
	  
	  if ( found.size() < v.size() ) {
	    std::set<std::string>::const_iterator vv = v.begin();
	    while ( vv != v.end() ) {
	      if ( found.find( *vv ) == found.end() )
		logger << "  *** warning - could not find " << *vv << "\n";
	      ++vv;
	    }
	  }		  
	  
	  if ( gfound.size() < vgrps.size() ) {
	      std::set<std::string>::const_iterator vv = vgrps.begin();
	      while ( vv != vgrps.end() ) {
		if ( gfound.find( *vv ) == gfound.end() )
		  logger << "  *** warning - could not find " << *vv << "\n";
		++vv;
	      }
	  }		  	  
	}


      if ( param.has( "Z" ) || param.has( "Zg" ) )
	{	  
	  cvs.clear();
	  std::set<std::string> v     = param.strset( "Z" );
	  std::set<std::string> vgrps = param.strset( "Zg" );
	  std::set<std::string> found;
	  std::set<std::string> gfound;
	  const int nv = vars.size();
	  for (int j=0; j<nv; j++) {
	    if ( v.find( vars[j] ) != v.end() ) {
	      cvs.push_back( j );
	      found.insert( vars[j] );
	    }
	    
	    if ( vgrps.find( var2group[ vars[j] ] ) != vgrps.end() ) {
	      cvs.push_back( j );
	      gfound.insert( vars[j] );
	    }
	  }
	  
	  if ( found.size() < v.size() ) {
	    std::set<std::string>::const_iterator vv = v.begin();
	    while ( vv != v.end() ) {
	      if ( found.find( *vv ) == found.end() )
		logger << "  *** warning -  could not find " << *vv << "\n";
	      ++vv;
	    }
	  }		  
	  
	  if ( gfound.size() < vgrps.size() ) {
	    std::set<std::string>::const_iterator vv = vgrps.begin();
	    while ( vv != vgrps.end() ) {
	      if ( gfound.find( *vv ) == gfound.end() )
		logger << "  *** warning - could not find " << *vv << "\n";
	      ++vv;
	    }
	  }		  
	  
	}

      
      //
      // Y is everything else that is left
      //  (i.e. post any prior col selection)
      //  unless Yg is specified too: (to select Y 
      //  based on groups
      //
      
      const std::set<std::string> ygroups = param.strset( "Yg" );
      const std::set<std::string> yvars = param.strset( "Y" );
      const bool has_yspecified = ygroups.size() != 0 || yvars.size() != 0;
      
      //
      // if 'all-by-all' added, then also set X == Y
      //
      
      const bool all_by_all = param.has( "all-by-all" );
      
      if ( all_by_all && ( param.has( "X" ) || param.has( "Xg" ) ) )
	Helper::halt( "cannot specify X (Xg) and all-by-all together" );

      
      std::set<int> v1 = Helper::vec2set( ivs );
      std::set<int> v2 = Helper::vec2set( cvs );
      
      const int nv = vars.size();      

      if ( all_by_all )
	{
	  for (int j=0; j<vars.size(); j++)
            if ( v2.find( j ) == v2.end() ) // not a covariate
	      {
		if ( ( ! has_yspecified )
		     || yvars.find( vars[j] ) != yvars.end()
		     || ygroups.find( var2group[ vars[j] ] ) != ygroups.end() )
		  {
		    ivs.push_back( j );
		    dvs.push_back( j );
		  }
	      }
	  
	  logger << "  selected " 
		 << cvs.size() << " Z vars, implying "
		 << dvs.size() << " X and Y vars (given 'all-by-all')\n";
	  
	}
      
      else // standard, assuming that X has been specified (and is now in v1)
	{
	  for (int j=0; j<vars.size(); j++)
	    if ( v1.find( j ) == v1.end() && v2.find( j ) == v2.end() )
	      {
		if ( ( ! has_yspecified )
                     ||	yvars.find( vars[j] ) != yvars.end()
                     ||	ygroups.find( var2group[ vars[j] ] ) !=	ygroups.end() )
		  {
		    dvs.push_back( j );
		  }
	      }
	  logger << "  selected " 
		 << ivs.size() << " X vars & "
		 << cvs.size() << " Z vars, implying "
		 << dvs.size() << " Y vars\n";
	  
	}

      //
      // internal/debug option: set all NA variables to 0? (e.g. in test mode) 
      //

      if ( param.has( "force-zeros" ) )
	{
	  for (int i=0; i<X.rows(); i++)
	    for (int j=0; j<X.cols(); j++)
	      if ( std::isnan( X(i,j) ) ) X(i,j) = 0;	  	  
	}

      //
      // kNN imputation of missing points (for DVs only)? 
      //
      
      if ( param.has( "knn" ) )
	knn_imputation( param );
      
      //
      // drop any null cols in this (potentially subsetted) dataset
      //

      if ( ! retain_cols ) 
	drop_null_columns();
      
      
      //
      // now do QC (on DVs only) - i.e. keeps 
      //

      double winsor_th = -9; // default = N
      if ( param.has( "winsor" ) )
	{
	  if ( param.value( "winsor" ) == "F"
	       || param.value( "winsor" ) == "N"
	       || param.value( "winsor" ) == "0" ) 
	    winsor_th = -9; // i.e. none
	  else
	    {
	      winsor_th = param.requires_dbl( "winsor" ) ;
	      if ( winsor_th < 0 || winsor_th > 0.2 )
		Helper::halt( "winsor must be set between 0 and 0.2" ); 
	    }
	}
      
      // can skip QC with qc=F option
      if ( ( ! param.has( "qc" ) ) || param.yesno( "qc" ) )
	qc( winsor_th , param.has( "stats" ) );


      
      
      //
      // optionally dump variables and/or manifest?
      //

      if ( param.has( "dump" ) )
	{
	  dump();
	  return;
	}
      
      if ( param.has( "manifest" ) )
	{
	  manifest();
	  return;
	}

      if ( param.has( "stats" ) )
	{

	  if ( param.has( "strata" ) )
	    {
	      std::vector<std::string> sx = param.strvector( "strata" );
	      stats( &sx );
	    }
	  else
	    stats();
	  
	  return;
	}
      
      if ( param.has( "summary" ) || param.has( "summarize" ) || param.has( "desc" ) )
	{
	  summarize();
	  return;
	}

    
      //
      // drop any stratifying variables (note: currently, s_incgrps() is
      //  not used, can be dropped; likely these were already dropped - i.e.
      //  due to invariance (i.e. post selection)... but double-check here
      //
      
      if ( s_incvars.size() )
	{
	  bool to_drop = false;
	  std::vector<int> retain;
	  const int nv = X.cols();
	  for (int j=0; j<nv; j++)
	    {
	      if ( s_incvars.find( vars[j] ) == s_incvars.end() )
		retain.push_back( j );
	      else
		to_drop = true; 
	    }	  

	  if ( to_drop )
	    {
	      Eigen::MatrixXd X1 = X( Eigen::all , retain );
	      X = X1;
	      logger << "  dropped " << s_incvars.size() << " stratifying variables\n";
	    }
	}
            
      //
      // run association tests
      //
           
      const bool run_assoc =  dvs.size() != 0 && ivs.size() != 0 ; 

      if ( ! run_assoc ) 
	logger << "  not fitting any association models, as no dependent, or no predictor vars specified\n";
      
      if ( run_assoc )
	{	  

	  if ( X.rows() <= 2 )
	    {
	      logger << "  *** only " << X.rows() << " left, cannot fit linear models\n";
	      return ;
	    }
	  
	  // options: permutations?
	  nreps = param.has( "nreps" ) ? param.requires_int( "nreps" ) : 0 ;
	  if ( nreps < 0 ) Helper::halt( "nreps must be positive" );

	  show_progress = param.has( "progress" ) ? param.yesno( "progress" ) : true ; 
	  
	  // level of multiple-test correction
	  correct_all_X = param.has( "adj-all-X" ) ? param.yesno( "adj-all-X" ) : false;

	  // adjustments
	  // BONF (=Bonferroni single-step adjusted p-values),
	  // HOLM (=Holm (1979) step-down adjusted p-values),
	  // FDR_BH (=Benjamini & Hochberg (1995) step-up FDR control),
	  // FDR_BY (=Benjamini & Yekutieli (2001) step-up FDR control).

	  // for B&H FDR by default
	  adj_fdr_bh = param.has( "fdr" ) ? param.yesno( "fdr" ) : true; 
	  adj_holm = param.has( "holm" );
	  adj_bonf = param.has( "bonf" );
	  adj_fdr_by = param.has( "fdr-by" );
	  if ( param.has( "adj") )
	    adj_fdr_bh = adj_bonf = adj_holm = adj_fdr_by = true;	  
	  adj_any = adj_holm || adj_fdr_bh || adj_bonf || adj_fdr_by ;
	  
	  // # tests requested
	  int ntests = 0;
	  for (int i=0; i<ivs.size(); i++)
	    for (int j=0; j<dvs.size(); j++)
	      if ( ivs[i] != dvs[j] ) ++ntests;


	  // optional: comparison against a template

	  if ( param.has( "comp" ) )
	    {
	      if ( param.empty( "comp" ) )
		Helper::halt( "comp requires a filename" );
	      comps.load( param.value( "comp" ) );
	      do_comparisons = true;
	    }
	  else do_comparisons = false;
	  
	  
	  //
	  // report tests/Ns (and flag potential problems)
	  //
	  
	  logger << "  " << ntests << " total tests specified\n";
	  
	  if ( ni_orig2 != ni_orig1 )
	    {
	      logger << "  analysis of " << X.rows() << " from " << ni_orig2 << " selected from " << ni_orig1 << " original individuals\n";
	      if ( X.rows() < ni_orig2 / 2 ) logger << "  *** warning - dropping over half of selected individuals\n"; 
	    }
	  else
	    {
	      logger << "  analysis of " << X.rows() << " of " << ni_orig1 << " individuals\n";
	      if ( X.rows() < ni_orig1 / 2 ) logger << "  *** warning - dropping over half of all individuals\n"; 
	    }
	  
	  if ( X.rows() < 10 ) 
	    logger << "  *** warning - running w/ a very small case N - perhaps check for missing values?\n";
	  

	  //
	  // do the actual work 	  
	  //
	  
	  if ( correct_all_X ) 
	    {
	      logger << "  adjusting for multiple tests across all X variables\n";
	      if ( nreps ) 
		logger << "  performing association tests w/ "
		       << nreps << " permutations... (may take a while)\n"; 
	      run();
	    }
	  else
	    {

	      logger << "  adjusting for multiple tests only within each X variable ('adj-all-X' to adjust across all X)\n";

	      if ( nreps ) 
		logger << "  performing association tests w/ "
		       << nreps << " permutations... (may take a while)\n";	  
	      run1X();
	    }
	  
	}
    
    }

}


void gpa_t::prep()
{

  std::map<std::string,int> id2slot;
  std::map<std::string,int> var2slot;
  std::map<int,std::string> slot2var;

  faclvl.clear();
  basevar.clear();
  var2group.clear();

  logger << "\n  preparing inputs...\n";
  
  std::map<int,std::map<int,double> > D; // nb. D[ var ][ ind ] --> value
    
  // iterate over all files
  std::map<std::string,std::set<std::string> >::const_iterator ff = infiles.begin();
  while ( ff != infiles.end() )
    {
      
      // was this group excluded?
      if ( incgrps.size() && incgrps.find( file2group[ ff->first ] ) == incgrps.end() )
	{
	  logger << "  -- " << ff->first << ": skipping due to grps requirement\n";
	  ++ff;
	  continue;
	}

      if ( excgrps.find( file2group[ ff->first ] ) != excgrps.end() )
	{
	  logger << "  -- " << ff->first << " skipping due to xgrps requirement\n";
	  ++ff;
	  continue;
	}

      
      // file-specific counts
      std::set<std::string> file_ids, file_bvars, file_evars;
            
      // named factors for this file
      const std::set<std::string> & facs = ff->second;

      // any file-specific variable inc/exc lists
      const std::set<std::string> & fincvars = file2incvars[ ff->first ];
      const std::set<std::string> & fexcvars = file2excvars[ ff->first ];

      // any file-specific fixed fac/lvls?
      const std::map<std::string,std::string> & fixed = file2fixed[ ff->first ];

      // any file-specific aliasing?
      const std::map<std::string,std::string> & aliases = file2var2alias[ ff->first ];

      // any file-specific/variable-specific mapping of strings -> numeric
      const std::map<std::string,std::map<std::string,double> > & mappings = file2var2mapping[ ff->first ];
      
      
      //
      // if including/excluding based on presence of factors, we can
      // decide here whether to skip the entire file or no
      //
      
      // if requiring factors, they must all be present
      if ( incfacs.size() )
	{
	  const int req_matched =  Helper::nmatches( incfacs , facs );
	  const int obs_matched =  Helper::nmatches( facs , incfacs );
	  if ( req_matched != facs.size() || req_matched != obs_matched || req_matched != incfacs.size() )
	    {
	      logger << "  -- " << ff->first << ": skipping due to facs requirement\n";
	      ++ff;
	      continue;
	    }
	}

      if ( excfacs.size() )
	{
          const int req_matched =  Helper::nmatches( excfacs , facs );
          if ( req_matched == excfacs.size() && req_matched == facs.size() )
            {
	      logger <<	"  -- " << ff->first << ": skipping due to xfacs requirement\n";
              ++ff;
              continue;
            }
        }
      
      
      //
      // read int
      //

      if ( ! Helper::fileExists( Helper::expand( ff->first ) ) )
	{
	  logger << "  -- " << ff->first << ": skipping, could not open file\n";
	  ++ff;
	  continue;
	}
      
      
      std::ifstream IN1( Helper::expand( ff->first ).c_str() , std::ios::in );

      //
      // get header
      //

      std::string hdr;      
      Helper::safe_getline( IN1 , hdr );

      // bad/empty file?
      if ( IN1.eof() || hdr == "" )
        {
          IN1.close();
	  logger << "  -- " << ff->first << ": skipping, bad/empty file\n";
	  ++ff;
	  continue;
        }

      // assume tab-delimited
      std::vector<std::string> tok = Helper::parse( hdr , "\t" );
      if ( tok.size() < 2 ) continue;
      
      int id_col = -1;
      std::vector<bool> col( tok.size() , false );
      std::map<std::string,int> fac2slot;
      
      std::vector<std::string> tok2;
      for (int j=0; j<tok.size(); j++)
	{
	  
	  // aliases?
	  if ( aliases.find( tok[j] ) != aliases.end() )
	    tok[j] = aliases.find( tok[j] )->second;

	  // ID col?
	  if ( tok[j] == "ID" )
	    {
	      id_col = j;
	      continue;
	    }

	  // factor?
	  if ( facs.find( tok[j] ) != facs.end() )
	    {
	      fac2slot[ tok[j] ] = j;
	      //allfacs.insert( tok[j] );
	      continue;
	    }
	  
	  // skip? (global options)
	  if ( incvars.size() != 0 && incvars.find( tok[j] ) == incvars.end() ) continue;
	  if ( excvars.find( tok[j] ) != excvars.end() ) continue;

	  // skip? (file-local options)
          if ( fincvars.size() != 0 && fincvars.find( tok[j] ) == fincvars.end() ) continue;
          if ( fexcvars.find( tok[j] ) != fexcvars.end() ) continue;
	  
	  // if still here, means we want to read this var
	  tok2.push_back( tok[j] );
	  col[j] = true;
	  file_bvars.insert( tok[j] );
	}

      // no ID col?
      if ( id_col == -1 )
	Helper::halt( "no ID column for " + ff->first );

      // no vars to read?
      if ( tok2.size() == 0 )
	{	  
	  logger << "  -- " << ff->first << ": skipping, no selected (non-factor) variables\n";
	  ++ff;
	  continue;
	}
      
      // not all factors found?
      if ( facs.size() != fac2slot.size() )
	Helper::halt( "not all factors found for " + ff->first );

      //
      // build any fixed faclvl string once
      //

      std::string fixed_str;
      std::map<std::string,std::string>::const_iterator xx = fixed.begin();
      while ( xx != fixed.end() )
	{
	  fixed_str += "_" + xx->first + "_" + xx->second;
	  ++xx;
	}


      //
      // read rows ( to find unique IDs and unique VAR+FACLVL combos )
      //

      
      while ( 1 )
	{

	  std::string dat;

	  Helper::safe_getline( IN1 , dat );
	  
	  // all done?
	  if ( IN1.eof() || dat == "" )
	    {
	      IN1.close();
	      break;
	    }

	  std::vector<std::string> dtok = Helper::parse( dat , "\t" );
	  if ( dtok.size() != tok.size() )
	    {
	      logger << "  *** reading file " << ff->first << "\n";
	      logger << "  *** expecting " << tok.size() << " fields, but found " << dtok.size() << "\n";
	      logger << "  *** line [" << dat << "]\n";
	      Helper::halt( "bad line - col # doesn't match header" );
	    }

	  // a new ID?
	  const std::string & id = dtok[ id_col ] ;
	  file_ids.insert( id );
	  if ( id2slot.find( id ) == id2slot.end() )
	    {
	      const int n1 = id2slot.size(); 
	      id2slot[ id ] = n1;
	    }
	  
	  // construct the faclvl for this row, initiating w/ any fixed fac/lvls
	  std::string fl = fixed_str;
	  std::map<std::string,std::string> ffll = fixed;
	  if ( facs.size() )
	    {
	      // at least one faclvl
	      std::map<std::string,int>::const_iterator ff = fac2slot.begin();
	      while ( ff != fac2slot.end() )
		{
		  fl += "_" + tok[ ff->second ] + "_" + dtok[ ff->second ];
		  // store in fac->lvl map  
		  ffll[ tok[ ff->second ] ] = dtok[ ff->second ]; 
		  ++ff;
		}	      
	    }

	  // register each var (w/ unique faclvl if specified)
	  // and add values
	  for (int j=0; j<tok.size(); j++)
	    if ( col[j] )
	      {
		std::string expand_vname = tok[j] + fl;
		file_evars.insert( expand_vname );
		// have we seen this before? if not, add
		if ( var2slot.find( expand_vname ) == var2slot.end() )
		  {
		    const int n1 = var2slot.size();
		    var2slot[ expand_vname ] = n1;
		    slot2var[ n1 ] = expand_vname;
		    
		    faclvl[ expand_vname ] = ffll; // and store
		    basevar[ expand_vname ] = tok[j];
		    var2group[ expand_vname ] = file2group[ ff->first ];
		  }
		

		// check we haven't already seen this:		
		std::map<int,std::map<int,double> >::const_iterator dd = D.find( var2slot[ expand_vname ] );
		if ( dd != D.end() )
		  {
		    const std::map<int,double> & D2 = dd->second;
		    if ( D2.find( id2slot[ id ] ) != D2.end() )
		      logger << "  *** warning - repeated instances of " << expand_vname << " for " << id
			     << " (" << D2.find( id2slot[ id ] )->second << " and " << dtok[j] << ")\n";
		  }
		
		// store actual (numeric, non-NaN) value
		double val;
		if ( Helper::str2dbl( dtok[j] , &val ) ) 
		  {		    
		    D[ var2slot[ expand_vname ] ][ id2slot[ id ] ] = val;
		  }
		else
		  {
		    // do we have a mapping object for this file/var?
		    std::map<std::string,std::map<std::string,double> >::const_iterator mm = mappings.find( tok[j] );
		    if ( mm != mappings.end() )
		      {			
			const std::map<std::string,double> & mp = mm->second;
			std::map<std::string,double>::const_iterator kk = mp.find( dtok[j] );
			if ( kk != mp.end() )
			  {			   
			    D[ var2slot[ expand_vname ] ][ id2slot[ id ] ] = kk->second ;
			  }
		      }
		    else
		      {
			// special case of sex M/F? -> 1/0
			// otherwise, **assume** T/F or Y/N -> 1/0 

			char f1 = std::toupper( dtok[j][0] );

			if ( Helper::iequals( "sex" , expand_vname ) || Helper::iequals( "gender" , expand_vname ) ) 
			  {
			    if ( f1 == 'F' ) D[ var2slot[ expand_vname ] ][ id2slot[ id ] ] = 0;
			    else if ( f1 == 'M' ) D[ var2slot[ expand_vname ] ][ id2slot[ id ] ] = 1;
			  }
			else // try default T/F and Y/N mappings
			  {

			    if ( f1 == 'F' || f1 == 'N' )
			      {
				// check doesn't match NA or NAN
				std::string nn = Helper::toupper( dtok[j] );
				if ( nn != "NA" && nn != "NAN" ) 
				  D[ var2slot[ expand_vname ] ][ id2slot[ id ] ] = 0;
			      }
			    else if ( f1 == 'T' || f1 == 'Y' )
			      {			   
				D[ var2slot[ expand_vname ] ][ id2slot[ id ] ] = 1;
			      }
			  }
		      }
		  }
	      }
	  
	  // next data row
	}
      
      // all done with this file      
      IN1.close();

      logger << "  ++ " << ff->first << ": "
	     << "read " << file_ids.size() << " indivs, "
	     << file_bvars.size() << " base vars --> "
	     << file_evars.size() << " expanded vars\n";
	           
      ++ff;
    }

  logger << "\n";
  
  //
  // ------- all inputs read now
  //

  const int ni = id2slot.size();
  const int nv = var2slot.size();
  X = Eigen::MatrixXd::Constant( ni , nv , std::numeric_limits<double>::quiet_NaN() );


  // slot2var / var2slot gives the order as read, which tracks factors, but we want to organize outputs
  // base the group and base var combined
  std::set<std::pair<std::string,std::string> > gbvars;
  std::map<std::string,std::string>::const_iterator bb = basevar.begin();
  while ( bb != basevar.end() )
    {
      gbvars.insert( std::make_pair( var2group[ bb->first ] , bb->second ) );
      ++bb;
    }
  
  // make final, ordered vars list
  vars.clear();
  std::map<std::string,int> final_var2slot;
  int cidx = 0;
  std::set<std::pair<std::string,std::string> >::const_iterator qq = gbvars.begin();
  while ( qq != gbvars.end() )
    {
      // use ordered list (slot2var)      
      std::map<int,std::string>::const_iterator ss = slot2var.begin();
      while ( ss != slot2var.end() )
	{
	  // matching group & basevar?
	  if ( var2group[ ss->second ] == qq->first && basevar[ ss->second ] == qq->second )
	    {
	      vars.push_back( ss->second );
	      final_var2slot[ ss->second ] = cidx++;
	    }
	  ++ss;
	}
      
      ++qq;
    }

  //
  // now add all data
  //

  ids.clear();
  ids.resize( id2slot.size() );
  
  std::map<std::string,int>::const_iterator ii = id2slot.begin();
  while ( ii != id2slot.end() )
    {
      // IDs
      ids[ ii->second ] = ii->first ;
      
      // data
      for (int j=0; j< vars.size(); j++)
	{
	  int final_slot = final_var2slot[ vars[j] ];
	  int slot = var2slot[ vars[j] ];
	  
	  if ( D[ slot  ].find( ii->second ) != D[ slot ].end() )
	    X( ii->second , final_slot ) = D[ slot ][ ii->second ];	    
	}
      
      ++ii;
    }

  //
  // -------- drop any completely empty cols
  //
  
  if ( ! retain_cols ) 
    drop_null_columns();
  
  //
  // -------- write to binary file
  //

  logger << "\n";
  bfile_t bf( bfile );
  bf.write( ids, vars , var2group, basevar, faclvl , X );
  logger << "  ...done\n";

  //
  // --------- show manifest (or dump)
  //
  
  if ( ! dump_file )
    manifest();
  else
    dump();
   
}


void gpa_t::read()
{
  bfile_t bf( bfile );
  bf.read( incvars , excvars,
	   inclvars, exclvars,
	   incnums, excnums,
	   incfacs, excfacs, incfaclvls, excfaclvls,
	   incgrps, excgrps,
	   xz_incvars, xz_incgrps,
	   y_incvars, y_incgrps,
	   s_incvars, s_incgrps,
	   &ids, &vars , &var2group, &basevar, &faclvl , &X );
  logger << "  read " << ids.size() << " individuals and "
	 << vars.size() << " variables from " << bfile << "\n";
}


void gpa_t::run()
{
  
  if ( ivs.size() == 0 || dvs.size() == 0 ) return;
  std::vector<std::string> xvars, yvars;
  for (int j=0; j<dvs.size(); j++) yvars.push_back( vars[ dvs[j] ] );
  for (int j=0; j<ivs.size(); j++) xvars.push_back( vars[ ivs[j] ] );
  
  // assume at this point that all dvs, ivs and cvs actually exist

  // run linmod_t separately for each X column,
  // but include all Z cols and all Y cols for correction

  const Eigen::MatrixXd & Y  = X( Eigen::all , dvs );
  const Eigen::MatrixXd & Z  = X( Eigen::all , cvs );
  
  // initiate (single X swapped in below)
  // correction done across all Y, but fixed within a single X (for now)
  linmod_t lm( Y , yvars, X( Eigen::all , ivs ), xvars, Z );

  // run
  linmod_results_t results = lm.run( nreps , show_progress , do_comparisons ? &comps : NULL );

  // outputs

  // manifest details
  const int ni = X.rows();
  const int nv = X.cols();
  std::set<std::string> allfacs;  
  // var -> fac -> lvl  
  std::map<std::string,std::map<std::string,std::string> >::const_iterator ii = faclvl.begin();
  while ( ii != faclvl.end() )
    {
      const std::map<std::string,std::string> & xx = ii->second;
      std::map<std::string,std::string>::const_iterator ff = xx.begin();
      while ( ff != xx.end() )
	{
	  allfacs.insert( ff->first );
	  ++ff;
	}      
      ++ii;
    }
  
  int count_p05 = 0 , count_padj05 = 0 , count_fdr05 = 0, count_all = 0;
  
  // iterate over X
  
  for (int j=0; j<ivs.size(); j++)
    {
      
      const std::string & xvar = vars[ ivs[j] ];
      writer.level( xvar , "X" );
      
      // over Y
      bool shown_y = false;
      for (int y=0; y<dvs.size(); y++)
	{
	  const std::string & var = vars[ dvs[y] ];
	  
	  const bool self = xvar == var; 
	  
	  shown_y = true; 
	  writer.level( var , "Y" );
	  
	  if ( ! self ) // keep as NA if self (Y == X)
	    {

	      writer.value( "B"  , results.beta[ xvar ][ var ] );
	      writer.value( "T"  , results.t[xvar][ var ] );
	      writer.value( "N" , (int)X.rows() );
	      
	      const bool valid_test = results.p[xvar][ var ] > -1 ; // -9 code means not valid
	      
	      if ( valid_test ) 
		{
		  writer.value( "P" , results.p[xvar][ var ] );
		  writer.value( "P_FDR"  , results.fdr_bh( xvar , var ) );
		}

	      if ( nreps != 0 )
		{
		  writer.value( "EMP" , results.emp[xvar][ var ] );
		  writer.value( "EMPADJ" , results.emp_corrected[xvar][ var ] );
		  if ( results.emp[xvar][ var ] < 0.05 ) count_p05++;
		  if ( results.emp_corrected[xvar][ var ] < 0.05 ) count_padj05++;
		  if ( valid_test && results.fdr_bh( xvar , var ) < 0.05 ) count_fdr05++;
		}
	      else // else summaary based on asymptotic & FDR  
		{
		  if ( valid_test )
		    {
		      if ( results.p[xvar][ var ] < 0.05 ) count_p05++;
		      if ( results.fdr_bh( xvar , var ) < 0.05 ) count_fdr05++;
		    }
		}
	      
	      if ( valid_test )
		{
		  if ( adj_fdr_by )
		    writer.value( "P_FDR_BY"  , results.fdr_by( xvar , var ) );
		  if ( adj_bonf )
		    writer.value( "P_BONF"  , results.bonf( xvar , var ) );
		  if ( adj_holm )
		    writer.value( "P_HOLM"  , results.holm( xvar , var ) );
		}

	      count_all++;
	    }
	  
	  // manifest details
	  writer.value( "GROUP" ,  var2group[ var ] );
	  writer.value( "BASE"  , basevar[ var ] );
	  
	  std::string x;
	  std::set<std::string>::const_iterator gg = allfacs.begin();
	  while ( gg != allfacs.end() )
	    {
	      std::map<std::string,std::map<std::string,std::string> >::const_iterator kk = faclvl.find( var );
	      if ( kk != faclvl.end() )
		{
		  if ( kk->second.find( *gg ) != kk->second.end() )
		    {
		      if ( x != "" ) x += ";";
		      x += *gg + "=" + kk->second.find( *gg )->second;
		    }
		}
	      ++gg;
	    }
	  writer.value( "STRAT" , x );
	  
	  // optional X-variable details too?
	  if ( show_xfacs )
	    {
	      writer.value( "XGROUP" ,  var2group[ xvar ] );
	      writer.value( "XBASE"  , basevar[ xvar ] );
	      std::string x;
	      std::set<std::string>::const_iterator gg = allfacs.begin();
	      while ( gg != allfacs.end() )
		{		      
		  std::map<std::string,std::map<std::string,std::string> >::const_iterator kk = faclvl.find( xvar );
		  if ( kk != faclvl.end() )
		    {
		      if ( kk->second.find( *gg ) != kk->second.end() )
			{
			  if ( x != "" ) x += ";";
			  x += *gg + "=" + kk->second.find( *gg )->second;
			}
		    }		      
		  ++gg;
		}	  
	      writer.value( "XSTRAT" , x );
	    }
	}
      if ( shown_y )
	writer.unlevel( "Y" );
      
    }
  
  writer.unlevel( "X" );
  
  logger << "  " << count_p05 << " (prop = " << count_p05 / (double)count_all << ") "
	 << "significant at nominal p < 0.05\n";
  
  logger << "  " << count_fdr05 << " (prop = " << count_fdr05 / (double)count_all << ") "
	 << "significant at FDR p < 0.05\n";
  
  if ( nreps != 0 ) 
    logger << "  " << count_padj05 << " (prop = " << count_padj05 / (double)count_all << ") "
	   << "significant after empirical family-wise type I error control p < 0.05\n";


  // optional comparison test
  if ( do_comparisons )
    {

      for (int j=0; j<ivs.size(); j++)
	{	  
	  const std::string & xvar = vars[ ivs[j] ];
	  writer.level( xvar , "X" );

	  linmod_comp_t * comp = comps.comp( xvar );
    	  
	  comp->evaluate();
	  writer.value( "MATCH_OBS" , comp->stat() );
	  writer.value( "MATCH_EXP" , comp->expected() );
	  writer.value( "Z" , comp->zstat() );
	  writer.value( "P" , comp->pvalue() );
	  writer.value( "P2" , comp->pvalue2() );

	}
      
      // iterate over all X values
      writer.unlevel( "X" );
      
    }
  
  
}




void gpa_t::run1X() // correction within X
{
  
  if ( ivs.size() == 0 || dvs.size() == 0 ) return;
  std::vector<std::string> xvars, yvars;
  for (int j=0; j<dvs.size(); j++) yvars.push_back( vars[ dvs[j] ] );
  for (int j=0; j<ivs.size(); j++) xvars.push_back( vars[ ivs[j] ] );
 
  
  const Eigen::MatrixXd & Y  = X( Eigen::all , dvs );
  const Eigen::MatrixXd & Z  = X( Eigen::all , cvs );
  
  // initiate (single X swapped in below)
  
  linmod_t lm( X( Eigen::all , dvs ) , yvars, // Y
	       X( Eigen::all , ivs ), xvars,  // X
	       X( Eigen::all , cvs ) );       // Z

  // manifest details
  const int ni = X.rows();
  const int nv = X.cols();

  std::set<std::string> allfacs;  
  // var -> fac -> lvl  
  std::map<std::string,std::map<std::string,std::string> >::const_iterator ii = faclvl.begin();
  while ( ii != faclvl.end() )
    {
      const std::map<std::string,std::string> & xx = ii->second;
      std::map<std::string,std::string>::const_iterator ff = xx.begin();
      while ( ff != xx.end() )
	{
	  allfacs.insert( ff->first );
	  ++ff;
	}      
      ++ii;
    }


  int count_p05 = 0 , count_padj05 = 0 , count_fdr05 = 0, count_all = 0;

  // iterate over X

  for (int j=0; j<ivs.size(); j++)
    {
      
      const std::string & xvar = vars[ ivs[j] ];
      writer.level( xvar , "X" );

      // swap in X
      lm.set_IV( X.col( ivs[j] ) , xvar );
      
      // run
      linmod_results_t results = lm.run( nreps , show_progress , do_comparisons ? &comps : NULL );

      // output
      bool shown_y = false;
      for (int y=0; y<dvs.size(); y++)
	{
	  const std::string & var = vars[ dvs[y] ];
	  
	  const bool self = xvar == var;
	  shown_y = true; 
	  writer.level( var , "Y" );
	  
	  if ( ! self )
	    {
	      writer.value( "B"  , results.beta[ xvar ][ var ] );
	      writer.value( "T"  , results.t[ xvar ][ var ] );
	      writer.value( "N" , (int)X.rows() );
	      
	      const bool valid_test = results.p[xvar][ var ] > -1; // valid p-value?
	      
	      if ( valid_test )
		{
		  writer.value( "P" , results.p[xvar][ var ] );
		  writer.value( "P_FDR"  , results.fdr_bh( xvar , var ) );
		}

	      if ( nreps != 0 )
		{
		  writer.value( "EMP" , results.emp[xvar][ var ] );
		  writer.value( "EMPADJ" , results.emp_corrected[xvar][ var ] );

		  if ( results.emp[xvar][ var ] < 0.05 ) count_p05++;
		  if ( results.emp_corrected[xvar][ var ] < 0.05 ) count_padj05++;
		  if ( valid_test && results.fdr_bh( xvar , var ) < 0.05 ) count_fdr05++;
		}
	      else // else summaary based on asymptotic & FDR  
		{
		  if ( valid_test )
		    {
		      if ( results.p[xvar][ var ] < 0.05 ) count_p05++;
		      if ( results.fdr_bh( xvar , var ) < 0.05 ) count_fdr05++;
		    }	   
		}   	      
	      
	      if ( valid_test )
		{
		  if ( adj_fdr_by)
		    writer.value( "P_FDR_BY"  , results.fdr_by( xvar , var ) );
		  if ( adj_bonf )
		    writer.value( "P_BONF"  , results.bonf( xvar , var ) );
		  if ( adj_holm )
		    writer.value( "P_HOLM"  , results.holm( xvar , var ) );
		}

	      count_all++;
	      
	    }

	  
	  // manifest details
	  writer.value( "GROUP" ,  var2group[ var ] );
	  writer.value( "BASE"  , basevar[ var ] );
	  
	  std::string x;
	  std::set<std::string>::const_iterator gg = allfacs.begin();
	  while ( gg != allfacs.end() )
	    {
	      std::map<std::string,std::map<std::string,std::string> >::const_iterator kk = faclvl.find( var );
	      if ( kk != faclvl.end() )
		{
		  if ( kk->second.find( *gg ) != kk->second.end() )
		    {			  
		      if ( x != "" ) x += ";";
		      x += *gg + "=" + kk->second.find( *gg )->second;
		    }
		}		  
	      ++gg;
	    }
	  writer.value( "STRAT" , x );
	  
	  // optional X var manifest
	  if ( show_xfacs )
	    {
	      writer.value( "XGROUP" ,  var2group[ xvar ] );
	      writer.value( "XBASE"  , basevar[ xvar ] );
	      
	      std::string x;
	      std::set<std::string>::const_iterator gg = allfacs.begin();
	      while ( gg != allfacs.end() )
		{	      
		  std::map<std::string,std::map<std::string,std::string> >::const_iterator kk = faclvl.find( xvar );
		  if ( kk != faclvl.end() )
		    {
		      if ( kk->second.find( *gg ) != kk->second.end() )
			{
			  if ( x != "" ) x += ";";
			  x += *gg + "=" + kk->second.find( *gg )->second;
			}
		    }
		  ++gg;
		}
	      writer.value( "XSTRAT" , x );
	    }	  
	}
      if ( shown_y )
	writer.unlevel( "Y" );
      
    }  
  writer.unlevel( "X" );
  

  logger << "  " << count_p05 << " (prop = " << count_p05 / (double)count_all << ") "
	 << "significant at nominal p < 0.05\n";
  
  logger << "  " << count_fdr05 << " (prop = " << count_fdr05 / (double)count_all << ") "
	 << "significant at FDR p < 0.05\n";
  
  if ( nreps != 0 ) 
    logger << "  " << count_padj05 << " (prop = " << count_padj05 / (double)count_all << ") "
	   << "significant after empirical family-wise type I error control p < 0.05\n";

  // optional comparison test
  
  if ( do_comparisons )
    {
      for (int j=0; j<ivs.size(); j++)
	{	  
	  const std::string & xvar = vars[ ivs[j] ];
	  writer.level( xvar , "X" );

	  std::cout << "xvar = " << xvar << "\n";
	  linmod_comp_t * comp = comps.comp( xvar );
	  
	  std::cout << "comp = " << ( comp == NULL ) << "\n";

	  comp->evaluate();
	  writer.value( "MATCH_OBS" , comp->stat() );
	  writer.value( "MATCH_EXP" , comp->expected() );
	  writer.value( "Z" , comp->zstat() );
	  writer.value( "P" , comp->pvalue() );
	  writer.value( "P2" , comp->pvalue2() );
	}
      
      // iterate over all X values
      writer.unlevel( "X" );

    }
  
}



// ------------------------------------------------------------
//
// bfile reader/writer
//
// ------------------------------------------------------------

bool bfile_t::write( const std::vector<std::string> & ids ,
		     const std::vector<std::string> & vars ,
		     const std::map<std::string,std::string> & var2group,
		     const std::map<std::string,std::string> & basevar,		     
		     const std::map<std::string,std::map<std::string,std::string> > & faclvl ,
		     const Eigen::MatrixXd & X )
{
  ni = X.rows();
  nv = X.cols();

  if ( vars.size() != nv || ids.size() != ni )
    Helper::halt( "internal error in bfile_t::write()" );
  
  logger << "  writing binary data (" << ni
	 << " obs, " << nv
	 << " variables) to " << name << "\n";

  
  // compile all unique factors (write all)
  std::set<std::string> facs;
  
  // var -> fac -> lvl  
  std::map<std::string,std::map<std::string,std::string> >::const_iterator ii = faclvl.begin();
  while ( ii != faclvl.end() )
    {
      const std::map<std::string,std::string> & xx = ii->second;
      std::map<std::string,std::string>::const_iterator ff = xx.begin();
      while ( ff != xx.end() )
	{
	  facs.insert( ff->first );
	  ++ff;
	}      
      ++ii;
    }
  
  const int nf = facs.size();
  
  // ID, ne, nf { features , row-major }   
  std::ofstream OUT1( Helper::expand( name ).c_str() , std::ios::binary | std::ios::out );

  // N's
  bwrite( OUT1, ni ) ;
  bwrite( OUT1, nv ) ;
  bwrite( OUT1, nf ) ;
  
  // IDs
  for (int i=0; i<ni; i++) 
    bwrite( OUT1, ids[i] ) ;

  // headers: expanded vars
  for (int i=0; i<nv; i++) 
    bwrite( OUT1, vars[i] ) ;

  // var-groups
  for (int i=0; i<nv; i++) 
    bwrite( OUT1, var2group.find( vars[i] )->second ) ;
  
  // base-vars
  for (int i=0; i<nv; i++) 
    bwrite( OUT1, basevar.find( vars[i] )->second ) ;
  
  // for each fac: (by alpha-order)
  std::set<std::string>::const_iterator ff = facs.begin();
  while ( ff != facs.end() )
    {
      // name of factors
      bwrite( OUT1, *ff ) ;
      
      // and then value for each variable (or '.' if N/A)
      
      for (int i=0; i<nv; i++)
	{
	  
	  // defailt (-->N/A)
	  std::string lvl = ".";
	  
	  // check var is in the faclvl map 
	  std::map<std::string,std::map<std::string,std::string> >::const_iterator vv = faclvl.find( vars[i] );
	  
	  if ( vv != faclvl.end() )
	    {
	      // does this var have this factor associated?
	      std::map<std::string,std::string>::const_iterator ll = vv->second.find( *ff );
	      if ( ll != vv->second.end() )
		lvl = ll->second;
	    }
	  	  
	  // may be '.' still
	  bwrite( OUT1, lvl );

	}

      // next factor
      ++ff;
    }

  // now all numeric data  [ var x ind for easier subsetting of vars on reads ] 
  for (int j=0; j<nv; j++)
    for (int i=0; i<ni; i++)    
      bwrite( OUT1, X(i,j) );
  
  // done
  OUT1.close();

  return true;
}

 
bool bfile_t::read( const std::set<std::string> & incvars ,
		    const std::set<std::string> & excvars ,
		    const std::set<std::string> & inclvars,
		    const std::set<std::string> & exclvars,
		    const std::vector<std::pair<int,int> > & incnums,
		    const std::vector<std::pair<int,int> > & excnums, 
		    const std::set<std::string> & incfacs,
		    const std::set<std::string> & excfacs,
		    const std::map<std::string,std::set<std::string> > & incfaclvls,
		    const std::map<std::string,std::set<std::string> > & excfaclvls, 		    
		    const std::set<std::string> & incgrps,
		    const std::set<std::string> & excgrps,		    
		    const std::set<std::string> & xz_incvars,
		    const std::set<std::string> & xz_incgrps,		    
		    const std::set<std::string> & y_incvars,
		    const std::set<std::string> & y_incgrps,		    
		    const std::set<std::string> & s_incvars,
		    const std::set<std::string> & s_incgrps,		    
		    std::vector<std::string> * ids ,
		    std::vector<std::string> * vars ,
		    std::map<std::string,std::string> * var2group ,
		    std::map<std::string,std::string> * basevar,		    
		    std::map<std::string,std::map<std::string,std::string> > * faclvl ,
		    Eigen::MatrixXd * X )
{


  if ( ! Helper::fileExists( Helper::expand( name ) ) )
    Helper::halt( "could not open " + name );
      
  std::ifstream IN1( Helper::expand( name ).c_str() , std::ios::binary | std::ios::in );
  
  // size of data (in file) 
  const int ni = bread_int( IN1 );
  const int nv = bread_int( IN1 );
  const int nf = bread_int( IN1 );

  //
  // initial IDs and full (file) vars
  //
  
  // IDs
  ids->resize( ni );
  for (int i=0; i<ni; i++)
    (*ids)[i] = bread_str( IN1 );
  
    
  // variables read in (which may be subsetted)
  std::vector<std::string> all_vars( nv );
  std::vector<std::string> all_groups( nv );
  std::vector<std::string> all_basevars( nv );
  std::map<std::string,std::map<std::string,std::string > > all_faclvl;
  
  // variables (expanded)
  for (int j=0; j<nv; j++)
    all_vars[j] = bread_str( IN1 );
  
  // group variable
  for (int j=0; j<nv; j++)
    all_groups[j] = bread_str( IN1 );

  // base variable
  for (int j=0; j<nv; j++)
    all_basevars[j] = bread_str( IN1 );

  // fac/lvls  
  std::vector<std::string> facs( nf );
  for (int k=0; k<nf; k++)
    {
      facs[k] = bread_str( IN1 );
      
      for (int j=0; j<nv; j++)
	{
	  const std::string lvl = bread_str( IN1 );
	  if ( lvl != "." )
	    all_faclvl[ all_vars[j] ][ facs[k] ] = lvl;
	}
    }

  
  //
  // restrict inputs?
  //
  
  std::vector<bool> readvar( nv , true );

  // explicitly track exclusions, to handle xyz_includes
  std::vector<bool> excluded( nv ,  false );
  
  if ( incvars.size() || excvars.size()
       || inclvars.size() || exclvars.size()
       || incnums.size() || excnums.size()
       || incgrps.size() || excgrps.size() 
       || incfacs.size() || excfacs.size()
       || incfaclvls.size() || excfaclvls.size()
       || xz_incvars.size() || xz_incgrps.size()
       || s_incvars.size() || s_incgrps.size() 
       || y_incvars.size() || y_incgrps.size() )
    {
      
      // includes: within a type (e.g. vars) apply OR logic
      //           across types apply AND logic (for subsetting)
      
      // excludes: apply OR logic across all (i.e. 

      // XYZ includes: these are added (irrespective of other include AND-logic)
      //  unless they are explicitly excluded; do this on a second round: always
      //  added, unless we flag that we've explicitly excluded this (i.e. versus
      //  a failure to include)

      // Y includes handled separately from XZ includes:
      //   i.e. not being in an XZ include does not imply we don't want to read (as may be a Y)
      //        but if we have Y includes, then not being part of that implies we don't want to
      //        read (unless in XZ)
      //   i.e. this arises as the Y variables are otherwise just what is left over
      //        after reading all and removing X/Z
      
      std::set<int> incnums_set;
      for (int i=0; i<incnums.size(); i++  )
	{
	  int s1 = incnums[i].first < incnums[i].second ? incnums[i].first : incnums[i].second ;
          int s2 = incnums[i].first < incnums[i].second ? incnums[i].second : incnums[i].first ;
	  for (int j=s1; j<=s2; j++)
	    incnums_set.insert( j );
	}

      std::set<int> excnums_set;      
      for (int i=0; i<excnums.size(); i++  )
	{
	  int s1 = excnums[i].first < excnums[i].second ? excnums[i].first : excnums[i].second ;
          int s2 = excnums[i].first < excnums[i].second ? excnums[i].second : excnums[i].first ;
	  for (int j=s1; j<=s2; j++)
	    excnums_set.insert( j );
	}
      
      const bool has_incvars = incvars.size();
      const bool has_inclvars = inclvars.size();
      const bool has_incnums = incnums.size();
      const bool has_incgrps = incgrps.size();
      const bool has_incfacs   = incfacs.size();
      const bool has_incfaclvls = incfaclvls.size();

      const bool has_xz_vars = xz_incvars.size();
      const bool has_xz_grps = xz_incgrps.size();

      const bool has_y_vars = y_incvars.size();
      const bool has_y_grps = y_incgrps.size();

      const bool has_s_vars = s_incvars.size();
      const bool has_s_grps = s_incgrps.size();
      
      const bool has_excvars = excvars.size();
      const bool has_exclvars = exclvars.size();
      const bool has_excnums = excnums.size();
      const bool has_excgrps = excgrps.size();
      const bool has_excfacs = excfacs.size();
      const bool has_excfaclvls = excfaclvls.size();

      //
      // cannot both have Y or Yg as well as if we've selected 
      //
      
      //
      // iterate over each variable
      //
      
           
      for (int j=0; j<nv; j++)
	{

	  // exclusions first
	  if ( has_excnums )
	    {
	      if ( excnums_set.find( j ) != excnums_set.end() )
		{
		  readvar[j] = false;
		  excluded[j] = true;
		  continue;
		}
	    }
	  
	  // base-var exc?
	  if ( has_excvars ) 
	    {
	      if ( excvars.find( all_basevars[j] ) != excvars.end() )
		{
		  readvar[j] = false;
		  excluded[j] = true;
		  continue;
		}
	    }
      
	  // base-var exc?
	  if ( has_exclvars ) 
	    {
	      if ( exclvars.find( all_vars[j] ) != exclvars.end() )
		{
                  readvar[j] = false;
		  excluded[j] = true;
                  continue;
		}
	    }
	  
	  // var-group exc?
	  if ( has_excgrps ) 
	    {
	      if ( excgrps.find( all_groups[j] ) != excgrps.end() )
		{
		  readvar[j] = false;
		  excluded[j] = true;
		  continue;
		}
	    }
	  
	  // excfacs
	  if ( has_excfacs )
	    {

	      const std::map<std::string,std::string> & fl = all_faclvl[ all_vars[j] ];
	      
	      int c=0;
	      std::set<std::string>::const_iterator kk = excfacs.begin();
	      while ( kk != excfacs.end() )
		{
		  if ( fl.find( *kk ) != fl.end() ) ++c;
		  ++kk;
		}

	      // exclude if has /all/ these factors
	      if ( c == excfacs.size() && c == fl.size() )
		{
		  readvar[j] = false;
		  excluded[j] = true;
		  continue;
		}
	    }
	  
	  // excfaclvls      
	  if ( has_excfaclvls )
	    {
	      
	      // only include/exclude based on factor levels
	      // irrespective of whether the full fac-set for that
	      // variable matches;
	      //  i.e. if faclvls=CH/CZ|FZ,B=SIGMA|GAMMA
	      //       then only pull CH==CZ whether fac is just CH, or CH+B, or CH+F, etc
	      //  do includes first, then excludes
	      
	      // the above, has an AND logic between CH and B
	      // but an OR logic between CZ vs FZ, or SIGMA vs GAMMA
	  		  
	      const std::map<std::string,std::string> & fl = all_faclvl[ all_vars[j] ];				  
		  
	      bool consistent = true;
		      
	      // consider each factor;  okay for it not to exist, but if it does, then
	      // we must *NOT* match one of the listed factors
	      std::map<std::string,std::set<std::string> >::const_iterator ll = excfaclvls.begin();
	      while ( ll != excfaclvls.end() )
		{
		  // we have this factor...
		  if ( fl.find( ll->first ) != fl.end() )
		    {
		      // ... do we have an unacceptable level?
		      const std::string & lvl = fl.find( ll->first )->second;
		      
		      if ( ll->second.find( lvl ) != ll->second.end() ) // yes
			{
			  // then we are a no-go
			  consistent = false;
			  break;
			}
		    }
		  ++ll;
		}
	      
	      if ( ! consistent )
		{
		  readvar[j] = false;		    
		  excluded[j] = true;
		  continue;
		}	      
	    }


	  //
	  // now, if still here, process any inclusions
	  //  -- i.e. set to missing if *not* a match for any
	  //          one category, which implies an OR match
	  //          across types for inclusions
	  //
	  
          if ( has_incnums )
            {
              if ( incnums_set.find( j ) == incnums_set.end() )
                {
                  readvar[j] = false;
                  continue;
                }
            }

	  // base-var inc?
	  if ( has_incvars ) 
	    {
	      if ( incvars.find( all_basevars[j] ) == incvars.end() )
		{
		  readvar[j] = false;
		  continue;
		}
	    }
      
	  // l-var inc?
	  if ( has_inclvars ) 
	    {
	      if ( inclvars.find( all_vars[j] ) == inclvars.end() )
		{
		  readvar[j] = false;
		  continue;		  
		}
	    }
	  
	  // var-group inc?
	  if ( has_incgrps )
	    {
	      if ( incgrps.find( all_groups[j] ) == incgrps.end() )
		{
                  readvar[j] = false;
                  continue;
		}
	    }
	  	  
	  // infacs
	  if ( has_incfacs )
	    {
	      // exclude if does not contain all these factors
	      
	      bool consistent = true;
	      
	      // infacs: include if has all these factors
	      const std::map<std::string,std::string> & fl = all_faclvl[ all_vars[j] ];
	      
	      if ( fl.size() != incfacs.size() )
		consistent = false;
	      else
		{
		  std::set<std::string>::const_iterator kk = incfacs.begin();
		  while ( kk != incfacs.end() )
		    {
		      if ( fl.find( *kk ) == fl.end() )
			{
			  consistent = false;
			  break;
			}
		      ++kk;
		    }
		}

	      // matched?
	      if ( ! consistent )
		{
		  readvar[j] = false;
		  continue;
		}
	    }
	      
	  
	  // infaclvls 
	  
	  if ( has_incfaclvls )
	    {
	      
	      // only include/exclude based on factor levels
	      // irrespective of whether the full fac-set for that
	      // variable matches;
	      //  i.e. if faclvls=CH/CZ|FZ,B=SIGMA|GAMMA
	      //       then only pull CH==CZ whether fac is just CH, or CH+B, or CH+F, etc
	      //  do includes first, then excludes
	      
	      // the above, has an AND logic between CH and B
	      // but an OR logic between CZ vs FZ, or SIGMA vs GAMMA
	      
	      // if the variable doesn't have /any/ relevant factors, then we leave as is
	      
	      const std::map<std::string,std::string> & fl = all_faclvl[ all_vars[j] ];
	      
	      bool consistent = true; 
	      bool has_facs = false;
	      
	      // consider each factor;  okay for it not to exist, but if it does, then
	      // we must match one of the listed factors
	      std::map<std::string,std::set<std::string> >::const_iterator ll = incfaclvls.begin();
	      while ( ll != incfaclvls.end() )
		{
		  // we have this factor...
		  if ( fl.find( ll->first ) != fl.end() )
		    {
		      
		      has_facs = true;
		      
		      // ... do we have an acceptable level?
		      const std::string & lvl = fl.find( ll->first )->second;
		      
		      if ( ll->second.find( lvl ) == ll->second.end() ) // no
			{
			  // no... then we are a no-go
			  consistent = false;
			  break;
			}
		    }
		  ++ll;
		}
	      
	      if ( ! has_facs )
		{
		  readvar[j] = false;
		  continue;
		}

	      if ( ! consistent )
		{
		  readvar[j] = false;
                  continue;
		}
		  
	    }

	  //
	  // next variable to process for inc/excl status
	  //

	}

      
      //
      // now do xyz-includes as a second loop
      //

      if ( has_xz_vars || has_xz_grps || has_y_vars || has_y_grps | has_s_vars || has_s_grps )
	{
	  for (int j=0; j<nv; j++)
	    {
	      // if already *explicitly* excluded, do not try to add in
	      if( ! excluded[j] ) 
		{
		  
		  if ( has_xz_vars && 			
		       xz_incvars.find( all_basevars[j] ) != xz_incvars.end() )
		    readvar[j] = true;		      
		  else if ( has_xz_grps && 			
			    xz_incgrps.find( all_groups[j] ) != xz_incgrps.end() )
		    readvar[j] = true;
		  else if ( has_s_vars &&
			    s_incvars.find( all_basevars[j] ) != s_incvars.end() )
		    readvar[j] = true;
		  else if ( has_s_grps &&
			    s_incgrps.find( all_groups[j] ) != s_incgrps.end() )
		    readvar[j] = true;
		  else
		    {
		      // if not XZ included, *and* we have a Y include, then
		      // if we don't find it, set to exclude
		      
		      if ( has_y_vars || has_y_grps )
			{
			  bool found = false;
			  if ( has_y_vars &&
			       y_incvars.find( all_basevars[j] ) != y_incvars.end() )
			    found = true;
			  
			  if ( has_y_grps &&
			       y_incgrps.find( all_groups[j] ) != y_incgrps.end() )
			    found = true;
			  
			  if ( found ) readvar[j] = true;
			  else readvar[j] = false;
			}
		      
		    }		
		}
	    }
	}
      
    }


  //
  // copy subset into file return values
  //
  
  vars->clear();
  basevar->clear();
  var2group->clear();
  
  for (int i=0; i<nv; i++)
    if ( readvar[i] )
      {
	vars->push_back( all_vars[i] );
	(*var2group)[ all_vars[i] ] = all_groups[i];
	(*basevar)[ all_vars[i] ] = all_basevars[i];	
      }
  
  const int nv2 = vars->size();;
  
  logger << "  reading " << nv2 << " of " << nv << " vars on " << ni << " indivs\n";
    
  // copy selected factors over
  faclvl->clear();
  // for (int k=0; k<nf; k++)
  //   {
  //     for (int j=0; j<nv; j++)
  // 	if ( readvar[j] )
  // 	  (*faclvl)[ all_vars[j] ][ facs[k] ] = all_faclvl[ all_vars[j] ][ facs[k] ];
  //   }
  
  for (int j=0; j<nv; j++)
    if ( readvar[j] )
      (*faclvl)[ all_vars[j] ] = all_faclvl[ all_vars[j] ];
  
    
  
  // data [ var x ind ] --> X[ ind x var ] 

  *X = Eigen::MatrixXd::Zero( ni , nv2 );
  
  int cidx = 0;
  for (int j=0; j<nv; j++)
    {
      if ( readvar[j] ) 
	{
	  for (int i=0; i<ni; i++)
	    (*X)(i,cidx) = bread_dbl( IN1 );
	  ++cidx;
	}
      else
	bskip_dbl( IN1 , ni ); // skip variable
    }
  
  // all done
  IN1.close();
  return true;
}


void gpa_t::dump()
{
  const int ni = X.rows();
  const int nv = X.cols();

  // header
  std::cout << "ID";
  for (int j=0; j<nv; j++)
    std::cout << "\t" << vars[j];
  std::cout << "\n";

  // data
  for (int i=0; i<ni; i++)
    {
      std::cout << ids[i] ;
      for (int j=0; j<nv; j++)
	std::cout << "\t" << X(i,j);
      std::cout << "\n";
    }
  
}


// make strata (for stats())
std::vector<std::map<std::string,std::string> > gpa_t::make_strata( const std::vector<std::string> & s ,
								    std::vector<std::vector<int>> * x )
{

  // given a set of variables in s[], find all N>0 combinations of values (i.e. we're assuming
  // categorical/integer values passed) and populate two vectors (each element = one combination)
  //   x - the rows included
  //   retval - the labels 
  
  const int ns = s.size();
  const int ni = X.rows();
  
  // helper for lookups
  std::map<std::string,int> v2n;
  const int nv = X.cols();
  for (int j=0; j<nv; j++)
    v2n[ vars[j] ] = j;
  
  // get slots for each var
  std::vector<int> cols( ns );  
  for (int i=0; i<s.size(); i++)
    {
      if ( v2n.find( s[i] ) == v2n.end() )
	Helper::halt( "could not find " + s[i] );
      cols[i] = v2n[ s[i] ];
    }

  // get all unique values of each stratifying variable
  // as well as who has that combo
  std::vector<std::map<std::string,std::set<int> > > sxs(ns);  
  std::vector<std::vector<std::string> > lxs (ns);  
  std::vector<int> nl( ns );
  
  for (int j=0; j<cols.size(); j++)
    {
      const int slot = cols[j];
      std::map<std::string,std::set<int> > sx;
      for (int i=0; i<ni; i++)
	sx[ Helper::dbl2str( X(i,slot) ) ].insert( i );
      sxs[j] = sx ;
      
      // track rows
      std::map<std::string,std::set<int> >::const_iterator ii = sx.begin();
      while ( ii != sx.end() ) { lxs[j].push_back( ii->first ) ; ii++ ; } 
      nl[j] = sx.size();
    }

  // now sxs tracks for each stratifying var, what the levels are, and for each
  // level, who the obs w/ that value are
  //  sxs[ var ][ level ][ (inds...) ]
  
  // want to populate this: std::vector<std::vector<int>> for each
  // unique combo of levels across vars

  // use this to go over all combos of levels, i.e. 0 to nl[j]-1 
  std::vector<int> n1( ns , 0 );

  // return vals
  x->clear();
  std::vector<std::map<std::string,std::string> > res;
  
  do {
    
    // populate levels
    std::vector<std::string> labels( ns );
    std::map<std::string,std::string> f2l;    
    std::map<int,int> includes;
    
    for (int i=0; i<ns; i++)
      {
	const std::map<std::string,std::set<int> > & sx = sxs[ i ];	
	labels[i] = lxs[ i ][ n1[i] ] ;
		
	// get people w/ this level (should always be good)
	const std::set<int> & p = sx.find( labels[i] )->second;
	//std::cout << "  " << s[i] << " = " << labels[i] << " n = " << p.size() << "\n";
	f2l[ s[i]  ] = labels[i] ;
	
	std::set<int>::const_iterator pp = p.begin();
	while ( pp != p.end() )
	  {
	    includes[ *pp ]++;
	    ++pp;
	  }	
	
      }

    // now find people who match all vars
    std::vector<int> pinc;
    std::map<int,int>::const_iterator ii = includes.begin();
    while ( ii != includes.end() )
      {
	if ( ii->second == ns )
	  pinc.push_back( ii->first );
	++ii;
      }
        
    // add people list to return, if non-missing
    const int nsi = pinc.size();
    if ( nsi > 0 )
      {
	x->push_back( pinc );
	res.push_back( f2l ) ;
      }
    
  } while ( next_combo( n1 , nl ) );
    
  // return strata defns
  return res;
    
}

// helper: mixed radix counter
bool gpa_t::next_combo(std::vector<int>& a, const std::vector<int>& L) {
  const size_t n = a.size();
  for (size_t i = 0; i < n; ++i) {          // i = 0 is least significant
    if (++a[i] < L[i]) return true;         // carry stops here
    a[i] = 0;                               // wrap and carry to next digit
  }
  return false;                             // overflow: we’re done
}


void gpa_t::stats( const std::vector<std::string> * s )
{

  // total
  const int ni = X.rows();
  const int nv = X.cols();

  // by strata?
  const bool by_strata = s != NULL && s->size() != 0 ; 
  
  // which individuals for which strata?
  std::vector<std::vector<int>> sindiv;

  // strata labels
  std::vector<std::map<std::string,std::string> > slabels;

  // strata variable lookup
  std::set<std::string> sset;
  
  // make strata
  if ( by_strata )
    {
      slabels = make_strata( *s , &sindiv );
      for (int i=0; i<s->size(); i++) sset.insert( (*s)[i] );
    }

  // number of strata
  const int ns = by_strata ? slabels.size() : 0 ; 

  // report for each variable
  for (int j=0; j<nv; j++)
    {
      
      // except skip any stratifing variables
      if ( sset.find( vars[j] ) != sset.end() )
	continue;
      
      writer.level( vars[j] , globals::var_strat );
      
      // get variable
      const Eigen::VectorXd & col = X.col(j);

      // calc/report for unstratified version
      if ( ! by_strata )
	{
	  // actual number of non-missing obs
	  Eigen::VectorXd x = eigen_ops::removeNaNs( col );
	  int na = x.size();
	  writer.value( "NOBS" , na );
	  
	  if ( na > 0 ) {	    
	    const double mean = x.mean();
	    writer.value( "MEAN" , mean );
	    if ( na > 1 ) { 
	      const double sd = eigen_ops::sdev( x );
	      writer.value( "SD" , sd );
	    }
	  }
	}

      // else, stratified report

      if ( by_strata )
	{
	  // stratum count
	  int si = 0; 

	  while ( 1 )
	    {
	      // get indivs for this stratum
	      const std::vector<int> & idx = sindiv[si];
	      Eigen::Map<const Eigen::VectorXi> ei(idx.data(), idx.size());
	      const Eigen::VectorXd & col2 = col(ei);
	      int ni2 = col2.size();
	      
	      writer.level( Helper::kv_print( slabels[si] ) , "STRATUM" );	      
	      
	      // actual number of non-missing obs
	      Eigen::VectorXd x = eigen_ops::removeNaNs( col2 );
	      int nobs = x.size();

	      writer.value( "NOBS" , nobs );
	      if ( nobs > 0 )
		{
		  const double mean = x.mean();
		  writer.value( "MEAN" , mean );
		  if ( nobs > 1 )
		    {
		      const double sd = eigen_ops::sdev( x ) ; 
		      writer.value( "SD" , sd );
		    }
		}
	      
	      // all done?
	      ++si;
	      if ( si == ns ) break;
	    }
	  writer.unlevel( "STRATUM" );
	}

    } // next variable
  
  writer.unlevel( globals::var_strat );
  
}


void gpa_t::summarize()
{

  const int ni = X.rows();
  const int nv = X.cols();


  logger << "\n"
	 << "  ------------------------------------------------------------\n"
	 << "  Summary for (the selected subset of) " << bfile << "\n\n";

  logger << "  # individuals (rows)  = " << ni << "\n"
	 << "  # variables (columns) = " << nv << "\n";

  // get all factors left in 
  std::set<std::string> allfacs;  
  // var -> fac -> lvl  
  std::map<std::string,std::map<std::string,std::string> >::const_iterator aa = faclvl.begin();
  while ( aa != faclvl.end() )
    {
      const std::map<std::string,std::string> & xx = aa->second;
      std::map<std::string,std::string>::const_iterator ff = xx.begin();
      while ( ff != xx.end() )
	{
	  allfacs.insert( ff->first );
	  ++ff;
	}      
      ++aa;
    }
  
  // grp -> # of vars
  std::map<std::string,int> grp2vars;
  
  // grp -> base-vars
  std::map<std::string,std::set<std::string> > grp2basevars;
  
  // grp -> # fac -> lvls -> # vars
  std::map<std::string,std::map<std::string,std::map<std::string,int> > > grp2fac2lvl2vars;
  
  // basevar -> # of vars
  std::map<std::string,int> basevar2vars;
  

  // vars
  for (int j=0; j<nv; j++)
    {

      const std::string g = var2group[ vars[j] ];

      grp2vars[ g ]++;
      grp2basevars[ g ].insert( basevar[ vars[j] ] );
      basevar2vars[ basevar[ vars[j] ] ]++;

      std::set<std::string>::const_iterator gg = allfacs.begin();
      while ( gg != allfacs.end() )
	{
	  std::map<std::string,std::map<std::string,std::string> >::const_iterator kk = faclvl.find( vars[j] );
	  if ( kk != faclvl.end() )
	    {
	      if ( kk->second.find( *gg ) != kk->second.end() )
		grp2fac2lvl2vars[ g ][ *gg ][ kk->second.find( *gg )->second ]++;
	    }
	  ++gg;
	}

    }  


  //
  // summary by basevars
  // 
  
  logger << "\n"
	 << "  ------------------------------------------------------------\n"
	 << "  " << basevar2vars.size() << " base variable(s):\n\n";
  
  std::map<std::string,int>::const_iterator bb =  basevar2vars.begin();
  while ( bb != basevar2vars.end() )
    {
      logger << "    " << bb->first << " --> " << bb->second << " expanded variable(s)\n";
      ++bb;
    }

  
  //
  // summary by groups
  //
  
  logger << "\n"
	 << "  ------------------------------------------------------------\n"
         << "  " << grp2vars.size() << " variable groups:\n\n";

  std::map<std::string,int>::const_iterator ii =  grp2vars.begin();
  while	( ii != grp2vars.end() )
    {

      const std::set<std::string> & ref = grp2basevars[ ii->first ];

      logger << "    " << ii->first << ":\n"
	     << "         " << ref.size() << " base variable(s) (";
      
      int c = 0;
      std::set<std::string>::const_iterator rr = ref.begin();
      while ( rr != ref.end() ) {
	logger << " " << *rr;
	++rr;

	// enough?
	++c;	
	if ( c == 6 && ref.size() > 5 )
	  {
	    logger << " ...";
	    break;
	  }
      }
      logger << " )\n";

      logger << "         " << ii->second << " expanded variable(s)\n";


      // grp -> # faclvls
      const std::map<std::string,std::map<std::string,int> > & ref2 = grp2fac2lvl2vars[ ii->first ];

      if ( ref2.size() == 0 )
	logger << "         " << "no stratifying factors\n";
      else
	logger << "         " << ref2.size() << " unique factors:\n";

      std::map<std::string,std::map<std::string,int> >::const_iterator ff = ref2.begin();
      while ( ff != ref2.end() )
	{
	  const std::map<std::string,int> & ref3 = ff->second;
	  logger << "           " << ff->first << " -> " << ref3.size() << " level(s):\n";
	  std::map<std::string,int>::const_iterator rr = ref3.begin();

	  int c = 0;
	  while ( rr != ref3.end() )
	    {
	      logger << "             " << ff->first << " = " << rr->first << " --> " << rr->second << " expanded variable(s)\n";	      
	      ++rr;

	      // enough?
	      ++c;
	      if ( c == 6 && ref3.size() > 5 )
		{
		  logger << "             ...\n";
		  break;
		}
	    }
	  ++ff;
	}
      logger << "\n";
      
      ++ii;
    }

  logger << "  ------------------------------------------------------------\n\n";
  
}


void gpa_t::manifest()
{

  const int ni = X.rows();
  const int nv = X.cols();
  
  //
  // verbose listing to output
  //

  // first, get all factors
  
  std::set<std::string> allfacs;
  
  // var -> fac -> lvl  
  std::map<std::string,std::map<std::string,std::string> >::const_iterator ii = faclvl.begin();
  while ( ii != faclvl.end() )
    {
      const std::map<std::string,std::string> & xx = ii->second;
      std::map<std::string,std::string>::const_iterator ff = xx.begin();
      while ( ff != xx.end() )
	{
	  //std::cout << " adding " << ii->first << " " << ff->first << "\n";
	  allfacs.insert( ff->first );
	  ++ff;
	}      
      ++ii;
    }
  
  std::cout << "NV\t"
	    << "VAR\t"
	    << "NI\t"
	    << "GRP\t"
	    << "BASE";

  std::set<std::string>::const_iterator gg = allfacs.begin();
  while ( gg != allfacs.end() )
    {
      std::cout << "\t" << *gg ;
      ++gg;
    }  
  std::cout << "\n";

  // all
  std::cout << "0\t"
            << "ID\t"
            << ni << "\t"
	    << ".\t"
            << "ID";
  
  gg = allfacs.begin();
  while ( gg != allfacs.end() )
    {
      std::cout << "\t." ;
      ++gg;
    }
  std::cout << "\n";

  
  // vars
  for (int j=0; j<nv; j++)
    {
      int na = 0;
      const Eigen::VectorXd & col = X.col(j);
      for (int i=0; i<ni; i++) if ( std::isnan( col[i] ) ) ++na;
      
      std::cout << j+1 << "\t"
		<< vars[j] << "\t"
		<< ni - na << "\t"
		<< var2group[ vars[j] ] << "\t"
		<< basevar[ vars[j] ] ;

      //std::cout << " | \n-->";
		
      std::set<std::string>::const_iterator gg = allfacs.begin();
      while ( gg != allfacs.end() )
	{
	  std::string x = ".";
	  std::map<std::string,std::map<std::string,std::string> >::const_iterator kk = faclvl.find( vars[j] );
	  if ( kk != faclvl.end() )
	    {
	      if ( kk->second.find( *gg ) != kk->second.end() )
		{
		  x = kk->second.find( *gg )->second;
		  //std::cout << " *gg " << *gg << " x = " << x << "\n";		  
		}
	    }
	  std::cout << "\t" << x ;
	  ++gg;
	}
      std::cout << "\n";
      
    } // next var
  
}


// subset rows
void gpa_t::subset( const std::set<int> & rows , const std::map<int,bool> & cols )
{
  int ni = X.rows();
  
  const bool id_subsetting = rows.size() != 0 && rows.size() != ni;
  const bool col_subsetting = cols.size() != 0;
  
  // nothing to do
  if ( ! ( id_subsetting || col_subsetting ) ) return;
  

  std::vector<bool> included( ni , true );
  
  // first, select any named IDs (if any specified, else retain all) 
  if ( id_subsetting )
    {
      for (int i=0; i<ni; i++)
	if ( rows.find( i ) == rows.end() )
	  included[i] = false;
    }
  
  // second, find other rows to select based on matching col conditions
  // if flip == T then reverse selection
  // if multiple cols, implies must match on all
  std::stringstream ss;
  if ( col_subsetting )
    {            
      std::map<int,bool>::const_iterator cc = cols.begin();
      while ( cc != cols.end() )
	{
	  const bool pos_match = cc->second;      
	  const Eigen::VectorXd & v = X.col( cc->first );

	  if ( cc != cols.begin() )
	    ss << ", ";
	  ss << vars[cc->first];
	  
	  for (int i=0;i<ni; i++)
	    {
	      // 0 or NaN is a 'negative' value for a subsetting col
	      const bool pos_value = ! ( std::isnan( v[i] ) || fabs( v[i] ) < 1e-4 );
	      if ( pos_match != pos_value ) included[i] = false;
	    }
	  // next term
	  ++cc;
	}
    }

  std::vector<int> r;
  for (int i=0;i<ni;i++)
    if ( included[i] ) r.push_back(i);

  // nothing to drop?
  if ( r.size() == ni ) return;

  // copy in case of aliasing etc
  Eigen::MatrixXd X1 = X( r, Eigen::all );
  X = X1;

  // also, updated IDs
  std::vector<std::string> ids2 = ids;
  ids.clear();
  for (int i=0;i<ni;i++)
    if ( included[i] ) ids.push_back( ids2[i] );

  // report
  const int nvsub = cols.size();
  
  std::map<int,bool>::const_iterator cc = cols.begin();
  
  logger << "  subsetted X to " << X.rows() << "/" << ni << " rows based on " ;
  if ( id_subsetting )
    {
      logger << "row exclusion";
      if ( col_subsetting ) logger << ", ";
    }
  if ( col_subsetting )
    logger << nvsub << " variable(s) (" << ss.str() << ")";
  logger << "\n";
}

// drop bad cols
void gpa_t::drop_null_columns()
{

  logger << "  checking for too much missing data ('retain-cols' to skip; 'verbose' to list)\n";

  // nothing to do
  if ( n_req == 0 && n_prop < 1e-6 )
    {
      logger << "  nothing to check (n-rep and n-prop set to 0)\n";
      return; 
    }

  logger << "  requiring at least n-req=" << n_req << " non-missing obs "
	 << "(or proportion n-prop=" << n_prop << ")\n";
    
  // uses n_req and n_prop
  
  std::vector<int> good_cols;

  const int ni = X.rows();
  const int nv = X.cols();

  // n_req cannot be larger than ni
  if ( n_req > ni ) n_req = ni;

  for (int j=0; j<nv; j++)
    {
      int na = 0;
      for (int i=0; i<ni; i++) if ( std::isnan( X(i,j) ) ) ++na;
      int ng = ni - na;
      if ( ng >= n_req && ng / (double)(ni) >= n_prop )
	good_cols.push_back( j );
      else if ( verbose ) 
	logger << "  *** dropping " << vars[j] 
	       << " due to missing values: " << ng << "/" << ni << " = " << ng / (double)(ni) 
	       << " good values, given n-req=" << n_req << " and n-prop=" << n_prop << " required\n";
    }

  // any to remove? 
  if ( good_cols.size() < nv )
    {
      
      Eigen::MatrixXd X1 = X( Eigen::all , good_cols );
      X = X1;
      
      // now need to update vars, var2group, basevar
      // okay to keep faclvl as is

      // var labels
      std::vector<std::string> vars2 = vars;
      std::map<std::string,std::map<std::string,std::string> > faclvl2 = faclvl;
      std::map<std::string,std::string> basevar2 = basevar;
      std::map<std::string,std::string> var2group2 = var2group;;

      std::set<int> dvs2 = Helper::vec2set( dvs );
      std::set<int> ivs2 = Helper::vec2set( ivs );
      std::set<int> cvs2 = Helper::vec2set( cvs );
      
      vars.clear();
      faclvl.clear();
      basevar.clear();
      var2group.clear();
      dvs.clear();
      ivs.clear();
      cvs.clear();

      for (int j=0; j<good_cols.size(); j++)
	{	  
	  const std::string & v = vars2[ good_cols[j] ];	  
	  vars.push_back( v );
	  basevar[ v ] = basevar2[ v ];
	  var2group[ v ] = var2group2[ v ];
	  faclvl[ v ] = faclvl2[ v ];

	  if ( dvs2.find( good_cols[j] ) != dvs2.end() )
	    dvs.push_back( j );
	  if ( ivs2.find( good_cols[j] ) != ivs2.end() )
	    ivs.push_back( j );
	  if ( cvs2.find( good_cols[j] ) != cvs2.end() )
	    cvs.push_back( j );
	  
	}

      logger << "  dropped " << nv - good_cols.size()  << " vars with too many NA values\n";
    }
  else
    logger << "  no vars dropped based on missing-value requirements\n";
}


// QC matrix
void gpa_t::qc( const double winsor , const bool stats_mode )
{

  logger << "  running QC (add 'qc=F' to skip)";
  if ( winsor > 0 ) logger << " with winsor=" << winsor << "\n";
  else logger << " without winsorization (to set, e.g. 'winsor=0.05')\n";

  // 1) case-wise deletion

  std::vector<int> retained;
  const int ni = X.rows();
  const int nv = X.cols();

  
  if ( ! retain_rows )
    {
      for (int i=0;i<ni;i++)
	{
	  int num_missing = 0;
	  const Eigen::VectorXd & row = X.row(i);
	  for (int j=0; j<nv; j++)
	    if ( std::isnan( row[j] ) ) ++num_missing;
	  if ( num_missing == 0 ) retained.push_back( i );
	  else if ( verbose )
	    logger << "  dropping indiv. " << ids[i] << " due to missing values (case-wise deletion)\n";
	}
      
      if ( retained.size() < ni )
	{
	  Eigen::MatrixXd X1 = X( retained, Eigen::all );
	  X = X1;
	  
	  std::vector<std::string> ids2 = ids;
	  ids.clear();
	  for (int i=0;i<retained.size();i++)
	    ids.push_back( ids2[retained[i]] );
	  
	  logger << "  case-wise deletion subsetted X from " << ni
		 << " to " << X.rows() << " indivs (add 'verbose' to list)\n";
	  
	}
      else
	logger << "  retained all observations following case-wise deletion screen\n";
    }
  
  // nothing left?
  if ( X.rows() == 0 ) return;

  // if dumping means/SDs, then we don't want to normalize
  if ( stats_mode ) return;

  
  // 2) robust norm & winsorize
  //    but only for DVs
  
  std::set<int> s1 = Helper::vec2set( ivs );
  std::set<int> s2 = Helper::vec2set( cvs );
  s1.insert( s2.begin(), s2.end() );
  if ( s1.size() != ivs.size() + cvs.size() )
    Helper::halt( "overlapping terms in X and Z" );  
  std::vector<int> v1 = Helper::set2vec( s1 );
  
  // copy X and Z vars (will retain originals) 
  Eigen::MatrixXd XZ = X( Eigen::all , v1 );
  
  std::vector<int> zeros;

  eigen_ops::robust_scale( X ,
			   true , // mean center
			   true , // normalize
			   winsor ,
			   true , // second-round norm (post winsorizing)
			   true , // ignore invariants
			   & zeros ); // but track them here

  // swap X and Z values back
  for (int j=0; j<v1.size(); j++)
    X.col( v1[j] ) = XZ.col(j);
  

  // remove any dead cols
  if ( zeros.size() && ! retain_rows )
    {
      std::set<int> z = Helper::vec2set( zeros );
      std::vector<int> nonzeros;
      const int nv = X.cols();
      for (int j=0; j<nv; j++)
	{
	  if ( z.find( j ) == z.end() )
	    nonzeros.push_back( j );
	  else
	    logger << "  dropping " << vars[j] << " due to invariance\n";
	}
      
      Eigen::MatrixXd X1 = X( Eigen::all , nonzeros );
      X = X1;
      
      // now need to update vars, var2group, basevar
      // okay to keep faclvl as is

      // var labels
      std::vector<std::string> vars2 = vars;
      std::map<std::string,std::map<std::string,std::string> > faclvl2 = faclvl;
      std::map<std::string,std::string> basevar2 = basevar;
      std::map<std::string,std::string> var2group2 = var2group;;

      std::set<int> dvs2 = Helper::vec2set( dvs );
      std::set<int> ivs2 = Helper::vec2set( ivs );
      std::set<int> cvs2 = Helper::vec2set( cvs );
      
      vars.clear();
      faclvl.clear();
      basevar.clear();
      var2group.clear();
      dvs.clear();
      ivs.clear();
      cvs.clear();
      

      for (int j=0; j<nonzeros.size(); j++)
	{
	  
	  const std::string & v = vars2[ nonzeros[j] ];	  
	  vars.push_back( v );
	  basevar[ v ] = basevar2[ v ];
	  var2group[ v ] = var2group2[ v ];
	  faclvl[ v ] = faclvl2[ v ];

	  if ( dvs2.find( nonzeros[j] ) != dvs2.end() )
	    dvs.push_back( j );
	  if ( ivs2.find( nonzeros[j] ) != ivs2.end() )
	    ivs.push_back( j );
	  if ( cvs2.find( nonzeros[j] ) != cvs2.end() )
	    cvs.push_back( j );
	  
	}
      
      logger << "  reduced data from " << nv << " to " << X.cols() << " vars ('retain-rows' to skip) \n"; 
      
    }
  
  logger << "  standardized and winsorized all Y variables\n";
  
}


// ------------------------------------------------------------
//
//  linmod_t implementation
//   effectively same as cpt_t, but w/out clustering
//
// ------------------------------------------------------------


void linmod_t::set_DV( const Eigen::MatrixXd & Y_ )
{
  
  Y = Y_;

  if ( ni != 0 && ni != Y.rows() )
    Helper::halt( "unequal number of rows" );
  else
    ni = Y.rows();
  ny = Y.cols();
  
}

void linmod_t::set_IV( const Eigen::MatrixXd & X_ )
{
  
  X = X_;
  
  if ( ni != 0 && ni != X.rows() )
    Helper::halt( "unequal number of rows" );
  else
    ni = X.rows();
  
}

void linmod_t::set_IV( const Eigen::VectorXd & X_ , const std::string & n )
{
  
  X = X_;
  
  if ( ni != 0 && ni != X.rows() )
    Helper::halt( "unequal number of rows" );
  else
    ni = X.rows();

  // fix a vector of X=1
  nx = 1;
  xname.resize(1);
  xname[0] = n;
  
}

void linmod_t::set_Z( const Eigen::MatrixXd & Z_ )
{
  Z = Z_;
  
  if ( ni != 0 && ni != Z.rows() )
    Helper::halt( "unequal number of rows" );
  else
    ni = Z.rows();
  
  nz = Z.cols();
  
}




linmod_results_t linmod_t::run( const int nreps , const bool show_progress , linmod_comps_t * comps )
{
  
  linmod_results_t results;

  ni = Y.rows();
  ny = Y.cols();
  nx = X.cols();
  nz = Z.cols();
  
  if ( ni == 0 || nx == 0 || ny == 0 ) Helper::halt( "linmod_t has no obs, or not X/Y vars" );
  

  //
  // Permutation matrix
  //

  Eigen::MatrixXd P = Eigen::MatrixXd::Identity(ni,ni);
  

  
  //
  // Freedman–Lane (following Winkler et al 2014) 
  //

  // intercept
      
  Eigen::VectorXd II = Eigen::VectorXd::Ones( ni );
      
  //
  // Z is intercept plus nuissance variables
  //
  
  Eigen::MatrixXd ZZ( ni , 1 + nz );  
  if ( nz > 0 ) 
    ZZ << II , Z;
  else
    ZZ << II;
  
  //
  // Rz = 1 - ZZ+
  //
  
  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cqr( ZZ );
  Eigen::MatrixXd Zinv = cqr.pseudoInverse();      
  Eigen::MatrixXd Hz = ZZ * Zinv;
  Eigen::MatrixXd Rz = Eigen::MatrixXd::Identity( Hz.rows() , Hz.cols() ) - Hz ;   

  
  //
  // Initiate permutation counters
  //
  
  // family-wise (init. @ 1); uncorrect U below w/in x-loop
  Eigen::ArrayXXd F = Eigen::ArrayXXd::Ones( nx, ny );  
  
  // re-use perms across different X to maintain dependence
  std::vector<std::vector<int> > pord( nreps );
  for (int r=0; r<nreps; r++)
    {
      std::vector<int> a( ni );
      CRandom::random_draw( a );
      pord[r] = a;
    }
  
  // track max stats per replicate
  std::vector<double> max_t( nreps , 0 );
  
  //
  // iterate over each X
  //

  for (int x=0; x<nx; x++)
    {

      // comps? add a new one for this X var
      linmod_comp_t * comp = comps != NULL ? comps->add( xname[x] ) : NULL ; 
      
      //
      // M is full model: intercept + nuissance + design (1 X term)
      //
  
      Eigen::MatrixXd MM( ni , 1 + nz + 1 );
      MM << ZZ , X.col(x) ; 
      Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cqrM( MM );
      Eigen::MatrixXd Minv = cqrM.pseudoInverse();
      Eigen::MatrixXd Hm = MM * Minv;
      Eigen::MatrixXd Rm = Eigen::MatrixXd::Identity( Hm.rows() , Hm.cols() ) - Hm ;

      //
      // Get get observed statistics
      //

      Eigen::MatrixXd YZres = Rz * Y;
      Eigen::MatrixXd B = Minv * YZres;
      Eigen::MatrixXd Yres = Rm * YZres;
      Eigen::MatrixXd VX = ( MM.transpose() * MM ).inverse() ;
      const int nterms = 1 + nz + 1 ; // intercept + covariates + single IV 
      const int idx = nterms - 1;
      
      // also get asymptotic p-values from the original
      Eigen::VectorXd Pasym;
      Eigen::VectorXd T = get_tstats( B.row(idx) , Yres , VX(idx,idx) , ni - nterms , &Pasym );
      
      Eigen::ArrayXd U = Eigen::ArrayXd::Ones( ny );

      // optional, template comparisons
      if ( comp != NULL )
	comp->compare_observed( vname , B.row(idx) ) ;
            
      //
      // begin perms
      //
      
      if ( nreps != 0 && show_progress ) 
	logger << "  ";
      
      for (int r=0; r<nreps; r++)
	{

	  if ( show_progress )
	    {
	      logger << ".";
	      if ( (r+1) % 10 == 0 ) logger << " ";
	      if ( (r+1) % 50 == 0 ) logger << " " << r+1
					    << " perms\n"
					    << ( r+1 == nreps ? "" : "  " ) ;
	    }
	  
	  // shuffle
	  // std::vector<int> pord( ni );
	  // CRandom::random_draw( pord );
      
	  // permutation matrix
	  Eigen::MatrixXd P = Eigen::MatrixXd::Zero( ni , ni );
	  for (int i=0; i<ni; i++) P(i,pord[r][i]) = 1;
      
	  // shuffle RHS
	  Eigen::MatrixXd MM_perm = P * MM;
	  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cqrM( MM_perm );
	  Eigen::MatrixXd Minv_perm = cqrM.pseudoInverse();
	  Eigen::MatrixXd Hm_perm = MM_perm * Minv_perm;
	  Eigen::MatrixXd Rm_perm = Eigen::MatrixXd::Identity( Hm_perm.rows() , Hm_perm.cols() ) - Hm_perm ;
	  
	  // get permuted statistics
	  Eigen::MatrixXd B_perm = Minv_perm * YZres;
	  Eigen::MatrixXd Yres_perm = Rm_perm * YZres;      
	  Eigen::VectorXd T_perm = get_tstats( B_perm.row(idx) , Yres_perm , VX(idx,idx) , ni - nterms );

	  // optional, template comparisons
	  if ( comp != NULL )
	    comp->compare_null( vname , B_perm.row(idx) ) ;

	  // accumulate
	  //double max_t = 0;
	  for (int y=0; y<ny; y++)
	    {
	      double abs_t = fabs( T_perm[y] );

	      // store point-wise emp-P now
	      if ( abs_t >= fabs( T[y] ) ) ++U[y];

	      // and track max per repl. for family-wise correction when done all X
	      // unless is a self-test (Y==X), i.e. ignoring here means we do not
	      // pay price for extra tests (that are never reported)
	      if ( abs_t > max_t[r] )
		if ( xname[x] != vname[y] )
		  max_t[r] = abs_t ;
	    }
	  
		  
	} // next replicate

      if ( nreps != 0 && show_progress ) 
	logger << "\n";
      
      //
      // Get point-wise p-values
      //
      
      U /= (double)(nreps+1);

      Eigen::MatrixXd R( ny , 4 );

      R << B.row(idx).transpose() , T , U , Pasym;
      
      //
      // Point-wise results
      //

      for (int y=0; y<ny; y++)
	{

	  if (  xname[x] != vname[y] ) // ignore any self (Y == X) tests
	    {
	      results.beta[ xname[x] ][ vname[y] ] = R(y,0);
	      results.t[ xname[x] ][ vname[y] ] = R(y,1);
	      results.emp[ xname[x] ][ vname[y] ] = R(y,2);
	      results.p[ xname[x] ][ vname[y] ] = R(y,3);
	    }
	}
            
      // next X
    }

  
  //
  // family-wise corrected P
  //
  
  for (int x=0; x<nx; x++)
    for (int y=0; y<ny; y++)
      {
	const double t = fabs( results.t[ xname[x] ][ vname[y] ] );
	for (int r=0; r<nreps; r++)
	  if ( max_t[r] >= t ) ++F(x,y);
      }
  
  F /= (double)(nreps+1);
  
  for (int x=0; x<nx; x++)
    for (int y=0; y<ny; y++)
      if (  xname[x] != vname[y] ) // ignore self tests
	results.emp_corrected[ xname[x] ][ vname[y] ] = F(x,y);
  
  // and make other adjusted stats
  results.make_corrected( xname , vname );
  
  // all done
  return results;
}



Eigen::VectorXd linmod_t::get_tstats( const Eigen::VectorXd & B ,
				      const Eigen::MatrixXd & Yres ,
				      const double vx ,
				      const int denom ,
				      Eigen::VectorXd * pvalues )
{

  const int n = B.rows();
  
  Eigen::VectorXd T = Eigen::VectorXd::Zero( n );
  for (int i=0; i<n; i++)
    T[i] = Yres.col(i).transpose() * Yres.col(i); 
  
  // OR Eigen::VectorXd T = Yres.transpose() * Yres;
  
  for (int i=0; i<n; i++)
    T[i] = B[i] / sqrt( vx * T[i] / (double)denom ) ;
  
  // residual degrees of freedom = ni - 2 - nz  = denom
  
  // double sigmaSq =  ee[0]  / ( ni - 1 - 1 - nz ) ;
  // double se = sqrt( sigmaSq * vx );

  // for (int i=0; i<n; i++)
  //   std::cout << " p = " << T[i] << " " << MiscMath::pT( T[i] , denom ) << "\n";

  if ( pvalues != NULL )
    {
      *pvalues = Eigen::VectorXd::Zero( n );
      for (int i=0; i<n; i++)
	(*pvalues)[i] = MiscMath::pT( T[i] , denom );
    }
  
  return T;
}


void gpa_t::parse( const std::string & pfile )
{

  std::string pfile1 = Helper::expand( pfile );
  logger << "  parsing " << pfile1 << "\n";

  if ( ! Helper::fileExists( pfile1 ) ) Helper::halt( "could not open " + pfile );
  
  std::ifstream IN1( pfile1.c_str() , std::ios::in );
    
  //json doc{json::parse(IN1)};

  try 
    {
      json doc = json::parse(IN1);
      
      bool has_inputs = doc.count( "inputs" ); 
      bool has_specs  = doc.count( "specs" );
      
      // general specs: only supports a limited set on prep-mode
      
      // vars -> incvars
      // xvars -> excvars  
      // facs -> incfacs
      // xfacs -> excfacs
      // grps -> incgrps
      // xgrps -> excgrps
      
      
      // ** these not supported in prep-mode (or JSON input)  
      // nvars -> incnums
      // xnvars -> excnums
      // lvars -> inclvars (long var names)
      // xlvars -> exclvars 
      // faclvls -> incfaclvls
      // xfaclvls -> excfaclvls.clear();
      
      if ( has_specs )
	{
	  
	  logger << "  reading general specifications ('specs') from " << pfile << "\n";
	  
	  json s = doc[ "specs" ];
	  
	  std::vector<std::string> tok;
	  
	  if ( s.contains( "vars" ) )
	    {
	      json x = s[ "vars" ];
	      if ( x.is_string() ) x = std::vector<std::string>(1,x);
	      if ( x.is_array() ) incvars = Helper::vec2set( x.get<std::vector<std::string>>() );
	    }
	  
	  if ( s.contains( "xvars" ) )
	    {
	      json x = s[ "xvars" ];
	      if ( x.is_string() ) x = std::vector<std::string>(1,x);
	      if ( x.is_array() ) excvars = Helper::vec2set( x.get<std::vector<std::string>>() );
	    }
	  
	  if ( s.contains( "facs" ) )
	    {
	      json x = s[ "facs" ];
	      if ( x.is_string() ) x = std::vector<std::string>(1,x);
	      if ( x.is_array() ) incfacs = Helper::vec2set( x.get<std::vector<std::string>>() );
	    }
	  
	  if ( s.contains( "xfacs" ) )
	    {
	      json x = s[ "xfacs" ];
	      if ( x.is_string() ) x = std::vector<std::string>(1,x);
	      if ( x.is_array() ) excfacs = Helper::vec2set( x.get<std::vector<std::string>>() );
	    }
	  
	  if ( s.contains( "grps" ) )
	    {
	      json x = s[ "grps" ];
	      if ( x.is_string() ) x = std::vector<std::string>(1,x);
	      if ( x.is_array() ) incgrps = Helper::vec2set( x.get<std::vector<std::string>>() );
	    }
	  
	  if ( s.contains( "xgrps" ) )
	    {
	      json x = s[ "xgrps" ];
	      if ( x.is_string() ) x = std::vector<std::string>(1,x);
	      if ( x.is_array() ) excgrps = Helper::vec2set( x.get<std::vector<std::string>>() );
	    }
	  
	}
      
      
      //
      // parse inputs
      //
      
      if ( has_inputs )
	{
	  
	  logger << "  reading file specifications ('inputs') for " << doc[ "inputs" ].size() << " files:\n";
	  
	  for (auto & item : doc[ "inputs" ] ) {
	    
	    // expect minimally, group and file
	    if ( ! item.contains( "group" ) ) Helper::halt( "expecting 'group' key for all inputs in " + pfile );
	    if ( ! item.contains( "file" ) ) Helper::halt( "expecting 'file' key for all inputs in " + pfile );
	    
	    std::string file_name = item[ "file" ];
	    std::string file_group = item[ "group" ];
	    
	    logger << "   " << file_name << " ( group = " << file_group << " ): ";
	    
	    //
	    // and also
	    //   1) define factors for this file 
	    //      (mirroring command line form: inputs=file|group|fac1|fac2
	    //
	    //   2) allow for vars, xvars ( --> file specific versions) for prep only
	    //
	    //   3) mappings (for a given variable, label --> numeric) 
	    //
	    //   4) allowing a fixed faclvl to be added to a file
	    //
	    //   5) allowing a file to be designated as IV only (i.e. will not be added as a DV by default) 
	    //
	    
	    //
	    // vars
	    //
	    
	    if ( item.contains( "vars" ) )
	      {
		
		json x = item[ "vars" ] ;
		
		if ( x.is_string() )
		  x = std::vector<std::string>( 1 , x );
		
		if ( ! x.is_array() )
		  {
		    if ( ! x.is_null() ) 
		      logger << "  *** expecting vars: [ array ] in " << pfile << ", skipping...\n";
		  }
		else
		  {		
		    for (auto v : x )
		      {
			if ( v.is_string() )
			  {
			    std::string var = v;
			    file2incvars[ file_name ].insert( var );
			  }
			else if ( v.is_object() )
			  {
			    for ( auto & vv : v.items() )
			      {
				std::string var = vv.key();
				std::string val = vv.value();
				file2var2alias[ file_name ][ var ] = val;
				file2incvars[ file_name ].insert( val ); // nb. is new, aliased name
			      }
			  }
		      }
		    
		  }
	      }
	    
	    
	    //
	    // xvars (no aliasing)
	    //
	    
	    if ( item.contains( "xvars" ) )
	      {
		
		json x = item[ "xvars" ] ;
		if ( x.is_string() )
		  x = std::vector<std::string>( 1 , x );
		
		if ( ! x.is_array() )
		  {
		    if ( ! x.is_null() )
		      logger << "  *** expecting xvars: [ array ] in " << pfile << ", skipping...\n";
		  }
		else
		  {		
		    for (auto v : x )
		      {
			if ( v.is_string() )
			  {
			    std::string var = v;
			    file2excvars[ file_name ].insert( var );
			  }
		      }
		  }		
	      }
	    
	    
	    //
	    // specify facs (i.e. matching inputs=[] --> infiles
	    //   ;;; and allow aliasing here too
	    
	    std::set<std::string> file_facs;
	    
	    if ( item.contains( "facs" ) )
	      {
		json x = item[ "facs" ] ;
		if ( x.is_string() )
		  x = std::vector<std::string>( 1 , x );
		
		if ( ! x.is_array() )
		  {
		    if ( ! x.is_null() )
		      logger << "  *** expecting facs: [ array ] in " << pfile << ", skipping...\n";
		  }
		else
		  {		
		    for (auto v : x )
		      {
			if ( v.is_string() )
			  {
			    std::string var = v;
			    file_facs.insert( var );
			  }
			else if ( v.is_object() )
			  {
			    for ( auto & vv : v.items() )
			      {
				std::string var = vv.key();
				std::string val = vv.value();
				file2var2alias[ file_name ][ var ] = val; // add generic alias
				file_facs.insert( val ); // nb. is new, aliased name to be searched for                                                              
			      }
			    
			  }		
			
		      }
		  }
	      }
	    
	    
	    //
	    // fixed faclvls
	    //
	    
	    std::map<std::string,std::string> file_fixed;
	    
	    if ( item.contains( "fixed" ) )
	      {

		json x = item[ "fixed" ] ;
		if ( x.is_string() )
		  x = std::vector<std::string>( 1 , x );
		
		if ( ! x.is_array() )
		  {
		    if ( ! x.is_null() ) 
		      logger << "  *** expecting fixed: [ array ] in " << pfile << ", skipping...\n";
		  }
		else
		  {		
		    for (auto v : x )
		      {
			// expecting [ {  str : str } , ... ]
			
			if ( v.is_object() )
			  {			
			    for ( auto & vv : v.items() )
			      {
				std::string fac = vv.key();
				std::string lvl = vv.value();				
				if ( file_facs.find( fac ) != file_facs.end() )
				  Helper::halt( "cannot specify a fixed factor that is also a named factor in the file" );
				
				// store
				file_fixed[ fac ] = lvl; 
			      
			      }
			  }
		      }
		    
		  }
	      }
	    
		    
	    //
	    // mappings
	    //
	    
	    if ( item.contains( "mappings" ) )
	      {
		
		json x = item[ "mappings" ] ;
		if ( x.is_string() )
		  x = std::vector<std::string>( 1 , x );

		if ( ! x.is_array() )
		  {
		    if ( ! x.is_null() ) 
		      logger << "  *** expecting mappings: [ array ] in " << pfile << ", skipping...\n";
		  }
		else
		  {		
		    for (auto v : x )
		      {
			// expecting [ {  var : { str , num } } , ... ]
			
			// if using aliases, should be in new alias form
			if ( v.is_object() )
			  {			
			    for ( auto & vv : v.items() )
			      {
				std::string var = vv.key();
				if ( vv.value().is_object() && vv.value().size() == 1 )
				  {
				    for ( auto & vvv : vv.value().items() )
				      {
					std::string str = vvv.key();
					double num = vvv.value();					
					// store
					file2var2mapping[ file_name ][ var ][ str ] = num;
				      }
				  }
			      }
			  }
		      }
		    
		  }
	      }
	    
	
	    
	    //
	    // add in
	    //
	    
	    infiles[ file_name ] = file_facs;
	    file2group[ file_name ] = file_group;
	    file2fixed[ file_name ] = file_fixed;
	    
	    logger << "\n    expecting " << file_facs.size() << " factors";
	    if  (file_fixed.size() ) logger << " and " << file_fixed.size() << " fixed factors";
	    if ( file2incvars[ file_name ].size() ) logger << ", extracting " << file2incvars[ file_name ].size() << " var(s)";
	    if ( file2excvars[ file_name ].size() ) logger << ", ignoring " << file2excvars[ file_name ].size() << " var(s)";	
	    if ( file2incvars[ file_name ].size() == 0 && file2excvars[ file_name ].size() == 0 ) logger << ", reading all var(s)";
	    logger << "\n";
	
	
	  }
	
	}
  
    }
  catch(const std::exception& e)
    {
      Helper::halt( "problem parsing JSON file " + pfile1 + ":\n --> " + e.what() );
    }

}

// kNN imputation of missing points
void gpa_t::knn_imputation( param_t & param )
{
  
  int k = param.requires_int( "knn" );
  
  if ( k == 0 ) return;
  
  // only do this for DVs; but copy whole X (all vars) for simplicity;
  // but only do updates for dvs[]
  
  const int ni = X.rows();
  const int nv = X.cols();
  const int ndv = dvs.size();
  
  //
  // count missing data per indiv.
  //

  std::vector< std::set<int> > nmissing( ni );
  
  for (int i=0; i<ni; i++)
    {
      for (int j=0; j<ndv; j++)
	{
	  if ( std::isnan( X(i,dvs[j]) ) )
	    nmissing[i].insert(dvs[j]);
	}
    }


  //
  // Get a normalized version of X (special case handling existing NaNs)
  //   --> must be room for improvement here, but brief searching around
  //       suggests something like this needed for Eigen

  
  // extract only DVs (no missingness in IV or covariates allowed)
  // so just have all non-DV variables as 0, i.e. so they won't feature
  // in the neighbour calcs
  
  Eigen::MatrixXd Z = Eigen::MatrixXd::Zero( ni , nv ); 
  
  for (int j0=0; j0<ndv; j0++)
    {
      
      // get X index
      const int j = dvs[j0];

      // add column in (otherwise we leave it blank)
      Z.col(j) = X.col(j);

      // normalize
      Eigen::VectorXd col = Z.col(j);
      
      // get indices of missing elements
      Eigen::ArrayXi idx = col.array().isNaN().cast<int>();
      
      const int non_missing = ni - idx.count() ;

      if ( non_missing > 1 )
	{
	  
	  // make a std::vector<int> mask      
	  std::vector<int> mask( non_missing );	  
	  int zi = 0;
	  for (int z = 0; z < ni; z++ ) 
	    if (!idx(z)) mask[zi++] = z;
	  
	  // normalize non-NaN values
	  double mean = Z.col(j)(mask).mean();
	  double sd = eigen_ops::sdev( Z.col(j)(mask) );

	  Z.col(j)(mask)= (Z.col(j)(mask).array() - mean ) / ( sd > 0 ? sd : 1 );;

	}

    }
   
  //
  // process each indiv
  //


  int nans_imputed = 0 , nans_remaining = 0;
  std::set<int> indivs_with_nans;
  
  for (int i=0; i<ni; i++)
    {

      int m = nmissing[i].size();
      
      if ( m == 0 ) continue;

      // require at least 50% of variables to be present
      // or else leave missing vars as is
      if ( m >= nv / 2 ) continue;
      
      // get data (each col is normalized)
      const Eigen::VectorXd & f1 = Z.row(i);

      // rank by Euclidean distance to all others
      std::map<double,int> nearest;
      for (int k=0; k<ni; k++)
	{
	  if ( k == i ) continue;
	  
	  const Eigen::VectorXd & f2 = Z.row(k);	  
	  double dst = (f1 - f2).array().isNaN().select(0,f1-f2).squaredNorm();
	  int n = ni - (f1 - f2).array().isNaN().count() - 1;
	  // weighted Euclidean distance
	  nearest[ sqrt( n / (double)(ni-1) * dst ) ] = k;
	}

      // now impute each missing spot

      const std::set<int> & mvars = nmissing[i];

      std::set<int>::const_iterator mm = mvars.begin();
      while ( mm != mvars.end() )
	{
	  // get 'k' people to fill this
	  // may have to skip some if missing
	  int c = 0;

	  double imp = 0;
	  
	  // X(i,*mm) should be NaN -- check
	  //std::cout << " X(i,*mm) should be missing = " << X(i,*mm)  << "\n";
	  
	  std::map<double,int>::const_iterator nn = nearest.begin();
	  while ( nn != nearest.end() )
	    {

	      // neighbour missing
	      if ( std::isnan( X( nn->second , *mm ) ) )
		{
		  ++nn;
		  continue;
		}

	      // addup
	      imp += X( nn->second , *mm );
	      
	      // done?
	      if ( ++c == k ) break;

	      ++nn;
	    }

	  // update, if possible
	  if ( c > 0 ) 
	    {
	      nans_imputed++;
	      X(i,*mm) = imp / (double)c ; 
	    }
	  else
	    {
	      nans_remaining++;
	      indivs_with_nans.insert( i );
	    }
	  // next variable for this person
	  ++mm;
	}
      
      // next person
    }

  logger << "  finished kNN imputation: imputed " << nans_imputed << " missing values\n";

  if ( nans_remaining )
    logger << "  " << nans_remaining << " missing values remaining, from "
	   << "  " << indivs_with_nans.size() << " individuals\n";
  else
    logger << "  all missing values imputed\n";

}
  
Eigen::MatrixXd linmod_t::correct( const Eigen::VectorXd & p )
{

  // 4 corrected values: Bonf, Holm, BH, BY

  // Benjamini, Y., and Hochberg, Y. (1995). Controlling the false
  // discovery rate: a practical and powerful approach to multiple
  // testing. Journal of the Royal Statistical Society Series B, 57,
  // 289–300.

  // Benjamini, Y., and Yekutieli, D. (2001). The control of the false
  // discovery rate in multiple testing under dependency. Annals of
  // Statistics 29, 1165–1188.

  //  Yekutieli and Benjamini (1999) developed another FDR controlling
  //  procedure, which can be used under general dependence by
  //  resampling the null distribution.
  
  if ( p.size() == 0 )
    return Eigen::MatrixXd::Zero(0,4);
  
  if ( p.size() == 1 ) 
    {
      Eigen::MatrixXd R = Eigen::MatrixXd::Zero(1,4);
      R(0,0) = R(0,1) = R(0,2) = R(0,3) = p[0];
      return R;
    }

  struct pair_t {
    double p;
    int l;
    bool operator< (const pair_t & p2) const { return ( p < p2.p ); }
  };

  std::vector<pair_t> sp;

  const int npv = p.size();
  
  for (int l=0; l<npv; l++)
    {      
      if ( p[l] > -1) 
	{
	  pair_t pt;
	  pt.p = p[l];
	  pt.l = l;      
	  sp.push_back(pt);
	}      
    }
  
  // sort p-values
  std::sort(sp.begin(),sp.end());

  double t = (double)sp.size();
  int ti = sp.size();
  
  // Consider each test
  Eigen::VectorXd pv_holm = Eigen::VectorXd::Zero( ti );
  Eigen::VectorXd pv_BH = Eigen::VectorXd::Zero( ti );
  Eigen::VectorXd pv_BY = Eigen::VectorXd::Zero( ti );
  
  // check that we have >1 valid p-value
  if ( ti == 1 )
    {
      pv_holm[0] = pv_BH[0] = pv_BY[0] = sp[0].p;
    } 
  else if ( ti > 1 )
    {  
      // Holm 
      pv_holm[0] = sp[0].p*t > 1 ? 1 : sp[0].p*t;
      for (int i=1;i<ti;i++)
	{
	  double x = (ti-i)*sp[i].p < 1 ? (ti-i)*sp[i].p : 1;
	  pv_holm[i] = pv_holm[i-1] > x ? pv_holm[i-1] : x;
	}
  
      // // Sidak SS
      // for (int i=0;i<ti;i++)
      //   pv_sidakSS[i] = 1 - pow( 1 - sp[i].p , t );
      //
      // // Sidak SD
      // pv_sidakSD[0] = 1 - pow( 1 - sp[0].p , t );
      // for (int i=1;i<ti;i++)
      //   {
      //     double x = 1 - pow( 1 - sp[i].p , t - i  );
      //     pv_sidakSD[i] = pv_sidakSD[i-1] > x ? pv_sidakSD[i-1] : x ; 
      //   }
  
      // BH
      pv_BH[ti-1] = sp[ti-1].p;
      for (int i=ti-2;i>=0;i--)
	{
	  double x = (t/(double)(i+1))*sp[i].p < 1 ? (t/(double)(i+1))*sp[i].p : 1 ;
	  pv_BH[i] = pv_BH[i+1] < x ? pv_BH[i+1] : x;
	}
      
      // BY
      double a = 0;
      for (double i=1; i<=t; i++)
	a += 1/i;
      pv_BY[ti-1] = a * sp[ti-1].p < 1 ? a * sp[ti-1].p : 1 ; 
      for (int i=ti-2;i>=0;i--)
	{      
	  double x = ((t*a)/(double)(i+1))*sp[i].p < 1 ? ((t*a)/(double)(i+1))*sp[i].p : 1 ;
	  pv_BY[i] = pv_BY[i+1] < x ? pv_BY[i+1] : x;
	}
      
    }
  
  // 
  // update original p-value list
  //
  
  Eigen::MatrixXd res = Eigen::MatrixXd::Constant( npv , 4 , -9 );
  
  for (int l=0; l<ti; l++)
    {
      const int idx = sp[l].l;
      
      // bonferroni
      res(idx,0) = sp[l].p*t > 1 ? 1 : sp[l].p*t;
      // others
      res(idx,1) = pv_holm[l];
      res(idx,2) = pv_BH[l]; // fdr
      res(idx,3) = pv_BY[l];
      
    }

  return res;
}



// generate corrected results
void linmod_results_t::make_corrected( const std::vector<std::string> & xvars ,
				       const std::vector<std::string> & yvars )
  
{
  
  // index results (will always be all-by-all)
  int nx = xvars.size();
  int ny = yvars.size();
  
  index.clear();
  int cnt = 0;
  
  for (int i=0; i<xvars.size(); i++)
    for (int j=0; j<yvars.size(); j++)
      index[ xvars[i] ][ yvars[j] ] = cnt++;
    
  // make vector of asymptotic p-values
  
  Eigen::VectorXd P = Eigen::VectorXd::Zero( cnt );

  cnt = 0;
  for (int i=0; i<xvars.size(); i++)
    for	(int j=0; j<yvars.size(); j++)
      P[ cnt++ ] = p[ xvars[i] ][ yvars[j] ];
    
  corr = linmod_t::correct( P );
  
}

double linmod_results_t::bonf(const std::string & x , const std::string & y )
  { return corr( index[ x ][ y ] , 0 ); }

// double linmod_results_t::sidak_ss(const std::string & x , const std::string & y )
//   { return corr( index[ x ][ y ] , 1 ); }

// double linmod_results_t::sidak_sd(const std::string & x , const std::string & y )
//   { return corr( index[ x ][ y ] , 2 ); }

double linmod_results_t::holm(const std::string & x , const std::string & y )
  { return corr( index[ x ][ y ] , 1 ); }

double linmod_results_t::fdr_bh(const std::string & x , const std::string & y )
  { return corr( index[ x ][ y ] , 2 ); }

double linmod_results_t::fdr_by(const std::string & x , const std::string & y )
  { return corr( index[ x ][ y ] , 3 ); }



//
// comparison mode
//

linmod_comp_t * linmod_comps_t::add( const std::string & xvar )
{
  // insert if absent
  comp_map.try_emplace(xvar, this);
  return comp( xvar );
}

void linmod_comps_t::load( const std::string & file )
{

  deff.clear();
  
  std::string filename = Helper::expand( file );

  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not load " + filename );

  std::ifstream IN1( filename.c_str() , std::ios::in );

  int t = 0;
  
  std::string line;
  while ( Helper::safe_getline( IN1 , line ) )
    {
      if ( line.empty() ) continue;
      
      // whitespce delimited  :   XVAR YVAR DEFF
      std::vector<std::string> tok = Helper::parse( line, " \t" );

      double d;
      
      if ( tok.size() != 3 )	
	Helper::halt( "bad line: " + line );
      
      if ( ! Helper::str2dbl( tok[2] , &d ) )
	Helper::halt( "bad line: " + line );	
      
      if ( d != 0 )
	{
	  deff[ tok[1] ] = d > 0 ? 1 : -1;
	  ++t;
	}
    }
  
  logger << "  read " << t << " direction-of-effect values from " << filename << "\n"; 
  
  IN1.close();
  
}


void linmod_comp_t::compare_observed( const std::vector<std::string> & yvars , const Eigen::VectorXd & b )
{
  // also set # of tests actually found
  stat1 = compare( yvars, b , &nt );
}

void linmod_comp_t::compare_null( const std::vector<std::string> & yvars , const Eigen::VectorXd & b )
{
  double s0 = compare( yvars, b );
  stat0.push_back( s0 );
}

double linmod_comp_t::compare( const std::vector<std::string> & yvars , const Eigen::VectorXd & b , int * nt )
{
  double s = 0;
  
  const int n = yvars.size();
  int cnt = 0;
  
  for (int i=0; i<n; i++)
    {      
      std::map<std::string,int>::const_iterator yy = parent->deff.find( yvars[i] );      
      if ( yy != parent->deff.end() )
	{
	  // found
	  ++cnt;
	  
	  // passed
	  int sign1 = b[i] > 0 ? +1 : -1;
	  bool zero1 = b[i] == 0 ;
	  
	  // template
	  int sign2 = yy->second  > 0 ? +1 : -1;

	  // no weighting etc
	  if ( sign1 == sign2 && ! zero1 )
	    s++;  

	}      
    }
  
  // number of tests found 
  if ( nt != NULL ) *nt = cnt;
  
  // return statistic 
  return s;
}

void linmod_comp_t::evaluate()
{

  int r = 0;
  int n = stat0.size();
  if ( n == 0 ) return;
  for (int i=0; i<n; i++)
    if ( stat0[i] >= stat1 ) ++r;
  pval = (r+1) / double(n+1);

  xval = MiscMath::mean( stat0 );
  double sd = MiscMath::sdev( stat0 , xval );
  zval = ( stat1 - xval ) / sd;
  
  // two-sided p-value: statistic is difference from expectation (which we get from stat0)
  const double t1 = fabs( stat1 - xval );
  
  r = 0;
  for (int i=0; i<n; i++)
    if ( fabs( stat0[i] - xval ) >= t1 ) ++r;

  pval2 = (r+1) / double(n+1);

}

double linmod_comp_t::stat() const
{
  return stat1;
}

double linmod_comp_t::expected() const
{
  return xval;
}

double linmod_comp_t::zstat() const
{
  return zval;
}

double linmod_comp_t::pvalue() const
{
  return pval;
}

double linmod_comp_t::pvalue2() const
{
  return pval2;
}





