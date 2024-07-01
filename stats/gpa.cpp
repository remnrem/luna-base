
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

#include "helper/logger.h"
#include "db/db.h"
#include "stats/eigen_ops.h"
#include "miscmath/crandom.h"

extern logger_t logger;
extern writer_t writer;


gpa_t::gpa_t( param_t & param , const bool prep_mode )
{

  // invoked in one of two modes:
  //   1) prepare a binary data file  [ --prep-gpa ] --> prep_mode == T 
  //   2) apply associaiton models    [ --gpa ]

  // prep:
  //   inputs=file1|grp|X|Y|Z,file2|grp,file|grp|CH|F , i.e. always group name, then any factors
  //   dat=<bfile name>
  //   vars=<set of vars to include only>
  //   xvars=<set of vars to exclude from inputs>

  bfile = param.requires( "dat" );

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
  
  if ( prep_mode )
    {

      std::set<std::string> files = param.strset( "inputs" );
      if ( files.size() == 0 ) return;
      
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
	      for (int j=2;j<tok.size(); j++) 
		infiles[ tok[0] ].insert( tok[j] );
	    }
	  
	  ++ff;
	}
    }

  
  //
  // include/exclude vars
  //

  if ( param.has( "vars" ) )
    incvars = param.strset( "vars" );
  if ( param.has( "xvars" ) )
    excvars = param.strset( "xvars" );

  //
  // include vars based on presence/absence of a factor
  //

  if ( param.has( "facs" ) )
    incfacs = param.strset( "facs" );
  if ( param.has( "xfacs" ) )
    excfacs = param.strset( "xfacs" );


  //
  // include/exclude based on (file) group
  //
  
  if ( param.has( "grps" ) )
    incgrps = param.strset( "grps" );
  if ( param.has( "xgrps" ) )
    excgrps = param.strset( "xgrps" );

  //
  // inc/exc based on fac/lvl pairs (for now, can only do on read)
  //

  if ( prep_mode )
    {
      if (  param.has( "faclvls" ) ||  param.has( "xfaclvls" ) )
        Helper::halt( "cannot specify faclvls/xfaclvls with --gpa-prep" );
      if (  param.has( "grps" ) ||  param.has( "xgrps" ) )
        Helper::halt( "cannot specify grps/xgrps with --gpa-prep" );
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
  
  
  // nvars/xnvars only allowed in read/run mode
  if ( prep_mode )
    {
      if (  param.has( "nvars" ) ||  param.has( "xnvars" ) )
	Helper::halt( "cannot specify nvars/xnvars with --gpa-prep" );
    }

  
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
  // read data / make binary file? 
  //

  if ( prep_mode )
    prep();
  else
    {
      logger << "  reading binary data from " << bfile << "\n";

      // read data in  ( and this does filtering of columns) 
      read();

      // secondarily, subset to a smaller # of rows (e.g. case-only analysis)
      if ( param.has( "subset" ) || param.has( "ids" ) )
	{
	  std::set<std::string> sub_ids, sub_cols;
	  
	  if ( param.has( "ids" ) ) sub_ids = param.strset( "ids" );
	  if ( param.has( "subset" ) ) sub_cols = param.strset( "subset" );
	  std::set<int> rows;
	  std::map<int,bool> cols;
	  const int ni = X.rows();
	  const int nv = X.cols();
	  // build ID rows
	  for (int i=0; i<ni; i++)
	    if ( sub_ids.find( ids[i] ) != sub_ids.end() )
	      rows.insert( i );

	  // build cols (i.e. searching for non-null value to include, not NaN or missing)
	  // allows for each term to be a pos or neg match subset=-MALE implies MALE == 0 (-->F)
	  //  where subset=MALE or subset=+MALE implies --> M
	  
	  for (int j=0; j<nv; j++)
	    {	      
	      if ( sub_cols.find( vars[j] ) != sub_cols.end() )
		cols[j] = true;
	      else if ( sub_cols.find( "+" + vars[j] ) != sub_cols.end() )
		cols[j] = true;
	      else if ( sub_cols.find( "-" + vars[j] ) != sub_cols.end() )
		cols[j] = false;

	    }
	  subset( rows, cols );
	}

      
      //
      // select predictors (X) and covariates (Z) - assume everything else is
      // a DV (i.e. sleep metric) unless explicitly told so
      //
      
      if ( param.has( "X" ) )
	{
	  std::set<std::string> v = param.strset( "X" );
	  std::set<std::string> found;
	  const int nv = vars.size();
	  ivs.clear();
	  for (int j=0; j<nv; j++)
	    if ( v.find( vars[j] ) != v.end() )
	      {
		ivs.push_back( j );
		found.insert( vars[j] );
	      }
	  if ( found.size() < v.size() )
	    {
	      std::set<std::string>::const_iterator vv = v.begin();
	      while ( vv != v.end() )
		{
		  if ( found.find( *vv ) == found.end() )
		    logger << "  *** warning, could not find " << *vv << "\n";
		  ++vv;
		}
	    }		  
	}
      
      if ( param.has( "Z" ) )
	{
	  std::set<std::string> v = param.strset( "Z" );
	  std::set<std::string> found;
	  const int nv = vars.size();
	  cvs.clear();
	  for (int j=0; j<nv; j++)
	    {
	      if ( v.find( vars[j] ) != v.end() )
		{
		  cvs.push_back( j );
		  found.insert( vars[j] );
		}
	    }
	  if ( found.size() < v.size() )
	    {
	      std::set<std::string>::const_iterator vv = v.begin();
	      while ( vv != v.end() )
		{
		  if ( found.find( *vv ) == found.end() )
		    logger << "  *** warning, could not find " << *vv << "\n";
		  ++vv;
		}
	    }
	}
      
      // Y is everything else that is left
      //  (i.e. post any prior col selection)
      std::set<int> v1 = Helper::vec2set( ivs );
      std::set<int> v2 = Helper::vec2set( cvs );      
      
      const int nv = vars.size();      
      for (int j=0; j<vars.size(); j++)
	if ( v1.find( j ) == v1.end() && v2.find( j ) == v2.end() )
	  dvs.push_back( j );

      logger << "  selected " 
	     << ivs.size() << " X vars & "
	     << cvs.size() << " Z vars, implying "
	     << dvs.size() << " Y vars\n";
      
      //
      // now do QC (on DVs only) - i.e. keeps 
      //

      double winsor_th = param.has( "winsor" ) ? param.requires_dbl( "winsor" ) : 0.01 ; 
      if ( winsor_th < 0 || winsor_th > 0.2 )
	Helper::halt( "winsor must be set between 0 and 0.2" ); 

      // can skip QC with qc=F option
      if ( ( ! param.has( "qc" ) ) || param.yesno( "qc" ) )
	qc( winsor_th );
      

      //
      // optionally dump variables and/or manifest?
      //
      
      if ( param.has( "dump" ) )
	dump();
      
      if ( param.has( "manifest" ) )
	manifest();

      //
      // run association tests
      //

      if ( dvs.size() != 0 && ivs.size() != 0 )
	{	  
	  // options
	  nreps = param.requires_int( "nreps" ) ;
	  if ( nreps < 0 ) Helper::halt( "nreps must be positive" );

	  pthresh = param.has( "p" ) ? param.requires_dbl( "p" ) : 99 ;
	  pthresh_adj = param.has( "padj" ) ? param.requires_dbl( "padj" ) : 99 ;

	  correct_all_X = param.has( "correct-all-X" ) ? param.yesno( "correct-all-X" ) : false;

	  if ( correct_all_X ) 
	    {
	      logger << "  adjusting for multiple tests across all X variables\n";
	      logger << "  performing association tests w/ " << nreps << " permutations... (may take a while)\n"; 
	      run();
	    }
	  else
	    {
	      logger << "  adjusting for multiple tests only within each X variable\n";
	      logger << "  performing association tests w/ " << nreps << " permutations... (may take a while)\n";	  
	      run1X();
	    }
	  
	  logger << "  ...done\n";
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
  
  std::set<std::string> allfacs;
  
  std::map<int,std::map<int,double> > D; // nb. D[ var ][ ind ] --> value

  // track lvl order over all files
  std::map<std::string,std::map<std::string,int> >  faclvlcnt;
	
  // iterate over all files
  std::map<std::string,std::set<std::string> >::const_iterator ff = infiles.begin();
  while ( ff != infiles.end() )
    {

      // file-specific counts
      std::set<std::string> file_ids, file_bvars, file_evars;
            
      // named factors for this file
      const std::set<std::string> & facs = ff->second;

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
	      logger << "  skipping " << ff->first << " due to facs requirement\n";
	      ++ff;
	      continue;
	    }
	}

      if ( excfacs.size() )
	{
          const int req_matched =  Helper::nmatches( excfacs , facs );
          if ( req_matched == excfacs.size() && req_matched == facs.size() )
            {
	      logger <<	"  skipping " << ff->first << " due to xfacs requirement\n";
              ++ff;
              continue;
            }
        }

      
      //
      // read int
      //
      
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
	  logger << "  *** skipping " << ff->first << " ... empty file\n";
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
	      allfacs.insert( tok[j] );
	      continue;
	    }
	      
	  // skip?
	  if ( incvars.size() != 0 && incvars.find( tok[j] ) == incvars.end() ) continue;
	  if ( excvars.find( tok[j] ) != excvars.end() ) continue;
	  
	  
	  // read this var
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
	  ++ff;
	  logger << "  *** skipping " << ff->first << ", no selected (non-factor) variables\n";
	  continue;
	}
      
      // not all factors found?
      if ( facs.size() != fac2slot.size() )
	Helper::halt( "not all factors found for " + ff->first );

      
      //
      // read rows ( to find unique IDs and unique VAR+FACLVL combos
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
	    Helper::halt( "bad line - col # doesn't match header" );

	  // a new ID?
	  const std::string & id = dtok[ id_col ] ;
	  file_ids.insert( id );
	  if ( id2slot.find( id ) == id2slot.end() )
	    {
	      const int n1 = id2slot.size(); 
	      id2slot[ id ] = n1;
	    }

	  // new faclvl?
	  std::string fl = "";
	  std::map<std::string,std::map<int, std::string> > ffll;
	  if ( facs.size() )
	    {
	      // at least one faclvl
	      std::map<std::string,int>::const_iterator ff = fac2slot.begin();
	      while ( ff != fac2slot.end() )
		{
		  fl += "_" + tok[ ff->second ] + "_" + dtok[ ff->second ];
		  
		  // count how many unique levels for this factor (in this file at least...)
		  if ( faclvlcnt[ tok[ ff->second ] ].find( dtok[ ff->second ] ) == faclvlcnt[ tok[ ff->second ] ].end() )
		    {
		      int n1 = faclvlcnt[ tok[ ff->second ] ].size();
		      faclvlcnt[ tok[ ff->second ] ][ dtok[ ff->second ] ] = n1;
		    }

		  // curr number of levels for this factor
		  std::map<int,std::string> dt;
		  int nl = faclvlcnt[ tok[ ff->second ] ][ dtok[ ff->second ] ];
		  dt[ nl ] = dtok[ ff->second ];
		  //std::cout << " nl = " << nl << " " << dtok[ ff->second ] << "\n";
		  ffll[ tok[ ff->second ] ] = dt; 
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
		
		
		// store actual value
		double val;
		if ( Helper::str2dbl( dtok[j] , &val ) ) 
		  D[ var2slot[ expand_vname ] ][ id2slot[ id ] ] = val;
		  
	      }
	  
	  // next data row
	}
      
      // all done with this file      
      IN1.close();

      logger << "  read " << file_ids.size() << " indivs, "
	     << file_bvars.size() << " base vars & "
	     << file_evars.size() << " expanded vars from "
	     << ff->first << "\n";
      
      ++ff;
    }


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
  
  std::map<std::string,int>::const_iterator ii = id2slot.begin();
  while ( ii != id2slot.end() )
    {
      // IDs
      ids.push_back( ii->first );

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
  // -------- write to binary file
  //

  logger << "\n";
  bfile_t bf( bfile );
  bf.write( ids, vars , var2group, basevar, faclvl , X );
  
  logger << "  ...done\n";

  //
  // --------- show manifest
  //

  manifest();
  
  
}


void gpa_t::read()
{
  bfile_t bf( bfile );
  bf.read( incvars , excvars, incnums, excnums,
	   incfacs, excfacs, incfaclvls, excfaclvls,
	   incgrps, excgrps, 
	   &ids, &vars , &var2group, &basevar, &faclvl , &X );
  logger << "  read " << ids.size() << " individuals and " << vars.size() << " variables from " << bfile << "\n";
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
  linmod_results_t results = lm.run( nreps );

  // outputs

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
	  
	  if ( results.emp[xvar][ var ] < pthresh && results.emp_corrected[xvar][ var ] < pthresh_adj )
	    {
	      shown_y = true; 
	      writer.level( var , "Y" );
	      writer.value( "B"  , results.beta[ xvar ][ var ] );
	      writer.value( "T"  , results.t[xvar][ var ] );
	      writer.value( "P" , results.emp[xvar][ var ] );
	      writer.value( "PADJ" , results.emp_corrected[xvar][ var ] );	  
	    }
	}
      if ( shown_y )
	writer.unlevel( "Y" );
            
    }
  
    writer.unlevel( "X" );

  
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

    // iterate over X

  for (int j=0; j<ivs.size(); j++)
    {
      
      const std::string & xvar = vars[ ivs[j] ];
      writer.level( xvar , "X" );

      // swap in X
      lm.set_IV( X.col( ivs[j] ) , xvar );
      
      // run
      linmod_results_t results = lm.run( nreps );

      // output
      bool shown_y = false;
      for (int y=0; y<dvs.size(); y++)
	{
	  const std::string & var = vars[ dvs[y] ];
	  
	  if ( results.emp[xvar][ var ] < pthresh && results.emp_corrected[xvar][ var ] < pthresh_adj )
	    {
	      shown_y = true; 
	      writer.level( var , "Y" );
	      writer.value( "B"  , results.beta[ xvar ][ var ] );
	      writer.value( "T"  , results.t[xvar][ var ] );
	      writer.value( "P" , results.emp[xvar][ var ] );
	      writer.value( "PADJ" , results.emp_corrected[xvar][ var ] );	  
	    }
	}
      if ( shown_y )
	writer.unlevel( "Y" );
            
    }  
  writer.unlevel( "X" );
  
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
		     const std::map<std::string,std::map<std::string,std::map<int,std::string> > > & faclvl ,
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
  std::map<std::string,std::map<std::string,std::map<int,std::string> > >::const_iterator ii = faclvl.begin();
  while ( ii != faclvl.end() )
    {
      const std::map<std::string,std::map<int,std::string> > & xx = ii->second;
      std::map<std::string,std::map<int,std::string> >::const_iterator ff = xx.begin();
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
      
      // and then value for each variable
      for (int i=0; i<nv; i++)
	{
	  std::map<int,std::string> vl;
	  std::map<std::string,std::map<std::string,std::map<int,std::string> > >::const_iterator vv = faclvl.find( vars[i] );
	  if ( vv != faclvl.end() )
	    {
	      std::map<std::string,std::map<int,std::string> >::const_iterator ll = vv->second.find( *ff );
	      if ( ll != vv->second.end() )
		vl = ll->second;
	    }
	  
	  bwrite( OUT1, vl.begin()->second ) ;
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
		    const std::vector<std::pair<int,int> > & incnums,
		    const std::vector<std::pair<int,int> > & excnums, 
		    const std::set<std::string> & incfacs,
		    const std::set<std::string> & excfacs,
		    const std::map<std::string,std::set<std::string> > & incfaclvls,
		    const std::map<std::string,std::set<std::string> > & excfaclvls, 		    
		    const std::set<std::string> & incgrps,
		    const std::set<std::string> & excgrps,		    
		    std::vector<std::string> * ids ,
		    std::vector<std::string> * vars ,
		    std::map<std::string,std::string> * var2group ,
		    std::map<std::string,std::string> * basevar,		    
		    std::map<std::string,std::map<std::string,std::map<int,std::string> > > * faclvl ,
		    Eigen::MatrixXd * X )
{
  
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
  
    
  // variables (which may be subsetted)
  std::vector<std::string> all_vars( nv );
  std::vector<std::string> all_groups( nv );
  std::vector<std::string> all_basevars( nv );
  
  std::map<std::string,std::map<std::string,std::map<int,std::string> > > all_faclvl;
  
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
	  if ( lvl != "" )
	    {
	      // i.e. keep ordering as per file
	      std::map<int,std::string> dt;
	      const int nl = (*faclvl)[ all_vars[j] ][ facs[k] ].size();
	      dt[ nl ] = lvl;
	      all_faclvl[ all_vars[j] ][ facs[k] ] = dt;
	    }
	}
    }

  
  //
  // restrict inputs?
  //
  
  std::vector<bool> readvar( nv , true );
  
  if ( incvars.size() || excvars.size()
       || incnums.size() || excnums.size()
       || incgrps.size() || excgrps.size() 
       || incfacs.size() || excfacs.size()
       || incfaclvls.size() || excfaclvls.size() )
    {

      const bool has_incvars = incvars.size();
      const bool has_incnums = incnums.size();
      const bool has_incgrps = incgrps.size();
      
      // if include lists specified, then set all to F initially
      if ( has_incvars || has_incnums || has_incgrps )
	{
	  readvar.clear();
	  readvar.resize( nv , false );
	  
	  for (int i=0; i<incnums.size(); i++  )
	    {
	      int s1 = incnums[i].first < incnums[i].second ? incnums[i].first : incnums[i].second ;
	      int s2 = incnums[i].first < incnums[i].second ? incnums[i].second : incnums[i].first ;
	      for (int j=s1; j<=s2; j++)
		readvar[j] = true;
	    }

	  // base-var inc?
	  if ( has_incvars ) 
	    for (int j=0; j<nv; j++)
	      {
		if ( incvars.find( all_basevars[j] ) != incvars.end() )
		  readvar[j] = true;
	      }

	  // var-group inc?
	  if ( has_incgrps )
	    for (int j=0; j<nv; j++)
              {
                if ( incgrps.find( all_groups[j] ) != incgrps.end() )
                  readvar[j] = true;
              }
	}

      // now do any exclusions
      const bool has_excvars = excvars.size();
      const bool has_excnums = excnums.size();
      const bool has_excgrps = excgrps.size();

      // if include lists specified, then set all to F initially
      if ( has_excvars || has_excnums || has_excgrps )
	{
	  
	  for (int i=0; i<excnums.size(); i++ )
	    {
	      int s1 = excnums[i].first < excnums[i].second ? excnums[i].first : excnums[i].second ;
	      int s2 = excnums[i].first < excnums[i].second ? excnums[i].second : excnums[i].first ;
	      for (int j=s1; j<=s2; j++)
		readvar[j] = false;
	    }

	  // base-var exc?
	  if ( has_excvars ) 
	    for (int j=0; j<nv; j++)
	      {
		if ( excvars.find( all_basevars[j] ) != excvars.end() )
		  readvar[j] = false;
	      }

	  // var-group exc?
	  if ( has_excgrps ) 
	    for (int j=0; j<nv; j++)
	      {
		if ( excgrps.find( all_groups[j] ) != excgrps.end() )
		  readvar[j] = false;
	      }
	  
	}


      //
      // finally, inc/exc based on facs or faclvls
      //

      // infacs

      if ( incfacs.size() || excfacs.size() )
	{
	  // exclude if does not contain all these factors
	  for (int j=0; j<nv; j++)
	    {
	      if ( readvar[j] )
		{
		  const std::map<std::string,std::map<int,std::string> > & fl = all_faclvl[ all_vars[j] ];

		  // incfacs || excfacs

		  if ( incfacs.size() )
		    {
		      if ( fl.size() != incfacs.size() ) readvar[j] = false;
		      else
			{
			  std::set<std::string>::const_iterator kk = incfacs.begin();
			  while ( kk != incfacs.end() )
			    {
			      if ( fl.find( *kk ) == fl.end() )
				{
				  readvar[j] = false;
				  break;
				}
			      ++kk;
			    }
			}
		    }

		  if ( excfacs.size() )
		    {
		      int c=0;
		      std::set<std::string>::const_iterator kk = excfacs.begin();
		      while ( kk != excfacs.end() )
			{
			  if ( fl.find( *kk ) != fl.end() ) ++c;
			  ++kk;
			}
		      if ( c == excfacs.size() && c == fl.size() ) readvar[j] = false;
		    }
		}
	    }
	}
    
      
      // infaclvls | excfaclvls
      
      if ( incfaclvls.size() || excfaclvls.size() )
	{

	  // only include/exclude based on factor levels
	  // irrespective of whether the full fac-set for that
	  // variable matches;
	  //  i.e. if faclvls=CH/CZ|FZ,B=SIGMA|GAMMA
	  //       then only pull CH==CZ whether fac is just CH, or CH+B, or CH+F, etc
	  //  do includes first, then excludes

	  // the above, has an AND logic between CH and B
	  // but an OR logic between CZ vs FZ, or SIGMA vs GAMMA

	  for (int j=0; j<nv; j++)
	    {
	      if ( readvar[j] )
		{
		  
		  const std::map<std::string,std::map<int,std::string> > & fl = all_faclvl[ all_vars[j] ];

		  // incfaclvls
		  if ( incfaclvls.size() )
		    {

		      bool match = true; 

		      // consider each factor;  okay for it not to exist, but if it does, then
		      // we must match one of the listed factors
		      std::map<std::string,std::set<std::string> >::const_iterator ll = incfaclvls.begin();
		      while ( ll != incfaclvls.end() )
			{
			  // we have this factor...
			  if ( fl.find( ll->first ) != fl.end() )
			    {
			      // ... do we have an acceptable level?
			      const std::string & lvl = fl.find( ll->first )->second.begin()->second;

			      if ( ll->second.find( lvl ) == ll->second.end() ) // no
				{
				  // no... then we are a no-go
				  match = false;
				  break;
				}
			    }
			  ++ll;
			}

		      if ( ! match )
			readvar[j] = false;

		    }



		  // excfaclvls
		  if ( excfaclvls.size() )
		    {
		      
		      bool match = true; 
		      
		      // consider each factor;  okay for it not to exist, but if it does, then
		      // we must *NOT* match one of the listed factors
		      std::map<std::string,std::set<std::string> >::const_iterator ll = excfaclvls.begin();
		      while ( ll != excfaclvls.end() )
			{
			  // we have this factor...
			  if ( fl.find( ll->first ) != fl.end() )
			    {
			      // ... do we have an unacceptable level?
			      const std::string & lvl = fl.find( ll->first )->second.begin()->second;

			      if ( ll->second.find( lvl ) != ll->second.end() ) // yes
				{
				  // then we are a no-go
				  match = false;
				  break;
				}
			    }
			  ++ll;
			}

		      if ( ! match )
			readvar[j] = false;

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
    
  // factors
  faclvl->clear();
  for (int k=0; k<nf; k++)
    {
      for (int j=0; j<nv; j++)
	if ( readvar[j] )
	  (*faclvl)[ all_vars[j] ][ facs[k] ] = all_faclvl[ all_vars[j] ][ facs[k] ];
    }
  
  
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
  std::map<std::string,std::map<std::string,std::map<int,std::string> > >::const_iterator ii = faclvl.begin();
  while ( ii != faclvl.end() )
    {
      const std::map<std::string,std::map<int,std::string> > & xx = ii->second;
      std::map<std::string,std::map<int,std::string> >::const_iterator ff = xx.begin();
      while ( ff != xx.end() )
	{
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
      std::cout << j+1 << "\t"
		<< vars[j] << "\t"
		<< ni - X.col(j).array().isNaN().sum() << "\t"
		<< var2group[ vars[j] ] << "\t"
		<< basevar[ vars[j] ] ;
		
      std::set<std::string>::const_iterator gg = allfacs.begin();
      while ( gg != allfacs.end() )
	{
	  std::string x = ".";
	  std::map<std::string,std::map<std::string,std::map<int,std::string> > >::const_iterator kk = faclvl.find( vars[j] );
	  if ( kk != faclvl.end() )
	    {
	      if ( kk->second.find( *gg ) != kk->second.end() )
		{
		  std::map<int,std::string> dt = kk->second.find( *gg )->second;
		  x = dt.begin()->second; // should always be exactly one value 
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

  const bool id_subsetting = rows.size() != 0 ;
  const bool col_subsetting = cols.size() != 0;
  
  // nothing to do
  if ( ! ( id_subsetting || col_subsetting ) ) return;
  
  int ni = X.rows();
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
  if ( col_subsetting )
    {            
      std::map<int,bool>::const_iterator cc = cols.begin();
      while ( cc != cols.end() )
	{
	  const bool pos_match = cc->second;      
	  const Eigen::VectorXd & v = X.col( cc->first );
	  
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
  
  logger << "  subsetted X from " << ni << " to " << X.rows() << " indivs\n";
}

// QC matrix
void gpa_t::qc( const double winsor )
{

  // 1) case-wise deletion
  std::vector<int> retained;
  const int ni = X.rows();
  for (int i=0;i<ni;i++)
    {
      const int num_missing = X.row(i).array().isNaN().sum();
      if ( num_missing == 0 ) retained.push_back( i );
    }

  if ( retained.size() < ni )
    {
      Eigen::MatrixXd X1 = X( retained, Eigen::all );
      X = X1;
      
      std::vector<std::string> ids2 = ids;
      ids.clear();
      for (int i=0;i<retained.size();i++)
	ids.push_back( ids2[retained[i]] );
     
      logger << "  case-wise deletion subsetted X from " << ni << " to " << X.rows() << " indivs\n";

    }

  
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
  if ( zeros.size() )
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
      std::map<std::string,std::map<std::string,std::map<int,std::string> > > faclvl2 = faclvl;
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

      logger << "  reduced data from " << nv << " to " << X.cols() << " vars\n"; 
    }
  
  
  
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


linmod_results_t linmod_t::run( int nreps )
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
      Eigen::VectorXd T = get_tstats( B.row(idx) , Yres , VX(idx,idx) , ni - nterms );

      Eigen::ArrayXd U = Eigen::ArrayXd::Ones( ny );
                  
      //logger << "  ";
      for (int r=0; r<nreps; r++)
	{
	  // logger << ".";
	  // if ( (r+1) % 10 == 0 ) logger << " ";
	  // if ( (r+1) % 50 == 0 ) logger << " " << r+1 << " perms\n" << ( r+1 == nreps ? "" : "  " ) ;

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

	  // accumulate
	  //double max_t = 0;
	  for (int y=0; y<ny; y++)
	    {
	      double abs_t = fabs( T_perm[y] );
	      // store point-wise emp-P now
	      if ( abs_t >= fabs( T[y] ) ) ++U[y];
	      // and track max per repl. for family-wise correction when done all X
	      if ( abs_t > max_t[r] ) max_t[r] = abs_t ;
	    }
	  
	  
	} // next replicate
      
      //
      // Get point-wise p-values
      //
      
      U /= (double)(nreps+1);

      Eigen::MatrixXd R( ny , 3 );

      R << B.row(idx).transpose() , T , U ;
      
      //
      // Point-wise results
      //
      
      for (int y=0; y<ny; y++)
	{      
	  results.beta[ xname[x] ][ vname[y] ] = R(y,0);
	  results.t[ xname[x] ][ vname[y] ] = R(y,1); 
	  results.emp[ xname[x] ][ vname[y] ] = R(y,2);
	  //	  results.emp_corrected[ xname[x] ][ vname[y] ] = R(y,3);
	}
            
      // next X
    }


  //
  // family-wise corrected P
  //
  
  for (int r=0; r<nreps; r++)
    for (int x=0; x<nx; x++)
      for (int y=0; y<ny; y++)
	if ( max_t[r] >= fabs( results.t[ xname[x] ][ vname[y] ] ) ) ++F(x,y);
  
  F /= (double)(nreps+1);
  
  for (int x=0; x<nx; x++)
    for (int y=0; y<ny; y++)
      results.emp_corrected[ xname[x] ][ vname[y] ] = F(x,y);
  
  return results;
}



Eigen::VectorXd linmod_t::get_tstats( const Eigen::VectorXd & B ,
				      const Eigen::MatrixXd & Yres ,
				      const double vx ,
				      const int denom )
{
  
  const int n = B.rows();
  
  Eigen::VectorXd T = Eigen::VectorXd::Zero( n );
  for (int i=0; i<n; i++)
    T[i] = Yres.col(i).transpose() * Yres.col(i); 
  
  // OR Eigen::VectorXd T = Yres.transpose() * Yres;
  
  for (int i=0; i<n; i++)
    T[i] = B[i] / sqrt( vx * T[i] / (double)denom ) ;
  
  // double sigmaSq =  ee[0]  / ( ni - 1 - 1 - nz ) ;
  // double se = sqrt( sigmaSq * vx );
  
  return T;
}


