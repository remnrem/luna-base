
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

#include "stats/cpt.h"
#include "eval.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "miscmath/miscmath.h"
#include "stats/eigen_ops.h"
#include "clocs/clocs.h"
#include "miscmath/crandom.h"

extern writer_t writer;
extern logger_t logger;

// notes: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4010955/
// TODO:   two-sided tests
// TODO:   more efficient adjacency detection

void cpt_wrapper( param_t & param )
{
  
  //
  // Similar input structure as for PSC; currently, 
  //   - specify a single IV 
  //   - F, CH or CH1+CH2 are the valid stratifiers currently
  //   - an ID-linked covariate table that includes a DV
  //   - a model specification (linear / logistic) 
  //   - adjacent frequencies defined in the obvious manner
  //   - adjacent spatially defined based on clocs, based on a distance threshold
  //

  //
  // Adjacency definition for EEG channel neighbours; likewise for frequency domains
  //

  double spatial_threshold = param.has( "th-spatial" ) ? param.requires_dbl( "th-spatial" ) : 0 ;
  
  double freq_threshold = param.has( "th-freq" ) ? param.requires_dbl( "th-freq" ) : 0 ; 


  //
  // Clustering threshold  (-ve means no clustering done)
  //
  
  double cl_threshold = param.has( "th-cluster" ) ? param.requires_dbl( "th-cluster" ) : -1 ; 

  bool one_sided_test = param.has( "1-sided" ) || param.has( "one-sided" );

  bool no_clustering = cl_threshold < 0 || param.has( "no-clustering" ); 


  
  //
  // Covariate and IV : assume a single file
  //
  
  // file containing IV and covariates
  std::string iv_file = param.requires( "iv-file" );
  
  // IV, e.g. iv=DIS
  std::string iv = param.requires( "iv" );

  // covariates, e.g. AGE,SEX [ optional ] 
  std::set<std::string> covars;
  if ( param.has( "covar" ) ) covars = param.strset( "covar" );
  
  //
  // Sleep metrics (DVs) : potentially many files
  //  (may be multiple files concatenated, but fine, and we can skip rows that are 'ID' ) 
  //

  if ( ! param.has( "dv-file" ) ) Helper::halt( "no dv-file option" );
  std::vector<std::string> dv_files = param.strvector( "dv-file" );

  // // e.g. 'COH' which may be stratified by F, CH, CH1 or CH2
  // std::set<std::string> dvars = param.strset( "dv" );
  
  
  //
  // IV outlier removal: assume that DV and covariates have already been assessed here
  //

  bool remove_Y_outliers = param.has( "th" );
  double outlier_threshold = remove_Y_outliers ? param.requires_dbl( "th" ) : 0 ; 

  bool winsorize = param.has( "winsor" );
  double winsor_pct = winsorize ? param.requires_dbl( "winsor" ) : 0 ; 

  if ( remove_Y_outliers && winsorize )
    Helper::halt( "cannot specify both th and winsor" );
  
  //
  // Number of replicates 
  //
  
  const int nreps = param.requires_int( "nreps" );


  //
  // Channel locations
  //

  std::string clocs_file = param.has( "clocs" ) ? param.value( "clocs" ) : "" ;
  
  
  //
  // Attach sleep metrics, and define adjacent points 
  //
  
  
  //
  // which variables to pull from these?
  //

  if ( ! param.has( "dv" ) ) Helper::halt( "no dv=var1,var2 specified" );
  std::set<std::string> dvars = param.strset( "dv" );

  //
  // which channels to include (if not all)
  //

  std::set<std::string> chs;
  if ( param.has( "ch" ) )    
    {
      chs = param.strset( "ch" );
      logger << "  expecting to retrain only " << chs.size() << " channels\n";
    }
  
  
  //
  // individual include/exclude lists? (e.g. skip outliers)
  // nb. can use @{exapansion} for file reading
  //

  std::set<std::string> id_includes, id_excludes;
  if ( param.has( "inc-ids" ) ) id_includes = param.strset( "inc-ids" );
  if ( param.has( "ex-ids" ) ) id_excludes = param.strset( "ex-ids" );
  
  if ( id_includes.size() > 0 && id_excludes.size() > 0 ) 
    Helper::halt( "cannot specify both inc-ids and ex-ids lists" );

  if ( id_includes.size() > 0 ) logger << "  read " << id_includes.size() << " IDs to include\n";
  if ( id_excludes.size() > 0 ) logger << "  read " << id_excludes.size() << " IDs to exclude\n";

  //
  // Require all obs, versus case-wise dropping
  //
  
  bool drop_incomplete_rows = ! param.has( "complete-obs" );
  
  //
  // make any variables absolute values
  //

  std::set<std::string> toabs;
  if ( param.has( "abs" ) ) toabs = param.strset( "abs" );

  //
  // log any variables? (dB)
  //

  std::set<std::string> tolog;
  if ( param.has( "dB" ) ) tolog = param.strset( "dB" );

  //
  // Outliers / normalization
  //

  std::vector<double> th;
  if ( param.has( "th" ) ) th = param.dblvector( "th" );
  

  //
  // any frequencies limits?
  //

  double flwr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 0 ;
  double fupr = param.has( "f-upr" ) ? param.requires_dbl( "f-upr" ) : 0 ;


  
  //
  // Attach covariates, define main IV: note, this sets the total
  // number of individuals
  //

  iv_file = Helper::expand( iv_file );
		  
  if ( ! Helper::fileExists( iv_file ) )
    Helper::halt( "could not load " + iv_file );
  
  std::ifstream IN1( iv_file.c_str() , std::ios::in );

  // parse headers
  std::vector<std::string> iv_header;
  std::string hline; 

  Helper::safe_getline( IN1 , hline );
  if ( IN1.eof() ) Helper::halt( "problem reading from " + iv_file );

  std::vector<std::string> tok = Helper::parse( hline , "\t" );
  std::map<std::string,int> iv_cols;
  for (int i=0; i<tok.size(); i++) iv_cols[ tok[i] ] = i;

  // check we have requested IV and covaraites
  if ( iv_cols.find( iv ) == iv_cols.end() ) Helper::halt( "could not find variable " + iv + " in " + iv_file );
  int iv_col = iv_cols[ iv ];

  // get ID col
  if ( iv_cols.find( "ID" ) == iv_cols.end() ) Helper::halt( "could not find ID column (case-sensitive) in " + iv_file );
  int id_col = iv_cols[ "ID" ];

  // get covar columns
  std::vector<int> covar_col;
  std::vector<std::string> covar_label;
  std::set<std::string>::const_iterator cc = covars.begin();
  while ( cc != covars.end() )
    {
      if ( iv_cols.find( *cc ) == iv_cols.end() ) Helper::halt( "could not find variable " + *cc + " in " + iv_file );
      covar_col.push_back( iv_cols[ *cc ] );
      covar_label.push_back( *cc );
      ++cc;
    }

  const int iv_coln = tok.size();

  
  //
  // now read rest of data : one row per individual
  //

  int row_cnt = 0;

  std::vector<std::string> ids;
  std::map<std::string,int> ids_map; // ID --> row
  std::map<std::string,std::vector<double> > ivdata;
  while ( ! IN1.eof() )
    {
      std::string dline;
      Helper::safe_getline( IN1 , dline );      
      if ( IN1.eof() ) break;
      if ( dline == "" ) continue;
      if ( dline[0] == '%' || dline[0] == '#' ) continue; 

      ++row_cnt;

      std::vector<std::string> tok = Helper::parse( dline , "\t" );
      if ( tok.size() != iv_coln ) Helper::halt( "bad number of columns in " + iv_file + "\n" + dline );

      // add ID
      const std::string this_id = tok[ id_col ] ;
      if ( id_excludes.size() != 0 && id_excludes.find( this_id ) != id_excludes.end() ) continue;
      if ( id_includes.size() != 0 && id_includes.find( this_id ) == id_includes.end() ) continue;

      const int this_idn = ids_map.size();
      ids.push_back( tok[ id_col ] );
      ids_map[ tok[ id_col ] ] = this_idn;

      // add IV
      double iv_num;
      if ( ! Helper::str2dbl( tok[ iv_col ] , &iv_num ) )
	Helper::halt( "problem with numeric value" + tok[ iv_col ] );
      ivdata[ iv ].push_back( iv_num );

      // add covariates
      for (int c=0; c<covar_col.size(); c++)
	{
	  double cov_num;
	  if ( ! Helper::str2dbl( tok[ covar_col[c] ] , &cov_num ) )
	    Helper::halt( "problem with numeric value" + tok[ covar_col[c] ] );
	  ivdata[ covar_label[c] ].push_back( cov_num );
	}

      // next row
      
    }
  IN1.close();

  logger << "  read " << ids.size() << " people from " << iv_file << " (of total " << row_cnt << " data rows)\n";

  //
  // Number of people in the IV file (might be > than for the sleep metrics) 
  //

  const int ni_ivfile = ids.size();


  //
  // Read in sleep metrics : multiple metrics, multiple stratifications (F and/or CH, or CH1+CH2)
  // Same approach as for --psc command
  //

  // can stratify by channel, frequency (columns)
  // expect input 
  // can also have coherence (CH1/CH2) values too
  // expect file to contain multiple individuals
  
  // in output, one row per individual; columns are CH x F (xCH1/CH2)
  
  // id -> ch -> f -> var -> value



  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > i2c2f2v;
  
  for (int i=0; i<dv_files.size(); i++)
    {
      std::string infile = Helper::expand( dv_files[i] );

      logger << "  reading metrics from " << infile << "\n";

      if ( ! Helper::fileExists( infile ) ) 
	Helper::halt( "could not find " + infile );

      std::ifstream IN1( infile.c_str() , std::ios::in );

      // header row

      std::string hline;
      
      Helper::safe_getline( IN1 , hline );

      if ( IN1.eof() )
	{
	  IN1.close();
	  continue;
	}

      //
      // expects a tab-delimited file
      //
      
      std::vector<std::string> tok = Helper::parse( hline , "\t" );
      std::set<std::string> cols;
      for (int i=0; i<tok.size(); i++) cols.insert( tok[i] );

      //
      // Check for key column names
      //

      if ( cols.find( "ID" ) == cols.end() ) Helper::halt( "no ID column in " + infile );
      bool ch1 = cols.find( "CH" ) != cols.end();
      bool ch2 = cols.find( "CH1" ) != cols.end() && cols.find( "CH2" ) != cols.end();
      bool frq = cols.find( "F" ) != cols.end() ;
      
      // channel and frequency are optional

      int id_slot = -1;      
      int ch_slot = -1;
      int ch1_slot = -1;
      int ch2_slot = -1;
      int f_slot = -1;      
      
      std::map<int,std::string> slot2var;
      for (int i=0;i<tok.size();i++)
	{
	  if ( tok[i] == "ID" ) id_slot = i;
	  if ( tok[i] == "F" ) f_slot = i;
	  if ( tok[i] == "CH" ) ch_slot = i;
	  if ( tok[i] == "CH1" ) ch1_slot = i;
	  if ( tok[i] == "CH2" ) ch2_slot = i;	  
	  if ( dvars.find( tok[i] ) != dvars.end() ) slot2var[i] = tok[i];
	}
      
      if ( slot2var.size() == 0 ) 
	Helper::halt( "no variables dvars=<...> in " + infile );

      const int ncols = tok.size();
      
      //
      // all set, now start reading
      //

      while ( ! IN1.eof() )
	{
	  // looking for ID F CH         --> PSD 
	  //  OR         ID F CH1 CH2    --> LCOH (default)  
	  //  OR         ID              --> HYPNO

	  std::string line;

	  Helper::safe_getline( IN1 , line );

	  std::vector<std::string> tok = Helper::parse( line , "\t" );

	  if ( IN1.eof() || tok.size() == 0 ) continue;
	  
	  if ( tok.size() != ncols ) Helper::halt( "incorrect number of columns in " + infile );
	  
	  std::string id = tok[ id_slot ] ;
	    
	  // skip if person was not included in the IV dataset
	  if ( ids_map.find( id ) == ids_map.end() ) continue;
	  	  
	  // channels requested (if channels present in this file)
          if ( chs.size() != 0 && ( ch1 || ch2 ) )
            {
              bool okay = true;
              if ( ch_slot == -1 )
		okay = chs.find( tok[ ch1_slot ] ) != chs.end() && chs.find( tok[ ch2_slot ] ) != chs.end() ;
              else
                okay = chs.find( tok[ ch_slot ] ) != chs.end() ;
              // skip if channel(s) not found
	      
              if ( ! okay ) continue;
            }

	  std::string ch = "";
	  
	  if ( ch_slot != -1 )
	    ch = tok[ ch_slot ];
	  else if ( ch1_slot != -1 && ch2_slot != -1 )
	    ch = tok[ ch1_slot ] + "." + tok[ ch2_slot ];
	  else
	    ch = "-"; // or '-' if no CH vars

	  // store as string and also numeric (for output) [ or '0' if no F variable ] 
	  std::string f = f_slot != -1 ? tok[ f_slot ] : "0" ; 

	  // filter on frequency, if present
	  if ( f_slot != -1 && ( flwr > 0 || fupr > 0 ) ) 
	    {
	      double fn;
	      if ( ! Helper::str2dbl( f , &fn ) )
		Helper::halt( "problem with frequency value: " + f );
	      // skip this frequency?
	      if ( flwr > 0 && fn < flwr ) continue;
	      if ( fupr > 0 && fn > fupr ) continue;
	    }

	  
	  // Get values into the map 
	  std::map<int,std::string>::const_iterator ii = slot2var.begin();
	  while ( ii != slot2var.end() )
	    {
	      
	      // is this missing? 
	      if ( Helper::iequals( tok[ ii->first ] , "NA" ) 
		   || Helper::iequals( tok[ ii->first ] , "nan" ) 
		   || Helper::iequals( tok[ ii->first ] , "inf" ) ) 
		{
		  ++ii;
		  continue;
		}
		   
	      
	      double x;
	      if ( ! Helper::str2dbl( tok[ ii->first ] , &x ) ) 
		Helper::halt( "bad value in " + infile + "\n" 
			      + ii->second + " --> [" + tok[ ii->first ] + "]" );

	      // make absolute?
	      if ( toabs.find( ii->second ) != toabs.end() ) x = fabs( x );
	      
	      // take log?
	      if ( tolog.find( ii->second ) != tolog.end() ) x = 10 * log10( x ) ;

	      // save
	      i2c2f2v[ id ][ ch ][ f ][ ii->second ] = x;
	      ++ii;
	    }
	  
	}
      IN1.close();
      
      // next input file
    }


  //
  // Construct data matrix
  //

  logger << "  converting input files to a single matrix\n";
    
  std::vector<std::string> vname;  
  std::map<std::string,std::map<std::string,std::map<std::string,int> > > slot;
  std::map<std::string,std::string> col2ch, col2var;
  std::map<std::string,std::string> col2ch1, col2ch2; 
  std::map<std::string,double> col2f;
  std::set<std::string> rows, cols;
  std::set<std::string> drop_indivs;
  std::vector<std::string> id;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > >::const_iterator ii1 = i2c2f2v.begin();
  while ( ii1 != i2c2f2v.end() )
    {
      // unique individual IDs
      rows.insert( ii1->first );
      id.push_back( ii1->first );

      // make col names
      std::map<std::string,std::map<std::string,std::map<std::string,double> > >::const_iterator ii2 = ii1->second.begin();
      while ( ii2 != ii1->second.end() )
	{
	  std::map<std::string,std::map<std::string,double> >::const_iterator ii3 = ii2->second.begin();
	  while ( ii3 != ii2->second.end() )
	    {
	      std::map<std::string,double>::const_iterator ii4 = ii3->second.begin();
	      while ( ii4 != ii3->second.end() )
		{
		  std::string col_name = ii2->first + "~" + ii3->first + "~" + ii4->first ;

		  col2ch[ col_name ] = ii2->first;
		  std::vector<std::string> ctok = Helper::parse( ii2->first , "." );
		  col2ch1[ col_name ] = ctok[0] ;
		  col2ch2[ col_name ] = ctok.size() == 2 ? ctok[1] : "." ; 
		  
		  double ff;
		  if ( ! Helper::str2dbl( ii3->first , &ff ) ) Helper::halt( "problem with F non-numeric value" );
		  col2f[ col_name ] = ff;
		  col2var[ col_name ] = ii4->first;

		  if ( cols.find( col_name ) == cols.end() )
		    {
		      cols.insert( col_name );
		      vname.push_back( col_name );
		      slot[ ii2->first ][ ii3->first ][ ii4->first ] = cols.size() - 1 ;
		    }
		  ++ii4;
		}
	      ++ii3;
	    }
	  ++ii2;
	}
      ++ii1;
    }

  logger << "  found " << rows.size() << " rows (individuals) and " << cols.size() << " columns (features)\n";

  if ( rows.size() == 0 || cols.size() == 0 ) 
    return;


  //
  // Find individuals to drop (prior to populating the matrix
  //
  
  std::map<std::string,std::map<std::string,std::map<std::string,int> > >::const_iterator ss1 = slot.begin();
  while ( ss1 != slot.end() )
    {
      std::map<std::string,std::map<std::string,int> >::const_iterator ss2 = ss1->second.begin();
      while ( ss2 != ss1->second.end() )
	{
	  std::map<std::string,int>::const_iterator ss3 = ss2->second.begin();
	  while ( ss3 != ss2->second.end() )
	    {

	      // get all individuals for this column?	      
	      std::set<std::string>::const_iterator ii = rows.begin();
	      while ( ii != rows.end() )
		{
		 
		  // find channel?
		  const std::map<std::string,std::map<std::string,std::map<std::string,double> > > & dat = i2c2f2v.find( *ii )->second;
		  if ( dat.find( ss1->first ) == dat.end() ) 
		    {
		      if ( drop_incomplete_rows )
			{
			  drop_indivs.insert( *ii );
			  ++ii; continue; 
			}
		      Helper::halt( "no channel " + ss1->first + " for individual " + *ii );
		    }

		  // find frequency?
		  const std::map<std::string,std::map<std::string,double> > & dat2 = dat.find( ss1->first )->second;
		  if ( dat2.find( ss2->first ) == dat2.end() ) 
		    {
		      if ( drop_incomplete_rows )
			{
			  drop_indivs.insert( *ii );
			  ++ii; continue; 
			}
		      Helper::halt( "no frequency " + ss2->first + " for individual " + *ii );
		    }

		  // find variable?
		  const std::map<std::string,double> & dat3 = dat2.find( ss2->first )->second;
		  if ( dat3.find( ss3->first ) == dat3.end() ) 
		    {
		      if ( drop_incomplete_rows )
			{
			  drop_indivs.insert( *ii );
			  ++ii; continue; 
			}
		      Helper::halt( "no variable " + ss3->first + " for individual " + *ii );
		    }
		  
		  ++ii;
		}
	      
	      ++ss3;
	    }
	  ++ss2;
	}
      ++ss1;
    }


  if ( drop_incomplete_rows )
    logger << "  identified " << drop_indivs.size() << " of " << rows.size() << " individuals with at least some missing data\n";
  
  if ( rows.size() - drop_indivs.size() <= 2 ) 
    Helper::halt( "not enough observationns for CPT analysis" );

  // clean up rows
  if ( drop_indivs.size() > 0 ) 
    {
      std::set<std::string>::const_iterator dd = drop_indivs.begin();
      while ( dd != drop_indivs.end() )
	{
	  rows.erase( rows.find( *dd ) );
	  ++dd;
	}
      
      // re-populate ID list
      id.clear();
      dd = rows.begin();
      while ( dd != rows.end() )
	{
	  id.push_back( *dd );
	  ++dd;
	}
    }


  //
  // Populate matrix
  //
  
  Eigen::MatrixXd Y( rows.size() , cols.size() );
  
  //std::map<std::string,std::map<std::string,std::map<std::string,int> > >::const_iterator ss1 (above)
  ss1 = slot.begin();
  while ( ss1 != slot.end() )
    {
      std::map<std::string,std::map<std::string,int> >::const_iterator ss2 = ss1->second.begin();
      while ( ss2 != ss1->second.end() )
	{
	  std::map<std::string,int>::const_iterator ss3 = ss2->second.begin();
	  while ( ss3 != ss2->second.end() )
	    {

	      // get all individuals for this column?
	      int row = 0;
	      std::set<std::string>::const_iterator ii = rows.begin();
	      while ( ii != rows.end() )
		{
		  
		  // find channel
		  const std::map<std::string,std::map<std::string,std::map<std::string,double> > > & dat = i2c2f2v.find( *ii )->second;
		  
		  // find frequency
		  const std::map<std::string,std::map<std::string,double> > & dat2 = dat.find( ss1->first )->second;
		  
		  // find variable
		  const std::map<std::string,double> & dat3 = dat2.find( ss2->first )->second;
		  
		  // these have all been checked now: should be okay to add to the store
		  Y( row, ss3->second ) = dat3.find( ss3->first)->second;
		  ++row;
		  ++ii;
		}

	      ++ss3;
	    }
	  ++ss2;
	}
      ++ss1;
    }

  logger << "  finished making regular data matrix on " << rows.size()  << " individuals\n";



  //
  // free up main memory store
  //

  i2c2f2v.clear();
  

  //
  // Check for invariant columns
  //

  const double EPS = 1e-6;
  int N = Y.rows();

  Eigen::Array<double, 1, Eigen::Dynamic> means = Y.colwise().mean();
  Eigen::Array<double, 1, Eigen::Dynamic> sds = ((Y.array().rowwise() - means ).square().colwise().sum()/(N-1)).sqrt();
  for (int i=0; i < sds.size(); i++)
    if ( sds[i] < EPS ) 
      Helper::halt( "at least one invariant column in input; first = : " + vname[i] );

  int ni = rows.size();
  
  //
  // Case-wise outlier removal?
  //


  if ( remove_Y_outliers )
    {
      
      const int nv = cols.size();
      
      std::vector<bool> inc( ni , true );
      
      for (int t=0; t<th.size(); t++)
	{
	  // copy current:
	  std::vector<bool> prior = inc;
	  
	  for (int j=0; j<nv; j++)
	    {
	      // urgh, for now copy vector from Eigen to muse outliers() function... 
	      // this is not a big timesink in the flow of the PSC option, so 
	      // should not matter... we will clean up later converting all matrix/;vector helpers
	      // to assume Eigen objects
	      
	      // this sets 'inc' values to missing, but uses the same prior for all channels
	      
	      std::vector<double> tmp(ni);
	      Eigen::VectorXd::Map( &tmp[0], ni ) = Y.col(j);
	      
	      int removed = MiscMath::outliers( &tmp , th[t] , &inc , &prior);
	      
	      //	  logger << "  removing " << removed << " var " << j << " round " << t << "\n";
	    }
	  
	}
  
      //
      // Remove rows from input U
      //
      
      Eigen::MatrixXd Y2 = Y;      
      std::vector<std::string> id2 = id;
      
      Y.resize(0,0);
      id.clear();
      
      ni = 0;
      for (int i=0;i<inc.size();i++)
	if ( inc[i] ) ++ni;
      
      logger << "  after outlier removal, " << ni << " individuals remaining\n";
      
      Y.resize( ni , nv );
      id.resize( ni );
      int r = 0;
      for (int i=0;i<Y2.rows();i++)
	if ( inc[i] )
	  {
	    for (int j=0;j<nv;j++) 
	      Y(r,j) = Y2(i,j);
	    id[r] = id2[i];
	    ++r;
	  }
    }


  //
  // Winsorization?
  //

  if ( winsorize )
    {
      
      if ( ! eigen_ops::robust_scale( Y , winsor_pct ) )
	Helper::halt( "one or more features with no variability... quitting\n(better error message hopefully forthcoming..." );      

    }

  
  
  //
  // Create X and Z matrices that match Y
  //

  const int nz = covar_label.size();
  
  Eigen::VectorXd X = Eigen::VectorXd::Zero( ni );
  Eigen::MatrixXd Z = Eigen::MatrixXd::Zero( ni , nz );

  for (int i=0; i<ni; i++)
    {
      const int idx = ids_map[ id[i] ];
      
      // main IV
      X[i] = ivdata[ iv ][ idx ];
	
      // Covariates
      for (int c=0; c<covar_label.size(); c++)
	Z(i,c) = ivdata[ covar_label[c] ][ idx ];
      
    }
  

  logger << "  final datasets contains " << Y.cols() << " DVs on " << X.rows() << " individuals, "
	 << X.cols() << " primary IV(s), and " << Z.cols() << " covariate(s)\n";
  
  
  //
  // Perform CPT 
  //
  
  cpt_t cpt( Y , X , Z ) ;


  //
  // Channel locations
  //
  
  clocs_t clocs;
  if ( clocs_file != "" )
    clocs.load_cart( clocs_file );

  //
  // Verbose output for adjacencies
  //

  const bool verbose = param.has( "dump-adj" ) || param.has( "verbose" );
  
  //
  // Define adjacencies
  //

  if ( ! no_clustering ) 
    {
      logger << "  defining adjacent variables...\n";
      
      cpt.calc_adjacencies( vname , col2var , col2f , col2ch1 , col2ch2 ,
			    freq_threshold ,
			    clocs_file == "" ? NULL : &clocs ,
			    spatial_threshold , 
			    verbose ) ;
    }
  else
    {
      // but need to store these
      cpt.vname = vname;
    }
  
  //
  // Run permutations
  //

  logger << "  running permutations, assuming a " << ( one_sided_test ? "one" : "two" ) << "-sided test...\n";
  
  cpt_results_t results = cpt.run( nreps , cl_threshold , ! one_sided_test , verbose );

  logger << "  all done.\n";
  
  //
  // Report results
  //

  for (int y=0; y<vname.size(); y++)
    {
      const std::string & var =  vname[y] ;
      writer.level( var , globals::var_strat );
      writer.value( "B"  , results.beta[ var ] );
      writer.value( "T"  , results.t[ var ] );
      writer.value( "PU" , results.emp[ var ] );
      writer.value( "PC" , results.emp_corrected[ var ] );      
      writer.value( "CLST" , results.inclst[ var ] ); // 0 if not in a cluster

      // extra information on CH, F (or CH1, CH2)
      if ( col2ch2[ var ] != "." )
	{
	  writer.value( "CH1" ,  col2ch1[ var ] );
	  writer.value( "CH2" ,  col2ch2[ var ] );		  
	}
      else
	{
	  writer.value( "CH" ,  col2ch[ var ] );
	}

      if ( col2f[ var ] > 0 )
	writer.value( "F" ,  col2f[ var ] );
      
    }
  writer.unlevel( globals::var_strat );

  if ( ! no_clustering ) 
    {
      // clusters
      int cln = 0;
      std::map<std::string,double>::const_iterator qq = results.cluster_emp.begin();
      while ( qq != results.cluster_emp.end() )
	{
	  const std::set<std::string> & members = results.cluster_members.find( qq->first )->second;
	  
	  writer.level( ++cln , globals::cluster_strat );
	  
	  writer.value( "SEED" , qq->first );
	  writer.value( "P" ,  qq->second );
	  writer.value( "N" , (int)members.size() );
	  
	  // members
	  int memn = 0;
	  std::set<std::string>::const_iterator mm = members.begin();
	  while ( mm != members.end() )
	    {
	      writer.level( ++memn , "M" );
	      writer.value( "VAR" , *mm );
	      ++mm;
	    }
	  writer.unlevel( "M" );
	  ++qq;
	}
      writer.unlevel( globals::cluster_strat );
    }


  //
  // All done
  //

  

}
  


void cpt_t::set_DV( const Eigen::MatrixXd & Y_ )
{
  
  Y = Y_;

  if ( ni != 0 && ni != Y.rows() )
    Helper::halt( "unequal number of rows" );
  else
    ni = Y.rows();
  ny = Y.cols();
  
}

void cpt_t::set_IV( const Eigen::VectorXd & X_ )
{
  
  X = X_;
  
  if ( ni != 0 && ni != X.rows() )
    Helper::halt( "unequal number of rows" );
  else
    ni = X.rows();
  
}

void cpt_t::set_Z( const Eigen::MatrixXd & Z_ )
{
  Z = Z_;
  
  if ( ni != 0 && ni != Z.rows() )
    Helper::halt( "unequal number of rows" );
  else
    ni = Z.rows();
  
  nz = Z.cols();
  
}


void cpt_t::calc_adjacencies( const std::vector<std::string> & vname_ , 
			      const std::map<std::string,std::string> & col2var,
			      const std::map<std::string,double> & col2f,
			      const std::map<std::string,std::string> & col2ch1,
			      const std::map<std::string,std::string> & col2ch2,
			      double fth ,
			      clocs_t * clocs ,
			      double sth , 
			      bool dump_adj )
{

  // store these for output
  vname = vname_;
  
  // std::map<int,std::set<int> > adjacencies;
  adjacencies.clear();

  // rules:
  //  - assume uniform grid of frequencies for each channel / channel pair
  //  - different root variables (e.g. PSD and COH) cannot be adjacent (i.e. clusters only defined temporally and spatially within root variable)
  //     ( we still want the multiple test correction to potentially operate over multiple variables though, thus the inclusion of multiple root vars)
  //  - variables have 0 or 1 associated frequency
  //  - variables either have 0, 1 or 2 associated channels; adjacencies are only defined within these groups
  //  -  for 2 channel variables, adjacncy means both pairs are adjacent: but can be flipped,
  //  -   i.e. for pairs A1-B1   and A2-B2       adj if : A1 adj to A2 AND B1 adj to B2 ... OR  A1 adj B2 AND A2 adj B1 
  //  -  two variables are adjacent only if both channel AND frequency terms are adjacent (assuming those are defined) 

  // use 'vname' as the key for look up 

  const int nv = vname.size();
  if ( nv != Y.cols() && Y.cols() > 0 )
    Helper::halt( "variable definitions do not match Y matrix # of cols" );

  std::map<std::string,int> var2col;
  for (int i=0; i<nv; i++) var2col[ vname[i] ] = i;

  std::vector<std::string> var( nv );
  std::vector<double> freq( nv , -1 ); // -1 means 'not defined'
  std::vector<std::string> ch1( nv , "." ); // '.' means 'not defined'
  std::vector<std::string> ch2( nv , "." ); // '.' means 'not defined'

  // to speed up seaerch, store all frequencyes, and so we only look at other variables where we
  // know the frequency may match
  // use string for F key (numerical precision)
  std::map<std::string,std::set<int> > f2slot;
  std::map<std::string,double> f2num;

  std::set<std::string> chs;

  for (int i=0; i<nv; i++)
    {      
      freq[i] = col2f.find( vname[i] )->second <= 0 ? -1 : col2f.find( vname[i] )->second;
      ch1[i] = col2ch1.find( vname[i] )->second ;
      ch2[i] = col2ch2.find( vname[i] )->second ;      
      var[i] = col2var.find( vname[i] )->second ;
      
      // for speeding up search below
      std::string s = Helper::dbl2str( freq[i] );
      f2slot[ s ].insert( i );
      f2num[ s ] = freq[i];

      // track channels
      chs.insert( ch1[i] );
      chs.insert( ch2[i] ); 
    }
 

  // pre-calculate distance matrix
  std::map<std::string,std::map<std::string,double> > dist_matrix;
  if ( clocs != NULL )
    {
      std::set<std::string>::const_iterator cc1 = chs.begin();
      while ( cc1 != chs.end() )
	{
	  if ( *cc1 != "." )
	    {
	      std::set<std::string>::const_iterator cc2 = chs.begin();
	      while ( cc2 != chs.end() )
		{
		  if ( *cc2 != "." ) 
		    dist_matrix[ *cc1 ][ *cc2 ] = clocs->distance( *cc1 , *cc2 , 2 );
		  ++cc2;
		}
	    }
	  ++cc1;
	}
    }
  

  //
  // main pairwise comparisons of all variables
  //
  
  for (int i=0; i<nv; i++)
    {
      
      std::map<std::string,double>::const_iterator ff = f2num.begin();
      
      while ( ff != f2num.end() )
	{
	  
	  // freq adjacent?
	  bool freq_adjacent = false;
	  if ( freq[i] <= 0 && ff->second <= 0 ) freq_adjacent = true;  // if F not defined, then they are 'adjacent'
	  else if ( freq[i] <= 0 != ff->second <= 0 ) freq_adjacent = false; // but both variables must be freq agnostic
	  else if ( fabs( freq[i]  - ff->second ) <= fth ) freq_adjacent = true;
	  	  	  
	  if ( freq_adjacent )
	    {

	      // we need to check the following variables
	      
	      const std::set<int> & tocheck = f2slot[ ff->first ];
	      
	      // rather than all pairs
	      //for (int j=i+1;j<nv; j++)
	      //

	      std::set<int>::const_iterator jj = tocheck.begin();

	      while ( jj != tocheck.end() )
		{

		  // same variable?
		  if ( i == *jj ) { ++jj; continue; } 
		  
		  // different root variable? 
		  if ( var[i] != var[*jj] ) { ++jj; continue; } 

		  // we know this will be freq adjacent 
		  // // freq adjacent?
		  // bool freq_adjacent = false;
		  // if      ( freq[i] <= 0 && freq[j] <= 0 ) freq_adjacent = true;  // if F not defined, then they are 'adjacent'
		  // else if ( freq[i] <= 0 != freq[j] <= 0 ) freq_adjacent = false; // but both variables must be freq agnostic
		  // else if ( fabs( freq[i]  - freq[j] ) <= fth ) freq_adjacent = true;  	
	
		  // spatial match
		  bool spatial_adjacent = false;
		  int ci = ch1[i] == "." ? 0 : ( ch2[i] == "." ? 1 : 2 );
		  int cj = ch1[*jj] == "." ? 0 : ( ch2[*jj] == "." ? 1 : 2 );

		  if ( ci != cj ) spatial_adjacent = false; // different domain
 		  else if ( ci == 0 ) spatial_adjacent = true; // no spatial info
		  else
		    {
		      if ( clocs != NULL )
			{		

			  // single channels matching
			  if ( ci == 1 )
			    {			      			      
			      if ( dist_matrix[ ch1[i] ][ ch1[ *jj ] ] < sth )
				spatial_adjacent = true;
			    }
			  else // or this is for a pair of channels
			    {			      
			      // scenario 1:
			      double d11 = dist_matrix[ ch1[i] ][ ch1[*jj] ];
			      double d22 = dist_matrix[ ch2[i] ][ ch2[*jj] ];
			      
			      // scenario 2: (flipped)
			      double d12 = dist_matrix[ ch1[i] ][ ch2[*jj] ];
			      double d21 = dist_matrix[ ch2[i] ][ ch1[*jj] ];

			      // i.e. for A1-B1 and A2-B2, then
			      //  either A1/A2 AND B1/B2 must match
			      //      OR A1/B2 AND A2/B1

			      if ( ( d11 < sth && d22 < sth ) || ( d12 < sth || d21 < sth ) ) 
				spatial_adjacent = true;
			    }
			}
		      else
			Helper::halt( "no clocs supplied" );
		    }
		  
		  // so... are these two variables adjacent or no?
		  
		  if ( spatial_adjacent && freq_adjacent )
		    {
		      adjacencies[ i ].insert( *jj );
		      adjacencies[ *jj ].insert( i );	    
		    }
		  
		  ++jj;
		} // next set of 'j' vars
	    } // end 'if freq adj'
	  
	  ++ff; 
	} // next freq bin 
    } // outer loop, next 'i' variable
  

  //
  // dump adjacencies
  //
  
  double mean_adjn = 0;
  int cnt = 0;

  std::map<int,std::set<int> >::const_iterator ss = adjacencies.begin();
  while ( ss != adjacencies.end() )
    {
      const std::set<int> & adj = ss->second;
      if ( dump_adj ) 
	std::cout << vname[ ss->first ] << "\t" 
		  << adj.size() << "\n";
      
      mean_adjn += adj.size();

      if ( dump_adj )
	{
	  std::set<int>::const_iterator qq = adj.begin();
	  while ( qq != adj.end() )
	    {
	      std::cout << " -> " << vname[ *qq ] ; 
	      ++qq;
	    }
	}
      if ( dump_adj ) std::cout << "\n";
      ++ss;
    }
  
  logger << "  on average, each variable has " << mean_adjn / (double)vname.size() << " adjacencies\n"; 

  logger << "  " << vname.size() - adjacencies.size() << " variable(s) have no adjacencies\n";
  
  
  // all done
}


cpt_results_t cpt_t::run( int nreps , double cl_threshold , bool two_sided_test , bool verbose )
{

  cpt_results_t results;

  ni = Y.rows();
  ny = Y.cols();
  //nx = X.cols(); // has to be 1 currently
  nz = Z.cols();

  if ( X.cols() != 1 ) Helper::halt( "cpt_t not set up yet for multiple X" );


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
  // M is full model: intercept + nuissance + design (1 term currently)
  //
  
  Eigen::MatrixXd MM( ni , 1 + nz + 1 );
  MM << ZZ , X ; 
  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cqrM( MM );
  Eigen::MatrixXd Minv = cqrM.pseudoInverse();
  Eigen::MatrixXd Hm = MM * Minv;
  Eigen::MatrixXd Rm = Eigen::MatrixXd::Identity( Hm.rows() , Hm.cols() ) - Hm ;
  
  //
  // Get get obserevd statistics
  //

  Eigen::MatrixXd YZres = Rz * Y;
  Eigen::MatrixXd B = Minv * YZres;
  Eigen::MatrixXd Yres = Rm * YZres;
  Eigen::MatrixXd VX = ( MM.transpose() * MM ).inverse() ;
  const int nterms = 1 + nz + 1 ; // intercept + covariates + IV 
  const int idx = nterms - 1;
  Eigen::VectorXd T = get_tstats( B.row(idx) , Yres , VX(idx,idx) , ni - nterms );


  //
  // Get clusters
  //

  cpt_clusters_t clusters( T , cl_threshold , adjacencies , two_sided_test , verbose , &vname ); 
  
  logger << "  found " << clusters.clusters.size()
	 << " clusters, maximum statistic is "
	 << clusters.max_stat << "\n";
  
  //
  // Initiate permutation counters
  //

  // uncorrected (initiate at 1 to include the observed data)
  Eigen::ArrayXd U = Eigen::ArrayXd::Ones( ny );  

  // family-wise
  Eigen::ArrayXd F = Eigen::ArrayXd::Ones( ny );  

  logger << "  ";
  for (int r=0; r<nreps; r++)
    {
      logger << ".";
      if ( (r+1) % 10 == 0 ) logger << " ";
      if ( (r+1) % 50 == 0 ) logger << " " << r+1 << " perms\n" << ( r+1 == nreps ? "" : "  " ) ;

      // shuffle
      std::vector<int> pord( ni );
      CRandom::random_draw( pord );

      // permutation matrix
      Eigen::MatrixXd P = Eigen::MatrixXd::Zero( ni , ni );
      for (int i=0; i<ni; i++) P(i,pord[i]) = 1;

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
      double max_t = 0;
      for (int y=0; y<ny; y++)
	{
	  double abs_t = fabs( T_perm[y] );
	  if ( abs_t >= fabs( T[y] ) ) ++U[y];
	  if ( abs_t > max_t ) max_t = abs_t ;
	}
      
      for (int y=0; y<ny; y++)
	if ( max_t >= fabs( T[y] ) ) ++F[y];


      //
      // clustering
      //
      
      cpt_clusters_t perm_clusters( T_perm , cl_threshold , adjacencies , two_sided_test );
      
      clusters.update( perm_clusters.max_stat );

    }


  //
  // Get point-wise p-values
  //

  U /= (double)(nreps+1);
  F /= (double)(nreps+1);
  
  Eigen::MatrixXd R( ny , 4 );
  R << B.row(idx).transpose() , T , U , F ;

  //
  // Store results
  //

  
  
  //
  // Report significant clusters 
  //

  std::map<int,int> inclst;
    
  int cnt = 0 , pos = 0;
  std::set<cpt_cluster_t>::const_iterator cc = clusters.clusters.begin();
  while ( cc != clusters.clusters.end() )
    {

      // get cluster-corrected p-value
      clusters.perm[ cnt ] /= (double)( nreps+1 ) ;
      
      // only track # of significant clusters
      if ( clusters.perm[ cnt ] <= 0.05 )
	{
	  ++pos;	      
	  
	  results.cluster_emp[ vname[ cc->seed ] ] = clusters.perm[ cnt ] ;
	  
	  std::set<int>::const_iterator ii = cc->members.begin();
	  while ( ii != cc->members.end() )
	    {	      
	      inclst[ *ii ] = pos;	      
	      results.cluster_members[ vname[ cc->seed ] ].insert( vname[ *ii ] );
	      ++ii;
	    }
	}
      
      ++cnt;
      ++cc;
    }

  logger << "  " << pos << " clusters significant at corrected empirical P<0.05\n";
  

  
  //
  // Point-wise results
  //
  
  for (int y=0; y<ny; y++)
    {
      
      results.beta[ vname[y] ] = R(y,0);
      results.t[ vname[y] ] = R(y,1);		  
      results.emp[ vname[y] ] = R(y,2);
      results.emp_corrected[ vname[y] ] = R(y,3);
      // cluster membership?
      if ( inclst.find( y ) != inclst.end() )
	results.inclst[ vname[y] ] =  inclst[y];      
    }
  
  return results;
}



Eigen::VectorXd cpt_t::get_tstats( const Eigen::VectorXd & B ,
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
  


//
// Find clusters
//

struct cpt_sorter_t {
  cpt_sorter_t( double s , int v ) : stat( fabs(s) ) , v(v) { } 
  double stat;
  int v;
  bool operator< ( const cpt_sorter_t & rhs ) const
  {
    if ( stat > rhs.stat ) return true;  // i.e. largest first
    if ( stat < rhs.stat ) return false;
    return v < rhs.v;
  }
};


cpt_clusters_t::cpt_clusters_t( const Eigen::VectorXd & T ,
				double threshold ,
				const std::map<int,std::set<int> > & adj , 
				bool two_sided , 
				bool verbose , 
				const std::vector<std::string> * labels )
{
  
  max_stat = 0;
  clusters.clear();
  
  // no clustering if threhsold is neg.
  if ( threshold < 0 ) return;

  // get ordered statstics
  
  std::set<cpt_sorter_t> o;
  const int n = T.size();
  for (int i=0; i<n; i++)
    o.insert( cpt_sorter_t( T[i] , i ) );

  // start looking through 
  std::set<cpt_sorter_t>::const_iterator oo = o.begin();
  std::set<int> clustered;
  while ( oo != o.end() )
    {
    
      
      // if below thresjold, all done
      if ( oo->stat < threshold ) break;
      
      if ( verbose ) std::cout << "  flagging " << (*labels)[ oo->v ] << " = " << oo->stat << " > threshold = " << threshold << "\n";

      // is this variable already spoken for? go to next 
      if ( clustered.find( oo->v ) != clustered.end() )
	{
	  if ( verbose ) std::cout << "  --- already spoken for\n";  
	  ++oo;
	  continue;
	}

      // if here, we have a new seed:
      cpt_cluster_t cluster;
      cluster.seed = oo->v;
      cluster.members.insert( cluster.seed );
      cluster.stat = oo->stat; // start sum
      clustered.insert( cluster.seed );
      
      // add all friends -- and then friends of friends
      if ( adj.find( cluster.seed ) != adj.end() )
	{
	  // individuals to test
	  std::set<int> friends = adj.find( cluster.seed )->second ;
	  
	  while ( friends.size() ) 
	    {
	      
	      std::set<int> newfriends;

	      if ( verbose ) 
		std::cout << "  --- considering " << friends.size() << " friends\n";
	  
	      std::set<int>::const_iterator ff = friends.begin();
	      while ( ff != friends.end() )
		{
		  
		  if ( verbose ) std::cout << "     -- " << (*labels)[ *ff ] << " ";
		  
		  if ( clustered.find( *ff ) == clustered.end() &&
		       fabs( T[ *ff ] ) >= threshold )
		    {
		      
		      if ( verbose ) std::cout << " above threshold ";
		      
		      // if two-sided, only cluster groups that have same direction of effect
		      if ( (!two_sided) || ( T[ *ff ] <= 0 == T[ cluster.seed ] <= 0 ) )
			{
			  if ( verbose ) std::cout << "  adding ";
			  cluster.members.insert( *ff );
			  cluster.stat += fabs( T[ *ff ] );
			  clustered.insert( *ff );
			  
			  // loop in new friends
			  const std::set<int> nf = adj.find( *ff )->second;
			  std::set<int>::const_iterator nn = nf.begin();
			  while ( nn != nf.end() )
			    {
			      if ( clustered.find( *nn ) == clustered.end() )
				{
				  bool above_threshold = fabs( T[ *nn ] ) >= threshold ;
				  bool direction = (!two_sided) || ( T[ *nn ] <= 0 == T[ cluster.seed ] <= 0 );
				  if ( above_threshold && direction ) 
				    newfriends.insert( *nn );
				}
			      ++nn;
			    }
			}
		    }
		  
		  if ( verbose ) std::cout << "\n";
		  
		  ++ff;
		}	  
	      
	      // go back and loop in friends of friends
	      
	      friends = newfriends;
	      
	    }

	}
      
      // save this cluster 
      clusters.insert( cluster );
      
      // carry on down the list
      ++oo;
    }
  
  // get maximum cluster statistic

  std::set<cpt_cluster_t>::const_iterator ii = clusters.begin();
  while ( ii != clusters.end() )
    {
      if ( ii->stat > max_stat )
	max_stat = ii->stat;
      ++ii;
    }

  // initiate vector for counting empirical results

  perm.resize( clusters.size() , 1 ); // 1 to include original 
  
  
}


