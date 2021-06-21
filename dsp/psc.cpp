
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


#include "dsp/psc.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"

std::vector<std::string> psc_t::vname;
Eigen::Array<double, 1, Eigen::Dynamic> psc_t::means;
Eigen::Array<double, 1, Eigen::Dynamic> psc_t::sds;
Eigen::VectorXd psc_t::W;
Eigen::MatrixXd psc_t::V;   

extern writer_t writer;


void psc_t::construct( param_t & param )
{

  //
  // Always require input from 1+ output files (i.e. PSD or COH)
  //
  
  if ( ! param.has( "spectra" ) ) Helper::halt( "no spectra=<files> specified" );
  std::vector<std::string> infiles = param.strvector( "spectra" );

  //
  // which variables to pull from these?
  //

  if ( ! param.has( "v" ) ) Helper::halt( "no v=<variables> specified" );
  std::set<std::string> vars = param.strset( "v" );

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
  
  bool drop_incomplete_rows = param.yesno( "drop-incomplete-rows" );
  
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
  // any frequencies limits?
  //

  double flwr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 0 ;
  double fupr = param.has( "f-upr" ) ? param.requires_dbl( "f-upr" ) : 0 ;
  
  //
  // Epoch level input?  (in which case,  ID ==>  ID:E internally
  //

  bool epoch = param.yesno( "epoch" );

  //
  // Notes that all features are 1) signed, and 2) pairwise (CH1xCH2)
  //

  bool signed_stats = param.yesno( "signed-pairwise" );
  
  //
  // Save projection?
  //
  
  std::string projection = param.has( "proj" ) ? param.value( "proj" ) : "" ;

    
  //
  // Dump component definitions and raw data in a separate text file
  //

  std::string vdump = param.has( "dump" ) ? param.value( "dump" ) : "" ;
  
  
  //
  // Report input variables by quantile of each PSC? (for v-dump only)
  //
  
  int q = param.has( "q" ) ? param.requires_int( "q" ) : 0;
  
  // if ( q && vdump == "" ) 
  //   Helper::halt( "can only report feature quantile-stratified medians in v-matrix=<file> mode" );
  
  if ( q < 0 || q > 10 ) 
    Helper::halt( "q should be between 1 and 10" );
  
  
  //
  // PSC parameters
  //
  
  int nc = param.has( "nc" ) ? param.requires_int( "nc" ) : 10 ;

  std::vector<double> th;
  if ( param.has( "th" ) ) th = param.dblvector( "th" );
  
  const bool standardize_inputs = param.yesno( "norm" ) ;

		 
  //
  // Read spectra into matrix X
  //


  // can stratify by channel, frequency (columns)
  // expect input 
  // can also have coherence (CH1/CH2) values too
  // expect file to contain multiple individuals
  
  // in output, one row per individual; columns are CH x F (xCH1/CH2)
  
  // id -> ch -> f -> var -> value
  
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > i2c2f2v;

  for (int i=0; i<infiles.size(); i++)
    {
      std::string infile = Helper::expand( infiles[i] );

      logger << "  reading spectra from " << infile << "\n";

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
      if ( epoch && cols.find( "E" ) == cols.end() ) Helper::halt( "no E column in " + infile );
      bool ch1 = cols.find( "CH" ) != cols.end();
      bool ch2 = cols.find( "CH1" ) != cols.end() && cols.find( "CH2" ) != cols.end();
      bool frq = cols.find( "F" ) != cols.end() ;
      
      // channel and frequency are optional
      // if ( ch1 == ch2 ) Helper::halt( "require either CH or CH1 & CH2 in " + infile );
      // if ( cols.find( "F" ) == cols.end() ) Helper::halt( "no F column in " + infile );

      int id_slot = -1;
      int e_slot = -1;
      int ch_slot = -1;
      int ch1_slot = -1;
      int ch2_slot = -1;
      int f_slot = -1;      

      std::map<int,std::string> slot2var;
      for (int i=0;i<tok.size();i++)
	{
	  if ( tok[i] == "ID" ) id_slot = i;
	  if ( tok[i] == "E" ) e_slot = i;
	  if ( tok[i] == "F" ) f_slot = i;
	  if ( tok[i] == "CH" ) ch_slot = i;
	  if ( tok[i] == "CH1" ) ch1_slot = i;
	  if ( tok[i] == "CH2" ) ch2_slot = i;	  
	  if ( vars.find( tok[i] ) != vars.end() ) slot2var[i] = tok[i];
	}
      
      if ( slot2var.size() == 0 ) 
	Helper::halt( "no variables v=<...> in " + infile );

      const int ncols = tok.size();
      
      //
      // all set, now start reading
      //

      while ( ! IN1.eof() )
	{
	  // looking for ID F CH         --> PSD 
	  //  OR         ID F CH1 CH2    --> LCOH (default)  
	  //  OR         ID              --> HYPNO

	  // in epoch mode, looking for ID and E... we merge these to ID = ID:E

	  std::string line;

	  Helper::safe_getline( IN1 , line );

	  std::vector<std::string> tok = Helper::parse( line , "\t" );

	  if ( IN1.eof() || tok.size() == 0 ) continue;
	  
	  if ( tok.size() != ncols ) Helper::halt( "incorrect number of columns in " + infile );
	  
	  std::string id = epoch ? 
	    tok[ id_slot ] + ":" + tok[ e_slot ] : 
	    tok[ id_slot ] ;
	  
	  // skip if person not on an include list
	  if ( id_includes.size() != 0 && id_includes.find( id ) == id_includes.end() ) continue;

	  // skip if person is on an exclude list
	  if ( id_excludes.size() != 0 && id_excludes.find( id ) != id_excludes.end() ) continue;	  
	  
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
  // output on transformations
  //

  if ( param.has( "dB" ) )
    logger << "  taking 10log10(X) of " << param.value( "dB" ) << "\n";

  if ( flwr > 0 || fupr > 0 )
    {
      logger << "  restricting to ";
      if ( flwr > 0 ) logger << flwr << " <= ";
      logger << "F";
      if ( fupr > 0 ) logger << " <= " << fupr;
      logger << "\n";
    }
  
  
  //
  // Construct data matrix
  //

  logger << "  converting input spectra to a matrix\n";
  
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
    Helper::halt( "not enough observationns for PSC analysis" );

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
  
  Eigen::MatrixXd U( rows.size() , cols.size() );
  
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
		 		  
		  // find channel?
		  const std::map<std::string,std::map<std::string,std::map<std::string,double> > > & dat = i2c2f2v.find( *ii )->second;
		  //if ( dat.find( ss1->first ) == dat.end() ) Helper::halt( "no channel " + ss1->first + " for individual " + *ii );

		  // find frequency?
		  const std::map<std::string,std::map<std::string,double> > & dat2 = dat.find( ss1->first )->second;
		  //if ( dat2.find( ss2->first ) == dat2.end() ) Helper::halt( "no frequency " + ss2->first + " for individual " + *ii );
		    		    
		  // find variable?
		  const std::map<std::string,double> & dat3 = dat2.find( ss2->first )->second;
		  //if ( dat3.find( ss3->first ) == dat3.end() ) Helper::halt( "no variable " + ss3->first + " for individual " + *ii );
		  
		  // these have all been checked now: should be okay to add to the store
		  U( row, ss3->second ) = dat3.find( ss3->first)->second;
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

  int N = U.rows();
  means = U.colwise().mean();
  sds = ((U.array().rowwise() - means ).square().colwise().sum()/(N-1)).sqrt();
  for (int i=0; i < sds.size(); i++)
    if ( sds[i] < EPS ) 
      {
	Helper::halt( "at least one invariant column in input; first = : " + vname[i] );
      }


  //
  // Outliers?
  //
  
  int ni = rows.size();
  
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
	  Eigen::VectorXd::Map( &tmp[0], ni ) = U.col(j);

	  int removed = MiscMath::outliers( &tmp , th[t] , &inc , &prior);
	  
	  //logger << "  removing " << removed << " var " << j << " round " << t << "\n";
	}
      
    }

  
  //
  // Remove rows from input U
  //

  Eigen::MatrixXd U2 = U;

  std::vector<std::string> id2 = id;

  U.resize(0,0);
  id.clear();

  ni = 0;
  for (int i=0;i<inc.size();i++)
    if ( inc[i] ) ++ni;
  
  logger << "  after outlier removal, " << ni << " individuals remaining\n";
  
  U.resize( ni , nv );
  id.resize( ni );
  int r = 0;
  for (int i=0;i<U2.rows();i++)
    if ( inc[i] )
      {
	for (int j=0;j<nv;j++) 
	  U(r,j) = U2(i,j);
	id[r] = id2[i];
	++r;
      }

  //
  // Track new means / SDs of the nv input variables
  //
  
  N = ni;
  means = U.colwise().mean();
  sds = ((U.array().rowwise() - means ).square().colwise().sum()/(N-1)).sqrt();

  //
  // for directional (coherence) terms, also optionally track means separately for positive and 
  // negative components;  0 if not defined.
  //

  std::vector<double> pos_means, neg_means;

  if ( signed_stats )
    {
      pos_means.resize( nv , 0 );
      neg_means.resize( nv , 0 );
      
      for (int k=0; k<nv; k++)
	{
	  double psum = 0 , nsum = 0;
	  int    pcnt = 0 , ncnt = 0;

	  const Eigen::VectorXd & col = U.col( k );
	  
	  for (int i=0; i<ni; i++)
	    {
	      if      ( col(i) > 0 ) { psum += col(i) ; ++pcnt; }
	      else if ( col(i) < 0 ) { nsum += col(i) ; ++ncnt; }
	    }
	  
	  pos_means[k] = pcnt == 0 ? 0 : psum / (double)pcnt;
	  neg_means[k] = ncnt == 0 ? 0 : nsum / (double)ncnt;

	}
      
    }


  //
  // Get copy of original input data, before normalization
  //

  Eigen::MatrixXd O = U;


  
  //
  // SVD
  //

  if ( nc > nv )
    {
      logger << "  reducing nc to the number of features, " << nv << "\n";
      nc = nv;
    }

  //
  // First standardize columns? 
  //
  
  if ( standardize_inputs )
    {
      logger << "  standardizing data matrix\n";      
      U.array().rowwise() -= means;
      U.array().rowwise() /= sds;
    }
  else
    {
      logger << "  mean-centering data matrix\n";
      U.array().rowwise() -= means;      
    }


  //
  // need to change 'nc' ? 
  //

  int maxk = U.rows() < U.cols() ? U.rows() : U.cols();
  if ( nc > maxk )
    {
      logger << "  reducing nc to " << maxk << "\n";
      nc = maxk;
    }
      
  //
  // SVD
  //

  logger << "  about to perform SVD...\n";

  Eigen::BDCSVD<Eigen::MatrixXd> svd( U , Eigen::ComputeThinU | Eigen::ComputeThinV );
  U = svd.matrixU();
  V = svd.matrixV();
  W = svd.singularValues();

  logger << "  done... now writing output\n";
  
  // components are already sorted in decreasing order, so so just take the first 'nc' components
  //  assuming ni << nv , and so ni components returned, of which we take nc (which is less than ni )

  //  U[ ni x ni ]  -->   U[ ni x nc ]
  //  V[ nv x ni ]  -->   V[ nv x nc ]
  //  W[ ni ]       -->   W[ nc ]

  //
  // Output
  //

  writer.id( "." , "." );

  // U ( x ID )

  for (int i=0;i<ni;i++)    
    {
      writer.id( id[i] , "." );
      for (int j=0;j<nc;j++)
	{
	  writer.level( j+1 , "PSC" );
	  // either as factor or variable?  col2ch col2f col2var
	  writer.value( "U" , U(i,j) );
	}
      writer.unlevel( "PSC" );
    }


  //
  // PSC frequency peaks, and channel means: 
  //  a) based on MEAN, SD, then PSCs
  //    sum over channels/channel pairs
  // to get the frequencies of maximal variability for that PSC
  //

  if ( 0 ) 
    {
      double s=0 , s_sd = 0 ;
      std::map<double,double> sfrqs, sfrqs_sd;
      std::map<std::string,double> schs, schs_sd;
      for (int k=0; k<nv; k++)
	{
	  
	  sfrqs[ col2f[ vname[k] ] ] += means[k] * means[k];
	  sfrqs_sd[ col2f[ vname[k] ] ] += sds[k] * sds[k];
	  
	  const std::string & ch1 = col2ch1[ vname[k] ];
	  schs[ ch1 ] += means[k] * means[k];
	  schs_sd[ ch1 ] += sds[k] * sds[k];
	  
	  const std::string & ch2 = col2ch2[ vname[k] ];
	  if ( ch2 != "." ) 
	    {
	      schs[ ch2 ] += means[k] * means[k];
	      schs_sd[ ch2 ] += sds[k] * sds[k];
	    }
	}
      
      // normalize each to 1.0
      double sf = 0 , sf_sd = 0 , sc = 0 , sc_sd = 0;
      
      // report by frequency
      std::map<double,double>::iterator kk = sfrqs.begin();
      {
	sf += kk->second;
	sf_sd += sfrqs_sd[ kk->first ];
	++kk;
      }
      
      kk = sfrqs.begin();
      while ( kk != sfrqs.end() )
	{
	  writer.level( kk->first , globals::freq_strat );
	  writer.value( "FW" , kk->second / sf );
	  writer.value( "FW_SD" , kk->second / sf_sd );
	  ++kk;
	}
      writer.unlevel( globals::freq_strat );
      
      // report by channels
      std::map<std::string,double>::iterator cc = schs.begin();
      while ( cc != schs.end() )
	{
	  sc += cc->second;
	  sc_sd += schs_sd[ cc->first ];
	  ++cc;
	}
      
      cc = schs.begin();
      while ( cc != schs.end() )
	{
	  writer.level( cc->first , globals::signal_strat );
	  writer.value( "FW" , cc->second / sc );
	  writer.value( "FW_SD" , cc->second / sc_sd );
	  ++cc;
	}
          

  //
  // Same by PSC
  //
  
  for (int j=0;j<nc;j++)
    {
      writer.level( j+1 , "PSC" );
      
      const Eigen::VectorXd & vec = V.col(j);
      
      // add squared term across channels for each freq
      double s =0;
      std::map<double,double> sfrqs;
      for (int k=0; k<nv; k++)
	{
	  s += vec(k) * vec(k);
	  sfrqs[ col2f[ vname[k] ] ] += vec(k) * vec(k);
	}

      std::map<double,double>::const_iterator kk = sfrqs.begin();
      while ( kk != sfrqs.end() )
	{
	  writer.level( kk->first , globals::freq_strat );
	  writer.value( "FW" , kk->second / s );
	  ++kk;
	}
      writer.unlevel( globals::freq_strat );
      
    }
  writer.unlevel( "PSC" );


  //
  // Repeat for PSC by channel summaries
  //

  for (int j=0;j<nc;j++)
    {
      writer.level( j+1 , "PSC" );
      
      const Eigen::VectorXd & vec = V.col(j);
      
      std::map<std::string,double> schs;
      for (int k=0; k<nv; k++)
	{
	  const std::string & ch1 = col2ch1[ vname[k] ];
	  schs[ ch1 ] += vec(k) * vec(k);
	  
	  const std::string & ch2 = col2ch2[ vname[k] ];
	  if ( ch2 != "." )
	    schs[ ch2 ] += vec(k) * vec(k);	  
	}
      
      double sum = 0;
      std::map<std::string,double>::const_iterator cc = schs.begin();
      while ( cc != schs.end() )
	{
	  sum += cc->second;
	  ++cc;
	}
      
      cc = schs.begin();
      while ( cc != schs.end() )
	{
	  writer.level( cc->first , globals::signal_strat );
	  writer.value( "FW" , cc->second / sum );
	  ++cc;
	}
      writer.unlevel( globals::signal_strat );
      
    }
  writer.unlevel( "PSC" );
  
  
    }


  //
  // Make quantile summaries
  //

  // variable --> PSC --> QUANTILE --> mean
  // if this is populated below, will be output in vdump
  std::map<std::string,std::map<int,std::map<int,double> > > qsumms;

  if ( q && ! signed_stats )
    {

      // loop over each PSC
      //  for each quantile, - get mean value of features
      
      // col2ch[ vname[k] ] 

      //
      // For each PSC
      //

      for (int j=0; j<nc; j++)
	{
	  
	  //
	  // 1) get PSC across all individuals
	  //
	  
	  const Eigen::VectorXd & psc = U.col(j);
	  	  
	  // helper to get quantiles;
	  std::set<psc_sort_t> sp;
	  for (int i=0; i<ni; i++)
	    sp.insert( psc_sort_t( i , psc(i) ) );
	  	  
	  std::vector<int> qt = psc_sort_t::quantile( sp , q );

	  //
	  // 2) for each quantile, then:
	  //    get means across all quantiles
	  //

	  for (int qq=0; qq<q; qq++)
	    {
	      
	      int nq = 0;
	      for (int i=0;i<ni; i++) 
		if ( qt[i] == qq ) ++nq;
	      
	      // get means for this quantile
	      std::vector<double> xx( nv , 0 );
	      for (int i=0; i<ni; i++)
		{
		  if ( qt[i] == qq ) 
		    for (int k=0; k<nv; k++) xx[k] += O(i,k);
		}	      
	      for (int k=0; k<nv; k++) 
		xx[k] /= (double)nq;
	      
	      // output
	      for (int k=0; k<nv; k++) 
		qsumms[ vname[ k ] ][ j ][ qq ] = xx[k] ;
	      
	    } // next quantile
	  	  
	} // next PSC     
      
    }


  //
  // Special case: Signed stats & coherence: make POS and NEG channel-level summaries
  //

  if ( q && signed_stats )
    {
      // loop over each PSC
      //  for each quantile, - get median values of features (pairwise)
      //                     - collapse to channel-level summaries (POS/NEG separately)
      
      // this assumes that ALL features of CH1 x CH2 
      // get these channels now:
      
      std::map<std::string,int> ch2slot;
      std::map<int,std::string> slot2ch;

      std::vector<int> ch1( nv ), ch2( nv );
      std::vector<double> frq( nv );

      for (int k=0; k<nv; k++)
	{
	  std::vector<std::string> ctok = Helper::parse( col2ch[ vname[k] ] , "." ) ;
	  
	  if ( ctok.size() != 2 ) 
	    Helper::halt( "q & signed-stats requires that all features are pairwise stats: CH1 x CH2" );

	  if ( ch2slot.find( ctok[0] ) == ch2slot.end() ) 
	    {
	      int sz = ch2slot.size();
	      ch2slot[ ctok[0] ] = sz;
	      slot2ch[ sz ] = ctok[0];
	    }
	  
	  if ( ch2slot.find( ctok[1] ) == ch2slot.end() )
            {
	      int sz = ch2slot.size();
              ch2slot[ ctok[1] ] = sz;
              slot2ch[ sz ] = ctok[1];
            }
	  
	  ch1[ k ] = ch2slot[ ctok[0] ];
	  ch2[ k ] = ch2slot[ ctok[1] ];	  

	  // store just for quick lookup below (prob. no diff/not needed)
	  frq[ k ] = col2f[ vname[k] ];
	  
	}
      
      const int nch = ch2slot.size();

      //
      // For each PSC
      //

      for (int j=0; j<nc; j++)
	{
	  
	  writer.level( j+1 , "PSC" );

	  //
	  // 1) get PSC across all individuals
	  //

	  const Eigen::VectorXd & psc = U.col(j);

	  // helper to get quantiles;
	  std::set<psc_sort_t> sp;
	  for (int i=0; i<ni; i++)
	    sp.insert( psc_sort_t( i , psc(i) ) );
	  
	  std::vector<int> qt = psc_sort_t::quantile( sp , q );
	  
	  	  
	  //
	  // 2) for each quantile, then:
	  //
	  // 2a) get means across all quantiles
	  // 2b) Collapse to frequency-specific, channel-level summaries
	  // 2c) Output the channel level summaries
	  //
	  
	  for (int qq=0; qq<q; qq++)
	    {
	      
	      writer.level( qq+1 , "Q"  );
	      
	      int nq = 0;
	      for (int i=0;i<ni; i++) 
		if ( qt[i] == qq ) ++nq;
	      
	      // get means for this quantile
	      std::vector<double> xx( nv , 0 );
	      for (int i=0; i<ni; i++)
		{
		  if ( qt[i] == qq ) 
		    for (int k=0; k<nv; k++) 
		      xx[k] += O(i,k);
		}	      

	      for (int k=0; k<nv; k++) 
		xx[k] /= (double)nq;
	      
	      // collapse to [ch][freq]
	      std::vector<std::map<double,double> > pos( nch );
	      std::vector<std::map<double,double> > neg( nch );
	      for (int k=0; k<nv; k++)
		{
		  if ( xx[k] > 0 ) 
		    {
		      // +ve:   A -> B   means A=POS, B=NEG
		      pos[ ch1[k] ][ frq[k] ] += xx[k]; 
		      neg[ ch2[k] ][ frq[k] ] -= xx[k]; 
		    }
		  else
		    {
		      // -ve: A -> means A=NEG, B=POS
		      neg[ch1[k]][ frq[k] ] += xx[k]; // i.e. adds a negative value for A
		      pos[ch2[k]][ frq[k] ] -= xx[k]; //      adds a positive value for B		      
		    }
		}
	      
	      //
	      // normalize within ch/freq combo
	      //

	      for (int c=0;c<nch;c++) 
		{
		  writer.level( slot2ch[ c ] , globals::signal_strat );  

		  int nn = pos[c].size();
		  
		  std::map<double,double>::const_iterator ff = pos[c].begin();
		  while ( ff != pos[c].end() )
		    {		      
		      writer.level( ff->first , globals::freq_strat ); 		      
		      writer.value( "POS" , ff->second / (double)nn );
		      writer.value( "NEG" , neg[ c ][ ff->first ] / (double)nn  );		  
		      ++ff;
		    }
		  writer.unlevel( globals::freq_strat );
		}
	      writer.unlevel( globals::signal_strat );
	    }
	  writer.unlevel( "Q" );
	}
      writer.unlevel( "PSC" );
    }
  


  //
  // Output V to a file, along with variable information (channel, variable, freq)
  //

  if ( vdump != "" ) 
    {
      
      logger << "  dumping V and meta-information to file: " << vdump << ".vars \n";

      std::ofstream V1( (vdump+".vars").c_str() , std::ios::out );
      
      // V : first nc component
      
      V1 << "VAR\tCH\tCH1\tCH2\tF";

      V1 << "\tMN\tSD";
      
      // components
      for (int c=0;c<nc; c++ )
	V1 << "\tV" << c+1 ;      
            
      // quantiles for each variable?
      if ( q > 0 ) 
	for (int c=0;c<nc; c++ )
	  for (int qq=0;qq<q; qq++ )
	    V1 << "\tV" << c+1 << ".Q" << qq+1 ;
      
      V1 << "\n";
      
      //
      // Data rows
      //

      for (int k=0;k<nv;k++)
	{
	  // ch ~ f ~ var
	  // ch1.ch2 ~ f ~ var

	  // variable name
	  V1 << col2var[ vname[k] ];

	  // channel(s)
	  std::vector<std::string> ctok = Helper::parse( col2ch[ vname[k] ] , "." ) ;
	  if ( ctok.size() == 1 ) 
	    V1 << "\t" << ctok[0] << "\t.\t.";
	  else if ( ctok.size() == 2 )  
	    V1 << "\t." << "\t" << ctok[0] << "\t" << ctok[1] ;
	  else 
	    Helper::halt( "bad format in vname/channel");
	  
	  // freq
	  V1 << "\t" << col2f[ vname[k] ];
	  
	  // means / SD from raw data 
	  V1 << "\t" << means[ k ]
	     << "\t" << sds[ k ];
	  
	  // V coefficients
	  for (int c=0;c<nc; c++ )
	    V1 << "\t" << V(k,c);
	  	  
	  // Quantile means?
	  if ( q > 0 ) 
	    for (int c=0;c<nc; c++ )
	      for (int qq=0;qq<q; qq++ )
		V1 << "\t" << qsumms[ vname[k] ][ c ][ qq ] ;

	  // done
	  V1 << "\n";
	}
      V1.close();

      

      //
      // Input matrix (and PSC)
      //
      
      logger << "  dumping U and PSCs to file: " << vdump << ".data\n";
      
      std::ofstream V2( (vdump+".data").c_str() , std::ios::out );
      
      V2 << "ID";
      
      // components
      for (int c=0;c<nc; c++ )
	V2 << "\tPSC" << c+1 ;      
      
      // variables
      for (int c=0; c<nv; c++)
	V2 << "\t" << vname[ c ] ; 
      
      V2 << "\n";
      
      //
      // Data rows
      //
           
      for (int i=0;i<ni;i++)
	{
	  
	  V2 << id[i] ;
	  
	  // components
	  for (int j=0;j<nc;j++)
	    V2 << "\t" << U(i,j);
	  
	  // variables  (O is original data)
	  for (int c=0;c<nv; c++ )
	    V2 << "\t" << O(i,c);	      

	  V2 << "\n";
	}
      
      V2.close();
        
    }
  

  //
  // Output to standard DB
  //
  
  // W
  
  double wsumsq = W.array().square().sum();
  
  double cve = 0;
  
  for (int j=0;j< W.size(); j++)
    {
      writer.level( j+1 , "I" );
      writer.value( "W" , W(j) );
      double ve = ( W(j) * W(j) ) / wsumsq;
      cve += ve;
      writer.value( "VE" , ve );
      writer.value( "CVE" , cve );	  
      writer.value( "INC" , j < nc ? 1 : 0 );
    }
  
  writer.unlevel( "W" );
  
  
  // V (transpose) (only first 'nc' components, as U/V etc are sorted 
  // in decreasing order
  
  if ( 1 || vdump == "" ) 
    {
      for (int j=0;j<nc;j++)
	{
	  writer.level( j+1 , "I" );
	  for (int k=0;k<nv;k++)
	    {
	      writer.level( vname[k] , "J" );
	      writer.value( "V" , V(k,j) );	// nb transpose
	    }
	  writer.unlevel( "J" );
	}
      writer.unlevel( "I" );
      
      // VARS
      for (int j=0;j<nv;j++)
	{
	  writer.level( vname[j] , "J" );
	  writer.value( "CH" , col2ch[ vname[j] ] );
	  writer.value( "F" , col2f[ vname[j] ] );
	  writer.value( "VAR" , col2var[ vname[j] ] );
	}
      writer.unlevel( "J" );      
    }
  
  
  
  //
  // Output projection to a separate file?
  // How to keep track of CH F and VAR to make sure that new projections
  // will line up?  Need to output that too
  //
  
  if ( projection != "" ) 
    {

      logger << "  writing projection to " << projection << "\n";
      
      std::ofstream OUT1( projection.c_str() , std::ios::out );

      // variables (with mean/SD in reference population)
      OUT1 << "NV: " << nv ;
      for (int j=0;j<nv;j++)
	OUT1 << " " << vname[j]
	     << " " << means[j]
	     << " " << sds[j];
      OUT1 << "\n";

      // means/SDs in reference population
      
      // components
      OUT1 << "NC: " << nc << "\n";
      
      OUT1 << "W:";
      for (int i=0;i<nc;i++)
	OUT1 << " " << W(i);
      OUT1 << "\n";
      
      OUT1 << "V:";
      for (int i=0;i<nv;i++)
	for (int j=0;j<nc;j++)
	  OUT1 << " " << V(i,j);
      OUT1 << "\n";   
      
      OUT1.close();
    }
  
}


void psc_t::attach( param_t & param )
{

  //
  // check if already attached
  //

  if ( W.size() != 0 ) return;

  //
  // Read W and V matrices in this file
  //

  std::string infile = param.requires( "proj" );

  if ( ! Helper::fileExists( infile ) )
    Helper::halt( "could not find " + infile );
  
  logger << "  reading projection from " << infile << "\n";
  
  std::ifstream IN1( infile.c_str() , std::ios::in );


  
  //
  // variables
  //

  std::string dummy;
  int nv;
  
  IN1 >> dummy >> nv;
  vname.resize( nv );
  means.resize( nv );
  sds.resize( nv );
  
  for (int j=0;j<nv;j++)
    IN1 >> vname[j] >> means(j) >> sds(j);
  
  // components
  IN1 >> dummy >> nc;
        
  W.resize( nc );
  V.resize( nv , nc );
    
  
  // W 
  IN1 >> dummy;
  for (int i=0;i<nc;i++)
    IN1 >> W(i);
  
  // V
  IN1 >> dummy;
  for (int i=0;i<nv;i++)
    for (int j=0;j<nc;j++)
      IN1 >> V(i,j);
  
  IN1.close();


  //
  // reformat of W for projection
  //
  
  W = W.cwiseInverse();

  //
  // Set some components to zero potential?
  //

  //
  // all PSC or a subset (via 'nc=' or 'drop/keep=')
  //
  
  logger << "  found " << nc << " PSCs based on " << nv << " variables\n";

  if ( param.has( "nc" ) )
    {
      int k = param.requires_int( "nc" ) ;

      if ( k > nc )
	Helper::halt( "requested " + Helper::int2str( k )
		      + " PSCs but only " + Helper::int2str( nc ) + " present" );

      // set to zero
      logger << "  subsetting to the first " << k << " of " << nc << " PSCs\n";
      for (int i=k;i<nc;i++) W[i] = 0;
    }

  std::vector<int> drop, keep;
  if ( param.has( "drop" ) ) drop = param.intvector( "drop" );
  if ( param.has( "keep" ) ) keep = param.intvector( "keep" );

  if ( drop.size() != 0 && keep.size() != 0 )
    Helper::halt( "cannot specify both drop and keep" );

  // if 'keep', start with all set to drop (T)
  // else if 'drop' start with none set to drop (F)
  std::vector<int> to0( nc , keep.size() != 0 );
  for (int i=0; i<drop.size(); i++)
    {
      if ( drop[i] < 1 || drop[i] > nc ) Helper::halt( "drop parameter out of range" );
      to0[ drop[i] - 1 ] = true; // nb. convert 1-based input to 0-base
    }
  
  for (int i=0; i<keep.size(); i++)
    {
      if ( keep[i] < 1 || keep[i] > nc ) Helper::halt( "keep parameter out of range" );
      to0[ keep[i] - 1 ] = false; // nb. convert 1-based input to 0-base 
    }  

  if ( drop.size() ) logger << "  dropping " << drop.size() << " of " << nc << " components\n";
  if ( keep.size() ) logger << "  retainging only " << keep.size() << " of " << nc << " components\n";
  
  if ( drop.size() + keep.size() != 0 )
    for (int i=0; i<to0.size(); i++)
      if ( to0[i] ) W(i) = 0;
  
  //  std::cout << "1/W = " << W << "\n";  
    
}


void psc_t::project( edf_t & edf , param_t & param )
{
  
  //
  // Cache w/ spectral data for this individual 
  //

  // i.e. expecting to be given a cache from PSD (or similar)
  // that can populate the vnames above
  //
  // PSC vnames are in the format  CH_F_VAR
  // which would correspond to cache[ VAR ] , ckey( "CH" , "F" )
  //
  // COH vnames are in the format  CH1.CH2_F_VAR  [not yet implemented]
  // which would correspond to cache[ COH ] , ckey( "CH1" ,"CH2" ,  "F" )
  // etc...

  std::string cache_name = param.requires( "cache" );

  if ( ! edf.timeline.cache.has_num( cache_name ) )
    Helper::halt( "cache not found for this individual: " + cache_name );

  cache_t<double> * cache = edf.timeline.cache.find_num( cache_name );

  const bool norm = param.yesno( "norm" );
  
  // see which variables exist, i.e. psc_t::vname[] 

  int nv = vname.size();

  Eigen::VectorXd X( nv );
  
  for (int i=0; i<nv; i++)
    {
      std::vector<std::string> tok = Helper::parse( vname[i] , "~" );
      if ( tok.size() != 3 ) 
	Helper::halt( "bad format for PSC vnames, expecting 3 fields, '~'-delimited" );


      // can we find this value in the cache?  

      // variable =  tok[2]

      ckey_t key( tok[2] );

      // channel(s) = tok[0]
      
      std::vector<std::string> tokch = Helper::parse( tok[0] , "." );

      if ( tokch.size() == 2 )
	{
	  key.add( "CH1" , tokch[0] );
	  key.add( "CH2" , tokch[1] );
	}
      else if ( tokch.size() == 1 )
	{
	  key.add( "CH" , tok[0] );
	}
      else
	Helper::halt( "bad format for PSC vname: ch " + tok[0] );
      
      // frequency = tok[1] 

      if ( tok[1] != "0" ) 
	{
	  double f;
	  if ( ! Helper::str2dbl( tok[1] , &f ) ) 
	    Helper::halt( "bad frequency value in PSC vname" );
	  
	  key.add( "F" , f );
	}

      // fetch 

      std::vector<double> cx = cache->fetch( key );
      
      if ( cx.size() != 1 )
	Helper::halt( "could not find cached variable: " + vname[i] );
      
      X(i) = cx[0];
    }
  
  logger << "  all " << nv << " features found in the cache\n";


  //
  // Mean center data
  //

  // remove original mean

  X.array() -= means;

  // optionally, scale by original population SD

  if ( norm ) 
    X.array() /= sds;

    
  //
  // Project given W and V to get U, PSCs for this individual
  //
  

  Eigen::MatrixXd U_proj = X.transpose() * V * W.asDiagonal();

  
  //
  // Output
  //

  for (int i=0; i<nc; i++)
    {
      writer.level( i+1 , "PSC" );
      writer.value( "U" , U_proj(0,i) );
    }
  writer.unlevel( "PSC" );
  
}
