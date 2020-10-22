
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


Data::Vector<double> psc_t::W;
Data::Matrix<double> psc_t::DW;
Data::Matrix<double> psc_t::V;

std::vector<std::string> psc_t::vname;
Data::Vector<double> psc_t::vmean;
Data::Vector<double> psc_t::vsd;

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
  // Save projection?
  //

  std::string projection = param.has( "proj" ) ? param.value( "proj" ) : "" ;

  //
  // Output input matrix too?
  //

  bool output_input = param.has( "output-input" );
  
  //
  // PSC parameters
  //
  
  int nc = param.has( "nc" ) ? param.requires_int( "nc" ) : 10 ;

  std::vector<double> th;
  if ( param.has( "th" ) ) th = param.dblvector( "th" );
  
  const int q = param.has( "q" ) ? param.requires_int( "q" ) : 5 ; 

  const bool standardize_inputs = param.has( "norm" ) ;
		 
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

      logger << "  reading from " << infile << "\n";
      
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
      if ( ch1 == ch2 ) Helper::halt( "require either CH or CH1 & CH2 in " + infile );
      if ( cols.find( "F" ) == cols.end() ) Helper::halt( "no F column in " + infile );

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
	  if ( vars.find( tok[i] ) != vars.end() ) slot2var[i] = tok[i];
	}
      
      if ( slot2var.size() == 0 ) Helper::halt( "no variables v=<...> in " + infile );

      const int ncols = tok.size();
      
      //
      // all set, not start reading
      //

      while ( ! IN1.eof() )
	{
	  // looking for ID F CH         --> PSD 
	  //  OR         ID F CH1 CH2    --> LCOH (default)  

	  std::string line;

	  Helper::safe_getline( IN1 , line );

	  std::vector<std::string> tok = Helper::parse( line , "\t" );

	  if ( IN1.eof() || tok.size() == 0 ) continue;
	  
	  if ( tok.size() != ncols ) Helper::halt( "incorrect number of columns in " + infile );

	  std::string id = tok[ id_slot ];

	  // skip if person not on an include list
	  if ( id_includes.size() != 0 && id_includes.find( id ) == id_includes.end() ) continue;

	  // skip if person is on an exclude list
	  if ( id_excludes.size() != 0 && id_excludes.find( id ) != id_excludes.end() ) continue;	  

	  std::string ch = "";

	  if ( ch_slot == -1 )
	    ch = tok[ ch1_slot ] + "." + tok[ ch2_slot ];
	  else
	    ch = tok[ ch_slot ];

	  // store as string and also numeric (for output)
	  std::string f = tok[ f_slot ] ; 

	  if ( flwr > 0 || fupr > 0 )
	    {
	      double fn;
	      if ( ! Helper::str2dbl( f , &fn ) )
		Helper::halt( "problem with frequency value: " + f );
	      // skip this frequency?
	      if ( flwr > 0 && fn < flwr ) continue;
	      if ( fupr > 0 && fn > fupr ) continue;
	    }
		  
	  std::map<int,std::string>::const_iterator ii = slot2var.begin();
	  while ( ii != slot2var.end() )
	    {
	      double x;
	      if ( ! Helper::str2dbl( tok[ ii->first ] , &x ) ) Helper::halt( "bad value in " + infile );

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

  logger << "  converting input spectra to a matrix\n";
  
  std::map<std::string,std::map<std::string,std::map<std::string,int> > > slot;
  std::map<std::string,std::string> col2ch, col2var;
  std::map<std::string,double> col2f;
  std::set<std::string> rows, cols;
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

  
  //
  // Populate matrix
  //
  
  Data::Matrix<double> U( rows.size() , cols.size() );
  
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
	      int row = 0;
	      std::set<std::string>::const_iterator ii = rows.begin();
	      while ( ii != rows.end() )
		{

		  // find channel?
		  const std::map<std::string,std::map<std::string,std::map<std::string,double> > > & dat = i2c2f2v.find( *ii )->second;
		  if ( dat.find( ss1->first ) == dat.end() ) Helper::halt( "no channel " + ss1->first + " for individual " + *ii );

		  // find frequency?
		  const std::map<std::string,std::map<std::string,double> > & dat2 = dat.find( ss1->first )->second;
		  if ( dat2.find( ss2->first ) == dat2.end() ) Helper::halt( "no frequency " + ss2->first + " for individual " + *ii );
		    		    
		  // find variable?
		  const std::map<std::string,double> & dat3 = dat2.find( ss2->first )->second;
		  if ( dat3.find( ss3->first ) == dat3.end() ) Helper::halt( "no variable " + ss3->first + " for individual " + *ii );
		  
		  // all okay, store
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

  logger << "  good, all expected observations found, no missing data\n";

  //
  // Check for invariant columns
  //

  Data::Vector<double> colvars = Statistics::variance( U );

  for (int i=0;i<colvars.size(); i++)
    if ( colvars[i] < EPS ) Helper::halt( "invariant column in input\n" );
					    
  
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
	  // this sets 'inc' values to missing, but uses the same prior for all channels
	  int removed = MiscMath::outliers( U.col(j).data_pointer() , th[t] , &inc , &prior);
	  //logger << "  removing " << removed << " var " << j << " round " << t << "\n";
	}
      
    }

  
  Data::Matrix<double> U2 = U;
  std::vector<std::string> id2 = id;

  U.clear();
  id.clear();

  ni = 0;
  for (int i=0;i<inc.size();i++)
    if ( inc[i] ) ++ni;

  logger << "  after outlier removal, " << ni << " individuals remainingin\n";

  U.resize( ni , nv );
  id.resize( ni );
  int r = 0;
  for (int i=0;i<U2.dim1();i++)
    if ( inc[i] )
      {
	for (int j=0;j<nv;j++) U(r,j) = U2(i,j);
	id[r] = id2[i];
	++r;
      }


  //
  // Output input data ?
  //

  Data::Matrix<double> X;

  // track means / SDs
  vmean = Statistics::mean( U );
  vsd = Statistics::sdev( U , vmean );
  
  if ( output_input ) X = U;
  
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
      Statistics::standardize( U );
    }
  else
    {
      Statistics::mean_center_cols( U );
    }

  
  Data::Vector<double> mm = Statistics::mean( U );
  
  W.clear(); V.clear();
  W.resize( nv );
  V.resize( nv , nv );

  bool okay = Statistics::svdcmp( U , W , V );
  if ( ! okay ) Helper::halt( "problem with SVD" );
  
  int rank = Statistics::orderSVD( U , W , V );
  if ( rank == 0 ) Helper::halt( "problem with input data, rank 0" );


  //
  // Output
  //

  // W

  double wsumsq = 0;
  for (int j=0;j<nv;j++) wsumsq += W[j] * W[j];

  writer.id( "." , "." );
  
  for (int j=0;j<nv;j++)
    {
      writer.level( j+1 , "I" );
      writer.value( "W" , W[j] );
      writer.value( "VE" , (W[j]*W[j])/wsumsq );
    }
  writer.unlevel( "W" );
  
  // V (transpose) (only first 'nc' components
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

  // data
  if ( output_input ) 
    {
      for (int i=0;i<ni;i++)
	{
	  writer.id( id[i] , "." );
	  for (int j=0;j<nv;j++)
	    {
	      writer.level( j+1 , "VAR" );
	      writer.value( "X" , X(i,j) );
	    }
	  writer.unlevel( "VAR" );
	}
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
	     << " " << vmean[j]
	     << " " << vsd[j];
      OUT1 << "\n";

      // means/SDs in reference population
      
      // components
      OUT1 << "NC: " << nc << "\n";
      
      OUT1 << "W:";
      for (int i=0;i<nc;i++)
	OUT1 << " " << W[i];
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

  if ( W.dim1() != 0 ) return;

  //
  // Read W and V matrices in this file
  //

  std::string infile = param.requires( "proj" );

  if ( ! Helper::fileExists( infile ) )
    Helper::halt( "could not find " + infile );
  
  logger << "  reading projection from " << infile << "\n";
  
  std::ifstream IN1( infile.c_str() , std::ios::in );

  // variables
  std::string dummy;
  int nv;
  
  IN1 >> dummy >> nv;
  vname.resize( nv );
  vmean.resize( nv );
  vsd.resize( nv );
  
  for (int j=0;j<nv;j++)
    IN1 >> vname[j] >> vmean[j] >> vsd[j];
    
  // components
  IN1 >> dummy >> nc;
      

  // all PSC (0) or a subset?

  logger << "  found " << nc << " PSCs based on " << nv << " variables\n";
  
  int k = param.has( "nc" ) ? param.requires_int( "nc" ) : 0 ;   

  if ( k >= 1 && k < nc )
    logger << "  subsetting to the first " << k << " of " << nc << " PSCs\n";
  else if ( k > nc )
    {
      logger << "  requested " << k << " PSCs but only " << nc << " present\n";
      k = nc;
    }
  else
    k = nc;
  
  W.resize( k );
  V.resize( nv , k );

  // buffer
  double inp;

  // W 
  IN1 >> dummy;
  for (int i=0;i<nc;i++)
    {
      IN1 >> inp;
      if ( i < k ) W[i] = inp;
    }
  
  // V
  IN1 >> dummy;
  for (int i=0;i<nv;i++)
    for (int j=0;j<nc;j++)
      {
	IN1 >> inp;
	if ( j < k ) V(i,j) = inp;
      }
  
  IN1.close();


  // set actual 'nc' to the desired value (i.e. if less than nc in file)
  nc = k;
  
  // reformat of W for projection
  DW.resize( nc , nc );
  for (int i=0; i<nc; i++)
    DW(i,i) = 1.0 / W[i];

  
  if ( 0 )
    {
      std::cout << "W\n" << W.print() << "\n";
      std::cout << "V\n" << V.print() << "\n";
    }
  
  
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

  // TODO: allow multiple cache here?
  std::string cache_name = param.requires( "cache" );

  if ( ! edf.timeline.cache.has_num( cache_name ) )
    Helper::halt( "cache not found for this individual: " + cache_name );

  cache_t<double> * cache = edf.timeline.cache.find_num( cache_name );

  // see which variables exist, i.e. psc_t::vname[] 

  int nv = vname.size();

  Data::Matrix<double> X( 1 , nv );
  
  for (int i=0; i<nv; i++)
    {
      std::vector<std::string> tok = Helper::parse( vname[i] , "~" );
      if ( tok.size() != 3 ) Helper::halt( "bad format for PSC vnames, expecting 3 fields, '~'-delimited" );
      std::string ch = tok[0];
      double f;
      if ( ! Helper::str2dbl( tok[1] , &f ) ) Helper::halt( "bad frequency value in PSC vname" );
      std::string var = tok[2];
      
      // can we find this value in the cache?

      ckey_t key( var );
      key.add( "CH" , ch );
      key.add( "F" , f );

      std::vector<double> cx = cache->fetch( key );

      if ( cx.size() != 1 )
	Helper::halt( "could not find cached variable: " + vname[i] );

      X(0,i) = cx[0];
    }
  
  logger << "  all " << nv << " features found in the cache\n";

  
  //
  // Mean center data
  //
  
  Statistics::subtract_cols( X , vmean );
  
  
  //
  // Project given W and V to get U, PSCs for this individual
  //
  
  Data::Matrix<double> U_proj = X * V * DW;


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
