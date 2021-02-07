
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

#include "ica-wrapper.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "eval.h"
#include "db/db.h"

#include <map>
#include <set>

extern writer_t writer;

extern logger_t logger;

void dsptools::ica_wrapper( edf_t & edf , param_t & param )
{

  std::string signal_label = param.requires( "sig" );

  std::string component_tag = param.has( "tag" ) ? param.value( "tag" ) : "IC_";

  std::string Aout = param.requires( "A" );

  bool write_matrix = param.has( "file" );
  
  bool write_S_matrix = param.has( "S" );
  
  std::string matrix_fileroot = write_matrix ? param.value( "file" ) : "xxx";

  bool original_signals = param.has( "original-signals" );

  bool do_not_add_channels = param.has( "no-new-channels" );

  //
  // Fetch all data signals from 'sig'
  //

  const bool no_annotations = true;
  
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );  

  const int ns = signals.size();


  //
  // Check sample rates
  //
  
  if ( ns < 2 ) return;
  
  const int sr = edf.header.sampling_freq( signals(0) );
  
  for (int i=1;i<ns;i++)
    {      
      if ( edf.header.sampling_freq( signals(i) ) != sr ) 
	Helper::halt( "all signals must have similar SR for ICA" );
    }


  //
  // Fetch sample matrix
  //
  
  eigen_matslice_t mslice( edf , signals , edf.timeline.wholetrace() );
  
  // nb. fastICA() this will modify X/mslice...
  Eigen::MatrixXd & X = mslice.nonconst_data_ref();

  const int rows = X.rows();
  const int cols = X.cols();

  //
  // Number of components 'nc' parameter
  //

  int nc = param.has( "nc" ) ? param.requires_int( "nc" ) : ns ;


  //
  // ICA: note, this alters input 'X'
  //
  
  eigen_ica_t ica( X , nc );


  //
  // Add new signals
  //
  
  if ( ! do_not_add_channels ) 
    {
      logger << "  adding " << nc << " new signals to EDF:";
    
      for (int c=0;c<nc;c++)
	{
	  std::vector<double> copy( rows );
	  Eigen::VectorXd::Map( &copy[0], rows ) = ica.S.col(c);
	  logger << " " << component_tag + Helper::int2str( c+1 ) ;
	  edf.add_signal( component_tag + Helper::int2str( c+1 ) , sr , copy );	  
	}
      logger << "\n";
    }
  
  
  
  //
  // Output mixing matrix, etc
  //

  // K : cols x nc
  // A : nc x cols
  // W : nc x nc
  
  for (int i=0;i<nc;i++)
    {
      writer.level( i+1 , "IC1" );
      for (int j=0;j<nc;j++)
	{
	  writer.level( j+1 , "IC2" );	  
	  writer.value( "W" , ica.W(i,j) );
	}
      writer.unlevel( "IC2" );
    }
  writer.unlevel( "IC1" );


  for (int i=0;i<nc;i++)
    {
      writer.level( i+1 , "IC" );	  
      
      for (int j=0;j<cols;j++)
	{
	  writer.level( signals.label(j) , globals::signal_strat );
	  writer.value( "A" , ica.A(i,j) );
	}
      writer.unlevel( globals::signal_strat );      
    }
  writer.unlevel( "IC" );
  

  // use different labels, so that K and A are in separate strata
  for (int i=0;i<cols;i++)
    {
      writer.level( signals.label(i) , "K" + globals::signal_strat );
      for (int j=0;j<nc;j++)
	{
	  writer.level( j+1 , "KIC" );	  
	  writer.value( "K" , ica.K(i,j) );
	}
      writer.unlevel( "KIC" );
    }
  writer.unlevel( "K" + globals::signal_strat );
  

  //
  // File-based output
  //
  
  if ( write_matrix )
    {
      

      if ( write_S_matrix )
	{
	  	  
	  std::ofstream S( (matrix_fileroot + "S").c_str() , std::ios::out );
	  for (int j=0;j<nc;j++) S << ( j ? "\t" : "" ) << "S" << j+1;
	  S << "\n";  
	  for (int i=0;i<rows;i++)
	    {
	      for (int j=0;j<nc;j++) S << ( j ? "\t" : "" ) << ica.S(i,j) ;      
	      S << "\n";
	    }
	  S.close();
	}

      if ( original_signals ) 
	{
	  std::ofstream F( (matrix_fileroot + "X").c_str() , std::ios::out );
	  for (int j=0;j<cols;j++) F << ( j ? "\t" : "" ) << "X" << j+1;
	  F << "\n";
	  for (int i=0;i<rows;i++)
	    {
	      for (int j=0;j<cols;j++) F << ( j ? "\t" : "" ) << X(i,j);
	      F << "\n";
	    }
	  F.close();
	}

      //
      // other matrices
      //

      // K : cols x nc   pre-whitening matrix: data --> principal components
      // A : nc x col    un-mixing matrix
      // W : nc x nc     mixing matrix
      // S : rows x cols : IC time courses
      
      std::ofstream K( (matrix_fileroot + "K").c_str() , std::ios::out );
      for (int i=0;i<cols;i++)
	{
	  for (int j=0;j<nc;j++) K << ( j ? "\t" : "" ) << ica.K(i,j);
	  K << "\n";
	}
      K.close();
      
      std::ofstream W( (matrix_fileroot + "W").c_str() , std::ios::out );
      for (int i=0;i<nc;i++)
	{
	  for (int j=0;j<nc;j++) W << ( j ? "\t" : "" ) << ica.W(i,j);
	  W << "\n";
	}
      W.close();
            
      std::ofstream A( (matrix_fileroot + "A").c_str() , std::ios::out );
      for (int i=0;i<nc;i++)
	{
	  for (int j=0;j<cols;j++) A << ( j ? "\t" : "" ) << ica.A(i,j);
	  A << "\n";
	}
      A.close();

    }


  std::ofstream A( Aout.c_str() , std::ios::out );
  for (int i=0;i<nc;i++)
    for (int j=0;j<cols;j++) 
      A << component_tag + Helper::int2str( i+1 )  << "\t"
	<< signals.label(j) << "\t"
	<< ica.A(i,j) << "\n";
  A.close();
  
}



void dsptools::ica_adjust( edf_t & edf , param_t & param )
{
  
  // requires an un-mixing matrix A 
  std::string Af = param.requires( "A" );
  if ( ! Helper::fileExists( Af ) ) 
    Helper::halt( "could not find matrix A: " + Af );
  
  std::map<std::string,std::map<std::string,double> > Am;
  int cnt = 0;
  std::set<std::string> ics, chs;
  std::ifstream Ain( Af.c_str() , std::ios::in );
  while ( ! Ain.eof() )
    {
      std::string ic , ch;
      double value;
      Ain >> ic >> ch >> value;
      if ( Ain.eof() ) break;
      Am[ ic ][ ch ] = value;
      ics.insert( ic );
      chs.insert( ch );
      ++cnt;
    }
  Ain.eof();  

  logger << "  read " << ics.size() << " ICs, " << chs.size() << " channels from " << Af << "\n";

  // check squared off matrix
  if ( ics.size() * chs.size() != cnt )
    Helper::halt( "problem with format of " + Af );

  // requires IC tag, i.e. IC_1, IC_2 etc by default
  std::string tag = param.has( "tag" ) ? param.value( "tag" ) : "IC_" ; 

  // signals to be adjusted
  signal_list_t signals = edf.header.signal_list( param.value( "sig" ) );
  
  // (putative) ICs to adjust for (i.e. subtract) i.e. IC_1, IC_2, 
  signal_list_t adjs = edf.header.signal_list( param.value( "adj" ) );

  // nothing to do
  if ( signals.size() == 0 || adjs.size() == 0 ) return;

  const int ns = signals.size();
  const int na = adjs.size();

  // check all requested signals exist in A
  for (int s=0; s<ns; s++)
    if ( chs.find( signals.label(s) ) == chs.end() ) 
      Helper::halt( "could not find " + signals.label(s) + " in " + Af );

  // check all requested IC exist i Af
  for (int a=0; a<na; a++)
    if ( ics.find( adjs.label(a) ) == ics.end() )
      Helper::halt( "could not find " + adjs.label(a) + " in " + Af );

  // track SR
  const int sigsr = (int)edf.header.sampling_freq( signals(0) );

  //
  // correlative factors?
  //

  bool corr_criteria = param.has( "corr-sig" );
  std::map<std::string,double> corr;
  signal_list_t corrsigs;
  if ( corr_criteria )
    {
      if ( ! param.has( "corr-th" ) ) Helper::halt( "requires corr-th with corr-sig" );
      std::vector<std::string> tok = param.strvector( "corr-sig" );
      std::vector<double> tokd = param.dblvector( "corr-th" );
      if ( tok.size() != tokd.size() )
	Helper::halt( "corr-sig and corr-th not of similar lengths" );
      for (int i=0; i<tok.size(); i++) corr[ tok[i] ] = tokd[i] ;
      corrsigs = edf.header.signal_list( param.value( "corr-sig" ) );      
    }

  //
  // fetch adjusting factors
  //
    
  std::vector<std::vector<double> > adjdata(na); // sig x points

  for (int a=0;a<na;a++)
    if ( edf.header.n_samples[ adjs(a) ] != edf.header.n_samples[ signals(0) ] )
      Helper::halt( "different SRs, not to RESAMPLE first" );

  int rec = edf.timeline.first_record();
  while ( rec != -1 )
    {      
      edf.ensure_loaded( rec );
      edf_record_t & record = edf.records.find(rec)->second;            
      // get data
      for (int a=0;a<na;a++) 		
	{	  
	  std::vector<double> t = record.get_pdata( adjs(a) );	  
	  for (int j=0;j<t.size();j++)
	    adjdata[a].push_back( t[j] );
	}
      // next record
      rec = edf.timeline.next_record(rec);       
    }

    
  //
  // fetch any correlation factors
  //

  const int nc = corrsigs.size();

  std::vector<std::vector<double> > corrdata(nc);

  
  for (int c=0;c<nc;c++)
    if ( edf.header.n_samples[ corrsigs(c) ] != edf.header.n_samples[ signals(0) ] )
      Helper::halt( "different SRs, not to RESAMPLE first" );
  
  rec = edf.timeline.first_record();
  while ( rec != -1 )
    {
      edf.ensure_loaded( rec );
      edf_record_t & record = edf.records.find(rec)->second;      
      // get data
      for (int c=0;c<nc;c++)
	{
	  std::vector<double> t = record.get_pdata( corrsigs(c) );
	  for (int j=0;j<t.size();j++)
	    corrdata[c].push_back( t[j] );
	}
      // next record
      rec = edf.timeline.next_record(rec);       
    }
  
  
  //
  // fetch target (to-be-adjusted) signal data
  //
  
  std::vector<std::vector<double> > sigdata(ns);
  
  for (int s=0;s<ns;s++)
    if ( edf.header.n_samples[ signals(s) ] != edf.header.n_samples[ signals(0) ] )
      Helper::halt( "different SRs, not to RESAMPLE first" );
  
  rec = edf.timeline.first_record();
  while ( rec != -1 )
    {
      edf.ensure_loaded( rec );
      edf_record_t & record = edf.records.find(rec)->second;      
      // get data
      for (int s=0;s<ns;s++)
	{
	  std::vector<double> t = record.get_pdata( signals(s) );
	  for (int j=0;j<t.size();j++)
	    sigdata[s].push_back( t[j] );
	}
      // next record
      rec = edf.timeline.next_record(rec);       
    }

    
  //
  // now we have all data in:
  //   sigdata
  //   adjdata
  //   corrdata
  //  go through signals and make any adjustments needed
  //
  
  logger << "  adjusting " << ns << " signals based on " << na << " adjustment-signals\n";

  //
  // Identify only adjustment factors that are signficiantly (time-domain)
  // correlated with one of more other channels
  //
  
  std::set<std::string> excludes;
  if ( corrsigs.size() != 0 )
    {
      for (int a=0;a<na; a++)
	{
	  bool flagged = false;
	  for (int c=0;c<nc; c++)
	    {
	      double correl = Statistics::correlation( adjdata[a] , corrdata[c] );
	      // threshold for this corr signal
	      std::string label = corrsigs.label(c);
	      double threshold = corr[ label ];
	      std::cout << "checking " << adjs.label(a) << "  -- " << corrsigs.label(c) << " = " << correl << " " << threshold << "\n";
	      if ( fabs( correl ) > threshold )
		{
		  logger << "   including " << adjs.label(a) << " based on its absolute correlation with "
			 << corrsigs.label(c) << ", r = " << correl << "\n";
		  flagged = true;
		  break;
		}
	    }
	  // was this factor flagged? (i.e. based on strength of correlation with at least one corr-sig ? (above corr-th)
	  if ( ! flagged )
	    excludes.insert( adjs.label(a) );
	}
    }

  
  if ( nc > 0 )
    logger << "  " << na-excludes.size() << " adjustment-signals retained based on correlations with " << nc << " correlative-signals\n";


    
  //
  // Process each signal
  //
  
  for (int s=0; s<ns; s++)
    {
      
      // check SR
      if ( (int)edf.header.sampling_freq( signals(s) ) != sigsr )
	Helper::halt( "must have similar sampling rates for ADJUST - use RESAMPLE first" );

      // reference to the original data
      std::vector<double> & d = sigdata[s];
      const int np = d.size();
      
            
      // adjust for each adjusting signal
      for (int a=0; a<na; a++)
	{
	  // skip?
	  if ( excludes.find( adjs.label(a) ) != excludes.end() )
	    continue;
	  
	  std::vector<double> & adj = adjdata[a];
	  if ( adj.size() != np )
	    Helper::halt( "internal error in edf_t::adjust()" );
	  
	  // A matrix value
	  double aval = Am[ adjs.label(a) ][  signals.label(s) ];

	  // adjust D <- D - S_i A_i
	  for (int p=0; p<np; p++)
	    d[p] -= adj[p] * aval;
	  
	}
      
      // update signal                                                                                                                                
      edf.update_signal( signals(s) , &d );

    } // next signal

  // all done
}


