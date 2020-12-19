
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

#include "dsp/microstates.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/logger.h"
#include "stats/kmeans.h"


extern writer_t writer;
extern logger_t logger;

void dsptools::microstates( edf_t & edf , param_t & param )
{

  const bool no_annotations = true;
  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , no_annotations );  

  const int ns = signals.size();

  // Check sample rates
  
  if ( ns < 2 ) return;
  
  const int sr = edf.header.sampling_freq( signals(0) );
  
  for (int i=1;i<ns;i++)
    {      
      if ( edf.header.sampling_freq( signals(i) ) != sr ) 
	Helper::halt( "all signals must have similar SR for MICROSTATES" );
    }

  //
  // Microstate analysis parameters
  //
  
  std::vector<int> ks = param.intvector( "k" );

  
  //
  // Fetch sample matrix
  //
  
  matslice_t mslice( edf , signals , edf.timeline.wholetrace() );

  const Data::Matrix<double> & X = mslice.data_ref();
  
  //
  // Perform basic microstate analysis
  //
  
  microstates_t mstates( X , signals , ks );

  
  //
  // Get solutions
  //
  
  
}

  

microstates_t::microstates_t( const Data::Matrix<double> & X , 
			      const signal_list_t & signals , 
			      const std::vector<int> & ks )
{

  //
  // Normalize within channel
  //

  Data::Matrix<double> Z = X;
  
  std::cout << " channel SDs " << Statistics::sdev( X , Statistics::mean(X) ).print() << "\n";
  
  Statistics::standardize( Z );

  //
  // Global field power 
  //

  const int np = X.dim1();
  const int nc = X.dim2();

  logger << "  calculating GFP for sample\n";

  Data::Vector<double> GFP( np );

  for (int i=0; i<np; i++)
    {

      // get time-points across channels
      Data::Vector<double> p = Z.row( i );

      // ignore polarity
      for (int j=0; j<nc; j++) p[j] = fabs( p[j] );

      GFP[i] = sqrt( Statistics::variance( p ) );
      
      //      std::cout << GFP[i] << "\n";
   }
  
  logger << "  finding GFP peaks\n";

  //
  // Find peaks in GFP
  //

  std::vector<bool> peak( np , false );
  std::vector<int> peak_idx;
  int n_peaks = 0;
  for (int i=1; i<(np-1); i++)
    {
      if ( GFP[i] > GFP[i-1] && GFP[i] > GFP[i+1] ) 
	{
	  peak[i] = true;
	  peak_idx.push_back(i);
	  ++n_peaks;
	}
    }

  //
  // Copy subset of data for clustering (nb. ignores signal polarity again here)
  //

  Data::Matrix<double> P( n_peaks , nc );
  for (int r=0; r<n_peaks; r++)
    for (int c=0; c<nc; c++)
      P( r , c ) = fabs( Z( peak_idx[r] , c ) );

  logger << "  extracted " << n_peaks << " peaks from " << np << " samples ("
	 << round( 100 * ( n_peaks / (double)np ) ) << "%)\n";

  std::cout << "P\n" << P.print() << "\n";

  //
  // K-means clustering on peaks
  //

  for (int ki = 0 ; ki < ks.size(); ki ++ )
    {

      const int k = ks[ki];

      std::cout << "K = " << k << "\n";

      kmeans_t kmeans;
      std::vector<int> sol;
      Data::Matrix<double> means = kmeans.kmeans( P , k , &sol );
      
      std::map<int,int> cnts;
      for (int i=0; i<sol.size(); i++) cnts[ sol[i] ]++;
      std::map<int,int>::const_iterator cc = cnts.begin();
      while ( cc != cnts.end() )
	{
	  std::cout << "cnts = " << cc->first << "\t" << cc->second << "\t"
		    << cc->second / (double)n_peaks << "\n";
	  ++cc;
	}

      for (int s=0;s<signals.size();s++) std::cout << signals.label(s) << "\n";

      means = Statistics::transpose( means );
      Statistics::standardize( means );
      std::cout << "M\n" << means.print() << "\n\n";
      
    }
   
}

     
Data::Vector<int> microstates_t::solution( const int k )
{
  Data::Vector<int> retval;
  return retval;
}


