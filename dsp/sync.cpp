
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

#include "dsp/sync.h"

#include "dsp/hilbert.h"
#include "dsp/wrappers.h"
#include "stats/eigen_ops.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger;

void dsptools::sync( edf_t & edf , param_t & param )
{

  //
  // Get signals
  //

  std::string signal_label = param.requires( "sig" );
  
  const bool no_annotations = true;
  
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );  
  
  const int ns = signals.size();

  if ( ns < 2 ) return;

  //
  // Get frequencies
  //

  std::vector<double> lwr, upr;

  if ( param.has( "f-lwr" ) && param.has( "f-upr" ) && param.has( "w" ) && param.has( "r" ) )
    {
      double w = param.requires_dbl( "w" );
      double r = param.requires_dbl( "r" );
      double fl = param.requires_dbl( "f-lwr" );
      double fu = param.requires_dbl( "f-upr" );
      for (double ff = fl ; ff <= (fu+0.5*r) ; ff += r )
	{
	  if ( ff - w/2.0 > 0 )
	    {
	      lwr.push_back( ff - w/2.0 );
	      upr.push_back( ff + w/2.0 );
	    }
	}
    }
  else if ( param.has( "f-lwr" ) && param.has( "f-upr" ) )
    {
      lwr = param.dblvector( "f-lwr" );
      upr = param.dblvector( "f-upr" );
      if ( lwr.size() != upr.size() ) 
	Helper::halt( "f-lwr and f-upr have different sizes" );
      for (int i=0;i<lwr.size();i++)
	if ( lwr[i] >= upr[i] ) Helper::halt( "f-lwr >= f-upr" );
    }
  else if ( param.has( "f" ) )
    {
      lwr = upr = param.dblvector( "f" );
      double w = param.has( "w" ) ? param.requires_dbl( "w" ) : 3 ; // plus/minus 3 Hz by default
      for (int i=0; i<lwr.size(); i++) 
	{
	  lwr[i] -= w;
	  upr[i] += w;
	  if ( lwr[i] <= 0 ) 
	    Helper::halt( "frequency below 0 Hz specified" );
	}
    }
  

  bool has_freqs = lwr.size() > 0 ;

  
  //
  // Check sample rates
  //  
  
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

  const Eigen::MatrixXd & X = mslice.data_ref();  

  const int rows = X.rows();
  const int cols = X.cols();


  const int nf = lwr.size();

  Eigen::MatrixXd kop = Eigen::MatrixXd::Zero( rows , nf );

  for (int f=0; f<nf; f++)
    {

      //
      // filter-Hilbert 
      //
      
      double tw = param.has( "tw" ) ? param.requires_dbl( "tw" ) : 1 ;

      double ripple = param.has( "ripple" ) ? param.requires_dbl( "ripple" ) : 0.025;
      
      Eigen::MatrixXd phase = Eigen::MatrixXd::Zero( rows , cols );
      
      for (int s = 0; s < cols; s++)
	{
	  
	  const std::vector<double> d = eigen_ops::copy_vector( X.col(s) );
	  
	  std::vector<double> p;
	  
	  run_hilbert( d , sr , lwr[f] , upr[f] , ripple , tw , NULL , &p , NULL , NULL );
	  
	  phase.col(s) = Eigen::Map<Eigen::VectorXd>(p.data(), p.size() );
	  
	}
     
      
      //
      // Kuramoto order parameter
      //
      
      for (int i=0; i<rows; i++)
	{
	  dcomp k( 0, 0 );
	  for (int j=0; j<cols; j++)
	    k += exp( dcomp ( 0 , phase(i,j) ) );
	  k /= (double)cols;
	  kop(i,f) = abs( k );
	}

    }
  
  std::cout << "KOP\n" << kop << "\n";

  //
  // Group calcs  (region vars) 
  //

  //
  // Permuation / time-shuffle for original values 
  //


  //
  // Statistics on KOP
  //

  
}  

