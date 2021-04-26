
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
#include "fftw/fftwrap.h"
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
  // Options
  //

  bool do_not_add_channels = param.has( "no-new-channels" );

  std::string kop_tag = param.has( "tag" ) ? param.value( "tag" ) : "KOP" ; 

  double fmin = param.has( "min" ) ? param.requires_dbl( "min" ) : 0.5; 
  double fmax = param.has( "max" ) ? param.requires_dbl( "max" ) : 20 ; 
  

  //
  // Epoch-wise FFT or full HT? 
  //
  
  bool do_ht = param.has( "w" );

  bool do_fft = param.has( "fft" );
  
  
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
  // First pass FFT (epoch-wise) 
  //
  
  if ( do_fft )
    {
      edf.timeline.first_epoch();
 
      // set up FFT
      slice_t slice( edf , signals(0) , edf.timeline.epoch( 0 ) );
      const std::vector<double> * d = slice.pdata();
      int index_length = d->size();
      real_FFT fftseg( index_length , index_length , sr , WINDOW_NONE );
      fftseg.apply( &((*d)[0]) , index_length );
      int my_N = fftseg.cutoff;
      
      while ( 1 ) 
	{
	  int epoch = edf.timeline.next_epoch();
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
	  Eigen::MatrixXd ekop = Eigen::MatrixXd::Zero( my_N , ns );
	  
	  writer.epoch( edf.timeline.display_epoch( epoch ) );
	  
	  // all channels
	  
	  for ( int s=0; s<ns; s++ )
	    {
	      
	      slice_t slice( edf , signals(s) , interval );
	      
	      const std::vector<double> * d = slice.pdata();
	      
	      if ( index_length != d->size() ) Helper::halt( "internal error in sync() " );
	      
	      fftseg.apply( &((*d)[0]) , index_length );
	      
	      // Extract the raw transform
	      std::vector<std::complex<double> > t = fftseg.transform();

	      if ( my_N != fftseg.cutoff ) Helper::halt( "internal error in sync() " );
	      
	      // Extract the raw transform scaled by 1/n
	      //	  std::vector<std::complex<double> > t2 = fftseg.scaled_transform();
	      
	      for (int f=0;f<my_N;f++)
		ekop( f , s ) = std::arg( t[f] );
	      
	      // next signal
	    }
      
	  // get KOP for this 
	  Eigen::ArrayXd skop = Eigen::ArrayXd::Zero( my_N );
      
	  for (int f=0; f<my_N; f++)
	    {
	      
	      if ( fftseg.frq[f] >= fmin && fftseg.frq[f] <= fmax )
		{
		  dcomp k( 0, 0 );
		  
		  for (int s=0; s<ns; s++)
		    k += exp( dcomp ( 0 , ekop(f,s) ) );
		  k /= (double)ns;
		  skop[f] = abs( k );
		  
		  writer.level( fftseg.frq[f] , globals::freq_strat );
		  writer.value( "KOP" , skop[f] );
		}
	    }
	  writer.unlevel( globals::freq_strat );
	  
	  // next epoch
	}
      writer.unepoch();
      
    }


  //
  // filter-Hilbert approach 
  //


  if ( ! do_ht ) return;


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
  else
    Helper::halt( "no frequency bins specified" );



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
      
      logger << "  filter-Hilbert " << lwr[f] << "-" << upr[f] << "Hz:";

      //
      // filter-Hilbert 
      //
      
      double tw = param.has( "tw" ) ? param.requires_dbl( "tw" ) : 1 ;

      double ripple = param.has( "ripple" ) ? param.requires_dbl( "ripple" ) : 0.025;
      
      Eigen::MatrixXd phase = Eigen::MatrixXd::Zero( rows , cols );
      
      for (int s = 0; s < cols; s++)
	{
	  logger << ".";
	  
	  const std::vector<double> d = eigen_ops::copy_vector( X.col(s) );
	  
	  std::vector<double> p;
	  
	  run_hilbert( d , sr , lwr[f] , upr[f] , ripple , tw , NULL , &p , NULL , NULL );
	  
	  phase.col(s) = Eigen::Map<Eigen::VectorXd>(p.data(), p.size() );
	  
	}
     
      logger << "\n";
      
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
  
  //
  // Add new signals
  //
  
  if ( ! do_not_add_channels )
    {
      logger << "  adding " << nf << " new signals to EDF:";
      for (int c=0;c<nf;c++)
	{
	  std::vector<double> copy( rows );
	  Eigen::VectorXd::Map( &copy[0], rows ) = kop.col(c);
          logger << " " << kop_tag + Helper::int2str( c+1 ) ;
          edf.add_signal( kop_tag + Helper::int2str( c+1 ) , sr , copy );
        }
      logger << "\n";
    }

  //
  // Epoch-level summaries
  //

  edf.timeline.ensure_epoched();

  const uint64_t epoch_sp = sr * edf.timeline.epoch_length();
  
  const int ne = edf.timeline.num_epochs();

  const int expected_ne = rows / epoch_sp;
  
  if ( ne != expected_ne ) 
    logger << "  warning : expecting " << expected_ne << " but found " << ne << "\n";

  
  uint64_t pos = 0;
  int epoch = 0;
  for (int pos = 0 ; pos < rows ; pos += epoch_sp )
    {
      
      writer.epoch( edf.timeline.display_epoch( epoch ) );
      
      for (int f=0; f<nf; f++)
	{
	  
	  double m = 0 ; 
	  uint64_t pos2 = pos + epoch_sp ; 
	  for (int p = pos ; p < pos2; p++ )
	    m += kop(p,f);
	  m /= epoch_sp;
	  
	  writer.level( f+1 , globals::freq_strat );
	  writer.value( "KOP" , m );
	}
      writer.unlevel( globals::freq_strat );
      ++epoch;
    }
  writer.unepoch();
  
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

