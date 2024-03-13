
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

#include "slice.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "eval.h"

#include "timeline/timeline.h"
#include "fftw/fftwrap.h"
#include "defs/defs.h"
#include "edf.h"
#include "intervals/intervals.h"



//
// slice_t
//

interval_t slice_t::duration() const 
{ 
  // interval is 1-past end of interval
  return interval_t( time_points[0] , time_points[ time_points.size()-1 ] + 1LLU );
}


slice_t::slice_t( edf_t & edf , 
		  int signal ,
		  const interval_t & interval ,
		  int downsample , 
		  bool digital ,
		  bool get_smps )
  : edf(edf) , signal(signal) , interval(interval) , downsample(downsample) 
{

  //
  // Extract data for a single signal from an EDF structure, given a
  // mask (in timeline)
  //

  data.clear();
  time_points.clear();
  records.clear();
  dig_data.clear();
  smps.clear();
  
  //
  // Empty?
  //

  if ( interval.empty() ) return;
  
  //
  // Get signals and labels
  //
  
  if ( signal < 0 || signal >= edf.header.ns ) 
    {
      Helper::halt( "problem in slice(), bad signal requested: " 
		    + Helper::int2str(signal) 
		    + " of " + Helper::int2str( edf.header.ns ) );
    }
      

  // 
  // Populate data matrix
  //
  
  const uint64_t duration = interval.duration();
  
  //
  // use fixed channel/signal sampling rate (i.e. array can be ragged)
  //

  if ( ! digital ) 
    data = edf.fixedrate_signal( interval.start , 
				 interval.stop , 
				 signal , 
				 downsample , 
				 &time_points , 
				 &records ,
				 get_smps ? &smps : NULL , 
				 NULL );
  
  else
    data = edf.fixedrate_signal( interval.start , 
				 interval.stop , 
				 signal , 
				 downsample , 
				 &time_points , 
				 &records , 
				 get_smps ? &smps : NULL ,
				 &dig_data ); // requests digital value; here 'data' will be empty   
  

}
 


//
// mslice_t
//

mslice_t::mslice_t( edf_t & edf , 
		    const signal_list_t & signals , 
		    const interval_t & interval , 
		    int    downsample )
{
  channel.clear();
  labels.clear();
  const int ns = signals.size();
  for (int s=0;s<ns;s++)
    {
      slice_t * slice = new slice_t(edf , signals(s) , interval , downsample );
      channel.push_back( slice ); 
      labels.push_back( signals.label( s ) );
    }
}



Data::Matrix<double> mslice_t::extract()
{
  const int nr = channel[0]->size(); 
  const int nc = channel.size();
  
  Data::Matrix<double> d;
  for (int c=0;c<nc;c++)
    {
      if ( nr != channel[c]->size() ) 
	Helper::halt( "internal error in mslice, SRs different" );
      d.add_col( *channel[c]->pdata() );
    }
  return d;
}




//
// matslice_t : only for equal SR data, returns a single matrix
//


matslice_t::matslice_t( edf_t & edf , 
			const signal_list_t & signals , 
			const interval_t & interval )
{

  data.clear();
  
  time_points.clear();
  
  labels.clear();
  
  //
  // Empty?
  //

  const int ns = signals.size();
  
  if ( ns == 0 ) return;
  
  if ( interval.empty() ) return;
  
  
  //
  // Get labels and check SR
  //

  const int Fs = edf.header.n_samples[ signals(0) ];

  labels.push_back( signals.label(0) );
  
  for (int s=1;s<ns;s++)
    {
      if ( edf.header.n_samples[ signals(s) ] != Fs ) Helper::halt( "unequal sample rates in matslice_t: use RESAMPLE" );
      labels.push_back( signals.label(s) );
    }


  
  //
  // use fixed channel/signal sampling rate (i.e. array can be ragged), and populate time-points
  //
  
  data.add_col( edf.fixedrate_signal( interval.start , 
				      interval.stop , 
				      signals(0) ,  // first channel
				      1 ,  // no downsampling
				      &time_points ,  // get TPs for first channel only
				      NULL , // no recs
				      NULL , // no smps
				      NULL ) ); // not digital data
  
  
  //
  // get all other channels (w/out time-points)
  //
  
  for (int s=1;s<ns;s++) 
    data.add_col( edf.fixedrate_signal( interval.start , interval.stop , signals(s) , 1 , NULL , NULL , NULL , NULL ) );  
  
}
  



//
// eigen_matslice_t : only for equal SR data, returns a single matrix
//


eigen_matslice_t::eigen_matslice_t( edf_t & edf , 
				    const signal_list_t & signals , 
				    const interval_t & interval )
{

  
  data.resize(0,0);
  
  time_points.clear();
  
  labels.clear();

  //
  // Empty?
  //

  const int ns = signals.size();
  
  if ( ns == 0 ) return;
  
  if ( interval.empty() ) return;
  
  
  //
  // Get labels and check SR
  //

  const int Fs = edf.header.n_samples[ signals(0) ];

  labels.push_back( signals.label(0) );
  
  for (int s=1;s<ns;s++)
    {
      if ( edf.header.n_samples[ signals(s) ] != Fs )
	Helper::halt( "unequal sample rates in matslice_t: use RESAMPLE" );
      labels.push_back( signals.label(s) );
    }


  
  //
  // use fixed channel/signal sampling rate (i.e. array can be ragged), and populate time-points
  //

  std::vector<double> ch1 = edf.fixedrate_signal( interval.start , 
						  interval.stop , 
						  signals(0) ,    // first channel
						  1 ,             // no downsampling
						  &time_points ,  // get TPs for first channel only
						  NULL, NULL, NULL ) ;        // no records, smps or dig. data

  const int nr = ch1.size();

  // size output matrix, and assign first row:
  data.resize( nr , ns );

  data.col(0) = Eigen::VectorXd::Map( &ch1[0] , nr );
  
    
  
  //
  // get all other channels (w/out time-points)
  //
  
  for (int s=1;s<ns;s++) 
    {
      // TODO: add a version of fixedrate_signal() that directly outputs to Eigen matrix class
      std::vector<double> tmp = edf.fixedrate_signal( interval.start , interval.stop , signals(s) , 1 , NULL , NULL , NULL , NULL ) ;
      data.col(s) = Eigen::VectorXd::Map( &tmp[0] , nr );
    }
  
  
}



//
// Non-epoch based slicer
//


void edf_t::slicer( const std::set<interval_t> & intervals1 , param_t & param , int extract )
{
  
  // nothing to do
  if ( intervals1.size() == 0 ) return;
  
  // options
  bool flatten_overlap = ! param.has( "allow-overlap" );
  
  // only a single signal allowed
  std::string signal_label = param.requires( "sig" );   

  // optionally dump the raw signal?
  bool raw_signal = param.has( "dump-signal" );

  // window around segment? (i.e. 2 seconds, means 2 secs before, and after)
  double window = param.has("window") ? param.requires_dbl("window") : 0 ;
  uint64_t window_tp = window * globals::tp_1sec;
  
  // and define the number of windows (i.e. 2 * 2second)
  double n_window = param.has("n-window" ) ? param.requires_int( "n-window" ) : 1 ;
  uint64_t total_window_tp = window_tp * n_window ;

  // max timeline.last_time_point_tp;

  // Short-time FFT parameters

  const double fft_segment_sec = 1.0; // analyse in 1-second bins
  const double fft_segment_inc = 0.5; // overlap by 50% (0.5 seconds)

  // Attach signal

  signal_list_t signals = header.signal_list( signal_label );  
  if ( signals.size() > 1 ) Helper::halt( "only a single signal allowed for SLICE" );
  if ( header.is_annotation_channel(signals(0)) ) return;
  
  const int Fs = header.sampling_freq( signals(0) );

  int window_samples = window * Fs; 
  
  // check each interval will be within bounds
  std::set<interval_t> intervals0;
  std::set<interval_t>::const_iterator ii0 = intervals1.begin();

  if ( timeline.last_time_point_tp < total_window_tp ) return;

  uint64_t mx = timeline.last_time_point_tp - total_window_tp;
  while ( ii0 != intervals1.end() )
    {
      if ( ii0->start > total_window_tp && ii0->stop <= mx ) intervals0.insert( *ii0 );
      //else std::cout << "dropping " << ii0->start << " " << ii0->stop << "\n";
      ++ii0;
    }

  // process intervals

  std::set<interval_t> intervals;

  uint64_t expected = total_window_tp;

  if ( flatten_overlap )
    {

      interval_t w = *intervals0.begin();      
      w.expand( total_window_tp );
      
      // check for overlap
      uint64_t furthest = w.stop;
      
      // add first element
      intervals.insert( w );
      
      std::set<interval_t>::const_iterator ii = intervals0.begin();
      ++ii; // first one already added
      while ( ii != intervals0.end() )
	{
	  
	  interval_t w = *ii;	  
	  w.expand( total_window_tp );

	  if ( w.stop <= furthest ) { ++ii; continue; }
	  if ( w.start <= furthest ) 
	    { 
	      if ( w.stop > furthest+1 )
		intervals.insert( interval_t( furthest+1 , w.stop ) ); 
	    }
	  else 
	    intervals.insert( w );
	  
	  furthest = w.stop;
	  ++ii;
	}
    }
  else
    {

      std::set<interval_t>::const_iterator ii = intervals0.begin();
      while ( ii != intervals0.end() )
	{
	  interval_t w = *ii;
	  w.expand( total_window_tp );
	  if ( w.duration() > expected ) 
	    intervals.insert( w );
	  ++ii;
	}
    }
   
  std::vector<uint64_t> tp;
  std::vector<double> d;
  std::vector<int> segment_start;
  std::vector<int> segment_stop;

  // use slice_t to get all the data that we wish to keep
  // for now, only do this for 'extract' modes

  std::set<interval_t>::const_iterator ii = intervals.begin();
  while ( ii != intervals.end() )
    {
      
      slice_t slice( *this , signals(0) , *ii );
      
      const std::vector<double> * d1 = slice.pdata();      
      const std::vector<uint64_t> * tp1 = slice.ptimepoints();
      const int n = tp1->size();

      // did we managage to extract the entire segment?  if not, by
      // default, just drop
      
//       std::cout << "dets " << ( ii->start != (*tp1)[0]  ) << "\t" << ii->start << " " << (*tp1)[0] << "\n";
//       std::cout << "dets " << ( ii->stop  != (*tp1)[ tp1->size()-1 ]  ) << "\t" << ii->stop << " " << (*tp1)[ tp1->size()-1 ] << "\n";

      if ( n == 0 ) { ++ii; continue; }
//       if ( ii->start != (*tp1)[0] ) { ++ii; continue; }
//       if ( ii->stop  != (*tp1)[ tp1->size()-1 ] ) { ++ii; continue; }
      

      // okay, if here, we got the full complement...
      
      // track each individual segment
      segment_start.push_back( tp.size() );
      
      for (int i=0;i<n;i++)
	{	  
	  tp.push_back( (*tp1)[i] );
	  d.push_back( (*d1)[i] );
	}
      
      segment_stop.push_back( tp.size()-1 );
      
      ++ii;
    }


  //
  // Report on these segments
  //
  
  const int nsegs = segment_start.size();

  for (int s=0;s<nsegs;s++)
    {
      
      // Complete [ window(s) -- segment -- windows(s) ]  indices
      int j = segment_start[s];
      int k = segment_stop[s];
      int tot_length = k-j+1;
      
      // dump the raw signal?
      if ( raw_signal )
	for (int i=j;i<k;i++)
	  std::cout << "SLICE\t"
		    << id << "\t"
		    << "[" << globals::current_tag << "]\t"
		    << signals.label(0) << "\t"		  
		    << s << "\t"
		    << tp[i] << "\t" << d[i] << "\n";
      
      // if j and k define the start and stop in this space, then
      // figure out the actual window and the intervals (to the nearest bin)
      
      int index_start = j + window_samples * n_window;
      int index_length = (k-j+1) - ( 2 * n_window * window_samples);
      int index_stop  = index_start + index_length - 1;
      
      uint64_t tp_start = tp[ j ];
      uint64_t tp_stop  = tp[ k ];

           
      // std::cout << "slice: " << j << "/" << k << " (" << tot_length << " samples)\n";
      // std::cout << "tp   : " << tp_start << "/" << tp_stop << "\n";
      // std::cout << "index: " << index_start << "/" << index_stop << " (" << index_length << ")\n";
      // std::cout << "w-smp: " << window_samples << " (window length in samples)\n";

      
      // all     : j , tot_length
      // segment : index_start + index_length
      // window 'w' win_start[w] + window_samples

      std::vector<int> win_start;
      // before segment
      for (int w=0;w<n_window;w++)
	win_start.push_back( j + w*window_samples );
      // after segment
      for (int w=0;w<n_window;w++)
	win_start.push_back( index_stop + 1 + w*window_samples );

      //
      // FFT
      //
      
      //      std::cout << "index = " << index_length << "\n";
      FFT fftseg( index_length , MiscMath::nextpow2( index_length ) , Fs , FFT_FORWARD , WINDOW_HANN );
      
      fftseg.apply( &(d[index_start]) , index_length );

      double delta_lwr = globals::freq_band[ DELTA ].first;
      double delta_upr = globals::freq_band[ DELTA ].second;
      double delta = 0;

      double theta_lwr = globals::freq_band[ THETA ].first;
      double theta_upr = globals::freq_band[ THETA ].second;
      double theta = 0;

      double alpha_lwr = globals::freq_band[ ALPHA ].first;
      double alpha_upr = globals::freq_band[ ALPHA ].second;
      double alpha = 0;

      double sigma_lwr = globals::freq_band[ SIGMA ].first;
      double sigma_upr = globals::freq_band[ SIGMA ].second;
      double sigma = 0;

      double beta_lwr = globals::freq_band[ BETA ].first;
      double beta_upr = globals::freq_band[ BETA ].second;
      double beta = 0;

      double gamma_lwr = globals::freq_band[ GAMMA ].first;
      double gamma_upr = globals::freq_band[ GAMMA ].second;
      double gamma = 0;

      int N = fftseg.cutoff;
      for (int f=0;f<N;f++)
	{
	  if ( fftseg.frq[f] >= delta_lwr && fftseg.frq[f] < delta_upr ) delta += fftseg.X[f];
	  if ( fftseg.frq[f] >= theta_lwr && fftseg.frq[f] < theta_upr ) theta += fftseg.X[f];
	  if ( fftseg.frq[f] >= alpha_lwr && fftseg.frq[f] < alpha_upr ) alpha += fftseg.X[f];
	  if ( fftseg.frq[f] >= sigma_lwr && fftseg.frq[f] < sigma_upr ) sigma += fftseg.X[f];
	  if ( fftseg.frq[f] >= beta_lwr  && fftseg.frq[f] < beta_upr )  beta  += fftseg.X[f];
	  if ( fftseg.frq[f] >= gamma_lwr  && fftseg.frq[f] < gamma_upr ) gamma  += fftseg.X[f];
	}

      delta /= globals::band_width( DELTA );
      theta /= globals::band_width( THETA );
      alpha /= globals::band_width( ALPHA );
      sigma /= globals::band_width( SIGMA );
      beta  /= globals::band_width( BETA );
      gamma /= globals::band_width( GAMMA );

      std::cout << "SLICER" << "\t"
		<< id << "\t"
		<< "[" << globals::current_tag << "]\t"
		<< tp_start << ".." << tp_stop << "\t"
		<< 0 << "\t"
		<< delta << "\t"
		<< theta << "\t"
		<< alpha << "\t"
		<< sigma << "\t"
		<< beta << "\t"
		<< gamma << "\n";
      

      //
      // Whole interval
      //

      FFT fftall( tot_length , MiscMath::nextpow2( tot_length ) , Fs , FFT_FORWARD , WINDOW_HANN );
      fftall.apply( &(d[j]) , tot_length );
      
      // reset
      delta = theta = alpha = sigma = beta = gamma = 0;
      
      N = fftall.cutoff;
      
      for (int f=0;f<N;f++)
	{
	  if ( fftall.frq[f] >= delta_lwr && fftall.frq[f] < delta_upr ) delta += fftall.X[f];
	  if ( fftall.frq[f] >= theta_lwr && fftall.frq[f] < theta_upr ) theta += fftall.X[f];
	  if ( fftall.frq[f] >= alpha_lwr && fftall.frq[f] < alpha_upr ) alpha += fftall.X[f];
	  if ( fftall.frq[f] >= sigma_lwr && fftall.frq[f] < sigma_upr ) sigma += fftall.X[f];
	  if ( fftall.frq[f] >= beta_lwr  && fftall.frq[f] < beta_upr )  beta  += fftall.X[f];
	  if ( fftall.frq[f] >= gamma_lwr  && fftall.frq[f] < gamma_upr ) gamma  += fftall.X[f];
	}
      
      delta /= globals::band_width( DELTA );
      theta /= globals::band_width( THETA );
      alpha /= globals::band_width( ALPHA );
      sigma /= globals::band_width( SIGMA );
      beta  /= globals::band_width( BETA );
      gamma /= globals::band_width( GAMMA );
      
      std::cout << "SLICER" << "\t"
		<< id << "\t"
		<< "[" << globals::current_tag << "]\t"
		<< tp_start << ".." << tp_stop << "\t"
		<< "A" << "\t"
		<< delta << "\t"
		<< theta << "\t"
		<< alpha << "\t"
		<< sigma << "\t"
		<< beta << "\t"
		<< gamma << "\n";	   


      const bool spindle_bins = true;

      if ( spindle_bins )
	{
	  const double lwr = 8;
	  const double upr = 16;
	  const double thr = 2;
	  
	  const int    nbins = ( upr - lwr ) / thr;
	  //	  std::vector<double> s( bins );
	  	  
	  for (int b=0;b<nbins;b++)
	    {
	      const double l0 = lwr + b * thr;
	      const double u0 = l0 + thr;
	      
	      double sx = 0;
	      int nx = 0;

	      for (int f=0;f<N;f++)
		{
		  if ( fftall.frq[f] >= l0 && fftall.frq[f] < u0 ) 
		    {
		      sx += fftall.X[f];
		      nx++;
		    }
		}

	      // get mean
	      sx /= (double)nx;
	      
	      std::cout << "SLICE\tSPINDLE\t" 
			<< id << "\t"
			<< "[" << globals::current_tag << "]\t"
			<< tp_start << ".." << tp_stop << "\t"	      
			<< s << "\t" 
			<< l0 << "\t" 
			<< u0 << "\t"
			<< sx << "\n";
	    }
	}
      
      
    
      std::vector<std::string> win_label;
      for (int w=0;w<n_window;w++) win_label.push_back(  Helper::dbl2str( - n_window + w ) );
      for (int w=0;w<n_window;w++) win_label.push_back(  Helper::int2str( w+1 ) );
      
      //
      // Windows
      //
    
      for (int w=0;w<win_start.size();w++)
	{
	  delta = theta = alpha = sigma = beta = gamma = 0;
	  
	  FFT fftwin( window_samples , MiscMath::nextpow2( window_samples ) , Fs , FFT_FORWARD , WINDOW_HANN );      
	  fftwin.apply( &(d[win_start[w]]) , window_samples );
	  
	  
	  int N = fftwin.cutoff;
	  
	  for (int f=0;f<N;f++)
	    {
	      if ( fftwin.frq[f] >= delta_lwr && fftwin.frq[f] < delta_upr ) delta += fftwin.X[f];
	      if ( fftwin.frq[f] >= theta_lwr && fftwin.frq[f] < theta_upr ) theta += fftwin.X[f];
	      if ( fftwin.frq[f] >= alpha_lwr && fftwin.frq[f] < alpha_upr ) alpha += fftwin.X[f];
	      if ( fftwin.frq[f] >= sigma_lwr && fftwin.frq[f] < sigma_upr ) sigma += fftwin.X[f];
	      if ( fftwin.frq[f] >= beta_lwr  && fftwin.frq[f] < beta_upr )  beta  += fftwin.X[f];
	      if ( fftwin.frq[f] >= gamma_lwr  && fftwin.frq[f] < gamma_upr ) gamma  += fftwin.X[f];
	    }
	   
	   delta /= globals::band_width( DELTA );
	   theta /= globals::band_width( THETA );
	   alpha /= globals::band_width( ALPHA );
	   sigma /= globals::band_width( SIGMA );
	   beta  /= globals::band_width( BETA );
	   gamma /= globals::band_width( GAMMA );
	   
	   std::cout << "SLICER" << "\t"
		     << id << "\t"
		     << "[" << globals::current_tag << "]\t"
		     << tp_start << ".." << tp_stop << "\t"
		     << win_label[w] << "\t"		     
		     << delta << "\t"
		     << theta << "\t"
		     << alpha << "\t"
		     << sigma << "\t"
		     << beta << "\t"
		     << gamma << "\n";	   
	   
	   


	}
      
      

      // 	std::cout << "w " << win_start[w] << " " << win_start[w] + window_samples - 1 << "\n";
      
    }
  
  
}


