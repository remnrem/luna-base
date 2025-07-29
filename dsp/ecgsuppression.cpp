
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


#include "ecgsuppression.h"

#include "resample.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "fir.h"
#include "dsp/spline.h"
#include "spectral/welch.h"
#include "fftw/fftwrap.h"

#include "eval.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

#include <iomanip>
#include <cmath>
#include <map>

extern writer_t writer;

extern logger_t logger;

// assumes broad range of 'normal' sleeping heart beat is 40-100 

std::vector<uint64_t> rpeaks_t::beats( interval_t & interval ) const
{  
  std::vector<uint64_t> ret;
  for (int i=0;i<R_t.size();i++)
    {
      if ( R_t[i] > interval.stop ) break;
      if ( R_t[i] >= interval.start ) ret.push_back( R_t[i] );
    }
  return ret;
}

double rpeaks_t::bpm( interval_t & interval , double lwr , double upr ) 
{
  
  std::set<int> trk;
  double b = 0;
  double s = interval.duration() * globals::tp_duration ;
  for (int i=0;i<R_t.size();i++)
    {
      if ( R_t[i] > interval.stop ) break;
      if ( R_t[i] >= interval.start ) { trk.insert(i); ++b; }      
    }  


  double ret = (b/s)*60.0;
  
  // zero-out bad epochs
  
  if ( lwr != 0 )
    {
      if ( ret < lwr || ret > upr )
	{
	  std::vector<uint64_t> c;
	  std::vector<uint64_t> ci;
	  for (int i=0;i<R_t.size();i++)
	    if ( trk.find( i ) == trk.end() ) 
	      {
		c.push_back( R_t[i] );
		ci.push_back( R_i[i] );
	      }
	  R_t = c;
	  R_i = ci;
	  
	}
    }

  // return BPM
  return ret;

}

int rpeaks_t::clean( double mn , double mx ) 
{
  
  uint64_t tpmin = mn * globals::tp_1sec;
  uint64_t tpmax = mx * globals::tp_1sec;

  if ( R_t.size() < 2 ) return 0;

  std::set<int> trk;

  const int n = R_t.size()-1;

  for (int i=1;i<n;i++)
    {
      double t12 = ( R_t[i] - R_t[i-1] ) * globals::tp_duration;
      double t23 = ( R_t[i+1] - R_t[i] ) * globals::tp_duration;
      
      if ( t12 < mn )
	{
	  trk.insert(i-1);
	  trk.insert(i);	  
	}
      
      if ( t23 < mn )
	{
	  trk.insert(i);
	  trk.insert(i+1);	  
	}
      
      // isolated beat ; just remove middle one
      if ( t12 > mx && t23 > mx )
	{
	  trk.insert(i);	  
	}
      
    }
  
  if ( trk.size() )
    {
      std::vector<uint64_t> c;
      std::vector<uint64_t> ci;
      for (int i=0; i<R_t.size(); i++)
	if ( trk.find(i) == trk.end() ) 
	  {
	    c.push_back( R_t[i] );
	    ci.push_back( R_i[i] );
	  }
      R_t = c;
      R_i = ci;
    }
  
  // report number of removed beats
  //  std::cout << "DET removed beats = " << trk.size() << "\n";
  return trk.size();
}

	 

void dsptools::ecgsuppression( edf_t & edf , param_t & param )
{
 
  // do not update signal (but remove epochs w/out ECG)
  // i.e. this is just for methods comparisons, in practice
  
  bool nosuppression = param.has( "no-suppress" ) ;
  
  // EEG channels (to be modified)
  std::string signal_label = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( signal_label );    
  const int ns = signals.size();
  
  // if not specified, SR is set to the first channel
  int sr = param.has( "sr" ) ? param.requires_int( "sr" ) : 0 ;
  if ( sr == 0 ) sr = edf.header.sampling_freq( signals(0) ) ;
  
  logger << " setting SR to " << sr << "\n";
  

  // default is to leave EEG in ECG-bad epochs 'as is'
  const bool mask_bad_epochs = param.has( "mask-bad-epochs" );
  
  //
  // ECG channel
  //

  std::string ecg_label = param.requires("ecg");

  int ecg_n = edf.header.signal( ecg_label );

  if ( ecg_n == -1 ) 
    {
      logger << "could not find ECG (label " << ecg_label << "), skipping ECG suppression\n";
      return;
    }

  //
  // check SR for each channel
  //
  
  for (int s=0;s<ns;s++)
    {
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	resample_channel( edf, signals(s) , sr );

    }

  // and for ECG...
  if ( edf.header.sampling_freq( ecg_n ) != sr ) 
    resample_channel( edf, ecg_n , sr );
  
  
  //
  // pull entire trace (assumes continuous/contiguous data)
  //
  
  interval_t interval = edf.timeline.wholetrace();
  
  
  //
  // ECG
  //
  
  slice_t slice1( edf , ecg_n , interval );    
  const std::vector<double> * ecg = slice1.pdata();
  const std::vector<uint64_t> * tp = slice1.ptimepoints();  

  //  logger << "  pulled " << tp->size() << " samples\n";

  //
  // find ECG peaks
  //
  
  rpeaks_t peaks = mpeakdetect( edf , ecg , tp , sr );

  logger << "  detected R peaks\n";

  //
  // clean beats
  //
  
  int removed = peaks.clean( 0.3 , 2 );  

  //  logger << "  removed " << removed << " bad beats\n";
  
  
  //
  // find bad ECG epochs (i.e. calculate implied HR)
  //
  
  int ne = edf.timeline.set_epoch( 30 , 30 );
  
  std::vector<double> epoch_bpm;
  std::map<int,double> epoch2bpm;
  
  int removed_epochs = 0;
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      
      
      if ( epoch == -1 ) break;
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      double bpm = peaks.bpm( interval , 40 , 100 );

      epoch2bpm[ epoch ] = bpm;

      if ( bpm > 40 && bpm < 100 ) 
	epoch_bpm.push_back( bpm );
      else 
	{
	  
	  // mask that epoch?
	  if ( mask_bad_epochs ) 
	    edf.timeline.set_epoch_mask( epoch );
	  else
	    {
	      // else, just remove all R peaks from that epoch, i.e. so no correction is made
	      
	    }
	  ++removed_epochs;
	}
      
      
            
    }


  //
  // Also remove outliers (note, above, we exclude non-biological values <40, >100 from calculations)
  //
  
  double mean_bpm = MiscMath::mean( epoch_bpm );
  double sd_bpm = MiscMath::sdev( epoch_bpm , mean_bpm );

  double lwr95 = mean_bpm - 2 * sd_bpm;
  double upr95 = mean_bpm + 2 * sd_bpm;
  
  writer.value( "BPM" , mean_bpm );
  writer.value( "BPM_L95" , lwr95 );
  writer.value( "BPM_U95" , upr95 );
  
  //
  // Output for each epoch
  //

  for (int epoch=0; epoch<ne; epoch++)
    {
      
      if ( epoch2bpm[ epoch ] < lwr95 || epoch2bpm[ epoch ] > upr95 ) 
	{
	  if ( ! edf.timeline.masked( epoch ) ) ++removed_epochs;
	  edf.timeline.set_epoch_mask( epoch );
	}
      
      writer.epoch( edf.timeline.display_epoch( epoch ) );
      writer.value( "BPM" , epoch2bpm[ epoch ] );
      writer.value( "BPM_MASK" , edf.timeline.masked( epoch ) ? 1 : 0  );
      
    }
  writer.unepoch();
  writer.value( "BPM_N_REMOVED" , removed_epochs );
  writer.value( "BPM_PCT_REMOVED" , removed_epochs/(double)ne );


  //
  // process each signal
  //

  for (int s = 0 ; s < ns ; s++ )
    {

      // skip annotations channels
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      // skip ECG if that is included
      if ( signals.label(s) == ecg_label ) continue;

      //
      // Output stratified by signal
      //

      writer.level( signals.label(s) , globals::signal_strat );
      
      // get average waverform/adjust EEG
      slice_t slice( edf , signals(s) , interval );    
      std::vector<double> * sig  = slice.nonconst_pdata();
      
      if ( ecg->size() != sig->size() ) 
	Helper::halt( "internal error, signals different length in coherence()");


      //
      // Mean centre signal first
      //
      
      MiscMath::centre( sig );
      
      
      //
      // Build up average profile around peak
      //
      
      const int npeaks = peaks.R_t.size();
      
      // i.e. averaged signal is only up to 2.0 seconds max
      const int mxdur = 2 * sr; 
      const int pnts_200ms = 0.2 * sr;
      
      // matrix is R/C time-points/beats
      // i.e. so we will average for each time-point, over all beats
      
      std::vector<std::vector<double> > mat( mxdur );
      for (int i=0;i<npeaks-1;i++)
	{
	  uint64_t p = peaks.R_i[i];
	  if ( p < pnts_200ms ) continue;	  
	  uint64_t p1 = p - pnts_200ms;
	  uint64_t p2 = peaks.R_i[i+1] - pnts_200ms;
	  
	  // ensure no larger than 2 seconds
	  if ( p2 - p1 >= mxdur ) p2 = p1 + mxdur - 1;	  
	  int c = 0;
	  for (int j=p1;j<=p2;j++) 
	    mat[c++].push_back( (*sig)[j] );	    
	  
	}
      
           
      
      //
      // Get average artifact matrix (i.e. in EEG signals) 
      //
      
      std::vector<double> art( mat.size() , 0 );
      
      double art_rms = 0;

      for (int i=0;i<mat.size();i++)
	{
	  
	  double denom = mat[i].size() / (double)npeaks;
	  
	  double m = MiscMath::mean( mat[i] );
	  
	  art[i] = m * denom;
	  
	  //
	  // Output correction matrix
	  //

	  writer.level( i , "SP" );
	  writer.value( "ART" , art[i] );

	  art_rms += art[i] * art[i];
	  
	}

      writer.unlevel( "SP" );

      // RMS of mean artifact signature, i.e. a simple index
      // of the extent of cardiac contamination

      art_rms /= (double)mat.size();
      art_rms = sqrt( art_rms );
      writer.value( "ART_RMS" , art_rms );

      
    
      // 
      // Subtract out around peaks (will leave non-peaks untouced, but 
      // really aberrant epochs should have been masked out anyway)
      //

      std::vector<double> nsig = *sig;

      for (int i=0;i<npeaks-1;i++)
	{
	  uint64_t p = peaks.R_i[i];
	  if ( p < pnts_200ms ) continue;	  
	  uint64_t p1 = p - pnts_200ms;
	  uint64_t p2 = peaks.R_i[i+1] - pnts_200ms;

	  // ensure no larger than 2 seconds
	  if ( p2 - p1 >= mxdur ) p2 = p1 + mxdur - 1;	  
	  int c = 0;
	  for (int j=p1;j<=p2;j++) 
	    nsig[j] -= art[c++];
	}
      
      //
      // Put EEG signal back
      //

      if ( ! nosuppression )
	{
	  logger << " updating ECG-corrected signal " << signals.label(s) << "\n";
	  edf.update_signal( signals(s) , &nsig );      
	}


      //
      // Next signal
      //

    }
  
  writer.unlevel( globals::signal_strat );
  
}

struct rpeak_opt_t {

  rpeak_opt_t()
  {
    ripple = 0.01;
    tw = 3;
    flwr = 5;
    fupr = 20;
    median_filter_window = 9;

    //smoothed Z
    lag_sec = 10;      
    influence = 0.001;
 
    th = 2;
    th2 = 1;
    max = 0; // no max.

    mindur_sec = 0.04; // 40 msec
    mindur2_sec = 0.04; // same
  }

  double ripple;
  double tw;
  double flwr;
  double fupr;
  int median_filter_window;

  // smoothedZ peak finding
  double lag_sec;
  double influence;

  // core
  double th;
  double max; 
  double mindur_sec; // 40 ms
  
  // flanking region
  double th2;
  double mindur2_sec; // 40 ms

      
  
};


rpeaks_t dsptools::mpeakdetect2( const std::vector<double> * d ,
				 const std::vector<uint64_t> * tp ,
				 int Fs ,
				 const rpeak_opt_t & opt )
{
  
  const int n = d->size();
  
  if ( tp->size() != n )
    Helper::halt("error in mpeakdetect2");
  
  rpeaks_t peaks;
  
  // bandpass 5 - 15 Hz
  std::vector<double> ripple( 1, opt.ripple ); // 40dB
  std::vector<double> tw( 1, opt.tw ); // transition width...can flex  
  std::vector<double> bpf = dsptools::apply_fir( *d , Fs , fir_t::BAND_PASS ,
						 1, // Kaiser window
						 ripple , tw , // ripple, TW
						 opt.flwr , opt.fupr ) ;  // f1, f2
  
  // t0    t1     t2     ...
  // X0    X1     X2     ...
  // .     d1-d0  d2-d1  ...

  // differentiate and sqr data
  std::vector<double> sq( n );
  sq[0] = 0;
  for (int i=1;i<n;i++) 
    sq[i] = ( bpf[i] - bpf[i-1] ) * ( bpf[i] - bpf[i-1] );

  // moving window integration:
  //  default to 7 points, higher for higher SRs
  const int ds = Fs > 256 ? 7 * Fs / (double)256  : 7 ; 
  std::vector<double> ss( n , 0 );
  for (int i=0; i<n; i++)
    for (int j=0; j < ( ds > i ? i : ds ) ; j++ )
      ss[i] += sq[i-j];
  
  // median filter size 9
  ss = MiscMath::median_filter( ss , opt.median_filter_window );
  
  std::vector<interval_t> regions;

  const bool only_positives = true;
  
  std::vector<int> pks = MiscMath::smoothedZ( ss ,
					      Fs * opt.lag_sec,
					      opt.th,
					      opt.influence,
					      Fs * opt.mindur_sec,
					      opt.max,
					      opt.th2,
					      Fs * opt.mindur2_sec,
					      only_positives,
					      &regions);  

  // within each region, look to bpf signal to get min/max  
  std::vector<int> maxloc, minloc;
  std::vector<double> maxval, minval;
  
  const int npks = regions.size();
  
  // get max point within each segment
  for (int i=0;i<regions.size();i++)
    {

      int p1 = regions[i].start;
      int p2 = regions[i].stop;
      
      double mx = bpf[p1];
      int mxi = p1;

      double mn = bpf[p1];
      int mni = p1;
      
      for (int j=p1+1;j<=p2;j++) 
	{
	  
	  if ( bpf[j] > mx ) 
	    {
	      mx = bpf[j];
	      mxi = j;
	    }
	  
	  if ( bpf[j] < mn ) 
	    {
	      mn = bpf[j];
	      mni = j;
	    }
	}
      
      maxloc.push_back( mxi );
      maxval.push_back( mx );
      
      minloc.push_back( mni );
      minval.push_back( mn );
    }

  // check for lead inversion
  //  --> do minima precede maxima more often than not?
    
  int inv1 = 0;
  for (int i=0;i<npks;i++)
    if ( minloc[ i ] < maxloc[ i ] ) ++inv1;
  const double p_inverted = inv1 / (double)minloc.size();  
  const bool inverted = p_inverted > 0.5;

  peaks.R_t.resize( npks );
  peaks.R_i.resize( npks );
  
  for (int i=0;i<npks;i++)
    {
      peaks.R_t[ i ] = (*tp)[ inverted ? minloc[i] : maxloc[i] ] ;
      peaks.R_i[ i ] = inverted ? minloc[i] : maxloc[i] ;
    }
  

  //
  // Output
  //

  peaks.npks = npks;
  peaks.p_inverted = p_inverted;
  peaks.inverted = inverted;
    
  return peaks;

  
}



rpeaks_t dsptools::mpeakdetect( const edf_t & edf , 
				const std::vector<double> * d , 
				const std::vector<uint64_t> * tp , 				
				int Fs , 
				const std::vector<double> * eeg , 
				bool * force )
{
  
  const int n = d->size();

  if ( tp->size() != n ) Helper::halt("error in mpeakdetect");

  rpeaks_t peaks;
  
  // mean-center ECG
  std::vector<double> x = MiscMath::centre( *d );
  
  // band-pass filter: 0.5-40Hz 

  logger << "  filtering ECG...\n";
  std::vector<double> ripple( 1, 0.01 ); // 40dB
  std::vector<double> tw( 1, 2 ); // transition width...can flex  
  std::vector<double> bpf = dsptools::apply_fir( x , Fs , fir_t::BAND_PASS ,
						 1, // Kaiser window
						 ripple , tw , // ripple, TW
						 0.5 , 40 ) ;  // f1, f2

  
  // differentiate and sqr data
  std::vector<double> sq( n - 1 );
  for (int i=0;i<n-1;i++) 
    sq[i] = ( bpf[i+1] - bpf[i] ) * ( bpf[i+1] - bpf[i] );
  
  // integrate over 7 points (i.e. sum) - but expanded for higher SRs  
  const int ds = Fs > 256 ? 7 * Fs / (double)256  : 7 ; 
  
  std::vector<double> ss( n - 1 , 0 );
  for (int i=0; i < n-1; i++)
    for (int j = 0 ; j < ds ; j++ )
      if ( i-j >= 0 ) ss[i] += sq[i-j];
  
  // median filter, window size of 10
  std::vector<double> mf = MiscMath::median_filter( ss , 10 );
  
  // remove filter delay
  int delay = ceil( ds / 2.0 );
  
  // i.e. skip first 'delay-1' elements
  std::vector<double> mdfint;
  for (int i=delay-1;i<ss.size();i++)
    mdfint.push_back( ss[i] );
  
  //  std::cout << " delay = " << delay << " " << mdfint.size() << " " << tp->size() << "\n";

  
  //
  // segment search 'sq'
  //
  
  int len = mdfint.size();
  
  // avoid issues, calculate max within each 30s EPOCH, taking median max value
  // over epochs
  int e30 = Fs * 30;

  // note, as we dropped some samples, need to ensure we count all
  // thus e <= ne below
  int ne = len / e30;
  
  std::vector<double> maxvals;
  for (int e = 0; e <= ne ; e++)
    {
      int s1 = e * e30;
      int s2 = s1 + e30 - 1 ;
      double max_h = 0;
      for (int i = s1; i <= s2 ; i++ ) { 
	if ( i < len ) 
	  if ( mdfint[i] > max_h ) max_h = mdfint[i];
      }
      maxvals.push_back( max_h );
    }

  double max_h = median_destroy( &maxvals[0] , maxvals.size() );
  double thresh = 0.2;  
  double th = max_h * thresh;

  // get segments
  std::vector<int> left;
  std::vector<int> right;
  bool inregion = false;
  for (int i=0;i<mdfint.size();i++)
    {
      const bool pk = mdfint[i] > th ;
      if ( pk && ! inregion ) left.push_back(i);
      if ( inregion && !pk ) right.push_back(i);
      inregion = pk;
    }
  if ( inregion ) right.push_back( mdfint.size()-1 );

  std::vector<int> maxloc, minloc;
  std::vector<double> maxval, minval;

  // get max point within each segment
  for (int i=0;i<left.size();i++)
    {

      int p1 = left[i];
      int p2 = right[i];
      
      // sanity check on size of region?
      double mx = bpf[p1];
      int mxi = p1;

      double mn = bpf[p1];
      int mni = p1;
      
      for (int j=p1;j<=p2;j++) 
	{

	  if ( bpf[j] > mx ) 
	    {
	      mx = bpf[j];
	      mxi = j;
	    }
	  
	  if ( bpf[j] < mn ) 
	    {
	      mn = bpf[j];
	      mni = j;
	    }
	}
      
      maxloc.push_back( mxi );
      maxval.push_back( mx );

      minloc.push_back( mni );
      minval.push_back( mn );
    }

  // check for lead inversion
  //  --> do minima precede maxima more often than not?
    
  int inv1 = 0;
  for (int i=0;i<minloc.size();i++)
    if ( minloc[ i ] < maxloc[ i ] ) ++inv1;
  const double p_inverted = inv1 / (double)minloc.size();  
  const bool inverted = p_inverted ; 
  
  // Or, swap in user-forced value (option not for general use...)
  //  if ( force != NULL ) inverted = *force;
  
  std::vector<uint64_t> R_t;  
  std::vector<uint64_t> R_i;    

  std::vector<uint64_t> t_R; // mirror... 
  std::vector<uint64_t> i_R; 
  
  for (int i=0;i<maxloc.size();i++)
    {
      uint64_t ht = (*tp)[ inverted ? minloc[i] : maxloc[i] ] ;
      R_t.push_back( ht );
      R_i.push_back( inverted ? minloc[i] : maxloc[i] );

      uint64_t th = (*tp)[ inverted ? maxloc[i] : minloc[i] ] ;
      t_R.push_back( th );
      i_R.push_back( inverted ? maxloc[i] : minloc[i] );
      
    }
  

  //
  // take the one with the smallest HRV  (threshold at 2)
  //
  
  // const int nb = R_i.size();
  
  // const uint64_t two_sec = globals::tp_1sec + globals::tp_1sec; 
  // const uint64_t point3_sec = 0.3 * globals::tp_1sec;
  
  // std::vector<double> v1, v2;
  // for (int i=1;i<nb;i++)
  //   {
  //     double i1 = R_t[i] - R_t[i-1];
  //     if ( i1 > point3_sec && i1 < two_sec ) v1.push_back( i1 * globals::tp_duration );

  //     double i2 = t_R[i] - t_R[i-1];
  //     if ( i2 > point3_sec && i2 < two_sec ) v2.push_back( i2 * globals::tp_duration ); 
  //   }
  // double var1 = MiscMath::variance( v1 );
  // double var2 = MiscMath::variance( v2 );

  //
  // Output
  //

  writer.value( "P_INV" , p_inverted );
  writer.value( "ECG_INVERTED" , inverted );

  
  // writer.value( "ECG_VAR1" , var1 );
  // writer.value( "ECG_VAR2" , var2 );
  // writer.value( "ECG_VAR_FLIPPED" , var1 < var2 ? 0 : 1 );
  
  //
  // Based on variance if inter-beat intervals, select the orientation that
  // minimizes this value
  //
  
  if ( p_inverted )
    {
      peaks.R_t = t_R;
      peaks.R_i = i_R;
    }
  else
    {
      peaks.R_t = R_t;
      peaks.R_i = R_i;
    }

  // total time
  double mx = (*tp)[ tp->size() - 1 ];
  
  return peaks;
  
}


void dsptools::bpm( edf_t & edf , param_t & param )
{

  std::string ecg_label = param.requires("ecg");
  int ecg_n = edf.header.signal( ecg_label );
  if ( ecg_n == -1 ) 
    {
      logger << "could not find ECG (label " << ecg_label << "), skipping ECG suppression\n";     
      return;
    }


  //
  // Sample rate
  //

  int sr = edf.header.sampling_freq( ecg_n );


  //
  // Pull entire trace 
  //
  
  interval_t interval = edf.timeline.wholetrace();
  
  //
  // ECG
  //
  
  slice_t slice1( edf , ecg_n , interval );    
  const std::vector<double> * ecg = slice1.pdata();
  const std::vector<uint64_t> * tp = slice1.ptimepoints();  

  //
  // find ECG peaks
  //

  rpeaks_t peaks = mpeakdetect( edf , ecg , tp , sr );

  
  //
  // clean beats
  //
  
  int removed = peaks.clean( 0.3 , 2 );  
    
  //
  // find bad ECG epochs w.r.t inferred HR
  //
  
  int ne = edf.timeline.set_epoch( 30 , 30 );
  
  std::vector<double> epoch_bpm;
  std::map<int,double> epoch2bpm;
  
  int removed_epochs = 0;
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      
      
      if ( epoch == -1 ) break;
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      double bpm = peaks.bpm( interval , 40 , 100 );
      
      epoch2bpm[ epoch ] = bpm;

      if ( bpm > 40 && bpm < 100 ) 
	epoch_bpm.push_back( bpm );
      else 
	{
	  edf.timeline.set_epoch_mask( epoch );
	  ++removed_epochs;
	}
            
    }


  //
  // Remove statistical outliers
  //
  
  double mean_bpm = MiscMath::mean( epoch_bpm );
  double sd_bpm = MiscMath::sdev( epoch_bpm , mean_bpm );

  writer.value( "BPM" , mean_bpm );

  double lwr95 = mean_bpm - 2 * sd_bpm;
  double upr95 = mean_bpm + 2 * sd_bpm;

  writer.value( "BPM_L95" , lwr95 );
  writer.value( "BPM_U95" , upr95 );

  double bpm2 = 0 ; 
  int cnt2 = 0;
  
  for (int epoch=0; epoch<ne; epoch++)
    {
      
      if ( epoch2bpm[ epoch ] < lwr95 || epoch2bpm[ epoch ] > upr95 ) 
	{
	  if ( ! edf.timeline.masked( epoch ) ) ++removed_epochs;
	  edf.timeline.set_epoch_mask( epoch );
	}
      
      bpm2 += epoch2bpm[ epoch ] ;
      ++cnt2;
      
      writer.epoch( edf.timeline.display_epoch( epoch ) );
      writer.value( "BPM_MASK" , edf.timeline.masked( epoch ) ? 1 : 0 );
      writer.value( "BPM" , epoch2bpm[ epoch ] );

    }
  writer.unepoch();

  if ( cnt2 != 0 ) 
    writer.value( "BPM2" , bpm2 / (double)cnt2 );
  
}





//
// HRV metrics (time-domain and frequency domain)
//


struct hrv_opt_t {
  bool freq_domain;
  bool time_domain;
  double rr_lwr;
  double rr_upr;
  int median_filter_width;
  int welch_nsamples;

  // for annotations
  edf_t * edf;
  bool annot_stratify;
  bool inst_stratify;
  std::vector<annot_t*> annots;  
};


void dsptools::hrv( edf_t & edf , param_t & param )
{
    
  //
  // Options
  //

  
  const bool annotate = param.has( "add-annot" ) || param.has( "add-annot-ch" );
  const bool annotate_ch = param.has( "add-annot-ch" );

  
  std::string alabel = "Rpk";
  if ( annotate )
    {
      if ( annotate_ch )
	alabel = param.empty( "add-annot-ch" ) ? "Rpk" :param.value( "add-annot-ch" ) ;
      else
	alabel = param.empty( "add-annot" ) ? "Rpk" :param.value( "add-annot" ) ;
    }

  
  //
  // R-peak detection options
  //

  rpeak_opt_t ropt;
  
  if ( param.has( "rp-lag" ) ) ropt.lag_sec = param.requires_dbl( "rp-lag" );
  if ( param.has( "rp-infl" ) ) ropt.influence = param.requires_dbl( "rp-infl" );
  if ( param.has( "rp-th" ) ) ropt.th = param.requires_dbl( "rp-th" );
  if ( param.has( "rp-th2" ) ) ropt.th2 = param.requires_dbl( "rp-th2" );
  if ( param.has( "rp-max" ) ) ropt.max = param.requires_dbl( "rp-max" );
  if ( param.has( "rp-dur" ) ) ropt.mindur_sec = param.requires_dbl( "rp-dur" );
  if ( param.has( "rp-dur2" ) ) ropt.mindur2_sec = param.requires_dbl( "rp-dur2" );

  if ( param.has( "rp-ripple" ) ) ropt.ripple = param.requires_dbl( "rp-ripple" );
  if ( param.has( "rp-tw" ) ) ropt.tw = param.requires_dbl( "rp-tw" );
  if ( param.has( "rp-flwr" ) ) ropt.flwr = param.requires_dbl( "rp-flwr" );
  if ( param.has( "rp-fupr" ) ) ropt.fupr = param.requires_dbl( "rp-fupr" );
  if ( param.has( "rp-w" ) ) ropt.median_filter_window = param.requires_dbl( "rp-w" );
        
  
  //
  // HRV analysis options
  //

  hrv_opt_t opt;

  opt.edf = &edf;

  opt.freq_domain = param.has( "freq-domain" ) ? param.yesno( "freq-domain" ) : true ;
  opt.time_domain = param.has( "time-domain" ) ? param.yesno( "time-domain" ) : true ; 
  opt.rr_lwr = param.has( "lwr" ) ? param.requires_dbl( "lwr" ) : 0.3;
  opt.rr_upr = param.has( "upr" ) ? param.requires_dbl( "upr" ) : 2;
  opt.median_filter_width = param.has( "w" ) ? param.requires_int( "w" ) : 5 ;
  opt.welch_nsamples = param.has( "ns" ) ? param.requires_int( "ns" ) : 512 ; 

  
  //
  // Get time-domain HRV values stratified by annot
  //

  opt.annot_stratify = param.has( "annot" );

  opt.inst_stratify = param.has( "by-instance" );
  
  const std::vector<std::string> annot_labels = param.strvector( "annot" );
  for (int a=0; a<annot_labels.size(); a++)
    {
      annot_t * annot = edf.annotations->find( annot_labels[a] );
      if ( annot != NULL ) opt.annots.push_back( annot );
    }
    

  //
  // Pre-processing for annot-stratified analyses
  //  --> a list of all annots/instance pairs (. if not instance)
  //

  std::map<std::pair<std::string,std::string> , std::set<interval_t> > atypes;
  
  if ( opt.annot_stratify )
    {
      // iterate over each annotation
      for (int a=0; a<opt.annots.size(); a++)
	{
	  annot_t * annot = opt.annots[a] ;
	  const std::string label = annot->name;	  
	  // iterate over all events
	  annot_map_t::const_iterator ii = annot->interval_events.begin();
	  while ( ii != annot->interval_events.end() )
	    {
	      const instance_idx_t & idx = ii->first;
	      const std::pair<std::string,std::string> key(  annot->name, opt.inst_stratify ? idx.id : "." );
	      atypes[ key ].insert( idx.interval );		  
	      ++ii;
	    }	  
	}

      // enumerate
      logger << "  stratifying analyses by " << atypes.size() << " annotation types\n";
      std::map<std::pair<std::string,std::string> , std::set<interval_t> >::const_iterator aa = atypes.begin();
      while ( aa != atypes.end() )
	{
	  logger << "   --> " << aa->first.first <<
	    ( opt.inst_stratify ? "/" + aa->first.second : "" ) << " "
		 << " n = " << aa->second.size() << " intervals\n";
	  ++aa;
	}
    }


  //
  // Do by epochs or no?
  //

  const bool by_epoch = edf.timeline.epoched();
  
  const bool epoch_output = by_epoch && param.has( "epoch" );

  if ( param.has( "epoch" ) && ! by_epoch )
    Helper::halt( "data are not epoched yet" );

  
  //
  // ECG signals
  //
  
  const std::string signal_label = param.requires( "sig" );

  const bool NO_ANNOTS = true; 

  signal_list_t signals = edf.header.signal_list( signal_label , NO_ANNOTS );

  const int ns = signals.size();

  if ( ns == 0 ) return ;

  
  //
  // Get each ECG signal
  //

  for (int s=0; s<ns; s++)
    {
      
      writer.level( signals.label(s) , globals::signal_strat );
      
      const int sr = edf.header.sampling_freq( signals(s) );

      //
      // Add annotations?
      //

      annot_t * r_annot = NULL;

      if ( annotate_ch )
	r_annot = edf.annotations->add( alabel + "_"  + signals.label(s) );
      else if ( annotate ) 
	r_annot = edf.annotations->add( alabel );

      //
      // Optionally, track all all rpeaks for second-round annot-stratified
      // analyses
      //

      std::set<uint64_t> tps;
      
      //
      // Iterate over epochs
      //
      
      const int ne = by_epoch ? edf.timeline.first_epoch() : 1 ; 

      std::vector<hrv_res_t> results;
	
      while ( 1 )
	{
	  
	  int epoch = by_epoch ? edf.timeline.next_epoch() : 1 ; 

	  if ( epoch == -1 ) break;

	  interval_t interval = by_epoch ? edf.timeline.epoch( epoch ) : edf.timeline.wholetrace(); 

	  if ( by_epoch ) 
	    writer.epoch( edf.timeline.display_epoch( epoch ) );
	  
	  //
	  // Get ECG signal and time-points
	  //

	  slice_t slice( edf , signals(s) , interval );
	  
	  const std::vector<double> * ecg = slice.pdata();

	  const std::vector<uint64_t> * tp = slice.ptimepoints();
      
	  //
	  // Detect R-peaks
	  //
	  
	  rpeaks_t peaks = mpeakdetect2( ecg , tp , sr , ropt );

	  //
	  // Track?
	  //

	  if ( opt.annot_stratify )
	    for (int i=0; i<peaks.R_t.size(); i++)
	      tps.insert( peaks.R_t[i] );
			
	  //
	  // Add annotations
	  //
	  
	  if ( annotate )
	    {
	      const int npeaks = peaks.R_t.size();
	      
	      for (int i=0; i<npeaks; i++)
		r_annot->add( "." , interval_t( peaks.R_t[i] , peaks.R_t[i] ) , signals.label(s) );
	    }
      
	  //
	  // Derive and resample RR intervals
	  //
	  
	  rr_intervals_t rr( peaks , opt );

	  //
	  // Store/report epochwise results
	  //
	  
	  if ( by_epoch )
	    {
	      results.push_back( rr.res );
	      
	      if ( epoch_output ) 
		rr.res.write( opt );
	    }
	  else
	    {
	      // write whole-trace results
	      rr.res.write( opt );		
	    }
	  
	  
	  //
	  // Next epoch
	  //
	  
	  if ( ! by_epoch ) break;
	  
	}


      //
      // Summarize over epochs
      //
      
      if ( by_epoch )
	{
	  writer.unepoch();

	  // get means over epochs
	  hrv_res_t means = hrv_res_t::summarize( results );

	  // output
	  means.write( opt );
	  
	}


      //
      // Second round : annotation-stratified analyses
      //

      if ( opt.annot_stratify )
	{

	  // 1) we have a list of all annots/intervals here
	  // std::map<std::pair<std::string,std::string> , std::set<interval_t> > atypes;
	  
	  // 2) we've stored all mask-aware rpeaks in tps

	  //  - intersect this with annotation events to make
	  //    new subsetted rpeaks_t objects to pass to
	  //    rr_intervals_t - for time-domain analyses only
	  //    as these will contain gaps
	  //  - we can remove those RR intervals rather than mean
	  //    impute now, if only doing time-domain stuff(?)

	  
	  hrv_opt_t opt2 = opt;
	  opt2.freq_domain = false;
	  
	  std::map<std::pair<std::string,std::string> , std::set<interval_t> >::const_iterator aa = atypes.begin();
	  while ( aa != atypes.end() )
	    {
	      // get subset of peaks : i.e. only tps that fall within an interval
	      //   tps        = std::set<uint64_t>
	      //   aa->second = std::set<interval_t>

	      rpeaks_t peaks2 = rpeaks_t::intersect( tps , aa->second );
	      
	      //
	      // ~not enough data?
	      //
	      
	      if ( peaks2.R_t.size() < 10 )
		{ ++aa;
		  continue;
		}
	      
	      //
	      // Handle outputs
	      //
	      
	      writer.level( aa->first.first , globals::annot_strat );

	      if ( opt.inst_stratify )
		writer.level( aa->first.second , globals::annot_instance_strat );

	      //
	      // HRV analysis: 
	      //

	      rr_intervals_t rr( peaks2 , opt2 );
	      
	      rr.res.write( opt, true );

	      //
	      // Next annot class
	      //
	      
	      ++aa;
	    }
	  
	  writer.unlevel( globals::annot_strat );

	  if ( opt.inst_stratify )
	    writer.unlevel( globals::annot_instance_strat );
	  
	}

      
      //
      // Next channel
      //
    }

  writer.unlevel( globals::signal_strat );
  
}


rr_intervals_t::rr_intervals_t( const rpeaks_t & pks ,
				const hrv_opt_t & opt )
{

  // given a set of peaks, get RR intervals (in ms units)
  //   RR[i] = R[i+1] - R[i]  
  
  // first, excluding any intervals that are too large or too small and
  // replace w/ the mean of the remaining intervals

  t.clear();
  tp.clear();
  rr.clear();
  
  const int np = pks.R_t.size();

  if ( np < 10 ) return;
  
  for (int i=0; i<np-1; i++)
    {
      const double rr1 = globals::tp_duration * ( pks.R_t[i+1] - pks.R_t[i] );
      if ( rr1 >= opt.rr_lwr && rr1 <= opt.rr_upr ) rr.push_back( rr1 );
    }
  
  const int len_imputed = rr.size();
  const double prop_imputed = len_imputed / (double)(np-1);

  if ( len_imputed < 10 ) 
    {
      logger << "  warning: epoch with <10 NN-intervals detected, skipping\n";
      return;
    }
  
  const double mean_rr = MiscMath::mean( rr );
  
  // go back and make a more continuous time series but putting mean
  //  RR in for unrealistic values i.e. this will reduce HRV
  //  presumably a bit, but should allow for spanning very gaps
  //  (i.e. over minutes)
  
  t.clear();
  rr.clear();
  tp.clear();
  
  for (int i=1; i<np; i++)
    {
      const double rr1 = globals::tp_duration * ( pks.R_t[i] - pks.R_t[i-1] );
      tp.push_back(  pks.R_t[i] );

      if ( rr1 >= opt.rr_lwr && rr1 <= opt.rr_upr )
	{
	  rr.push_back( 1000.0 * rr1 );
	  t.push_back( i == 1 ? 0 : t[i-2] + rr1 );
	}
      else
	{
	  rr.push_back( 1000.0 * mean_rr );
	  t.push_back( i == 1 ? 0 : t[i-2] + mean_rr );	  
	}      
    }

  

  //
  // Median filter 
  //

  if ( opt.median_filter_width != 0 ) 
    rr = MiscMath::median_filter( rr , opt.median_filter_width );
  
    
  //
  // Time-domain stats
  //

  res.NP = res.NP_TOT = pks.npks;
  
  // if not enough points in this epoch
  if ( pks.npks < 10 ) return;   
  
  res.P_INV = pks.p_inverted;
  res.INV = pks.inverted;
  
  res.RR = 1000.0 * mean_rr;
  res.HR = 60 / mean_rr;
  res.IMPUTED = prop_imputed;


  if ( opt.freq_domain )
    {
  
      //
      // Cubic spline
      //
      
      // t will now start at 0 and proceed w/out
      // major gaps (in seconds)
      
      tk::spline spline;
      spline.set_points( t, rr );
      
      const int    hrv_sr  = 4;
      
      const double tmax = t[ t.size() - 1 ];    
      const double tinc = 1.0 / (double)hrv_sr;
      
      std::vector<double> rri;
      std::vector<double> ti;
      for (double tt = 0 ; tt <= tmax ; tt += tinc )
	{
	  ti.push_back( tt ); 
	  rri.push_back( spline(tt) );
	}
      
      // resample intervals to create a uniform grid:
      //  cubic spline, 4 Hz  
      //  welch using e.g. 512 sample segment length 
      //  but reduce if segment too short (i.e. not enough RR intervals) 

      const int total_points    = rri.size();  
      const int segment_points  = opt.welch_nsamples > rri.size() ? rri.size() : opt.welch_nsamples ; // segment_sec * sr;
      const int noverlap_points = segment_points / 2 ; // overlap_sec * sr;
      
      const double segment_sec  = segment_points / (double)hrv_sr ; 
      const double overlap_sec  = segment_sec / 2.0 ; 
      
      int noverlap_segments = floor( ( total_points - noverlap_points) 
				     / (double)( segment_points - noverlap_points ) );      
      
      window_function_t window_function = WINDOW_TUKEY50;

      // remove mean
      double rri_mean = MiscMath::mean( rri );
      for (int i=0; i<rri.size(); i++)
	rri[i] -= rri_mean;      
            
      // for (int i=0; i<ti.size(); i++)
      // 	std::cout << "rr" << "\t"
      //  		  << ti[i] << "\t"
      //  		  << rri[i] << "\n";
      
      // Welch
      PWELCH pwelch( rri ,
		     hrv_sr , 
		     segment_sec ,
		     noverlap_segments ,
		     window_function );
      
      bin_t bin( 0 , 0.5 , 1 );
      bin.bin( pwelch.freq , pwelch.psd );
      
      bool bad = false;
      for ( int i = 0 ; i < bin.bfa.size() ; i++ )
	{
	  if ( bin.bspec[i] <= 0 )
	    {
	      bad  = true;		       
	      bin.bspec[i] = 1e-4 ; // set to -40dB as a fudge		   
	    }
	}
      
      // for ( int i = 0 ; i < bin.bfa.size() ; i++ )
      //  	std::cout << "HRV\t"
      //  		  << bin.bfa[i] << "\t" 
      //  		  << bin.bspec[i] << "\t"
      //  		  << 10*log10( bin.bspec[i] ) << "\n";
      
      // 7) Power ratios
      //  ULF	< 0.0033	Uncertain
      
      //  VLF	0.0033–0.04	Circadian/metabolic
      //    rhythms with periods between 25 and 300 secs
      
      //  LF	0.04–0.15	Sympathetic + parasympathetic
      //   rhythms with periods between 7 and 25 s and is affected by breathing from ~3 to 9 bpm
      //   within a 5 min sample, 12–45 complete periods of oscillation
      
      //  HF	0.15–0.4	Parasympathetic (vagal)
      //   breathing from 9 to 24 bpm
      
      //  LF/HF ratio	--	Autonomic balance index
      //   ratio between sympathetic nervous system (SNS) and
      //   parasympathetic nervous system (PNS) activity (under controlled conditions)

      //  Total power:
      //     sum of the energy in the ULF, VLF, LF, and HF bands for 24 h
      //     VLF, LF, and HF bands for short-term recordings


      double vlf = 0 , lf = 0 , hf = 0;
      int n_vlf = 0 , n_lf = 0 , n_hf = 0;
      int lf_peak = 0 , hf_peak = 0;
      double lf_max = 0  , hf_max = 0;

      // uniform freq. but still, just take mean..
      // [ to swap w/ triangulation in any case.]
      std::vector<double> vlf_fd, lf_fd, hf_fd;
      
      for ( int i = 0 ; i < bin.bfa.size() ; i++ )
	{
	  
	  // if ( bin.bfa[i] >= 0.0033 && bin.bfa[i] < 0.04 )
	  //   {
	  //     vlf += bin.bspec[i];
	  //     if ( i != 0 ) 
	  // 	vlf_fd.push_back( bin.bfa[i] - bin.bfa[i-1] );
	  //     ++n_vlf;
	  //   }

	  if ( bin.bfa[i] >= 0.04 && bin.bfa[i] < 0.15 )
	    {
	      lf += bin.bspec[i];
	      ++n_lf;
	      if ( bin.bspec[i] > lf_max ) {
		lf_max = bin.bspec[i];
		lf_peak = i;
	      }
	      lf_fd.push_back( bin.bfa[i] - bin.bfa[i-1] );
	    }
	  else if ( bin.bfa[i] >= 0.15 && bin.bfa[i] < 0.4 )
	    {
	      hf += bin.bspec[i];
	      ++n_hf;
	      if ( bin.bspec[i]	> hf_max ) {
		hf_max = bin.bspec[i];
                hf_peak = i;
	      }
	      hf_fd.push_back( bin.bfa[i] - bin.bfa[i-1] );
	    }	  
	}

      //vlf *= MiscMath::mean( vlf_fd );
      lf *= MiscMath::mean( lf_fd );
      hf *= MiscMath::mean( hf_fd );


      // store
      res.LF = lf;
      res.LF_N = lf / ( lf + hf );
      res.LF_PK = bin.bfa[ lf_peak ] ;
      res.HF = hf;
      res.HF_N = hf / ( lf + hf );
      res.HF_PK = bin.bfa[ hf_peak ] ;
      res.LF2HF = lf / hf ;
      
    }
   
  // req.
  // Accept range:  300 ms ≤ RR ≤ 2000 ms   (i.e., 30–200 bpm)
  // Some studies are stricter:  400 ms ≤ RR ≤ 1600 ms  (i.e., 37.5–150 bpm)
  // drop these prior to interpolation

  
  //

  // HRV on different time-scales
  // 24 h
  // short-term (~5 min),
  // ultra-short-term (<5 min)

  // NN - 'normal' RR intervals (i.e. as above) 
  //
  // Time domain HRV metrics
  //

  // SDNN: Standard deviation of NN intervals
  // SDANN: Standard deviation of the average NN intervals for each 5 min segment of a 24 h HRV recording
  // SDNNI: Mean of the standard deviations of all the NN intervals for each 5 min segment of a 24 h HRV recording
  // pNN50: Percentage of successive RR intervals that differ by more than 50 ms


  // RMSSD: Root mean square of successive RR interval differences
  //   Reflects Parasympathetic Activity
  //   RMSSD is sensitive to short-term variations in heart rate and is primarily influenced by the parasympathetic nervous system.
  //  Indicator of Recovery and Stress:
  //   Higher RMSSD values generally indicate better recovery and lower stress levels,
  //   while lower values may suggest stress or fatigue.


  const int npks = rr.size();
  
  res.SDNN = MiscMath::sdev( rr );
  res.SDNN_R = MiscMath::sdev_robust( rr );  
  res.pNN50 = 0;
  res.RMSSD = 0;

  std::vector<double> diffs;
  for (int i=1;i<npks;i++)
    {
      const double df = fabs( rr[i] - rr[i-1] );
      if ( df > 50 ) ++res.pNN50;
      res.RMSSD += df*df;
      diffs.push_back( df );
    }
  
  res.RMSSD /= (double)(npks-1.0);
  res.RMSSD = sqrt( res.RMSSD );  
  res.RMSSD_R = MiscMath::sdev_robust( diffs );  
  res.pNN50 /= (double)npks;

    
}


hrv_res_t hrv_res_t::summarize( const std::vector<hrv_res_t> & x )
{
  const int n = x.size();
  hrv_res_t res;
  
  for (int i=0; i<n; i++)
    {
      res.IMPUTED += x[i].IMPUTED;
      res.P_INV += x[i].P_INV;
      res.INV += x[i].INV;
      res.NP += x[i].NP;
      res.NP_TOT += x[i].NP_TOT;
      
      res.RR += x[i].RR;
      res.HR += x[i].HR;
      
      res.SDNN += x[i].SDNN;
      res.SDNN_R += x[i].SDNN_R;
      res.RMSSD += x[i].RMSSD;
      res.RMSSD_R += x[i].RMSSD_R;
      res.pNN50 += x[i].pNN50;
      
      res.LF += x[i].LF;
      res.HF += x[i].HF;
      res.LF_N += x[i].LF_N;
      res.HF_N += x[i].HF_N;
      res.LF_PK += x[i].LF_PK;
      res.HF_PK += x[i].HF_PK;
      res.LF2HF += x[i].LF2HF;	
    }
  
  res.IMPUTED /= (double)n;
  res.P_INV /= (double)n;
  res.INV /= (double)n;

  res.NP /= (double)n;
  // i.e. do not divide NP_TOT by n
  
  res.RR  /= (double)n;
  res.HR  /= (double)n;
  
  res.SDNN  /= (double)n;
  res.SDNN_R  /= (double)n;
  res.RMSSD  /= (double)n;
  res.RMSSD_R  /= (double)n;
  res.pNN50  /= (double)n;
  
  res.LF  /= (double)n;
  res.LF_N  /= (double)n;
  res.LF_PK  /= (double)n;
  res.HF  /= (double)n;
  res.HF_N  /= (double)n;
  res.HF_PK  /= (double)n;
  res.LF2HF  /= (double)n;


  // **assumes** that all strata have the same # 
  // std::map<std::string,double> strat_RR;
  // std::map<std::string,double> strat_SDNN_R;
  // std::map<std::string,double> strat_RMSSD_R;
  // std::map<std::string,double> strat_pNN50;
  
  return res;
}

  
void hrv_res_t::write( const hrv_opt_t & opt, const bool reduced ) const
{
  // reduced == T for annot-stratifed analyses:
  //   i.e. we'll drop freq. domain and some other metrics
  
  if ( ! reduced ) {  
    writer.value( "IMPUTED" , IMPUTED);
    writer.value( "P_INV" , P_INV);
    writer.value( "INV" , INV);
  }
  
  writer.value( "NP" , NP);

  // only for epoch-summarized view
  if ( NP_TOT > NP ) 
    writer.value( "NP_TOT" , NP_TOT);  

  writer.value( "RR" , RR );
  writer.value( "HR" , HR );

  if ( opt.time_domain ) { 
    writer.value( "SDNN" ,SDNN );
    writer.value( "SDNN_R" ,SDNN_R );
    writer.value( "RMSSD" ,RMSSD );
    writer.value( "RMSSD_R" ,RMSSD_R );
    writer.value( "pNN50" ,pNN50 );
  }
  
  if ( ! reduced )
    {
      if ( opt.freq_domain )
	{
	  writer.value( "LF" ,LF );
	  writer.value( "HF" ,HF );
	  
	  if ( HF + LF > 0 ) { 
	    writer.value( "LF_N" ,LF_N );
	    writer.value( "HF_N" ,HF_N );
	  }
	  
	  if ( LF > 0 ) 
	    writer.value( "LF_PK" ,LF_PK );
	  
	  if ( HF > 0 ) {
	    writer.value( "HF_PK" ,HF_PK );
	    writer.value( "LF2HF" ,LF2HF );
	  }
	}
    }
  
}


rpeaks_t rpeaks_t::intersect( const std::set<uint64_t> & x ,
			      const std::set<interval_t> & y )
{
  rpeaks_t res;
  for (const auto& interval : y) {
    auto it = x.lower_bound(interval.start);  
    while (it != x.end() && *it < interval.stop) {
      res.R_t.push_back( *it );
      ++it;
    }
  }
  res.npks = res.R_t.size();
  return res;
}

 
