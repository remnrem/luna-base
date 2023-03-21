
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

#include "eval.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

#include <iomanip>
#include <cmath>

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

  //
  // find ECG peaks
  //

  rpeaks_t peaks = mpeakdetect( edf , ecg , tp , sr );
  

  //
  // clean beats
  //
  
  int removed = peaks.clean( 0.3 , 2 );  
    
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
 

rpeaks_t dsptools::mpeakdetect( const edf_t & edf , 
				const std::vector<double> * d , 
				const std::vector<uint64_t> * tp , 				
				int Fs , 
				const std::vector<double> * eeg , 
				bool * force )
{
  
  const int n = d->size();

  if ( tp->size() != n ) Helper::halt("error in mpeakdetect");

  rpeaks_t peaks( n );
  
  // mean-center ECG
  std::vector<double> x = MiscMath::centre( *d );
  
  // band-pass filter: 0.05-40Hz 
    
  std::vector<double> bpf = dsptools::apply_fir( x , Fs , fir_t::BAND_PASS ,
						 1, // Kaiser window
						 0.02 , 0.5 , // ripple, TW
						 0.5 , 40 ) ;  // f1, f2

  // differentiate and sqr data
  std::vector<double> sq( n - 1 );
  for (int i=0;i<n-1;i++) 
    sq[i] = ( bpf[i+1] - bpf[i] ) * ( bpf[i+1] - bpf[i] );

  // integrate over 7 points (i.e. sum)  
  // d=[1 1 1 1 1 1 1]; % window size - intialise
  // x = filter(d,1,sqr);
  
  std::vector<double> ss( n - 1 , 0 );
  for (int i=0; i < n-1; i++)
    for (int j = 0 ; j < 7 ; j++ )
      if ( i-j >= 0 ) ss[i] += sq[i-j];

  // median filter, window size of 10
  std::vector<double> mf = MiscMath::median_filter( ss , 10 );

  // remove filter delay
  int delay = ceil( 7 / 2.0 );
  
  // i.e. skip first 'delay-1' elements, by starting at element 'delay-1', 
  // put back into 'sq'
  
  std::vector<double> mdfint;
  for (int i=delay-1;i<ss.size();i++) mdfint.push_back( ss[i] );

  //
  // segment search area (on 'sq')
  //
  
  //  max_h = max (mdfint(round(len/4):round(3*len/4)));
  int len = mdfint.size();
  
  // int s1 = round( len/4.0 );
  // int s2 = round( 3*len/4.0 );
  
  // avoid issues, calculate max within each EPOCH, take median max value 
  int e30 = Fs * 30;
  int ne = len / e30;
  
  std::vector<double> maxvals;
  for (int e = 0; e < ne ; e++)
    {
      int s1 = e * e30;
      int s2 = s1 + e30 - 1 ;  
      double max_h = 0;
      for (int i = s1; i <= s2 ; i++ ) 
	if ( mdfint[i] > max_h ) max_h = mdfint[i];
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
      bool pk = mdfint[i] > th ;
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
  //  do minima precede maxima?
    
  int inv1 = 0;
  for (int i=0;i<minloc.size();i++)
    if ( minloc[ i ] < maxloc[ i ] ) ++inv1;


  //  bool inverted = minloc[ minloc.size()-1 ] < maxloc[ maxloc.size() - 1 ];

  bool inverted = inv1 / (double)minloc.size() > 0.5;   
  
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
      
      //std::cout <<std::setprecision(16) << "HB " << (double)ht << "\n"; 
      //std::cout << "HB " << (double)ht * globals:: tp_duration << " secs\n"; 
      //std::cout << "xx " << ( inverted ? minloc[i] : maxloc[i] ) << "\n"; 
    }
  
  // take the one with the smallest HRV  (threshold at 2)
  const int nb = R_i.size();
  
  const uint64_t two_sec = globals::tp_1sec + globals::tp_1sec; 
  const uint64_t point3_sec = 0.3 * globals::tp_1sec;
  
  std::vector<double> v1, v2;
  for (int i=1;i<nb;i++)
    {
      double i1 = R_t[i] - R_t[i-1];
      if ( i1 > point3_sec && i1 < two_sec ) v1.push_back( i1 * globals::tp_duration );

      double i2 = t_R[i] - t_R[i-1];
      if ( i2 > point3_sec && i2 < two_sec ) v2.push_back( i2 * globals::tp_duration ); 
    }
  
  double var1 = MiscMath::variance( v1 );
  double var2 = MiscMath::variance( v2 );

  //
  // Output
  //

  writer.value( "ECG_INVERTED" , inverted );
  writer.value( "ECG_VAR1" , var1 );
  writer.value( "ECG_VAR2" , var2 );
  writer.value( "ECG_VAR_FLIPPED" , var1 < var2 ? 0 : 1 );

  //
  // Based on variance if inter-beat intervals, select the orientation that
  // minimizes this value
  //
  
  if ( var1 < var2 )
    {
      // swap in (best guess)
      peaks.R_t = R_t;
      peaks.R_i = R_i;
    }
  else
    {
      peaks.R_t = t_R;
      peaks.R_i = i_R;
    }

  // mirror
  // mirror->R_t = t_R;
  // mirror->R_i = i_R;

  // summarize
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

int rpeaks_t::strip( const std::vector<interval_t> & bad_epochs )
{
  // return # of peaks removed

  int removed = 0;
  
  

  return removed;
}
