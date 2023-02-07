
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

#include "dsp/ripples.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "annot/annot.h"

#include "dsp/hilbert.h"

#include "helper/logger.h"
#include "db/db.h"

extern logger_t logger;
extern writer_t writer;

void dsptools::ripple_wrapper( edf_t & edf , param_t & param )
{
  
  // Signals

  const bool no_annotations = true;

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , no_annotations );
  
  if ( signals.size() == 0 ) return;

  const int ns = signals.size();

  //
  // Analysis parameters
  //
  
  const double flwr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 70 ;

  const double fupr = param.has( "f-upr" ) ? param.requires_dbl( "f-upr" ) : 150 ;

  const double kwin_ripple = param.has( "ripple" ) ? param.requires_dbl( "ripple" ) : 0.02 ;

  const double kwin_tw = param.has( "tw" ) ? param.requires_dbl( "tw" ) : 5 ;
    
  const double th = param.has( "th" ) ? param.requires_dbl( "th" ) : 90 ;

  if ( th <= 0 || th >= 100 ) Helper::halt( "expecting th percentile between 0 - 100% " );
  
  const double req_msec = param.has( "msec" ) ? param.requires_dbl( "msec" ) : 6.0; 

  const int req_peaks_flt = param.has( "peaks" ) ? param.requires_int( "peaks" ) : 6;

  const int req_peaks_raw = param.has( "peaks-raw" ) ? param.requires_int( "peaks-raw" ) : req_peaks_flt ; 

  const double req_raw_p2p_prop = param.has( "peaks-raw-prop" ) ? param.requires_dbl( "peaks-raw-prop" ) : 0.01 ; 

  const double max_amp_thresh_abs = param.has( "max-abs" ) ? param.requires_dbl( "max-abs" ) : -1 ;

  const double max_amp_thresh_pct = param.has( "max-pct" ) ? param.requires_dbl( "max-pct" ) : -1 ; 
  
  const int hfbands = param.has( "bands" ) ? param.requires_int( "bands" ) : 1 ;

  const double edge_secs = param.has( "edges" ) ? param.requires_dbl( "edges" ) : 1.0 ; 
  
  const std::set<std::string> excludes = param.strset( "exclude" ); // annots to exclude, e.g. IEDs
  
  const double combine_msec = param.has( "combine" ) ? param.requires_dbl( "combine" ) : 10.0;


  //
  // Check excludes
  //

  std::set<std::string>::const_iterator ee = excludes.begin();
  while ( ee != excludes.end() )
    {
      annot_t * annot = edf.timeline.annotations( *ee );
      if ( annot == NULL )
	Helper::halt( "could not find annotation " + *ee );
      ++ee;
    }

  
  //
  // Outputs
  //

  const std::string annot_label = param.has( "annot" ) ? param.value( "annot" ) : "" ;
  
  annot_t * annots = annot_label != "" ? edf.timeline.annotations.add( annot_label ) :  NULL  ; 

  const bool verbose = param.has( "verbose" );

  const bool otsu = param.has( "otsu" );

  const int otsu_k = otsu ? 100 : -1;

  
  //
  // Check samples rates 
  //

  const double min_nyquist = fupr / 2.0;
  std::vector<double> Fs = edf.header.sampling_freq( signals );  

  int sr = Fs[0];
  for (int s=0;s<ns;s++)
    {
      if ( Fs[ s ] < min_nyquist )
	Helper::halt( "sample rate not sufficient for f-upr" );
      
      if ( Fs[ s ] != sr )
	Helper::halt( "all sampling rates must be similar for RIPPLES" );
    }

  std::vector<uint64_t> * tp = NULL ;
  
  for (int s=0; s<ns; s++)
    {

      // skip any annotation channels
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      // output
      writer.level( signals.label(s) , globals::signal_strat );
      
      // get data, whole signals
      slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );

      const std::vector<double> * d = slice.pdata();

      //if ( tp == NULL )
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      
      // detect ripples

      logger << "\n  processing " << signals.label(s) << "...\n";
      
      ripples_t ripples( *d , *tp , sr , flwr, fupr, kwin_ripple, kwin_tw, verbose ,
			 hfbands, th, req_msec, req_peaks_flt, req_peaks_raw , req_raw_p2p_prop,
			 max_amp_thresh_abs , max_amp_thresh_pct , 
			 combine_msec, edge_secs ,
			 excludes.size() == 0 ? NULL : &edf , excludes,
			 otsu_k );
      

      // stats & annots
      if ( ! otsu )
	{
	  
	  ripples.output( true );
	  
	  if ( annots ) 
	    ripples.annotate( annots , signals.label(s) );
	  
	}
      
    }

  writer.unlevel( globals::signal_strat );

}



ripples_t::ripples_t( const std::vector<double> & x ,
		      const std::vector<uint64_t> & tp ,
		      const int sr_ ,
		      const double flwr ,
		      const double fupr ,
		      const double kwin_ripple ,
		      const double kwin_tw ,
		      const bool verbose_ , 
		      const int hfbands , 
		      const double th ,		      
		      const double req_msec ,
		      const int req_peaks_flt ,
		      const int req_peaks_raw ,
		      const double req_raw_p2p_prop , 
		      const double max_amp_thresh_abs ,
		      const double max_amp_thresh_pct ,
		      const double combine_msec , 
		      const double edge_secs ,
		      edf_t * edf , 
		      const std::set<std::string> & excludes , 
		      const int otsu_k )
{

  sr = sr_;
  verbose = verbose_; 
  const bool otsu = otsu_k != -1;
  const uint64_t combine_tp = ( combine_msec / 1000.0 ) *  globals::tp_1sec ;
  
  logger << "  excluding edges of segments, for " << edge_secs << " seconds\n"
	 << "  requiring ripples to be at least " << req_msec << " msec\n"
	 << "  combining ripples nearer than " << combine_msec << " msec\n"
	 << "  requiring at least " << req_peaks_flt << " peaks in the filtered signal, " << req_peaks_raw << " in the raw signal\n"
	 << "  splitting range into " << hfbands << " equal bands\n";
  logger << "  FIR tw = " << kwin_tw << ", ripple = " << kwin_ripple << "\n";



  //
  // set up 
  //

  const int n = x.size();
  
  if ( tp.size() != n )
    Helper::halt( "internal error in ripples, #1" );

  ripples.clear();

  
  //
  // Duration (calculate channel specific, but will be common to all)
  //
  
  totdur_mins = ( n / (double)sr ) ;
  
  
  //
  // filter-Hilbert
  //
  
  std::vector<double> xf = dsptools::apply_fir( x , sr , fir_t::BAND_PASS ,
						1 , // Kaiser window
						kwin_ripple , kwin_tw ,
						flwr ,  fupr ,
						0 ,  // order/ignored
						fir_t::HAMMING , // ignored for KW      
						true ) ; // use FFT
  
  
  std::vector<double> mag( n , 0 );
  
  const double fwin = ( fupr - flwr ) / (double)hfbands;
  
  for (int b=0; b<hfbands; b++)
    {
      
      const double f1 = flwr + b * fwin ;
      const double f2 = flwr + (b+1) * fwin ;

      logger << "  filtering " << f1 << " Hz -- " << f2 << " Hz\n";
      
      std::vector<double> txf =
	hfbands == 1 ? xf :	
	dsptools::apply_fir( x , sr , fir_t::BAND_PASS ,
			     1 , // Kaiser window
			     kwin_ripple , kwin_tw ,
			     f1, f2 , 
			     0 ,  // order/ignored
			     fir_t::HAMMING , // ignored for KW
			     true ) ; // use FFT
      
      if ( txf.size() != xf.size() )
	Helper::halt( "internal error in ripples_t(), flt length " );
      
      // Hilbert transform
      hilbert_t hilbert( txf );
      
      std::vector<double> tmag = *hilbert.magnitude();
      
      // normalize
      tmag = MiscMath::Z( tmag );
      
      // aggregate
      for (int i=0; i<n; i++) mag[i] += tmag[i];
      
    }

  //
  // aggregate & unit scale
  //

  for (int i=0; i<n; i++) mag[i] /= (double)hfbands;
  
  // unit scale
  
  double mmin, mmax;
  MiscMath::minmax( mag , &mmin, &mmax );
  double mrng = mmax - mmin == 0 ? 1 : mmax - mmin ;  
  for (int i=0;i<n;i++)
    {
      if ( mag[i] <= mmin ) mag[i] = 0;
      else if ( mag[i] >= mmax ) mag[i] = 1;
      else mag[i] = ( mag[i] - mmin ) / mrng;
    }
  
  //
  // threshold
  //
  
  const double thx = MiscMath::percentile( mag , th / 100.0 );
  
  logger << "  thresholding at percentile = " << th << " (" << thx << ")\n";


  //
  // absolute amplitude threshold exclusion
  //

  double th_amp = max_amp_thresh_abs > 0 ? max_amp_thresh_abs : -1 ;
  
  if ( max_amp_thresh_pct > 0 ) 
    {

      std::vector<double> ax( x.size() );

      for (int i=0; i<x.size(); i++) ax[i] = fabs( x[i] );

      const double tpct = MiscMath::percentile( ax , max_amp_thresh_pct / 100.0 );

      logger << "  " << max_amp_thresh_pct << " percentile = " << tpct << "\n";

      if ( max_amp_thresh_abs > 0 )
	{
	  if ( tpct < th_amp ) th_amp = tpct ; 
	}
      else
	{
	  th_amp = tpct;
	}
    }

  if ( max_amp_thresh_pct > 0 || max_amp_thresh_abs > 0 )
    logger << "  using amplitude threshold of " << th_amp << "\n";
  
  //
  // empirical threshold determination 
  //

  if ( otsu_k != -1 )
    {
      
      if ( verbose )
	{
	  for (int i=0;i<mag.size();i++)
	    std::cout << "raw\t" << mag[i] << "\n";
	}
      
      std::map<double,double> tvals, fvals;

      double empf;
      const double otsu_est = MiscMath::threshold2( mag , &empf, otsu_k , &fvals, &tvals );
      

      logger << "  Otsu threshold estimate : unit scale (raw) = " << otsu_est << "\n";
      logger << "                          : percentile       = " << 100 * empf << "\n";

      writer.value( "EMPTH" , otsu_est );
      writer.value( "EMPF" , empf );
      
      std::map<double,double>::const_iterator tt =  tvals.begin();
      while ( tt != tvals.end() )
	{
	  writer.level( tt->first , "TH" );
	  writer.value("SIGMAB" , tt->second );
	  writer.value("F" , fvals[ tt->first ] );
	  ++tt;
	}
      writer.unlevel( "TH" );
      
            
      return ;
    }  
  
  
  
  //
  // detect putative ripples
  //
  
  std::vector<ripple_t> all_ripples;
  
  int start = -1;
  bool in_ripple = false;

  //
  // to detect discontinuities & edges
  //

  uint64_t dt = 1.5 * ( globals::tp_1sec / sr ) ; // time in tp-units  
  int edge_sp = edge_secs * sr ; 
  std::vector<bool> edge( n , false );

  n_segments = 1;

  for (int i=0; i<n; i++)
    {
      if ( i == 0 || i == n-1 || ( i != 0 && ( tp[i] - tp[i-1] ) > dt ) ) 
	{
	  // gap between i-1 and i ?  (ignore start/stop here)
	  if ( i != 0 && i != n-1 ) 
	    ++n_segments;

	  // start
	  if ( i == 0 )
	    {
	      for (int j=0; j<edge_sp; j++)
		{
		  if ( j == n ) break;
		  edge[j] = true;
		}	  
	    }
	  else if ( i == n-1 ) // stop
	    {
	      int p = i;
	      for (int j=0; j<edge_sp; j++)
		{
		  if ( p < 0 ) break;
		  edge[p--] = true;	      
		}
	    }
	  else // boundary (start & stop )
	    {
	      int p = i-1;
	      for (int j=0; j<edge_sp; j++)
		{
		  if ( p < 0 ) break;
		  edge[p--] = true;	      
		}
	      
	      // forward
	      p = i;
	      for (int j=0; j<edge_sp; j++)
		{
		  if ( p == n ) break;
		  edge[p++] = true;
		}	  
	    }
	}
    }
  

  // get interval region, excluding 'edges'
  incdur_mins = 0;
  for (int i=0; i<n; i++)
    {
      //      std::cout << " det " << floor( i / 1024.0 ) << " " << i << " " << edge[i] << "\n";
      if ( ! edge[i] ) ++incdur_mins;
    }
  
  // sp -> seconds (then -> mins below)
  incdur_mins /= sr;

  writer.value( "MINS_TOT" , totdur_mins / 60.0 );
  writer.value( "MINS" , incdur_mins / 60.0 );
  writer.value( "SECS_TOT" , totdur_mins );
  writer.value( "SECS" , incdur_mins );
  writer.value( "NSEG" , n_segments );
  
  //
  // iterate over signal
  //

  int fail_dur = 0, fail_amp = 0 , fail_flt_hw = 0 , fail_raw_hw = 0, fail_exc_annot = 0; 
  
  
  for (int i=0; i<n; i++)
    {
      
      if ( ! in_ripple ) // start a new ripple ?
	{
	  if ( mag[i] >= thx )
	    {
	      in_ripple = true ;
	      start = i;
	    }
	}
      else // end one? 
	{
	  
	  // end of a ripple?  (discontinuity/edge, end of data, or below threshold signal)
	  if ( edge[i] || i == n-1 || mag[i] < thx )
	    {
	      
	      // if end, add +1
	      if ( i == n-1 ) ++i;

	      // nb. i is one past end
	      const int & stop = i;

	      // meets duration criterion?
	      const int len_sp = stop - start ; 
	      const double len_msec = ( len_sp / (double)sr ) * 1000.0 ; 
	      
	      bool okay = len_msec >= req_msec ;

	      if ( ! okay ) ++fail_dur;
	      
	      // amplitude threshold?
	      if ( okay && th_amp > 0 )
		{
		  for (int s=start; s<stop; s++)
		    if ( fabs( x[s] ) > th_amp )
		      {
			okay = false;
			++fail_amp;
			break;
		      }
		}
	      
	      // peak count okay?
	      if ( okay && req_peaks_flt )
		{

		  // for (int s=start; s<stop; s++)
		  //   std::cout << " " << xf[s] << "\n";
		  // std::cout << "\n";
		  
		  int flt_peakn = 0;
		  for (int s=start+1; s<stop-1; s++)
		    {
		      if ( xf[s] > xf[s-1] && xf[s] > xf[s+1] ) 
			++flt_peakn;
		      else if ( xf[s] < xf[s-1] && xf[s] < xf[s+1] )
			++flt_peakn;
		    }

		  //		  std::cout << " flt_peakn = " << flt_peakn << ", dur = " << len_sp << " " << len_msec << "\n";

		  if ( flt_peakn < req_peaks_flt )
		    {
		      ++fail_flt_hw;
		      okay = false;
		    }
		}
	      
	      // same, in raw data -- and track sizes
	      if ( okay && req_peaks_raw )
		{
		  
		  std::vector<int> pk;
		  for (int s=start+1; s<stop-1; s++)		    
		    {
		      if ( x[s] >= x[s-1] && x[s] > x[s+1] )	
			pk.push_back( s );
                      else if (	x[s] <= x[s-1] && x[s] < x[s+1] )
			pk.push_back( s );		      
		    }
		  
		  std::vector<double> p2p;
		  std::vector<bool> pdir;
		  double max_p2p = 0;
		  for (int p=1;p<pk.size();p++)
		    {
		      double t = fabs( x[pk[p-1]] ) - fabs( x[pk[p]] ) ;
		      if ( t > max_p2p ) max_p2p = t;
		      p2p.push_back( t );
		      pdir.push_back( x[pk[p-1]] >  x[pk[p]] );
		    }
		  
		  // require these to be at least e.g. 20% of the max p2p to be counted
		  // and require at least req_peaks_raw / 2 of each

		  int raw_peakn_pos = 0;
		  int raw_peakn_neg = 0;
		  for (int p=0; p<p2p.size(); p++)
		    {		      
		      if ( fabs( p2p[p] ) >= max_p2p * req_raw_p2p_prop )
			{
			  if ( pdir[p] ) 
			    ++raw_peakn_pos;
			  else
			    ++raw_peakn_neg;
			}
		    }
		  
		  if ( raw_peakn_pos < req_peaks_raw/2.0 || raw_peakn_neg < req_peaks_raw/2.0 )
		    {
		      ++fail_raw_hw;
		      okay = false;		   
		    }
		}

	      // any annotation exclusions?
	      if ( okay && edf != NULL && excludes.size() != 0 )
		{
		  std::set<std::string>::const_iterator ee = excludes.begin();
		  while ( ee != excludes.end() )
		    {
		      // tp[stop] is +1 end already
		      
		      interval_t interval( tp[start] , tp[stop] );

		      // nb. as extract() is a terrible function and currently
		      // only extracts things that are completely within the interval
		      // here we will make a kludge and expand this region to +/- 0.5 secs
		      
		      interval.expand( 0.5 * globals::tp_1sec );
		      
                      annot_t * annot = edf->timeline.annotations( *ee );
		      if ( annot == NULL ) Helper::halt( "could not find annotation " + *ee );
                      annot_map_t events = annot->extract( interval );
                      const bool has_annot = events.size() ;
                      if ( has_annot )
			{
			  okay = false;
			  ++fail_exc_annot;
			  break;
			}
		      ++ee;
		    }
		}
	      
	      // add?
	      if ( okay )
		{		  
		  // tp[i] will be +1 end, so good
		  ripple_t ripple( tp[start] , tp[stop] , start , stop ) ;
		  all_ripples.push_back( ripple );
		}
	      
	      // mark end of ripple either way
	      in_ripple = false;
	      
	    }
	  
	}
	  
    }

  logger << "  " << all_ripples.size() << " ripples retained, "
	 << " failed N: dur = " << fail_dur
	 << ", amp = " << fail_amp	 
	 << ", half-waves (filtered) " << fail_flt_hw
	 << ", half-waves (raw) = " << fail_raw_hw
	 << ", annotated exclusions = " << fail_exc_annot 
	 << "\n";

  if ( all_ripples.size() == 0 ) return;
  
  //
  // merge nearby
  //
  
  ripples.clear();
  
  const int n0 = all_ripples.size();

  int prev = 0;
  
  for (int i=1; i<n0; i++)
    {
            
      // gap between the previous and this one?
      const bool gap = ( all_ripples[ i ].pos.start - all_ripples[ i - 1 ].pos.stop + 1 ) >= combine_tp ;

      if ( gap )
	{
	  // add prior
	  // typically, prev == i, but this handles the case of adding multiple adjoined regions
	  //  (each separated by < 10ms )
	  ripples.push_back( ripple_t( all_ripples[ prev ].pos.start , all_ripples[ i - 1 ].pos.stop ,
                                       all_ripples[ prev ].start_sp , all_ripples[ i - 1 ].stop_sp  ) );

	  // track that the current is now in the previous, to-be-added set
	  prev = i;
	}
      
    }
  
  // add last set
  ripples.push_back( ripple_t( all_ripples[ prev ].pos.start , all_ripples[ n0 - 1 ].pos.stop ,
			       all_ripples[ prev ].start_sp , all_ripples[ n0 - 1 ].stop_sp  ) );
    
  
  logger << "  found " << all_ripples.size() << " ripples, merged to " << ripples.size() << "\n";
  
  
  //
  // add ripple meta-data: magnitude, frequency, mid-point
  //

  for (int i=0; i< ripples.size(); i++)
    {

      ripple_t & rip = ripples[i];

      //
      // number of sample points
      //

      rip.n = rip.stop_sp - rip.start_sp;
      
      //
      // mean magnitude
      //
      
      rip.x = 0;
      
      for (int j= rip.start_sp; j < rip.stop_sp; j++)
	rip.x += mag[j];
      
      rip.x /= (double)(rip.stop_sp - rip.start_sp);
            
      //
      // frequency (from full-range filtered signal, xf) half-wave sample points
      //

      std::vector<double> xx;
      for (int j = rip.start_sp ; j < rip.stop_sp; j++)
	xx.push_back( xf[j] );
      
      // copy 
      std::vector<double> rawx = xx;
      
      // ensure locally mean-centered
      MiscMath::centre( &xx );

      // count zero-crossings
      // (w/ linear interpolation) 
      std::vector<double> hwsp;

      double last = -9;
      double dt = 1.0 / (double) sr;
      const int nxx = xx.size();

      // track which are pos2neg or not
      std::vector<bool> pos2neg;
      std::vector<int> zc_idx;
      
      for (int j=1; j<nxx; j++)
	{
	  
	  const bool neg2pos_zc = xx[j-1] <= 0 && xx[j] > 0 ;
	  const bool pos2neg_zc = xx[j-1] > 0  && xx[j] <= 0 ;
	  
	  if ( neg2pos_zc || pos2neg_zc ) 
	    {
	      // track for determination of ripple peak/middle
	      pos2neg.push_back( pos2neg_zc );
	      zc_idx.push_back( j );

	      // get fractional ZC point
	      double wj = fabs( xx[j] ) / ( fabs( xx[j-1] ) + fabs( xx[j] ) ) ; 
	      double jj = dt * ( j * wj + (j-1)*(1-wj) );
	      
	      if ( last > -1 )
		hwsp.push_back( jj - last );
	      
	      last = jj;
	    }
	}
      
      
      rip.frq = 1.0 / ( 2.0 * MiscMath::mean( hwsp ) ) ; 

      //
      // verbose of (complete) half-waves ... i.e. may be less than # of peaks
      //
      
      rip.nhw = hwsp.size() ;
      
      //
      // skew & kurtosis based on raw signals
      //

      double mean   = MiscMath::mean( rawx );
      double sd     = MiscMath::sdev( rawx , mean );
      rip.skew   = MiscMath::skewness( rawx , mean , sd );
      rip.kurt   = MiscMath::kurtosis( rawx , mean );
      
      //
      // max peak-to-peak amplitude (based on neg-to-pos)
      //
      
      rip.p2pamp = 0;

      double max_neg = 0 , max_pos = 0;
      bool seen_pos = false; bool seen_neg = false;
      
      for (int i=1; i<zc_idx.size(); i++)
	{
	  
	  // from ZC previous to this one
	  const bool pos_hw = pos2neg[i]; // this is end of this HW

	  if ( pos_hw )
	    {
	      seen_pos = true;
	      for (int j=zc_idx[i] ; j>=zc_idx[i-1]; j-- )
		if ( xx[j] > max_pos ) max_pos = xx[j];
	    }
	  else
	    {
	      seen_neg = true;
	      for (int j=zc_idx[i] ; j>=zc_idx[i-1]; j-- )
                if ( xx[j] < max_neg ) max_neg = xx[j];
	    }
	  
	  // eval only neg-peak --> pos-peak
	  
	  if ( seen_neg && seen_pos )
	    {
	      double p2p = max_pos - max_neg;	      
	      if ( p2p > rip.p2pamp )
		rip.p2pamp = p2p;

	      // keep most recent HW peak; clear one before so we can compare
	      // (i-1) vs i,  and now i vs (i+1)
	      if ( pos_hw )
		{
		  seen_neg = false;
		  max_neg = 0;
		}
	      else
		{
		  seen_pos = false;
		  max_pos = 0;
		}
	      
	    }
	  
	}

      
      //
      // get most central negative peak      
      //

      int didx = 0;
      int dbest = -1;
      
      int mid = nxx / 2;

      for (int i=1; i<zc_idx.size(); i++)
	{
	  // only consider negative halfwaves
	  if ( pos2neg[i] ) continue;

	  // mid-point of this negative halfwave (in sp)
	  int j = ( zc_idx[i] + zc_idx[i-1] ) / 2 ;
	  
	  // distance (in sp) from middle of ripple
	  int d = abs( j - mid );
	  
	  if ( dbest == -1 )
	    {
	      dbest = d;
	      didx = 1;	// of second (neg-to-pos) ZC
	    }
	  else
	    {
	      if ( d < dbest )
		{
		  dbest = d;
		  didx = i; // of second (neg-to-pos) ZC
		}	      
	    }
	}

      // mid-point = peak of negative halfwave (in sp) closet to ripple middle
      int mid_sp = ( zc_idx[didx] + zc_idx[didx-1] ) / 2 ;
      double lowest = 999;

      for (int j = zc_idx[didx-1]; j <= zc_idx[didx]; j++)
	if ( xx[ j ] < lowest )
	  {
	    mid_sp = j;
	    lowest = xx[j] ;
	  }
      
      // convert to time-points from EDF start
      
      double mid_fraction = mid_sp / (double) nxx ;
      
      // in time-points (from EDF start)
      rip.midp = rip.pos.start + ( rip.pos.stop - rip.pos.start ) * mid_fraction ;
      
      // std::cout << " mid_sp = " << rip.mid_sp << "\n"
      // 		<< didx << " " << zc_idx[didx] << " " << zc_idx[didx-1] << "\n";
      
      // for (int z=0; z<zc_idx.size(); z++)
      // 	std::cout << "z = " << z << " " << zc_idx[z] << " " << pos2neg[z] << "\n";
      // std::cout << "\n";
      
       // for (int j = 0 ; j < nxx  ; j++)
       // 	std::cout << "  " << j << "\t" << j / (double)sr << "\t" << rip.start_sp + j << "\t" << ( rip.start_sp + j == rip.mid_sp ) << "\t" << xx[j] << "\n";
    
      
      
    }
  
  //
  // score 
  //

  // add wgt/score for each ripple
  //  go from 'th' percentile to '100'
  
  std::map<double,int> cnts;
  for (int i=0;i<ripples.size();i++)
    cnts[ ripples[i].x ]++;

  int cum = 0;
  std::map<double,int>::iterator cc = cnts.begin();
  while ( cc != cnts.end() )
    {
      cc->second += cum ;
      //std::cout << " cum = " << cc->first << " -> " << cc->second << "\n";
      cum = cc->second ;
      ++cc;
    }
 
  
  // i.e. 10, if th = 90
  const double th2 = 100 - th;
  const double nt = ripples.size();
  //  std::cout << "cum = " << cum << " " << nt << "\n" ; 

  for (int i=0;i<ripples.size();i++)
    {
      ripple_t & rip = ripples[i];
      
      double t = cnts[ rip.x ];
      
      rip.wgt = th + th2 * ( t / (double)nt );
      
    }
  
  
  // nb...
  
  // https://www.nature.com/articles/s41467-019-11444-x#Sec10

  // Ripples: We utilized automatic IED detection algorithms to detect
  // prominent discharges, however, this approach does not capture
  // below-threshold epileptiform or artifactual activity in the
  // HFB. Therefore, we imposed additional measures to isolate true
  // ripples from artifactual ripples, which can easily be mistaken
  // when data is analyzed in the frequency domain or after applying a
  // band-pass filter (Fig. 4b). We first narrowed down our search
  // window for ripples, to epochs during cortical spindle events
  // (±1 s) during NREM sleep when no simultaneous IED had been
  // detected. We time-locked the HFB traces relative to those spindle
  // events and determined the strongest HFB peak during any given,
  // artifact-free spindle event. Then we extracted the number of
  // peaks in the raw signal in a 40-ms window around this peak. If
  // equal or more than three peaks were detected, the event was
  // classified as a true ripple, while events that only exhibited one
  // or two peaks were classified as artifactual. Ripples were peak
  // locked (Supplementary Fig. 4c). Note that for group averages
  // (Fig. 4a, b), we followed the convention that ripples were nested
  // in the upwards swing of the sharp wave and hence, we adjusted
  // ripple polarity, which was less informative given the bipolar
  // referencing scheme (Supplementary Fig. 4c). Resulting ripples
  // were then visually inspected and we carefully assessed their
  // spatiotemporal profile relative to other oscillatory events, such
  // as slow oscillations and spindles (Fig. 4d–g).

  
  // for (int r=0; r<2; r++)
  //   {
      
  //     int s1 = ripples[r].start_sp;
  //     int s2 = ripples[r].stop_sp;

  //     s1 -= 400;
  //     s2 += 400;

  //     if ( s1 < 0 ) s1 = 0;
      
  //     for (int s=s1; s<s2; s++)
  // 	std::cout << "r" << r+1 << "\t" << s << "\t" << mag[s] << "\t" << x[s] << "\t" << xf[s] << "\n";
  //     std::cout << "\n";
  //   }
  
}


void ripples_t::output( const bool verbose )
{
  
  writer.value( "N" , (int)ripples.size() );
  writer.value( "DENS" , ripples.size() / (double) ( totdur_mins / 60.0 ) );

  if ( verbose )
    {
      for (int i=0; i<ripples.size(); i++)
	{

	  const ripple_t & ripple = ripples[i];
	  
	  writer.level( i+1 , globals::count_strat );
	  writer.value( "START" , globals::tp_duration * ripple.pos.start );
	  writer.value( "STOP" , globals::tp_duration * ripple.pos.stop );
	  writer.value( "MID" , globals::tp_duration * ripple.midp );
	  
	  writer.value( "START_SP" , ripple.start_sp );
	  writer.value( "STOP_SP" , ripple.stop_sp );
	  
	  
	  writer.value( "PCT" , ripple.wgt );
	  writer.value( "FRQ" , ripple.frq );	  
	  writer.value( "MAG" , ripple.x );
	  writer.value( "SP" , ripple.n );
	  writer.value( "NHW" , ripple.nhw );
	  writer.value( "AMP" , ripple.p2pamp );
	  writer.value( "DUR" , ripple.n / (double)sr );
	  writer.value( "SKEW" , ripple.skew );
	  writer.value( "KURT" , ripple.kurt );
	}

      writer.unlevel( globals::count_strat );
      
    }
    
}

  


void ripples_t::annotate( annot_t * a , const std::string & ch )
{
  
  
  for (int i=0;i<ripples.size();i++)
    {
      const ripple_t & ripple = ripples[i];
      instance_t * instance = a->add( "." , ripple.pos , ch );
      instance->set( "pct" , ripple.wgt );
      instance->set( "frq" , ripple.frq );
      instance->set( "n"   , ripple.n );
      instance->set( "nhw" , ripple.nhw );
      instance->set( "amp" , ripple.p2pamp );
      instance->set( "mag" , ripple.x );
      instance->set( "skew" , ripple.skew );
      instance->set( "kurt" , ripple.kurt );
      
      std::string mid_tp = "tp:" + Helper::int2str( ripple.midp );
      instance->set( "mid" , mid_tp );
    }
}


std::vector<double> ripples_t::percentile( const std::vector<double> & x) 
{

  const int n = x.size();
  
  std::map<double,int> cnts;
  for (int i=0; i<n; i++)
    cnts[ x[i] ]++;
  
  int cum = 0;
  std::map<double,int>::iterator cc = cnts.begin();
  while ( cc != cnts.end() )
    {
      cc->second += cum ;      
      cum = cc->second ;
      ++cc;
    }
  
  std::vector<double> y( n );
  for (int i=0;i<n; i++)
    y[i] = 100 * ( cnts[ x[i] ] / (double)n );
  
  return y;
  
}
