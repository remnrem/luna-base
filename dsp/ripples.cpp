
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

#include "dsp/hilbert.h"

#include "helper/logger.h"
#include "db/db.h"

extern logger_t logger;
extern writer_t writer;

void dsptools::ripple_wrapper( edf_t & edf , param_t & param )
{
  
  // Signals
  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  
  if ( signals.size() == 0 ) return;

  const int ns = signals.size();

  //
  // Analysis parameters
  //
  
  const double flwr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 70 ;

  const double fupr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 150 ;

  const double kwin_ripple = param.has( "ripple" ) ? param.requires_dbl( "ripple" ) : 0.02 ;

  const double kwin_tw = param.has( "tw" ) ? param.requires_dbl( "tw" ) : 5 ;
    
  const double th = param.has( "th" ) ? param.requires_dbl( "th" ) : 90 ;

  if ( th <= 0 || th >= 100 ) Helper::halt( "expecting th percentile between 0 - 100% " );
  
  const double req_msec = param.has( "msec" ) ? param.requires_dbl( "msec" ) : 6.0; 

  const int req_peaks_flt = param.has( "peaks" ) ? param.requires_int( "peaks" ) : 6;

  const int req_peaks_raw = param.has( "peaks-raw" ) ? param.requires_int( "peaks-raw" ) : req_peaks_flt ; 

  const int hfbands = param.has( "bands" ) ? param.requires_int( "bands" ) : 1 ;
  

  //
  // Outputs
  //

  const std::string annot_label = param.has( "annot" ) ? param.value( "annot" ) : "" ;
  
  //
  // Check samples rates 
  //

  const double min_nyquist = fupr / 2.0;
  std::vector<double> Fs = edf.header.sampling_freq( signals );  
  int sr = Fs[0];
  for (int s=0;s<ns;s++)
    {
      if ( Fs[s] < min_nyquist )
	Helper::halt( "sample rate not sufficient for f-upr" );
      if ( Fs[s] != sr )
	Helper::halt( "all sampling rates must be similar for PSI" );
    }

  std::vector<uint64_t> * tp = NULL ;
  
  for (int s=0; s<ns; s++)
    {

      // skip any annotation channels
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      // get data, whole signals
      slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );

      const std::vector<double> * d = slice.pdata();

      //if ( tp == NULL )
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      
      // detect ripples
      ripples_t ripples( *d , *tp , sr , flwr, fupr, kwin_ripple, kwin_tw, hfbands, 
			 th, req_msec, req_peaks_flt, req_peaks_raw );
      
    }
  
}



ripples_t::ripples_t( const std::vector<double> & x ,
		      const std::vector<uint64_t> & tp ,
		      const int sr ,
		      const double flwr ,
		      const double fupr ,
		      const double kwin_ripple ,
		      const double kwin_tw ,
		      const int hfbands , 
		      const double th ,		      
		      const double req_msec ,
		      const int req_peaks_flt ,
		      const int req_peaks_raw )
{

  const bool verbose = true ; 
  const bool otsu = false;
  const double combine_msec = 10.0;
  const uint64_t combine_tp = ( combine_msec / 1000.0 ) *  globals::tp_1sec ;
  
    
  //
  // set up 
  //

  const int n = x.size();

  if ( tp.size() != n )
    Helper::halt( "internal error in ripples, #1" );

  ripples.clear();

  
  //
  // filter-Hilert
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
      logger << "  filtering " << f1 << " -- " << f2 << "\n";
      
      std::vector<double> txf = dsptools::apply_fir( x , sr , fir_t::BAND_PASS ,
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
      std::cout << "tmag n = " << tmag.size() << "\n";
      
      // normalize
      tmag = MiscMath::Z( tmag );
      
      // aggregate
      for (int i=0; i<n; i++) mag[i] += tmag[i];
      
    }

  // aggrgetae
  for (int i=0; i<n; i++) mag[i] /= (double)hfbands;

  // unit scale
  double mmin, mmax;
  MiscMath::minmax( mag , &mmin, &mmax );
  const double mrng = mmax - mmin == 0 ? 1 : mmax - mmin ;  
  
  for (int i=0;i<n;i++)
    {
      if ( mag[i] <= mmin ) mag[i] = 0;
      else if ( mag[i] >= mmax ) mag[i] = 1;
      else mag[i] = ( mag[i] - mmin ) / mrng;
    }
  
  //
  // threshold
  //

  if ( otsu )
    {
      
      // scale by median
      const double med = MiscMath::median( mag );
      
      std::vector<double> adj_vals(  mag.size() );
      for (int i=0;i<adj_vals.size();i++)
	adj_vals[i] = mag[i] / med; 
      
      std::map<double,double> tvals;
      
      double thf = 0;
      
      const double otsu_th = MiscMath::threshold( adj_vals , 1 , 5 , 0.2 , &thf, &tvals );
      
      logger << "  estimated empirical thresholds as " << otsu_th << "\n";
      
      writer.var( "EMPTH" , "Empirical threshold" );
      writer.value( "EMPTH" , otsu_th );

      writer.var( "EMPF" , "Empirical threshold frequency" );
      writer.value( "EMPF" , thf );

      std::map<double,double>::const_iterator tt =  tvals.begin();
      while ( tt != tvals.end() )
	{
	  std::cout << " tt " << tt->first << " = " << tt->second << "\n";
	  ++tt;
	}

      const double thxr = MiscMath::percentile(adj_vals , th / 100.0 );
      logger << " percentile (rel) th = " << th << " = " << thxr << "\n";

    }
  
  
  const double thx = MiscMath::percentile( mag , th / 100.0 );
  
  logger << " percentile th = " << th << " = " << thx << "\n";

  
  //
  // detect putative ripples
  //
  
  std::vector<ripple_t> all_ripples;
  
  int start = -1;
  bool in_ripple = false;

  // to detect discontinuities
  uint64_t dt = 1.5 * ( globals::tp_1sec / sr ) ; // time in tp-units  


  //
  // edges
  //

  std::vector<bool> edge( n , false );
  edge[0] = edge[n-1] = true;

  // discontinuities
  for (int i=1; i<n; i++)
    {
      if ( ( tp[i] - tp[i-1] ) > dt )
	{
	  // fill in edges...
	  edge[i] = true;
	}
    }
  
    

  //
  // iterate over signal
  //

  for (int i=0; i<n; i++)
    {
      //std::cout << i << "\t" <<  (*mag)[i] << "\t" << thx << "\t" << in_ripple << "\n";
      
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

	  // can be i==0 here, so okay for i-1 idx
	  const bool disc = ( tp[i] - tp[i-1] ) > dt ; 
	  //std::cout << " disc = " << disc << "\t" << tp[i] - tp[i-1]  << "\t" << dt << "\n";
	  
	  // discontinuity, or end of sequence, 
	  if ( disc || i == n-1 || mag[i] < thx )
	    {

	      // if end, add +1
	      if ( i == n-1 ) ++i;

	      // nb. i is one past end
	      const int & stop = i;

	      // meets duration criterion?
	      const int len_sp = stop - start ; 
	      const double len_msec = ( len_sp / (double)sr ) * 1000.0 ; 
	      bool okay = len_msec >= req_msec ;
	      
	      // peak count okay?
	      if ( okay && req_peaks_flt )
		{
		  int flt_peakn = 0;
		  for (int s=start+1; s<stop-1; s++)
		    if ( fabs( xf[s] ) > fabs( xf[s-1] )  && fabs( xf[s] ) > fabs( xf[s+1] ) )
		      ++flt_peakn;
		  if ( flt_peakn < req_peaks_flt )
		    okay = false;
		}
	      
	      // same, in raw data
	      if ( okay && req_peaks_raw )
		{
		  int raw_peakn = 0;
		  for (int s=start+1; s<stop-1; s++)
		    if ( fabs(x[s]) > fabs(x[s-1]) && fabs(x[s]) > fabs(x[s+1]) )
		      ++raw_peakn;
		  if ( raw_peakn < req_peaks_raw )
		    okay = false;
		  
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


  logger << " found " << all_ripples.size() << " ripples, merged to " << ripples.size() << "\n";


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


  
  for (int r=0; r<2; r++)
    {
      
      int s1 = ripples[r].start_sp;
      int s2 = ripples[r].stop_sp;

      s1 -= 400;
      s2 += 400;

      if ( s1 < 0 ) s1 = 0;
      
      for (int s=s1; s<s2; s++)
	std::cout << "r" << r+1 << "\t" << s << "\t" << mag[s] << "\t" << x[s] << "\t" << xf[s] << "\n";
      std::cout << "\n";
    }
  
}


