
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
  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  
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

  const int hfbands = param.has( "bands" ) ? param.requires_int( "bands" ) : 1 ;

  const double edge_secs = param.has( "edges" ) ? param.requires_dbl( "edges" ) : 1.0 ; 

  const double combine_msec = param.has( "combine" ) ? param.requires_dbl( "combine" ) : 10.0;
  
  
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

      // output
      writer.level( signals.label(s) , globals::signal_strat );
      
      // get data, whole signals
      slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );

      const std::vector<double> * d = slice.pdata();

      //if ( tp == NULL )
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      
      // detect ripples
      ripples_t ripples( *d , *tp , sr , flwr, fupr, kwin_ripple, kwin_tw, verbose ,
			 hfbands, th, req_msec, req_peaks_flt, req_peaks_raw , combine_msec,			 
			 edge_secs , otsu_k );

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
		      const double combine_msec , 
		      const double edge_secs ,
		      const int otsu_k )
{

  sr = sr_;
  verbose = verbose_; 
  const bool otsu = otsu_k != -1;
  const uint64_t combine_tp = ( combine_msec / 1000.0 ) *  globals::tp_1sec ;
      
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

  
  logger << "  found " << all_ripples.size() << " ripples, merged to " << ripples.size() << "\n";
  

  //
  // add magnitudes
  //

  for (int i=0; i< ripples.size(); i++)
    {

      ripple_t & rip = ripples[i];
      
      rip.n = rip.stop_sp - rip.start_sp;
      
      //
      // mean magnitude
      //
      
      rip.x = 0;
      
      for (int j= rip.start_sp; j < rip.stop_sp; j++)
	rip.x += mag[j];
      
      rip.x /= (double)(rip.stop_sp - rip.start_sp);
            
      //
      // frequency (from filtered signal, xf) half-wave sample points
      //
      
      std::vector<int> hwsp;

      // based on ZC:

      int last = -1;
      
      for (int j = rip.start_sp == 0 ? 1 : rip.start_sp ; j < rip.stop_sp; j++)
	{
	  const bool zc = ( xf[j-1] <= 0 && xf[j] > 0 ) || ( xf[j-1] > 0 && xf[j] <= 0 );
	  
	  if ( zc )
	    {

	      //std::cout << " ZC = " << j << "\n";

	      if ( last != -1 )
		hwsp.push_back( j - last + 1 );	      
	      last = j;
	    }
	}
     
      
      rip.frq = sr / (double)( 2.0 * MiscMath::mean( hwsp ) ) ;

      //      std::cout	<< "frq\t" << rip.frq << " (p=" << hwsp.size() << ") ";

      //
      // based on Peaks
      //
      
      hwsp.clear();
      last = -1;
      for (int j = rip.start_sp + 1 ; j < rip.stop_sp - 1 ; j++)
	{
	  
	  const bool pos_peak = xf[j] >= xf[j-1] && xf[j] > xf[j+1] ;
	  const bool neg_peak = xf[j] <= xf[j-1] && xf[j] < xf[j+1] ;	  
	  
          if ( pos_peak || neg_peak )
            {
	      //std::cout << " PEAK = " << j << "\n";
	      if ( last != -1 )
		hwsp.push_back( j - last + 1 );
              last = j;
            }
        }


      rip.frqp2p = sr / (double)( 2.0 * MiscMath::mean( hwsp ) ) ;

      
      // std::vector<double> xx;
      // for (int j = rip.start_sp ; j < rip.stop_sp  ; j++)
      // 	xx.push_back( xf[j] );
      
      // std::cout << rip.frqp2p << " (p=" <<  hwsp.size() << ")  ; mean = " << MiscMath::mean(xx) << "\n";
      
      // for (int j = rip.start_sp ; j < rip.stop_sp  ; j++)
      // 	{
      // 	  std::cout << "  " << j << "\t" << ( j - rip.start_sp ) / (double)sr << "\t" << xf[j] << "\n";
      // 	}
      
      
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

	  writer.value( "START_SP" , ripple.start_sp );
	  writer.value( "STOP_SP" , ripple.stop_sp );

	  writer.value( "PCT" , ripple.wgt );
	  writer.value( "FRQ" , ripple.frq );
	  writer.value( "FRQ2" , ripple.frqp2p );
	  writer.value( "MAG" , ripple.x );
	  writer.value( "SP" , ripple.n );
	  writer.value( "DUR" , ripple.n / (double)sr );
	  
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
      instance->set( "frq2" , ripple.frqp2p );
      instance->set( "n" , ripple.n );
      instance->set( "mag" , ripple.x );      
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
