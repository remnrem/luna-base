
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

#include "coherence.h"

#include "resample.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "param.h"
#include "db/db.h"
#include "fftw/cohfft.h"
#include "helper/logger.h"
#include "helper/helper.h"


// notes:
//   http://doc.ml.tu-berlin.de/causality/
//   https://mne.tools/stable/generated/mne.connectivity.spectral_connectivity.html
//    [1] Nolte et al. “Robustly Estimating the Flow Direction of Information in
//     Complex Physical Systems”, Physical Review Letters, vol. 100, no. 23, pp. 1-4, Jun. 2008.


extern writer_t writer;

extern logger_t logger;

#include <cmath>


void dsptools::coherence( edf_t & edf , param_t & param )
{


  std::string signal_label1 , signal_label2 ; 

  bool all_by_all = false;
  
  const bool mirror_outputs = param.has( "mirror" );

  if ( param.has( "sig1" ) )
    {
      signal_label1 =  param.requires( "sig1" );
      signal_label2 =  param.requires( "sig2" );
    }
  else
    {
      all_by_all = true;
      signal_label1 = signal_label2 = param.requires( "sig" );
    }
  
  
  signal_list_t signals1 = edf.header.signal_list( signal_label1 );  

  signal_list_t signals2 = edf.header.signal_list( signal_label2 );  

  // accumulate these below
  int ns1 = 0 , ns2 = 0;

  for (int s=0;s<signals1.size();s++)
    if ( ! edf.header.is_annotation_channel( signals1(s) ) ) ++ns1;
  
  for (int s=0;s<signals2.size();s++)
    if ( ! edf.header.is_annotation_channel( signals2(s) ) ) ++ns2;

  //
  // Epochs or whole signal?  Other output param
  //

  bool epoched = param.has("epoch") || param.has( "epoch-spectrum" ) ;

  bool show_spectrum = param.has( "spectrum" ) || param.has( "epoch-spectrum" );
  
  bool show_epoch_spectrum = param.has( "epoch-spectrum" );
  
  double upper_freq = param.has( "max" ) ? param.requires_dbl( "max" ) : 20 ; 
  double lower_freq = param.has( "min" ) ? param.requires_dbl( "min" ) : -1  ; 

  int sr = param.has( "sr" ) ? param.requires_int( "sr" ) : 0 ;

  //
  // FFT parameters
  //

  const double segment_sec = param.has( "segment-sec" ) ? param.requires_dbl( "segment-sec" ) : 4; // was 5 

  const double overlap_sec = param.has( "segment-inc" ) ? param.requires_dbl( "segment-inc" ) : 2; // was 0 

  window_function_t window = WINDOW_TUKEY50;
  if ( param.has( "hann" ) ) window = WINDOW_HANN;

  const bool average_adj = false;

  const bool detrend = false;
     

  //
  // Get all implied signals
  //
  
  std::vector<int> sigs;
  std::map<int,std::string> sigset;

  for (int s=0;s<ns1;s++) 
    {
      if ( sigset.find( signals1(s) ) == sigset.end() )
	{
	  if ( edf.header.is_annotation_channel( signals1(s) ) ) continue;
	  sigset[ signals1(s) ] =  signals1.label(s);
	  sigs.push_back( signals1(s) );
	}
    }

  
  for (int s=0;s<ns2;s++) 
    {
      if ( sigset.find( signals2(s) ) == sigset.end() )
	{
	  if ( edf.header.is_annotation_channel( signals2(s) ) ) continue;
	  sigset[ signals2(s) ] =  signals2.label(s);
	  sigs.push_back( signals2(s) );
	}
    }

  const int ns = sigs.size();


  //
  // check/adjust all SRs now
  //

  if ( sr )
    {
      for (int s=0;s<ns;s++)
	{
	  
	  if ( edf.header.is_annotation_channel( sigs[s] ) ) continue;
	  
	  if ( edf.header.sampling_freq( sigs[s] ) != sr ) 
	    {
	      logger << "  resampling channel " << sigset[ sigs[s] ] 
		     << " from " << edf.header.sampling_freq( sigs[s] )
		     << " to " << sr << "\n";
	      resample_channel( edf, sigs[s] , sr );
	    }
	}
    }
  else
    {
      int sr_check = 0;
      for (int s=0;s<ns;s++)
	{
	  if ( edf.header.is_annotation_channel( sigs[s] ) ) continue;
	  if ( sr_check == 0 ) sr_check = edf.header.sampling_freq( sigs[s] );
	  else if ( sr_check != edf.header.sampling_freq( sigs[s] ) )
	    Helper::halt( "all signals must have similar sampling rates (add 'sr' option)" );
	}
    }

  //
  // Nothing to do
  //

  if ( sigs.size() == 0 ) return;
  
  // track what SR is 
  if ( ! sr ) sr = edf.header.sampling_freq( sigs[0] );
      
  
  
  //
  // Number of pairwise comparisons
  //

  const int np = all_by_all ? (ns*(ns-1))/2 : ns1 * ns2;
  


  //
  // Initiate main coherence_t object  
  //

  const int Fs = sr;

  interval_t interval;

  if ( epoched )
    {      
      edf.timeline.first_epoch();
      int epoch = edf.timeline.next_epoch();
      if ( epoch == -1 ) Helper::halt( "no epochs to analyse" );
      interval = edf.timeline.epoch( epoch );      
    }
  else
    {
      interval = edf.timeline.wholetrace();
    }
  
  slice_t slice1( edf , sigs[0] , interval );      
  const std::vector<double> * d1 = slice1.pdata();
  const int total_sample_points = d1->size();  

  //
  // nb. this assumes epochs have similar sizes (this will have been checked in eval.cpp before calling)
  //

  coherence_t coherence( total_sample_points, Fs, segment_sec, overlap_sec, window , average_adj , detrend );


  //
  // Iterate over epochs (potentially) 
  //
  
  int epoch0 = epoched ? edf.timeline.first_epoch() : -9 ; 

  logger << "  calculating coherence for " << np << " channel pairs; ";

  if ( epoched ) 
    logger << "  iterating over epochs\n";
  else
    logger << "  for entire (unmasked) interval\n";
  
  //
  // Track epoch level results here; alternatively, drop
  // whole-trace results here for unitary output
  //
  
  std::map<int,std::map<int,coh_t> > coh;
  
  while ( 1 ) 
    {

      interval_t interval;

      int epoch = -9;
      
      // either a single epoch, or the entire signal

      if ( epoched )
	{
	  epoch = edf.timeline.next_epoch();      

	  if ( epoch == -1 ) break;
	  
	  interval = edf.timeline.epoch( epoch );
	}
      else
	interval = edf.timeline.wholetrace();


      //
      // Clear cache storing univariate spectra
      //
      
      coherence.clear();

      //
      // (Re)populate cache with single-channel spectra
      //

      for (int i=0;i<ns;i++)
	{
	  if ( edf.header.is_annotation_channel( sigs[i] ) ) continue;
	  dsptools::coherence_prepare( edf , sigs[i] , interval , &coherence );
	}  

      //
      // Iterate over pairs of channels
      //
      
      for (int i=0;i<ns1;i++)
	{
	  
	  if ( edf.header.is_annotation_channel( signals1(i) ) ) continue;
	  
	  writer.level( signals1.label(i) , "CH1" );
	  
	  //
	  // Second channel
	  //

	  for (int j=0;j<ns2;j++)
	    {
	      

	      if ( edf.header.is_annotation_channel( signals2(j) ) ) continue;
	      
	      // if all-by-all, do not do redundant channels	      
	      if ( all_by_all &&  j <= i ) continue;

	      // good to evaluate this channel pair
	      writer.level( signals2.label(j) , "CH2" );
	  
	      //
	      // Calculate coherence metrics
	      //
	     
	      scoh_t scoh = dsptools::coherence_do( &coherence , signals1(i) , signals2(j) );
	      

	      //
	      // Track Sxx, Syy, Sxy across epochs, for this pair (i,j)
	      //
	      
	      coh[i][j].add( scoh );
	      

	      //
	      // Output? (+/- full spectra)
	      //
	      
	      if ( epoched )
		{
		  writer.epoch( edf.timeline.display_epoch( epoch ) );
		  
		  scoh.proc_and_output( coherence , true, 
					show_epoch_spectrum ? upper_freq : -1 , 
					lower_freq );
		  
		}
	      
	      
	    } // second channel, j

	} // first channel, i


      //
      // If we were not iterating over epochs, now all done
      //
      
      if ( ! epoched ) break;


    } // next epoch

  
  if ( epoched )
    writer.unepoch();
  

  //
  // Output whole-signal coherence (this may optionally contain an average over all coh_t
  // if we previously analyzed epochs)
  //
  
  logger << "  calculating overall coherence statistics\n";

  std::map<int,std::map<int,coh_t> >::iterator ii = coh.begin();
  while ( ii != coh.end() )
    {
      
      writer.level( signals1.label(ii->first) , "CH1" );
      
      std::map<int,coh_t>::iterator jj = ii->second.begin();
      while ( jj != ii->second.end() )
	{
	  writer.level( signals2.label(jj->first) , "CH2" );
	  
	  // output (and calc stats averaged over epochs)	
	  jj->second.calc_stats( coherence , show_spectrum ? upper_freq : -1 , lower_freq );
	  
	  ++jj;
	}

      ++ii;
    }

  //
  // all done
  //

  writer.unlevel( "CH1" );  
  writer.unlevel( "CH2" );
  
}



void dsptools::coherence_prepare( edf_t & edf , const int signal1 , const interval_t & interval , coherence_t * ptr_coherence )
{
  slice_t slice1( edf , signal1 , interval );      
  const std::vector<double> * d1 = slice1.pdata();
  ptr_coherence->prepare( signal1 , *d1 );
}

scoh_t dsptools::coherence_do( coherence_t * coherence , const int signal1 , const int signal2 )
{  

  //
  // compute coherence 
  //
  
  coherence->process( signal1 , signal2 );
  
  return coherence->res;
    
}





void coh_t::calc_stats( const coherence_t & coherence , const double upper_freq , const double lower_freq ) 
{
  
  
  const int ne = epochs.size();

  // nothing to do
  if ( ne == 0 ) return;


  // if only a single epoch, then just use scoh_t to output
  if ( ne == 1 )
    {
      epochs[0].proc_and_output( coherence , true , upper_freq , lower_freq );
      return;
   }


  //
  // otherwise, average over epochs and calc all overall stats
  //
  
  const std::vector<double> & f = coherence.frq();

  const int nf = f.size();
  
  scoh_t scoh ( nf );
  
  for (int i=0;i<nf;i++)
    {
      
      double s = 0;
      
      std::complex<double> plv0;

      for (int e=0;e<ne;e++)
	{
	  
	  scoh.sxx[i] += epochs[e].sxx[i];
	  scoh.syy[i] += epochs[e].syy[i];
	  scoh.sxy[i] += epochs[e].sxy[i];
	  
	  
	  double Im = std::imag( scoh.sxy[i] );
	  if ( Im < 0 ) --s; else if ( Im > 0 )  ++s;

	  //std::cout << " s = " << Im << "\t" << s << " " << ne << "\n";
	  plv0 += scoh.sxy[i] / std::abs( scoh.sxy[i] );
	  //std::cout << "plv0 " << plv0 << "\n";
	}
      
      scoh.sxx[i] /= (double)ne;
      scoh.syy[i] /= (double)ne;
      scoh.sxy[i] /= (double)ne;
      
      plv0 /= (double)ne;

      s /= (double)ne;
      
      
      //
      // connectivity statistics
      //

      double & Sxx = scoh.sxx[i];
      double & Syy = scoh.syy[i];

      std::complex<double> & Sxy = scoh.sxy[i];

      double Re = std::real( Sxy );

      double Im = std::imag( Sxy );
      
      double phi = abs( Sxy ) ;      

      double phi2 = phi * phi;

      //
      // Phase Slope Index (PSI)
      // http://doc.ml.tu-berlin.de/causality/
      //
            

      //
      // Final stats
      //

      // plv : phase-locking value      
      // pli : phase lag index
      // wpli : weight phase lag index
      

      // coh = phi2 / ( Sxx * Syy ); 
      
      // icoh = Im / sqrt( Sxx * Syy );

      // lcoh = Im / sqrt( Sxx * Syy - ( Re * Re ) );

      // plv = abs( plv0 );

      // pli = fabs( s );
            
      
     
    }
  

  //
  // Output
  //
  
  scoh.proc_and_output( coherence , true, upper_freq , lower_freq );
  
  
}






void scoh_t::proc_and_output( const coherence_t & coherence , 
			      const bool output , 
			      const double upper_freq ,  // -1 means no epoch output
			      const double lower_freq  )
{
  
  //
  // Band-level summaries
  //

  bcoh.clear();
  bicoh.clear();
  blcoh.clear();
  bn.clear();

  std::vector<frequency_band_t> bands;
  bands.push_back( SLOW );
  bands.push_back( DELTA );
  bands.push_back( THETA );
  bands.push_back( ALPHA );
  bands.push_back( SIGMA );
  bands.push_back( LOW_SIGMA );
  bands.push_back( HIGH_SIGMA );
  bands.push_back( BETA );
  bands.push_back( GAMMA );

  // track output status
  bool any = false;

    
  //
  // iterate over frequencies
  //

  const std::vector<double> & frq = coherence.frq();

  for ( int k=0; k< frq.size() ; k++)
    {
      
      if ( bad[k] ) continue;
      
      //
      // main stats
      //
      
      const double & Sxx = sxx[k];
      const double & Syy = syy[k];
      const std::complex<double> & Sxy = sxy[k];
      
      double Re = std::real( Sxy );
      double Im = std::imag( Sxy );
      double phi = abs( Sxy ) ;
      double phi2 = phi * phi;
      double coh = phi2 / ( Sxx * Syy );
      double icoh = Im / sqrt( Sxx * Syy );
      double lcoh = Im / sqrt( Sxx * Syy - ( Re * Re ) );

      double cross_spectra_dB = 5.0 * log10( phi2 );
      
      
      // if ( frq[k] < 20 ) 
      // 	std::cout << "dets = " << frq[k] << "\t" << Sxx << " " << Syy << " " << " " << Re << " " << Im << " " << coh << " " << icoh << "\n";
      

      //
      // band-level summaries
      //
      
      std::vector<frequency_band_t>::const_iterator bb = bands.begin();
      while ( bb != bands.end() )
	{

	  if ( frq[k] >= globals::freq_band[ *bb ].first && frq[k] < globals::freq_band[ *bb ].second ) 
	    {
	      
	      if ( Helper::realnum( coh ) && Helper::realnum( icoh ) && Helper::realnum( lcoh ) ) 
		{
		  bcoh[ *bb ] += coh;
		  bicoh[ *bb ] += icoh;
		  blcoh[ *bb ] += lcoh;
		  ++bn[ *bb ];	      
		}
	    }
	  
	  ++bb;
	}

      
      //
      // frequency-bin output (uf = -1 setting turns off output , used for epoch-output)
      //

      if ( output ) 
	if ( frq[k] >= lower_freq && upper_freq > 0 && frq[k] <= upper_freq ) 
	  {
	    
	    any = true;
	    
	    writer.level( frq[k] , globals::freq_strat );
	    
	    if ( Helper::realnum( coh ) )
	      writer.value( "COH" , coh );
	    if ( Helper::realnum( icoh ) )
	      writer.value( "ICOH" , icoh );
	    if ( Helper::realnum( lcoh ) )
	      writer.value( "LCOH" , lcoh );
	    if ( Helper::realnum( cross_spectra_dB ) )
	      writer.value( "CSPEC" , cross_spectra_dB );
	  }      
      
    } // next frequency 'k'

  
  //
  // clear any frequency-strata, if set
  //
  
  if ( output && any )
    writer.unlevel( globals::freq_strat );
  
  
  //
  // band summary/output
  //

  any = false;

  std::vector<frequency_band_t>::const_iterator bb = bands.begin();
  while ( bb != bands.end() )
    {
      
      //std::cerr << " b " << globals::band( *bb ) << " " << bn[ *bb ] << "\n";
	    
      if ( bn[ *bb ] )
	{

	  bcoh[ *bb ] /= (double)bn[ *bb ];
	  bicoh[ *bb ] /= (double)bn[ *bb ];
	  blcoh[ *bb ] /= (double)bn[ *bb ];

	  if ( output ) 
	    {
	      any = true;
	      writer.level( globals::band( *bb ) , globals::band_strat );	  
	      
	      if ( Helper::realnum( bcoh[ *bb ] ) )
		writer.value( "COH" , bcoh[ *bb ] );
	      
	      if ( Helper::realnum( bicoh[ *bb ] ) )
		writer.value( "ICOH" , bicoh[ *bb ] );
	      
	      if ( Helper::realnum( blcoh[ *bb ] ) )
		writer.value( "LCOH" , blcoh[ *bb ] );	  
	    }      
	}

      ++bb;
    }
  
  if ( output && any )
    writer.unlevel( globals::band_strat );
  
}


