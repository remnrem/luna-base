
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

#include "dsp/wrappers.h"
#include "edf/edf.h"
#include "eval.h"
#include "edf/slice.h"
#include "cwt/cwt.h"
#include "dsp/hilbert.h"

#include <vector>
  
void dsptools::cwt( edf_t & edf , param_t & param )
{

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  
  const int ns = signals.size();

  double fc = param.requires_dbl( "fc" );

  bool alt_spec = param.has( "fwhm" );

  double fwhm = alt_spec ? param.requires_dbl( "fwhm" ) : 0 ;
  
  int num_cycles = alt_spec ? 0 : param.requires_int( "cycles" );

  double timelength = alt_spec ? ( param.has( "len" ) ? param.requires_dbl( "len" ) : 20 ) : 0 ; 
  
  bool return_phase = param.has( "phase" );
  
  bool wrapped_wavelet = param.has( "wrapped" );

  std::string tag = param.has( "tag" ) ? "_" + param.value( "tag" ) : "" ; 

  for (int s=0;s<ns;s++)
    {

      if ( edf.header.is_annotation_channel( signals(s) ) ) 
	continue;
      
      const int Fs = edf.header.sampling_freq( signals(s) );
      
      interval_t interval = edf.timeline.wholetrace();
      
      slice_t slice( edf , signals(s) , interval );

      const std::vector<double> * d = slice.pdata();
      
      std::vector<double> mag , phase;

      if ( alt_spec )
	alt_run_cwt( *d , Fs , fc , fwhm , timelength , wrapped_wavelet , &mag , return_phase ? &phase : NULL );
      else
	run_cwt( *d , Fs , fc , num_cycles , &mag , return_phase ? &phase : NULL );
      
      std::string new_mag_label = signals.label(s) + tag + "_cwt_" + Helper::dbl2str(fc) + "_" + Helper::int2str( num_cycles ) + "_mag";
      std::string new_phase_label = signals.label(s) + tag + "_cwt_" + Helper::dbl2str(fc) + "_" + Helper::int2str( num_cycles ) + "_phase";
      if ( alt_spec )
	{
	  new_mag_label = signals.label(s) + tag + "_cwt_" + Helper::dbl2str(fc) + "_fwhm" + Helper::dbl2str( fwhm ) + "_mag";
	  new_phase_label = signals.label(s) + tag + "_cwt_" + Helper::dbl2str(fc) + "_fwhm" + Helper::dbl2str( fwhm ) + "_phase";
	}

      logger << " CWT for " << signals.label(s) << " --> " << new_mag_label ;      
      if ( return_phase ) 
	logger << ", " << new_phase_label ;
      logger << "\n";

      edf.add_signal( new_mag_label , Fs , mag );
      
      if ( return_phase ) 
	edf.add_signal( new_phase_label , Fs , phase );
      
    }  
  
}


void dsptools::hilbert( edf_t & edf , param_t & param )
{

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  
  const int ns = signals.size();
  
  std::vector<double> frqs = param.dblvector( "f" );

  bool no_filter = frqs.size() != 2 ;

  if ( no_filter )
    {
      // this is checked in run_hilbert() for impossible freqs
      frqs.resize( 2 , -1 );
    }

  double ripple = 0, tw = 0;;

  if ( ! no_filter )
    {
      ripple = param.requires_dbl( "ripple" );
      tw = param.requires_dbl( "tw" );
    }
  
  if ( no_filter && param.has( "ripple" ) ) Helper::halt( "cannot specify ripple with no filter" );
  if ( no_filter && param.has( "tw" ) ) Helper::halt( "cannot specify tw with no filter" );
    
  bool return_phase = param.has( "phase" ) || param.has( "angle" );

  bool return_angle = param.has( "angle" );

  bool return_ifrq = param.has( "ifrq" );
  
  std::string tag = param.has( "tag" ) ? "_" + param.value( "tag" ) : "" ; 

  for (int s=0;s<ns;s++)
    {

      if ( edf.header.is_annotation_channel( signals(s) ) ) 
	continue;
      
      const int Fs = edf.header.sampling_freq( signals(s) );
      
      interval_t interval = edf.timeline.wholetrace();
      
      slice_t slice( edf , signals(s) , interval );

      const std::vector<double> * d = slice.pdata();
      
      std::vector<double> mag , phase, ifrq;

      run_hilbert( *d , Fs , frqs[0] , frqs[1] , ripple , tw , &mag ,
		   return_phase ? &phase : NULL ,
		   return_angle ? &phase : NULL ,
		   return_ifrq ? &ifrq : NULL );
            
      std::string new_mag_label = signals.label(s) + tag + "_hilbert_"   + Helper::dbl2str(frqs[0]) + "_" + Helper::dbl2str( frqs[1] ) + "_mag";
      std::string new_phase_label = signals.label(s) + tag + "_hilbert_" + Helper::dbl2str(frqs[0]) + "_" + Helper::dbl2str( frqs[1] ) + "_phase";
      std::string new_angle_label = signals.label(s) + tag + "_hilbert_" + Helper::dbl2str(frqs[0]) + "_" + Helper::dbl2str( frqs[1] ) + "_angle";
      std::string new_ifrq_label = signals.label(s) + tag + "_hilbert_"  + Helper::dbl2str(frqs[0]) + "_" + Helper::dbl2str( frqs[1] ) + "_ifrq";

      if ( no_filter )
	{
	  new_mag_label = signals.label(s) + tag + "_hilbert_mag";
	  new_phase_label = signals.label(s) + tag + "_hilbert_phase";
	  new_angle_label = signals.label(s) + tag + "_hilbert_angle";
	  new_ifrq_label = signals.label(s) + tag + "_hilbert_ifrq";
	}

      logger << " Hilbert transform for " << signals.label(s) << " --> " << new_mag_label ;      

      if ( return_phase ) 
	logger << ", " << new_phase_label ;

      if ( return_angle )
	logger << ", " << new_angle_label ;

     if ( return_ifrq ) 
	logger << ", " << new_ifrq_label ;

      logger << "\n";

      edf.add_signal( new_mag_label , Fs , mag );
      
      if ( return_phase ) 
	edf.add_signal( new_phase_label , Fs , phase );
      
      if ( return_angle )
	edf.add_signal( new_angle_label , Fs , phase );

      if ( return_ifrq ) 
	{
	  // this returns n-1 estimates, i.e. based on the derivative of the phase
	  // add 0 to the end as a null marker, so it can be placed back in the EDF
	  ifrq.push_back(0);
	  edf.add_signal( new_ifrq_label , Fs , ifrq );
	}
    }  

}

  
void dsptools::alt_run_cwt( const std::vector<double> & data ,
			    const int Fs, 
			    const double fc ,
			    const double FWHM ,
			    const double tlen , 
			    const bool wrapped , 
			    std::vector<double> * mag , 
			    std::vector<double> * phase )
{
  
  CWT cwt;
  
  cwt.set_sampling_rate( Fs );
  
  cwt.set_timeframe( 50.0 / tlen );

  cwt.alt_add_wavelet( fc , FWHM , tlen );
  
  cwt.store_real_imag_vectors( true );

  cwt.load( &data );

  if ( wrapped ) 
    cwt.run_wrapped();
  else
    cwt.run();
  
  *mag = cwt.results(0);
  
  if ( phase != NULL ) 
    *phase = cwt.phase(0);
  
}



void dsptools::run_cwt( const std::vector<double> & data ,
			const int Fs, 
			const double fc ,
			const int num_cycles , 
			std::vector<double> * mag , 
			std::vector<double> * phase )
{
  
  CWT cwt;
  
  cwt.set_sampling_rate( Fs );
  
  cwt.add_wavelet( fc , num_cycles ); 
  
  cwt.load( &data );
  
  cwt.run();
  
  *mag = cwt.results(0);
  
  if ( phase != NULL ) 
    *phase = cwt.phase(0);
  
}


void dsptools::run_hilbert( const std::vector<double> & data , const int Fs , 
			    const double flwr , const double fupr , const double ripple , const double tw , 
			    std::vector<double> * mag , 
			    std::vector<double> * phase ,
			    std::vector<double> * angle ,
			    std::vector<double> * ifrq )
{
  
  if ( flwr < 0 )
    {
      // straight Hilbert , no filter
      hilbert_t hilbert( data );

      *mag = *(hilbert.magnitude());

      if ( phase != NULL ) *phase = *(hilbert.phase());
      
      if ( angle != NULL )
	{
	  *angle = *phase;
	  // convert to degrees with 0 as pos-to-neg crossing
	  for (int i=0;i<angle->size();i++) (*angle)[i] = MiscMath::as_angle_0_pos2neg( (*angle)[i] );	  
	}
      
      if ( ifrq != NULL )  *ifrq = hilbert.instantaneous_frequency( Fs );    

    }
  else
    {

      // filter-Hilbert 
      hilbert_t hilbert( data , Fs , flwr , fupr , tw , ripple );
      *mag = *(hilbert.magnitude());
      if ( phase != NULL ) *phase = *(hilbert.phase());

      if ( angle != NULL )
	{
	  *angle = *phase;
	  // convert to degrees with 0 as pos-to-neg crossing
	  for (int i=0;i<angle->size();i++) (*angle)[i] = MiscMath::as_angle_0_pos2neg( (*angle)[i] );	  
	}
      
      if ( ifrq != NULL )  *ifrq = hilbert.instantaneous_frequency( Fs );    
    }

}


