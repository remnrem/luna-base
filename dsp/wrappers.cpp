
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
  
  int num_cycles = param.requires_int( "cycles" );
  
  bool return_phase = param.has( "phase" );
  
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
      
      run_cwt( *d , Fs , fc , num_cycles , &mag , return_phase ? &phase : NULL );
      
      std::string new_mag_label = signals.label(s) + tag + "_cwt_" + Helper::dbl2str(fc) + "_" + Helper::int2str( num_cycles ) + "_mag";
      std::string new_phase_label = signals.label(s) + tag + "_cwt_" + Helper::dbl2str(fc) + "_" + Helper::int2str( num_cycles ) + "_phase";
      
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

  if ( frqs.size() != 2 ) Helper::halt( "expecting f=lower,upper" );

  double ripple = param.requires_dbl( "ripple" );
  
  double tw = param.requires_dbl( "tw" );
  
  bool return_phase = param.has( "phase" );

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
      
      run_hilbert( *d , Fs , frqs[0] , frqs[1] , ripple , tw , &mag , return_phase ? &phase : NULL , return_ifrq ? &ifrq : NULL );
            
      std::string new_mag_label = signals.label(s) + tag + "_hilbert_"   + Helper::dbl2str(frqs[0]) + "_" + Helper::dbl2str( frqs[1] ) + "_mag";
      std::string new_phase_label = signals.label(s) + tag + "_hilbert_" + Helper::dbl2str(frqs[0]) + "_" + Helper::dbl2str( frqs[1] ) + "_phase";
      std::string new_ifrq_label = signals.label(s) + tag + "_hilbert_"  + Helper::dbl2str(frqs[0]) + "_" + Helper::dbl2str( frqs[1] ) + "_ifrq";
      
      logger << " Hilbert transform for " << signals.label(s) << " --> " << new_mag_label ;      

      if ( return_phase ) 
	logger << ", " << new_phase_label ;

      if ( return_ifrq ) 
	logger << ", " << new_ifrq_label ;

      logger << "\n";

      edf.add_signal( new_mag_label , Fs , mag );
      
      if ( return_phase ) 
	edf.add_signal( new_phase_label , Fs , phase );
      
      if ( return_ifrq ) 
	{
	  // this returns n-1 estimates, i.e. based on the derivative of the phase
	  // add 0 to the end as a null marker, so it can be placed back in the EDF
	  ifrq.push_back(0);
	  edf.add_signal( new_ifrq_label , Fs , ifrq );
	}
    }  

}

  
void dsptools::run_cwt( const std::vector<double> & data , const int Fs, 
			const double fc , const int num_cycles , 
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
			    std::vector<double> * ifrq )
{

  hilbert_t hilbert( data , Fs , flwr , fupr , tw , ripple );
  
  *mag = *(hilbert.magnitude());

  if ( phase != NULL )
    *phase = *(hilbert.phase());

  if ( ifrq != NULL )   
      *ifrq = hilbert.instantaneous_frequency( Fs );    
}


