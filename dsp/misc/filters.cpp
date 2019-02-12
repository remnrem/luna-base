
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


#include "filters.h"

#include "../edf/edf.h"
#include "../main.h"


std::vector<double> band_pass_filter( const std::vector<double> & input , 
				      double Fs , 
				      int num_taps , 
				      double lwr , 
				      double upr )
{  
  MyFilter f;
  f.band_pass( num_taps , Fs , lwr , upr );            
  return f.filter( input );
}


void band_pass_filter( edf_t & edf , 
		       param_t & param )
{
  
  //
  // Signals
  //

  std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );      

  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  bool verbose = param.has( "verbose" );
  
  const int ns = signals.size();

  //
  // Filter type
  // 

  int num_taps = -1;
  if ( param.has( "num_taps" ) ) num_taps = param.requires_int( "num_taps" );
  
  double lwr = param.requires_dbl( "lower" );
  double upr = param.requires_dbl( "upper" );

  //
  // Filter each 
  //

  interval_t interval = edf.timeline.wholetrace();

  for (int s=0; s<ns; s++)
    {

      //
      // skip annotation channels
      //

      if ( edf.header.is_annotation_channel(s) ) continue;

      //
      // Filter data
      //
      
      std::cerr << " filtering channel " << edf.header.label[ signals(s) ] << "\n";

      //
      // Pull entire signals out
      //
      
      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();
      
      
      //
      // Filter order
      //
      
      int filter_order = num_taps == -1 ? 3 * ( Fs[s] / lwr ) : num_taps ; 
      
      //      std::cerr << " filter order is " << filter_order << " " << Fs[s] << " " << lwr << "\n";
	
      //
      // Filter each signal independently
      //
      
      MyFilter f;
      f.band_pass( filter_order , Fs[s] , lwr , upr );            
      std::vector<double> filtered = f.filter( *d );
      
      //
      // Place back
      //

      edf.update_signal( signals(s) , &filtered );
      
    } // next signal
}
