
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

#include "dsp/siggen.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include <vector>

#include <cmath>

void dsptools::siggen( edf_t & edf , param_t & param )
{

  //
  // Add or update a signal
  //
  
  const std::string siglab = param.requires( "sig" );
  if ( siglab == "*" ) Helper::halt( "must specify a single signal 'sig'" );
  if ( siglab.find( "," ) != std::string::npos ) Helper::halt( "must specify a single signal 'sig'" );

  // channel already exists?
  const bool update_existing_channel = edf.header.has_signal( siglab );
  const bool add_to_existing = param.has( "add" );

  // SR of a new channel?
  const int fs = update_existing_channel ? -1 : param.requires_int( "sr" );  
  if ( (!update_existing_channel) && fs < 1 )
    Helper::halt( "requires postive sample rate specified with 'sr'" );

  
  //
  // Options
  //

  // sine, square, saw, triangular

  bool sine_wave = param.has( "sine" );

  std::vector<double> sine_param;

  if ( sine_wave )
    {
      sine_param = param.dblvector( "sine" );
      if ( sine_param.size() == 2 ) sine_param.resize(3,0);
      else if ( sine_param.size() != 3 ) Helper::halt( "expecting sine=frq,amp{,phase}" );
      if ( sine_param[0] <= 0 ) Helper::halt( "frq must be positive" );
      if ( sine_param[0] >= fs / 2.0 ) Helper::halt( "frq not under Nyquist frequency, given sample rate" );
      if ( sine_param[1] <= 0 ) Helper::halt( "amp should be positive, non-zero" );
    }


  //
  // simple impulses
  //

  // impulse=T,A,D,T,A,D,...
  bool impulses = param.has( "impulse" ) ;
  std::vector<double> impulse_t, impulse_a;
  std::vector<int> impulse_d; // in samples
  std::vector<double> impulse_all;
  if ( impulses ) impulse_all = param.dblvector( "impulse" );
  if ( impulse_all.size() % 3 != 0 ) Helper::halt( "need impulse=T,A,D,T,A,D,..." );

  if ( impulses )
    {
      int n = impulse_all.size() ;
      for (int i=0; i<n; i+=3)
	{
	  impulse_t.push_back( impulse_all[i] );
	  impulse_a.push_back( impulse_all[i+1] );
	  impulse_d.push_back( impulse_all[i+2] );
	}
    }
  
  
  //
  // make synthetic signal
  //

  const int np = edf.header.record_duration * edf.header.nr * fs ;
    
  std::vector<double> d( np , 0 );
  

  //
  // sine waves
  //

  if ( sine_wave )
    {
      for ( int p=0 ; p<np; p++ )
	{
	  // time in seconds
	  double t = p /(double)fs;
	  
	  // add a sine wave?
	  if ( sine_wave )
	    d[p] += sine_param[1] * sin( 2 * M_PI * sine_param[0] * t + sine_param[2] );
	  // others to go here..

	}  
    }

  
  //
  // add impulse?
  //
  
  if ( impulses )
    {
      for (int i=0; i<impulse_t.size(); i++)
	{
	  int start = impulse_t[i] * np;
	  int end = start + impulse_d[i];
	  end = end >= np ? np-1 : end ;

	  for (int j=start; j<end; j++)
	    d[j] += impulse_a[i];
	}
    }
  
  //
  // Create/update signal
  //
  
  if ( update_existing_channel )
    {
      
      const int slot = edf.header.signal( siglab );
      
      if ( edf.header.is_annotation_channel( slot ) )
	Helper::halt( "cannot modify an EDF Annotation channel" );
      
      slice_t slice( edf , slot , edf.timeline.wholetrace() );
      
      const std::vector<double> * d1 = slice.pdata();
      
      if ( d1->size() != d.size() )
	Helper::halt( "internal error in siggen()" );
      
      if ( add_to_existing )
	{
	  // add existing to simulated signal
	  for (int i=0; i<d1->size(); i++)
	    d[i] += (*d1)[i];
	}
            
      // now update the channel       
      logger << "  updating " << siglab << "...\n";
      edf.update_signal( edf.header.signal( siglab ) , &d );

    }
  else
    {
      
      logger << "  creating new channel " << siglab << "...\n";
      edf.add_signal( siglab , fs , d );
    }
  
  //
  // all done
  //
  
  
}

