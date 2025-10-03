
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

#include "ed.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "stats/matrix.h"
#include "stats/statistics.h"

#include "param.h"
#include "db/db.h"

extern writer_t writer;

// Alschuler et al (2014) Identifying electrode bridging from
// electrical distance distributions: a survey of publicly-available
// EEG data using a new method. Clin Neurophysiol. 2014 Mar; 125(3):
// 484–490.

// http://psychophysiology.cpmc.columbia.edu/mmedia/SPR2000/iHjorth.pdf

void dsptools::elec_distance( edf_t & edf , param_t & param )
{
  
  std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  // drop annot channels
  edf.header.drop_annots_from_signal_list( &signals );

  const int ns = signals.size();

  bool epoch = param.has( "epoch" );
  
  // check all same SR
  int sr = 0;
  for (int s=0;s<ns;s++)
    {
      
      if ( sr == 0 ) sr = edf.header.sampling_freq( signals(s) ) ;
      if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	Helper::halt( "requires all signals to have similar sampling rate" );      
    }
  
  
  //
  // Step through each epoch/channel
  //
  
  int ne = edf.timeline.first_epoch();
  
  std::vector<std::vector<std::vector<double> > > ED( ne );

  std::vector<double> EDmedian;

  while ( 1 ) 
     {

       int epoch = edf.timeline.next_epoch();      

       if ( epoch == -1 ) break;

       interval_t interval = edf.timeline.epoch( epoch );
       
       mslice_t mslice( edf , signals , interval );

       Data::Matrix<double> D = mslice.extract();
       
       //
       // consider each channel pair
       //

       std::vector<std::vector<double> > & eED = ED[ epoch ];
       
       eED.resize( ns );
       
       for (int s1=0;s1<ns;s1++)
	 {
	   
	   for (int s2=s1;s2<ns;s2++)
	     {
	       if ( s1==s2 ) continue;
	       
	       const int nr = D.dim1();
	       Data::Vector<double> ed( nr );
	       for (int i=0;i<nr;i++) ed[i] = D(i,s1) - D(i,s2);
	       double _ED = Statistics::variance( ed );
	       eED[s1].push_back( _ED );
	       EDmedian.push_back( _ED );
	     }
	 }
     }




  //
  // Scaling factor = 100/median
  //

  double fac = 100 / MiscMath::median( EDmedian );
  
  EDmedian.clear();
  
  for (int e=0;e<ne;e++)
    for (int s1=0;s1<ns;s1++)
      for (int s2=s1+1;s2<ns;s2++)
	ED[e][s1][s2] *= fac ; 

  
  //
  //
  //



}
