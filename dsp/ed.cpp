
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
#include "../edf/edf.h"
#include "../main.h"
#include "../db/db.h"

extern writer_t writer;

// Alschuler et al (2014) Identifying electrode bridging from
// electrical distance distributions: a survey of publicly-available
// EEG data using a new method. Clin Neurophysiol. 2014 Mar; 125(3):
// 484–490.

// http://psychophysiology.cpmc.columbia.edu/mmedia/SPR2000/iHjorth.pdf

void dsptools::elec_distance( edf_t & edf , param_t & param )
{
  
  std::string signal_label = param.requires( "signal" );
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
  
//   int ne = edf.timeline.first_epoch();
  
//   std::vector<std::vector<std::vector<double> > > ED( ns );
  
//   while ( 1 ) 
//     {
      
//       int epoch = edf.timeline.next_epoch();      
       
//       if ( epoch == -1 ) break;
      
//       interval_t interval = edf.timeline.epoch( epoch );
      
//       mslice_t mslice( edf , signals , interval );
      
//       Data::Matrix<double> D = mslice.extract();
 
//       Data::Matrix<double> ED( ns , ns );

//       //
//       // consider each channel pair
//       //
      
//       for (int s1=0;s1<ns;s1++)
// 	for (int s2=s1;s2<ns;s2++)
// 	  {
// 	    if ( s1==s2 ) continue;
	    
// 	    const int nr = D.dim1();
// 	    Data::Vector<double> ed( nr );
// 	    for (int i=0;i<nr;i++) ed[i] = D(i,s1) - D(i,s2);
// 	    double _ED = Statistics::variance( ed );
// 	    ED(s1,s2) = _ED;
// 	  }

//       //
//       // Next epoch, so store the ch x ch matrix
//       //
      
//       epochED.push_back( ED );
      
//     }

  
//   //
//   // Scale by the median
//   //

//   Data::Matrix<double> m(ns,ns);

//   const int ne = epochED.size();
  
//   for (int s1=0;s1<ns;s1++)
//     for (int s2=s1;s2<ns;s2++)
//       {
// 	if ( s1==s2 ) continue;
	
// 	std::vector<double> x;
// 	for (int 
	  
//       writer.unlevel( globals::signal_strat );
      
//     }
//   writer.unepoch();
  



//       writer.epoch( edf.timeline.display_epoch( epoch ) );

// 	    // write
// 	    // writer.level( signals.label(s) , globals::signal_strat );
// 	    // 	  writer.value( "R" , r ); 
//       writer.unlevel( globals::signal_strat );
      
//     }
//   writer.unepoch();

  
// }

}
