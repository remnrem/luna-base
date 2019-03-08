
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


#include "staging.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "eval.h"
#include "miscmath/miscmath.h"
#include "fftw/fftwrap.h"
#include "db/db.h"

#include <iostream>
#include <string>

#include "helper/helper.h"
#include "helper/logger.h"

extern writer_t writer;

extern logger_t logger;

void zratio_t::calc( edf_t & edf , const std::string & signal_label )
{

  //
  // Attach signals
  //

  signal_list_t signals = edf.header.signal_list( signal_label );  
  std::vector<double> Fs = edf.header.sampling_freq( signals );

  const int ns = signals.size();

  //
  // FFT parameters
  //

  double fft_segment_size    = 2;
  double fft_segment_overlap = 0; 
  
  freq_range_t zr_delta( 0.5 , 2.0 );
  freq_range_t zr_theta( 2.5 , 7.5 );
  freq_range_t zr_alpha( 8.0 , 12.5 );
  freq_range_t zr_beta( 13.0 , 30 );
  
  //
  // Assume is as 30-second epochs
  //
  
  if ( ! edf.timeline.epoched() ) Helper::halt( "require epoched data" );
  int ne30 = edf.timeline.num_total_epochs();
    
  //
  // Set 2-second epochs
  //
  
  double saved_epoch_length = edf.timeline.epoch_length(); 
  double saved_epoch_inc = edf.timeline.epoch_inc(); 
  
  if ( fabs( saved_epoch_length - 30.0 ) > 0.001 ) 
    Helper::halt( "require 30-second epochs initially" );
  
  int ne = edf.timeline.set_epoch( 2 , 2 ); 
  

  //
  // Loop over signals
  //

  for (int s = 0 ; s < ns ; s++ )
    {

      // only consider data tracks
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      // stratify output by signal
      
      writer.level( signals.label(s) , globals::signal_strat );
      
      // Set first epoch
      
      edf.timeline.first_epoch();

      zr2.clear();
      zr30.clear();


      //
      // for each each epoch 
      //

      while ( 1 ) 
	{
	  
	  int epoch = edf.timeline.next_epoch();      
	  
	  if ( epoch == -1 ) break;
	  
	  writer.level( epoch , "E2" );
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  const std::vector<double> * d = slice.pdata();
	  
	  const int total_points = d->size();
	  
	  FFT fft( total_points , Fs[s] , FFT_FORWARD , WINDOW_NONE );

	  fft.apply( *d );
	  
	  int N = fft.cutoff;
	  
	  double delta = 0 , theta = 0 , alpha = 0 , beta = 0 ;

	  for (int f=0;f<N;f++)
	    {
	      //std::cout << "zr\t det " << f << "\t" << fft.frq[f] << "\t" << fft.X[f] << "\n";
	      if ( fft.frq[f] >= zr_delta.first && fft.frq[f] < zr_delta.second ) delta += fft.X[f];
	      if ( fft.frq[f] >= zr_theta.first && fft.frq[f] < zr_theta.second ) theta += fft.X[f];
	      if ( fft.frq[f] >= zr_alpha.first && fft.frq[f] < zr_alpha.second ) alpha += fft.X[f];
	      if ( fft.frq[f] >= zr_beta.first  && fft.frq[f] < zr_beta.second ) beta += fft.X[f];
	    }
	  
// 	  delta /= zr_delta.second - zr_delta.first ;
// 	  theta /= zr_theta.second - zr_theta.first ;
// 	  alpha /= zr_alpha.second - zr_alpha.first ;
// 	  beta /= zr_beta.second - zr_beta.first ;
	  
// 	  delta = log( delta );
// 	  theta = log( theta );
// 	  alpha = log( alpha );
// 	  beta = log( beta );

	  // calculate z-ratio, catching flat/clipped signals
	  
	  double denom =  delta + theta + alpha + beta ;
	  
	  double zr = denom > 0 ? ( ( delta + theta ) - ( alpha + beta ) ) / ( delta + theta + alpha + beta ) : -9;

	  if ( zr >= -1 ) 
	    writer.value( "ZR2" , zr );
	  
	  zr2.push_back( zr );
	  
	} // next epoch
      
      
      writer.unlevel( "E2" );
      
      //
      // compile into 30-second level stats epoch 
      //

      int ec2 = 0;
      int e2max = zr2.size();

      for (int e = 0; e < ne30; e++ )
	{
	  int elast = ec2 + 15 >= e2max ? e2max - 1 : ec2 + 15 - 1 ;
	  int nbins = elast - ec2 + 1;
	  
	  // TODO... need to check for -9 bad values and skip those. 
	  
	  // Zpage : average of the majority
	  if ( false )
	    {
	      int pos = 0 , neg = 0;
	      double mpos = 0 , mneg = 0;
	      for (int e2 = ec2 ; e2 <= elast ; e2++ )
		{
		  if   ( zr2[ e2 ] <= 0 ) { mneg += zr2[ e2 ]; ++neg; }
		  else if ( zr2[ e2 ] > 0 ) { mpos += zr2[ e2 ]; ++pos; }	      
		}
	      
	      double acc_zr;
	      
	      if ( pos > neg ) acc_zr = mpos / (double)pos;
	      else acc_zr = mneg / (double)neg;
	    }

	  // simple average
	  double mval = 0 , mcnt = 0;
	  for (int e2 = ec2 ; e2 <= elast ; e2++ )
	    {
	      mval += zr2[ e2 ];
	      mcnt++;
	    }
	  
	  double acc_zr = mval / mcnt ; 
	  
	  std::cout << "Z " << nbins << "\t" 
		    << e << "\t" << ec2 << ".." << elast <<"\t" << acc_zr << "\n";
	  
	  ec2 += 15;
	  zr30.push_back( acc_zr );
	}     
  
      //
      // Output
      //

      for (int epoch = 0 ; epoch < zr30.size() ; epoch++)
	{

	  writer.level( epoch + 1 , "E30" );
	  
	  writer.value( "ZR30" , zr30[ epoch ] );
	  
	}
      
      writer.unlevel( "E30" );
      
      } // next signal
	
	writer.unlevel( globals::signal_strat );
   
  //
  // reset epochs
  //

  edf.timeline.set_epoch( saved_epoch_length , saved_epoch_inc ); 


}

 
