
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

#include "mtm.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "eval.h"
#include "fftw/fftwrap.h"

#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "nrutil.h"

extern writer_t writer;
extern logger_t logger; 

void mtm::wrapper( edf_t & edf , param_t & param )
{
  
  std::string signal_label = param.requires( "sig" );

  signal_list_t signals = edf.header.signal_list( signal_label );    

  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  const int ns = signals.size();

  // other parameters
  bool whole_signal_analysis = true;
  bool epoch_level_output = param.has( "epoch" );
  if ( param.has( "epoch-only" ) )
    {
      whole_signal_analysis = false;
      epoch_level_output = true;
    }

  // MTM parameters
  double npi = param.has( "nw" ) ? param.requires_dbl( "nw" ) : 3 ;
  int nwin = param.has( "t" ) ? param.requires_int( "t" ) : 2*npi-1 ;
  
  double max_f = param.has( "max" ) ?  param.requires_dbl( "max" ) : 20;  // default up to 20 Hz
  double wid_f = param.has( "bin" ) ? param.requires_dbl( "bin" ) : 0.5 ; // default 0.5 Hz bins


  // output
  
  bool dB = param.has( "dB" );
  
  if ( param.has( "full-spectrum" ) ) wid_f = 0 ; // return full spectrum, no binning.
  
  
  //
  // Whole signal analyses
  //

  if ( whole_signal_analysis )
    {

      interval_t interval = edf.timeline.wholetrace();
      
      //
      // Get each signal
      //
      
      for (int s = 0 ; s < ns; s++ )
	{
	  
	  //
	  // only consider data tracks
	  //
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) )
	    continue;
	  
	  //
	  // Stratify output by channel
	  //
	  
	  writer.level( signals.label(s) , globals::signal_strat );
	  
	  //
	  // Get data
	  //
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  const std::vector<double> * d = slice.pdata();	   
	  
	  //
	  // call MTM
	  //
	  
	  mtm_t mtm( npi , nwin );
	  
	  mtm.dB = dB;

	  mtm.apply( d , Fs[s] );
	  
	  // for 1-sided spectrum, scale by x2
	  // (skipping DC and NQ)
	  
	  if ( wid_f > 0 )
	    {
	      
	      // 'x' Hz bins
	      bin_t bin( wid_f , max_f , Fs[s] );
	      
	      bin.bin( mtm.f , mtm.spec );
	      
	      // output
	      for ( int i = 0 ; i < bin.bfa.size() ; i++ ) 
		{
		  //writer.level( Helper::dbl2str( bin.bfa[i] ) + "-" + Helper::dbl2str( bin.bfb[i] ) ,  globals::freq_strat  );
		  writer.level( ( bin.bfa[i] + bin.bfb[i] ) / 2.0 , globals::freq_strat );
		  writer.value( "MTM" , bin.bspec[i] );
		}
	      writer.unlevel( globals::freq_strat );
	      
	    }
	  
	  // otherwise, original entire spectrum
	  else
	    {
	      
	      for ( int i = 0 ; i < mtm.f.size() ; i++ ) 
		{
		  if ( mtm.f[i] <= max_f ) 
		    {
		      writer.level( mtm.f[i] , globals::freq_strat  );
		      writer.value( "MTM" , mtm.spec[i] );
		    }
		}
	      writer.unlevel( globals::freq_strat );
	      
	    }
	  
	} // next signal
      
      writer.unlevel( globals::signal_strat );
      
    }

  
  //
  // Epoch-wise analyses
  //
  

  if ( ! epoch_level_output )
    return;

  
  edf.timeline.first_epoch();

  
  //
  // for each each epoch 
  //

  
  int total_epochs = 0;
  
  while ( 1 ) 
     {
       
       int epoch = edf.timeline.next_epoch();      
       
       if ( epoch == -1 ) break;              
  
       ++total_epochs;
       
       interval_t interval = edf.timeline.epoch( epoch );
       
       // stratify output by epoch?
       writer.epoch( edf.timeline.display_epoch( epoch ) );
       
       //
       // Get each signal
       //
       
       for (int s = 0 ; s < ns; s++ )
	 {
	   
	   //
	   // only consider data tracks
	   //

	   if ( edf.header.is_annotation_channel( signals(s) ) )
	     continue;
	   
	   //
	   // Stratify output by channel
	   //
	   
	   writer.level( signals.label(s) , globals::signal_strat );

	   //
	   // Get data
	   //
	   
	   slice_t slice( edf , signals(s) , interval );
	   
	   const std::vector<double> * d = slice.pdata();	   
	   
	   // call MTM
	   
	   mtm_t mtm( npi , nwin );
	   
	   mtm.dB = dB;
	   
	   mtm.apply( d , Fs[s] );
	   	   
	   if ( wid_f > 0 )  // binned output
	     {	       

	       // wid_f (default 1) Hz bins
	       bin_t bin( wid_f , max_f , Fs[s] );
	       
	       bin.bin( mtm.f , mtm.spec );

	       // output	       
	       for ( int i = 0 ; i < bin.bfa.size() ; i++ ) 
		 {
		   //writer.level( Helper::dbl2str( bin.bfa[i] ) + "-" + Helper::dbl2str( bin.bfb[i] ) ,  globals::freq_strat  );
		   writer.level(  ( bin.bfa[i] + bin.bfb[i] ) / 2.0 , globals::freq_strat );
		   writer.value( "MTM" , bin.bspec[i] );
		 }
	       writer.unlevel( globals::freq_strat );
	       
	     }

	   else // full-spectrum output
	     {
	       
	       for ( int i = 0 ; i < mtm.f.size() ; i++ ) 
		 {
		   if ( mtm.f[i] <= max_f ) 
		     {
		       writer.level( mtm.f[i] , globals::freq_strat  );
		       writer.value( "MTM" , mtm.spec[i] );
		     }
		 }
	      
	       writer.unlevel( globals::freq_strat );
	       
	     }
	   
	 } // next signal
       
       writer.unlevel( globals::signal_strat );
       
     } // next epoch

  writer.unepoch();
  
}



mtm_t::mtm_t( const double npi , const int nwin ) : npi(npi) , nwin(nwin) 
{
  // by default, set to use 'adaptive weights' (2)
  kind = 2 ; 
  
  // set to use 1/(N.Fs) weights (4)
  inorm = 4 ; 
  
  display_tapers = false;
}


// void mtm_t::bin( double w , double mx_f , double fs )
// {
  
//   if ( f.size() < 2 ) return;
  
//   bfa.clear();
//   bfb.clear();  
//   bspec.clear();
  
//   // DC component  
//   bspec.push_back( spec[0] );
//   bfa.push_back( 0 );
//   bfb.push_back( 0 );

//   // other frequencies:
//   // 0-0
//   // 0 < x <= t1
//   // t1 < x <= t2
//   // etc
  
//   int num_freqs = f.size();
  
//   double df = f[1] - f[0];
//   // i.e. 2*nyquist/klen;

//   int freqwin = (int) ( w / df ) ;      

//   for (int i = 1; i < num_freqs ; i += freqwin)
//     {
      
//       double tem = 0.0;

//       int k = 0;
      
//       for (int j = i ; j < i + freqwin - 1 ; j++) 
// 	{
	  
// 	  if (j > 0 && j < num_freqs - 1) // skip DC and Nyquist
// 	    {	      	      
// 	      if ( f[j] <= mx_f )
// 		{
// 		  tem += spec[j];
// 		  k++;
// 		}
// 	    }
// 	}  
      
//       if ( k > 0 ) 
// 	{	  
// 	  bspec.push_back( tem/(double)k );
// 	  bfa.push_back( f[i-1] ); // less than 
// 	  bfb.push_back( f[i+k] ); // greater than or equal to
// 	}
      
//     }      
    
// }


// void mtm_t::smooth( double w , double fs )
// {
    
//   int num_freqs = f.size();
  
//   double df = f[1] - f[0];
  
//   int freqwin = (int) ( w / df) /2 ;      
  
//   for (int i = 0; i < num_freqs ; i++)
//     {
//       double tem = 0.0;
//       int k = 0;
//       for (int j = i - freqwin; j <= i + freqwin; j++) {
// 	if (j > 0 && j < num_freqs - 1) {
// 	  tem += spec[j];
// 	  k++;
// 	}
//       }      
//       if(k>0) { spec[i] = tem/(double)k; } // else spec[i] = spec[i];       
//     }
  
// }


void mtm_t::apply( const std::vector<double> * d , const int fs )
{
  
  std::vector<double> d2 = *d;
  
  double * data = (double*)&d2[0];
  
  // Fs is samples per second
  
  double dt = 1.0/(double)fs;
  
  int num_points = d->size();
  
  double fWidth =  npi/((double)num_points*dt);
  
  int K = (int) 2*num_points*fWidth*dt;
  
  double nyquist = 0.5/dt;
  
  int klen = mtm::get_pow_2( num_points );
  
  //  logger << "  running MTM based on " << klen << "-point FFT\n";
  
  double df = 2*nyquist/klen;  
  
  int num_freqs = 1+klen/2;
  
  int npoints = num_points;
  
  int k = 1;
  
  // mean-center
  
  if ( 0 ) 
    {
      double mean = mtm::remove_mean( data, npoints ); 
    }

  spec.resize( klen ,  0 );  
  
  std::vector<double> dof( klen );
  std::vector<double> Fvalues( klen );
  
  mtm::do_mtap_spec(&(data)[0], npoints, kind,  nwin,  npi, inorm, dt,
		    &(spec)[0], &(dof)[0], &(Fvalues)[0], klen , display_tapers );
  
  // shrink to positive spectrum 
  // and scale x2 for 
  spec.resize( num_freqs );
  
  f.resize( num_freqs , 0 );
  
  for (int i = 0; i < num_freqs; i++)
    {
      
      f[i] = df*i;
      
      if ( i > 0 && i < num_freqs - 1 ) 
	spec[i] *= 2;
      
      // report dB?
      if ( dB ) spec[i] = 10 * log10( spec[i] );
      
      // std::cerr << i << "\t" 
      //        		<< f[i] << "\t" 
      // 	 	<< spec[i] << "\t"
      //        		<< 10*log10(spec[i]) << "\n";
      
    }  
  
}

