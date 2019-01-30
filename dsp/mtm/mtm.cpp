
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

#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "nrutil.h"

extern writer_t writer;
extern logger_t logger; 

void mtm::wrapper( edf_t & edf , param_t & param )
{
  
  std::string signal_label = param.requires( "signal" );

  signal_list_t signals = edf.header.signal_list( signal_label );    

  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  const int ns = signals.size();

  // other parameters
  bool epoch_level_output = param.has( "epoch" );
  bool spectrum = param.has( "spectrum" );
  
  // MTM parameters
  int npi = param.has( "nw" ) ? param.requires_int( "nw" ) : 3 ;
  int nwin = param.has( "t" ) ? param.requires_int( "t" ) : 2*npi-1 ;
  
  double max_f = param.has( "max" ) ?  param.requires_dbl( "max" ) : 30;
  double wid_f = param.has( "bin" ) ? param.requires_dbl( "bin" ) : ( spectrum ? 0 : 1 ) ;
  
  
  
  //
  // Whole signal analyses
  //


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
      
      // call MTM
      
      mtm_t mtm( npi , nwin );

      mtm.apply( d , Fs[s] );
	   	   
      if ( wid_f > 0 )
	{
	  
	  // 1 Hz bins
	  mtm.bin( wid_f , max_f , Fs[s] );
	  
	  // output
	  for ( int i = 0 ; i < mtm.bfa.size() ; i++ ) 
	    {
	      writer.level( Helper::dbl2str( mtm.bfa[i] ) + "-" + Helper::dbl2str( mtm.bfb[i] ) ,  globals::freq_strat  );
	      writer.value( "MIDF" , ( mtm.bfa[i] + mtm.bfb[i] ) / 2.0 );
	      writer.value( "MTM" , mtm.bspec[i] );
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
	   
	   mtm.apply( d , Fs[s] );
	   
	   
	   if ( wid_f > 0 )  // binned output
	     {	       

	       // wid_f (default 1) Hz bins
	       mtm.bin( wid_f , max_f , Fs[s] );
	       
	       // output	       
	       for ( int i = 0 ; i < mtm.bfa.size() ; i++ ) 
		 {
		   writer.level( Helper::dbl2str( mtm.bfa[i] ) + "-" + Helper::dbl2str( mtm.bfb[i] ) ,  globals::freq_strat  );
		   writer.value( "MIDF" , ( mtm.bfa[i] + mtm.bfb[i] ) / 2.0 );
		   writer.value( "MTM" , mtm.bspec[i] );
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



mtm_t::mtm_t( const int npi , const int nwin ) : npi(npi) , nwin(nwin) 
{
  
  // by default, set to use 'adaptive weights' (2)
  kind = 2 ; 
  
  // by default, set to use 1/N weights 
  inorm = 0 ; 

  display_tapers = true;
}

void mtm_t::bin( double w , double mx_f , double fs )
{

  if ( f.size() < 2 ) return;

  bfa.clear();
  bfb.clear();  
  bspec.clear();
  
  int num_freqs = f.size();
  
  double df = f[1] - f[0];
  // i.e. 2*nyquist/klen;

  int freqwin = (int) ( w / df ) ;      

  for (int i = 1; i < num_freqs ; i += freqwin)
    {

      double tem = 0.0;
      int k = 0;
      
      for (int j = i ; j < i + freqwin - 1 ; j++) 
	{
	  if (j > 0 && j < num_freqs - 1) 
	    {	      
	      if ( f[j] <= mx_f )
		{
		  tem += spec[j];
		  k++;
		}
	    }
	}  
      
      if ( k > 0 ) 
	{	  
	  bspec.push_back( tem/(double)k );
	  bfa.push_back( f[i-1] );
	  bfb.push_back( f[i+k] );
	}
      
    }      
    
}


void mtm_t::smooth( double w , double fs )
{
    
  int num_freqs = f.size();
  
  double df = f[1] - f[0];
  
  int freqwin = (int) ( w / df) /2 ;      
  
  for (int i = 0; i < num_freqs ; i++)
    {
      double tem = 0.0;
      int k = 0;
      for (int j = i - freqwin; j <= i + freqwin; j++) {
	if (j > 0 && j < num_freqs - 1) {
	  tem += spec[j];
	  k++;
	}
      }      
      if(k>0) { spec[i] = tem/(double)k; } // else spec[i] = spec[i];       
    }

}


void mtm_t::apply( const std::vector<double> * d , const int fs )
{

  std::vector<double> d2 = *d;
  
  double * data = (double*)&d2[0];
  
  // Fs is samples per second
  
  double dt = 1.0/(double)fs;
  
  int num_points = d->size();

  double fWidth =  npi/((double)num_points*dt);
  
  //  std::cout << "fWidth = " << fWidth << "\n";

  int K = (int) 2*num_points*fWidth*dt;
  
  double nyquist = 0.5/dt;
  
  int klen = mtm::get_pow_2( num_points );

  double df = 2*nyquist/klen;  

  int num_freqs = 1+klen/2;
  
  int npoints = num_points;

  int k = 1;

  // mean-center

  double mean = mtm::remove_mean( data, npoints ); 
  
  spec.resize( klen ,  0 );  
  
  std::vector<double> dof( klen );
  std::vector<double> Fvalues( klen );
  
  mtm::do_mtap_spec(&(data)[0], npoints, kind,  nwin,  npi, inorm, dt,
		    &(spec)[0], &(dof)[0], &(Fvalues)[0], klen , display_tapers );
  
  // shrink to positive spectrum 
  spec.resize( num_freqs );
  
  f.resize( num_freqs , 0 );
  
  for (int i = 0; i < num_freqs; i++)
    {
      
      f[i] = df*i;
      
      // std::cerr << i << "\t" 
      //        		<< f[i] << "\t" 
      // 	 	<< spec[i] << "\t"
      //        		<< 10*log10(spec[i]) << "\n";
      
    }
    

}



void mtm_t::apply2( const std::vector<double> * d , const int fs )
{

  std::vector<double> d2 = *d;
  
  double * data = (double*)&d2[0];
  
  // Fs is samples per second
  
  double dt = 1.0/(double)fs;
  
  int num_points = d->size();

  // time track
  
  std::vector<double> ex( num_points );
  
  for (int i=0;i<num_points;i++) ex[i] = i * dt;
  
  double fWidth =  npi/((double)num_points*dt);

  int K = (int) 2*num_points*fWidth*dt;
  
  double nyquist = 0.5/dt;
  
  int klen = mtm::get_pow_2( num_points );

  //  std::cout << "klen = " << klen << "\n";

  int num_freqs = 1+klen/2;
  
  int npoints = num_points;

  int k = 1;

  if ( 0 ) 
    {
      double mean = mtm::remove_mean( data, npoints ); 
    }

  //
  // simple (naive) periodogram
  //
  
  std::vector<double> naive_spec( num_freqs );
  
  std::vector<double> dtemp( klen ); 
  
  // 10% cosine taper

  std::cout << "lookig at taper\n";
  
  for (int i = 0; i < num_points; i++)
    {
      double vwin = mtm::get_cos_taper(num_points, i, .05); 
      dtemp[i] = vwin*data[i];
      
      if ( i < 10 || i > num_points - 10 ) 
	std::cout << i << "\t" << vwin <<  "\t" << dtemp[i] << "\t" << data[i] << "\n" ;
    }
  


  double anrm = num_points;
  
  switch (inorm)
    {
    case 0:
      anrm = 1.;
      break;
    
    case 1:
      anrm = num_points;
      break;

    case 2:
      anrm = 1 / dt;
      break;
    case 3:
      anrm = sqrt((double) num_points);
      break;
    default:
      anrm = 1.;
      break;
    }

  
  double norm = 1./(anrm*anrm);
  
  std::cout << "NORM = " << norm << " (inorm " << inorm << ")\n";

  mtm::zero_pad(&(dtemp)[0], num_points, klen);

  int isign = 1;
  mtm::jrealft(&(dtemp)[0]-1, (unsigned long) klen, isign);
  
  for(int i=1; i<num_freqs-1; i++)
    naive_spec[i] = norm*(SQR(dtemp[2*i+1])+SQR(dtemp[2*i]));
  
  naive_spec[0] = norm*SQR(fabs(dtemp[0]));

  naive_spec[num_freqs-1] = norm*SQR(fabs(dtemp[1]));

  double df = 2*nyquist/klen;

  int freqwin = (int) ( fWidth/df) /2 ;      
  

#if 1

  // smooth the periodogram 

  fprintf(stderr, "smooth the periodogram 4, freqwin=%d\n", freqwin);
  
  for (int i = 0; i < num_freqs ; i++)
    {
      double tem = 0.0;
      k = 0;
      for (int j = i - freqwin; j <= i + freqwin; j++) {
	if (j > 0 && j < num_freqs - 1) {
	  tem += naive_spec[j];
	  k++;
	}
      }
      
      if(k>0) { dtemp[i] = tem/(double)k; } else  dtemp[i] = naive_spec[i]; 
      
    }
  
  /*for (i = 1; i < num_freqs - 1; i++) naive_spec[i] = dtemp[i];*/
  
#endif


  
  for (int i = 0; i < num_freqs ; i++)
    {
      
      if(  naive_spec[i] < 0.0 || dtemp[i] < 0.0 )
	{
	  fprintf(stderr,"negative or zero spectrum: %d\n",i);
	  fprintf(stderr,"%g  %g\n", naive_spec[i], dtemp[i]);
	  exit(0);
	}
      
      naive_spec[i] = 10.*log10(naive_spec[i]);
      dtemp[i] = 10.*log10(dtemp[i]);
      
    }
  
  /**********************************************/

  std::vector<double> spec( klen );
  std::vector<double> dof( klen );
  std::vector<double> Fvalues( klen );

  std::cout << "\n\nentering MTM\n\n";

  mtm::do_mtap_spec(&(data)[0], npoints, kind,  nwin,  npi, inorm, dt,
		    &(spec)[0], &(dof)[0], &(Fvalues)[0], klen , true );

  std::cerr << " done with do_mtap_spec: " << num_freqs << "\n";
  
  for (int i = 0; i < num_freqs; i++)
    {
      
      double frq1 =  df*i;
      
      std::cout << i << "\t" 
		<< frq1 << "\t" 
		<< spec[i] << "\t"
		<< 10*log10(spec[i]) << "\n";
		
// 		<< "\t" << naive_spec[i] << "\t"
// 		<< dtemp[i] << "\t" << dof[i] << "\t"
// 		<< Fvalues[i] << "\n";
      
    }
  
  
}

