
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


#include "pac.h"
#include "param.h"

#include "cwt/cwt.h"
#include "miscmath/crandom.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/helper.h"

#include <iostream>

extern writer_t writer;

void dsptools::pac( edf_t & edf , param_t & param )
{

  std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );    
  
  const int ns = signals.size();

  // get freq for phase, then for power
  std::vector<std::string> sf4p = param.strvector( "ph" );
  std::vector<std::string> sf4a = param.strvector( "amp" );
    
  std::vector<double> f4p, f4a;
  for (int i=0;i<sf4p.size();i++) 
    {
      double x;
      if ( ! Helper::str2dbl( sf4p[i] , &x ) ) Helper::halt( "not valid freq for ph=" );
      f4p.push_back( x );
    }		
			       
  for (int i=0;i<sf4a.size();i++) 
    {
      double x;
      if ( ! Helper::str2dbl( sf4a[i] , &x ) ) Helper::halt( "not valid freq for amp=" );
      f4a.push_back( x );
    }		
  
  if ( f4p.size() == 0 || f4a.size() == 0 ) Helper::halt( "requires 'ph=' and 'amp=' parameters" );

  // by default, 1000 replicates for permutation
  const int nreps = param.has("nreps") ? param.requires_int( "nreps" ) : 1000;

  // using epochs or not?
  bool epoched = param.has( "epoch" );
  
  // for each signal
  for (int s=0;s<ns;s++)
    {
      
      const int srate = edf.header.sampling_freq( signals(s) ) ; 

      //
      // Output      
      //

      writer.level( signals.label(s) , globals::signal_strat );

      
      logger << "  running PAC...\n";
      
      //
      // either for each epoch, or for entire trace
      //

      while ( 1 ) 
	{
	  
	  //
	  // fetch either an EPOCH or the entire timeline
	  //
	  
	  interval_t interval;

	  int epoch = -1;
	  
	  if ( epoched )
	    {
	      epoch = edf.timeline.next_epoch();
	      if ( epoch == -1 ) break;
	      interval = edf.timeline.epoch( epoch );
	      
	      writer.epoch( edf.timeline.display_epoch( epoch ) );

	    }
	  else
	    {
	      interval = edf.timeline.wholetrace();
	    }
	  
	  //
	  // Fetch data slice
	  //
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  const std::vector<double> * signal = slice.pdata();


	  //
	  // Calculate PAC
	  //
	  
	  pac_t pac( signal , f4p , f4a , srate , nreps );
	  
	  bool okay = pac.calc();
	  
	  if ( ! okay )
	    Helper::halt( "problem in PAC calculation" );


	  //
	  // Output, looping over frequencies
	  //
	  
	  for (int a=0;a<f4p.size();a++)
	    for (int b=0;b<f4a.size();b++)
	      {
		writer.level( Helper::dbl2str( f4p[a] ) + "x" + Helper::dbl2str( f4a[b] ) , "FRQS" );
		writer.value( "ZPAC" , pac.zpac(a,b) );
		writer.value( "PPAC" , pac.ppac(a,b) );
	      }
	  
	  writer.unlevel( "FRQS" );

	  //
	  // Done for this signal ?
	  //

	  if ( ! epoched ) break;
	}
      
      if ( epoched ) writer.unepoch();
      
      // next signal      
    }

  

  writer.unlevel( globals::signal_strat );

  //
  // all done
  //

}


bool pac_t::calc() 
{
  
  for (int fa=0;fa<frq4phase.size();fa++)
    for (int fb=0;fb<frq4pow.size();fb++)
      {	
	
	// for wavelets, for high temporal resolution, set a relatively 
	// small number of cycles
	
	//	const int n_cycles = 4;
	const int n_cycles = 7;
	
	//
	// wavelet for phase
	//

	CWT phase_cwt;	
	phase_cwt.set_sampling_rate( srate );  
	phase_cwt.add_wavelet( frq4phase[fa] , n_cycles );  
	phase_cwt.load( data );
	phase_cwt.run();
	
	//
	// wavelet for power
	//
	
	CWT pow_cwt;	
	pow_cwt.set_sampling_rate( srate );  
	pow_cwt.add_wavelet( frq4pow[fb] , n_cycles );  
	pow_cwt.load( data );
	pow_cwt.run();
	
	//
	// Extract
	//

	std::vector<double> angle = phase_cwt.phase(0);
	std::vector<double> pwr   = pow_cwt.results(0); // get 'raw' power back
	
	const int n = angle.size();

	// std::cout << "sizes = " << angle.size() << " " << pwr.size() << "\n";
	// for (int i = 0 ; i < angle.size() ; i++ ) 
	//   std::cout << "PDET\t" << i << "\t"
	// 	    << angle[i] << "\t"
	// 	    << pwr[i] << "\n";
	
	// Note:   pwr is different compared to example
	// Note:   angle is offset by 1 compared to example
	
	// Note: various differences:  time resolution selected, in examples of using different N versus nextPow2 for convolution size, etc
	// In broad terms, things line up

	//
	// Precalculate x and exp(y) 
	//

	std::vector<std::complex<double> > x(n);
 
	std::vector<std::complex<double> > y(n);

	for (int i=0;i<n;i++)
	  {
	    x[i] = std::complex<double>( pwr[i] , 0 );
	    y[i] = exp( std::complex<double>( 0 , angle[i] ) );
	  }
	
	//
	// Calculate PAC
	//
	
	// obsPAC = abs(mean(pwr.*exp(1i*phase)));
	
	dcomp sm = 0;
	
	for (int i=0; i<n; i++)
	  sm += x[i] * y[i] ;
	
	double pac = abs( sm / (double)n );
	
	//
	// Permute
	// 
	
	std::vector<double> ppac( nreps , 0 ); // permuted PACs
	double p = 1;
	
	for (int r=0; r<nreps ; r++ ) 
	  {

	    // from 0..n, select a random time-point (from within 10-90% of signal)

	    int pp = n * 0.1 + CRandom::rand( int( n * 0.8 ) );
	    
	    dcomp sm = 0;
	
	    for (int i=0; i<n; i++)
	      {

		sm += x[pp] * y[i] ;
		
		// advance permuted position, 
		++pp;

		// wrapping around when needed
		if ( pp == n ) pp = 0;
	      }
	    
	    ppac[r] = abs( sm / (double)n );
	    
	    if ( ppac[r] >= pac ) ++p;
	    
	  }
	
	p /= (double)(nreps+1);
	
	//
	// Z-transform observed PAC
	//

	double pacmean = MiscMath::mean( ppac );
	double pacsd   = MiscMath::sdev( ppac );
	double pacZ    = ( pac - pacmean ) / pacsd;
	
	// store
	z[fa][fb] = pacZ;
	pval[fa][fb] = p;
      }
  
  return true;
  
}
