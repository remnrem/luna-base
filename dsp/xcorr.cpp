
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


#include "dsp/xcorr.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "db/db.h"
#include "helper/logger.h"
#include "fftw/fftwrap.h"

extern writer_t writer;
extern logger_t logger;


// wrapper
void dsptools::xcorr( edf_t & edf , param_t & param )
{

  //
  // get signals
  //

  std::string signal_label = param.requires( "sig" );

  //
  // epoch-level or whole-trace
  //

  const bool epoch_output = param.has( "epoch" ) ;
  
  //
  // analysis/output options
  //

  // window (in seconds ) 0 = none
  const double w_sec = param.has( "w" ) ? param.requires_dbl( "w" ) : 0 ;

  // center (in seconds) 0 = symm / usual

  // report / find max only for  c-w to c+w
  const double c_sec = param.has( "c" ) ? param.requires_dbl( "c" ) : 0 ;
  
  const bool verbose = param.has( "verbose" );


  //
  // get signals
  //
      
  const bool no_annotations = true;

  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );

  const int ns = signals.size();

  if ( ns == 0 )
    {
      logger << "  *** none of the requested signals found... bailing\n";
      return;
    }

  //
  // sample rates
  //

  std::vector<double> FsV = edf.header.sampling_freq( signals );
  
  const int Fs = FsV[0];
  for (int i=0; i<signals.size(); i++)
    if ( FsV[i] != Fs )
      Helper::halt( "unequal sampling frequencies" );
  
  
  //
  // iterate over epochs
  //
  
  int ne = edf.timeline.first_epoch();
  
  std::map<int,std::map<int,std::map<int,double> > > xcorr;
  std::map<int,std::map<int,int> > delay;
  std::map<int,std::map<int,std::vector<double> > > delay_tracker;
  
  int cnt_epoch = 0 ;

  int mxlag = Fs * w_sec;
  int cent = Fs * c_sec;
  
  while ( 1 )
    {
      
      int epoch = edf.timeline.next_epoch();
      
      if ( epoch == -1 ) break;
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      ++cnt_epoch;
      
      if ( epoch_output )
	writer.epoch( edf.timeline.display_epoch( epoch ) );		
      
      // get data      
      matslice_t mslice( edf , signals , interval );
      
      // need to make modifiable versions
      std::vector<std::vector<double> > D( ns );
      for (int i=0; i<ns; i++)
	D[i] =  *mslice.col(i);

      
      // consider all pairs      
      for (int s1=0; s1 < ns-1; s1++)
	{

	  // epoch-level outputs?
	  if ( epoch_output )
	    writer.level( signals.label(s1) , globals::signal1_strat );
	  
	  for (int s2=s1+1; s2 < ns ; s2++)
	    {
	      
	      if ( epoch_output )
		writer.level( signals.label(s2) , globals::signal2_strat );

	      // calculate cross-correl
	      xcorr_t xc( D[s1], D[s2] , mxlag , cent );
	      
	      // track delay (and mirror) [ for means ] 
	      delay[s1][s2] += xc.lags[xc.mx]; 
	      delay[s2][s1] -= xc.lags[xc.mx];
	      
	      // for medians
	      delay_tracker[s1][s2].push_back( xc.lags[xc.mx] );
	      delay_tracker[s2][s1].push_back( -xc.lags[xc.mx] );
	      
	      // output?
	      if ( epoch_output)
		writer.value( "D" , xc.lags[xc.mx] / (double)Fs );
	      
	      // track Xcorr
	      for (int i=0; i<xc.lags.size(); i++)
		{
		  // in reporting range? this will fill in symm pair also 
		  if ( xc.lags[i] >= cent - mxlag && xc.lags[i] <= cent + mxlag )
		  {
		      xcorr[s1][s2][ xc.lags[i] ] += xc.C[i];
		      xcorr[s2][s1][ - xc.lags[i] ] += xc.C[i];
		    }
		}
	      
	    }	  
	  if ( epoch_output ) writer.unlevel( globals::signal2_strat );	  
	}
      if ( epoch_output ) writer.unlevel( globals::signal1_strat );
      
      // next epoch
    } 
  if ( epoch_output ) writer.unepoch();
  

  
  //
  // Report averaged results
  //
      
  for (int s1=0; s1<ns; s1++)
    {
      
      writer.level( signals.label(s1) , globals::signal1_strat );
      
      for (int s2=0; s2<ns; s2++)
	{

	  if ( s1 == s2 ) continue;
	  
	  writer.level( signals.label(s2) , globals::signal2_strat );

	  // D based on average of epochs
	  writer.value( "D_MN" , ( delay[ s1 ][ s2 ] / (double)cnt_epoch ) / (double)Fs );
	  writer.value( "S_MN" , ( delay[ s1 ][ s2 ] / (double)cnt_epoch ) );

	  // also based on medians
	  const double med_d = MiscMath::median( delay_tracker[ s1 ][ s2 ]  );
	  writer.value( "D_MD" , med_d / (double)Fs );
	  writer.value( "S_MD" , med_d );
	  
	  // also get D based on the average of lagged XCs
	  double mx = 0;
	  int mxi = 0;
	  bool first = true;
	  std::map<int,double>::const_iterator ff = xcorr[s1][s2].begin();	  
	  while( ff != xcorr[s1][s2].end() )
	    {
	      if ( first )
		{
		  mx = fabs( ff->second );
		  mxi = ff->first;
		  first = false;
		}
	      else if ( fabs( ff->second ) > mx )		
		{
		  mx = fabs( ff->second );
		  mxi = ff->first;
		}
	      ++ff;
	    }

	  // delay from averaged Xcs
          writer.value( "D" , mxi / (double)Fs );
	  writer.value( "S" , mxi  );	  
	  
	  if ( verbose )
	    {
	      std::map<int,double>::const_iterator ff = xcorr[s1][s2].begin();
	      
	      while( ff != xcorr[s1][s2].end() )
		{
		  // delay
		  writer.level( ff->first , "D" );
		  writer.value( "T" , ff->first / (double)Fs );
		  writer.value( "XCORR" , ff->second / (double)cnt_epoch );
		  ++ff;
		}
	      writer.unlevel( "D" );
	    }	      
	}
      writer.unlevel( globals::signal2_strat );
    }
  writer.unlevel( globals::signal1_strat );
  
}




// primary function : pairwise xcorr via FFT
xcorr_t::xcorr_t( std::vector<double> a ,
		  std::vector<double> b ,
		  const int mxlag ,
		  const int cent )
{
  
  const int na = a.size();
  const int nb = b.size();

  const int nm = na > nb ? na : nb ;
  
  // need to zero padded?

  if ( na != nm ) a.resize( nm , 0 );
  if ( nb != nm ) b.resize( nm , 0 );
  
  // transform both vectors
  const int Fs = 100; // arb;

  long int np2 = MiscMath::nextpow2( 2 * nm - 1 );
  
  FFT ffta( nm , np2 , Fs , FFT_FORWARD, WINDOW_NONE );
  FFT fftb( nm , np2 , Fs , FFT_FORWARD, WINDOW_NONE );
  ffta.apply( a );
  fftb.apply( b );
  
  // Extract the raw transform
  std::vector<std::complex<double> > At = ffta.transform();
  std::vector<std::complex<double> > Bt = fftb.transform();

  // std::cout << "A(t)\n";
  // for (int i=0; i<At.size() ; i++)
  //   std::cout << " A " << At[i] << "\n";
  
  // std::cout << "B(t)\n";
  // for (int i=0; i<Bt.size() ; i++)
  //   std::cout << " B " << Bt[i] << "\n";
  
  // compute cross-correlation
  //c = ifft(X.*conj(Y));
  const int nt = At.size();
  
  std::vector<std::complex<double> > I( nt );
  for (int i=0; i<nt; i++)
    I[i] = At[i] * std::conj( Bt[i] );
  
  
  // std::cout << "I(t)\n";
  // for (int i=0; i<I.size() ; i++)
  //   std::cout << " I " << I[i]	<< "\n";
  
  // iFFT
  //  FFT ifft( n_conv_pow2 , n_conv_pow2 , 1 , FFT_INVERSE );
  //  std::cout << " I sx " << I.size() << " " << np2 << "\n";
  FFT ifft( np2, np2 , Fs , FFT_INVERSE );
  ifft.apply( I );
  
  std::vector<double> C0 = ifft.inverse();
  
  // std::cout << "C0(t)\n";
  // for (int i=0; i<C0.size() ; i++)
  //   std::cout << " C0 " << C0[i]  << "\n";


  
  // lags
  const int maxlag = nm - 1 ; // M-1
  const int nl = maxlag * 2 + 1 ;
  const int nc = C0.size();

  //std::cout << " maxlag nl nm nc " << maxlag << " " << nl << " " << nm << " " << nc << "\n";

  // reorder
  C.resize( nl );
  lags.resize( nl );

  //c = [c(end-maxlag+1:end,:);c(1:maxlag+1,:)];

  // for (int k = 0 ; k < nc ; k++)
  //   std::cout << k << "\t" << C[k] << "\n";
  
  double mc = 0;
  mx = 0;
  
  int lag = -maxlag;
  int idx = 0;
  for (int k = nc - maxlag ; k < nc ; k++)
    {
      C[idx] = C0[k];

      if ( mxlag == 0 || ( lag >= cent - mxlag && lag <= cent + mxlag ) )
	{
	  if ( fabs( C[idx] ) > mc )
	    {
	      mc = fabs( C[idx] );
	      mx = idx;
	    }
	}
      
      lags[idx] = lag;
      ++lag;
      ++idx;
    }
  
  for (int k = 0 ; k < maxlag + 1 ; k++)
    {
      C[idx] = C0[k] ;

      if ( mxlag == 0 || ( lag >= cent - mxlag && lag <= cent + mxlag ) )
	{
	  if ( fabs( C[idx] ) > mc )
	    {
	      mc = fabs( C[idx] );
	      mx = idx;
	    }
	}
      
      lags[idx] = lag;
      ++lag;
      ++idx;
    }
  
  //  std::cout << " mx " << mx << " " << lags[mx] << " " << C[mx] << "\n";
}
