
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

#include "dsp/linedenoiser.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "fftw/fftwrap.h"

void dsptools::line_denoiser( edf_t & edf , param_t & param )
{
  
  // signal(s)
  signal_list_t signals = edf.header.signal_list( param.value( "sig" ) );
  const int ns = signals.size();

  // target frequencies
  if ( ! param.has( "f" ) ) Helper::halt( "no 'f' option for LINE-DENOISE" );
  std::vector<double> f = param.dblvector( "f" );

  // do by epoch?
  bool epoch = param.has( "epoch" );
  if ( epoch ) logger << "  iterating over epochs\n";
  else logger << "  correcting for entire signal\n";
	 
  // options
  double w_noise = 1;
  double w_neigh = 1;
  if ( param.has( "w" ) )
    {
      std::vector<double> w = param.dblvector("w");
      if ( w.size() != 2 )
	Helper::halt( "requires 'w' to be a two-element vector" );
      w_noise = w[0];
      w_neigh = w[1];
    }
  
  logger << "  running line denoiser for " << f.size() << " target frequencies\n"
	 << "  noise/neighbour band width " << w_noise << " and " << w_neigh << " Hz respectively\n";

  // process each signal
  for (int s=0; s<ns; s++)
    {
      // only process data channels
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      // get sample rate
      int sr = edf.header.sampling_freq( signals(s) );

      // get whole signal  (although we update epoch-by-epoch)
      // we need to this edit and return back 

      slice_t slice0( edf , signals(s) , edf.timeline.wholetrace() );
      std::vector<double> orig = * slice0.pdata();
  
      // do by epoch
      int ne = epoch ? edf.timeline.first_epoch() : 1 ; 
	
      // iterate over each epoch
      
      // compile new signal (by epochs)
      std::vector<std::vector<double> > filt;
      
      while ( 1 )
	{
	  
	  // next epoch
	  int epoch = epoch ? edf.timeline.next_epoch() : 1 ; 

	  // all done?
	  if ( epoch == -1 ) break;
	  
	  // get data
	  interval_t interval = epoch ? edf.timeline.epoch( epoch ) : edf.timeline.wholetrace();

	  slice_t slice( edf , signals(s) , interval );

	  const std::vector<double> * data = slice.pdata();
	  
	  // do procedure, and store
	  filt.push_back( line_denosier( data , sr , f , w_noise , w_neigh ) );
	  
	  // next epoch
	}
      
      // update signal
      int p = 0;
      const int n = orig.size();
      for (int e=0; e<filt.size(); e++)
	for (int i=0; i<filt[e].size(); i++)
	  orig[p++] = filt[e][i];
      
      logger << "  updating " << signals.label(s) << "\n";

      // and send back to the EDF
      edf.update_signal( signals(s) , &orig );

      if ( ! epoch ) break;
    }
			    
  
}

std::vector<double> line_denosier( const std::vector<double> * x ,
				   const int Fs ,
				   const std::vector<double> & fl ,
				   const double w_noise ,
				   const double w_neigh )
{
  
  const int n = x->size();
  // DFT of data
  real_FFT eegfft;
  eegfft.init( n , n , Fs );
  eegfft.apply( *x );
  std::vector<dcomp> eegfftX = eegfft.transform();
  
  // frqs: eegfft.apply
  std::vector<double> frq = eegfft.frq;
  int cutoff = eegfft.cutoff;
  int nfreqs = frq.size();

  double fmax = frq[ frq.size() - 1 ] ; 

  // interpolate each region
  const int nf = fl.size();
  for (int i=0;i<nf;i++)
    {

      // get windows      
      double flwr = fl[i] - w_noise;
      double fupr = fl[i] + w_noise;

      if ( flwr < 0 || fupr > fmax ) continue;
	   
      int lwr_idx = 0 , upr_idx = 0;
      for ( int j=0; j<cutoff; j++)
	if ( frq[j] >= flwr ) { lwr_idx = j; break; }
      for ( int j=lwr_idx; j<cutoff; j++)
	if ( frq[j] > fupr ) { upr_idx = j-1; break; }

      // neighbiours :  left --> flwr
      //                fupr --> right
      double left = fl[i] - w_noise - w_neigh;
      double right= fl[i] + w_noise + w_neigh;;

      if ( left < 0 ) left = 0;
      if ( right > fmax ) right = fmax;
      
      int left_idx = 0 , right_idx = 0;
      for ( int j=0; j<cutoff; j++)
	if ( frq[j] >= left ) { left_idx = j; break; }
      for ( int j=upr_idx; j<cutoff; j++)
	if ( frq[j] > right ) { right_idx = j-1; break; }
      
      // get left and right means
      double mn = 0, nn = 0;
      // left
      for ( int j=left_idx; j<lwr_idx; j++)
	{
	  mn += eegfft.mag[j];
	  ++nn;
	}
      
      // right
      for ( int j=upr_idx + 1 ; j<=right_idx; j++)
	{
	  mn += eegfft.mag[j];
	  ++nn;
	}
      
      mn /= nn;

      //% Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
      //data_fft(:,smpl2int) = bsxfun(@times, exp(bsxfun(@times,angle(data_fft(:,smpl2int)),1i)), mns4int);

      // do interpolation (simple mean) of amplitude, followed by Euler's formula
      // to replace complex values
      for ( int j=lwr_idx; j<=upr_idx; j++)
	eegfftX[j] = exp( dcomp( 0 , std::arg( eegfftX[j] ) ) ) * mn ; 

    }
  
  //% complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
  //% to ensure a real valued signal after iFFT
  //filt = ifft(data_fft,[],2,'symmetric');

  real_iFFT ifft( n, n , Fs );
  ifft.apply( eegfftX );
    
  std::vector<double> f = ifft.inverse();
  return f;

}
