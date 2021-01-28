
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

#include "dsp/psi.h"
#include "miscmath/miscmath.h"
#include <iostream>
#include "defs/defs.h"
#include "fftw/fftwrap.h"
#include "stats/matrix.h"
#include "stats/statistics.h"

#include "db/db.h"
#include "helper/logger.h"
#include "edf/edf.h"
#include "edf/slice.h"

// This implements the phase slope index (PSI) as described here:
//    Nolte G, Ziehe A, Nikulin VV, Schl\"ogl A, Kr\"amer N, Brismar T, M\"uller KR. 
//    Robustly estimating the flow direction of information in complex physical systems. 
//    Physical Review Letters. To appear. 
// This is based on their Matlab implementation, available: http://doc.ml.tu-berlin.de/causality/

extern writer_t writer;
extern logger_t logger;


void dsptools::psi_wrapper( edf_t & edf , param_t & param )
{

  //
  // get signals
  //

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );

  if ( signals.size() < 2 ) return;

  const int ns = signals.size();

  //
  // check sample rates
  //

  std::vector<double> Fs = edf.header.sampling_freq( signals );

  int sr = Fs[0];
  for (int s=1;s<ns;s++)
    if ( Fs[s] != sr )
      Helper::halt( "all sampling rates must be similar for PSI" );
  
  //
  // Epoch or whole trace?
  //

  bool by_epoch = param.has( "epoch" );
  
  int ne = edf.timeline.first_epoch();

  logger << "  running within " << ne << " " << edf.timeline.epoch_length() << " second epochs\n";

  //
  // Frequency bins?
  //

  std::vector<int> lwr = param.intvector( "f-lwr" );
  std::vector<int> upr = param.intvector( "f-upr" );
  if ( lwr.size() != upr.size() ) 
    Helper::halt( "f-lwr and f-upr have different sizes" );
  for (int i=0;i<lwr.size();i++)
    if ( lwr[i] >= upr[i] ) Helper::halt( "f-lwr >= f-upr" );
  bool has_freqs = lwr.size() > 0 ;


  //
  // PSI epoch/segment lengths 
  //

  // default 4 seconds epochs, with 2-second segments (w/ 50% overlap) 
  int eplen = param.has( "eplen" ) ? sr * param.requires_dbl( "eplen" ) : sr * 4 ;
  
  // by default, segment length is 1/2 of epoch length (and w/ 50% overlap of segments) 
  int seglen = param.has( "seglen" ) ? sr * param.requires_dbl( "seglen" ) :  eplen / 2 ;
  
  logger << "  running PSI with " << eplen << " samples per epoch, " << seglen << " per segment\n";

  //
  // Get data
  //

  if ( ! by_epoch )
    {

      logger << "  running across entire trace\n";
      
      matslice_t mslice( edf , signals , edf.timeline.wholetrace() );
      
      const Data::Matrix<double> & X = mslice.data_ref();

      psi_t psi( &X , eplen , seglen );
      
      if ( has_freqs )
	for (int i=0;i<lwr.size();i++)
	  psi.add_freqbin( lwr[i] , upr[i] );

      psi.calc();
  
      psi.report( signals );

      return; 
    }


  //
  // By epoch
  //

  
  while ( 1 ) 
    {

      
      int epoch = edf.timeline.next_epoch();      
      
      if ( epoch == -1 ) break;
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      matslice_t mslice( edf , signals , interval );
      
      const Data::Matrix<double> & X = mslice.data_ref();
      
      psi_t psi( &X , eplen , seglen );

      if ( has_freqs )
	for (int i=0;i<lwr.size();i++)
	  psi.add_freqbin( lwr[i] , upr[i] );
      
      psi.calc();
      
      writer.epoch( edf.timeline.display_epoch( epoch ) );
      
      psi.report( signals );

    }

  writer.unepoch();
    
}



void psi_t::report( const signal_list_t & signals )
{
  const double EPS = 1e-8;

  // note: need to tidy for text-table output

  if ( n_models == 0 ) return;
  
  for (int m=0; m<n_models; m++)
    {
      writer.level( m+1 , "M" );

      writer.value( "F1" , freqbins[m][0] );
      writer.value( "F2" , freqbins[m][ freqbins[m].size()-1 ] );
      
      // channel sums
      for (int i=0;i<nchan;i++)
	{
	  writer.level( signals.label(i) , globals::signal_strat );
	  writer.value( "PSI_RAW" , psi_sum[m][i] );
	  writer.value( "STD" , std_psi_sum[m][i] );
	  writer.value( "PSI" , psi_sum[m][i]  / ( EPS + std_psi_sum[m][i] ) );
	}
      writer.unlevel( globals::signal_strat );
      
      // channel pairs
      for (int i=0;i<nchan;i++)
	{
	  writer.level( signals.label(i) , "CH1" );
		    
	  for (int j=0;j<nchan;j++)
	    {
	      if ( i == j ) continue;
	      writer.level( signals.label(j) , "CH2" );

	      writer.value( "PSI_RAW" , psi[m](i,j) );
	      writer.value( "STD" , std_psi[m](i,j) );
	      writer.value( "PSI" , psi[m](i,j)  / ( EPS + std_psi[m](i,j) ) );
	    }
	  writer.unlevel("CH2");
	}
      writer.unlevel("CH1");
    }
  writer.unlevel( "M" );

  
}

void psi_t::calc()
{

  const int ndat = data->dim1();
  
  nchan = data->dim2();  // struct member

  bool jackknife = true;
  segshift = seglen / 2;
  const int epjack = eplen;

  // no jackknife is eplen was set to 0
  if ( eplen == 0 ) 
    {
      jackknife = false;
      eplen = ndat;	
    }

  // if no specified frequency bins, get highest frequency based on segment length
  //  (seglen = seglen in sample points)
  
  int maxfreqbin = floor( seglen/2 ) + 1;
  
  if ( freqbins.size() == 0 )
    add_freqbin( 1 , maxfreqbin );
  else
    maxfreqbin = max_freq();
    
  // store
  n_models = freqbins.size();

  int nepoch = floor( ndat / eplen );

  int nepochjack =  epjack > 0 ? floor(ndat/epjack) : 2;

  // para.segave=1;
  // para.subave=0;

  //[cs,nave]=data2cs_event(data,segleng,segshift,epleng,maxfreqbin,para);
  // fbin / MxM 
  std::vector<Data::Matrix<std::complex<double> > > cs = data2cs_event( data , maxfreqbin );
  
   // for (int i=0;i<cs.size();i++)
   //   {
   //     Data::Matrix<dcomp> r = cs[i];
   //     std::cout << "CS " << i << " " << r.print() << "\n";
   //   }

  const int nm = freqbins.size();

  const int nf = freqbins[0].size();
  // nb. this assumes that if there are multiple freqbins, they
  // all need to be the same size.

  // psall=zeros(nchan,nchan,nm);
  std::vector<Data::Matrix<double> > psall( nm );
  for (int i=0;i<nm;i++) psall[i].resize( nchan , nchan );

  // pssumall=zeros(nchan,nm);
  std::vector<Data::Vector<double> > pssumall( nm );
  for (int i=0;i<nm;i++) pssumall[i].resize( nchan );

  for (int ii=0; ii<nm; ii++)
    {
      const int nf = freqbins[ii].size();

      //psall[ii] = cs2ps(cs(:,:,freqbins(ii,:)));
      std::vector<Data::Matrix<std::complex<double> > > tt(nf);
      for (int f = 0 ; f < nf ; f++) tt[f] = cs[ freqbins[ii][f] - 1 ] ;

      psall[ii] = cs2ps( tt );
      // can redo the above to avoid copying

      // pssumall(:,ii)=sum(psall(:,:,ii),2);
      for (int i=0;i<nchan;i++)
	for (int j=0;j<nchan;j++)
	  pssumall[ii][i] += psall[ii](i,j);
      
    }

  
  //
  //psisum=squeeze(pssumall);  [ squeeze not needed ] 
  //  std::vector<Data::Vector<double> > psisum = pssumall;


  //
  // bootstrap
  //
  
  std::vector<Data::Matrix<std::complex<double> > > csall = cs;

  // nepochjack x nm x MxM
  //   psloc=zeros(nchan,nchan,nepochjack,nm);
  std::vector<std::vector<Data::Matrix<double> > > psloc( nepochjack );
  for (int b=0;b<nepochjack;b++)
    {
      psloc[b].resize( nm );
      for (int ii=0; ii<nm; ii++) psloc[b][ii].resize( nchan , nchan );
    }

  
  //   pssumloc=zeros(nchan,nepochjack,nm);  
  std::vector<std::vector<Data::Vector<double> > > pssumloc( nepochjack );
  for (int b=0;b<nepochjack;b++)
    {
      pssumloc[b].resize( nm );
      for (int ii=0; ii<nm; ii++) pssumloc[b][ii].resize( nchan );
    }
  
  //  if ( epjack > 0 )

  for (int b=0; b<nepochjack; b++)
    {

      //dataloc=data((b-1)*epjack+1:b*epjack,:);
      Data::Matrix<double> dataloc( epjack , nchan );

      // std::cout << "boot strap " << b << "\t"
      // 		<<  b*epjack << "\t"
      // 		<< (b+1)*epjack << "\n";

      int r = 0;
      for (int i = b*epjack ; i < (b+1)*epjack ; i++ )
	{
	  for (int j=0; j<nchan; j++)
	    dataloc(r,j) = (*data)(i,j);
	  ++r;
	}
      
      //      csloc=data2cs_event(dataloc,segleng,segshift,epleng,maxfreqbin,para);

      std::vector<Data::Matrix<std::complex<double> > > csloc = data2cs_event( &dataloc , 
									       maxfreqbin );
      
      
      //cs=(nepochjack*csall-csloc)/(nepochjack+1);
      dcomp nepochjack_cmp( nepochjack , 0 );
      dcomp nepochjack_cmp_p1( nepochjack+1 , 0 );
      std::vector<Data::Matrix<std::complex<double> > > cs = csall;
      for (int f=0;f<maxfreqbin; f++)
	for (int i=0;i<nchan;i++)
	  for (int j=0;j<nchan;j++)
	    {
	      //	      std::cout << "CSLOC = " << csloc[f](i,j) << "\n";
	      cs[f](i,j) = ( nepochjack_cmp * csall[f](i,j) - csloc[f](i,j) ) / nepochjack_cmp_p1 ; 
	      //	      std::cout << " >> cs = " << csall[f](i,j) << " " << csloc[f](i,j)  << " " <<  cs[f](i,j)  << "\n";
	    }
      
      for (int ii=0; ii<nm; ii++)
	{

	  //psloc(:,:,b,ii)=cs2ps(cs(:,:,freqbins(ii,:)));
	  std::vector<Data::Matrix<std::complex<double> > > tt(nf);
	  for (int f = 0 ; f < nf ; f++) tt[f] = cs[ freqbins[ii][f] - 1 ] ;
	  psloc[b][ii] = cs2ps( tt );
	  
	  //pssumloc(:,b,ii)=sum(psloc(:,:,b,ii),2);
	  for (int i=0;i<nchan;i++)
	    for (int j=0;j<nchan;j++)
	      pssumloc[b][ii][i] += psloc[b][ii](i,j);
	  
	  
	}
    }

  
  // store results:
  // vector of MxM PSI matrices [ freqbins / nm ] 
  
  psi = psall;

  // vector of vectors of total flux [ M ] 
  psi_sum = pssumall;

  // bootstreap SD deviations
  
  std_psi.clear();
  std_psi_sum.clear();
  
  for (int ii=0;ii<nm;ii++)
    {
      // PSI
      Data::Matrix<double> S( nchan , nchan ); 
      for (int i=0;i<nchan; i++)
	for (int j=0;j<nchan; j++)
	  {
	    std::vector<double> xx;
	    for (int b=0;b<nepochjack;b++) xx.push_back( psloc[b][ii](i,j) );
	    S(i,j) = MiscMath::sdev( xx ) * sqrt( nepochjack ) ;
	  }
      std_psi.push_back( S );

      // PSI SUM
      Data::Vector<double> V( nchan );
      for (int i=0;i<nchan; i++)
	{
	  std::vector<double> xx;
	  for (int b=0;b<nepochjack;b++)
	    {
	      xx.push_back( pssumloc[b][ii](i) );
	      //	      std::cout << " xx = " << pssumloc[b][ii](i) << "\n";
	    }
	      V(i) = MiscMath::sdev( xx ) * sqrt( nepochjack ) ;
	}
      //      std::cout << "V std = " << V.print() << "\n";
      std_psi_sum.push_back( V );
      
    }
        
  return;
  
}


// expects nf x [ M x M matrices ] 
// returns a M x M matrix
Data::Matrix<double> psi_t::cs2ps( const std::vector<Data::Matrix<std::complex<double> > > & cs )
{
  int df = 1;
    
  // number of frequencies (within this band)
  int nf = cs.size();
  
  // copy
  std::vector<Data::Matrix<std::complex<double> > > pp = cs;

  // for each 'f'
  for (int f=0; f<nf ; f++)
    {
      // divide each element by 
      //pp(:,:,f)=cs(:,:,f)./sqrt(diag(cs(:,:,f))*diag(cs(:,:,f))');
      // std::cout << "f = " << f << "\n";
      // std::cout << "cs size = " << cs.size() << "\n";
      
      const Data::Matrix<std::complex<double> > & cc = cs[f];
      Data::Matrix<std::complex<double> > & ppf = pp[f];

      for (int i=0;i<nchan;i++)
	for (int j=0;j<nchan;j++)
	  ppf(i,j) = cc(i,j) / sqrt( cc(i,i) * conj( cc(j,j) ) ) ;
      
      //pp(:,:,f)=cs(:,:,f)./sqrt(diag(cs(:,:,f))*diag(cs(:,:,f))');
    }

  Data::Matrix<double> ps( nchan , nchan );

  // sum of nf-1   f / f+df comparisons
  for (int f=1;f<nf;f++)
    {
      const Data::Matrix<std::complex<double> > & pp0 = pp[f-1];
      const Data::Matrix<std::complex<double> > & pp1 = pp[f];
	    
      for (int i=0;i<nchan;i++)
	for (int j=0;j<nchan;j++)
	  ps(i,j) += std::imag( conj( pp0(i,j) ) * pp1(i,j) );
      
      // ps=sum( imag(  conj( pp(:,:,1:end-df) ) .* pp(:,:,1+df:end)  ) ,3);
    }
  
  return ps;
}



// calculates cross-spectra from data for event-related measurement

//  requires:
// data: ndat times nchan matrix each colum is the time-series in one channel;
// segleng: length of each segment in bins, e.g. segleng=1000;  
// segshift: numer of bins by which neighboring segments are shifted;
//           e.g. segshift=segleng/2 makes overlapping segments
// epleng: length of each epoch
// maxfreqbin: max frequency in bins

// para: optional structure:
//       para.segave=0     -> no averaging across segments 
//       para.segave neq 0 -> averaging across segments (default is 0)% \
//       para.subave = 1 subtracts the average across epochs,  
//       para.subave ~= 1 -> no subtraction (default is 1) 
//       IMPORTANT: if you just one epoch (e.g. for continuous data)
//         set para.subave=0          
//       -> averaging across segments (default is 0)

//       para.proj must be a set of vector in channel space,  
//       if it exists then the output raw contains the single trial 
//       Fourier-transform in that channel   
     
// output: 
// cs: nchan by chan by maxfreqbin by nseg tensor cs(:,:,f,i) contains 
//     the cross-spectrum at frequency f and segment i
// nave: number of averages

std::vector<Data::Matrix<std::complex<double> > > psi_t::data2cs_event( const Data::Matrix<double> * mydata , 
									int maxfreqbin )
{


  

  // from main: 
  //  para.segave=1;   YES : AVERAGE ACROSS SEGMENTS
  //  para.subave=0;   NO  : SUBTRACTION OF AVERAGE ACROSS EPOCHS
  
  const bool segave = true;   
  const bool subave = false; 

  // ensure maxfreqbin is valid
  if ( maxfreqbin > floor(seglen/2)+1 )
    maxfreqbin =  floor(seglen/2)+1;
  
  const bool mydetrend = false; // no detrending  
  
  const int ndat = mydata->dim1();

  const int nchan = mydata->dim2();

  const int nep = floor(ndat/eplen);

  //  std::cout << "eplen seglen segshift " << eplen << " " << seglen << " " << segshift << "\n";
  
  const int nseg = floor((eplen-seglen)/segshift)+1; //total number of segments

  std::vector<Data::Matrix<std::complex<double> > > cs( maxfreqbin );
  for (int f=0;f<maxfreqbin;f++) cs[f].resize( nchan , nchan );

  // not needed(?) for now, if no epoch subtraction
  // av=zeros(nchan,maxfreqbin);
  
  // get Hanning window of seglen, in mywindow
  std::vector<double> window = MiscMath::hanning_window( seglen );

  // for (int i=0;i<window.size(); i++)
  //   std::cout << "window " << i << " " << window[i] << "\n";
  
  // for each epoch
  int nave = 0;

  // nb. ML useds j in 1-based encoding
  for (int j=0; j<nep; j++) 
    {
      // get epoch slice of 'mydata'  -- horrible recopying but keeping for now
      //     dataep=data((j-1)*epleng+1:j*epleng,:);
      // rows : j * epleng
      //  to  : (j+1)*epleng-1 
      Data::Matrix<double> dataep( eplen , nchan );
      int start = j*eplen;
      int stop =  (j+1)*eplen-1;
      //      std::cout << "eplen " << eplen << " from start , stop = " << start << " " << stop << "\n";
      int r = 0;
      for (int i=start;i<=stop;i++)
	{
	  for (int j=0;j<nchan;j++)
	    dataep(r,j) = (*mydata)(i,j) ;
	  ++r;
	}

      // average over each segment 
      for (int i=0; i<nseg; i++) 
	{

	  // get data
	  //	 dataloc=dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
	  int start = i*segshift;
	  int stop =  i*segshift+seglen-1;
	  //std::cout << "seglen " << seglen << " from start , stop = " << start << " " << stop << "\n";
	  
	  Data::Matrix<double> dataloc( seglen , nchan );
	  int r = 0;
	  for (int i=start;i<=stop;i++)
	    {
	      for (int j=0;j<nchan;j++)
		dataloc(r,j) = dataep(i,j) * window[r];
	      ++r;
	    }

	 // get FFT (w/ window on data)
	 // datalocfft=fft(dataloc.*mywindow);

	  // have to do channels separately... windowing done above
	  std::vector<std::vector<dcomp> > datalocfft;
	  for (int j=0; j<nchan; j++)
	    {

	      int my_Fs = 400; // arbitrary
	      int index_length = seglen;
	      int index_start = 0;

	      FFT fftseg( index_length , index_length , my_Fs , FFT_FORWARD , WINDOW_NONE );

	      const std::vector<double> * pd = dataloc.col(j).data_pointer();
	      //fftseg.apply( &(x[index_start]) , index_length );
	      fftseg.apply( &((*pd)[0]) , index_length );

	      // Extract the raw transform
	      datalocfft.push_back( fftseg.transform() );
	      
	      
	    }

	  // std::cout << "FFT\n";
	  // for (int f=0; f<maxfreqbin; f++)
	  //   for (int j=0; j<nchan;j++)
	  //     std::cout << "FFT\t" << f+1 << "\t" << j+1 << "\t" << datalocfft[j][f] << "\n";

	  
	 // OR.. optionally, detrend segment
	 //datalocfft=fft(detrend(dataloc,0).*mywindow);

	 // for each frequency (from 1 Hz to MX)
         for (int f=0; f<maxfreqbin; f++) 
	   {

	     // cs(:,:,f)= cs(:,:,f)+conj(datalocfft(f,:)'*datalocfft(f,:));
	     for (int i=0; i<nchan; i++)
	       for (int j=0; j<nchan; j++)
	      	 cs[f](i,j) += conj( conj(datalocfft[i][f] ) * datalocfft[j][f] ) ;
	     

	     // av(:,f)=av(:,f)+conj(datalocfft(f,:)');
	   }
	 
       } // next segment

     nave=nave+1;

   } // next epoch


 // averageing 
 nave = nave * nseg;  
 
 // std::cout << " nave = " << nave << "\n";

 for (int f=0;f<maxfreqbin;f++)   
   {
     for (int i=0;i<nchan;i++)
       {
	 //av[f][i] /= (double)nave; // not needed if not epoch averaging
	 for (int j=0;j<nchan;j++)
	   cs[f](i,j) /= (double)nave;
       }
   }

 return cs;

}

