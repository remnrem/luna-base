
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

#include "simul.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include "eval.h"
#include "fftw/fftwrap.h"
#include "miscmath/crandom.h"
#include "dsp/spline.h"

#include "db/db.h"

extern writer_t writer;


void dsptools::simul( edf_t & edf , param_t & param )
{

  //
  // Update or create a signal
  //
 
  std::string siglab = param.requires( "sig" );

  const bool channel_mode = siglab != "*" ;

  if ( ! channel_mode ) Helper::halt( "need to specify a single channel" ) ;
  
  const bool update_existing_channel = edf.header.has_signal( siglab );

  // add to an existing signal (versus replace)?

  const bool add_to_existing = param.has( "add" ) ;

  if ( add_to_existing && ! update_existing_channel )
    Helper::halt( "specified 'add' to modify an existing signal, but it does not exist" );
  
  //
  // baseline signal from a specified PSD, either:
  //
  //  file    : from a file
  //  alpha   : 1/f^a
  //  freq    : freq=10,20,30  amp=1,1,1

  //
  // Also allow transients (pulses)
  //
  //  here we simply mask out everything *not* in a pulse, i.e. at the end of creating it
  //  this makes sense with the 'add' command, i.e. to inject pulses intp a real signal / simulated
  //
  
  const bool pulses = param.has( "pulses" );

  
  // 
  // peaks specification
  // 

  const bool simple = param.has( "frq" );
  
  
  //
  // M/f^a slope
  //

  const bool functional = param.has( "alpha" );

  const bool from_file = param.has( "file" );

  if ( from_file && ( functional || simple ) )
    Helper::halt( "cannot specify alpha/frq as well as file" );


  //
  // impulses (i.e. for impulse and step functions)
  //

  const bool impulses = param.has( "impulse" );
  
  //
  // sample rate & frequency resolution
  //

  int fs = 0;

  if ( update_existing_channel )
    {

      const int slot = edf.header.signal( siglab );
      
      if ( edf.header.is_annotation_channel( slot ) )
	Helper::halt( "cannot modify an EDF Annotation channel" );
      
      fs = edf.header.sampling_freq( slot );
      
      if ( param.has( "sr" ) )
	{	  
	  int fs0 = param.requires_int( "sr" );	  
	  if ( fs0 != fs )
	    Helper::halt( "cannot specify a different 'sr' if updating an existing channel (which has Fs = " + Helper::int2str( fs ) );
	}

    }
  else // we requre one to be specified 
    {
      fs = param.requires_int( "sr" );
      logger << "  using sample rate " << fs << "\n";
    }
  
  const double fmax = fs / 2.0;      
  
  const int n = edf.header.record_duration * edf.header.nr * fs ;       
  
  const int m = floor( n / 2 ) + 1;

  const double df = fmax / (double)(m-1) ; 

  
  //
  // Generate PSD to be simulated
  //

  
  std::vector<double> frqs , psds;

  if ( ! from_file )
    {
      frqs.resize( m , 0 );
      psds.resize( m , 0 );
    }
  
  //
  // Read from a file?
  //

  if ( from_file )
    {

      // assume RAW PSD as per Luna output format:  cols F and PSD
      // *IF* we see a negative PSD value, then assume dB scaled
      
      const std::string psd_file = Helper::expand( param.requires( "file" ) );

      if ( ! Helper::fileExists( psd_file ) )
	Helper::halt( "cannot read PSD from " + psd_file );
      
      int f_slot = -1 , psd_slot = -1;
      int header_size;
      bool has_neg = false;
      
      std::ifstream IN1( psd_file.c_str() , std::ios::in );
      while ( ! IN1.eof() )
	{
	  std::string line;
	  Helper::safe_getline( IN1 , line );
	  if ( IN1.eof() ) break;
	  if ( line == "" ) continue;
	  std::vector<std::string> tok = Helper::parse( line , "\t" );
	  
	  // look for headers?
	  if ( f_slot == -1 )
	    {
	      for (int i=0; i<tok.size(); i++)
		{
		  if ( tok[i] == "F" ) f_slot = i;
		  else if ( tok[i] == "PSD" ) psd_slot = i;
		}
	      if ( f_slot == -1 || psd_slot == -1 )
		Helper::halt( "did not find F and PSD in the header" );
	      header_size = tok.size();
	    }
	  else
	    {
	      double f = 0;
	      double psd = 0;
	      
	      if ( tok.size() != header_size )
		Helper::halt( "bad number of columns: " + line );
	      
	      if ( ! Helper::str2dbl( tok[ f_slot ] , &f ) )
		Helper::halt( "bad numeric value for F = " + tok[ f_slot ] );
	      
	      if ( ! Helper::str2dbl( tok[ psd_slot ] , &psd ) )
		Helper::halt( "bad numeric value for PSD = " + tok[ psd_slot ] );

	      // negatives?  assume dB scale if so (allow for rounding)
	      if ( psd < 0 ) has_neg = true;
	      
	      // store
	      frqs.push_back( f );
	      psds.push_back( psd );
	    }
	  
	}
      IN1.close();

      //
      // Set on raw scale if needed, assuming 10log10(X) 
      //

      if ( has_neg )
	{
	  for (int i=0; i<psds.size(); i++)
	    psds[i] = pow( 10.0 , psds[i] / (double)10 );
	}
      
      //
      // Resample to given length (cubic spline interpolation)
      //
      
      if ( frqs.size() != m )
	{
	  std::vector<double> f0 = frqs;
	  std::vector<double> p0 = psds;
	  
	  tk::spline spline;
	  spline.set_points( f0 , p0 );
	  	  
	  frqs.resize( m );
	  psds.resize( m );

	  for (int i=0; i<m; i++)
	    {
	      double f = i * df;
	      frqs[ i ] = f ;
	      psds[ i ] = spline( f );
	      if ( psds[ i ] < 0 ) psds[i] = 0;
	      //std::cout << " up sampled " << frqs[i] << "\t" << psds[i] << "\n";
	    }
	}

    }

  
  //
  // If functional form 
  //

  if ( functional )
    {
      const double alpha = param.requires_dbl( "alpha" );
      const double intercept = param.requires_dbl( "intercept" );
      
      for (int i=0; i<m; i++)
	{
	  frqs[i] = i * df;
	  psds[i] = i == 0 ? 0 : intercept / pow( frqs[i] , alpha );
	  //std::cout << " PSD " << i << "\t" << frqs[i] << "\t" << psds[i] << "\n";
	}
    }
  

  //
  // Simple peaks
  //

  if ( simple )
    {
	
      std::vector<double> f = param.dblvector( "frq" );
      std::vector<double> a = param.dblvector( "psd" );
      
      // also, assume a Gaussian distribution around each frequency, going to zero?
      const double pwidth = param.has( "w" ) ? param.requires_dbl( "w" ) : 0 ;
      if ( pwidth < 0 ) Helper::halt( "cannot have negative w");
      
      if ( f.size() != a.size() || f.size() == 0 ) Helper::halt( "bad frq=X,Y,Z psd=X,Y,Z specification" );
            
      // fixed points (no variance)
      
      if ( ! param.has( "w" ) )
	{
	  
	  for (int i=0; i<m; i++)
	    {
	      double frq = i * df;
	      frqs[i] = frq;
	      for (int j=0; j<f.size(); j++)
		{
		  if ( fabs( f[j] - frq ) < 1e-8 )
		    psds[i] = a[j];		      
		}
	    }
	}
      else
	{
	  for (int j=0; j<f.size(); j++)
	    {
	      std::vector<double> w( frqs.size() , 0 );
	      double mx = 0;
	      for (int i=0; i<m; i++)
		{
		  const double frq = i * df;
		  frqs[i] = frq;
		  w[i] = Statistics::normden( frq , f[j] , pwidth );
		  if ( w[i] > mx ) mx = w[i];
		  if ( w[i] < 0 ) w[i] = 0;
		}
	      
	      // nb. add to existing (i.e. if overlapping peaks)
	      const double y = a[j] / mx;
	      for (int i=0; i<m; i++)
		psds[i] += w[i] * y;	      
	    }
	}
    }
  
  
  //
  // Nothing specified?
  //

  if ( frqs.size() == 0 ) Helper::halt( "no PSD specified" );


  //
  // Verbose report?
  //


  if ( param.has( "verbose" ) )
    {

      for (int i=0; i < psds.size(); i++)
	{
	  writer.level( frqs[i] , globals::freq_strat );

	  if ( frqs[i] > 0 )
	    writer.value( "LF" , log( frqs[i] ) );
	  
	  writer.value( "P" , psds[i] );
	  
	  if ( psds[i] >= 0 ) 
	    {
	      writer.value( "LP" , log( psds[i] ) );
	      writer.value( "DB" , 10 * log10( psds[i] ) );
	    }
	  
	}
      writer.unlevel( globals::freq_strat );
    }
  

  //
  // Resample PSD at the required number of points
  //
  
  // take PSD, rescale --> amplitude
  //  set phases as random (no dependency between frequencies); fixed DC & Nyquist
  //  make complex variable Z,  then apply inverse FFT to get the final signal

  //
  // power --> amplitude
  //

  std::vector<double> amps = psds;
  for (int i=0; i<amps.size(); i++)
    amps[i] = sqrt( 0.5 * psds[i] * n * fs );

  
  // randomise phases (uniform between 0 and 2PI)
  // skip DC & Nyquist
  
  std::vector<double> phases( amps.size() );
  for (int i=1; i<phases.size()-1; i++) 
    phases[i] = 2 * M_PI * CRandom::rand();
  
  // make freq. domain signal Z = amp x e^(i phi), where i = sqrt(-1)
  
  std::vector<dcomp> z( amps.size() );
  for (int i=0; i<z.size(); i++)
    z[i] = amps[i] * exp( dcomp( 0 , phases[i] ) ); 
  
  // mirror negative freqs (complex conjugate)
  const bool even = n % 2 == 0 ;
  for (int p = even ? amps.size() - 2 : amps.size() - 1 ; p != 0;  p-- )
    z.push_back( std::conj( z[p] ) );
  
  if ( z.size() != n )
    Helper::halt( "internal error in dsptools::simul()" );
  
  // inverse FFT to get time-domain signal     
  real_iFFT ifft( n , n , fs );
  ifft.apply( z );
  std::vector<double> rdat = ifft.inverse();

  // check back: do FFT on this result

  real_FFT fft( n , n , fs );
  fft.apply( rdat );
  std::vector<dcomp> ft = fft.transform();
  
  // for (int j=0;j<fft.frq.size();j++)
  //   std::cout << " RET\t" << fft.frq[j] << "\t" << ft[j] << "\t z = " << z[j] << " PSD " << psds[j] << "  vs  " << fft.X[j] << "\n";
  

  //
  // Simulate pulses (i.e. blank out everything that is *not* a signal
  //  for now, just have on/off pulses, applied in the time domain
  //   i.e. not specific to a given frequency
  //

  if ( pulses )
    {
      // pulses=N,T
      //  N = number of pulses
      //  T = duration (in seconds)
      // by default, pulses do not overlap

      std::vector<double> popt = param.dblvector( "pulses" );
      if ( popt.size() != 2 )
	Helper::halt( "expecting pulses=N,T" );

      const int pn = popt[0];
      const double pt = popt[1];

      if ( pn < 0 ) Helper::halt( "cannot specify a -ve number of pulses" );
      if ( pt < 0 ) Helper::halt( "cannot specify a -ve pulse duration" );

      // rough sanity check.... cannot have most of signal as a 'pulse'
      const double totsec = rdat.size() / (double)fs;
      if ( ( pn * pt ) / totsec > 0.8 )
	Helper::halt( "cannot specify >80% of the signal as expected to be a pulse" );

      logger << "  applying " << pn << " pulses, each of " << pt << " seconds\n";
      
      const int mxtry = 1000;
      const int plen  = fs * pt;
      
      std::vector<bool> pmask( n , true );

      for (int i=0; i<pn; i++)
	{

	  int try1 = 0;

	  while ( 1 ) 
	    {
	      ++try1;
	      
	      if ( try1 > mxtry )
		Helper::halt( "could not apply all pulses (w/out overlap)... reduce pulse number or duration" );
	      
	      int p0 = CRandom::rand( n );
	      if ( ! pmask[p0] ) continue;
		
	      bool okay = true;
	      int p = p0;
	      for (int j=0; j<plen; j++)
		{
		  ++p;
		  if ( p == n ) { okay = false; break; }
		  if ( ! pmask[p] ) { okay = false; break; }		  
		}
	      //std::cout << " okay = " << okay << "\n";

	      // failed. try again
	      if ( ! okay ) continue;
	      
	      // std::cout << " Setting... p = " << p << " " << try1 << "\n";

	      // else, mark this pulse
	      for (int j=0; j<plen; j++)
		pmask[ p0 + j ] = false;
	      
	      // place next pulse 
	      break;
	    }
	}

      // now blank out non-pulses
      for (int i=0; i<n; i++)
	if ( pmask[i] ) rdat[i] = 0;
      
    }
  
  
  //
  // Create / update signal
  //
  
  if ( update_existing_channel )
    {
      
      const int slot = edf.header.signal( siglab );
      
      if ( edf.header.is_annotation_channel( slot ) )
	Helper::halt( "cannot modify an EDF Annotation channel" );

      // whether adding, or completely changing, we need to make sure we have 'read'
      // this channel (i.e. to appropriately populate the edf_record_t stores.  otherwise,
      // update_signal() will try to modify something that has not been read.   i.e.
      // previously, we always *ASSumed* that if we were to update something, we would
      // have read it first, but this is not the case if add_to_existing == false 
      // (unless we specifically read it)
      
      slice_t slice( edf , slot , edf.timeline.wholetrace() );
      
      const std::vector<double> * d = slice.pdata();
      
      if ( d->size() != rdat.size() )
	Helper::halt( "internal error in simul()" );
      
      if ( add_to_existing )
	{
	  // add existing to simulated signal:
	  for (int i=0; i<d->size(); i++)
	    rdat[i] += (*d)[i];	  
	}

      // now update the channel
      logger << "  updating " << siglab << "...\n";
      edf.update_signal( edf.header.signal( siglab ) , &rdat );

    }
  else
    {      
      logger << "  creating new channel " << siglab << "...\n";
      edf.add_signal( siglab , fs , rdat );
    }
  

  
}

