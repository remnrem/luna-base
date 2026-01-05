
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

#include "dsp/wrappers.h"
#include "edf/edf.h"
#include "param.h"

#include "edf/slice.h"
#include "cwt/cwt.h"
#include "dsp/hilbert.h"
#include "fftw/fftwrap.h"
#include "dynamics/qdynam.h"
#include "defs/defs.h"

#include "db/db.h"
#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger;

#include <vector>
  
void dsptools::cwt( edf_t & edf , param_t & param )
{

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  
  const int ns = signals.size();

  std::vector<double> fc;
  if ( param.has( "fc-inc" ) )
    {
      std::vector<double> f = param.dblvector( "fc-inc" );
      if ( f.size() != 3 ) Helper::halt( "expecting fc-inc=start,stop,inc" );
      for (double ff = f[0] ; ff <= f[1] ; ff += f[2] ) fc.push_back(ff);
    }
  else
    {
      if ( ! param.has( "fc" ) ) Helper::halt( "no fc for CWT" );
      
      fc = param.dblvector( "fc" );
    }

  bool alt_spec = param.has( "fwhm" );

  double fwhm = alt_spec ? param.requires_dbl( "fwhm" ) : 0 ;
  
  int num_cycles = alt_spec ? 0 : param.requires_int( "cycles" );

  double timelength = alt_spec ? ( param.has( "len" ) ? param.requires_dbl( "len" ) : 20 ) : 0 ; 
  
  bool return_phase = param.has( "phase" );
  
  bool wrapped_wavelet = param.has( "wrapped" );

  std::string tag = param.has( "tag" ) ? "_" + param.value( "tag" ) : "" ; 

  for (int s=0;s<ns;s++)
    {

      if ( edf.header.is_annotation_channel( signals(s) ) ) 
	continue;
      
      const int Fs = edf.header.sampling_freq( signals(s) );
      
      interval_t interval = edf.timeline.wholetrace();
      
      slice_t slice( edf , signals(s) , interval );

      const std::vector<double> * d = slice.pdata();

      for (int fi=0; fi<fc.size(); fi++)
	{
    
	  std::vector<double> mag , phase;
	  
	  if ( alt_spec )
	    alt_run_cwt( *d , Fs , fc[fi] , fwhm , timelength , wrapped_wavelet , &mag , return_phase ? &phase : NULL );
	  else
	    run_cwt( *d , Fs , fc[fi] , num_cycles , &mag , return_phase ? &phase : NULL );
      
	  std::string new_mag_label = signals.label(s) + tag + "_cwt_mag";
	  std::string new_phase_label = signals.label(s) + tag + "_cwt_ph";

	  if ( fc.size() > 1 ) 
	    {
	      new_mag_label += "_" + Helper::int2str( fi+1 );
	      new_phase_label += "_" + Helper::int2str( fi+1 );
	    }
      
	  logger << "  CWT, Fc = " << fc[fi] << ", for " << signals.label(s) << " --> " << new_mag_label ;      
	  
	  if ( return_phase ) 
	    logger << ", " << new_phase_label ;
	  logger << "\n";
	  
	  edf.add_signal( new_mag_label , Fs , mag );
	  
	  if ( return_phase ) 
	    edf.add_signal( new_phase_label , Fs , phase );
	  
	} // next Fc
      
    } // next signal
  
}


void dsptools::hilbert( edf_t & edf , param_t & param )
{

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  
  const int ns = signals.size();
  
  bool use_kaiser = param.has( "tw" );

  bool use_file = param.has( "file" );

  bool use_fixed = param.has( "order" );

  bool no_filter = ! ( use_kaiser || use_file || use_fixed );

  if ( ! no_filter ) 
    if ( ! ( param.has( "f" ) || param.has( "bandpass" ) ) )
      Helper::halt( "requires 'f' or 'bandpass'" );
  
  std::vector<double> frqs = param.dblvector( param.has( "f" ) ? "f" : "bandpass" );

  double ripple = use_kaiser ? param.requires_dbl( "ripple" ) : 0 ;

  double tw = use_kaiser ? param.requires_dbl( "tw" ) : 0 ;

  int order = use_fixed ? param.requires_int( "order" ) : 0;

  fir_t::windowType window = fir_t::HAMMING;
  if ( param.has( "rectangular" ) ) window = fir_t::RECTANGULAR;
  else if ( param.has( "bartlett" ) ) window = fir_t::BARTLETT;
  else if ( param.has( "hann" ) ) window = fir_t::HANN;
  else if ( param.has( "blackman" ) ) window = fir_t::BLACKMAN;

  std::string fir_file = use_file ? param.value( "file" ) : "";
  
  bool return_phase = param.has( "phase" ) || param.has( "angle" );

  bool return_angle = param.has( "angle" );

  bool return_ifrq = param.has( "ifrq" );
  
  std::string tag = param.has( "tag" ) ? "_" + param.value( "tag" ) : "" ; 

  for (int s=0;s<ns;s++)
    {

      if ( edf.header.is_annotation_channel( signals(s) ) ) 
	continue;
      
      const int Fs = edf.header.sampling_freq( signals(s) );
      
      interval_t interval = edf.timeline.wholetrace();
      
      slice_t slice( edf , signals(s) , interval );

      const std::vector<double> * d = slice.pdata();
      
      std::vector<double> mag , phase, ifrq;

      if ( use_kaiser ) 
	run_hilbert( *d , Fs , frqs[0] , frqs[1] , ripple , tw ,
		     &mag ,
		     return_phase ? &phase : NULL ,
		     return_angle ? &phase : NULL ,
		     return_ifrq ? &ifrq : NULL );
      else if ( use_fixed )
	run_hilbert( *d , Fs , frqs[0] , frqs[1] , order , window ,
		     &mag ,
		     return_phase ? &phase : NULL ,
		     return_angle ? &phase : NULL ,
		     return_ifrq ? &ifrq : NULL );
      else if ( use_file )
	run_hilbert( *d , Fs , fir_file , 
		     &mag ,
		     return_phase ? &phase : NULL ,
		     return_angle ? &phase : NULL ,
		     return_ifrq ? &ifrq : NULL );
      else
	run_hilbert( *d , Fs , 
		     &mag ,
		     return_phase ? &phase : NULL ,
		     return_angle ? &phase : NULL ,
		     return_ifrq ? &ifrq : NULL );

      //
      // labels for new EDF channel(s)
      //

      std::string new_mag_label = signals.label(s) + tag + "_ht_mag";
      std::string new_phase_label = signals.label(s) + tag + "_ht_ph";
      std::string new_angle_label = signals.label(s) + tag + "_ht_ang";
      std::string new_ifrq_label = signals.label(s) + tag + "_ht_ifrq";

      logger << " Hilbert transform for " << signals.label(s) << " --> " << new_mag_label ;      

      if ( return_phase && ! return_angle ) 
	logger << ", " << new_phase_label ;
      
      if ( return_angle )
	logger << ", " << new_angle_label ;

     if ( return_ifrq ) 
	logger << ", " << new_ifrq_label ;

      logger << "\n";

      edf.add_signal( new_mag_label , Fs , mag );
      
      if ( return_phase && ! return_angle ) 
	edf.add_signal( new_phase_label , Fs , phase );
      
      if ( return_angle )
	edf.add_signal( new_angle_label , Fs , phase );

      if ( return_ifrq ) 
	{
	  // this returns n-1 estimates, i.e. based on the derivative of the phase
	  // add 0 to the end as a null marker, so it can be placed back in the EDF
	  ifrq.push_back(0);
	  edf.add_signal( new_ifrq_label , Fs , ifrq );
	}
    }  

}

  
void dsptools::alt_run_cwt( const std::vector<double> & data ,
			    const int Fs, 
			    const double fc ,
			    const double FWHM ,
			    const double tlen , 
			    const bool wrapped , 
			    std::vector<double> * mag , 
			    std::vector<double> * phase )
{
  
  CWT cwt;
  
  cwt.set_sampling_rate( Fs );
  
  cwt.set_timeframe( 50.0 / tlen );

  cwt.alt_add_wavelet( fc , FWHM , tlen );
  
  cwt.store_real_imag_vectors( true );

  cwt.load( &data );

  if ( wrapped ) 
    cwt.run_wrapped();
  else
    cwt.run();
  
  *mag = cwt.results(0);
  
  if ( phase != NULL ) 
    *phase = cwt.phase(0);
  
}



void dsptools::run_cwt( const std::vector<double> & data ,
			const int Fs, 
			const double fc ,
			const int num_cycles , 
			std::vector<double> * mag , 
			std::vector<double> * phase )
{
  
  CWT cwt;
  
  cwt.set_sampling_rate( Fs );
  
  cwt.add_wavelet( fc , num_cycles ); 
  
  cwt.load( &data );
  
  cwt.run();
  
  *mag = cwt.results(0);
  
  if ( phase != NULL ) 
    *phase = cwt.phase(0);
  
}




void dsptools::run_hilbert( const std::vector<double> & data , const int Fs , 
			    std::vector<double> * mag , 
			    std::vector<double> * phase ,
			    std::vector<double> * angle ,
			    std::vector<double> * ifrq )
{
  
  // straight Hilbert , no filter
  hilbert_t hilbert( data );
  
  if ( mag != NULL ) *mag = *(hilbert.magnitude());
  
  if ( phase != NULL ) *phase = *(hilbert.phase());
  
  if ( angle != NULL )
    {
      *angle = *phase;
      // convert to degrees with 0 as pos-to-neg crossing
      for (int i=0;i<angle->size();i++) (*angle)[i] = MiscMath::as_angle_0_pos2neg( (*angle)[i] );	  
    }
  
  if ( ifrq != NULL )  *ifrq = hilbert.instantaneous_frequency( Fs );    
  
}



// filter-Hilbert, Kaiser window
void dsptools::run_hilbert( const std::vector<double> & data , const int Fs , 
			    const double flwr , const double fupr , const double ripple , const double tw , 
			    std::vector<double> * mag , 
			    std::vector<double> * phase ,
			    std::vector<double> * angle ,
			    std::vector<double> * ifrq )
{
  
  // filter-Hilbert 
  hilbert_t hilbert( data , Fs , flwr , fupr , ripple , tw );

  if ( mag != NULL ) *mag = *(hilbert.magnitude());

  if ( phase != NULL ) *phase = *(hilbert.phase());
  
  if ( angle != NULL )
    {
      *angle = *phase;
      // convert to degrees with 0 as pos-to-neg crossing
      for (int i=0;i<angle->size();i++) (*angle)[i] = MiscMath::as_angle_0_pos2neg( (*angle)[i] );	  
    }
  
  if ( ifrq != NULL )  *ifrq = hilbert.instantaneous_frequency( Fs );    

}

// filter-Hilbert, from file
void dsptools::run_hilbert( const std::vector<double> & data , const int Fs , 
			    const std::string & fir_file , 
			    std::vector<double> * mag , 
			    std::vector<double> * phase ,
			    std::vector<double> * angle ,
			    std::vector<double> * ifrq )
{
  
  hilbert_t hilbert( data , Fs , fir_file );

  if ( mag != NULL ) *mag = *(hilbert.magnitude());

  if ( phase != NULL ) *phase = *(hilbert.phase());
  
  if ( angle != NULL )
    {
      *angle = *phase;
      // convert to degrees with 0 as pos-to-neg crossing
      for (int i=0;i<angle->size();i++) (*angle)[i] = MiscMath::as_angle_0_pos2neg( (*angle)[i] );	  
    }
  
  if ( ifrq != NULL )  *ifrq = hilbert.instantaneous_frequency( Fs );    

}


// filter-Hilbert, fixed order
void dsptools::run_hilbert( const std::vector<double> & data , const int Fs , 
			    const double flwr , const double fupr , const int order , const fir_t::windowType window , 
			    std::vector<double> * mag , 
			    std::vector<double> * phase ,
			    std::vector<double> * angle ,
			    std::vector<double> * ifrq )
{
  
  // filter-Hilbert 
  hilbert_t hilbert( data , Fs , flwr , fupr , order , window );

  if ( mag != NULL ) *mag = *(hilbert.magnitude());

  if ( phase != NULL ) *phase = *(hilbert.phase());
  
  if ( angle != NULL )
    {
      *angle = *phase;
      // convert to degrees with 0 as pos-to-neg crossing
      for (int i=0;i<angle->size();i++) (*angle)[i] = MiscMath::as_angle_0_pos2neg( (*angle)[i] );	  
    }
  
  if ( ifrq != NULL )  *ifrq = hilbert.instantaneous_frequency( Fs );    

}

//
// Welch
//

bool dsptools::welch( const std::vector<double> & x ,
		      double fs,
		      std::vector<double>& freqs_hz,
		      std::vector<double>& psd,		      
		      const double segment_sec ,
		      const double overlap_sec ,
		      double upr )		      
{

  // const double segment_sec  = 4.0;
  // const double overlap_sec = 2.0;
  const int total_points = x.size();
  const int segment_points = segment_sec * fs;
  const int noverlap_points  = overlap_sec * fs;

  // too short
  if ( total_points < fs * segment_sec ) return false;
  
  // return up to this Hz, or Nyquist if not defined
  if ( upr < 0 ) upr = fs / 2.0;
  
  // implied number of segments
  int noverlap_segments = floor( ( total_points - noverlap_points) 
				 / (double)( segment_points - noverlap_points ) );
  
  // Welch uses mean (not median) over segments
  const bool use_median = false; 
  
  window_function_t window_function = WINDOW_TUKEY50;

  PWELCH pwelch( x , 
		 fs , 
		 segment_sec , 
		 noverlap_segments , 
		 window_function ,
		 use_median ); 
  
  freqs_hz.clear();
  psd.clear();
  
  for ( int i = 0 ; i < pwelch.freq.size(); i++) 
    {
      if ( pwelch.freq[i] > upr ) break;
      freqs_hz.push_back( pwelch.freq[i] );
      psd.push_back( pwelch.psd[i] );
    }

  return true;
}



//
// FFT
//

void dsptools::fft( edf_t & edf , param_t & param )
{

  //
  // whole signal FFT
  //

  // show real/imaginary frequency-domain values?
  const bool verbose = param.has( "verbose" );
  
  //
  // iterate over signals
  //

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  
  const int ns = signals.size();

  logger << "  calculating DFT:";

  for (int s=0; s<ns; s++)
    {
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) 
	continue;
      
      writer.level( signals.label(s) , globals::signal_strat );

      logger << " " << signals.label(s) ;
      
      const int Fs = edf.header.sampling_freq( signals(s) );
      
      interval_t interval = edf.timeline.wholetrace();
      
      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();
      
      run_fft( *d , Fs , verbose );
      
      writer.unlevel( globals::signal_strat );
    }
  logger << "\n";
}


std::vector<double> dsptools::readcin()
{
  std::vector<double> x;

  int cnt = 0;
  while ( ! std::cin.eof() )
    {
      double xx;
      std::cin >> xx;
      if ( std::cin.bad() ) Helper::halt( "bad input" );
      if ( std::cin.eof() ) break;	  
      x.push_back( xx );	  
      if ( ++cnt % 100000  == 0 )
	logger << " line " << cnt << "\n";
    }
  logger << x.size() << " values read\n";
  return x;
}


void dsptools::cmdline_fft( param_t & param )
{
  
  std::vector<double> x = dsptools::readcin();

  const int sr = param.has( "sr" ) ? param.requires_int( "sr" ) : 100 ; 

  logger << "  setting sr = " << sr << "\n";
  
  const bool verbose = param.has( "verbose" );
  
  dsptools::run_fft( x , sr , verbose );
  
}

void dsptools::run_fft( const std::vector<double> & x , const int Fs , const bool verbose )
{
  
  int index_length = x.size();
  
  FFT fftseg( index_length , index_length , Fs , FFT_FORWARD , WINDOW_NONE );
  
  fftseg.apply( &(x[0]) , index_length );
      
  // Extract the raw transform
  std::vector<std::complex<double> > t = fftseg.transform();
  
  // Extract the raw transform scaled by 1/n
  std::vector<std::complex<double> > t2 = fftseg.scaled_transform();
  
  int my_N = fftseg.cutoff;      
  
  for (int f=0;f<my_N;f++)
    {
      writer.level( fftseg.frq[f] , globals::freq_strat );
      
      if ( verbose )
	{
	  writer.value( "RE" , std::real( t[f] ) );
	  writer.value( "IM" , std::imag( t[f] ) );
	  writer.value( "UNNORM_AMP" , fftseg.mag[f] );
	  writer.value( "NORM_AMP" , ( f == 0 ? 1 : 2 ) * fftseg.mag[f] / (double)index_length );
	}
      
      writer.value( "PSD" , fftseg.X[f] );
      
      if ( fftseg.X[f] > 0 ) 
	writer.value( "DB" , log10( fftseg.X[f] )  );
    }
  writer.unlevel( globals::freq_strat );
    
}



//
// Otsu
//
  
void dsptools::otsu( edf_t & edf , param_t & param )
{
  
  const int k = param.has( "k" ) ? param.requires_int( "k" ) : 100 ; 

  const bool verbose = param.has( "verbose" ) ; 
  
  // iterate over signals
  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  
  const int ns = signals.size();

  logger << "  evaluating Otsu thresholds:";
  
  for (int s=0; s<ns; s++)
    {
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) 
	continue;
      
      writer.level( signals.label(s) , globals::signal_strat );
      
      logger << " " << signals.label(s) ;
      
      interval_t interval = edf.timeline.wholetrace();
      
      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();
      
      run_otsu( *d , k );
      
      writer.unlevel( globals::signal_strat );
    }
  
  logger << "\n";
  
}

void dsptools::cmdline_otsu( param_t & param )
{

  std::vector<double> x = dsptools::readcin();
  
  const int k = param.has( "k" ) ? param.requires_int( "k" ) : 100 ; 

  dsptools::run_otsu( x , k );
  
}
  
void dsptools::run_otsu( const std::vector<double> & x , const int k )
{  
  std::map<double,double> tvals, fvals;      

  double f;  
  double th = MiscMath::threshold2( x , &f, k , &fvals , &tvals );
  
  logger << "  Otsu threshold = " << th << " percentile = " << f << "\n";

  writer.value( "EMPTH" , th );
  writer.value( "EMPF" , f );
  
  std::map<double,double>::const_iterator tt =  tvals.begin();
  while ( tt != tvals.end() )
    {
      writer.level( tt->first , "TH" );
      writer.value("SIGMAB" , tt->second );
      writer.value("F" , fvals[ tt->first ] );
      ++tt;
    }
  writer.unlevel( "TH" );
    
}
  

//
// qdynam
//

void dsptools::qdynam( edf_t & edf , param_t & param )
{
  
  qdynam_t qd;

  qd.init( edf , param );
  
  const int ne = edf.timeline.first_epoch();

  // assume an ID field to select rows
  const bool ignore_id = param.has( "no-id" ) ? param.yesno( "no-id" ) : false ; 
  
  // get input(s)
  // match on ID
  // look for 'E'
  // pull all vars
  
  const std::set<std::string> vars = param.strset( "vars" );

  // facs must be the same across all input files (although vars can be different)
  const std::set<std::string> facs = param.strset( "facs" );
  
  std::vector<std::string> inputs = param.strvector( "inputs" );
  
  for (int i=0; i<inputs.size(); i++)
    {

      const std::string filename = Helper::expand( inputs[i] );
      if ( ! Helper::fileExists( filename ) )
	{
	  logger << "  *** could not open " << inputs[i] << "\n";
	  continue;
	}
      
      int n = -1;
      int slot_e = -1;
      int slot_id = -1;
      
      std::map<std::string,int> var2slot;
      std::map<std::string,int> fac2slot;
      
      std::ifstream IN1( filename.c_str() , std::ios::in );

      //
      // header 
      //
      
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() ) break;
      if ( line == "" ) continue;
      std::vector<std::string> hdr = Helper::parse( line , "\t" );
      n = hdr.size();
      for (int j=0; j<n; j++)
	{
	  if ( hdr[j] == "E" ) slot_e = j;
	  else if ( hdr[j] == "ID" ) slot_id = j;
	  else if ( facs.find( hdr[j] ) != facs.end() ) fac2slot[ hdr[j] ] = j;	    
	  else if ( vars.size() == 0 || vars.find( hdr[j] ) != vars.end() ) var2slot[ hdr[j] ] = j;
	}
      
      if ( slot_e == -1 )
	{
	  logger << "  ** no E column in " << inputs[i] << "\n";
	  break;
	}

      if ( slot_id == -1 && ! ignore_id )
	{
	  logger << "  ** no ID column in " << inputs[i] << "\n";
	  break;
	}

      if ( fac2slot.size() != facs.size() )
	{
	  logger << "  ** not all specified factors found in " << inputs[i] << "\n";
          break;
	}
      

      //
      // data rows
      //

      bool processed = false;
      
      while ( 1 )
	{
	  
	  std::string line;
	  Helper::safe_getline( IN1 , line );
	  if ( IN1.eof() ) break;
	  if ( line == "" ) continue;
	  std::vector<std::string> row = Helper::parse( line , "\t" );
	  if ( n != row.size() ) Helper::halt( "bad format in " + inputs[i] + " - variable # of cols" );

	  // row does match required ID?
	  if ( ( ! ignore_id ) && row[ slot_id ] != edf.id ) continue;
	  
	  // expecting 1-based in input; but qdynam wants display epoch # -1 
	  int epoch;
	  if ( ! Helper::str2int( row[slot_e] , &epoch ) )
	    Helper::halt( "bad format in " + inputs[i] + " - invalid epoch code" );	  
	  
	  processed = true;
	  
	  // get fac/lvl pairs
	  std::map<std::string,int>::const_iterator ff = fac2slot.begin();
          while ( ff != fac2slot.end() )
	    {
	      writer.level( row[ff->second], ff->first );
	      ++ff;
	    }
	  
	  // store values
	  std::map<std::string,int>::const_iterator ii = var2slot.begin();
	  while ( ii != var2slot.end() )
	    {
	      double x;
	      if ( ! Helper::str2dbl( row[ ii->second ] , &x ) )
		Helper::halt( "bad numeric format for " + inputs[i] + "\n" + line );
	      qd.add( writer.faclvl_notime() , ii->first , epoch - 1  , x );
	      ++ii;
	    }
	  
	}
      IN1.close();

      // undo factors
      if ( processed )
	{
	  std::map<std::string,int>::const_iterator ff = fac2slot.begin();
	  while ( ff != fac2slot.end() )
	    {
	      writer.unlevel( ff->first );
	      ++ff;
	    }
	}
      

      //
      // report
      //

      qd.proc_all();

      // next dataset
    }
  

}




void dsptools::make_bands( edf_t & edf , param_t & param )
{

  // signals: as we are iteratively reading/adding channels, postpone
  // creation of the sample list, then reading

  std::vector<std::string> siglabels = param.strvector( "sig" );
    
  // new label (default {orig}"_"{band}
  // for env:           {orig}"_"{band}_ht_mag

  const std::string tag = param.has( "tag" ) ? param.value( "tag" ) : "_";
  
  // S --> S_1, S_2, S_3, ...
  // default: S --> S_SLOW, S_DELTA, S_THETA , ... 
  const bool numeric = param.has( "numeric" ) ; 

  // envelopes
  const bool flt = param.has( "filtered" ) ? param.yesno( "filtered" ) : true ;  
  const bool env = param.has( "envelope" ) ? param.yesno( "envelope" ) : false ; 

  // faster filter
  const bool butterworth = param.has( "butterworth" );
  const int  butterworth_order = butterworth && ! param.empty( "butterworth" ) ? param.requires_int( "butterworth" ) : 4 ; 
												       
  if ( ! ( flt || env ) ) return;
  const int nx = flt + env;
    
  // bands
  std::vector<freq_range_t> bands;
  std::vector<std::string> labels;  
  
  if ( param.has( "bands" ) )
    {
      std::vector<std::string> tok = param.strvector( "bands" );
      std::set<std::string> utok;

      for (const auto &s : tok) {
	std::string u = s;
	std::transform(u.begin(), u.end(), u.begin(),
		       [](unsigned char c){ return std::toupper(c); });
	utok.insert(u);
      }

      if ( utok.find( "SLOW" ) != utok.end() ) { bands.push_back( globals::freq_band[ SLOW ] ); labels.push_back( "SLOW" ) ; }
      if ( utok.find( "DELTA" ) != utok.end() ) { bands.push_back( globals::freq_band[ DELTA ] ); labels.push_back( "DELTA" ) ; }      
      if ( utok.find( "THETA" ) != utok.end() ) { bands.push_back( globals::freq_band[ THETA ] ); labels.push_back( "THETA" ) ; }
      if ( utok.find( "ALPHA" ) != utok.end() ) { bands.push_back( globals::freq_band[ ALPHA ] ); labels.push_back( "ALPHA" ) ; }
      if ( utok.find( "SIGMA" ) != utok.end() ) { bands.push_back( globals::freq_band[ SIGMA ] ); labels.push_back( "SIGMA" ) ; }
      if ( utok.find( "SLOW_SIGMA" ) != utok.end() ) { bands.push_back( globals::freq_band[ LOW_SIGMA ] ); labels.push_back( "SLOW_SIGMA" ) ; }
      if ( utok.find( "FAST_SIGMA" ) != utok.end() ) { bands.push_back( globals::freq_band[ HIGH_SIGMA ] ); labels.push_back( "FAST_SIGMA" ) ; }
      if ( utok.find( "BETA" ) != utok.end() ) { bands.push_back( globals::freq_band[ BETA ] ); labels.push_back( "BETA" ) ; }
      if ( utok.find( "GAMMA" ) != utok.end() ) { bands.push_back( globals::freq_band[ GAMMA ] ); labels.push_back( "GAMMA" ) ; }
    }
  else if ( param.has( "freqs" ) )
    {
      // expect a-b,c-d,...
      std::vector<std::string> tok = param.strvector( "freqs" );
      
      for (const auto& s : tok) {
        std::stringstream ss(s);
        double lwr, upr;
        char dash;

        if (!(ss >> lwr >> dash >> upr) || dash != '-' || !ss.eof()) {
	  throw std::runtime_error("Invalid range format: " + s);
        }
	
        if (lwr >= upr) {
	  throw std::runtime_error("Invalid range (lwr >= upr): " + s);
        }

	std::string label = Helper::int2str( (int)(bands.size() + 1 ) );
	labels.push_back( label );
	bands.push_back( std::pair<double,double>( lwr, upr ) );
      }
    }
  else // default
    {

      bands.push_back( globals::freq_band[ SLOW ] ); labels.push_back( "1_SLOW" ) ; 
      bands.push_back( globals::freq_band[ DELTA ] ); labels.push_back( "2_DELTA" ) ; 
      bands.push_back( globals::freq_band[ THETA ] ); labels.push_back( "3_THETA" ) ; 
      bands.push_back( globals::freq_band[ ALPHA ] ); labels.push_back( "4_ALPHA" ) ;

      bool split_sigma = param.has( "split-sigma" ) ? param.yesno( "split-sigma" ) : false;

      if ( split_sigma )
	{
	  bands.push_back( globals::freq_band[ LOW_SIGMA ] ); labels.push_back( "5_SLOW_SIGMA" ) ; 
	  bands.push_back( globals::freq_band[ HIGH_SIGMA ] ); labels.push_back( "6_FAST_SIGMA" ) ;
	  bands.push_back( globals::freq_band[ BETA ] ); labels.push_back( "7_BETA" ) ; 
	  bands.push_back( globals::freq_band[ GAMMA ] ); labels.push_back( "8_GAMMA" ) ; 
	}
      else	
	{
	  bands.push_back( globals::freq_band[ SIGMA ] ); labels.push_back( "5_SIGMA" ) ;
	  bands.push_back( globals::freq_band[ BETA ] ); labels.push_back( "6_BETA" ) ; 
	  bands.push_back( globals::freq_band[ GAMMA ] ); labels.push_back( "7_GAMMA" ) ; 
	}
            
    }


  
  //
  // now make the channels
  //

  const int ns = siglabels.size();

  const bool NO_ANNOTS = true;

  for (int s=0; s<ns; s++)
    {

      signal_list_t signals = edf.header.signal_list( siglabels[s] , NO_ANNOTS );    
      
      if ( signals.size() != 1 ) continue;
     
      std::string slab = signals.label(0);
      
      logger << "  " << slab << ":\n";

      for (int i=0;i<labels.size(); i++)
	{

	  double lwr = bands[i].first;
	  double upr = bands[i].second;
	  std::string lu = Helper::dbl2str( lwr ) + "," + Helper::dbl2str( upr ); 
	      
	  // heuristic to get transition widths
	  double tw_lwr = 0.25 > 0.20 * lwr ? 0.25 : 0.20 * lwr ;
	  double tw_upr = 2    > 0.10 * upr ? 2    : 0.10 * upr ;
	  std::string tw = Helper::dbl2str( tw_lwr ) + "," + Helper::dbl2str( tw_upr );
	  
	  std::string ripple = "0.001";	      
	  std::string slab2 = tag + labels[i] ;

	  logger << "  --> " << slab << slab2 << " ( " << bands[i].first << " - " << bands[i].second << " Hz)\n";
	  
	  // make copy
	  param_t copy_param;
	  copy_param.add( "sig" , slab );
	  copy_param.add( "tag" , slab2 );
	  copy_param.add( "silent" );
	  proc_copy_signal( edf , copy_param );	  
	  
	  // apply filter (Kaiser IIR or Butterworth 4-th order IIR)
	  param_t filter_param;
	  filter_param.add( "sig" , slab + slab2 );	  
	  filter_param.add( "silent" );
	  filter_param.add( "bandpass" , lu );
	  if ( butterworth ) {
	    filter_param.add( "butterworth" , Helper::int2str( butterworth_order ) );
	  }
	  else {
	    filter_param.add( "tw" , tw );
	    filter_param.add( "ripple" , ripple );
	  }
	  proc_filter( edf , filter_param );
	  
	  // apply hilbert transform?
	  if ( env )
	    {
	      param_t hilbert_param;
	      hilbert_param.add( "sig" , slab + slab2 );	      
	      dsptools::hilbert( edf , hilbert_param );
	    }

	  // clean up flt if we only want envelope
	  if ( ! flt )
	    {
	      int s1 = edf.header.signal( slab + slab2 );
	      if ( s1 != -1 ) 		
		edf.drop_signal( s1 );
	    }
	}
      logger << "\n";
    }
}
