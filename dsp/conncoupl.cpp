
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

#include "conncoupl.h"
#include "edf/edf.h"
#include "eval.h"

#include "db/db.h"
#include "helper/logger.h"
#include "dsp/hilbert.h"
#include "edf/slice.h"
#include "miscmath/crandom.h"
#include "cwt/cwt.h"
#include "stats/matrix.h"

extern writer_t writer;

extern logger_t logger;

//
// Wrapper
//

void dsptools::connectivity_coupling( edf_t & edf , param_t & param )
{
  
  //
  // Signals
  //

  std::string signal_label = param.requires( "sig" );

  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();

    
  //
  // Sampling rates
  //
  
  std::vector<double> Fs = edf.header.sampling_freq( signals );

  for ( int s=1; s< Fs.size(); s++)
    if ( Fs[s] != Fs[0] )
      Helper::halt( "all sampling rates must be equal across channels: use RESAMPLE" );


  //
  // Number of permutation replicates
  //

  const int nreps = param.has( "nreps" ) ? param.requires_int( "nreps" ) : 0;

  //
  // analysis split by epoch (and end results are averages over epochs)
  // or do whole signal in one go?  (in which case, epoch size is set to 0)
  //

  if ( param.has( "epoch" ) && ! edf.timeline.epoched() )
    Helper::halt( "no EPOCHs set" );

  // use 30 seconds if not otherwise specified
  const double epoch_sec = param.has( "epoch" ) ? edf.timeline.epoch_length() : 30 ;

  

  //
  // TF-decomposition approach
  //
  
  bool use_hilbert = param.has( "hilbert" );
   
  bool do_pac = param.has( "pac" );
  
  bool do_xch = param.has( "xch" );
  
  bool do_xpac = param.has( "xpac" );
    

      
  //
  // filter-Hilbert approach
  //

  if ( use_hilbert )
    {

      double w = param.requires_dbl( "w" );
	    
      std::vector<double> lf0 = param.dblvector( "lwr" );
      std::vector<freq_range_t> fb1;
      
      if ( lf0.size() < 1 || lf0.size() > 3 ) Helper::halt( "expecting lwr to have 1,2 or 3 values: start,end,step" );
      if ( lf0.size() == 1 ) lf0.push_back( lf0[0] );
      if ( lf0.size() == 2 ) lf0.push_back( 1 );
      if ( lf0[0] - w < 0 ) Helper::halt( "bad format for lwr, lower value too low given w" );
      for (double v = lf0[0] ; v <= lf0[1] ; v += lf0[2] )
	fb1.push_back( freq_range_t( v-w , v+w ) );
      
      // optional, second freq. range

      std::vector<freq_range_t> fb2; // used for cross-frequency (PAC) measures
      if ( param.has( "upr" ) )
	{
	  std::vector<double> uf0 = param.dblvector( "upr" );

	  if ( uf0.size() == 1 ) uf0.push_back( uf0[0] );
	  
	  if ( uf0.size() < 1 || uf0.size() > 3 ) Helper::halt( "expecting upr to have 1,2 or 3 values: start,end,step" );

	  if ( lf0[0] > lf0[1] ) Helper::halt( "bad format for lwr" );
	  
	  if ( uf0.size() == 2 ) uf0.push_back( 1 );
	  
	  if ( uf0[0] > uf0[1] ) Helper::halt( "bad format for upr" );
	  
	  if ( uf0[0] - w < 0 ) Helper::halt( "bad format for upr, lower value too low given w" );
	  
	  for (double v = uf0[0] ; v <= uf0[1] ; v += uf0[2] )
	    fb2.push_back( freq_range_t( v-w , v+w ) );
	}


      // require at least 4 Hz separation between 
      //      double a2c_frq_gap = param.has( "pac-gap" ) ? param.require_dbl( "pac-gap" ) : 0 ; 

      // filter parameters
            
      double ripple = param.has( "ripple" ) ? param.requires_dbl( "ripple" ) : 0.05;
      
      double tw = param.has( "tw" ) ? param.requires_dbl( "tw" ) : 2 ;

      // output
      bool epoch_level_output = ! param.has( "no-epoch-output" );
      
      // do it
      
      conncoupl_t cc( edf , signals , Fs[0] ,
		      fb1 , fb2 , 
		      ripple, tw , 
		      nreps ,
		      epoch_sec ,
		      do_pac , do_xch , do_xpac ,
		      epoch_level_output );


    }
  
  
  //
  // OR Wavelet parameters
  //

  if ( ! use_hilbert )
    {
    
      // for cross-channel connectivity
      // for now, use alternate/ 'Cox' specification
      
      std::vector<double> fc, fwhm;
      int num = 0 ;
      // F_C
      if ( param.has( "fc" ) ) fc = param.dblvector( "fc" );
      else if ( param.has( "fc-range" ) )
	{	  
	  fc = param.dblvector( "fc-range" );
	  if ( fc.size() != 2 ) Helper::halt( "expecting fc-range=min/max" );
	  num = param.requires_int( "num" );
	  fc = param.has( "linear" ) ? MiscMath::linspace( fc[0] , fc[1] , num ) : MiscMath::logspace( fc[0] , fc[1] , num );
	}

      // FWHM
      if ( param.has( "fwhm" ) ) fwhm = param.dblvector( "fwhm" );
      else if ( param.has( "fwhm-range" ) )
	{	  
	  fwhm = param.dblvector( "fwhm-range" );
	  if ( fwhm.size() != 2 ) Helper::halt( "expecting fwhm-range=min/max" );
	  num = param.requires_int( "num" );
	  fwhm = param.has( "linear" ) ? MiscMath::linspace( fwhm[0] , fwhm[1] , num ) : MiscMath::logspace( fwhm[0] , fwhm[1] , num );
	}

      // secondary frequency band(s)

      std::vector<double> fc2, fwhm2;
      int num2 = 0 ; 
      // F_C
      if ( param.has( "fc2" ) ) fc2 = param.dblvector( "fc2" );
      else if ( param.has( "fc-range2" ) )
	{	  
	  fc2 = param.dblvector( "fc-range2" );
	  if ( fc2.size() != 2 ) Helper::halt( "expecting fc-range2=min/max" );
	  num2 = param.requires_int( "num2" );
	  fc2 = param.has( "linear" ) ? MiscMath::linspace( fc2[0] , fc2[1] , num2 ) : MiscMath::logspace( fc2[0] , fc2[1] , num2 );
	}

      // FWHM
      if ( param.has( "fwhm2" ) ) fwhm2 = param.dblvector( "fwhm2" );
      else if ( param.has( "fwhm-range2" ) )
	{	  
	  fwhm2 = param.dblvector( "fwhm-range2" );
	  if ( fwhm2.size() != 2 ) Helper::halt( "expecting fwhm-range2=min/max" );
	  num2 = param.requires_int( "num2" );
	  fwhm2 = param.has( "linear" ) ? MiscMath::linspace( fwhm2[0] , fwhm2[1] , num2 ) : MiscMath::logspace( fwhm2[0] , fwhm2[1] , num2 );
	}
      
      // check consistency
      if ( param.has( "num" ) && ( param.has( "fc" ) || param.has( "fwhm" ) ) )
	Helper::halt( "cannot use fc/fwhm and num: use fc-range/fwhm-range" );
      if ( param.has( "num2" ) && ( param.has( "fc2" ) || param.has( "fwhm2" ) ) )
	Helper::halt( "cannot use fc2/fwhm2 and num2: use fc-range2/fwhm-range2" ) ;
      
      if ( fc.size() == 0 ) Helper::halt( "bad specification of fc/fwhm" );

      // if FWHM not specified at at, use defaiults
      if ( fwhm.size() == 0 )
	{
	  fwhm.resize( fc.size() );
	  for (int i=0;i<fc.size();i++) fwhm[i] = CWT::pick_fwhm( fc[i] );
	}

      if ( fwhm2.size() == 0 )
	{
          fwhm2.resize( fc2.size() );
          for (int i=0;i<fc2.size();i++) fwhm2[i] = CWT::pick_fwhm( fc2[i] );
        }
            
      if ( fc.size() != fwhm.size() ) Helper::halt( "bad specification of fc/fwhm" );
      if ( fc2.size() != fwhm2.size() ) Helper::halt( "bad specification of fc2/fwhm2" );

      // wavelet duration (in seconds)
      double tline = param.has( "length" ) ? param.requires_dbl( "length" ) : 20 ; 

      // 
      
      conncoupl_t cc( edf , signals , Fs[0] ,
		      fc , fwhm , num , 
		      fc2 , fwhm2 , num2 ,
		      tline , 
		      nreps ,
		      epoch_sec ,
		      param.has( "pac" ) ,
		      param.has( "xch" ) , 
		      param.has( "xpac" ) ,		      
		      ! param.has( "no-epoch-output") ,
		      param.has( "dump-wavelets" ) 		      
		      );

    }

  
}



void conncoupl_t::setup()
{

  //
  // offsets used by all shifts for a given replicate
  //

  es_pts = es * sr ; 
    
  offset.resize( nreps );
  
  for (int i=0;i<nreps;i++)
    offset[i] = CRandom::rand( es_pts );

  
  //
  // contrats to report:
  // i) CFC-within channel (PAC) [ but only for f1 < f2 ]
  // ii) cross-channel, within frequency
  // iii) cross-channel, cross-frequency (xpac) [ but only for f1 < f2 ] 
  //

  if ( do_xch )
    {
      for (int si1=0;si1<signals.size();si1++)
	for (int si2=si1+1;si2<signals.size();si2++)
	  //for (int si2=0;si2<signals.size();si2++)
	  {
	    
	    if ( use_hilbert )
	      for (int f=0;f<fint1.size();f++)
		{
		  s1.push_back( si1 ); s2.push_back( si2 );
		  f1.push_back( str( fint1[f] ) ); f2.push_back( str( fint1[f] ) ); 
		  cfc.push_back( false );
		  xch.push_back( true );
		}
	    else
	      for (int f=0;f<fc1.size();f++)
		{
		  s1.push_back( si1 ); s2.push_back( si2 );
		  f1.push_back( str( freq_range_t( fc1[f] , fwhm1[f] ) ) ); f2.push_back( str( freq_range_t( fc1[f] , fwhm1[f] ) ) ); 
		  disp_f1.push_back(  fc1[f] ) ; disp_f2.push_back(  fc1[f] ) ; 
		  cfc.push_back( false );
		  xch.push_back( true );
		}
	  }
    }

  // not used for now... e.g. could set to 4Hz or so to
  // nesure 'gaps' between P and A bands
  const double pac_gap = 0;
  
  // within-channel, across frequency
  if ( do_pac )
    {
      for (int si1=0;si1<signals.size();si1++)
	{
	  if ( use_hilbert )
	    for (int fi1=0;fi1<fint1.size();fi1++)
	      for (int fi2=0;fi2<fint2.size();fi2++)
		{
		  if ( fint1[fi1] < fint2[fi2] )
		    {
		      s1.push_back( si1 ); s2.push_back( si1 );
		      f1.push_back( str( fint1[fi1] ) ); f2.push_back( str( fint2[fi2] ) ); 
		      cfc.push_back( true ); xch.push_back( false );
		    }
		}
	  else
	    for (int fi1=0;fi1<fc1.size();fi1++)
	      for (int fi2=0;fi2<fc2.size();fi2++)
		{
		  if ( fc2[fi2] - fc1[fi1] > pac_gap )
		    {
		      s1.push_back( si1 ); s2.push_back( si1 );
		      f1.push_back( str( freq_range_t( fc1[fi1] , fwhm1[fi1] ) ) ); f2.push_back( str( freq_range_t( fc2[fi2] , fwhm2[fi2] ) ) ); 
		      disp_f1.push_back(  fc1[fi1] ) ; disp_f2.push_back(  fc2[fi2] ) ; 
		      cfc.push_back( true ); xch.push_back( false );
		    }
		}
	}
    }
  
  // cross-channel, cross-frequency PAC
  if ( do_xpac )
    {
      for (int si1=0;si1<signals.size();si1++)
	for (int si2=0;si2<signals.size();si2++) // need to do all pairs i-j as well as j-i 
	  {
	    if ( si1 == si2 ) continue; // alreadt done within-channel CFC
	    if ( use_hilbert )
	      for (int fi1=0;fi1<fint1.size();fi1++)
		for (int fi2=0;fi2<fint2.size();fi2++)
		  {
		    if ( fint1[fi1] < fint2[fi2] )
		      {
			s1.push_back( si1 ); s2.push_back( si2 );
			f1.push_back( str( fint1[fi1] ) ); f2.push_back( str( fint2[fi2] ) ); 
			cfc.push_back( true ); xch.push_back( true );
		      }
		  }
	    else
	      for (int fi1=0;fi1<fc1.size();fi1++)
		for (int fi2=0;fi2<fc2.size();fi2++)
		  {
		    if ( fc2[fi2] - fc1[fi1] > pac_gap )
		      {
			s1.push_back( si1 ); s2.push_back( si2 );
			f1.push_back( str( freq_range_t( fc1[fi1] , fwhm1[fi1] ) ) ); f2.push_back( str( freq_range_t( fc2[fi2] , fwhm2[fi2] ) ) ); 
			disp_f1.push_back(  fc1[fi1] ) ; disp_f2.push_back(  fc2[fi2] ) ; 
			cfc.push_back( true ); xch.push_back( true );
		      }
		  }
	  }
    }

  logger << "  registered " << s1.size() << " channel/frequency combinations to evaluate per epoch\n";

  if ( s1.size() == 0 ) Helper::halt( "no combinations specified: add pac, xch and/or xpac options" );
}
  


void conncoupl_t::pre_calc()
{

  const int ns = signals.size();
  
  // get all unique frequencies/intervals
  // (using strings to avoid floating point equalities)    

  if ( use_hilbert )
    {
      for ( int i = 0 ; i < fint1.size() ; i++ )
	fmap[ str( fint1[i] ) ] = fint1[i];
      for ( int i = 0 ; i < fint2.size() ; i++ )
	fmap[ str( fint2[i] ) ] = fint2[i];
    }
  else // wavelet (just use first slot of freq_range_t for Fc, second slot for FWHM)
    {
      for ( int i = 0 ; i < fc1.size() ; i++ )
	fmap[ str( freq_range_t( fc1[i] , fwhm1[i] ) ) ] = freq_range_t( fc1[i] , fwhm1[i] ) ;
      for ( int i = 0 ; i < fc2.size() ; i++ )
	fmap[ str( freq_range_t( fc2[i] , fwhm2[i] ) ) ] = freq_range_t( fc2[i] , fwhm2[i] ) ;
    }


  //
  // iterate over epochs (then channels, then frequencies)
  //

  const int ne = edf.timeline.first_epoch();
  
  a.resize( ne );  
  a_conj.resize( ne );
  for (int e=0; e<ne; e++)
    {
      a[e].resize( ns );
      a_conj[e].resize( ns );
    }
  
  //
  // populate amplitude/complex measures 
  //

  const int nf = fmap.size();
  
  
  //
  // function only applies to entire current trace
  //

  
  interval_t interval = edf.timeline.wholetrace();
  
  //
  // Iterate over channels
  //
  
  for (int s=0;s<ns;s++)
    {
      
      //
      // Only data-channels
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) )
	Helper::halt( "can only apply to data channels" );
      
      
      //
      // Get signal (whole signal)
      //
      
      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();
      

      //
      // Optionally, first dump wavelets
      //

      if ( dump_wavelets )
	{
	  logger << "  dumping " << fmap.size() << " wavelets\n";
	  	  	      
	  std::map<std::string,freq_range_t>::const_iterator ff = fmap.begin();
	  while ( ff != fmap.end() )
	    {
	      
	      double freq = ff->second.first;
	      double fwhm = ff->second.second;
	      CWT cwt;
	      cwt.set_sampling_rate( sr );
	      cwt.set_timeframe( 50.0 / tlen );
	      cwt.alt_add_wavelet( freq , fwhm , tlen );	      

	      std::vector<double> t = cwt.get_timeframe();	      
	      std::vector<dcomp> w = cwt.alt_wavelet(0);
	      	      
	      const int n = w.size();

	      writer.level( freq , "F" );
	      writer.level( fwhm , "FWHM" );
	      
	      for (int i=0;i<n;i++)
	       	{
	       	  writer.level( t[i] , "SEC" );
	       	  writer.value( "REAL" , std::real( w[i] ) );
	       	  writer.value( "IMAG" , std::imag( w[i] ) );
	       	}
	      writer.unlevel( "SEC" );
	      
	      ++ff;
	    }

	  writer.unlevel( "F" );
	  writer.unlevel( "FWHM" );
	  
	}
      
      //
      // get wavelets per channel across all frequencies
      //
      
      if ( ! use_hilbert )
	{

	  //
	  // Initialize CWT
	  //

	  logger << "  estimating " << fmap.size() << " wavelets for " << signals.label(s) << "\n";
	  
	  CWT cwt;

	  cwt.set_sampling_rate( sr );
	  	      
	  std::map<std::string,freq_range_t>::const_iterator ff = fmap.begin();
	  while ( ff != fmap.end() )
	    {
	      
	      std::string label = ff->first;
	      double freq = ff->second.first;
	      double fwhm = ff->second.second;
	      cwt.set_timeframe( 50.0 / tlen );
	      cwt.alt_add_wavelet( freq , fwhm , tlen );	      

	      ++ff;
	    }


	  //
	  // Run wavelets
	  //

	  cwt.store_real_imag_vectors( true );
	  cwt.load( d );  
	  cwt.run_wrapped();
	

	  //
	  // store as epochs
	  //

	  int fc = 0;

	  ff = fmap.begin();
	  while ( ff != fmap.end() )
	    {	      
	      std::string label = ff->first;
	      double freq = ff->second.first;
	      double fwhm = ff->second.second;
	      
	      std::vector<dcomp> res =  cwt.get_complex( fc++ );
	      
	      int ec = 0;
	      for (int e = 0 ; e < ne ; e++ )
		{		  
		  std::vector<dcomp> tmp( es_pts );
		  for (int i=0;i<es_pts;i++) tmp[i] = res[ ec++ ];		  
		  a[ e ][ s ][ label ] = tmp;
		  // also store conjugate, for speed up of permutation
		  for (int i=0;i<es_pts;i++) tmp[i] = conj( tmp[i] );
		  a_conj[ e ][ s ][ label ] = tmp;
		}
	      
	      ++ff;
	  
	    }
	  
	}
	  
      
      // 	      // hilbert
      // double freq_lwr = ff->second.first;
	  //	      double freq_upr = ff->second.second;	
      
	  // if ( use_hilbert )
	  //   {
	  //     hilbert_t hilbert( *d , sr , freq_lwr , freq_upr , ripple , tw , true );
	  //     // TODO: store
	  //     //a[ ec ][ s ][ label ] = hilbert.get_complex();	      
	  //   }

      
      
      //
      // next signal
      //
    }

  
}




//
// Caclulate metrics
//

void conncoupl_t::calc()
{

  // at this point, for each epoch/signal we have real/imag parts of the signal (from wavelet or filter-Hilbert)
  // for each signal for each frequency/interval in both f1 and f2
  
  // across channel, consider all f1-f1 cross-channel pairings (for wPLI)
  // within channel, consider all f1xf2 intra-channel pairsings (for PAC, PLV) 
  
  // for this target frequency, get wavelets for all channels

  // number of tests
  const int nt = s1.size();

  // number of epochs
  const int ne = a.size();

  
  //
  // Store statistics (per epoch)
  //
  

  results[ "wPLI" ] = conncoupl_res_t( ne , nt );
  results[ "dPAC" ] = conncoupl_res_t( ne , nt );
  //  results[ "dPAC_PHASE" ] = conncoupl_res_t( ne , nt );
  
  
  //
  // For each epoch
  //

  logger << "  iterating over " << ne << " epochs...\n  ";

  for (int e=0;e<ne;e++)
    {
      
      if ( e % 50 == 49 ) logger << ". " << e+1 << " epochs\n  ";
      else logger << ".";
      
      
      //
      // Surrogate time-series (pre-calculate, fixed for all tests in a given epoch)
      // although independent between epochs (not that this should matter...)
      //      
      
      std::vector<std::vector<int> > shuffle( nreps );
      for (int r=0;r<nreps;r++)
	{
	  int offset = CRandom::rand( es_pts );
	  shuffle[r].resize( es_pts );
	  for (int i=0;i<es_pts;i++)
	    {
	      if ( offset >= es_pts ) offset = 0;
	      shuffle[r][i] = offset++;
	    }
	}


      //
      // For each test 
      //
  
      for (int t=0; t<nt; t++)
	{
	  
	  const bool pac = cfc[t];
	  
	  //
	  // CFC
	  //

	  if ( pac )
	    {

	      //
	      // obtain the two frequencies
	      //
	      
	      const std::string & flabel1 = f1[ t ] ;
	      const std::string & flabel2 = f2[ t ] ;
	      
	      //
	      // extract two relevant vectors (cross-channel)
	      //
	      
	      const std::vector<dcomp> & x = a[ e ][ s1[t] ][ flabel1 ];
	      const std::vector<dcomp> & y = a[ e ][ s2[t] ][ flabel2 ];
	      	      
	      //	      
	      // calculate phase for 'x' and magnitude squared for 'y'cross-spectral density
	      //

	      dcomp debias_term(0,0);
	      
	      std::vector<dcomp> ph( es_pts ) , mag( es_pts );

	      for (int i=0; i<es_pts; i++)
		{		  
		  ph[i] = exp( dcomp( 0 , arg( x[i] ) ) );
		  mag[i] = dcomp( pow( abs( y[i] ) , 2 ) , 0 );
		  // accrue debias term (which is similar across all null perms)
		  debias_term += ph[i];
		}
	      
	      // get mean
	      debias_term /= dcomp( (double)es_pts , 0 );

	      	      
	      //
	      // dPAC 
	      //

	      dcomp dPAC(0,0);
	      
	      for (int i=0; i<es_pts; i++)
		dPAC += ( ph[i] - debias_term ) * mag[i] ;
	      
	      dPAC /= dcomp( es_pts , 0 );

	      double obs_dPAC = abs( dPAC );
	      //	      double obs_dPAC_PHASE = arg( dPAC );

	      results[ "dPAC" ].stats( e , t ) = obs_dPAC;
	      // results[ "dPAC_PHASE" ].stats( e , t ) = obs_dPAC_PHASE;

	      //
	      // Surrogates
	      //

	      std::vector<double> null_stats( nreps );
	      std::vector<double> null_stats2( nreps );
	      
	      for (int r=0; r<nreps; r++)
		{
		  
		  dcomp perm_dPAC(0,0);

		  const std::vector<int> & sh = shuffle[r];
		  
		  for (int i=0; i<es_pts; i++)
		    perm_dPAC += ( ph[i] - debias_term ) * mag[ sh[i] ] ;
		  
		  perm_dPAC /= dcomp( es_pts , 0 );
		  
		  null_stats[r] = abs( perm_dPAC );
		  null_stats2[r] = arg( perm_dPAC );
		  
		}

	      if ( nreps )
		{
		  double mean = MiscMath::mean( null_stats );
		  double sd   = MiscMath::sdev( null_stats , mean );
		  results[ "dPAC" ].emp_z( e , t ) = ( obs_dPAC - mean ) / sd;
		  //std::cout << "xx " << e << " " << s1[t] << " " << s2[t] << " " << flabel1 << " " << flabel2 << " obs " << obs_dPAC << "\t" << mean << "\t" << sd << "\t" << ( obs_dPAC - mean ) / sd << "\n";
		  // mean = MiscMath::mean( null_stats2 );
		  // sd   = MiscMath::sdev( null_stats2 , mean );
		  // results[ "dPAC_PHASE" ].emp_z( e , t ) = ( obs_dPAC_PHASE - mean ) / sd;
		  
		}
	      
	    }
	  

	  
	  //
	  // Connectivity
	  //
	  
	  if ( ! pac )
	    {
	      
	      const std::string & flabel = f1[ t ] ;
	      
	      //
	      // extract two relevant vectors (cross-channel)
	      //
	      
	      const std::vector<dcomp> & x = a[ e ][ s1[t] ][ flabel ];
	      const std::vector<dcomp> & y = a_conj[ e ][ s2[t] ][ flabel ];
	      const int np = x.size();
	      
	      
	      //	      
	      // cross-spectral density
	      //
	      
	      std::vector<dcomp> sxy( np );
	      std::vector<double> isxy( np );
	      for (int i=0; i<np; i++)
		{
		  sxy[i] = x[i] * y[i] ;
		  isxy[i] = std::imag( sxy[i] );
		}
	      
	      
	      //
	      // wPLI: abs( mean( imag(X) )  ) / mean( abs( imag(X) ) )
	      //
	      
	      double numer = 0 , denom = 0; 
	      for (int i=0; i<np; i++)
		{
		  numer += isxy[i] ;
		  denom += abs( isxy[i] );
		}
	      
	      double wPLI = abs( numer / (double)np ) / ( denom / (double)np );
	      
	      results[ "wPLI" ].stats( e , t ) = wPLI;

	    }

	  //
	  // Next test
	  //
	}

              
 

      //
      // Iterate over tests on null data
      //
      
      for (int t=0; t<nt; t++)
	{

	  // assume just one test (for now...)

	  std::vector<double> null_stats( nreps ) ;

	  std::vector<dcomp> sxy( es_pts );
	  std::vector<double> isxy( es_pts );
	
	  //
	  // For each null replicate 
	  //
	  
	  for (int r=0; r<nreps; r++)
	    {
	      
	      //
	      // phase-amplitude coupling test?
	      //

	      const bool pac = cfc[t];
	      

	      //
	      // Connectivity
	      //
	      
	      if ( ! pac )
		{
		  
		  const std::string & flabel = f1[ t ] ;
		  
		  //
		  // extract two relevant vectors (cross-channel), with one shuffled
		  //
		  
		  const std::vector<dcomp> & x = a[ e ][ s1[t] ][ flabel ];
		  const std::vector<dcomp> & y_conj = a_conj[ e ][ s2[t] ][ flabel ];		  
		  
		  if ( x.size() != es_pts ) Helper::halt( "probo" );

		  //	      
		  // cross-spectral density, except with one channel shuffled
		  //

		  double wPLI_numer = 0 , wPLI_denom = 0;
		  
		  for (int i=0; i< es_pts; i++)
		    {
		      //sxy[i] = x[i] * conj( y[ shuffle[r][i] ] );
		      sxy[i] = x[i] * y_conj[ shuffle[r][i] ];
		      isxy[i] = std::imag( sxy[i] );
		      
		  
		      //
		      // wPLI: abs( mean( imag(X) )  ) / mean( abs( imag(X) ) )
		      //
		  
		      wPLI_numer += isxy[i] ;
		      wPLI_denom += abs( isxy[i] );
		    }

		  //
		  // calc stats for this epoch 
		  //

		  double wPLI = abs( wPLI_numer / (double)es_pts ) / ( wPLI_denom / (double)es_pts );
		  
		  null_stats[ r ] = wPLI ;

		}
	      
	      //
	      // Next null replicate
	      //
	      
	    }
	  
	  //
	  // wPLI Z score
	  //

	  if ( nreps )
	    {
	      double mean = MiscMath::mean( null_stats );
	      double sd   = MiscMath::sdev( null_stats , mean );	  
	      results[ "wPLI" ].emp_z( e , t ) = ( results[ "wPLI" ].stats( e , t ) - mean ) / sd;
	    }
	  
	  //
	  // Next test
	  //

	}    


      //
      // Next epoch
      //
    }	      
	      
  logger << " done\n";

      
  //
  // Average results over epochs, and report results
  //

  Data::Vector<double> mean_s_wPLI = Statistics::mean( results[ "wPLI" ].stats );
  Data::Vector<double> mean_s_dPAC = Statistics::mean( results[ "dPAC" ].stats );

  Data::Vector<double> mean_z_wPLI, mean_z_dPAC;
  if ( nreps )
    {
      mean_z_wPLI = Statistics::mean( results[ "wPLI" ].emp_z );
      mean_z_dPAC = Statistics::mean( results[ "dPAC" ].emp_z );
    }

  //
  // calculate empirical p-value based on distrubution of mean STAT for each test, over epochs
  //
  
  for (int t=0;t<nt;t++)
    {

      writer.level( signals.label( s1[t] ) , "CH1" );
      writer.level( signals.label( s2[t] ) , "CH2" );
      
      if ( use_hilbert )
	{
	  writer.level( f1[t] , "F1" );
	  writer.level( f1[t] , "F2" );	      
	}
      else
	{
	  writer.level( disp_f1[t] , "F1" );
	  writer.level( disp_f2[t] , "F2" );
	}
      
      writer.value( "CFC" , cfc[t] );
      writer.value( "XCH" , xch[t] );      

      if ( xch[t] & ! cfc[t] ) {	
	writer.value( "wPLI" , mean_s_wPLI[t] );
	if ( nreps ) writer.value( "wPLI_Z" , mean_z_wPLI[t] );
      }
      
      if ( cfc[t] ) {	
	writer.value( "dPAC" , mean_s_dPAC[t] );
	if ( nreps ) writer.value( "dPAC_Z" , mean_z_dPAC[t] );
      }
      
    }
  
  writer.unlevel( "CH1" );
  writer.unlevel( "CH2" );
  writer.unlevel( "F1" );
  writer.unlevel( "F2" );

  //
  // epoch level stats
  //

  if ( epoch_level_output )
    {
      
      for (int e=0;e<ne;e++)
	for (int t=0;t<nt;t++)
	  {
	    writer.epoch( edf.timeline.display_epoch( e ) );
	    
	    writer.level( signals.label( s1[t] ) , "CH1" );
	    writer.level( signals.label( s2[t] ) , "CH2" );
	    
	    if ( use_hilbert )
	      {
		writer.level( f1[t] , "F1" );
		writer.level( f2[t] , "F2" );	
	      }
	    else
	      {
		writer.level( disp_f1[t] , "F1" );
		writer.level( disp_f2[t] , "F2" );	
	      }
	    
	    if ( xch[t] ) 
	      {
		writer.value( "wPLI" , results["wPLI"].stats(e,t) );
		if ( nreps )
		  writer.value( "wPLI_Z" , results["wPLI"].emp_z(e,t) );
	      }
	    
	    if ( cfc[t] ) 
	      {
		writer.value( "dPAC" , results["dPAC"].stats(e,t) );
		//	    writer.value( "dPAC_PHASE" , results["dPAC_PHASE"].stats(e,t) );
		if ( nreps )
		  {
		    writer.value( "dPAC_Z" , results["dPAC"].emp_z(e,t) );
		    //writer.value( "dPAC_Z_PHASE" , results["dPAC_PHASE"].emp_z(e,t) );
		  }
	      }
	  }
      
      writer.unepoch();
      writer.unlevel( "CH1" );
      writer.unlevel( "CH2" );
      writer.unlevel( "F1" );
      writer.unlevel( "F2" );
    }
  
  
  //
  // All done
  //

}
  




// // void phsyn_t::calc() 
// // {

// //   //
// //   // bin boundaries (0..360)
// //   //
  
// //   double bs = 360 / (double)nbins;  
  
// //   std::vector<double> bb( nbins );
  
// //   bb[0] = 0;
  
// //   for (int i=1;i<nbins;i++) bb[i] = bb[i-1] + bs;
    
// //   //
// //   // size output grids
// //   //
  
// //   obs.resize( nbins );
// //   perm.resize( nbins );
// //   pv.resize( nbins );
// //   z.resize( nbins );
// //   z2.resize( nbins );
  
// //   for (int b=0;b<nbins;b++) 
// //     {      
// //       obs[b].resize( nbins , 0 );
// //       perm[b].resize( nbins , 0 );
// //       pv[b].resize( nbins , 0 );
// //       z[b].resize( nbins , 0 );
// //       z2[b].resize( nbins , 0 );
// //     }
  

// //   //
// //   // get phase for each f1 and f2 frequency band  
// //   //
  
// //   std::map<freq_range_t,std::vector<double> > ph;  
  
// //   std::set<freq_range_t> frqs;
// //   for (int f=0;f<f1.size();f++) frqs.insert( f1[f] );
// //   for (int f=0;f<f2.size();f++) frqs.insert( f2[f] );
  
// //   std::set<freq_range_t>::const_iterator ff = frqs.begin();
  
// //   logger << "  considering " << frqs.size() << " frequency bands\n";

// //   //
// //   // get phase angles at each time point
// //   //  

// //   while ( ff != frqs.end() )
// //     {
      
// //       logger << "  hilbert " << ff->first << "-" << ff->second << "\n";

// //       // filter-Hilbert
// //       hilbert_t hilbert( x , sr , ff->first , ff->second , ripple , tw );
      
// //       // convert to degrees with 0 as pos-to-neg crossing
// //       std::vector<double> angle = *hilbert.phase();
// //       for (int i=0;i<angle.size();i++) angle[i] = MiscMath::as_angle_0_pos2neg( angle[i] );
// //       // store
// //       ph[ *ff ] = angle;            
// //       ++ff;
// //     }
  

// //   const int npoints = x.size();
  
// //   //
// //   // create observed frequency bin
// //   //
        
// //   for (int i1=0;i1<f1.size();i1++) 
// //     for (int i2=0;i2<f2.size();i2++) 
// //       {

// // 	logger << f1[i1].first << " and " << f2[i2].first << "\n";

// // 	//
// // 	// Clear working count matrix
// // 	//
	
// // 	for (int b1=0;b1<nbins;b1++)
// // 	  for (int b2=0;b2<nbins;b2++)
// // 	    obs[b1][b2] = perm[b1][b2] = z[b1][b2] = z2[b1][b2] = pv[b1][b2] = 0 ;
	
// // 	const std::vector<double> * ph1 = &ph[ f1[i1] ];
// // 	const std::vector<double> * ph2 = &ph[ f2[i2] ];
	
// // 	int b1 = 0, b2 = 0; 
// // 	for (int i=0;i<npoints;i++)
// // 	  {
// // 	    bin( (*ph1)[i] , &b1 , bb , nbins );
// // 	    bin( (*ph2)[i] , &b2 , bb , nbins );

// // 	    // increment by one unit (default, unweighted)
// // 	    obs[ b1 ][ b2 ]++;
	    
// // 	  }
		

// // 	//
// // 	// Observed summary statistic
// // 	//

// // 	double obs_stat = test_uniform( obs );

	
// // 	//
// // 	// Replicates
// // 	//
	
// // 	int emp_stat = 0;
// // 	std::vector<double> perm_stats;

// // 	for (int r = 0; r < nreps; r++) 
// // 	  {

// // // 	    if ( r % 10 == 0 ) logger << ".";
// // // 	    if ( r % 100 == 0 ) logger << "\n";
	    
// // 	    // clear working perm matrix
	    
// // 	    for (int b1=0;b1<nbins;b1++)
// // 	      for (int b2=0;b2<nbins;b2++)
// // 		perm[b1][b2] = 0;

// // 	    // initially, whole trace perm
	    
// // 	    const std::vector<double> * ph1 = &ph[ f1[i1] ];
// // 	    const std::vector<double> * ph2 = &ph[ f2[i2] ];
	
// // 	    int b1 = 0, b2 = 0; 

// // 	    // within-epoch perm
	    
// // 	    std::vector<int> permpos( npoints );
// // 	    int j = 0;

// // 	    if ( es ) 
// // 	      {
// // 		int num_epoch = npoints / es ; 
// // 		std::vector<int> shuffle;
// // 		for (int e=0;e<num_epoch;e++) shuffle.push_back( CRandom::rand( es ) );
// // 		//std::cout << "fig " << num_epoch << " epochs\n";
// // 		for (int i=0;i<npoints;i++)
// // 		  {		    
// // 		    int cur_epoch = i / es;
// // 		    if ( cur_epoch >= num_epoch ) continue; // do not permute last epoch, if partial
// // 		    int off = i - cur_epoch * es; 
// // 		    if ( off + shuffle[ cur_epoch ] >= es ) permpos[i] = i - es + shuffle[cur_epoch ];
// // 		    else permpos[i] = i + shuffle[cur_epoch ];
// // 		    //std::cout << "pp\t" << cur_epoch << "\t" << i << "\t" << permpos[i] << "\t" << off << " " << shuffle[cur_epoch ] << "\n";
// // 		  }
// // 	      }
// // 	    else  // otherwise whole trace perm
// // 	      j = CRandom::rand( npoints );
	    

// // 	    for (int i=0;i<npoints;i++)
// // 	      {
		
// // 		bin( (*ph1)[i] , &b1 , bb , nbins );
		
// // 		if ( es ) // within-epoch shuffle
// // 		  j = permpos[i];
// // 		else // whole trace shuffle
// // 		  {
// // 		    ++j;
// // 		    if ( j == npoints ) j = 0;
// // 		  }

// // 		bin( (*ph2)[j] , &b2 , bb , nbins );

// // 		// increment by one unit (default, unweighted)
// // 		perm[ b1 ][ b2 ]++;
// // 	      }

	    

// // 	    // 
// // 	    // summary stat for surrogate 
// // 	    //

// // 	    double perm_stat = test_uniform( perm );
	    
// // 	    perm_stats.push_back( perm_stat );

// // 	    if ( perm_stat >= obs_stat ) ++emp_stat;

// // 	    //
// // 	    // accumulate point level pv/z file
// // 	    //

// // 	    for (int b1=0;b1<nbins;b1++)
// // 	      for (int b2=0;b2<nbins;b2++)
// // 		{
// // 		  if ( perm[b1][b2] >= obs[b1][b2] ) pv[b1][b2]++;
// // 		  z[b1][b2] += perm[b1][b2];
// // 		  z2[b1][b2] += perm[b1][b2] * perm[b1][b2];
// // 		}
	    
// // 	    // next replicate
	  
// // 	  }
	

// // 	double z_stat_mean = MiscMath::mean( perm_stats );
// // 	double z_stat_sd   = MiscMath::sdev( perm_stats );
// // 	double z_stat = ( obs_stat - z_stat_mean ) / z_stat_sd ; 
	
// // 	std::cout << f1[i1].first << "-" << f1[i1].second << "\t"
// // 		  << f2[i2].first << "-" << f2[i2].second << "\t"

// // 		  << obs_stat << "\t"
// // 		  << z_stat_mean << "\t"
// // 		  << z_stat_sd << "\t" 
// // 		  << z_stat << "\t"
	  
// // 		  << (emp_stat+1)/double(nreps+1) << "\n";



// // 	if ( true )
// // 	  {
	    
// // 	    for (int b1=0;b1<nbins;b1++)
// // 	      for (int b2=0;b2<nbins;b2++)
// // 		{
// // 		  //std::cout << "ha " << z[b1][b2] << "  " << (double)nreps << "\n";
// // 		  double zmean = z[b1][b2]/(double)nreps;
// // 		  double zsd   = sqrt(   z2[b1][b2]/(double)nreps  - zmean * zmean );
		  
// // 		  std::cout << "res " << b1 << " " << b2 << " " 
// // 			    << obs[b1][b2] << " " 		
// // 			    << (pv[b1][b2] + 1 )/double(nreps+1) << " " 
// // 			    << zmean << " " 
// // 			    << zsd << " " 
// // 			    << (obs[b1][b2] - zmean ) / zsd << "\n";
// // 		}
// // 	  }



// // 	//
// // 	// Next pair of frequencies
// // 	//

// //       }
  
// // }

// // bool phsyn_t::bin( double d , int * b , const std::vector<double> & th , const int nbins )
// // {
// //   // expect 0 to 360 
// //   if ( d < 0 || d > 360 ) return false;
// //   if ( *b < 0 || *b >= nbins ) return false;
  
// //   // bin 0 is 0 degree
// //   // last bin is 360 
  
// //   while ( 1 ) 
// //     {      
// //       //      std::cerr << "b = " << *b << " " << th[*b] << " " << th[ (*b ) + 1]  << "\t" << d << "\n";

// //       if ( *b == nbins - 1 )
// // 	{
// // 	  if ( d >= th[ *b ] ) return true;
// // 	  *b = 0;
// // 	}
      
// //       if ( d >= th[*b] && d < th[ (*b ) + 1] ) return true;      

// //       ++(*b);

// //       if ( *b == nbins ) *b = 0;
// //     }
  
// //   return false;
// // }


 
// // double phsyn_t::test_uniform( const std::vector<std::vector<double> > & m )
// // {
// //   // stat = ( O - E )^2 

// //   const int bs = m.size();
  
// //   std::vector<double> rows( bs , 0 );
// //   std::vector<double> cols( bs , 0 );
// //   double tot = 0;
// //   for (int b1=0;b1<bs;b1++)
// //     for (int b2=0;b2<bs;b2++)
// //       {
// // 	rows[b1] += m[b1][b2];
// // 	cols[b2] += m[b1][b2];
// // 	tot += m[b1][b2];
// //       }

// //   double stat = 0; 
    
// //   for (int b1=0;b1<bs;b1++)
// //     for (int b2=0;b2<bs;b2++)
// //       {
// // 	double exp = ( rows[b1] * cols[b2] ) / tot;
// // 	stat += ( m[b1][b2] - exp ) * ( m[b1][b2] - exp );
// //       }
// //   return stat;
    
// // }
