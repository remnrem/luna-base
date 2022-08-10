
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

#include "dsp/lat.h"
#include "eval.h"
#include "edf/edf.h"
#include "db/db.h"
#include "stats/statistics.h"

extern logger_t logger;
extern writer_t writer;

lat_t::lat_t( edf_t & edf , param_t & param )
{

  //
  // Get cache details (expecting epoch level power)
  //

  const std::string cache_name = param.requires( "cache" );

  // by default, variable name is PSD;  allow this to be altered, e.g.
  // if passing power results in from IRASA with PER and APER components
  // then we can pick which here
  
  const std::string cache_var = param.has( "cache-var" ) ? param.value( "cache-var" ) : "PSD" ;
  
  if ( ! edf.timeline.cache.has_num( cache_name ) )
    Helper::halt( "cache not found for this individual: " + cache_name );
  
  cache_t<double> * cache = edf.timeline.cache.find_num( cache_name );
  
  
  //
  // Get sleep stages
  //

  const int ne = edf.timeline.first_epoch();
  
  edf.timeline.annotations.make_sleep_stage( edf.timeline );

  // passes any hypno-related parameters too;
  // false --> not verbose
  bool has_staging = edf.timeline.hypnogram.construct( &(edf.timeline) , param , false );
  
  if ( (!has_staging)  || ne != edf.timeline.hypnogram.stages.size() )
    Helper::halt( "problem extracting stage information for trainer" );

  S.resize( ne , ASYMM_SS_IGNORE );
  
  for (int e=0; e<ne; e++)
    {      
      if ( edf.timeline.hypnogram.stages[ e ] == WAKE )
        S[e] = ASYMM_SS_WAKE ; 
      else if ( edf.timeline.hypnogram.stages[ e ] == NREM1 )
        S[e] = ASYMM_SS_NREM ; 
      else if ( edf.timeline.hypnogram.stages[ e ] == NREM2 )
        S[e] = ASYMM_SS_NREM ; 
      else if ( edf.timeline.hypnogram.stages[ e ] == NREM3 )
        S[e] = ASYMM_SS_NREM ; 
      else if ( edf.timeline.hypnogram.stages[ e ] == NREM4 )
        S[e] = ASYMM_SS_NREM ; 
      else if ( edf.timeline.hypnogram.stages[ e ] == REM )
        S[e] = ASYMM_SS_REM ; 
    }
  
  
  //
  // Requires pairings
  //

  left = param.strvector( "left" );
  right = param.strvector( "right" );
  
  if ( left.size() == 0 || left.size() != right.size() )
    Helper::halt( "expecting equal left= and right= channel sets" );

  sum_parts.clear();
  std::set<std::string> req_chs;

  // left
  for (int i=0; i<left.size(); i++)
    {
      std::vector<std::string> tok = Helper::parse( left[i] , "+" );
      for (int j=0; j<tok.size(); j++)
	{
	  sum_parts[ left[i] ].insert( tok[j] ) ;
	  req_chs.insert( tok[j] );
	}
    }

  // right
  for (int i=0; i<right.size(); i++)
    {
      std::vector<std::string> tok = Helper::parse( right[i] , "+" );
      for (int j=0; j<tok.size(); j++)
	{
	  sum_parts[ right[i] ].insert( tok[j] ) ;
	  req_chs.insert( tok[j] );
	}
    }

  //
  // Misc args
  //

  epoch_level_output = param.has( "epoch" );
  
  //
  // Extract and map power 
  //
  
  std::set<ckey_t> ckeys = cache->keys( cache_var );
  
  if ( ckeys.size() == 0 )
    Helper::halt( "no cache entries found for " + cache_name + " :: " + cache_var );
  
  // main data stores:
  
  // freq -> epoch -> channel -> power
  f2e2ch2psd.clear();
  
  // band -> epoch -> channel -> power
  b2e2ch2psd.clear();

  //
  // cache'd epochs are 1-based (from strata outputs), but
  // here set to 0-based on first read
  //
  
  // iterate over keys: only consider cache_var w/ CH, E, and (B or F)
  
  std::set<ckey_t>::const_iterator cc = ckeys.begin();    
  
  while ( cc != ckeys.end() )
    {
      // no channel
      if ( cc->stratum.find( globals::signal_strat ) == cc->stratum.end() )
	{
	  ++cc;
	  continue;
	}

      // no epoch
      if ( cc->stratum.find( globals::epoch_strat ) == cc->stratum.end() )
	{
	  ++cc;
	  continue;
	}

      // no F or B?
      const bool is_band =  cc->stratum.find( globals::band_strat ) != cc->stratum.end();
      const bool is_freq =  cc->stratum.find( globals::freq_strat ) != cc->stratum.end();
      if ( ! ( is_band || is_freq ) )
	{
	  ++cc;
	  continue;
	}
      
      //
      // If here, this is a valid entry: expecting just a single value
      //
      
      std::vector<double> cx = cache->fetch( *cc );

      if ( cx.size() != 1 )
	Helper::halt( "internal error in lat_t:: expecting scalar cache values" );
      
      const double psd = cx[0];

      int epoch1 = -1;
      double freq = -1;
      std::string band = "";
      std::string ch = cc->stratum.find( globals::signal_strat )->second;
      
      if ( ! Helper::str2int( cc->stratum.find( globals::epoch_strat )->second , &epoch1 ) )
	Helper::halt( "internal error with epoch encoding" );

      if ( is_freq )
	{
	  if ( ! Helper::str2dbl( cc->stratum.find( globals::freq_strat )->second , &freq ) )
	    Helper::halt( "internal error with epoch encoding" );
	}
      else
	{
	  band = cc->stratum.find( globals::band_strat )->second;
	}

      // 
      // recode as 0-based epoch count
      //
      
      const int epoch = epoch1 - 1;
      
      if ( is_band )
	b2e2ch2psd[ band ][ epoch ][ ch ] = psd ; 
      else
	f2e2ch2psd[ freq ][ epoch ][ ch ] = psd ; 
      
      //
      // next stratum
      //
      
      ++cc;
    }

  
  //
  // Check for squared off data
  //

  const bool any_freqs = f2e2ch2psd.size() != 0 ;
  const bool any_bands = b2e2ch2psd.size() != 0 ;

  // *assume* that bands and freqs have same epoch list, if both present...

  std::set<int> all_epochs;
  
  if ( any_freqs )
    {
      int data_points = 0;
      
      std::set<double> fset;
      std::set<int> eset;
      std::set<std::string> cset;
      
      auto ff = f2e2ch2psd.begin();
      while ( ff != f2e2ch2psd.end() )
	{
	  fset.insert( ff->first );
	  auto ee = ff->second.begin();
	  while ( ee != ff->second.end() )
	    {
	      eset.insert( ee->first );
	      auto cc = ee->second.begin();
	      while ( cc != ee->second.end() )
		{
		  cset.insert( cc->first );
		  ++data_points;
		  ++cc;
		}	  
	      ++ee;
	    }
	  ++ff;
	}

      //
      // check that data are squared off
      //

      if ( fset.size() * cset.size() * eset.size() != data_points )
	{
	  logger << "  expecting " << data_points << " data points (given "
		 << fset.size() << " freq bins x " << cset.size() << " chs x " << eset.size() << " epochs)\n"
		 << "  but observed only " << data_points << "\n";	    
	  Helper::halt( "ASYMM requires squared data" );
	}

      //
      // check that all necessary channels are present
      //

      auto cc = req_chs.begin();
      while ( cc != req_chs.end() )
	{
	  if ( cset.find( *cc ) == cset.end() )
	    Helper::halt( "could not find requested channel " + *cc + " in the cache" );
	  ++cc;
	}
      
      //
      // all good
      //

      logger << "  for " << fset.size() << " frequency bins, cached data on "
	     << cset.size() << " channels, " << eset.size() << " epochs\n";

      // track
      all_epochs = eset;
      
    }

  
  //
  // bands
  //
  
  if ( any_bands )
    {
      int data_points = 0;
      
      std::set<std::string> bset;
      std::set<int> eset;
      std::set<std::string> cset;
      
      auto bb = b2e2ch2psd.begin();
      while ( bb != b2e2ch2psd.end() )
	{
	  bset.insert( bb->first );
	  auto ee = bb->second.begin();
	  while ( ee != bb->second.end() )
	    {
	      eset.insert( ee->first );
	      auto cc = ee->second.begin();
	      while ( cc != ee->second.end() )
		{
		  cset.insert( cc->first );
		  ++data_points;
		  ++cc;
		}	  
	      ++ee;
	    }
	  ++bb;
	}

      //
      // Check squared off data
      //

      if ( bset.size() * cset.size() * eset.size() != data_points )
	{
	  logger << "  expecting " << data_points << " data points (given "
		 << bset.size() << " bands x " << cset.size() << " chs x " << eset.size() << " epochs)\n"
		 << "  but observed only " << data_points << "\n";	    
	  Helper::halt( "ASYMM requires squared data" );
	}

      //
      // check that all necessary channels are present
      //

      auto cc = req_chs.begin();
      while ( cc != req_chs.end() )
	{
	  if ( cset.find( *cc ) == cset.end() )
	    Helper::halt( "could not find requested channel " + *cc + " in the cache" );
	  ++cc;
	}

      if ( any_freqs && all_epochs.size() != eset.size() )
	Helper::halt( "internal error in ASYMM: mismatch of epoch sizes between bands and freqs" );
      else
	all_epochs = eset;

      //
      // All good
      //

      logger << "  for " << bset.size() << " bands, cached data on "
	     << cset.size() << " channels, " << eset.size() << " epochs\n";
    }
  
  
  //
  // Splice out necessary stages
  //

  std::vector<stg_t> S2 = S;
  
  S.clear();
  E.clear();
  auto ee = all_epochs.begin();
  while ( ee != all_epochs.end() )
    {
      // we're now working w/ 0-based epochs
      if ( *ee >= ne )
	{
	  logger << " expected (max) = " << ne << "\n";
	  logger << " observed (0-based index) " << *ee << "\n";	    
	  Helper::halt( "unexpected epoch number found, greater than implied ne" );
	}
      S.push_back( S2[*ee] );
      E.push_back( *ee );
      ++ee;
    }

  if ( S.size() < S2.size() )
    logger << "  spliced out " << S.size() <<" of " << S2.size() << " stages\n";

  // do the actual work
  proc( edf , param ); 
  
}
  

void lat_t::proc( edf_t & edf , param_t & param )
{

  //
  // Create sleep cycles  
  //

  const int ne = S.size();  

  // epoch level output 
  const bool verbose = false;

  bool okay = edf.timeline.hypnogram.construct( &edf.timeline , param , verbose ); 
  if ( ! okay ) Helper::halt( "problem constructing the hypnogram" );

  const bool epoch_lvl_output = false;
  const std::string eannot = "";  
  // cycle annotation --> epoch annots
  const std::string cycle_annot = "NREMC";
  
  // add annotations (this also generates some output, output...
  edf.timeline.hypnogram.output( verbose, epoch_lvl_output , eannot , cycle_annot );

  // look up to six cycles
  const bool has_cycles = edf.timeline.epoch_annotation( "_NREMC_1" ) 
    || edf.timeline.epoch_annotation( "_NREMC_2" )
    || edf.timeline.epoch_annotation( "_NREMC_3" )
    || edf.timeline.epoch_annotation( "_NREMC_4" )
    || edf.timeline.epoch_annotation( "_NREMC_5" )
    || edf.timeline.epoch_annotation( "_NREMC_6" );

  C.resize( ne , 0 );

  num_cycles = 0;
  
  if ( has_cycles )
    {
      
      // get current epoch encoding
      
      const int ne2 = edf.timeline.first_epoch();

      if ( ne != ne2 )
	Helper::halt( "EDF has been restructed prior to ASYMM... epochs encoding off" ); 

      std::vector<int> epochs;
      while ( 1 )
        {
	  int epoch = edf.timeline.next_epoch();
          if ( epoch == -1 ) break;
	  epochs.push_back( epoch );
	}
      
      for (int e=0;e<ne; e++)
	{	  
	  // up to 6 cycles
	  if      ( edf.timeline.epoch_annotation( "_NREMC_1" , epochs[e] ) ) C[e] = 1;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_2" , epochs[e] ) ) C[e] = 2;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_3" , epochs[e] ) ) C[e] = 3;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_4" , epochs[e] ) ) C[e] = 4;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_5" , epochs[e] ) ) C[e] = 5;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_6" , epochs[e] ) ) C[e] = 6;

	  if ( C[e] > num_cycles ) num_cycles = C[e];

	}
    }

  logger << "  detected " << num_cycles << " NREM cycles\n";
  
  //
  // NR/R transitions
  //

  T_R2NR.resize( ne, 0 );
  T_NR2R.resize( ne, 0 );

  // default, 3 mins (6 epochs) each side of a transition 
  const int e_window = param.has( "trans" ) ? param.requires_int( "trans" ) : 6;

  //
  // NREM -> REM 
  //

  for (int e=1; e<ne; e++)
    {
      if ( S[e] == ASYMM_SS_REM && S[e-1] == ASYMM_SS_NREM )
	{
	  bool okay = true;
	  int cnt = 1;
	  int idx = e-1;
	  while ( 1 )
	    {
	      --idx;
	      if ( idx == -1 ) { okay = false; break; }
	      if ( S[idx] != ASYMM_SS_NREM ) { okay = false; break; }
	      ++cnt;
	      if ( cnt == e_window ) break;
	    }

	  // left okay
	  if ( okay )
	    {

	      int cnt = 1;
	      int idx = e;
	      while	( 1 )
		{
		  ++idx;
		  if ( idx == ne ) { okay = false; break; }
		  if ( S[idx] != ASYMM_SS_REM ) { okay = false; break; }
		  ++cnt;
		  if ( cnt == e_window ) break;
		}

	      // both sides okay?
	      if ( okay )
		{
		  const int start = e - e_window;
		  const int stop  = e + e_window - 1;
		  for (int i=start; i<=stop; i++)
		    {
		      // encoding: ..., -3, -2, -1, +1, +2, +3 ...
		      // 0 means no transition
		      if ( i < e )
			T_NR2R[ i ] = i - e;
		      else
			T_NR2R[ i ] = i - e + 1; 
		    }
		}
	    }
	  
	}
    } // next epoch

  //
  // REM -> NREM transitions
  //

  for (int e=1; e<ne; e++)
    {
      if ( S[e] == ASYMM_SS_NREM && S[e-1] == ASYMM_SS_REM )
	{
	  bool okay = true;
	  int cnt = 1;
	  int idx = e-1;
	  while ( 1 )
	    {
	      --idx;
	      if ( idx == -1 ) { okay = false; break; }
	      if ( S[idx] != ASYMM_SS_REM ) { okay = false; break; }
	      ++cnt;
	      if ( cnt == e_window ) break;
	    }

	  // left okay
	  if ( okay )
	    {

	      int cnt = 1;
	      int idx = e;
	      while	( 1 )
		{
		  ++idx;
		  if ( idx == ne ) { okay = false; break; }
		  if ( S[idx] != ASYMM_SS_NREM ) { okay = false; break; }
		  ++cnt;
		  if ( cnt == e_window ) break;
		}
	      
	      // both sides okay?
	      if ( okay )
		{
		  const int start = e - e_window;
		  const int stop  = e + e_window - 1;
		  for (int i=start; i<=stop; i++)
		    {
		      // encoding: ..., -3, -2, -1, +1, +2, +3 ...
		      // 0 means no transition
		      if ( i < e )
			T_R2NR[ i ] = i - e;
		      else
			T_R2NR[ i ] = i - e + 1; 
		    }
		}
	    }
	  
	}
    } // next epoch

  //
  // Check
  //

  if ( 1  )
    {
      for (int e=0; e<ne; e++)
	{
	  logger << E[e]+1 << "\t"
		 << S[e] << "\t"
		 << C[e] << "\t"
		 << T_NR2R[e] << "\t"
		 << T_R2NR[e] << "\n";
	}
    }
  
  //
  // Iterate over bands/bins and channel pairs to do ASYMM analysis
  //

  const bool any_bands = b2e2ch2psd.size() != 0;
  const bool any_freqs = f2e2ch2psd.size() != 0;

  //
  // bands
  //

  if ( any_bands )
    {
      
      auto bb = b2e2ch2psd.begin();
      while ( bb != b2e2ch2psd.end() )
        {
	  
	  writer.level( bb->first , globals::band_strat );

	  //
	  // iterate 
	  //

	  for (int p = 0; p < left.size(); p++)
	    {

	      writer.level( left[p] + ".." + right[p] , "CHS" );
	      
	      // get signals / signal groups  [ treat as group, even if N=1 ] 
	      const std::set<std::string> & sigs1 = sum_parts[ left[p] ];
	      const std::set<std::string> & sigs2 = sum_parts[ right[p] ];
	      if ( sigs1.size() == 0 || sigs2.size() == 0 )
		Helper::halt( "internal error in finding sigs" );

	      std::vector<double> Lacc, Racc;
	      
	      // iterate of epochs
	      for (int e=0; e<ne; e++)
		{
		  const int epoch = E[e];
		  //const stg_t stg = S[e];
		  
		  std::map<std::string,double> & D = bb->second[ epoch ];
		  
		  // left
		  double l = 0;
		  auto ss1 = sigs1.begin();
		  while ( ss1 != sigs1.end() )
		    {
		      l += D[ *ss1 ];
		      ++ss1;
		    }
		  
		  // right
		  double r = 0;
		  auto ss2 = sigs2.begin();
		  while ( ss2 != sigs2.end() )
		    {
		      r += D[ *ss2 ];
		      ++ss2;
		    }
		  
		  // accum
		  
		  Lacc.push_back( l );
		  Racc.push_back( r );
		  
		}
	      
	      //
	      // analysis
	      //
	      
	      lat_results_t res = analyse( Lacc , Racc );

	      
	      //
	      // next pair
	      //
	    }

	  writer.unlevel( "CHS" );
	  
	  //
	  // done w/ this band
	  //
	  
	  ++bb;
        }
      writer.unlevel( globals::band_strat );
      
    }


  //
  // Frequency-bin based analysis
  //

  if ( any_freqs )
    {
      
      auto ff = f2e2ch2psd.begin();
      while ( ff != f2e2ch2psd.end() )
        {
	  
	  writer.level( ff->first , globals::freq_strat );

	  //
	  // iterate 
	  //

	  for (int p = 0; p < left.size(); p++)
	    {
	      
	      writer.level( left[p] + ".." + right[p] , "CHS" );
	      
	      // get signals / signal groups  [ treat as group, even if N=1 ] 
	      const std::set<std::string> & sigs1 = sum_parts[ left[p] ];
	      const std::set<std::string> & sigs2 = sum_parts[ right[p] ];
	      if ( sigs1.size() == 0 || sigs2.size() == 0 )
		Helper::halt( "internal error in finding sigs" );
	      
	      std::vector<double> Lacc, Racc;
	      
	      // iterate of epochs
	      for (int e=0; e<ne; e++)
		{
		  const int epoch = E[e];
		  //const stg_t stg = S[e];
		  
		  std::map<std::string,double> & D = ff->second[ epoch ];
		  
		  // left
		  double l = 0;
		  auto ss1 = sigs1.begin();
		  while ( ss1 != sigs1.end() )
		    {
		      l += D[ *ss1 ];
		      ++ss1;
		    }
		  
		  // right
		  double r = 0;
		  auto ss2 = sigs2.begin();
		  while ( ss2 != sigs2.end() )
		    {
		      r += D[ *ss2 ];
		      ++ss2;
		    }
		  
		  // accum
		  
		  Lacc.push_back( l );
		  Racc.push_back( r );
		  
		}
	      
	      //
	      // analysis
	      //
	      
	      lat_results_t res = analyse( Lacc , Racc );

	      
	      //
	      // next pair
	      //
	    }
	  writer.unlevel( "CHS" );
	  
	  //
	  // done w/ this freq
	  //
	  
	  ++ff;
        }
      writer.unlevel( globals::freq_strat );
      
    }

  
}


lat_results_t lat_t::analyse( const std::vector<double> & L ,
			      const std::vector<double> & R )
{
  
  // S, C, T_* and E will match length of L and R
  lat_results_t res;
    
  const int ne = L.size();
  
  //
  // analysis parameters
  // 

  // min. # of REM epochs per cycle
  const int req_rem_epochs = 10;
  
  // min. # of NREM epochs per cycle
  const int req_nrem_epochs = 10;

  // leading/trailing NREM search interval
  const int flanking_epochs_search = 100; 

  // max # of flanking NREM epochs to take (of /100)
  const int flanking_epochs_max = 40;

  // min # of flanking NREM epochs
  const int flanking_epochs_min = 25;

  // min # of good NREM following outlier removal (within region)
  const int req_good_nrem = 10;
  
  // min # of good REM following outlier removal (within region)
  const int req_good_rem = 10;
  
  // outlier limit range
  const double outlier_ratio = 2.0;						

  // outlier th
  const double psd_th = 3;

  //
  // As we can see major changes in L/R ratios across the night, this analysis is focussed to
  // matched NREM-REM-NREM regions 
  //  

  //
  // Track overall outliers 
  //

  std::vector<bool> outlier( ne , false );

  //
  // get log2(L/R)
  //

  const double ASYMM_EPS = 1e-6 ;
  
  std::vector<double> log2lr( ne , 0 );
  for (int e=0; e<ne; e++)
    {
      if ( L[e] < ASYMM_EPS || R[e] < ASYMM_EPS )
	outlier[e] = true;
      else
	{
	  log2lr[e] = log2( L[e] / R[e] ) ;
	  if ( log2lr[e] < -outlier_ratio || log2lr[e] > outlier_ratio )
	    outlier[e] = true;
	}
    }

  //
  // Raw, epoch-level output 
  //

  if ( epoch_level_output )
    {
      for (int e=0; e<ne; e++)
	{
	  writer.epoch( E[e] + 1 );
	  writer.value( "L" , L[e] );
	  writer.value( "R" , R[e] );
	  writer.value( "LR" , log2lr[e] );
	  writer.value( "OUT" , (int)outlier[e] );
	  writer.value( "C" , C[e] );

	  std::string ss = "?";
	  if ( S[e] == ASYMM_SS_WAKE ) ss = "W";
	  else if ( S[e] == ASYMM_SS_REM ) ss = "R";
	  else if ( S[e] == ASYMM_SS_NREM ) ss = "NR";
	  writer.value( "SS" , ss );
	}
      writer.unepoch();
    }
  
  
  //
  // Cycle based analysis 
  //

  for (int c=1; c<=num_cycles; c++)
    {
      std::vector<int> rems, nrems;
      for (int e=0; e<ne; e++)
	{
	  if ( C[e] > c ) break;
	  if ( outlier[e] ) continue;
	  if ( C[e] == c )
	    {
	      if ( S[e] == ASYMM_SS_REM ) rems.push_back(e);
	      else if ( S[e] == ASYMM_SS_NREM ) nrems.push_back(e);
	    }
	}
      
      logger << "  cycle " << c << ", found " << nrems.size() << " and " << rems.size() << "NR/R epochs\n";

      //
      // skip if not enough NREM or REM epochs
      //

      if ( rems.size() < req_rem_epochs ) continue;
      if ( nrems.size() < req_nrem_epochs ) continue;
      
      // get REM ranges
      
      const int rem_start = rems[0];
      const int rem_stop  = rems[ rems.size() - 1 ];

      // define NREM search space
      const int lwr_epoch = rem_start - flanking_epochs_search < 0 ? 0 : rem_start - flanking_epochs_search ;
      const int upr_epoch = rem_stop  + flanking_epochs_search >= ne ? ne - 1 : rem_stop  + flanking_epochs_search ;


      // get flanking NREM in here:
      nrems.clear();

      // leading NREM
      int cnt = 0;      
      for (int e=rem_start-1; e >= lwr_epoch ; e--)
	{
	  if ( outlier[e] ) continue;
	  if ( S[e] == ASYMM_SS_NREM )
	    {
	      ++cnt;
	      if ( cnt > flanking_epochs_max ) break;
	      nrems.push_back( e );
	    }
	}
      
      // trailing NREM
      cnt = 0;      
      for (int e=rem_stop+1; e <= upr_epoch ; e++)
	{
	  if ( outlier[e] ) continue;
	  if ( S[e] == ASYMM_SS_NREM )
	    {
	      ++cnt;
	      if ( cnt > flanking_epochs_max ) break;
	      nrems.push_back( e );
	    }
	}

      // enough NREM overall? if not, skip to next cycle
      if ( nrems.size() < flanking_epochs_min )
	continue;
      
      //
      // Flag outliers (in original channels, log-scaled)
      //

      std::vector<double> rem_l, rem_r, nrem_l, nrem_r;
      
      for (int i=0; i<rems.size(); i++)
	{
	  rem_l.push_back( log( L[ rems[i] ] ) );
	  rem_r.push_back( log( R[ rems[i] ] ) );
	}
      
      for (int i=0; i<nrems.size(); i++)
	{
	  nrem_l.push_back( log( L[ nrems[i] ] ) );
	  nrem_r.push_back( log( R[ nrems[i] ] ) );
	}
      
      const double rem_l_mean = MiscMath::mean( rem_l );
      const double rem_r_mean = MiscMath::mean( rem_r );
      
      const double rem_l_sd = MiscMath::sdev( rem_l , rem_l_mean );
      const double rem_r_sd = MiscMath::sdev( rem_r , rem_r_mean );

      for (int i=0; i<rems.size(); i++)
	{
	  if ( rem_l[i] < rem_l_mean - psd_th * rem_l_sd || rem_l[i] > rem_l_mean + psd_th * rem_l_sd )
	    outlier[ rems[i] ] = true;

	  if ( rem_r[i] < rem_r_mean - psd_th * rem_r_sd || rem_r[i] > rem_r_mean + psd_th * rem_r_sd )
	    outlier[ rems[i] ] = true;
	}

      const double nrem_l_mean = MiscMath::mean( nrem_l );
      const double nrem_r_mean = MiscMath::mean( nrem_r );

      const double nrem_l_sd = MiscMath::sdev( nrem_l , nrem_l_mean );
      const double nrem_r_sd = MiscMath::sdev( nrem_r , nrem_r_mean );

      for (int i=0; i<nrems.size(); i++)
	{
	  if ( nrem_l[i] < nrem_l_mean - psd_th * nrem_l_sd || nrem_l[i] > nrem_l_mean + psd_th * nrem_l_sd )
	    outlier[ nrems[i] ] = true;
	  
	  if ( nrem_r[i] < nrem_r_mean - psd_th * nrem_r_sd || nrem_r[i] > nrem_r_mean + psd_th * nrem_r_sd )
	    outlier[ nrems[i] ] = true;
	}

      
      //
      // Compile valid ratios
      //

      std::vector<double> rem_lr, nrem_lr;
      
      for (int i=0; i<rems.size(); i++)
	if ( ! outlier[rems[i]] )
	  rem_lr.push_back( log2lr[ rems[i] ] );

      for (int i=0; i<nrems.size(); i++)
	if ( ! outlier[nrems[i]] )
	  nrem_lr.push_back( log2lr[ nrems[i] ] );
      
           
      //
      // Normalize by NREM mean/SD
      //
      
      const double rem_mean = MiscMath::mean( rem_lr );
      const double nrem_mean = MiscMath::mean( nrem_lr );

      const double rem_sd = MiscMath::sdev( rem_lr , rem_mean );
      const double nrem_sd = MiscMath::sdev( nrem_lr , nrem_mean );

      std::vector<double> zrem_lr( rem_lr.size() );
      std::vector<double> znrem_lr( nrem_lr.size() );
      
      for (int i=0; i<rem_lr.size(); i++)
	zrem_lr[i] = ( rem_lr[i] - nrem_mean ) / nrem_sd ; 
      
      for (int i=0; i<nrem_lr.size(); i++)
        znrem_lr[i] = ( nrem_lr[i] - nrem_mean ) / nrem_sd ;
      
      
      // REM mean, in NREM SD units
      const double zrem_mean = MiscMath::mean( zrem_lr );
      const double zrem_sd = MiscMath::sdev( zrem_lr , zrem_mean );

      // should be 0 and 1, but just to check
      const double znrem_mean = MiscMath::mean( znrem_lr );
      const double znrem_sd = MiscMath::sdev( znrem_lr , znrem_mean );

      // t-test 

      double pvalue = 1;

      // Welch's t-test, unequal samples, unequal variances
      bool okay = Statistics::t_test( zrem_mean , zrem_sd * zrem_sd , zrem_lr.size() ,
				      znrem_mean , znrem_sd * znrem_sd , znrem_lr.size() ,
				      &pvalue );

      if ( !okay ) continue;
      
      writer.level( c, globals::cycle_strat );
      
      writer.value( "LR_REM" , rem_mean );
      writer.value( "LR_NREM" , nrem_mean );
      writer.value( "Z_REM" , zrem_mean );
      writer.value( "Z_NREM" , znrem_mean );
      writer.value( "N_REM" , (int)zrem_lr.size() );
      writer.value( "N_NREM" , (int)znrem_lr.size() ); 
      writer.value( "P" , pvalue );
      writer.value( "Z" , zrem_mean * -log10( pvalue ) );
      
      //
      // next cycle
      //
    }
  writer.unlevel( globals::cycle_strat );


  

    
  return res;
}

