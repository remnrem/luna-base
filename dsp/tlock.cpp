
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

#include "dsp/tlock.h"
#include "timeline/timeline.h"
#include "stats/statistics.h"
#include "db/db.h"
#include "helper/logger.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "assoc/massoc.h"
#include "dsp/mt_spectrogram.h"

extern writer_t writer;
extern logger_t logger;


void dsptools::tlock( edf_t & edf , param_t & param )
{
  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );

  if ( signals.size() < 1 ) return;

  const int ns = signals.size();

  // check sample rates: must be similar across all 

  std::vector<double> Fs = edf.header.sampling_freq( signals );
  

  for (int s=1;s<Fs.size();s++) 
    if ( Fs[s] != Fs[0] ) 
      Helper::halt( "sample rates must be similar across signals for TLOCK" );

  //
  // Options for how to handle inputs:  take logs before means, or treat as 
  // angles (circular values);  here, rather than calculate a circular mean, or ITPC
  // etc, for now output a binned histogram of counts (column nornalized) by angle bins
  // e.g. 18 bins of 20 degrees
  //

  bool take_log = param.has( "tolog" ) ;

  int angle_bins = param.has( "phase" ) ? param.requires_int( "phase" ) : 0 ;

  if ( take_log && angle_bins ) Helper::halt( "cannot specify both tolog and phase" );


  
  //
  // get time-points: expecting a cache in sample-point units (i.e. 
  // to ensure same number of points per window; i.e. this is also why
  // all signals must have the same SR
  //


  //
  // only self-channels
  //

  // by default, gives sync'ed values for all channels in sig
  // even if the cache has a CH strata that is different (i.e. seeded on sCH = "C3" )
  //  to only give output for same channels, add 'same-channel=T'

  const bool same_channel = param.yesno( "same-channel" );

  // also allow to match e.g. 'same-channel'  if seed = C3
  //   same-channel=T  channel-postfix=_SIGMA
  //  will match sig = C3 /and/ C3_SIGMA for this seed

  const std::string channel_postfix = param.has( "channel-postfix" ) ? param.value( "channel-postfix" ) : ""; 

  if ( channel_postfix != "" && ! same_channel )
    Helper::halt( "cannot specify channel-postfix without same-channel=T" );
  
  //
  // output
  // 
  
  bool verbose = param.has( "verbose" );

  bool to_massoc = param.has( "export" );
  
  //
  // get window
  //
  
  double half_window = param.requires_dbl( "w" );
  if ( half_window <= 0 ) Helper::halt( "w must be a positive number" );
  
  int half_points = half_window * Fs[0] ;
 
  // ensure an nice multiple of sample rate
  half_window = half_points / Fs[0];
  half_points = half_window * Fs[0] ;  
  int points = 1 + 2 * half_points;
  
  std::vector<double> t;
  // nb. fudge for floating point issues (allowing tenth of inc for stop 'w')
  const double inc = 1.0/Fs[0];
  for ( double w = -half_window ; w <= half_window + inc/10.0   ; w += inc ) 
    t.push_back(w);

  if ( t.size() != points ) 
    Helper::halt( "internal error constructing window:" + Helper::int2str( points ) + " vs " + Helper::int2str( (int) t.size() )  );

  //
  // Normalisation: e.g. 0.2 means 20% + 20% of window, i.e. skip middle 60%
  //                defauly = 0 , no normalisation
  
  int norm_points = points * ( param.has( "np" ) ? param.requires_dbl( "np" ) : 0 );
  
  if ( norm_points < 0 || norm_points > points / 2 ) 
    Helper::halt( "expecting np between 0 and 0.5" );

  //
  // Spectrogram
  //

  const bool spectrogram = param.has( "spectrogram" );
  
  const double mtm_nw = param.has( "nw" ) ? param.requires_dbl( "nw" ) : 3 ;
  const int    mtm_t  = param.has( "t" ) ? param.requires_int( "t" ) : 5 ;
  const double mtm_seg = param.has( "segment-sec" ) ? param.requires_dbl( "segment-sec" ) : half_window ;
  const double mtm_step = param.has( "segment-inc" ) ? param.requires_dbl( "segment-inc" ) : mtm_seg / 4.0;
  const double mtm_fmin = param.has( "f-lwr" ) ? param.requires_int( "f-lwr" ) : 1 ;
  const double mtm_fmax = param.has( "f-upr" ) ? param.requires_int( "f-upr" ) : 30 ;
  const bool mtm_dB = param.has( "dB" ) ? true : param.yesno( "dB" );
  const bool mtm_center = param.has( "center" ) ? true : param.yesno( "center" );

  
  //
  // Get seed sample-points from cache
  //

  std::string cache_name = param.requires( "cache" );

  if ( ! edf.timeline.cache.has_int( cache_name ) )
    Helper::halt( "cache not found for this individual: " + cache_name );

  cache_t<int> * cache = edf.timeline.cache.find_int( cache_name );

  //  cache->dump();

  //
  // Get any/all keys associated with the 'points' intrnal name cache
  //

  std::set<ckey_t> ckeys = cache->keys( "points" );

  std::set<ckey_t>::const_iterator cc = ckeys.begin();

  // can use CACHE command explicitly...
  //logger << cache->print();


  while ( cc != ckeys.end() )
    {
      
      int scnt = 0;

      std::vector<int> cx = cache->fetch( *cc );

      //
      // do we have a channel specification, and must this match?
      //

      std::string seed_channel = "";
      if ( same_channel && cc->stratum.find( globals::signal_strat ) != cc->stratum.end() )
	{
	  seed_channel = cc->stratum.find(  globals::signal_strat )->second;
	}
      
      //
      // add output stratifiers based on this key
      //

      std::stringstream sstr;
      std::map<std::string,std::string>::const_iterator ss = cc->stratum.begin();
      while ( ss != cc->stratum.end() )
	{
	  writer.level( ss->second , "s" + ss->first );
	  sstr << ss->first << "=" << ss->second << ";" ;
	  ++ss;
	}

      //
      // perform separately for each signal
      //

      for (int s=0; s<ns; s++)
	{

	  //
	  // track which intervals are kept (i.e. to crossing discontinuities)
	  // (this is only used in MASSOC mode)
	  //
	  
	  std::vector<int> eidx_base1;

	  int cnt_valid_intervals = 0;
	  
	  writer.level( signals.label(s) , globals::signal_strat );

	  // skip this channel?

	  if ( same_channel && seed_channel != "" )
	    {
	      if ( seed_channel != signals.label(s) &&
		   seed_channel + channel_postfix != signals.label(s) )
		continue;
	    }
	  
	  // get data and TP information 
	  
	  slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );
	  
	  const std::vector<double> * d = slice.pdata();
	  
	  const std::vector<uint64_t> * tp = slice.ptimepoints();
	  
	  // build up TLOCK matrix
	  // nb. we still need X constructed in MASSOC mode, even
	  // if we do not write to outptu ('export' versus 'verbose')
	  
	  tlock_t tlock(t, norm_points );
	  tlock.verbose = verbose || to_massoc || spectrogram ; 

	  // loop time-points to sync/lock on
	  
	  for ( int i=0; i<cx.size(); i++)
	    {
	      
	      //std::cout << "CX " << cx[i] << "\n";
	      
	      int lower = cx[i] - half_points;
	      int upper = cx[i] + half_points;
	      
	      //std::cout << " C1\n";
	      // ??TODO: check for discontunuities (EDF+)	  
	      //std::cout << " tp-> " << tp->size() << " " << Fs[0] << " " << lower << " " << upper << "\t";
	      
	      // interval out-of-range
	      if ( lower < 0 || upper > d->size() ) 
		{
		  //std::cout << "XXX2\n";			       
		  continue;
		}
	      
		  
	      if ( edf.timeline.discontinuity( *tp , Fs[0] , lower , upper ) )
		{
		  //std::cout << "XXX\n";
		  continue;
		}
	      //std::cout << "YES\n";
	      
	      // otherwise, add this interval to the tlock 
	      
	      tlock.add( d , lower , upper , take_log , angle_bins );
	      
	      ++cnt_valid_intervals;

	      if ( to_massoc )
		eidx_base1.push_back( i + 1 ); 
				     
	      // next interval
	    }
	  
	  // report number of interval considered/accepted
	  writer.value( "N" , cnt_valid_intervals );
	  writer.value( "N_ALL" , (int)cx.size() );
	  	  
	  logger << "  included " << cnt_valid_intervals
		 << " of " << (int)cx.size()
		 << " intervals for strata " << ++scnt << " "
		 << sstr.str() << " for channel " << signals.label(s) << "\n";
	
	  //
	  // Report either as phase angles (assuming radians), otherwise take the mean
	  //

	  if ( angle_bins ) 
	    {
	      Data::Matrix<double> angbin = tlock.angles();
	      
	      if ( angbin.dim1() != angle_bins || angbin.dim2() != t.size() )
		Helper::halt( "internal error in tlock_t()" );

	      for (int i=0; i<angle_bins; i++)
		{
		  writer.level( i+1 , "PHASE" );
		  for (int j=0; j<t.size(); j++)
		    {
		      writer.level( t[j] , "SEC" );		      
		      writer.value( "M" , angbin(i,j) );
		    }
		}
	      writer.unlevel( "PHASE" );
	      writer.unlevel( "SEC" );
	      
	    }

	  else if ( spectrogram )
	    {

	      logger << "  calculating mean MT spectrogram...\n";

	      Data::Vector<double> means = tlock.average( );
	      
	      mt_spectrogram_t mtm( tlock.X ,
				    Fs[0] ,
				    mtm_nw , mtm_t,
				    mtm_seg , mtm_step,
				    mtm_fmin, mtm_fmax,
				    mtm_dB,
				    mtm_center );

	      const int nf = mtm.frq.size();
	      const int nt = mtm.t.size();

	      // mean-centered, normalized versions (by row/freq and by col/time)
	      
	      // double min = 999;
	      // for (int i=0; i<nf; i++)
	      // 	for (int j=0; j<nt; j++)
	      // 	  if ( mtm.Z(i,j) < min ) min = mtm.Z(i,j);

	      //Data::Matrix<double> Z0 = mtm.Z;
	      
	      // for (int i=0; i<nf; i++)
	      // 	for (int j=0; j<nt; j++)
	      // 	  Z0(i,j) += min;

	      std::vector<double> row_min( nf );
	      std::vector<double> row_max( nf );
	      
	      for (int i=0; i<nf; i++)
		{
		  row_min[i] = mtm.Z(i,0);
		  row_max[i] = mtm.Z(i,0);
		  
		  for (int j=1; j<nt; j++)
		    {
		      if ( mtm.Z(i,j) < row_min[i] ) row_min[i] = mtm.Z(i,j);
		      if ( mtm.Z(i,j) > row_max[i] ) row_max[i] = mtm.Z(i,j);		      
		    }
		}
	      
	      // output
	      for (int i=0; i<nf; i++)
		{
		  writer.level( mtm.frq[i] , globals::freq_strat );
		  for (int j=0; j<nt; j++)
		    {
		      writer.level( mtm.t[j] , "SEC" );
		      writer.value( "PSD" , mtm.Z(i,j) );
		      // normalized by F (row) min/max
		      writer.value( "PSD_F" , ( mtm.Z(i,j) - row_min[i] ) / ( row_max[i] - row_min[i] ) );
		      
		      // variance in PSD
		      writer.value( "VAR" , mtm.ZZ(i,j) );
		    }
		  writer.unlevel( "SEC" );
		}
	      writer.unlevel( globals::freq_strat );
	    }
	 
	  //
	  // report summaries: for regular values, get the average
	  // and normalize (by edges)
	  //
	  
	  else
	    {
	      Data::Vector<double> means = tlock.average( );
	  
	      if ( means.size() != t.size() ) 
		{
		  logger << "  means.size() = " << means.size() << " " << t.size() << "\n";
		  Helper::halt( "internal error in tlock_t()" );
		}	      
	      
	      for (int i=0; i<means.size(); i++) 
		{
		  writer.level( t[i] , "SEC" );
		  writer.value( "M" , means[i] );
		}

	      writer.unlevel( "SEC" );
	    }

	  //
	  // Dump for MASSOC?
	  //

	  if ( to_massoc )
	    {
#ifdef HAS_LGBM
	      
	      // filename = indiv-id + strata
	      // ID =
	      // variables = T

	      const int nrow = tlock.X.dim2(); // note - obs in cols here
	      const int ncol = tlock.X.dim1();
	      
	      // row-IDs : strata + N
	      //  EDFID_C3_11_1
	      //  EDFID_C3_11_2 ...
	      //  etc
	      
	      std::vector<std::string> rowids( nrow );
	      std::vector<std::string> eids( nrow );
	      
	      std::string rowbase = signals.label(s);
	      std::map<std::string,std::string>::const_iterator ss = cc->stratum.begin();
	      while ( ss != cc->stratum.end() )
		{
		  // nb. here 'CH' means seed channel;  SIG means the readout
		  rowbase += "_" + ss->second ;
		  ++ss;
		}

	      if ( eidx_base1.size() != nrow )
		Helper::halt( "internal error in TLOCK w/ eidx_base1 size" );
	      
	      for (int i=0; i<nrow; i++)
		{
		  rowids[i] = rowbase ;
		  eids[i] = Helper::int2str( eidx_base1[i] );
		}
	      
	      // col IDs: simply 1, 2, 3, etc
	      std::vector<std::string> colids( ncol );
	      for (int i=0; i<ncol; i++)
		colids[i] = Helper::int2str( i+1 );
	      
	      // i.e. expecting something like export=path/
	      //      or export=path/root
	      //        filename will append root_ID_ROWBASE
	      const std::string filename = param.requires( "export" ) + "_" + edf.id + "_" + rowbase;

	      //std::cout << " rowids.size() = " << rowids.size() <<"  " << colids.size() << " " << tlock.X.dim2() <<  "x" << tlock.X.dim1() << "\n";

	      // save
	      massoc_t massoc( edf.id , rowids , eids, colids , tlock.X , filename );
#else
	      Helper::halt( "LGBM support not compiled in" );
#endif
	      
	    }
	  
	 	  
	  //
	  // Verbose output? Show whole matrix...
	  //
	  
	  if ( verbose )
	    {
	      for (int i=0; i<tlock.X.dim1(); i++)
                {
                  writer.level( t[i] , "SEC" );

		  for (int j=0; j<tlock.X.dim2(); j++)
		    {
		      writer.level( j+1 , "N" );
		      writer.value( "V" , tlock.X(i,j) );
		    }
		  writer.unlevel( "N" );
		}
              
	      writer.unlevel( "SEC" );
	    }

	  //
	  // all done
	  //
	  
	} // next signal
      
      writer.unlevel( globals::signal_strat );

      //
      // clear key output stratifiers
      // (clear now versus after last loop, as the strata may differ between keys)
      //
      

      ss = cc->stratum.begin();
      while ( ss != cc->stratum.end() )
	{
	  writer.unlevel( "s" + ss->first );
	  ++ss;
	}


      //
      // Next strataum
      //

      ++cc;
    }
  
  

}


tlock_t::tlock_t( const std::vector<double> & t , const int norm_points ) 
  : t(t) , norm_points(norm_points)
{
  verbose = false;
}

void tlock_t::add( const std::vector<double> * x , 
		   const int lower , const int higher , 
		   const bool take_log , 
		   const int angle_bins )
{

  if ( higher - lower + 1 != t.size() )
    Helper::halt( "internal error");

  Data::Vector<double> d( t.size() ) ;


  //
  // is this an angle? (build histogram)
  //

  if ( angle_bins )
    {
      // assume angle as 0..360
      double wa = 360.0 / (double)angle_bins;

      int j = 0 ;

      for ( int i = lower ; i <= higher ; i++ )
	{
	  // expect radians, -PI .. +PI
	  double deg = ( M_PI + (*x)[i] ) * 180.0 / M_PI;
	  if ( deg < 0 || deg > 360 ) Helper::halt( "value not a valid angle" );
	  int ia = deg / wa;
	  if ( ia == angle_bins ) ia = 0;
	  // store index
	  d[j++] = ia;
	}
      
      // at start, need to set up X to be [ angle bins x epochs ] 
      
      if ( X.dim1() == 0 )
	X.resize( d.size() , angle_bins  );
     
      // accumlate counts
      for (int j=0;j<d.size(); j++)
	X(j,d[j])++;
            
      
    }

  //
  // or else a regular value (take means, after normalizing, and adding smallest value if negative before division )
  //
  
  else
    {

      int j = 0;
      
      if ( take_log ) 
	for ( int i = lower ; i <= higher ; i++ )
	  d[j++] = log( (*x)[i] );
      else
	for ( int i = lower ; i <= higher ; i++ )
	  d[j++] = (*x)[i] ;
      
      if ( verbose ) 
	{
	  if ( X.dim1() == 0 )
	    {
	      X.resize( d.size() , 1 );
	      for (int j=0;j<d.size(); j++)
		X(j,0) = d[j];
	    }
	  else
	    X.add_col ( d ) ;
	}
      else
	{
	  if ( means.size() == 0 )
            {
              means = d; 
	      count = 1;
	    }
          else
            {
              for (int j=0;j<d.size(); j++)
                means[j] += d[j];
	      ++count;
	    }
	} // end of verb/non-verb mode

    } // end of regular value accumulator
  
}


Data::Vector<double> tlock_t::average( ) const
{  
  
  Data::Vector<double> means1 = means;

  //
  // get/calculate means 
  //

  if ( verbose ) 
    {
      // return row means (transpose + col means ) 
      Data::Matrix<double> Xt = Statistics::transpose( X );
      means1 = Statistics::mean( Xt );  
    }
  else
    {      
      for (int i=0; i<means1.size(); i++)
	means1[i] /= (double)count;
    }

  //
  // normalize by window edges (e.g. np=0.2 --> default 20% either side)
  //
  
  if ( norm_points )
    {
      
      const int n = means1.size() ; // should be 'points'
      
      // 1) rescale to minimumm of 0.0
      
      double minval = means1[0];
      for (int i=0; i<n; i++) 
	if ( means1[i] < minval ) minval = means1[i];
      for (int i=0; i<n; i++)
	means1[i] -= minval;
      
      // 2) normalize to get value of 1.0 for baseline, based on edges
      
      double norm = 0;
      for (int i=0; i<norm_points; i++)
	{
	  norm += means1[i];
	  norm += means1[n-(i+1)];
	}
      norm /= 2.0 * norm_points;
      for (int i=0; i<n; i++) means1[i] /= norm;
    }
     
  return means1;

}

Data::Matrix<double> tlock_t::angles() const 
{
  // [ time-points/rows X bins/cols ] 
  // normalize each time points 
  Data::Matrix<double> C = Statistics::transpose( X );
  // make each row sum to 1.0
  Data::Vector<double> S = Statistics::col_sums( C );
  
  for (int i=0;i<C.dim1();i++)
    for(int j=0;j<C.dim2();j++)
      C(i,j) /= S[j];
  
  return C;
}
