
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
#include "spectral/mt_spectrogram.h"


extern writer_t writer;
extern logger_t logger;

// tlock_t can operate either by cache points, or now by epochs too (i.e. as these
// can be defined via annotations and therefore have fixed width

void dsptools::tlock( edf_t & edf , param_t & param )
{
  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , true );
  
  if ( signals.size() < 1 ) return;
  
  const int ns = signals.size();

  // check sample rates: must be similar across all 

  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  for (int s=1;s<Fs.size();s++) 
    if ( Fs[s] != Fs[0] ) 
      Helper::halt( "sample rates must be similar across signals for TLOCK" );
  
  
  //
  // by epoch, or by cache?
  //
  
  const bool by_epoch = param.has( "epoch" );

  const bool by_cache = param.has( "cache" );

  if ( by_epoch == by_cache ) Helper::halt( "must specify either epoch or cache" );

  //
  // -------------------- epochs options --------------------
  //

  // can set '0' to any point
  double emid = param.has( "mid" ) ? param.requires_dbl( "mid" ) : 0 ; 
  
  //
  // -------------------- cache options --------------------
  //
  
  //
  // get time-points: expecting a cache in sample-point units
  //   also implies signals must have the same SR
  //
  
  std::string cache_name = by_epoch ? "" : param.requires( "cache" );

  if ( ! by_epoch ) 
    if ( ! edf.timeline.cache.has_int( cache_name ) )
      Helper::halt( "cache not found for this individual: " + cache_name );
  
  cache_t<int> * cache = by_epoch ? NULL : edf.timeline.cache.find_int( cache_name );
  
  // only self-channels? by default, gives sync'ed values for all
  // channels in sig even if the cache has a CH strata that is
  // different (i.e. seeded on sCH = "C3" ) to only give output for
  // same channels, add 'same-channel=T'
  
  const bool same_channel = param.yesno( "same-channel" );
  
  // also allow to match e.g. 'same-channel'  if seed = C3
  //   same-channel=T  channel-postfix=_SIGMA
  //  will match sig = C3 /and/ C3_SIGMA for this seed

  const std::string channel_postfix = param.has( "channel-postfix" ) ? param.value( "channel-postfix" ) : ""; 
  if ( channel_postfix != "" && ! same_channel )
    Helper::halt( "cannot specify channel-postfix without same-channel=T" );


  //
  // -------------------- Input options --------------------
  //

  // take logs before means?

  bool take_log = param.has( "tolog" ) ;
  
  // treat as angles (circular values) 
  //   --> output a binned histogram of counts rather than circular means, etc
  
  int angle_bins = param.has( "phase" ) ? param.requires_int( "phase" ) : 0 ;
  
  if ( take_log && angle_bins ) Helper::halt( "cannot specify both tolog and phase" );
  
  
  //
  // -------------------- Output options --------------------
  // 
  
  bool verbose = param.has( "verbose" );
  
  bool to_massoc = param.has( "export" );
  
  bool massoc_iids = param.yesno( "unique-ids" );

  
  //
  // -------------------- Analysis/window options -----------
  //
  
  const double half_window = by_epoch ? 0 : param.requires_dbl( "w" );
  if ( (!by_epoch) && half_window <= 0 ) Helper::halt( "w must be a positive number" );

  
  
  //
  // Normalisation: e.g. 0.2 means 20% + 20% of window, i.e. skip middle 60%
  //                defauly = 0 (or negative) -> no normalisation
  
  const double norm_pct = param.has( "np" ) ? param.requires_dbl( "np" ) : -1 ;

  if ( norm_pct > 0.5 )
    Helper::halt( "expecting np between 0 and 0.5" );
  
  const bool zero_trace = param.has( "zero" ) ? param.yesno( "zero" ) : false ; 
  
  const double outlier_th = param.has( "th" ) ? param.requires_dbl( "th" ) : -1 ; 
  
  const double outlier_winsor = param.has( "win" ) ? param.requires_dbl( "win" ) : -1; 
  
  //
  // --------------- Spectrogram options --------------------------
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
  // Primary loop, iterate over sig-specified channel (CH outputs) 
  //

  
  for (int s=0; s<ns; s++)
    {

      //
      // only data channels
      //

      if ( edf.header.is_annotation_channel( signals(s) ) )
        continue;

      writer.level( signals.label(s) , globals::signal_strat );
      
      //
      // initiate tlock_t object
      //
    
      tlock_t tlock( edf, Fs[s] );
      
      //
      // set options
      //
      
      tlock.take_log = take_log;
      tlock.angle_bins = angle_bins;
      tlock.norm_pct = norm_pct;
      tlock.zero_trace = zero_trace;
      tlock.outlier_th = outlier_th;
      tlock.outlier_winsor = outlier_winsor;
      tlock.emid = emid; // only epoch-mode
      
      //
      // build run tlock(s)
      //   - epoch builder will just construct a single (epoch-based) X
      //   - caches containing multiple strata will iterate (and give multiple, stratified outputs w/in this function below)
      //
      
      
      if ( by_epoch )
	tlock.epoch_builder( signals(s) );
      else
	tlock.cache_builder( cache , half_window ,
			     signals(s) , signals.label(s) , same_channel, channel_postfix );
      
      // next signal
    }
  writer.unlevel( globals::signal_strat );
  
  // all done
}


// // report number of interval considered/accepted
//       writer.value( "N" , cnt_valid_intervals );
//       writer.value( "N_ALL" , (int)cx.size() );
	  	  
//       logger << "  included " << cnt_valid_intervals
// 	     << " of " << (int)cx.size()
// 	     << " intervals for strata " << 
// 	     << sstr.str() << " for channel " << signals.label(s) << "\n";
      
// 	  //
// 	  // Report either as phase angles (assuming radians), otherwise take the mean
// 	  //

// 	  if ( angle_bins ) 
// 	    {
// 	      Data::Matrix<double> angbin = tlock.angles();
	      
// 	      if ( angbin.dim1() != angle_bins || angbin.dim2() != t.size() )
// 		Helper::halt( "internal error in tlock_t()" );

// 	      for (int i=0; i<angle_bins; i++)
// 		{
// 		  writer.level( i+1 , "PHASE" );
// 		  for (int j=0; j<t.size(); j++)
// 		    {
// 		      writer.level( t[j] , "SEC" );		      
// 		      writer.value( "M" , angbin(i,j) );
// 		    }
// 		}
// 	      writer.unlevel( "PHASE" );
// 	      writer.unlevel( "SEC" );
	      
// 	    }

// 	  else if ( spectrogram )
// 	    {

// 	      logger << "  calculating mean MT spectrogram...\n";

// 	      Data::Vector<double> means = tlock.average( );
	      
// 	      mt_spectrogram_t mtm( tlock.X ,
// 				    Fs[0] ,
// 				    mtm_nw , mtm_t,
// 				    mtm_seg , mtm_step,
// 				    mtm_fmin, mtm_fmax,
// 				    mtm_dB,
// 				    mtm_center );

// 	      const int nf = mtm.frq.size();
// 	      const int nt = mtm.t.size();

// 	      // mean-centered, normalized versions (by row/freq and by col/time)
	      
// 	      // double min = 999;
// 	      // for (int i=0; i<nf; i++)
// 	      // 	for (int j=0; j<nt; j++)
// 	      // 	  if ( mtm.Z(i,j) < min ) min = mtm.Z(i,j);

// 	      //Data::Matrix<double> Z0 = mtm.Z;
	      
// 	      // for (int i=0; i<nf; i++)
// 	      // 	for (int j=0; j<nt; j++)
// 	      // 	  Z0(i,j) += min;

// 	      std::vector<double> row_min( nf );
// 	      std::vector<double> row_max( nf );
// 	      std::vector<double> row_mean( nf , 0 );

// 	      for (int i=0; i<nf; i++)
// 		{
// 		  row_min[i] = mtm.Z(i,0);
// 		  row_max[i] = mtm.Z(i,0);
// 		  row_mean[i] = mtm.Z(i,0);
		  
// 		  for (int j=1; j<nt; j++)
// 		    {
// 		      if ( mtm.Z(i,j) < row_min[i] ) row_min[i] = mtm.Z(i,j);
// 		      if ( mtm.Z(i,j) > row_max[i] ) row_max[i] = mtm.Z(i,j);		      
// 		      row_mean[i] += mtm.Z(i,j);
// 		    }
// 		  row_mean[i] /= (double)nt;
// 		}
	      
// 	      // output
// 	      for (int i=0; i<nf; i++)
// 		{
// 		  writer.level( mtm.frq[i] , globals::freq_strat );
// 		  for (int j=0; j<nt; j++)
// 		    {
// 		      writer.level( mtm.t[j] , "SEC" );
// 		      // mean
// 		      writer.value( "PSD" , mtm.Z(i,j) );

// 		      // median
// 		      writer.value( "PSD_MD" , mtm.Z_median(i,j) );		      
		      
// 		      // normalized by F (row) min/max
// 		      writer.value( "PSD_F" , ( mtm.Z(i,j) - row_min[i] ) / ( row_max[i] - row_min[i] ) );

// 		      // normalized by row mean
// 		      writer.value( "PSD_M" , mtm.Z(i,j) - row_mean[i] ) ;
		      
// 		      // variance in PSD
// 		      writer.value( "VAR" , mtm.ZZ(i,j) );
// 		    }
// 		  writer.unlevel( "SEC" );
// 		}
// 	      writer.unlevel( globals::freq_strat );
// 	    }
	 
// 	  //
// 	  // report summaries: for regular values, get the average
// 	  // and normalize (by edges)
// 	  //
	  
// 	  else
// 	    {
// 	      Data::Vector<double> means = tlock.average( );
	  
// 	      if ( means.size() != t.size() ) 
// 		{
// 		  logger << "  means.size() = " << means.size() << " " << t.size() << "\n";
// 		  Helper::halt( "internal error in tlock_t()" );
// 		}	      
	      
// 	      for (int i=0; i<means.size(); i++) 
// 		{
// 		  writer.level( t[i] , "SEC" );
// 		  writer.value( "M" , means[i] );
// 		}

// 	      writer.unlevel( "SEC" );
// 	    }

// 	  //
// 	  // Dump for MASSOC?
// 	  //

// 	  if ( to_massoc )
// 	    {
// #ifdef HAS_LGBM
	      
// 	      // filename = indiv-id + strata
// 	      // ID =
// 	      // variables = T

// 	      const int nrow = tlock.X.dim2(); // note - obs in cols here
// 	      const int ncol = tlock.X.dim1();
	      
// 	      // row-IDs : strata + N
// 	      //  EDFID_C3_11_1
// 	      //  EDFID_C3_11_2 ...
// 	      //  etc

// 	      // optionally, make indiv-IDs unique, i.e. so we can select individual
// 	      // events as trainers.  This effectively breaks the "indiv" hierarchy.
// 	      // (that is, we will not be able to manipulate all events belonging to
// 	      // one person)
	      
// 	      std::vector<std::string> iids( massoc_iids ? nrow : 0 );
// 	      std::vector<std::string> rowids( nrow );
// 	      std::vector<std::string> eids( nrow );
	      
// 	      std::string rowbase = signals.label(s);
// 	      std::map<std::string,std::string>::const_iterator ss = cc->stratum.begin();
// 	      while ( ss != cc->stratum.end() )
// 		{
// 		  // nb. here 'CH' means seed channel;  SIG means the readout
// 		  rowbase += "_" + ss->second ;
// 		  ++ss;
// 		}

// 	      if ( eidx_base1.size() != nrow )
// 		Helper::halt( "internal error in TLOCK w/ eidx_base1 size" );
	      
// 	      for (int i=0; i<nrow; i++)
// 		{
		  
// 		  rowids[i] = rowbase ;
// 		  eids[i] = Helper::int2str( eidx_base1[i] );
		  
// 		  // make unique IDs?
// 		  if ( massoc_iids ) iids[i] = "_ind_" + edf.id + "_row_" + rowids[i] + "_obs_" + eids[i];
// 		}
	      
// 	      // col IDs: simply 1, 2, 3, etc
// 	      std::vector<std::string> colids( ncol );
// 	      for (int i=0; i<ncol; i++)
// 		colids[i] = Helper::int2str( i+1 );
	      
// 	      // i.e. expecting something like export=path/
// 	      //      or export=path/root
// 	      //        filename will append root_ID_ROWBASE
// 	      const std::string filename = param.requires( "export" ) + "_" + edf.id + "_" + rowbase;

// 	      //std::cout << " rowids.size() = " << rowids.size() <<"  " << colids.size() << " " << tlock.X.dim2() <<  "x" << tlock.X.dim1() << "\n";

// 	      // save

// 	      if ( massoc_iids )
// 		{
// 		  logger << "  writing with unique indiv-level IDs\n";
// 		  massoc_t massoc( iids , rowids , eids, colids , tlock.X , filename );
// 		}
// 	      else
// 		{
// 		  logger << "  preserving indiv-level IDs\n";
// 		  massoc_t massoc( edf.id , rowids , eids, colids , tlock.X , filename );
// 		}
// #else
// 	      Helper::halt( "LGBM support not compiled in" );
// #endif
	      
// 	    }
	  
	 	  
// 	  //
// 	  // Verbose output? Show whole matrix...
// 	  //
	  
// 	  if ( verbose )
// 	    {
// 	      for (int i=0; i<tlock.X.dim1(); i++)
//                 {
//                   writer.level( t[i] , "SEC" );

// 		  for (int j=0; j<tlock.X.dim2(); j++)
// 		    {
// 		      writer.level( j+1 , "N" );
// 		      writer.value( "V" , tlock.X(i,j) );
// 		    }
// 		  writer.unlevel( "N" );
// 		}
              
// 	      writer.unlevel( "SEC" );
// 	    }

// 	  //
// 	  // all done
// 	  //
	  
// 	} // next signal
      
//       writer.unlevel( globals::signal_strat );

//       //
//       // clear key output stratifiers
//       // (clear now versus after last loop, as the strata may differ between keys)
//       //
      

//       ss = cc->stratum.begin();
//       while ( ss != cc->stratum.end() )
// 	{
// 	  writer.unlevel( "s" + ss->first );
// 	  ++ss;
// 	}


//       //
//       // Next strataum
//       //

//       ++cc;
//     }
  
  

// }



void tlock_t::clearX()
{
  // clear data, but keep time-track and other options as is
  //  i.e. when resetting for a new signal/strata
  
  X.clear();
  ni = 0;
  
}


// void tlock_t::norm_within_intervals( const int np1 )
// {
//   if ( np1 <= 0 ) return;

//   // note - actually wants to norm cols, as these are 

//   const int nrow = X.dim1();
//   const int ncol = X.dim2();
  
//   // rescale each row to get min = 0
//   // and edges mean == 1 
  
//   for (int j=0; j<ncol; j++)
//     {
      
//       // 1) rescale to minimumm of 0.0
      
//       double minval = X(0,j);
//       for (int i=0; i<nrow; i++) 
// 	if ( X(i,j) < minval ) 
// 	  minval = X(i,j);
//       for (int i=0; i<nrow; i++)
// 	X(i,j) -= minval;
      
//       // 2) normalize to get value of 1.0 for baseline, based on edges
      
//       double norm = 0;
//       for (int i=0; i<np1; i++)
// 	{
// 	  norm += X(i,j);
// 	  norm += X(nrow-(i+1),j);
// 	}
//       norm /= 2.0 * np1;
//       for (int i=0; i<nrow; i++) 
// 	X(i,j) /= norm;
//     }

// }


Data::Vector<double> tlock_t::average( const double th , const double winsor ) const 
{  

  //
  // get means
  //
  
  // return row means (transpose + col means ) 
  Data::Matrix<double> Xt = Statistics::transpose( X );
  
  if ( th > 0 || winsor > 0 )
    Xt = remove_outliers( Xt , th , winsor );

  Data::Vector<double> means1 = Statistics::mean( Xt );
  

  //
  // normalize mean values by window edges (e.g. np=0.1 --> default 10% either side)?
  //
  
  if ( norm_pct > 0 || zero_trace )
    edge_normalization( &means1 , norm_pct * np );
  
  
  return means1;
}


Data::Matrix<double> tlock_t::remove_outliers( const Data::Matrix<double> & Y , const double th , const double winsor ) const 
{
  
  // std::vector<bool> exclude( ni , false );
  

  // for (int i=0; i<np; i++)
  //   {
  //     const std::vector<double> & 

  // 	  int outliers( const std::vector<double> * x , double th , 
  // 		std::vector<bool> * inc , const std::vector<bool> * prior = NULL );  


  //   }
  // // remove 
  
  return Y;
}

Data::Vector<double> tlock_t::median( const double th , const double winsor ) const
{  

  //
  // get means
  //
  
  // return row means (transpose + col means ) 
  Data::Matrix<double> Xt = Statistics::transpose( X );

  if ( th > 0 || winsor > 0 ) 
    Xt = remove_outliers( Xt , th , winsor );
  
  Data::Vector<double> med1( np );
  for (int i=0; i<np; i++) 
    med1[i] = MiscMath::median( *Xt.col(i).data_pointer() );
  

  //
  // normalize mean values by window edges (e.g. np=0.1 --> default 10% either side)?
  //
  
  if ( norm_pct > 0 || zero_trace )
    edge_normalization( &med1 , norm_pct * np );
  
  
  return med1;
}



void tlock_t::edge_normalization( Data::Vector<double> * m , const int p ) const
{

  // invalid sizes
  const int n = m->size();
  if ( n == 0 || n < 2 * p ) return;

  // rescale to minimumm of 0.0  
  if ( zero_trace ) 
    {
      std::cout <<" zero-ing\n";
      double minval = (*m)[0];
      for (int i=0; i<n; i++) 
	if ( (*m)[i] < minval ) minval = (*m)[i];
      for (int i=0; i<n; i++)
	(*m)[i] -= minval;
    }

  // normalize to get value of 1.0 for baseline, based on edges
  if ( p > 0 ) 
    { 
      double norm = 0;
      for (int i=0; i<p; i++)
	{
	  norm += (*m)[i];
	  norm += (*m)[n-(i+1)];
	}
      norm /= 2.0 * p;
      for (int i=0; i<n; i++)
	(*m)[i] /= norm;
    }

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



//
// -------------------- New code --------------------
//

tlock_t::tlock_t( edf_t & data , const int Fs )
  : edf(data) 
{
  
  sr = Fs;
  
  // number of points undefined
  np = 0;
  ni = 0;

  // options
  take_log = false;
  angle_bins = 0;

  verbose = false;
}


void tlock_t::epoch_builder( const int slot )
{

  // initially, we do not know epoch size
  // **assume** it is constant - will flag an error if not
  np = 0;
  
  // get data and TP information   
  slice_t slice( edf , slot , edf.timeline.wholetrace() );  
  const std::vector<double> * d = slice.pdata();  
  const std::vector<uint64_t> * tp = slice.ptimepoints();
  
  // Start building the X matrix
  clearX();
      
  // iterate over epochs
  edf.timeline.first_epoch();
  while ( 1 )
    {

      int epoch = edf.timeline.next_epoch();      
          
      if ( epoch == -1 ) break;
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      slice_t slice( edf , slot , interval );

      // we know this will be continuous, by definition of an 'epoch'
      std::vector<double> * d = slice.nonconst_pdata();

      // set window 
      const int np1 = d->size();

      //      std::cout << " xx np1 = " << np1 << "\t" << interval.duration() << "\n";
      
      if ( np == 0 )
	{
	  np = np1;
	  // set_window takes 
	  set_window_epoch( np1 );
	}
      else
	{
	  //std::cout << " np = " << np << " np1 " << np1 << "\n";
	  if ( abs( np - np1 ) > 1 )
	    Helper::halt( "cannot have variable-sized epochs in TLOCK" );
	}
      
      // otherwise, add this interval to the tlock 	      
      add( d , 0 , np-1 );
      ++ni;
      
      // next epoch
    }

  // run all outputs
  outputs();
}


void tlock_t::cache_builder( cache_t<int> * cache ,
			    const double half_window ,
			    const int slot , 
			    const std::string label ,
			    const bool same_channel ,
			    const std::string channel_postfix )
{

  // build np x nw matrix  [ samples x windows ] 
  //  np - number of samples in each interval/window
  //  nw - number of windows
    
  // this sets both np = 1 + 2 * half_points

  np = set_window( half_window * sr  );

  int half_points = ( np - 1 ) / 2 ;
  
  // Iterate over keys associated with the 'points' internal name cache
  
  std::set<ckey_t> ckeys = cache->keys( "points" );
  
  std::set<ckey_t>::const_iterator cc = ckeys.begin();
  
  while ( cc != ckeys.end() )
    {
      
      std::vector<int> cx = cache->fetch( *cc );
      
      //
      // any intervals? if not, skip
      // 
      
      if ( cx.size() == 0 ) 
	{
	  ++cc;
	  continue;
	}
      
      //
      // do we have a channel specification, and must this match? otherwise, skip this stata
      //
      
      std::string seed_channel = "";
      if ( same_channel && cc->stratum.find( globals::signal_strat ) != cc->stratum.end() )
	{
	  seed_channel = cc->stratum.find(  globals::signal_strat )->second;
	}

      
      // skip this strata for this channel?
      
      if ( same_channel && seed_channel != "" )
	{
	  if ( seed_channel != label &&
	       seed_channel + channel_postfix != label )
	    continue;
	}
      
      // add output stratifiers based on this key
      
      std::stringstream sstr;
      std::map<std::string,std::string>::const_iterator ss = cc->stratum.begin();
      while ( ss != cc->stratum.end() )
	{
	  writer.level( ss->second , "s" + ss->first );
	  sstr << ss->first << "=" << ss->second << ";" ;
	  ++ss;
	}
      
      // get data and TP information 
      
      slice_t slice( edf , slot , edf.timeline.wholetrace() );
      
      const std::vector<double> * d = slice.pdata();
      
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      
      
      //
      // Start building the X matrix
      //
      
      clearX();
      
      //
      // iterate over time-points (cache seeds)
      //
      
      for ( int i=0; i<cx.size(); i++)
	{
	  
	  int lower = cx[i] - half_points;
	  int upper = cx[i] + half_points;
	  
	  // interval out-of-range
	  if ( lower < 0 || upper > d->size() ) 
	    continue;
	  
	  // discontinuity?
	  if ( edf.timeline.discontinuity( *tp , sr , lower , upper ) )
	    continue;
	  
	  // otherwise, add this interval to the tlock 	      
	  add( d , lower , upper );
	  
	  ++ni;
	  
	  // next interval
	}

      //
      // run all outputs
      //
      
      outputs();

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

      // next cache strata
      ++cc;
      
    }

}
  


void tlock_t::add( const std::vector<double> * x , 
		   const int lower , const int higher )
{

  // because of small floating issues in defining windows based on annotation times (e.g. from spindles)
  // we might have +/-1 extra sample points, i.e. depending on how the annotation is aligned over sample points
  // allow this, just shift/drop one 

  const int implied = higher - lower + 1 ;
  const int expected = t.size();

  if ( abs( implied - expected ) > 1 )
    Helper::halt( "problem aligning equally-sized windows given epoch lengths and sample rates - internal error, likely floating point issues" );
  
  // if ( higher - lower + 1 != t.size() )
  //   Helper::halt( "internal error");
  
  Data::Vector<double> d( expected );
  const int nmin = expected < implied ? expected : implied ;
    
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
  // add a regular value - potentially skipping for not looking for an extra sample point
  //
  
  else
    {
      
      int j = 0;
      
      if ( take_log ) 
	for ( int i = lower ; i <= higher ; i++ )
	  {
	    d[j++] = log( (*x)[i] );
	    if ( j == nmin ) break;
	  }
      else
	for ( int i = lower ; i <= higher ; i++ )
	  {
	    d[j++] = (*x)[i] ;
	    if ( j == nmin ) break;
	  }
      
      if ( X.dim1() == 0 )
	{
	  X.resize( expected , 1 );
	  for (int j=0;j<expected; j++)
	    X(j,0) = d[j];
	}
      else
	X.add_col( d ) ;
      
     

    } // end of regular value accumulator
  
}


void tlock_t::set_window_epoch( int np1 )
{
  // np1 is total number of sample points (from epoch)
  np = np1;
  //  std::cout << " set_window_epoch() setting np = " << np << "\n";
  const double totsec = np1 / sr;
  const double half_window = totsec / 2.0 ; 
  //  logger << "  window has " << np1 << " samples, " << totsec << " seconds\n";
    
  t.clear();
  // nb. fudge for floating point issues (allowing tenth of inc for stop 'w')
  const double inc = 1.0/sr;
  
  
  for (int i=0;i<np;i++)
    {
      t.push_back( i * inc - emid );
      //      std::cout << " w i " << i << " " << i * inc - emid << "\n";
    }
  // for ( double w = -half_window ; w <= half_window + inc/10.0   ; w += inc ) 
  //   t.push_back(w - emid );  // any t-scale adjustment goes here, if emid != 0 
  
  if ( t.size() != np ) 
    Helper::halt( "internal error constructing window:" + Helper::int2str( np ) + " vs " + Helper::int2str( (int) t.size() )  );
  
}

int tlock_t::set_window( int half_points )
{
  // need option for angle_bins
  
  // ensure an nice multiple of sample rate
  int half_window = half_points / sr ;
  half_points = half_window * sr ;
  int points = 1 + 2 * half_points;
  
  t.clear();
  // nb. fudge for floating point issues (allowing tenth of inc for stop 'w')
  const double inc = 1.0/sr;
  for ( double w = -half_window ; w <= half_window + inc/10.0   ; w += inc ) 
    t.push_back(w);
  
  if ( t.size() != points ) 
    Helper::halt( "internal error constructing window:" + Helper::int2str( points ) + " vs " + Helper::int2str( (int) t.size() )  );
  
  return points;
}

void tlock_t::outputs()
{

  //  logger << " outputs\n";

  
  // std::cout << "X = " << X.dim1() << " " << X.dim2() << "\n";
  // for (int i=0;i<X.dim1();i++)
  //   {
  //     for (int j=0; j<X.dim2(); j++)
  // 	std::cout << ( j ? "\t" : "" ) << X(i,j) ;
  //     std::cout << "\n";
  //   }
  // std::cout << "\n";
  
  // means
  Data::Vector<double> m = average( outlier_th , outlier_winsor );

  // medians
  Data::Vector<double> md = median( outlier_th , outlier_winsor );
  
  if ( m.size() != np )
    {
      logger << "  means.size() = " << m.size() << " np = " << np << "\n";
      Helper::halt( "internal error in tlock_t()" );
    }
  
  writer.value( "N" , ni );

  for (int i=0; i<np; i++)
    {
      writer.level( t[i] , "SEC" );
      writer.value( "M" , m[i] );
      writer.value( "MD" , md[i] );
    }  
  writer.unlevel( "SEC" );
    
}
