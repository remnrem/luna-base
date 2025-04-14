
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

#include "edf/insert.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/logger.h"
#include "helper/helper.h"
#include "stats/eigen_ops.h"
#include "dsp/tsync.h"
#include "dsp/xcorr.h"
#include "dsp/spline.h"

extern writer_t writer;
extern logger_t logger;

edf_inserter_t::edf_inserter_t( edf_t & edf , param_t & param )
{
  
  // base EDF, already attached
  // second EDF, with signals to be inserted

  // both EDFs must be continuous (or effectively contiuous, i.e. no gaps)
  
  // two modes
  //   1) given pairs of signals, find the lag for each based on xcorr() or euclidean distance
  //        edf= pairs= w= verbose

  //   2) given a lag value (in seconds, or tp-units) insert the second channel, but with an offset 
  //        edf= sig (from secondary EDF)   align=<tp-units>

  annotation_set_t annotations2;
  edf_t edf2( &annotations2) ;
  
  if ( ! edf2.attach( param.requires( "edf" ) , "." ) )
    Helper::halt( "problem attaching second EDF, edf=" + param.value( "edf" ) );
  
  // currently, both EDFs must be continuous - or at least no gaps
  if ( ! ( edf.header.continuous &&  edf2.header.continuous ) )
    if ( edf.is_actually_discontinuous() || edf2.is_actually_discontinuous() )
      Helper::halt( "neither EDF can be discontinous, with gaps" );

  // insert mode?
  if ( param.has( "offset" ) )
    {
      const std::string signal_label = param.value( "sig" );
      const double offset = param.requires_dbl( "offset" );
      const std::string annot_label = param.has( "annot" ) ? param.value( "annot" ) : "" ;

      // optionally, timestretch secondary signal?   drift sec per 'secs' secs
      const double stretch_denom = param.has( "secs" ) ? param.requires_dbl( "secs" ) : 1 ;
      const double stretch_shift = param.has( "drift" ) ? param.requires_dbl( "drift" ) : 0 ;
      
      // e.g. if offset shifts (constant) by -10 seconds in 8 hours
      //  secs=28800 drift=-10
      // where 28800 = 8 * 60 * 60 
      
      // if just do shift  implies per 1 second... but to avoid floating point issues, probably
      // better to give a reasonable denominator (with secs) 
      const bool timestretch = param.has( "drift" ) && stretch_denom > 0 ;
      const double fac = timestretch ? stretch_shift / stretch_denom : 0 ; 
      
      insert( edf , edf2 , signal_label , offset , timestretch ? &fac : NULL , annot_label );
      // all done
      return;
    }


  // method
  const bool euclidean = ! param.has( "xcorr" );
  
  // verbose outputs
  const bool verbose = param.has( "verbose" );


  // XCORR: range to search (seconds)
  const double tmax = param.has( "w" ) ? param.requires_dbl( "w" ) : -1 ;
  const double tcent = param.has( "c" ) ? param.requires_dbl( "c" ) : 0 ; 

  
  // get signal pairs from edf1
  std::vector<std::string> pairs = param.strvector( "pairs" );
  
  if ( pairs.size() == 0 || pairs.size() % 2 )
    Helper::halt( "expecting an even number of channels pairs=sig1,sig2,sig1,sig2,... (or run in insert mode with offset arg)" );
  
  const int np = pairs.size() / 2 ;

  // slots
  std::vector<int> slot1, slot2, srs;
  std::vector<std::string> slab1, slab2;

  if ( tmax > 0 ) 
    logger << "  estimating signal lag, up to a maximum of " << tmax << " seconds, centered on " << tcent << " seconds\n";
  else
    logger << "  estimating signal lag (unconstrained)\n";
  
  for (int s=0; s<np*2; s+=2)
    {
      const int s1 = edf.header.signal( pairs[s] );
      const int s2 = edf2.header.signal( pairs[s+1] );

      if ( s1 == -1 ) Helper::halt( "could not find " + pairs[s] + " in primary EDF" );
      if ( s2 == -1 ) Helper::halt( "could not find " + pairs[s+1] + " in secondary EDF" );
      
      slot1.push_back( s1 );
      slot2.push_back( s2 );

      slab1.push_back( pairs[s] );
      slab2.push_back( pairs[s+1] );
            
      // SRs must match
      int sr1 = edf.header.sampling_freq( s1 );
      int sr2 = edf2.header.sampling_freq( s2 );

      if ( sr1 != sr2 )
	Helper::halt( "sample rates must match for " + pairs[s] + " and " + pairs[s+1] );

      srs.push_back( sr1 );
      
    }

  //
  // Enfore similar SR required for all signals... easier...
  //

  int sr = srs[0];
  for (int p=0; p<np; p++)
    if ( srs[p] != sr )
      Helper::halt( "sample must match for all signals" );
  
  
  // ------------------------------------------------------------
  // 
  // XCORR method      
  //
  // ------------------------------------------------------------
  

  if ( ! euclidean )
    {

      //
      // Get data 
      //
      
      for (int p=0; p<np; p++)
	{
	  
	  slice_t slice1( edf  , slot1[p] , edf.timeline.wholetrace() );
	  slice_t slice2( edf2 , slot2[p] , edf2.timeline.wholetrace() );
	  
	  const std::vector<double> * dx = slice1.pdata();
	  const std::vector<double> * dy = slice2.pdata();
	  
	  const int tmax_sp = tmax * srs[p];
	  const int tcent_sp = tcent * srs[p];
	  
	  xcorr_t xcorr( *dx , *dy , tmax_sp , tcent_sp );
	  
	  double lag_sec = xcorr.lags[xcorr.mx] / (double)srs[p];
	  
	  logger << "  cross-correlation for " << slab1[p] << " x " << slab2[p] << " estimated lag " << lag_sec << "\n";
	  
	  writer.level( slab1[p] + ".." + slab2[p] , "CHS" );
	  
	  writer.value( "SR" , srs[p] );
	  writer.value( "L1" , (int)dx->size() );
	  writer.value( "L2" , (int)dy->size() );
	  writer.value( "LAG_SP" , xcorr.lags[xcorr.mx] );
	  writer.value( "LAG_SEC" , lag_sec );
	  writer.value( "MX" , xcorr.C[xcorr.mx] );
	  
	  if ( verbose && tmax > 0 )
	    {
	      bool output = false;
	      
	      for (int i=0; i< xcorr.lags.size(); i++)
		{
		  const double t = xcorr.lags[i] / (double)srs[p];
		  if ( t >= tcent-tmax && t <= tcent+tmax )
		    {
		      output = true;
		      writer.level( xcorr.lags[i] , "SP" );
		      writer.value( "T" , t );
		      writer.value( "XC" , xcorr.C[i] );
		    }
		}
	      
	      if ( output ) 
		writer.unlevel( "SP" );	  
	    }
	  
	  // next pair of channels
	}
      writer.unlevel( "CHS" );
    }

  
  // ------------------------------------------------------------
  //
  // Distance-based slide
  //
  // ------------------------------------------------------------

  if ( euclidean )
    {

      
      const double ystart_sec = param.requires_dbl( "start" );
      const double ylen_sec   = param.requires_dbl( "len" );

      const double yinc_sec   = param.has( "inc" ) ? param.requires_dbl( "inc" ) : 600; // 10 min jumps
      const int    ysteps     = param.has( "steps" ) ? param.requires_int( "steps" ) : 1;

      const bool constrained_offset = param.has( "offset-range" );

      std::vector<double> offsets; 
      if ( constrained_offset )
	{
	  offsets = param.dblvector( "offset-range" );
	  if ( offsets.size() != 2 || offsets[1] <= offsets[0] ) Helper::halt( "expecting offset-range=min,max" );
	}
      
      // nb. only really makses sense when steps == 1 
      if ( verbose && ysteps != 1 )
	Helper::halt( "do not advise 'verbose' with multiple steps" );
      
      int ystart = ystart_sec * sr;
      const int ylen   = ylen_sec * sr;
      const int yinc   = yinc_sec * sr;
      
      // store all data
      std::vector<Eigen::VectorXd> dX(np), dY(np);
      
      for (int p=0; p<np; p++)
	{
	  
          slice_t slice1( edf  , slot1[p] , edf.timeline.wholetrace() );
          slice_t slice2( edf2 , slot2[p] , edf2.timeline.wholetrace() );
	  
	  const std::vector<double> * dx = slice1.pdata();
          const std::vector<double> * dy = slice2.pdata();
	  
	  const int nx = dx->size();
	  const int ny = dy->size();
	  
	  Eigen::VectorXd xx = Eigen::VectorXd::Zero( nx );
	  Eigen::VectorXd yy = Eigen::VectorXd::Zero( ny );
	  
	  for (int i=0; i<nx; i++) xx[i] = (*dx)[i];
	  for (int i=0; i<ny; i++) yy[i] = (*dy)[i];
	    
	  dX[p] = xx;
	  dY[p] = yy;

	}
      
      
      // now we have all data in dX and dY
           
      // total number of samples in each dataset
      const int nx = dX[0].size();
      const int ny = dY[0].size();
      
      logger << "  based on " << ysteps << " " << ylen_sec << "s segment(s), starting "
	     << ystart_sec << "s past 2ndary EDF start, advancing " << yinc_sec << "s each step\n";
      
      // consider each possible starting point

      int steps = 0;

      while ( 1 )
	{
	  // std::cout << " \nBEGIN\n";
	  // std::cout << "  ystart = " << ystart << "\n";
	  
	  if ( ystart + ylen >= ny )
	    {
	      logger << " done, reached end of secondary signal\n";
	      break;
	    }

	  // enough steps?
	  ++steps;
	  if ( steps > ysteps ) break;

	  writer.level( ystart / (double)sr , "WIN" );
	  
	  // Get segment for yy
	  std::vector<Eigen::VectorXd> sY(np);
	  for (int p=0;p<np;p++)
	    sY[p] = dY[p].segment( ystart , ylen );
	  
	  // 0 1 2 3 4 5 6 7 8 9
	  // 1 2 3
	  //               1 2 3
	  
	  // e.g. 8 possible alignments, going from 0 to 7 
	  const int na = nx - ylen + 1 ;
	  
	  // by default
	  // align 'y' segment from 0 up to nx - ylen 
	  // given constrains, only consider alignments compatible with those

	  int mina = 0;
	  int maxa = na;

	  if ( constrained_offset )
	    {
	      // in sample points
	      int minoff = offsets[0] * sr;
	      int maxoff = offsets[1] * sr;
	      
	      // offset = ystart - a 
	      // offset if what we substract from Y
	      // i.e. -ve offset means  Y - offset --> starts at +ve point in X
	      
	      // a = ystart - offset
	      mina = ystart - minoff;
	      maxa = ystart - maxoff;

	      if ( mina > maxa )
		{
		  int aa = mina;
		  mina = maxa;
		  maxa = aa;
		}

	      if ( mina < 0 ) mina = 0;
	      if ( maxa < 0 ) maxa = 0;
	      if ( mina >= na ) mina = na;
	      if ( maxa >= na ) maxa = na;

	      std::cout << " considering offsets "
			<< ( ystart - mina ) / (double)sr << " to " << ( ystart - maxa ) / (double)sr << "\n";
	      std::cout << " based on a range " << mina << " - " << maxa << " ..... " << mina/(double)sr << " - "
			<< maxa / (double)sr << "\n";
	      
	    }
	  
	  // score each alignment
	  std::vector<double> st( na , 0 );
	  
	  for (int p=0; p<np; p++)
	    for (int a=mina; a <maxa  ; a++)
	      st[a] += ( dX[p].segment(a,ylen) - sY[p] ).norm() ;
	  
	  // get min.
	  int minidx = 0;
	  double minst = 0;
	  
	  for (int a=mina; a < maxa ; a++)
	    {
	      if ( a == mina || st[a] < minst )
		{
		  minst = st[a];
		  minidx = a;
		}
	    }

	  // console outputs
	  logger << "  for segment starting " << ystart / (double)sr << "s, optimal offset = " << ystart - minidx
		 << " (" << ( ystart - minidx ) / (double)sr << "s )\n";

	  // output offset	  
	  if ( verbose ) 
	    for (int a=0; a < na ; a++)
	      std::cout << st[a] << "\n";
	  
	  writer.value( "SP" , ystart - minidx );
	  writer.value( "SEC" , ( ystart - minidx ) / (double)sr );
	  
	  // advance to the next block
	  ystart += yinc;
	  
	}
      writer.unlevel( "WIN" );
        
    }
  
}



void edf_inserter_t::insert( edf_t & edf , edf_t & edf2 , const std::string & siglabel ,
			     const double offset , const double * fac, 
			     const std::string annot_label )
{
  
  // insert as much of the signals from edf2 into edf
  // assume both are (effectively) continuous
  // by default, align at 0, i.e. start of edf2 signal equals start of edf, and add as much as we can  

  // if an offset is specified, then we insert after adding an offset
  //   -ve offset implies EDF2 start is AFTER EDF start (i.e. it needs to be shifted backwards) so that starts align
  //   +ve offset implies the reverse:  EDF2 start is BEFORE EDF start, needs to go forward 

  //            S
  // EDF        |-----------------------|

  // EDF2       |-----------------------|      offset = 0
  // EDF2       |-----------------|000000      offset = 0 , pad with zeros and add annotation to indicate missing signal

  // EDF2       000|--------------------|XXX|  offset = -ve (i.e. EDF2 start is after EDF start): pad w/ zeros; truncate at end X as needed
  // EDF2   |XXX|--------------------|000      offset = +ve (i.e. EDF2 start is before EDF start), need to shift forwards 

  
  // optionally, we can add annotations to indicate where the signal is missing from edf2 

  // if fac (secs/shift) is non-null, then apply this timestretch factor to the inserted channel
  // i.e. this is to adjust for linear difference in clock rates.

  const bool timestretch = fac != NULL; 
  
  const bool no_annotations = true;

  signal_list_t signals = edf2.header.signal_list( siglabel , no_annotations );

  const int ns = signals.size();
  
  logger << "  inserting " << ns << " signals from " << edf2.filename << ", ";
  logger << "using an offset of " << offset << " seconds\n";

  if ( timestretch )
    {
      if ( *fac > 0 ) 
	logger << "  shrinking secondary signals by a rate of " <<  (*fac) << " sec per second\n";
      else
	logger << "  stretching secondary signals by a rate of " << -1 * (*fac) << " sec per second\n";
    }
  


  for (int s=0; s<ns; s++)
    {
      // get SR of second signal (EDF2)
      const int Fs = edf2.header.sampling_freq( signals(s) );
      
      // new signal, with length to match EDF  
      const int np = edf.header.nr * edf.header.record_duration * Fs;

      // putative new signal name
      std::string sig = signals.label(s);
      
      // check new signal name will be unique
      if ( edf.header.has_signal( sig ) )
	{
	  int j = 1;
	  while ( 1 )
	    {
	      std::string sig2 = sig + "." + Helper::int2str( j );
	      if ( ! edf.header.has_signal( sig2 ) )
		{
		  sig = sig2;
		  break;
		}
	      ++j;
	    }
	}

      // make a new vector, set to zero-pad
      std::vector<double> d1( np , 0 );
      
      // pull the EDF2 signal
      slice_t slice( edf2 , signals(s) , edf2.timeline.wholetrace() );      
      std::vector<double> d2 = *slice.pdata();
      
      // calculate best offset in sample points
      const int offset_sp = offset * Fs;

      // time-stretch?
      if ( timestretch )
	{

	  // if a time-stretch factor is defined, then
	  // use spline interpolation to stretch or shrink
	  // the secondary signal as needed (based on
	  // the ratio of secs and shift as specified on
	  // the INSERT command line
	  
	  const int n_orig = d2.size();
	  const int n_scaled = n_orig - n_orig * (*fac);
	  
	  if ( n_scaled <= 0 ) Helper::halt( "rescaled signal not defined" );

	  std::vector<double> t( n_orig );
	  for (int i=0; i<n_orig; i++) t[i] = i;
	  
	  tk::spline spline;
	  spline.set_points( t, d2 );
	  
	  d2.clear();
	  d2.resize( n_scaled , 0 );
	  
	  for (int i=0; i<n_scaled; i++)
	    d2[i] = spline( n_orig * ( i / (double)n_scaled ) );	  
	  
	}
            
      // console messages
      logger << "  inserting " << sig << " ( SR = " << Fs << " Hz, offset = " << offset_sp << " samples ) into primary EDF\n";

      // -ve offset : pad w/ new N zeros , skip last N 
      //   ||||||||     d1
      //     ||||||||   d2   <---- need to shift back, aka start early 
      //   00IIIIII
      
      // +ve offset : pad w/ new N zeros , skip first N 
      //   ||||||||
      // ||||||||          -----> need to shift forward, but means start late in 2ndary
      //    IIIIII00
      
      // pointers (sample points to both files)
      int p1 = 0 , p2 = offset_sp ;
      
      // signal lengths
      const int n1 = np;
      const int n2 = d2.size();

      // advance sp in 2ndary signal 
      while ( 1 )
	{

	  // all done?
	  if ( p1 == n1 ) break;
	  
	  // zero-pad if out of input
	  if ( p2 >= n2 ) 
	    d1[ p1 ] = 0 ;
	  else if ( p2 >= 0 ) // else add or 0-pad (i.e. not yet at start of EDF1)
	    d1[ p1 ] = d2[ p2 ];
	  
	  // advance
	  ++p1;
	  ++p2;	  	  
	}

      // add the new signal
      edf.add_signal( sig , Fs , d1 );
      
      // add any annotations 
      // -- todo --
      
    }
  
  
  
}



  
