
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

extern writer_t writer;
extern logger_t logger;

edf_inserter_t::edf_inserter_t( edf_t & edf , param_t & param )
{
  
  // base EDF, already attached
  // second EDF, with signals to be inserted

  // both EDFs must be continuous (or effectively contiuous, i.e. no gaps)
  
  // two modes
  //   1) given pairs of signals, find the lag for each based on xcorr()
  //        edf= pairs= w= verbose

  //   2) given a lag value (in seconds, or tp-units) insert the second channel, but with an offset 
  //        edf= sig (from secondary EDF)   align=<tp-units>

  
  edf_t edf2;
  
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
      insert( edf , edf2 , signal_label , offset , annot_label );
      // all done
      return;
    }
  
  
  // verbose outputs
  const bool verbose = param.has( "verbose" );

  // range to search (seconds)
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

      // redundant now
      // adjust search to within a window? [ reset.mx ] 
      // if ( tmax > 0 )
      // 	{
      // 	  // only need to do if max is outside this window
      // 	  if ( fabs( lag_sec ) > tmax )
      // 	    {

      // 	      double mc = 0;
      // 	      int midx = 0;
      // 	      for (int i=0; i< xcorr.lags.size(); i++)
      // 		{
      // 		  const double t = xcorr.lags[i] / (double)srs[p];
      // 		  if ( fabs( t ) <= tmax )
      // 		    {
		    
      // 		      if ( fabs( xcorr.C[i] ) > mc )
      // 			{
		
      // 			  mc = fabs( xcorr.C[i] );
      // 			  midx = i;
      // 			}
      // 		    }
      // 		}

      // 	      // update
      // 	      xcorr.mx = midx;
      // 	      lag_sec = xcorr.lags[xcorr.mx] / (double)srs[p];
      // 	    }	      
      // 	}
      
      
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

 

void edf_inserter_t::insert( edf_t & edf , edf_t & edf2 , const std::string & siglabel , const double offset , const std::string annot_label )
{
  
  // insert as much of the signals from edf2 into edf
  // assume both are (effectively) continuous
  // by default, align at 0, i.e. start of edf2 signal equals start of edf, and add as much as we can  
  // if an offset is specified, then we insert
  // optionally, we can add annotations to indicate where the signal is missing from edf2 

  const bool no_annotations = true;

  signal_list_t signals = edf2.header.signal_list( siglabel , no_annotations );

  const int ns = signals.size();
  
  logger << "  inserting " << ns << " signals from " << edf2.filename << ", ";
  logger << "using an offset of " << offset << " seconds\n";

  
  // EDF        |-----------------------|

  // EDF2       |-----------------------|      offset = 0
  // EDF2       |-----------------|000000      offset = 0 , pad with zeros and add annotation to indicate missing signal
  // EDF2       000|--------------------|XXX|  offset = -ve (i.e. EDF comes before EDF2): pad w/ zeros; truncate at end X); add annot.  
  // EDF2   |XXX|--------------------|000      offset = +ve (i.e. EDF comes after EDF2) 


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

      // make a new vector
      std::vector<double> d1( np , 0 );

      // pull the EDF2 signal
      slice_t slice( edf2 , signals(s) , edf2.timeline.wholetrace() );      
      const std::vector<double> * d2 = slice.pdata();
                
      // calculate best offset in sample points
      const int offset_sp = offset * Fs;
      
      // console messages
      logger << "  adding " << sig << " ( SR = " << Fs << " Hz, offset = " << offset_sp << " samples ) to primary EDF\n";

      // +ve offset : pad w/ new N zeros , skip first N 
      //    ||||||||
      //  ||||||||

      // -ve offset : pad w/ new N zeros , skip last N 
      //  ||||||||
      //    ||||||||

      // pointers (sample points to both files)
      int p1 = 0; 
      int p2 = offset_sp; 
      
      // signal lengths
      const int n1 = np;
      const int n2 = d2->size();

      // advance sp in 2ndary signal 
      while ( 1 )
	{

	  if ( p2 < 0 || p2 >= n2 )
	    d1[ p1 ] = 0 ;
	  else
	    d1[ p1 ] = (*d2)[ p2 ];
	  
	  // advance
	  ++p1;
	  ++p2;
	  
	  // all done?
	  if ( p1 == n1 ) break;
	  
	}

      // add the new signal
      edf.add_signal( sig , Fs , d1 );
      
      // add any annotations 

    }
  
}



  
