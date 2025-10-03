
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

#include "dsp/pcoupl.h"
#include "param.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "dsp/hilbert.h"
#include "db/db.h"

extern writer_t writer;
extern logger_t logger;

void dsptools::phase_coupling( edf_t & edf , param_t & param )
{
  
  // slow signals

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );  
  const int ns = signals.size();  
  if ( ns == 0 ) return;
  
  // events
  if ( ! param.has( "events" ) ) Helper::halt( "requires events to point to one or more annotation classes" );
  std::set<std::string> evts = param.strset( "events" );
  if ( evts.size() == 0 ) return;

  // anchor : 0, 1, 2:  start(default) mid end
  int anchor = 0;
  if ( param.has( "anchor" ) )
    {
      if ( param.value( "anchor" ) == "mid" || param.value( "anchor" ) == "middle" ) anchor = 1;
      else if ( param.value( "anchor" ) == "end" || param.value( "anchor" ) == "stop" ) anchor = 2;
    }
  
  // permutations
  const int nreps = param.has( "nreps" ) ? param.requires_int( "nreps" ) : 1000 ; 
  if ( nreps < 0 ) Helper::halt( "nreps must be non-negative" );

  // default is to permute within intervals
  const bool epoch_perm = param.has( "perm-whole-trace" ) ? false : true; 
  
  // epoch size
  double epoch_sec = 0;
  if ( epoch_perm )
    {
      // if epoch-based permutation w/ generic epochs, we need to define a fixed epoch length      
      if ( edf.timeline.generic_epochs() )
	{
	  if ( param.has( "fixed-epoch-dur" ) )
	    epoch_sec = param.requires_dbl( "fixed-epoch-dur" );
	  else
	    Helper::halt( "cannot run within-epoch permutation with generic epochs: add 'fixed-epoch-dur or perm-whole-trace" );
	}
      else
	{
	  if ( ! edf.timeline.epoched() ) 
	    edf.timeline.ensure_epoched();
	  epoch_sec = edf.timeline.epoch_length();	  	  
	}
      logger << "  using epoch duration of " << epoch_sec << "s for within-epoch shuffling\n";
    }
  

  //
  // slow phases
  //

  const double phase_lwr = param.requires_dbl( "lwr" );
  const double phase_upr = param.requires_dbl( "upr" );
  const double fir_ripple = param.has( "ripple" ) ? param.requires_dbl( "ripple" ) : 0.01 ;
  const double fir_tw = param.has( "tw" ) ? param.requires_dbl( "tw" ) : 0.5 ;
  
  
  //
  // for each channel
  //
  
  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  for (int s=0; s<ns; s++)
    {
      // skip annot channels
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      writer.level( signals.label(s) , globals::signal_strat );

      logger << "  processing " << signals.label(s) << "\n";

      // get all data
      slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );
      	
      // signals
      const std::vector<double> * slow = slice.pdata();
      
      // time-points
      const std::vector<uint64_t> * tps = slice.ptimepoints();
      
      // filter-Hilbert
      hilbert_t hilbert( *slow , Fs[s] , phase_lwr , phase_upr , fir_ripple , fir_tw );
      
      logger << "  done filter-Hilbert...\n";


      //
      // verbose output?
      //

      const bool verbose = param.has( "verbose" );
      
      
      //
      // for each event
      //

      std::set<std::string>::const_iterator ee = evts.begin();
      while ( ee != evts.end() )
	{
	  
	  //
	  // get time-points
	  //

	  annot_t * annot = edf.annotations->find( *ee );
	  if ( annot == NULL ) { ++ee; continue; }
	  
	  writer.level( *ee  , globals::annot_strat );

	  logger << "   - processing " << *ee << "\n";

	  // convert annotation time-points to sample-pointns for this signal;
	  // if the anchor is more than one sample-rate interval, implies it is 
	  // in a gap or off the grid, so ignore those.
	  
	  std::set<uint64_t> evt_tp;
	  
	  const annot_map_t & events = annot->interval_events;
	  
	  annot_map_t::const_iterator aa = events.begin();
	  while ( aa != events.end() )
	    {	      
	      const interval_t & interval = aa->first.interval;
	      
	      // get anchor
	      uint64_t atp = anchor == 0 ? interval.start : ( anchor == 2 ? interval.stop - 1LLU : interval.mid() ) ;
	      
	      // add to sorted list
	      evt_tp.insert( atp );
	      
	      ++aa;
	    }
	  
	  //
	  // convert tps -> sps;   this might be a discontinous EDF, so we can't just 
	  // scale by SR;  other ways to do, but given sorted just do something dumb here
	  //
	  
	  const int np = tps->size();
	  
	  const uint64_t delta = globals::tp_1sec / Fs[s] ; 
	  
	  std::vector<int> evt_sp;
	  
	  int idx = 0;
	  int cnt = 0;
	  std::set<uint64_t>::const_iterator ii = evt_tp.begin();
	  while ( ii != evt_tp.end() )
	    {
	      // all done?
	      if ( idx == np ) break;
	      
	      // not yet?
	      if ( *ii > (*tps)[idx] && *ii - delta > delta )
		{
		  ++idx; // advance candidate tp, stick w/ same event
		  continue;
		}
	      
	      // passed by ?
	      if ( (*tps)[idx] > *ii && (*tps)[idx] - *ii > delta ) 
		{
		  ++ii; // advance event
		  continue;
		}
	      	      
	      // else implies a match
	      evt_sp.push_back( idx );	      
	      ++cnt;

	      // track
	      if ( verbose ) 
		{
		  writer.level( cnt , globals::count_strat );
		  //writer.value( "T" , (*tps)[idx] );
		  //writer.value( "SP" , idx );		  
		}

	      // move to next event
	      ++ii;
	    }
	  
	  logger << "  mapped " << evt_sp.size() << " of " << evt_tp.size() << " events\n";
	  
	  if ( verbose ) 
	    writer.unlevel( globals::count_strat );	  
	  
	  //
	  // coupling analysis
	  //

	  itpc_t itpc = hilbert.phase_events( evt_sp , 
					      NULL , // no mask
					      nreps ,
					      Fs[s] ,
					      epoch_sec ,  // opt: within-epoch shuffle
					      true         // stratify by slow signal phase bin
					      );
	  
	  
	  //
	  // Outputs
	  //
	  
	  writer.value( "N" , (int)evt_sp.size() );
	  writer.value( "MAG" , itpc.itpc.obs );
	  writer.value( "MAG_EMP" , itpc.itpc.p );
	  writer.value( "MAG_NULL" , itpc.itpc.mean );
	  writer.value( "MAG_Z" , ( itpc.itpc.obs - itpc.itpc.mean ) / itpc.itpc.sd );
	  
	  if ( itpc.angle.obs > -9 )
	    writer.value( "ANGLE" , itpc.angle.obs );

	  //
	  // asymptotic significance of coupling test; under
	  // the null, give mean rate of 'significant'
	  // (P<0.05) coupling
	  //
	  
	  writer.value( "PV" , itpc.pv.obs ) ;
	  if ( nreps )
	    writer.value( "SIGPV_NULL" , itpc.sig.mean ) ;
	  	  
	  //
	  // phase bins
	  //
	  
	  if ( nreps )	       
	    {
	      const int nbins = 18;
	      for (int b = 0 ; b < nbins ; b++ ) 
		{
		  writer.level( b * 20 + 10 , "PHASE" );
		  writer.value( "OVERLAP"      , itpc.phasebin[b].obs );
		  writer.value( "OVERLAP_EXP"   , itpc.phasebin[b].mean );
		  writer.value( "OVERLAP_EMP"    , itpc.phasebin[b].p );
		  
		  if ( itpc.phasebin[b].sd > 0 ) 
		    {
		      double Z = ( itpc.phasebin[b].obs - itpc.phasebin[b].mean ) / itpc.phasebin[b].sd ; 
		      writer.value( "OVERLAP_Z" , Z );
		    }
		}
	      writer.unlevel( "PHASE" );
	    }
	  
	  ++ee;
	} // next annotation					  
      
      writer.unlevel( globals::annot_strat );
      
    } // next signal
  
  writer.unlevel( globals::signal_strat );
  
  // all done
}






