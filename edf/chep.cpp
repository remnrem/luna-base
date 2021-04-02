
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


#include "helper/helper.h"
#include "helper/logger.h"
#include "edf.h"
#include "eval.h"
#include "artifacts/artifacts.h"

extern logger_t logger;


void proc_chep_mask( edf_t & edf , param_t & param )
{

  // fixed things first (clipped signals, flat or max values)
  // if param does not contain any relevant options, this will just 
  // return

  chep_mask_fixed( edf , param );

  // iterative procedure for Hjorth parameters
  chep_mask( edf , param );

}

void timeline_t::proc_chep( edf_t & edf , param_t & param )
{


  //
  // get a (default) signals list 
  //

  signal_list_t signals = edf.header.signal_list( param.value( "sig" ) );

  //
  // requires epoched data 
  //

  if ( ! edf.timeline.epoched() ) 
    Helper::halt( "data not EPOCH'ed");

  //
  // reset
  //

  if ( param.has( "clear" ) ) 
    edf.timeline.clear_chep_mask();


  //
  // read
  //
  
  if ( param.has( "load" ) )
    {
      std::string f = param.value( "load" );
      logger << "  reading chep from " << f << "\n";
      edf.timeline.read_chep_file( f );
    }


  //
  // prune to include only epochs/channels that actually, currently exist
  //

  std::set<int> epochs;
  const int ne = edf.timeline.first_epoch();
  while ( 1 ) 
    {
      int epoch = edf.timeline.next_epoch_ignoring_mask();
      if ( epoch == -1 ) break;
      // track all display epochs (whether currently masked or not)
      epochs.insert( edf.timeline.display_epoch( epoch ) );
    }

  std::set<std::string> channels;
  for (int s=0;s<edf.header.label.size();s++)
    {
      if ( edf.header.is_data_channel( s ) )
	channels.insert( edf.header.label[s] );
    }

  // copy/clear/reset chep mask
  std::map<int,std::set<std::string> > copy = edf.timeline.chep;
  edf.timeline.clear_chep_mask();

  std::map<int,std::set<std::string> >::const_iterator ii = copy.begin();
  while ( ii != copy.end() )
    {

      // got this epoch?
    if ( epochs.find( ii->first ) != epochs.end() )
      {
	std::set<std::string>::const_iterator jj = ii->second.begin();
	while ( jj != ii->second.end() )
	  {
	    if ( channels.find( *jj ) != channels.end() )
	      edf.timeline.chep[ ii->first ].insert( *jj );
	    ++jj;
	  }
      }
    ++ii;
  }

  //
  // manually specify good/bad channels/epochs
  //

  if ( param.has( "bad-channels" ) ) 
    {

      std::string sigstr = param.requires( "bad-channels" );
      signal_list_t bad_signals = edf.header.signal_list( sigstr );

      const int ns = bad_signals.size();
      const int ne = edf.timeline.num_epochs();
      
      if ( ns > 0 ) logger << "  setting as bad channels:";
      for (int s=0;s<ns;s++)
	{
	  logger << " " << bad_signals.label(s);
	  for (int e=0;e<ne;e++)
	    edf.timeline.set_chep_mask( e , bad_signals.label(s) );
	}
      if ( ns > 0 ) logger << "\n";
    }

  
  //
  // collapse epochs & set epoch-mask (restricted to signals_list from 'sig' option)
  //

  if ( param.has( "epochs" ) ) 
    {
      
      // epochs=prop{,n}
      // pct=0 means remove any epoch w/ 1 or more (x > p ) 
      // n=0 means ( x >= n )
      
      // defaults
      std::vector<double> p = { 0.0 , 0.0 } ;
      
      if ( param.value( "epochs" ) != "T" ) // i.e. if a specific value given
	p = param.dblvector( "epochs" );

      if ( p.size() == 1 ) 
	edf.timeline.collapse_chep2epoch( signals , p[0] , 0 ); 
      else if ( p.size() == 2 ) 
	edf.timeline.collapse_chep2epoch( signals , p[0] , p[1] ); 
      else 
	Helper::halt( "expecting with 0, 1 or 2 args: epochs=prop{,n}" );
    }


  //
  // collapse channels (just alters CH/EP, OR drop channel from EDF)
  //

  if ( param.has( "channels" ) || param.has( "drop-channels" ) ) 
    {
      
      bool drop = param.has( "drop-channels" );

      // channels=pct{,n}
      //  pct=0  remove any channel with any bad data
      //  n=0    do not remove 
      
      // defaults
      std::vector<double> p = { 0.0 , 0.0 } ;
      
      if ( drop && param.value( "drop-channels" ) != "T" ) // i.e. if a specific value given
	p = param.dblvector( "drop-channels" );      
      else if ( (!drop) && param.value( "channels" ) != "T" ) // i.e. if a specific value given
	p = param.dblvector( "channels" );
      
      if ( p.size() == 1 ) p.push_back(0);
      else if ( p.size() != 2 ) Helper::halt( "expecting channels or drop-channels with 0, 1 or 2 args" );

      // just alter CHEP mask (i.e. blank out whole row)
      // OR, actually drop the channel from the EDF?
      
      if ( drop )
	{

	  signal_list_t drops = edf.timeline.collapse_chep2ch( signals , p[0] , 0 , true , false ); 
	  
	  if ( drops.size() > 0 ) logger << "  dropping channels:";

	  std::set<std::string> labels;
	  for ( int s=0; s < drops.size(); s++) labels.insert( drops.label(s) );
	  
	  std::set<std::string>::const_iterator dd = labels.begin();
	  while ( dd != labels.end() )
	    {
	      logger << " " << *dd;
	      int s = edf.header.signal( *dd );
	      if ( s != -1 ) edf.drop_signal( s );
	      ++dd;
	    }
	  
	  if ( drops.size() > 0 ) logger << "\n";
	}
      else 
	{
	  // just alter CHEP mask, retain signal (i.e. for INTERPOLATE)
	  
	  // by default:
	  // true, false implies -- set all bad channels to all bad:  YES
	  //                        set all good channels to all good: NO
	  
	  // but let this second option be modified
	  bool black_and_white = param.has( "black-and-white" );
	  edf.timeline.collapse_chep2ch( signals , p[0] , p[1] , true , black_and_white ); 
	}

    }


  //
  // dump to standard output mechanism
  //

  if ( param.has( "dump" ) ) 
    edf.timeline.dump_chep_mask( signals , true );  // info to console plus output
  else
    edf.timeline.dump_chep_mask( signals , false );   // info to console only
  

  //
  // write to a file
  //

  if ( param.has( "save" ) )
    {
      std::string f = param.value( "save" );
      logger << "  saving chep to " << f << "\n";
      edf.timeline.write_chep_file( f );
    }

  
}

