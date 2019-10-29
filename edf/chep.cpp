
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

extern logger_t logger;

void proc_chep( edf_t & edf , param_t & param )
{

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
  // manually specify good/bad channels/epochs
  //

  if ( param.has( "bad-channels" ) ) 
    {

      std::string sigstr = param.requires( "bad-channels" );
      signal_list_t bad_signals = edf.header.signal_list( sigstr );
      const int ns = bad_signals.size();
      const int ne = edf.timeline.num_epochs();

      for (int s=0;s<ns;s++)
	{
	  logger << "  setting CHEP mask for " << bad_signals.label(s) << "\n";
	  for (int e=0;e<ne;e++)
	    edf.timeline.set_chep_mask( e , bad_signals(s) );
	}
    }


  //
  // collapse epochs
  //

  if ( param.has( "epochs" ) ) 
    {

      std::string sigstr = param.requires( "sig" );
      signal_list_t signals = edf.header.signal_list( sigstr );

      std::vector<double> p = param.dblvector( "epochs" );
      if ( p.size() == 1 ) 
	edf.timeline.collapse_chep2epoch( signals , p[0] , 0 ); 
      else if ( p.size() == 2 ) 
	edf.timeline.collapse_chep2epoch( signals , p[0] , p[1] ); 
      
    }


  //
  // collapse channels
  //

  if ( param.has( "channels" ) ) 
    {

      std::string sigstr = param.requires( "sig" );
      signal_list_t signals = edf.header.signal_list( sigstr );
      
      std::vector<double> p = param.dblvector( "channels" );
      if ( p.size() == 1 ) 
	edf.timeline.collapse_chep2ch( signals , p[0] , 0 ); 
      else if ( p.size() == 2 ) 
	edf.timeline.collapse_chep2ch( signals , p[0] , p[1] ); 
      
    }
  

  //
  // dump to standard output mechanism
  //

  if ( param.has( "dump" ) ) 
    edf.timeline.dump_chep_mask();  

  
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

