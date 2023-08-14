
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

#include "annot/align-annots.h"

#include "edf/edf.h"
#include "annot/annot.h"

#include "db/db.h"
#include "helper/logger.h"
#include "helper/helper.h"

extern writer_t writer;
extern logger_t logger;

align_annots_t::align_annots_t( edf_t & edf , param_t & param )
{

  // full (original) EDF is attached
  // requires an existing ALIGN-EPOCHS solution based on this EDF
  //  take the epoch-mapping, and alter all annotations in the implied
  //  new EDF;  save 


  //
  // Get all annotations from original 
  //
  
  std::vector<interval_t> epochs;
  int ne = edf.timeline.first_epoch();
  while ( 1 ) {
    int epoch = edf.timeline.next_epoch();
    if ( epoch == -1 ) break;
    interval_t interval = edf.timeline.epoch( epoch );
    epochs.push_back( interval );
  }


  //
  // Get prior solution 
  //

  // prior solution: only scan for this ID
  const std::string solfile = Helper::expand( param.requires( "sol" ) );
  
  if ( ! Helper::fileExists( solfile ) ) Helper::halt( "could not open " + solfile );

  // assumes the following header structure: (but check) 
  // ID	E2	D	E1	NEXT	NEXT_E1	ORDERED

  std::ifstream IN1( solfile.c_str() , std::ios::in );
  
  std::string line;
  Helper::safe_getline( IN1 , line );
  if ( IN1.eof() ) Helper::halt( "problem reading " + solfile );
  
  std::vector<std::string> tok = Helper::parse( line ) ;
  bool okay = true;
  if ( tok.size() != 7 ) okay = false;
  else if ( tok[0] != "ID" ) okay = false;
  else if ( tok[1] != "E2" ) okay = false;
  else if ( tok[3] != "E1" ) okay = false;
  else if ( tok[6] != "ORDERED" ) okay = false;

  if ( ! okay )
      Helper::halt( "expecting 7-tab delim cols: ID, E2, D, E1, NEXT, NEXT_E1 & ORDERED" );

  // by default, warning if ORDERED is 0 
    
  std::map<int,int> emap;
  while ( 1 )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() ) break;
      if ( line == "" ) continue;
      std::vector<std::string> tok = Helper::parse( line ) ;
      if ( tok.size() != 7 ) Helper::halt( "bad line: " + line );

      const std::string id = tok[0];

      int e1, e2;
      if ( ! Helper::str2int( tok[1] , &e2 ) )
	Helper::halt( "expecting integer epoch codes: " + line );
      if ( ! Helper::str2int( tok[3] , &e1 ) )
	Helper::halt( "expecting integer epoch codes: " + line );

      if ( e1 > ne )
	Helper::halt( "expecting original epoch codes between 1 and " + Helper::int2str( ne ) );
      
      // from (edf) -> to (edf2)
      if ( e1 < 1 || e2 < 1 )
	Helper::halt( "expecting 1-based epoch codes: " + line );

      // save as 0-based epoch numbers
      emap[ e1 - 1 ] = e2 - 1 ;
      
    }

  IN1.close();

  logger << "  read " << emap.size() << " epochs to remap from " << solfile << "\n";
  

  //
  // nb. FOR NOW, assume SIMPLE epch codes (i.e. offset = 0 in particular...) 
  //

  edf_t edf2;

  // not that it really matters, but make dummy EDF the "size" of 30 * number of saved epochs
  const int nr = emap.size() * globals::default_epoch_len ; 
  const double rs = 1;
  const std::string startdate = "01.01.85";
  const std::string starttime = "00.00.00";
  
  okay = edf2.init_empty( "__dummy" , nr , rs , startdate , starttime );
  if ( ! okay ) Helper::halt( "internal problem generating EDF" );

  //
  // Make new annotations
  //

  std::vector<std::string> names = edf.timeline.annotations.names();

  uint64_t elen = globals::tp_1sec * globals::default_epoch_len ; 

  std::map<std::string,int> skipped;
  
  for (int a=0;a<names.size();a++)
    {
      
      annot_t * annot1 = edf.timeline.annotations.find( names[a] );
      annot_t * annot2 = edf2.timeline.annotations.add( names[a] );
      
      // get all events
      annot_map_t::const_iterator aa = annot1->interval_events.begin();
      while ( aa != annot1->interval_events.end() )
        {

	  instance_idx_t instance_idx = aa->first;
	  interval_t interval1 = instance_idx.interval ;
	  
	  // figure out which epochs start/end occurred in 
	  // do -1 for stop
	  
	  int e1 = interval1.start / elen;
	  int e2 = (interval1.stop-1LLU) / elen; 

	  bool f1 = emap.find( e1 ) != emap.end();
	  bool f2 = emap.find( e2 ) != emap.end();

	  if ( ! ( f1 && f2 ) )
	    {
	      ++skipped[ names[a] ];
	      ++aa;
	      continue;
	    }
	  
	  int n1 = emap[ e1 ];
	  int n2 = emap[ e2 ];
	  
	  // check same span
	  if ( n2 - n1 != e2 - e1 )
	    {
              ++skipped[ names[a] ];
              ++aa;
              continue;
            }
	  
	  // good to convert
	  //  ( can drop the -1 now) 
	  uint64_t offset1 = interval1.start - e1 * elen;
	  uint64_t offset2 = interval1.stop  - e2 * elen;

	  // new mapping
	  interval_t interval2( n1 * elen + offset1 , n2 * elen + offset2 );
	  
	  // add new annot - note, drops any meta-data currently - easy to add if needed
	  annot2->add( instance_idx.id , interval2 , instance_idx.ch_str );
	  
          ++aa;
        }
      
    } // next annot class

       
  logger << "  copied " << names.size() << " new annotation classes\n";

  
  //
  // Save new annotations
  //

  edf2.timeline.annotations.write( param.requires( "out" ) , param , edf2 );

  
}
