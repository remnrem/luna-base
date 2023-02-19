
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

#include "edf/freezer.h"
#include "edf/edf.h"
#include "helper/logger.h"
#include "eval.h"

extern logger_t logger;


void freezer_t::edf2edf( const edf_t & from , edf_t & to )
{

  // primary shallow copy
  to = from;
  
  // update edf_t pointers in all records
  //to.update_edf_pointers( & to );
  
  logger << "  copied " << from.header.nr << " records\n";
  
  // we'll keep annotations constant, as those are not changed.
  // this means if we add new annotations, those are kept across the whole recording

  // we will not be able to change record size -- i.e. as that sets 'problem' flag and
  // skips -- but this is a reasonable constraint.
  
  
}
  
void freezer_t::freeze( const std::string & s , edf_t & edf )
{
  
  logger << "  freezing state, with tag " << s << "\n";

  //
  // first, ensure that all records have been pulled from the file;
  // i.e. so we will not need to touch the disk again;  this way,
  // we don't care about tracking the file handles (FILE * or edfz * pointers)
  // 
  
  // this command should pull all unread records into memory, but
  // will 
  //  a) skip records that are no longer retained,
  //  b) not re-read already loaded records
  
  edf.read_records( 0 , edf.header.nr_all - 1 );

  //
  // Allocate freeze target
  //
  
  edf_t * edf2 = new edf_t;
  
  //
  // Do the copy (includes deep-copying of some items)
  //

  edf2edf( edf , *edf2 );

  //
  // store pointer to this edf_t
  //
  
  store[ s ] = edf2;

  //
  // some output
  //

  logger << "  currently " << store.size() << " freeze(s):";
  std::map<std::string,edf_t*>::const_iterator ss = store.begin();
  while ( ss != store.end() ) { logger << " " << ss->first; ++ss; }
  logger << "\n";
  
}
  
bool freezer_t::thaw( const std::string & s , edf_t * edf , bool also_clean )
{

  if ( store.find( s ) == store.end() )
    Helper::halt( "could not find frozen EDF " + s );
  
  logger << "  thawing previous freeze " << s << "\n";
  
  edf_t * new_edf = store[s];

  logger << "  old dataset   : "
	 << edf->header.nr << " records, "
	 << edf->header.ns << " signals, "
	 << edf->timeline.annotations.names().size() << " annotations\n";
  logger << "  thawed freeze : "
	 << new_edf->header.nr << " records, "
	 << new_edf->header.ns << " signals, "
	 << new_edf->timeline.annotations.names().size() << " annotations\n";
  
  // do the copy back
  edf2edf( *(store[s]) , *edf );
  
  // and delete also?
  if ( also_clean ) clean( s );

  // it is not possible to freeze an empty dataset,
  // therefore set globals::empty to false (i.e. is okay
  // if it was previously set), to allow processing to restart

  globals::empty = false;
  
  return true;
}

void freezer_t::clean( const std::string & s )
{
  if ( store.find( s ) != store.end() )
    {
      logger << "  cleaning up freeze " << s << "\n";
      edf_t * p = store[ s ];
      delete p; 
    }
}


