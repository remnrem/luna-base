
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
#include "db/db.h"

extern logger_t logger;

extern writer_t writer;

void freezer_t::edf2edf( const edf_t & from , edf_t & to , bool preserve_cache )
{

  caches_t cache;
  std::string c_num, c_int, c_str;
  
  // do not overwrite cache? 
  if ( preserve_cache )
    {
      logger << "  restoring the cache and any recording options\n";
      cache = to.timeline.cache;
      c_num = writer.cache_num_name();
      c_int = writer.cache_int_name();
      c_str = writer.cache_str_name();
      //std::cout << " " << c_num << "] [" << c_int << "] [" << c_str << "\n";
    }
  else
    {
      // should revisit this - cannot simply set, as otherwise FREEZE wipes the cache
      //   think it implies that CACHE record should be set prior to any freeze/thaw...
      // TODO... at some point, revisit and fix it up..
      // logger << "  wiping any cache variables\n";
      // writer.no_cache();
    }

  // std::cout << " edf2edf_t...\n";
  // std::cout << " from annots N = " << from.annotations->names().size() << "\n";
  // std::cout << " to annots N = " << to.annotations->names().size() << "\n\n";
  

  // primary shallow copy
  to = from;
  
  // swap original cache back in
  if ( preserve_cache )
    {
      // restore actual cache
      to.timeline.cache = cache;

      // restore any writer/db recordings
      if ( c_num != "" )
	{	  
	  cache_t<double> * cache = to.timeline.cache.find_num( c_num );	  
          writer.cache( cache );
	}

      if ( c_int != "" )
	{
	  cache_t<int> * cache = to.timeline.cache.find_int( c_int );
          writer.cache( cache );
	}

      if ( c_str != "" )
	{
	  cache_t<std::string> * cache = to.timeline.cache.find_str( c_str );
          writer.cache( cache );
	}

    }

	
  // ??unnecessary now we've ensured all records are pulled into memory...
  // update edf_t pointers in all records
  to.update_edf_pointers( & to );

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
  // we should now be able to shut off the file inpts safely, if any exist
  //
  
  edf.closeout_inputs();
  
  //
  // Allocate freeze target (but point to same annotations)
  //  (nb. the shallow copy will copy edf.annotations anyway)
  //
  
  edf_t * edf2 = new edf_t( edf.annotations );

  
  //
  // Make a copy of edf into edf2
  //

  edf2edf( edf , *edf2 , false );



  
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
  
bool freezer_t::thaw( const std::string & s , edf_t * edf , bool also_clean , bool preserve_cache )
{
  
  if ( store.find( s ) == store.end() )
    Helper::halt( "could not find frozen EDF " + s );
  
  logger << "  thawing previous freeze " << s << "\n";
  
  edf_t * new_edf = store[s];

  logger << "  old dataset   : "
	 << edf->header.nr << " records, "
	 << edf->header.ns << " signals, "
	 << edf->annotations->names().size() << " annotations\n";
  logger << "  thawed freeze : "
	 << new_edf->header.nr << " records, "
	 << new_edf->header.ns << " signals, "
	 << new_edf->annotations->names().size() << " annotations\n";
  
  // do the copy back
  edf2edf( *(store[s]) , *edf , preserve_cache );
  
  // and delete also?
  if ( also_clean ) clean( s , true );
  
  // it is not possible to freeze an empty dataset,
  // therefore set globals::empty to false (i.e. is okay
  // if it was previously set), to allow processing to restart

  globals::empty = false;
  
  return true;
}

void freezer_t::clean( edf_t * self )
{
  if ( store.size() == 0 )
    return;
  
  std::map<std::string,edf_t*>::const_iterator ss = store.begin();
  while ( ss != store.end() )
    {
      // make sure not to clean self
      if ( ss->second != self )
	clean( ss->first , false ); // keep store key intact while we iterate
      ++ss;
      
    }
  store.clear();

}
  
  
void freezer_t::clean( const std::string & s , bool clear_store )
{

  if ( store.find( s ) != store.end() )
    {
      logger << "  cleaning up freeze " << s << "\n";
      edf_t * p = store[ s ];
      delete p;
    }
  
  if ( clear_store )
    store.erase( s );
}



