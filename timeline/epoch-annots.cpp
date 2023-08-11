
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

#include "timeline.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger;

void timeline_t::annotate_epochs( const std::string & label , 
				  const std::string & annot_label , 
				  const std::set<std::string> & values )
{


  //
  // Take information from the annot_t class, and make a simple
  // per-epoch annotation this can be performed after a restructure,
  // but (total) # of epochs must match the file exactly
  //


  //
  // Point to first epoch, and get the 'total' number of epochs
  // (i.e. both masked and unmasked), as first_epoch() only returns
  // the unmasked counts; nb. first_epoch() actually sets the pointer
  // *before* the first epoch, so fine to use this whether actual
  // first epoch is masked or not
  //
  
  first_epoch();
  
  int ne = num_total_epochs();
  
  //
  // Populate epoch-annotation vectors to the appropriate size
  //
  
  eannots[ label ].clear();
  

  //
  // Get annotations
  //
  
  annot_t * annot = annotations( annot_label );

  // if not found, then all eannots are effectively false
  // (i.e. missing)

  if ( annot == NULL ) return;


  //
  // for each epoch 
  //
  
  while ( 1 ) 
    {
      
      //
      // Get next epoch
      //
      
      int e = next_epoch_ignoring_mask();      

      if ( e == -1 ) break;

      // use 'e' for look-up here; but need to use orginal EDF
      // encoding for eannots, internally i.e. the zero-based version
      
      int e0 = original_epoch( e ) ;

      if ( e0 == -1 ) 
	Helper::halt( "internal error in annotate_epochs()" );

      interval_t interval = epoch( e );
      
      annot_map_t events = annot->extract( interval );
      
      // search for a matching value (at least one)
      
      annot_map_t::const_iterator ii = events.begin();

      while ( ii != events.end() )
	{	
	  
	  const instance_idx_t & instance_idx = ii->first;
	  const instance_t * instance = ii->second;

	  if ( values.find( instance_idx.id ) != values.end() )
	    {	      
	      // nb. store w.r.t. original epoch encoding e0
	      eannots[ label ][ e0 ] = true;
	      break;
	    }	      
	  
	  ++ii;
	  
	}
      
    } // next epoch
}

// should be used with current 0..(ne-1) mapping, will
// be converted to original epoch mapping if needed

void timeline_t::annotate_epoch( const std::string & label , int e )
{
  
  // do we need to remap the epoch?
  
  if ( has_epoch_mapping() )
    {
      // off-the-grid  
      if ( epoch_curr2orig.find( e ) == epoch_curr2orig.end() ) 
	return;
      
      // convert query to the original mapping
      e = epoch_curr2orig.find( e )->second;
    }
  
  eannots[ label ][ e ] = true;
}

  
void timeline_t::clear_epoch_annotations()
{
  if ( eannots.size() > 0 ) 
    logger << "  clearing all epoch-annotations\n";
  eannots.clear();
}



// Return all epoch annotations

std::set<std::string> timeline_t::epoch_annotations() const
{
  std::set<std::string> r;
  std::map<std::string,std::map<int,bool> >::const_iterator ii = eannots.begin();
  while ( ii != eannots.end() )
    {
      r.insert( ii->first );
      ++ii;
    }
  return r;
}

// does annotation 'k' exist at all? 

bool timeline_t::epoch_annotation(const std::string & k ) const
{
  return eannots.find( k ) != eannots.end() ;
}

// does EPOCH 'e' contain annotation 'k'?
// where 'e' is in the current 0..ne epoch form, 
// and will be remapped if necessary

bool timeline_t::epoch_annotation(const std::string & k, int e) const
{
  
  // look up this annotation 'k'
  std::map<std::string,std::map<int,bool> >::const_iterator ii = eannots.find( k );
  
  // annotation k does not exist anywhere
  if ( ii == eannots.end() ) return false;
  
  // do we need to remap the epoch? 
  if ( has_epoch_mapping() ) 
    {
      // off-the-grid
      if ( epoch_curr2orig.find( e ) == epoch_curr2orig.end() ) return false;
      // convert query to the original mapping
      e = epoch_curr2orig.find( e )->second;
    }
  
  // now we have the correct original-EDF epoch number, do
  // we have an flag? if not, means FALSE
  if ( ii->second.find( e ) == ii->second.end() ) return false;
  
  // return the boolean value (i.e. could still be set FALSE explicitly)
  return ii->second.find( e )->second; 
}
