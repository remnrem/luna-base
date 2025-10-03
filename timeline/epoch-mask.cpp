
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


#include "timeline/timeline.h"
#include "param.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "miscmath/crandom.h"
#include "helper/token-eval.h"

#include "db/db.h"
#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger;

//
// Epoch Masks
//


bool timeline_t::masked_epoch( int e ) const
{
  if ( ! mask_set ) return false;
  if ( e < 0 || e >= mask.size() ) return true; // out-of-bounds, so 'masked'
  return mask[e];
}

// flip all values of a mask
// i.e. to /include/ artifactual epochs only

void timeline_t::flip_epoch_mask()
{
  
  if ( ! mask_set ) return;

  const int ne = epochs.size();
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  
  // flip all (i.e. every epoch will change)
  for (int e=0;e<ne;e++)
    {
      
      mask[e] = ! mask[e];
      
      if ( mask[e] ) ++cnt_mask_set;
      else ++cnt_mask_unset;

    }
  
  logger << "  flipped all epoch masks\n";
  logger << "  total of " << cnt_mask_unset << " of " << epochs.size() << " retained\n";

}


// behavior if annotation missing, i.e. an 'empty' mask

// void timeline_t::apply_empty_epoch_mask( const std::string & label , bool include )
// {
  
//   // this is requested if the annotation is missing
//   // i.e. returns match == F for every epoch; treat as specified by include and mask_mode

//   // include T/F   means whether a 'match' means having (T) versus not-having (F) the annotation
//   // mask_mode will already have been set
  
//   mask_set = true;
  
//   const int ne = epochs.size();
  
//   int cnt_mask_set = 0;
//   int cnt_mask_unset = 0;
//   int cnt_unchanged = 0;
//   int cnt_now_unmasked = 0;
//   int cnt_basic_match = 0;  // basic count of matches, whether changes mask or not
  
//   for (int e=0;e<ne;e++)
//     {
      
//       bool matches = false;
      
//       // set new potential mask, depending on match_mode
      
//       bool new_mask = mask[e];

//       if ( include ) 
// 	{
// 	  if      ( mask_mode == 0 ) new_mask = matches;   // mask-if
// 	  else if ( mask_mode == 1 ) new_mask = !matches;  // unmask-if
// 	  else if ( mask_mode == 2 ) new_mask = matches ;  // if
// 	}
//       else
// 	{
// 	  if      ( mask_mode == 0 ) new_mask = !matches;  // mask-ifnot
// 	  else if ( mask_mode == 1 ) new_mask = matches;   // unmask-ifnot
// 	  else if ( mask_mode == 2 ) new_mask = ! matches; // ifnot
// 	}

//       int mc = set_epoch_mask( e , new_mask );

//       if      ( mc == +1 ) ++cnt_mask_set;
//       else if ( mc == -1 ) ++cnt_mask_unset;
//       else                 ++cnt_unchanged;
      
//       if ( !mask[e] ) ++cnt_now_unmasked;
      
//     }

//   logger << "  based on " << label << " " << cnt_basic_match << " epochs match; ";
//   logger << cnt_mask_set << " newly masked, "   
// 	 << cnt_mask_unset << " unmasked, " 
// 	 << cnt_unchanged << " unchanged\n";
//   logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";
  
//   // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
//   writer.level( label  , "EMASK" );

//   writer.var( "N_MATCHES"    , "Number of matching epochs" );
//   writer.var( "N_MASK_SET"   , "Number of epochs newly masked" ); 
//   writer.var( "N_MASK_UNSET" , "Number of epochs newly unmasked" );
//   writer.var( "N_UNCHANGED"  , "Number of epochs unchanged by this mask" );
//   writer.var( "N_RETAINED"   , "Number of epochs retained for analysis" );
//   writer.var( "N_TOTAL"      , "Total number of epochs" );

//   writer.value( "N_MATCHES"    , cnt_basic_match  );
//   writer.value( "N_MASK_SET"   , cnt_mask_set     );
//   writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
//   writer.value( "N_UNCHANGED"  , cnt_unchanged    );
//   writer.value( "N_RETAINED"   , cnt_now_unmasked );
//   writer.value( "N_TOTAL"      , (int)epochs.size()    );

//   writer.unlevel( "EMASK" );

// }


// void timeline_t::apply_epoch_mask( annot_t * a , std::set<std::string> * values , bool include )
// {
  
//   // include T/F   means whether a 'match' means having (T) versus not-having (F) the annotation
  
//   // mask_mode will already have been set
  
//   // if 'values' is NULL, then we just use presence of an annotation,
//   // rather than looking at the instance ID
  
//   bool value_mask = values != NULL;
  
//   mask_set = true;

//   const int ne = epochs.size();
  
//   //
//   // We do not clear the mask here, as we want to allow multiple
//   // filters to be added on top of oneanther
//   //

//   int cnt_mask_set = 0;
//   int cnt_mask_unset = 0;
//   int cnt_unchanged = 0;
//   int cnt_now_unmasked = 0;
//   int cnt_basic_match = 0;  // basic count of matches, whether changes mask or not
  
//   for (int e=0;e<ne;e++)
//     {

//       interval_t interval = epoch( e );
      
//       annot_map_t events = a->extract( interval );
      
//       bool matches = false;
      
//       if ( value_mask ) 
// 	{
// 	  // do any of the instance IDs match any of the values?
// 	  annot_map_t::const_iterator ii = events.begin();
// 	  while ( ii != events.end() )
// 	    {		  
// 	      const instance_idx_t & instance_idx = ii->first;	      
// 	      if ( values->find( instance_idx.id ) != values->end() )
// 		{
// 		  matches = true;
// 		  break;
// 		}
// 	      ++ii;
// 	    }
// 	}
//       else 
// 	{	  
// 	  matches = events.size() > 0 ;
// 	}

//       // count basic matches

//       if ( matches ) ++cnt_basic_match;
      
//       // set new potential mask, depending on match_mode
      
//       bool new_mask = mask[e];

//       if ( include ) 
// 	{
// 	  if      ( mask_mode == 0 ) new_mask = matches;   // mask-if
// 	  else if ( mask_mode == 1 ) new_mask = !matches;  // unmask-if
// 	  else if ( mask_mode == 2 ) new_mask = matches ;  // if
// 	}
//       else
// 	{
// 	  if      ( mask_mode == 0 ) new_mask = !matches;  // mask-ifnot
// 	  else if ( mask_mode == 1 ) new_mask = matches;   // unmask-ifnot
// 	  else if ( mask_mode == 2 ) new_mask = ! matches; // ifnot
// 	}

//       int mc = set_epoch_mask( e , new_mask );

//       if      ( mc == +1 ) ++cnt_mask_set;
//       else if ( mc == -1 ) ++cnt_mask_unset;
//       else                 ++cnt_unchanged;
      
//       if ( !mask[e] ) ++cnt_now_unmasked;
      
//     }
  
//   logger << "  based on " << a->name << ( value_mask ? "[" + Helper::stringize( *values , "|" ) + "]" : "" )  
// 	 << " " << cnt_basic_match << " epochs match; ";

//   logger << cnt_mask_set << " newly masked, " 
// 	 << cnt_mask_unset << " unmasked, " 
// 	 << cnt_unchanged << " unchanged\n";
//   logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

  
//   // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
//   writer.level( a->name , "EMASK" );

//   writer.var( "N_MATCHES"    , "Number of matching epochs" );
//   writer.var( "N_MASK_SET"   , "Number of epochs newly masked" ); 
//   writer.var( "N_MASK_UNSET" , "Number of epochs newly unmasked" );
//   writer.var( "N_UNCHANGED"  , "Number of epochs unchanged by this mask" );
//   writer.var( "N_RETAINED"   , "Number of epochs retained for analysis" );
//   writer.var( "N_TOTAL"      , "Total number of epochs" );

//   writer.value( "N_MATCHES"    , cnt_basic_match  );
//   writer.value( "N_MASK_SET"   , cnt_mask_set     );
//   writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
//   writer.value( "N_UNCHANGED"  , cnt_unchanged    );
//   writer.value( "N_RETAINED"   , cnt_now_unmasked );
//   writer.value( "N_TOTAL"      , (int)epochs.size()    );

//   writer.unlevel( "EMASK" );
// }





void timeline_t::apply_epoch_mask2( const std::map<annot_t *,std::set<std::string> > & annots ,
				    const std::set<std::string> & fullspan , 
				    const std::string & alabel ,
				    bool or_match , bool include )
				    
{

  // annots[] contains annot pointers
  //  if the value is non-null, implies to match on instance IDs
  //  otherwise, match on the class ID

  //  if annot is null, implies an 'empty' mask - i.e. no instances of that requested annotation found
  //    (previously we handled this as a special case in apply_empty_epoch_mask()
  //    simply eval as "F", i.e. no match, and then treat as usual

  //  if multiple annotations specified, a match is defined as AT LEAST ONE - i.e. OR
  //   AND logic can be specified 

  //   or_match   T    epoch matches if AT LEAST ONE annotation matches (i.e. OR mode)
  //              F    epoch matches if ALL annotations match (i.e. AND mode)

  //   include    T    a match means having that annotation requirement
  //              F    a match means not having that annotation requirement
  
  // Note: we do not clear the mask here, as we want to allow multiple
  // filters to be added on top of oneanther

  //
  // Initialise
  //

  mask_set = true;
  
  const int ne = epochs.size();

  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  int cnt_basic_match = 0;  // basic count of matches, whether changes mask or not

  //
  // Iterate over epochs
  //

  for (int e=0;e<ne;e++)
    {
      
      // for and/or logic
      //   i.e. AND matches = n_matches == n_annots
      //        OR  matches = n_matches > 0 ;    
      
      // note that AND logic is at the class level: annotations are always
      // considered as OR matches 
      
      int n_annots = annots.size();
      int n_matches = 0;

      // get the epoch 
      interval_t interval = epoch( e );
      
      // consider all annots
      std::map<annot_t *,std::set<std::string> >::const_iterator aa = annots.begin();
      while ( aa != annots.end() )
	{

	  // get the annotation
	  annot_t * a = aa->first;
	  
	  // empty? i.e. no instances, implies this epoch will not match either
	  if ( a == NULL )
	    {
	      ++aa;
	      continue;
	    }
	  
	  // require this epoch fully-spanned (by a each annotation of this class) to be a match?
	  const bool full = fullspan.find( aa->first->name ) != fullspan.end();
	  
	  // otherwise, get overlapping annots
	  annot_map_t events = full ? a->extract_complete_overlap( interval ) : a->extract( interval );

	  const bool check_instance_value = aa->second.size() != 0; 

	  if ( check_instance_value ) 
	    {
	      // do any of the instance IDs match any of the values?
	      // check each of potentially multiple events, but otherwise stop
	      // at the first match
	      annot_map_t::const_iterator ii = events.begin();
	      while ( ii != events.end() )
		{		  
		  const instance_idx_t & instance_idx = ii->first;	      
		  if ( aa->second.find( instance_idx.id ) != aa->second.end() )
		    {
		      // at least one value matches: we can count this as a class-level
		      // match and move on to the next annot class
		      ++n_matches;
		      break;
		    }
		  ++ii;
		}
	    }
	  else 
	    {
	      // otherwise, simply note if we see one or 1 in this area: 
	      if ( events.size() != 0 ) ++n_matches;
	    }
	  
	  // early stopping if on OR mode?
	  if ( or_match && n_matches > 0 ) break;
	  
	  // next annotation
	  ++aa;
	}
      
      // do we get a match?
      //  can apply either OR or AND logic across multiple annot clases here
      
      const bool matches = or_match ? n_matches != 0 : n_matches == n_annots ; 
      
      // count basic matches
      
      if ( matches ) ++cnt_basic_match;
      
      // set new potential mask, depending on mask mode and include/exclde mode
      
      bool new_mask = mask[e];

      if ( include ) 
	{
	  if      ( mask_mode == 0 ) new_mask = matches;   // mask-if
	  else if ( mask_mode == 1 ) new_mask = !matches;  // unmask-if
	  else if ( mask_mode == 2 ) new_mask = matches ;  // if
	}
      else
	{
	  if      ( mask_mode == 0 ) new_mask = !matches;  // mask-ifnot
	  else if ( mask_mode == 1 ) new_mask = matches;   // unmask-ifnot
	  else if ( mask_mode == 2 ) new_mask = ! matches; // ifnot
	}
      
      int mc = set_epoch_mask( e , new_mask );
      
      if      ( mc == +1 ) ++cnt_mask_set;
      else if ( mc == -1 ) ++cnt_mask_unset;
      else                 ++cnt_unchanged;
      
      if ( !mask[e] ) ++cnt_now_unmasked;
      
    }

  logger << "  applied annotation mask for " << annots.size() << " annotation(s)";
  if ( annots.size() > 1 ) logger << " (using " << ( or_match ? "or" : "and") << "-matching across multiple annotations)";
  logger << "\n";

  logger << "  " << cnt_basic_match << " epochs match; ";
  
  logger << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

  
  // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
  writer.level( alabel , "EMASK" );
  
  writer.var( "N_MATCHES"    , "Number of matching epochs" );
  writer.var( "N_MASK_SET"   , "Number of epochs newly masked" ); 
  writer.var( "N_MASK_UNSET" , "Number of epochs newly unmasked" );
  writer.var( "N_UNCHANGED"  , "Number of epochs unchanged by this mask" );
  writer.var( "N_RETAINED"   , "Number of epochs retained for analysis" );
  writer.var( "N_TOTAL"      , "Total number of epochs" );
  
  writer.value( "N_MATCHES"    , cnt_basic_match  );
  writer.value( "N_MASK_SET"   , cnt_mask_set     );
  writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
  writer.value( "N_UNCHANGED"  , cnt_unchanged    );
  writer.value( "N_RETAINED"   , cnt_now_unmasked );
  writer.value( "N_TOTAL"      , (int)epochs.size()    );
  writer.value( "MATCH_LOGIC"  , or_match ? "OR" : "AND" );
  writer.value( "MATCH_TYPE"   , include ? "if" : "ifnot" );
  writer.value( "MASK_MODE"    , mask_mode == 2 ? "force" : ( mask_mode == 0 ? "mask" : "unmask" ) ) ; 
 		
  writer.unlevel( "EMASK" );

  
}




  
// convert from full annotation to epoch-level masks
    
// void timeline_t::apply_epoch_include_mask( annot_t * a , std::set<std::string> * values )
// {
//   apply_epoch_mask( a , values , true );
// }

// void timeline_t::apply_epoch_exclude_mask( annot_t * a , std::set<std::string> * values )
// {
//   apply_epoch_mask( a , values , false );
// }



  // eval-based mask
// eval-based mask
void timeline_t::apply_eval_mask( const std::string & str , int mask_mode , const bool verbose )
{

  // mask_mode   0   mask
  //             1   unmask
  //             2   force  (T = mask    &  F = unmask)   [ drop mode ] 
  //            -2   force  (T = unmask  &  F = mask)     [ keep mode ]
  

  // is this 'KEEP' mode? 

  bool flip = false;
  
  if ( mask_mode == - 2 ) { mask_mode = 2 ; flip = true; } 


  //
  // Set mask mode
  //

  if ( mask_mode > -1 ) 
    {
      set_epoch_mask_mode( mask_mode );  
      logger << "  set masking mode to " << ( mask_mode == 2 ? "'force'" : mask_mode == 1 ? "'unmask'" : "'mask' (default)" ) << "\n";
    }


  //
  // allow both " and # quoting of EVAL expressions
  //

  std::string expression = Helper::trim( Helper::unquote( str , '#' ) );
  

  //
  // Get all existing annotations (overkill...)
  //
  
  std::vector<std::string> names = annotations->names();


  //
  // Keep track of changes
  //
  
  mask_set = true;

  const int ne = epochs.size();
  

  //
  // We do not clear the mask here, as we want to allow multiple
  // filters to be added on top of oneanther
  //

  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  int cnt_basic_match = 0;  

  
  //
  // Iterate over epochs
  //

  first_epoch();
  
  int acc_total = 0 , acc_retval = 0 , acc_valid = 0; 
  
  while ( 1 ) 
    {
      
      int e = next_epoch_ignoring_mask() ;

      if ( e == -1 ) break;
      
      interval_t interval = epoch( e );
	  
      std::map<std::string,annot_map_t> inputs;
      
      // get each annotations
      for (int a=0;a<names.size();a++)
	{
	  
	  annot_t * annot = annotations->find( names[a] );
	  
	  // get overlapping annotations for this epoch
	  annot_map_t events = annot->extract( interval );
	  
	  // store
	  inputs[ names[a] ] = events;
	}

      //
      // create a dummy new instance for the output variables (not saved)
      //
      
      instance_t dummy;
      
      //
      // evaluate the expression, but note, this is set to not 
      // allow any assignments.... this makes it cleaner and easier 
      // to spot bad//undefined variables as errors.
      //

      const bool no_assignments = true;
      
      Eval tok( expression , no_assignments );

      tok.bind( inputs , &dummy );

      bool is_valid = tok.evaluate( verbose );
      
      bool matches;
      
      if ( ! tok.value( matches ) ) is_valid = false;



      //
      // Flip?
      //

      if ( flip ) matches = ! matches;
            
      
      //
      // A match must be a valid value
      //
      
      if ( ! is_valid ) matches = false;


      //
      // apply mask (or not)
      //
      
      acc_total++;

      acc_valid += is_valid;
      
      //
      // Only for valud results
      //

      if ( is_valid ) 	
	{
	  acc_retval += matches;
      
	  
	  // count basic matches
	  
	  if ( matches ) ++cnt_basic_match;
      
	  // set new potential mask, depending on match_mode
	  
	  bool new_mask = mask[e];
	  
	  if      ( mask_mode == 0 ) new_mask = matches;   // mask
	  else if ( mask_mode == 1 ) new_mask = !matches;  // unmask
	  else if ( mask_mode == 2 ) new_mask = matches ;  // mask/unmask
	  
	  int mc = set_epoch_mask( e , new_mask );
	  
	  if      ( mc == +1 ) ++cnt_mask_set;
	  else if ( mc == -1 ) ++cnt_mask_unset;
	  else                 ++cnt_unchanged;
	}
      else ++cnt_unchanged;

      // count current mask state

      if ( !mask[e] ) ++cnt_now_unmasked;
      
      // next epoch
    } 
  
  
  logger << "  based on eval expression [" << expression << "]\n"
	 << "  " << acc_retval << " true, " << acc_valid - acc_retval << " false and " 
	 << acc_total - acc_valid << " invalid return values\n"
	 << "  " << cnt_basic_match << " epochs match; " 
	 << cnt_mask_set << " newly masked, "
	 << cnt_mask_unset << " unmasked, "
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

  
  // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
  writer.level( expression , "EMASK" );
  
  writer.var( "N_MATCHES"    , "Number of matching epochs" );
  writer.var( "N_MASK_SET"   , "Number of epochs newly masked" ); 
  writer.var( "N_MASK_UNSET" , "Number of epochs newly unmasked" );
  writer.var( "N_UNCHANGED"  , "Number of epochs unchanged by this mask" );
  writer.var( "N_RETAINED"   , "Number of epochs retained for analysis" );
  writer.var( "N_TOTAL"      , "Total number of epochs" );

  writer.value( "N_MATCHES"    , cnt_basic_match  );
  writer.value( "N_MASK_SET"   , cnt_mask_set     );
  writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
  writer.value( "N_UNCHANGED"  , cnt_unchanged    );
  writer.value( "N_RETAINED"   , cnt_now_unmasked );
  writer.value( "N_TOTAL"      , (int)epochs.size()    );

  writer.unlevel( "EMASK" );


  // all done 
  return;

}


  // other masks : randomly select up to 'n' epochs from the current set 
  // other masks : randomly select up to 'n' epochs from the current set 
void timeline_t::select_epoch_randomly( int n )
{

  mask_set = true;
  
  // from the unmasked set, pick at random 'n' (or as many as possible)
  std::vector<int> unmasked;
  
  const int ne = epochs.size();
  
  for (int e=0;e<ne;e++)
    if ( ! mask[e] ) unmasked.push_back(e);
  
  std::set<int> selected;
  
  int s = 0;
  
  const int num_unmasked = unmasked.size();
  
  const int n_to_select = num_unmasked < n ? num_unmasked : n;

  while ( s < n_to_select )
    {
      int rnd = CRandom::rand( num_unmasked );
      int sel = unmasked[ rnd ];

      if ( selected.find( sel ) == selected.end() )
	{
	  selected.insert( sel );
	  ++s;
	}
    }
  
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;


  // mask everything that was /not/ in the selected set
  for (int e=0;e<ne;e++)
    {
      
      if ( selected.find(e) == selected.end() ) 
	{
	  int mc = set_epoch_mask( e , true );
	  if ( mc == +1 ) ++cnt_mask_set;
	  else if ( mc == -1 ) ++cnt_mask_unset;
	  else ++cnt_unchanged;
	}
      
      if ( ! mask[e] ) ++cnt_now_unmasked;
    }

  logger << "  randomly selected up to " << n << " epochs; ";

  logger << cnt_mask_set << " newly masked " 
	 << cnt_mask_unset << " unmasked and " 
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

}


// trim leading and trailing epochs (allow only n)
void timeline_t::trim_epochs( std::string & label , int n )
{

  // find leading and trailing epochs with this 'label' annotation (e.g. "W")
  // and allow only up to 'n'
  // i.e. trimming from the start/end of the epochs

  annot_t * annot = annotations->find( Helper::unquote( label ) );

  if ( annot == NULL ) return;

  mask_set = true;
  
  const int ne = epochs.size();

  std::vector<bool> x( ne , false );
  
  for (int e=0;e<ne;e++)
    {
      interval_t interval = epoch( e );
      annot_map_t events = annot->extract( interval );
      x[e] = events.size() > 0 ;
    }
  
  // find first non-matching epoch
  // -1 if no leading matches
  int leading_end = -1;
  for (int e=0;e<ne;e++)
    {      
      if ( ! x[e] ) 
	{
	  leading_end = e - 1;
	  break;
	}
    }
  
  // find last nonmatching
  int trailing_start = ne;
  // ne if no trailing matches
  for (int e=ne-1;e>=0;e--)
    {      
      if ( ! x[e] ) 
	{	  
	  trailing_start = e + 1;
	  break;
	}
    }
  
  // allow up to 'n' trailing/leading matches
  leading_end -= n;
  trailing_start += n;

  if ( leading_end > 0 ) logger << "  trimming from start to epoch " << leading_end + 1 << "\n";
  if ( trailing_start < ne-1 ) logger << "  trimming from epoch " << trailing_start + 1 << " to end\n";

  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  int cnt_basic_match = 0;  // basic count of matches, whether changes mask or not                                                                                
  // blank out any ones needed
  for ( int e=0; e<ne; e++)
    {
      if ( e <= leading_end || e >= trailing_start )
	{
	  ++cnt_basic_match;
	  
	  // set new potential mask, depending on match_mode
	  
	  bool new_mask = true;
	  
	  int mc = set_epoch_mask( e , new_mask );
	  
	  if      ( mc == +1 ) ++cnt_mask_set;
	  else if ( mc == -1 ) ++cnt_mask_unset;
	  else                 ++cnt_unchanged;
	}
    
      if ( !mask[e] ) ++cnt_now_unmasked;
    
    }
  
  logger << "  based on leading/trailing " << label << " (w/ up to " << n << " epochs) " 
	 << cnt_basic_match << " epochs match; ";
  
  logger << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";
  
  // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
  writer.level( label , "EMASK" );
  writer.value( "N_MATCHES"    , cnt_basic_match  );
  writer.value( "N_MASK_SET"   , cnt_mask_set     );
  writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
  writer.value( "N_UNCHANGED"  , cnt_unchanged    );
  writer.value( "N_RETAINED"   , cnt_now_unmasked );
  writer.value( "N_TOTAL"      , (int)epochs.size()    );

  writer.unlevel( "EMASK" );

}



void timeline_t::mask_leading_trailing( const std::vector<std::string> & anns ,
					const bool _leading , const bool _trailing, const bool _mask )
{

  const std::string label = Helper::stringize( anns );
  
  // versus trim (above)
  //   this only does leading or trailing currently (as called)
  //   you can mask on multiple annots
  //   you can either mask or unmask those regions
  //   you cannot specify 'allow N'

  // find first/last annot has /doesn't/ have one of 'anns'
  //  --> this will delimit the interval to mask/unmask

  const int ne = epochs.size();

  std::vector<bool> x( ne , false );
  
  for (int a=0; a<anns.size(); a++)
    {
      
      annot_t * annot = annotations->find( Helper::unquote( anns[a] ) );
      
      if ( annot == NULL ) continue;
      
      mask_set = true;      
      
      for (int e=0;e<ne;e++)
	{
	  if ( ! x[e] ) { 
	    interval_t interval = epoch( e );
	    annot_map_t events = annot->extract( interval );
	    if ( events.size() > 0 ) x[e] = true;
	  }
	}
      
    } // next annot


  // find first/last epoch that does not match
  
  // find first non-matching epoch
  // -1 if no leading matches
  int leading_end = -1;
  if ( _leading ) 
    for (int e=0;e<ne;e++)
      {      
	if ( ! x[e] ) 
	  {
	    leading_end = e - 1;
	    break;
	  }
      }
  
  // find last nonmatching
  int trailing_start = ne;
  // ne if no trailing matches
  if ( _trailing ) 
    for (int e=ne-1;e>=0;e--)
      {      
	if ( ! x[e] ) 
	  {	  
	    trailing_start = e + 1;
	    break;
	  }
      }
  
  if ( leading_end > 0 ) logger << "  trimming from start to epoch " << leading_end + 1 << "\n";
  if ( trailing_start < ne-1 ) logger << "  trimming from epoch " << trailing_start + 1 << " to end\n";
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  int cnt_basic_match = 0;  // basic count of matches, whether changes mask or not                                                                                
  // blank out any ones needed
  for ( int e=0; e<ne; e++)
    {
      // if 'mask' mode (default, set mask for leading/trailing regions
      // if not, retain those, & set mask for interim region
      const bool match = _mask
	? ( e <= leading_end || e >= trailing_start )
	: ( e > leading_end && e < trailing_start ) ;
      
      if ( match )
	{
	  ++cnt_basic_match;
	  
	  // set new potential mask, depending on match_mode
	  
	  bool new_mask = true;
	  
	  int mc = set_epoch_mask( e , new_mask );
	  
	  if      ( mc == +1 ) ++cnt_mask_set;
	  else if ( mc == -1 ) ++cnt_mask_unset;
	  else                 ++cnt_unchanged;
	}
    
      if ( !mask[e] ) ++cnt_now_unmasked;
    
    }
  
  logger << "  based on leading/trailing " << Helper::stringize( anns ) << " "
	 << cnt_basic_match << " epochs match; ";
  
  logger << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";
  
  // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
  writer.level( label , "EMASK" );
  writer.value( "N_MATCHES"    , cnt_basic_match  );
  writer.value( "N_MASK_SET"   , cnt_mask_set     );
  writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
  writer.value( "N_UNCHANGED"  , cnt_unchanged    );
  writer.value( "N_RETAINED"   , cnt_now_unmasked );
  writer.value( "N_TOTAL"      , (int)epochs.size()    );

  writer.unlevel( "EMASK" );


  
}

// retain contiguous stretch of epochs that include these (+/- n epochs on both sides)
void timeline_t::retain_epochs( const std::set<std::string> & labels )
{

  // find first and last examples epochs w/ one of these labels (e.g. a sleep label)
  // expand outwards n epochs (if possible) 
  // mask leading and trailing gaps
  //  if no matches, then blank all
  
  // if one label is +int then interpret as n epoch border, not as a label
  int n = 0;

  int first = -1;
  int last = -1;
  
  const int ne = epochs.size();

  std::set<std::string>::const_iterator aa = labels.begin();
  while ( aa != labels.end() )
    {
      
      // matches +n ?
      if ( aa->substr(0,1) == "+" )
	{
	  if ( Helper::str2int( aa->substr(1) , &n ) ) 
	    {
	      ++aa;
	      continue;
	    }
	}
      
      // treat as annot
      annot_t * annot = annotations->find( Helper::unquote( *aa ) );
  
      if ( annot == NULL ) { ++aa; continue; }
      
      // find first match
      // -1 if no leading matches                                                                                                                       
      int leading_end = -1;
      for (int e=0;e<ne;e++)
	{
	  // already gone past a prior mark?
	  if ( e >= first && first != -1 ) break;
	  
	  interval_t interval = epoch( e );
	  annot_map_t events = annot->extract( interval );
	  if ( events.size() > 0 ) 
	    {
	      first = e; 
	      break;
	    }
	}
      
      // find final match
      int trailing_start = ne;      
      for (int e=ne-1;e>=0;e--)
	{      
	  // already gone past a prior mark?                                                                                                       
          if ( e <= last && last != -1 ) break;
	  interval_t interval = epoch( e );
          annot_map_t events = annot->extract( interval );
          if ( events.size() > 0 )
            {
              last = e;
              break;
            }
	}

      // look to next annotation now
      ++aa;
    }
  
  //
  // any matches? if not, mask all epochs
  //

  const bool no_matches = first == -1;

  
  // else mark up to 'first' (not including) and after 'last' (not including)
  
  mask_set = true;
  
  // adjust by up to 'n' epochs
  
  first = first > n ? first - n : 0 ; 
  last = (last + n ) < ne ? last + ne : ne - 1 ;

  if ( no_matches ) logger << "  did not find any matching epochs, masking all epochs\n";
  else logger << "  retaining from epoch " << first 
	      << " to " << last 
	      << " ( " << last - first + 1 << " epochs ), given up to " << n << " epoch border\n";
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  int cnt_basic_match = 0;  // basic count of matches, whether changes mask or not                                                                                
  // blank out any ones needed
  for ( int e=0; e<ne; e++)
    {
      if ( no_matches || e < first || e > last )
	{
	  ++cnt_basic_match;
	  
	  // set new potential mask, depending on match_mode
	  
	  bool new_mask = true;
	  
	  int mc = set_epoch_mask( e , new_mask );
	  
	  if      ( mc == +1 ) ++cnt_mask_set;
	  else if ( mc == -1 ) ++cnt_mask_unset;
	  else                 ++cnt_unchanged;
	}
      
      if ( !mask[e] ) ++cnt_now_unmasked;
    
    }
  
  std::string lstr = Helper::stringize( labels );

  logger << "  based on retained region " << lstr << " (w/ up to " << n << " epochs) " 
	 << cnt_basic_match << " epochs match; ";
  
  logger << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";
  
  // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
  writer.level( lstr , "EMASK" );
  writer.value( "N_MATCHES"    , cnt_basic_match  );
  writer.value( "N_MASK_SET"   , cnt_mask_set     );
  writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
  writer.value( "N_UNCHANGED"  , cnt_unchanged    );
  writer.value( "N_RETAINED"   , cnt_now_unmasked );
  writer.value( "N_TOTAL"      , (int)epochs.size()    );

  writer.unlevel( "EMASK" );

}



void timeline_t::unmask_interior()
{
  // look back / forward y epochs;
  // require that at least x of y are 'good' on at least one side
  // if both sides fail criterion, set this epoch to be masked
  // nb. out-of-range counts as 'bad'  (i.e. not a 'good' epoch)
  
  const int ne = epochs.size();
  
  // find first and last unmasked epochs
  int first_masked = 0;
  int last_masked = ne-1;
  for (int e=0;e<ne;e++)
    {
      if ( ! mask[e] ) 
	{
	  first_masked = e;
	  break;
	}	
    }

  for (int e=ne-1;e>=0;e--)
    {
      if ( ! mask[e] )
	{
          last_masked = e;
          break;
        }
    }
  
  std::set<int> tounmask;
  for (int e=first_masked; e<=last_masked; e++)
    {
      if ( mask[e] ) 
	{
	  tounmask.insert(e);
	}
    }
 
  //
  // now we have list of epochs that need masking in 'tomask'
  //
  
  std::set<int>::const_iterator ee = tounmask.begin();
  while ( ee != tounmask.end() )
    {
      mask[*ee] = false; 
      ++ee;
    }
  
  int cnt_now_unmasked = 0;
  for (int e=0;e<ne;e++)
    if ( ! mask[e] ) ++cnt_now_unmasked; 
  
  logger << "  based on unmask-interior: ";
  logger << tounmask.size() << " newly unmasked epochs\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

}


// regional mask
void timeline_t::regional_mask( int x , int y )
{

  // look back / forward y epochs;
  // require that at least x of y are 'good' on at least one side
  // if both sides fail criterion, set this epoch to be masked
  // nb. out-of-range counts as 'bad'  (i.e. not a 'good' epoch)

  if ( y < 1 || x > y || x < 1 ) 
    Helper::halt( "invalid values for regional mask" );
  
  const int ne = epochs.size();
  
  std::set<int> tomask;
  std::vector<bool> putative_mask( ne , false );
  
  for (int e=0;e<ne;e++)
    {
      // nothing to do
      if ( mask[e] ) 
	{
	  putative_mask[e] = true;
	  continue;
	}

      int backward = 0 , forward = 0;
      
      int curr = e;
      for (int i=0;i<y;i++)
	{
	  --curr;
	  if ( curr < 0 ) break;
	  if ( ! mask[ curr ] ) ++backward;
	}

      curr = e;
      for (int i=0;i<y;i++)
	{
	  ++curr;
	  if ( curr == ne ) break;
	  if ( ! mask[ curr ] ) ++forward;
	}
      
      // a bad outcome?
      if ( ! ( forward >= x || backward >= x ) )
	{
	  tomask.insert( e );
	  putative_mask[e] = true;
	}

    }
  
  // additionally, do not allow any 'island' bad masks 
  // which are possible given the simple heuristic above
  // i.e. for any epoch flanked by two bad epochs, mask

  for (int e=0;e<ne;e++)
    {
      if ( putative_mask[e] ) continue;
      int bad = 0;
      if ( e == 0 || putative_mask[e-1] ) ++bad;
      if ( e == ne-1 || putative_mask[e+1] ) ++bad;
      if ( bad == 2 ) tomask.insert( e );
    }

  //
  // now we have list of epochs that need masking in 'tomask'
  //

  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;

  std::set<int>::const_iterator ee = tomask.begin();
  while ( ee != tomask.end() )
    {
      int mc = set_epoch_mask( *ee , true );
      if ( mc == +1 ) ++cnt_mask_set;
      else if ( mc == -1 ) ++cnt_mask_unset;
      else ++cnt_unchanged;
      ++ee;
    }
  
  for (int e=0;e<ne;e++)
    if ( ! mask[e] ) ++cnt_now_unmasked; 

  logger << "  based on regional smoothing ("<<x << "/" << y << " good), ";

  logger << cnt_mask_set << " newly masked " 
	 << cnt_mask_unset << " unmasked and " 
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";
  
}

// mask range of epochs from a to b inclusive
void timeline_t::select_epoch_range( int a , int b , bool include )
{
  std::set<int> e;
  if ( a > b ) 
    {
      int tmp = b;
      b = a;
      a = tmp;
    }

  for (int i=a; i<=b; i++) e.insert( i );

  if ( include )
    logger << "  selecting epochs from " << a << " to " << b << "; ";
  else
    logger << "  masking epochs from " << a << " to " << b << "; ";

  return select_epoch_range( e , include );

}


// other masks : select up to 'n' epochs from the start of the record
void timeline_t::select_epoch_first( int n )
{
  
  mask_set = true;
  
  // from the unmasked set, pick at random 'n' (or as many as possible)
  //  std::vector<int> unmasked;
  
  const int ne = epochs.size();
     
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;

  // mask everything that was /not/ in the selected set
  for (int e=0;e<ne;e++)
    {
      if ( e >= n )
	{
	  int mc = set_epoch_mask( e , true );
	  if ( mc == +1 ) ++cnt_mask_set;
	  else if ( mc == -1 ) ++cnt_mask_unset;
	  else ++cnt_unchanged;
	}
      
      if ( ! mask[e] ) ++cnt_now_unmasked;
    }

  logger << "  selecting up to " << n << " epochs for start; ";
  logger << cnt_mask_set << " newly masked, "  
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";
}


  // mask set of epochs
void timeline_t::select_epoch_range( const std::set<int> & specified_epochs , bool include )
{

  mask_set = true;

  const int ne = epochs.size();
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  
  // mask everything that was /not/ in the selected set
  for (int e=0;e<ne;e++)
    {
      // use base-1 coding of epochs
      const int epoch = e+1;

      bool inset = specified_epochs.find( epoch ) != specified_epochs.end() ;

      bool match = include ? 
	! inset : inset ;
      
      if ( match ) 
	{
	  int mc = set_epoch_mask( e , true );
	  if ( mc == +1 ) ++cnt_mask_set;
	  else if ( mc == -1 ) ++cnt_mask_unset;
	  else ++cnt_unchanged;
	}
      
      if ( ! mask[e] ) ++cnt_now_unmasked;
    }

  if ( include )
    logger << "  selecting";
  else
    logger << "  masking";

  logger << " from set of " << specified_epochs.size() << " epochs; ";
  
  logger << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";

  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

}


bool timeline_t::masked_interval( const interval_t & interval , bool all_masked , bool * start_masked ) const
{
  
  // if all_masked,   returns T if /all/ of interval falls within masked regions
  // if not,          returns T if interval falls in at least one masked region
  // allow 0-duration intervals
    
  if ( start_masked != NULL ) *start_masked = false;  
  
  if ( edf->header.continuous )
    {

      if ( ! mask_set ) 
	{
	  return false;
	}
      
      int eleft = MiscMath::position2leftepoch( interval.start , epoch_length_tp, epoch_inc_tp , mask.size() );

      // end of interval is one past end of region:
      int eright = MiscMath::position2rightepoch( interval.stop-1LLU , epoch_length_tp, epoch_inc_tp , mask.size() );
      
      //      std::cout << "e1e2 = " << eleft << "  " << eright << "\n";
      
      if ( start_masked != NULL )
	{
	  if ( eleft == -1 || mask[eleft] ) *start_masked = true;
	}

      // i.e. no spanning epochs defined -- same as being 'masked'
      if ( eleft == -1 || eright == -1 ) return true;
      
      // above functions return -1 if the tp is off the map
      // (or epochs are not overlapping/contiguous); here it is
      // effectively 'masked'
      
      for (int e=eleft;e<=eright;e++) 
	{	  
	  if ( all_masked && ! mask[e] ) return false;
	  if ( mask[e] && ! all_masked ) return true;
	}

    }

  else // for EDF+D
    {
      
      // if ( ! mask_set )
      // 	{
      //     return false;
      // 	}

          
      std::set<int> records = records_in_interval( interval );

      // falls off edge of the map
      if ( records.size() == 0 ) return true;
      
      std::set<int>::const_iterator rr = records.begin();
      while ( rr != records.end() )
	{

	  const std::set<int> & epochs = rec2epoch.find( *rr )->second;
	  
	  std::set<int>::const_iterator ee = epochs.begin();

	  // no epochs defined over these records - i.e. same as 'masked'
	  if ( epochs.size() == 0 ) return true;
	  
	  if ( start_masked != NULL )
	    {
	      if ( mask[ *ee ] ) *start_masked = true;
	    }
	  
	  while ( ee != epochs.end() )
	    {	      	      
	      if ( all_masked && ! mask[ *ee ] ) return false;
	      if ( mask[ *ee ] && ! all_masked ) return true;
	      ++ee;
	    }
	  ++rr;
	}      
    }
  
  
  if ( all_masked )
    return true;
  else
    return false;

}


// select all EPOCHs until we come across an EPOCH that does /not/ have the 'str' annotation
void timeline_t::select_epoch_until_isnot( const std::string & str )
{
  Helper::halt( "timeline_t::select_epoch_until_isnot() is no longer supported" );

  mask_set = true;
  
  const int ne = epochs.size();

  bool found = false;

  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;

  for (int e=0;e<ne;e++)
    {      
      bool a = epoch_annotation( str , e );	  
      if ( ! a ) found = true;	  
      
      int mc = set_epoch_mask( e , found );
      if ( mc == +1 ) ++cnt_mask_set;
      else if ( mc == -1 ) ++cnt_mask_unset;
      else ++cnt_unchanged;
      
      if ( ! mask[e] ) ++cnt_now_unmasked;

    }

  logger << "  based on " << str << " leading epochs; ";

  logger << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

}


  // select EPOCHs within a run of 'n' (either side) other similarly annotated epochs
// select only EPOCHs that are in contiguous runs of EPOCH /str/ (i.e. +1 means one either side)

void timeline_t::select_epoch_within_run( const std::string & label , int b )
{

  if ( b < 1 ) Helper::halt( "epoch border must be 1 or greater" );

  annot_t * annot = annotations->find( Helper::unquote( label ) );

  const bool no_annots = annot == NULL ;

  mask_set = true;
  
  // get epoch annots for this label:
  const int ne = epochs.size();
  std::vector<bool> x( ne , false );

  if ( ! no_annots )
    for (int e=0;e<ne;e++)
      {
	interval_t interval = epoch( e );
	annot_map_t events = annot->extract( interval );
	x[e] = events.size() > 0 ;
      }
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;

  for (int e=0;e<ne;e++)
    {  
      
      bool set_mask = false;

      // if ( ! epoch_annotation( str , e ) ) 
      // 	set_mask = true;
      if ( ! x[e] )
	set_mask = true;
      
      if ( ! set_mask )
	{	  
	  
	  int cnt = 0;
	  
	  int current = e;
	  for (int bwk=0;bwk<b;bwk++)
	    {
	      if ( current == 0 ) continue;
	      --current;	      
	      //if ( epoch_annotation( str , current ) ) ++cnt;
	      if ( x[ current] ) ++cnt;
	    }
	  
	  current = e;
	  for (int fwd=0;fwd<b;fwd++)
	    {
	      if ( current == ne-1 ) continue;
	      ++current;
	      //if ( epoch_annotation( str , current ) ) ++cnt;
	      if ( x[ current ] ) ++cnt;
	    }
	  	  
	  if ( cnt < b * 2 ) set_mask = true;

	}
      
      int mc = set_epoch_mask( e , set_mask );
      if ( mc == +1 ) ++cnt_mask_set;
      else if ( mc == -1 ) ++cnt_mask_unset;
      else ++cnt_unchanged;
      
      if ( ! mask[e] ) ++cnt_now_unmasked;
      
    }
  
  logger << "  based on " << label << " with " << b << " flanking epochs; ";
  logger << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

}


void timeline_t::stable_mask( const int x , const std::vector<std::string> & anns , const bool unique )
{

  // if unique: can only have 1 of anns
  //   only select epochs contain a member of 'anns' that don't have an epoch within
  //   x epochs containing some other member of anns

  // if not unique:
  //   looking for runs of /any/ of that type
  
  // e.g. 2,B,C  - 
  // AAAAABBBBBBBCCCCCCBBBBABBB
  //      XX|||XXXX||XX
  // i.e. requires a buffer of 2 epochs
  //      before hitting another epoch
  
  std::map<std::string,annot_t*> annots;
  for (int a=0; a<anns.size(); a++)
    {
      annot_t * annot = annotations->find( Helper::unquote( anns[a] ) );
      if ( annot != NULL ) annots[ anns[a] ] = annot;
    }
  
  mask_set = true;

  // get epoch annots for this label:

  // collapse to 'any annot' if not-unique mode
  const int na = unique ? anns.size() : 1 ; 
  int na_obs = 0;
  const int ne = epochs.size();
  
  std::vector<bool> xany( ne , false ); // for non-unique mode
  std::vector<std::vector<bool> > XA( na );
  for (int a=0; a<anns.size(); a++)
    {
      std::vector<bool> x1( ne , false ); // none found
      if ( annots.find( anns[a] ) != annots.end() )
	{
	  ++na_obs;
	  annot_t * annot = annots[ anns[a] ] ;
	  for (int e=0;e<ne;e++)
	    {
	      interval_t interval = epoch( e );
	      annot_map_t events = annot->extract( interval );
	      x1[e] = events.size() > 0 ;

	      // track 'at least one' for non-unique mode
	      if ( ! unique )
		if ( x1[e] )
		  xany[e] = true;
	    }
	}

      // collapse if not unique mode
      if ( unique ) 
	XA[ a ] = x1 ;
      
    }
  
  // if not-unique mode, update XA w/ xany
  if ( ! unique ) 
    XA[ 0 ] = xany;
  
  // we should now have a ne x na matrix, where T means we see an annot
  // now determine (given x) which 
  std::vector<bool> X( ne , false );

  for (int e=0; e<ne; e++)
    {
      // a) we require one and only one annot
      int okay = 0;
      int a1 = 0;

      for (int a=0; a<na; a++)
	{
	  if ( XA[a][e] ) {
	    ++okay;
	    a1 = a;
	  }
	}

      if ( okay == 1 ) 
	{
	  // now check flanking regions - require x up, x back
	  // for which we see only a1 and no others
	  
	  // 0 1 2 3 4 5 6 7   ne = 8
	  //     .     .
	  // e.g. if x = 2
	  // require at least x epochs before after
	  if ( e < x  || e > ne - 1 - x ) 
	    {
	      okay = 0;
	    }
	  
	  if ( okay != 0 )
	    {
	      // now check annots
	      for (int f=1; f<=x; f++)
		{
		  // must still be 'A1'
		  if ( ! ( XA[a1][e-f] && XA[a1][e+f] ) ) { okay = 0; break; }
		  // and not anything else
		  for (int a=0; a<na; a++)
		    {
		      if ( a != a1 )
			{
			  if ( XA[a][e-f] || XA[a][e+f] )
			    {
			      okay = 0;
			      break;
			    }
			}
		    }
		}	  
	    }
	}
      
      // still okay?
      if ( okay == 1 ) 
	X[e] = true;
      
    }
    
  //
  // now update masks given X[]
  //
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  int cnt_any_match = 0;
  for (int e=0;e<ne;e++)
    {        
      const bool set_mask = ! X[e] ;
      if ( ! set_mask ) ++cnt_any_match;
      int mc = set_epoch_mask( e , set_mask );
      if ( mc == +1 ) ++cnt_mask_set;
      else if ( mc == -1 ) ++cnt_mask_unset;
      else ++cnt_unchanged;      
      if ( ! mask[e] ) ++cnt_now_unmasked;      
    }
  
  logger << "  " << cnt_any_match << " epochs match, based on stable/"
	 << ( unique ? "unique" : "any" ) << " annotations (" << x
	 << " flanking), seen " << na_obs << " of " << anns.size() << "\n";
  logger << "  " << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";
  
}


void timeline_t::clear_epoch_mask( bool b )
{
  mask.clear();
  mask_set = b;  // i.e. if b==T, equivalent to masking all entries
  mask.resize( epochs.size() , b );
  if ( epoched() )
    logger << "  reset all " << epochs.size() << " epochs to be " << ( b ? "masked" : "included" ) << "\n";
}

int timeline_t::set_epoch_mask( const int e , const bool b ) 
{
  //  std::cout <<" set_epoch_mask = " << e << " " << b << "\n";
  
  mask_set = true;
  
  if ( e < 0 || e >= mask.size() ) Helper::halt( "internal error setting mask" );
  
  bool original = mask[e];
  
  // implement mask mode
  // only mask
  if      ( mask_mode == 0 ) 
    {
      if ( (!original) && b ) mask[e] = true;  // default
    }
  else if ( mask_mode == 1 ) // 'unmask' --> only unmask
    {
      if ( original && (!b) ) mask[e] = false;
    }
  else if ( mask_mode == 2 ) // 'force' --> set either way
    {
      mask[e] = b ; // force (default)   
    }
  
  // teturn 0 if no change;
  // return +1 if set a mask (N->Y)
  // return -1 if freed a mask (Y->N)
  
  if ( original == mask[e] ) return 0;
  return mask[e] ? 1 : -1 ;     
}



bool timeline_t::is_epoch_mask_set() const
{
  return mask_set;
}

void timeline_t::set_epoch_mask_mode( const int m ) 
  {
    mask_mode = m;
  }

int timeline_t::epoch_mask_mode() const
{
  return mask_mode;
} 

bool timeline_t::masked( const int e ) const
{ 
  return mask[e]; 
}


void timeline_t::add_mask_annot( const std::string & tag )
{

  if ( ! epoched() ) return;
  
  first_epoch();

  logger << "  adding annotation " << tag << " to mark unmasked (included) epochs\n";

  annot_t * a = annotations->add( tag );

  a->description = "Included (unmasked) epoch";

  while ( 1 )
    {
      
      int e = next_epoch();

      if ( e == -1 ) break;
      
      instance_t * instance = a->add( "." , epoch(e) , "." );
      
    }
  
}

void timeline_t::dumpmask( const param_t & param )
{

  // also dump an 
  const bool dump_annot = param.has( "annot" );

  const std::string annot_str = dump_annot ? param.value( "annot" ) : "" ; 

  // default is to make annot when an epoch is /masked/ (versus opposite)
  const bool annot_unmasked = param.yesno( "annot-unmasked" );
  
  annot_t * ann = dump_annot ? annotations->add( annot_str ) : NULL ; 
  
  // output?
  const bool output = param.has( "output" ) ? param.yesno( "output" ) : true ; 
  
  // no mask set: means all clear so display that
  // if ( ! mask_set ) return;

  first_epoch();

  if ( output ) 
    logger << "  dumping MASK\n";
  if ( dump_annot )
    logger << "  creating annotation " << annot_str << " to indicate "
	   << ( annot_unmasked ? "unmasked" : "masked" ) << " epochs\n";
  
  
  while ( 1 ) 
    {
      
      int e = next_epoch_ignoring_mask();      
      
      if ( e == -1 ) break;
      
      interval_t interval = epoch( e );
      
      // EPOCH_INTERVAL will already have been output by the EPOCH command
      writer.epoch( display_epoch( e ) );

      const bool m =  mask_set ? mask[e] : false ;

      if ( output ) 
	writer.value( "EMASK" , m );
     
      if ( ann )
	{
	  if ( annot_unmasked && ! m ) 
	    ann->add( "." , interval , "." );
	  else if ( m && ! annot_unmasked )
	    ann->add( "." , interval , "." );
	}
    }

  writer.unepoch();
  
}

// mask as simply 0/1 in a file, which must match # of epochs
void timeline_t::load_mask( const std::string & f , bool exclude )
{
  
  if ( ! epoched() ) 
    {
      int ne = set_epoch( globals::default_epoch_len , globals::default_epoch_len );
      logger << "  set epochs to default " << globals::default_epoch_len << " seconds, " << ne << " epochs\n";
    }
     
  if ( ! Helper::fileExists( f ) ) Helper::halt( "could not find " + f );
  
  logger << "  attaching mask file " << f << "\n";
  
  logger << "  currently, mask mode set to: ";
  int mm = epoch_mask_mode();
  if ( mm == 0 ) logger << " mask (default)\n";
  else if ( mm == 1 ) logger << " unmask\n";
  else if ( mm == 2 ) logger << " force\n";

  
  // load
  std::ifstream FIN( f.c_str() , std::ios::in );

  int cnt_total = num_total_epochs();
  int cnt_mask0 = 0;
  int cnt_mask1 = 0;
  
  int e = 0;
  
  while ( ! FIN.eof() )
    {
      
      int m = 0;
      
      FIN >> m;

      if ( FIN.eof() ) continue;
      
      if ( ( exclude && m == 1 ) || ( (!exclude) && m == 0 ) )
	{
	  if ( ! masked(e) ) ++cnt_mask1;
	  set_epoch_mask( e );
	  ++cnt_mask0;	  
	}
      
      ++e;

      if ( e > cnt_total )
	{
	  logger << e << " masks read, for " << cnt_total << " existing epochs\n";
	  Helper::halt( "too many epochs specified in " + f );	
	}
    }
  
  
  FIN.close();

  logger << "  processed " << e
	    << " lines, with "
	    << cnt_mask0 << " masked epochs\n";

  logger << "  changed mask for " << cnt_mask1
	    << " of " << cnt_total << " epochs\n";

  return;

}



// interval list corresponding to original EDF times
void timeline_t::load_interval_list_mask( const std::string & f , bool exclude )
{

  
  Helper::halt( "not supported" );

  // assume format  time1   time2    { meta-data ....  ignored }

  // if +time1 +time2  implies an offset from start of record
  // otherwise, assume implies real clocktime, based on edf.header.starttime
  // need to handle discontinuous EDFs here
  
  
  if ( ! Helper::fileExists( f ) ) Helper::halt( "could not find " + f );
  
  logger << "  reading intervals to " << ( exclude ? " exclude" : "retain" ) << " from " << f << "\n";
  
  logger << "  currently, mask mode set to: ";
  int mm = epoch_mask_mode();
  if      ( mm == 0 ) logger << " mask (default)\n";
  else if ( mm == 1 ) logger << " unmask\n";
  else if ( mm == 2 ) logger << " force\n";
  
  // load
  std::ifstream FIN( f.c_str() , std::ios::in );
  
  std::vector<interval_t> intervals;
  int cnt = 0;
  while ( ! FIN.eof() )
    {

      std::string line;

      Helper::safe_getline( FIN , line );
      
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      
      if ( tok.size() == 0 ) 
	continue;
      
      if ( tok.size() < 2 ) 
	Helper::halt( "bad format in " + f + ", expecting at least 2 tab-delimited time fields" ) ;
      
      clocktime_t t1( tok[0] );
      clocktime_t t2( tok[1] );
      
      if ( ! t1.valid ) Helper::halt( "invalid HH:MM:SS timestring: " + tok[0] );
      if ( ! t2.valid ) Helper::halt( "invalid HH:MM:SS timestring: " + tok[1] );
      ++cnt;
    }
  
  FIN.close();

  logger << "  processed " << cnt << " " << intervals.size() << " intervals\n";


  //
  // figure out start time of EDF... either from header, or from EDF itself, i.e. if it has been editted. 
  //

  
  //
  // Make sure that we have a time-track set
  //
  
  edf->add_time_track();
  
    
  return;
    
}


// directly specify epoch-level masks
// void timeline_t::apply_simple_epoch_mask( const std::set<std::string> & labels , const std::string & onelabel , bool include )
// {
  
//   // if 'ifnot', can only specify a single 
//   if ( labels.size() > 1 && ! include ) Helper::halt( "can only specify a single mask for 'ifnot'");

//   mask_set = true;
  
//   const int ne = epochs.size();
 
//   // Note: we do not clear the mask here, as we want to allow multiple
//   // filters to be added on top of oneanther
  
//   int cnt_mask_set = 0;
//   int cnt_mask_unset = 0;
//   int cnt_unchanged = 0;
//   int cnt_now_unmasked = 0;
//   int cnt_basic_match = 0;  // basic count of matches, whether changes mask or not
  

//   for (int e=0;e<ne;e++)
//     {

//       bool matches = false;
      
//       std::set<std::string>::const_iterator ii = labels.begin();
//       while ( ii != labels.end() )
// 	{
// 	  if ( epoch_annotation( *ii , e ) ) { matches = true; break; }
// 	  ++ii;
// 	}
      
//       // count basic matches
//       if ( matches ) ++cnt_basic_match;
      
//       // set new potential mask, depending on match_mode      
//       bool new_mask = mask[e];
      
//       if ( include ) 
// 	{
// 	  if      ( mask_mode == 0 ) new_mask = matches;   // mask-if
// 	  else if ( mask_mode == 1 ) new_mask = !matches;  // unmask-if
// 	  else if ( mask_mode == 2 ) new_mask = matches ;  // if
// 	}
//       else
// 	{
// 	  if      ( mask_mode == 0 ) new_mask = !matches;  // mask-ifnot
// 	  else if ( mask_mode == 1 ) new_mask = matches;   // unmask-ifnot
// 	  else if ( mask_mode == 2 ) new_mask = ! matches; // ifnot
// 	}

//       int mc = set_epoch_mask( e , new_mask );

//       if      ( mc == +1 ) ++cnt_mask_set;
//       else if ( mc == -1 ) ++cnt_mask_unset;
//       else                 ++cnt_unchanged;
      
//       if ( !mask[e] ) ++cnt_now_unmasked;
      
//     }
  
//   logger << "  based on " << onelabel << " " << cnt_basic_match << " epochs match; ";

//   logger << cnt_mask_set << " newly masked, "   
// 	 << cnt_mask_unset << " unmasked, " 
// 	 << cnt_unchanged << " unchanged\n";
//   logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

//   // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
//   writer.level( onelabel , "EMASK" );

//   writer.var( "N_MATCHES"    , "Number of matching epochs" );
//   writer.var( "N_MASK_SET"   , "Number of epochs newly masked" ); 
//   writer.var( "N_MASK_UNSET" , "Number of epochs newly unmasked" );
//   writer.var( "N_UNCHANGED"  , "Number of epochs unchanged by this mask" );
//   writer.var( "N_RETAINED"   , "Number of epochs retained for analysis" );
//   writer.var( "N_TOTAL"      , "Total number of epochs" );

//   writer.value( "N_MATCHES"    , cnt_basic_match  );
//   writer.value( "N_MASK_SET"   , cnt_mask_set     );
//   writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
//   writer.value( "N_UNCHANGED"  , cnt_unchanged    );
//   writer.value( "N_RETAINED"   , cnt_now_unmasked );
//   writer.value( "N_TOTAL"      , (int)epochs.size()    );

//   writer.unlevel( "EMASK" );

// }



bool timeline_t::elapsed_seconds_to_spanning_epochs( const double t1, const double t2a, int * e1 , int * e2 )
{
  
  // reduce t2 by 1 tp, i.e. so we get time up to but not including t2
  // thus 0..30  means just 0 <= x < 30 , i.e. the first epoch
  
  const double t2 = t2a - globals::tp_duration;

  if ( t1 < 0 || t2 < 0 ) return false;
  
  if ( standard_epochs )
    {
      // return 1-based values
      *e1 = 1 + floor( t1 / epoch_length() );
      *e2 = 1 + floor( t2 / epoch_length() );
      if ( *e1 > *e2 ) return false;
    }

  // nb. use t2a here, i.e. tp2 is just after end 
  uint64_t tp1 = t1 * globals::tp_1sec;
  uint64_t tp2 = t2a * globals::tp_1sec;
  
  // need to count manually
  if ( ! standard_epochs )
    {
      *e1 = *e2 = -1;
      
      for (int e=0; e<epochs.size(); e++)
	{
	  //       T1----------------T2
	  // E1---
	  //      E2---
	  //                          
	  if ( *e1 == -1 && epochs[e].stop > tp1 ) *e1 = 1 + e ;
	  if ( epochs[e].start < tp2 ) *e2 = 1 + e;
	}
      
      if ( *e1 > *e2 ) return false;
      if ( *e1 == -1 || *e2 == -1 ) return false;
    }
  
  return true;
  
}
