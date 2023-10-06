
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
#include "annot/annotate.h"

extern logger_t logger;

//
// MASK : apply a mask to zero-out parts of the signal
//

void proc_mask( edf_t & edf , param_t & param )
{

  //
  // How keep things safe, we only allow a single parmeter for a MASK command
  //
  
  // single() takes into account 'hidden' params (i.e. signal)
  
  if ( ! param.single() ) 
    Helper::halt( "MASK commands can only take a single parameter" );
  
    
  if ( ! edf.timeline.epoched() ) 
    {
      int ne = edf.timeline.set_epoch( globals::default_epoch_len , globals::default_epoch_len );
      logger << "  set epochs, to default length " << globals::default_epoch_len << ", " << ne << " epochs\n";
    }


  //
  // Priamry annotation-based include/exclude masks
  //

  //
  // Mode:   mask[default]     // Existing      Evaluated      mask    unmask    both
  //                              N             N              N       N         N
  //                              N             Y             [Y]      N        [Y]

  //                              Y             N              Y      [N]       [N]   
  //                              Y             Y              Y       Y         Y       
  //
  //   mask 0 mask (default)
  //        1 unmask 
  //        2 force (i.e. both mask and unmask)
  //


  //  mask-if          include   mask
  //  unmask-if        include   unmask
  //  if               include   force 
  
  //  mask-ifnot       exclude   mask     // if doesn't have X, then mask
  //  unmask-ifnot     exclude   unmask   // if doesn't have X, then unmask
  //  ifnot                      force 

  //  mask-if-all       include   mask
  //  unmask-if-all     include   unmask
  //  if-all            include   force 
  
  //  mask-ifnot       exclude   mask     // if doesn't have X, then mask
  //  unmask-ifnot     exclude   unmask   // if doesn't have X, then unmask
  //  ifnot                      force 

  
  int mask_mode = -1;
  int match_mode = -1; 
  bool match_logic_or = true; // defaults to OR logic for multiple annotations
  
  std::string condition;
  
  if      ( param.has( "mask-if" ) )     condition = param.value( "mask-if" );
  else if ( param.has( "mask-if-any" ) ) condition = param.value( "mask-if-any" );  // mask-if == mask-if-any
  else if ( param.has( "mask-if-all" ) ) condition = param.value( "mask-if-all" );

  else if ( param.has( "unmask-if" ) )     condition = param.value( "unmask-if" );
  else if ( param.has( "unmask-if-any" ) ) condition = param.value( "unmask-if-any" );
  else if ( param.has( "unmask-if-all" ) ) condition = param.value( "unmask-if-all" );

  else if ( param.has( "if" ) )     condition = param.value( "if" );
  else if ( param.has( "if-any" ) ) condition = param.value( "if-any" );
  else if ( param.has( "if-all" ) ) condition = param.value( "if-all" );
  
  else if ( param.has( "mask-ifnot" ) )     condition = param.value( "mask-ifnot" );
  else if ( param.has( "mask-ifnot-any" ) ) condition = param.value( "mask-ifnot-any" );
  else if ( param.has( "mask-ifnot-all" ) ) condition = param.value( "mask-ifnot-all" );

  else if ( param.has( "unmask-ifnot" ) )     condition = param.value( "unmask-ifnot" );
  else if ( param.has( "unmask-ifnot-any" ) ) condition = param.value( "unmask-ifnot-any" );
  else if ( param.has( "unmask-ifnot-all" ) ) condition = param.value( "unmask-ifnot-all" );

  else if ( param.has( "ifnot" ) )     condition = param.value( "ifnot" );
  else if ( param.has( "ifnot-any" ) ) condition = param.value( "ifnot-any" );
  else if ( param.has( "ifnot-all" ) ) condition = param.value( "ifnot-all" );
  
  // 1) if vs ifnot match mode?

  match_mode = 0;  // all 'ifnot'
  if      ( param.has( "if" )     || param.has( "mask-if" )     || param.has( "unmask-if" ) )     match_mode = 1 ;
  else if ( param.has( "if-any" ) || param.has( "mask-if-any" ) || param.has( "unmask-if-any" ) ) match_mode = 1 ;
  else if ( param.has( "if-all" ) || param.has( "mask-if-all" ) || param.has( "unmask-if-all" ) ) match_mode = 1 ;


  // 2) mask mode: force / unmask / mask ?
  mask_mode = 0; // 'mask' by default
  // unmask?
  if ( param.has( "unmask-if" )    || param.has( "unmask-if-any" )    || param.has( "unmask-if-all" ) ) mask_mode = 1;
  if ( param.has( "unmask-ifnot" ) || param.has( "unmask-ifnot-any" ) || param.has( "unmask-ifnot-all" ) ) mask_mode = 1;    
  // force?
  if ( param.has( "if" )    || param.has( "if-any" )    || param.has( "if-all" ) )   mask_mode = 2;
  if ( param.has( "ifnot" ) || param.has( "ifnot-any" ) || param.has( "ifnot-all" ) ) mask_mode = 2;


      
  if ( mask_mode > -1 ) 
    {
      edf.timeline.set_epoch_mask_mode( mask_mode );  
      logger << "  set masking mode to " << ( mask_mode == 2 ? "'force'" : mask_mode == 1 ? "'unmask'" : "'mask' (default)" ) << "\n";
    }

  // 3) multiple annot logic mode
  match_logic_or = true;
  // and? (i.e. 'all' )
  if ( param.has( "if-all" )    || param.has( "mask-if-all" )    || param.has( "unmask-if-all" ) ) match_logic_or = false;
  if ( param.has( "ifnot-all" ) || param.has( "mask-ifnot-all" ) || param.has( "unmask-ifnot-all" ) ) match_logic_or = false;


  // apply an annotation-based mask? 
  if ( condition != "" )
    {      
      
      // expand any wildcards in annotation list
      //   note: this will not work for instance-ID options: class IDs only
      
      // if condition is comma-delimited, expand out any root* matches;
      //  otherwise all 
      std::set<std::string> conditions = annotate_t::root_match( condition , edf.timeline.annotations.names() );

      // build primary input for epoch mask:
      std::map<annot_t*,std::set<std::string> > amask;

      // for console output
      std::map<std::string,std::set<std::string> > amask_labels;

      // require full span for this annot?
      std::set<std::string> fullspan;

      // iterate over annot classes
      std::set<std::string>::const_iterator aa = conditions.begin();
      while ( aa != conditions.end() )
	{
	  // are instance values specified?   annot[V1|V2]
	  // is a value specified?

	  const std::vector<std::string> tok = Helper::parse( *aa , "[]" );

	  if ( tok.size() == 0 )
	    {
	      ++aa;
	      continue;
	    }
	  
	  std::string annot_label = Helper::unquote( tok[0] );

	  const std::string annot_label_orig = annot_label ; 
	  
	  // special syntax to require complete overlap
	  //  (nb. this means will not be able to combine +annot* )	  
	  
	  if ( annot_label[0] == '+' )
	    {
	      // drop special leading char
	      annot_label = annot_label.substr( 1 );
	      // send to timeline when doing matching
	      fullspan.insert( annot_label );	      
	    }
	       	  
          annot_t * annot = edf.timeline.annotations( annot_label );
	  
	  std::set<std::string> values;
	  
	  // not found?  need to add as an 'empty' annotation (but we can
	  // ignore any values) - the match function will always return a
	  // non-match

	  if ( annot == NULL )
	    {
	      amask[ NULL ] = values ;  // an empty record
	      amask_labels[ annot_label_orig ] = values; // for console output below
	      ++aa;
	      continue;
	    }

	  // parse for any explicitly specified instance ID matches?
          const bool has_values = tok.size() != 1 ;
	  
	  if ( has_values )
	    {
	      if ( tok.size() != 2 )
		Helper::halt( "incorrectly specified annot[value(s)] -- expecting ann1, ann1[val1] or ann1[val1|val2]" );
	      if ( (*aa)[ aa->size() - 1 ] != ']' )
		Helper::halt( "incorrectly specified annot[value(s)] -- expecting ann1, ann1[val1] or ann1[val1|val2]" );
	      
	      // multiple pipe-delimited instance IDs?  annot[v1|v2]
	      const std::vector<std::string> tok2 = Helper::parse( tok[1] , "|" );              
              for (int v=0;v<tok2.size();v++)
		if ( tok2[v] != "" )
		  values.insert( tok2[v] );
	      
	    }

	  // add to map
	  amask[ annot ] = values;
	  amask_labels[ annot_label_orig ] = values; // for console output below                                                                       
	  // next annotation
	  ++aa;
	}

      //
      // some console blurb 
      //

      logger  << "  annots:";
      std::string alabel;
      std::map<std::string,std::set<std::string> >::const_iterator ll = amask_labels.begin();
      while ( ll != amask_labels.end() )
	{
	  if ( alabel != "" ) alabel += ",";
	  logger << " " << ll->first;
	  alabel += ll->first;
	  if ( ll->second.size() != 0 )
	    {
	      const std::string s1 = Helper::stringize( ll->second , "|" );
	      logger << "[" << s1 << "]";
	      alabel += "[" + s1 + "]";
	    }
	  ++ll;
	}
      
      //
      // apply the actual mask
      //
      
      edf.timeline.apply_epoch_mask2( amask , fullspan , alabel , match_logic_or , match_mode );

      
      // all done
      return;
    }
  
  //
  // Wipe entire mask, i.e. include all 
  //
  
  if ( param.has( "clear" ) || param.has( "include-all" ) || param.has( "none" ) )
    {
      edf.timeline.clear_epoch_mask( false ); // i.e. set all mask to F
      return;
    }
  
  //
  // Excldue all
  //
  
  if ( param.has("all") || param.has( "exclude-all" ) || param.has( "total" ) )
    {
      edf.timeline.clear_epoch_mask( true ); // i.e. set all mask to T
      return;
    }
  

  //
  // Eval expression masks
  //
  
  bool verbose = false;

  if ( param.has( "expr" ) ) 
    {
      // force, mask_mode == 2 
      edf.timeline.apply_eval_mask( param.value( "expr" ) ,  2 , verbose ) ; 
      return;
   }

  if ( param.has( "not-expr" ) ) 
    {
      // force, mask_mode == -2 
      edf.timeline.apply_eval_mask( param.value( "not-expr" ) ,  -2 , verbose ) ; 
      return;
   }

  if ( param.has( "mask-expr" ) )
    {
      // mask, mask_mode == 0 
      edf.timeline.apply_eval_mask( param.value( "mask-expr" ) , 0 , verbose ) ; 
      return;      
    }
  
  if ( param.has( "unmask-expr" ) )
    {
      // unmask , mask_mode == 1 
      // i.e. if true, unmask
      edf.timeline.apply_eval_mask( param.value( "unmask-expr" ) , 1 , verbose ) ; 
      return;            
    }
  

  //
  // 'Special' masks
  //
  
  if ( param.has( "random" ) )
    {
      int n = param.requires_int( "random" );
      if ( n < 1 ) Helper::halt( "random value must be >= 1" );
      edf.timeline.select_epoch_randomly(n);
    }

  if ( param.has( "flip" ) ) 
    {
      edf.timeline.flip_epoch_mask();
    }

  if ( param.has( "leading" ) )
    {
      edf.timeline.select_epoch_until_isnot( param.value( "leading" ) );
    }

  // to add
  // if ( param.has( "mask-leading" ) )
  //   {
  //     edf.timeline.mask_leading_trailing( param.value( "mask-leading" ) , true , false ) ;
  //   }
  
  // if ( param.has( "mask-trailing" ) )
  //   {
  //     edf.timeline.mask_leading_trailing( param.value( "mask-leading" ) , false , true ) ;
  //   }

  if ( param.has( "regional" ) ) 
    {
      std::vector<int> r = param.intvector( "regional" );
      if ( r.size() != 2 ) Helper::halt( "expecting regional=x,y" );
      // on both sides of epoch E, requires that at least x of y flanking epochs are good
      edf.timeline.regional_mask( r[0] , r[1] );
    }

  if ( param.has( "first" ) )
    {
      int n = param.requires_int( "first" );
      if ( n < 1 ) Helper::halt( "first value must be >= 1" );
      edf.timeline.select_epoch_first( n );
    }
  
  if ( param.has( "trim" ) )
    {
      std::vector<std::string> ss = param.strvector( "trim" );
      int n = 0;
      std::string label = "";
      if ( ss.size() == 1 ) label = ss[0];
      else if ( ss.size() == 2 ) 
	{
	  label = ss[0];
	  if ( ! Helper::str2int( ss[1] , &n ) ) 
	    Helper::halt( "expecting positive integer for trim" );	  
	}
      else Helper::halt( "bad syntax for trim" );

      if ( n < 0 )  Helper::halt( "trim value must be >= 0" );
      edf.timeline.trim_epochs( label , n );
    }
  
  if ( param.has( "epoch" ) || param.has( "mask-epoch" ) )
    {
      // epoch --> 'force' mode (i.e. set all)
      // mask-epoch --> 'mask' mode
      // if 'end' means last epoch 

      bool include_mode = param.has( "epoch" );
      
      std::string label = include_mode ? "epoch" : "mask-epoch" ; 
      
      // for epoch=222-end 
      const int end_epoch = edf.timeline.num_total_epochs();

      std::vector<std::string> toks = Helper::parse( param.value( label ) , "," );
      
      std::set<int> epochs;

      for (int t=0; t<toks.size(); t++) 
	{
	  std::vector<std::string> val = Helper::parse( toks[t] , "-" );
	  if ( val.size() == 1 ) 
	    {
	      int v1=-1;
	      if ( ! Helper::str2int( val[0] , &v1 ) ) Helper::halt( label + " value must be integer" );
	      if ( v1 < 1 ) Helper::halt( label + " value must be >= 1" );
	      epochs.insert( v1 );
	    }
	  else if ( val.size() == 2 )
	    {
	      int v1=-1 , v2 = -1;
	      if ( ! Helper::str2int( val[0] , &v1 ) ) Helper::halt( label + " value must be integer" );
	      
	      if ( val[1] == "end" ) 
		v2 = end_epoch;
	      else  
		if ( ! Helper::str2int( val[1] , &v2 ) ) Helper::halt( label + " value must be integer" );
	      
	      if ( v1 > v2 ) Helper::halt( label + "=a-b requires a <= b");
	      for (int i=v1; i<= v2; i++) epochs.insert( i );
	    }
	  else
	    Helper::halt( label + "=a-b-c is bad format");
	  
	}

      edf.timeline.select_epoch_range( epochs , include_mode );

      
    }
  

  if ( param.has( "flanked" ) )
    {
      std::vector<std::string> val = param.strvector( "flanked" );
      if ( val.size() != 2 ) Helper::halt("flanked={annot},N");
      int n = 0;
      if ( ! Helper::str2int( val[1] , &n ) ) Helper::halt("flanked={annot},N");
      edf.timeline.select_epoch_within_run( val[0] , n );
    }


  if ( param.has( "sec" ) ) 
    {
      // convert to nearest epochs 
      std::vector<std::string> tok = param.strvector( "sec" , "-" );
      double t1, t2;      
      if ( tok.size() != 2 || (!Helper::str2dbl( tok[0] , &t1 )) || (!Helper::str2dbl( tok[1] , &t2 ) ) ) 
	Helper::halt( "expecting sec=<time1>-<time2> where <time> is in seconds" );

      // reduce t2 by 1 tp, i.e. so we get time up to but not including t2
      // thus 0..30  means just 0 <= x < 30 , i.e. the first epoch

      int epoch1, epoch2;
      if ( edf.timeline.elapsed_seconds_to_spanning_epochs( t1, t2, &epoch1, &epoch2 ) )
	edf.timeline.select_epoch_range( epoch1 , epoch2 , true );
      else
	logger << "  bad time ranges given from MASK sec\n";
    }


  if ( param.has( "hms" ) ) 
    {
      // convert to nearest epochs
      std::vector<std::string> tok = param.strvector( "hms" , "-" );
      if ( tok.size() != 2 ) Helper::halt( "expecting sec=<time1>-<time2> where <time> is in hh:mm:ss format" );
      
      clocktime_t starttime( edf.header.starttime );
      bool invalid_hms = ! starttime.valid;
      if ( invalid_hms ) Helper::halt( "EDF does not have valid start-time in header" );
      clocktime_t t1( tok[0] );
      clocktime_t t2( tok[1] );
      
      double s1 = clocktime_t::ordered_difference_seconds( starttime , t1 ) ; 
      double s2 = clocktime_t::ordered_difference_seconds( starttime , t2 ) ;

      int epoch1, epoch2;
      if ( edf.timeline.elapsed_seconds_to_spanning_epochs( s1, s2, &epoch1, &epoch2 ) )
        edf.timeline.select_epoch_range( epoch1 , epoch2 , true );
      else
        logger << "  bad time ranges given from MASK hms\n";
      
    }
  
  
}
