
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
  //  unmask-ifnot                        // if doesn't have X, then unmask
  //  ifnot                      force 
   
  int mask_mode = -1;
  int match_mode = -1; 
  
  std::string condition;
  
  if ( param.has( "mask-if" ) )
    {
      condition = param.value( "mask-if" );
      mask_mode = 0;
      match_mode = 1;
    }
  else if ( param.has( "unmask-if" ) )
    {
      condition = param.value( "unmask-if" );
      mask_mode = 1;
      match_mode = 1;
    }
  else if ( param.has( "if" ) )
    {
      condition = param.value( "if" );
      mask_mode = 2;
      match_mode = 1;
    }
  else if ( param.has( "mask-ifnot" ) )
    {
      condition = param.value( "mask-ifnot" );
      mask_mode = 0;
      match_mode = 0;
    }
  else if ( param.has( "unmask-ifnot" ) )
    {
      condition = param.value( "unmask-ifnot" );
      mask_mode = 1;
      match_mode = 0;
    }
  else if ( param.has( "ifnot" ) )
    {
      condition = param.value( "ifnot" );
      mask_mode = 2;
      match_mode = 0;
    }

  
  if ( mask_mode > -1 ) 
    {
      edf.timeline.set_epoch_mask_mode( mask_mode );  
      logger << "  set masking mode to " << ( mask_mode == 2 ? "'force'" : mask_mode == 1 ? "'unmask'" : "'mask' (default)" ) << "\n";
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
      t2 -= globals::tp_duration;
      int epoch1 = 1 + floor( t1 / edf.timeline.epoch_length() );
      int epoch2 = 1 + floor( t2 / edf.timeline.epoch_length() );
      edf.timeline.select_epoch_range( epoch1 , epoch2 , true );
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
      s2 -= globals::tp_duration;
      
      // for the continuous case:
      if ( edf.header.continuous )
	{
	  int epoch1 = 1 + floor( s1 / edf.timeline.epoch_length() );
	  int epoch2 = 1 + floor( s2 / edf.timeline.epoch_length() );	  
	  if ( epoch1 > epoch2 ) Helper::halt( "misspecified hms times: implies start after end" );	  
	  edf.timeline.select_epoch_range( epoch1 , epoch2 , true );
	}

      // for the discontinuous case, select all epochs that are spanned by this interval
      if ( ! edf.header.continuous )
	{
	  int epoch1 = -1 , epoch2 = -1;

	  // look at each epoch
	  edf.timeline.first_epoch();	  
	  while ( 1 ) 
	    {
	      int epoch = edf.timeline.next_epoch();      
	      if ( epoch == -1 ) break;
	      interval_t interval = edf.timeline.epoch( epoch );

	      // std::cout << "epoch "
	      // 		<< interval.start << "\t"
	      // 		<< interval.stop << "\t"
	      // 		<< interval.start_sec() << "\t"
	      // 		<< interval.stop_sec_exact() << "\t"
	      // 		<< s1 << "\t"
	      // 		<< s2 << "\t"
	      // 		<< epoch1 << "\t"
	      // 		<< epoch2 << "\n";

	      // is this epoch spanned by this hms range? track first / last spanned
	      if ( interval.start_sec() <= s2 && interval.stop_sec_exact() >= s1 ) 
		{
		  std::cout << "  ADDING.....\n";
		  if ( epoch1 != -1 )
		    epoch2 = epoch;
		  else
		    epoch1 = epoch;
		}
	    }

	  // if only a single epoch flagged
	  if ( epoch2 == -1 ) epoch2 = epoch1;
	  
	  // nb. 1-based selection here
	  if ( epoch1 != -1 )
	    edf.timeline.select_epoch_range( 1+epoch1 , 1+epoch2 , true );

	}
      
    }
  




  //
  // Include/exclude/annotate masks
  //
  
  // nb. these now mutually exclusive
  bool has_imask = match_mode == 1 ;
  bool has_xmask = match_mode == 0 ;
 
  std::string imask_str = has_imask ? condition : "";
  std::string xmask_str = has_xmask ? condition : "" ;

  // add epoch-level 'annotations' instead of masking
  bool has_amask = param.has( "flag" );
  bool has_alabels = param.has( "label" );
  if ( has_amask != has_alabels ) 
    Helper::halt( "need to specify both flag and labels together" );
  std::string amask_str = has_amask ? param.value( "flag" ) : "" ;
  std::string alabel_str = has_alabels ? param.value( "label" ) : "" ;
  
  
  //
  // MASK include [ if ]
  //
  
  if ( has_imask )
    {
      std::vector<std::string> im0 = Helper::parse( imask_str , "," );      

      // if MASK if=A,B does not work, i.e. in force mode, this is the 
      // same as saying MASK if=B
 
      if ( im0.size() > 1 && mask_mode == 2 ) 
	Helper::halt( "cannot specify multiple annotations with an 'if' mask" );
      
      for (int i=0;i<im0.size();i++)
	{
	  // is a value specified?
	  bool has_values = false;
	  
	  std::vector<std::string> im = Helper::parse( im0[i] , "[]" );
	  if      ( im.size() == 1 ) has_values = false;
	  else if ( im.size() != 2 ) Helper::halt( "incorrectly specified annot[value(s)] -- expecting ann1, ann1[val1] or ann1[val1|val2]" );
	  else if ( im0[i][im0[i].size()-1] != ']' ) Helper::halt( "incorrectly specified annot[value(s)] -- expecting ann1, ann1[val1] or ann1[val1|val2]" );
	  else    has_values = true;
	  
	  const std::string annot_label = Helper::unquote( im[0] );
	  
	  annot_t * annot = edf.timeline.annotations( annot_label );
	  
	  if ( annot == NULL ) 
	    { 
	      // does not have mask -- 
	      edf.timeline.apply_empty_epoch_mask( annot_label , true );
	      continue; // do nothing
	    }	  
	  
	  // do we have values? 
	  if ( has_values )
	    {
	      std::vector<std::string> imask_val = Helper::parse( im[1] , "|" );
	      std::set<std::string> ss;
	      for (int v=0;v<imask_val.size();v++) ss.insert( imask_val[v] );      
	      edf.timeline.apply_epoch_include_mask( annot , &ss );
	    }
	  else
	    {	      
	      edf.timeline.apply_epoch_include_mask( annot );
	    }
	  
	}
    }
  
  

  //
  // MASK exclude [ ifnot ]
  //
  
  if ( has_xmask )
    {
      std::vector<std::string> xm0 = Helper::parse( xmask_str , "," );      
      if ( xm0.size() > 1 ) Helper::halt( "cannot specify multiple annotations with an 'ifnot' mask" );
      for (int i=0;i<xm0.size();i++)
	{
	  bool has_values = false;
	  std::vector<std::string> xm = Helper::parse( xm0[i] , "[]" );
	  if      ( xm.size() == 1 ) has_values = false;
	  else if ( xm.size() != 2 ) Helper::halt( "incorrectly specified annot[value(s)] -- expecting ann1, ann1[val1] or ann1[val1|val2]" );
	  else if ( xm0[i][xm0[i].size()-1] != ']' ) Helper::halt( "incorrectly specified annot[value(s)] -- expecting ann1, ann1[val1] or ann1[val1|val2]" );
	  else    has_values = true;
	  
	  const std::string annot_label = Helper::unquote( xm[0] );
	  
	  annot_t * annot = edf.timeline.annotations( annot_label );

	  if ( annot == NULL ) 
	    {
	      edf.timeline.apply_empty_epoch_mask( annot_label , false );
	      continue; 
	    }

	  // do we have values? 
	  if ( has_values )
	    {
	      std::vector<std::string> xmask_val = Helper::parse( xm[1] , "|" );
	      std::set<std::string> ss;
	      for (int v=0;v<xmask_val.size();v++) ss.insert( xmask_val[v] );      
	      edf.timeline.apply_epoch_exclude_mask( annot , &ss );
	    }
	  else
	    {	      
	      edf.timeline.apply_epoch_exclude_mask( annot );
	    }

	}
    }
  
}
