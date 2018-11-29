
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
#include "edf.h"
#include "main.h"

//
// MASK : apply a mask to zero-out parts of the signal
//

void proc_mask( edf_t & edf , param_t & param )
{

//   std::cerr << " currently, mask mode set to: ";
//   int mm = edf.timeline.epoch_mask_mode();
//   if ( mm == 0 ) std::cerr << " mask (default)\n";
//   else if ( mm == 1 ) std::cerr << " unmask\n";
//   else if ( mm == 2 ) std::cerr << " force\n";
  
  if ( ! edf.timeline.epoched() ) 
    Helper::halt( "no EPOCHs set, cannot apply an epoch mask" );
  
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
      std::cerr << " set masking mode to " << ( mask_mode == 2 ? "'force'" : mask_mode == 1 ? "'unmask'" : "'mask' (default)" ) << "\n";
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

  if ( param.has( "first" ) )
    {
      int n = param.requires_int( "first" );
      if ( n < 1 ) Helper::halt( "first value must be >= 1" );
      edf.timeline.select_epoch_first( n );
    }

  if ( param.has( "epoch" ) )
    {
      std::vector<int> val = param.intvector( "epoch" , "-" );
      if ( val.size() == 1 ) 
	{
	  if ( val[0] < 1 ) Helper::halt( "epoch value must be >= 1" );
	  edf.timeline.select_epoch_range( val[0] , val[0] );
      	}
      else 
	{
	  if ( val.size() != 2 ) Helper::halt("epoch=a-b");
	  if ( val[0] > val[1] ) Helper::halt("epoch=a-b requires a <= b");
	  edf.timeline.select_epoch_range( val[0] , val[1] );
	}
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
      edf.timeline.select_epoch_range( epoch1 , epoch2 );      
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
      
      double h1 = clocktime_t::difference( starttime , t1 );
      double h2 = clocktime_t::difference( starttime , t2 );
      
      h2 -= globals::tp_duration;
      
      int epoch1 = 1 + floor( ( h1 * 3600 ) / edf.timeline.epoch_length() );
      int epoch2 = 1 + floor( ( h2 * 3600 ) / edf.timeline.epoch_length() );
      
      if ( epoch1 > epoch2 ) Helper::halt( "misspecified hms times: implies start after end" );

      edf.timeline.select_epoch_range( epoch1 , epoch2 );      
    }


  //
  // Include/exclude/annotate masks
  //
  
  // nb. these now mutually exclusive
  bool has_imask = match_mode == 1 ;
  bool has_xmask = match_mode == 0 ;
  
  // add epoch-level 'annotations' instead of masking
  bool has_amask = param.has( "flag" );
  bool has_alabels = param.has( "label" );

  if ( has_amask != has_alabels ) 
    Helper::halt( "need to specify both flag and labels together" );
  
  std::string imask_str = has_imask ? condition : "";
  std::string xmask_str = has_xmask ? condition : "" ;
  std::string amask_str = has_amask ? param.value( "flag" ) : "" ;
  std::string alabel_str = has_alabels ? param.value( "label" ) : "" ;
  
  //
  // MASK include [ if ]
  //
  
  if ( has_imask )
    {
      std::vector<std::string> im0 = Helper::parse( imask_str , "," );      
      for (int i=0;i<im0.size();i++)
	{
	  std::vector<std::string> im = Helper::parse( im0[i] , "[]" );
	  if ( im.size() == 1 ) im.push_back("1");
	  if ( im.size() != 2 ) Helper::halt( "incorrectly specified include[value]" );
	  std::vector<std::string> imask_val = Helper::parse( im[1] , "," );
	  const std::string annot_label = Helper::unquote( im[0] );
	  
	  attach_annot( edf , annot_label );      
	  annot_t * annot = edf.timeline.annotations( annot_label );
	  if ( annot == NULL ) 
	    { 
	      // does not have mask -- 
	      edf.timeline.apply_empty_epoch_mask( annot_label , true );
	      continue; // do nothing
	    }

	  // TMP -- ignore values
	  // std::set<std::string> ss;
	  // for (int v=0;v<imask_val.size();v++) ss.insert( imask_val[v] );
	  // edf.timeline.apply_epoch_include_mask( annot , &ss );
	  
	  edf.timeline.apply_epoch_include_mask( annot );
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
	  std::vector<std::string> xm = Helper::parse( xm0[i] , "[]" );
	  if ( xm.size() == 1 ) xm.push_back("1");
	  if ( xm.size() != 2 ) Helper::halt( "incorrectly specified exclude[value]" );
	  std::vector<std::string> xmask_val = Helper::parse( xm[1] , "," );
	  const std::string annot_label = Helper::unquote( xm[0] );
	  
	  attach_annot( edf , annot_label );      
	  annot_t * annot = edf.timeline.annotations( annot_label );

	  if ( annot == NULL ) 
	    {
	      edf.timeline.apply_empty_epoch_mask( annot_label , false );
	      continue; 
	    }
	  // IGNORE values
	  // std::set<std::string> ss;
	  // for (int v=0;v<xmask_val.size();v++) ss.insert( xmask_val[v] );      
	  // edf.timeline.apply_epoch_exclude_mask( annot , &ss );

	  edf.timeline.apply_epoch_exclude_mask( annot );
	}
    }
  

  //
  // MASK annot
  //
  
  if ( has_amask )
    {

      std::vector<std::string> am0 = Helper::parse( amask_str , "," );      
      for (int i=0;i<am0.size();i++)
	{
	  
	  std::vector<std::string> am = Helper::parse( am0[i] , "[]" );
	  if ( am.size() == 1 ) am.push_back("1");
	  if ( am.size() != 2 ) Helper::halt( "incorrectly specified flag[value]" );
	  
	  std::vector<std::string> amask_val = Helper::parse( am[1] , "," );
	  const std::string annot_label = Helper::unquote( am[0] );
	  
	  attach_annot( edf , annot_label );      
	  
	  annot_t * annot = edf.timeline.annotations( annot_label );
      
	  if ( annot == NULL ) 
	    edf.timeline.annotate_epochs( alabel_str , false );
	  else
	    {
	      std::set<std::string> ss;
	      for (int v=0;v<amask_val.size();v++) ss.insert( amask_val[v] );
	      edf.timeline.annotate_epochs( alabel_str , annot_label , ss );
	    }
	  
	  std::cerr << " set flag annotation for [" << annot_label << "], labeled as [" << alabel_str << "]\n";
	}      
            
    }

}
