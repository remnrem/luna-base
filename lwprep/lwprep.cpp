
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

#include "lwprep.h"

#include "edf/edf.h"
#include "spindles/spectral.h"
#include "helper/helper.h"
#include "sstore/sstore.h"
#include "eval.h"
#include "db/db.h"
#include "dsp/tv.h"

#include <cstdlib>
#include <vector>
#include <map>


extern writer_t writer;

extern logger_t logger;


// step 0) for whole sample, initiate 'lw' folder, with samples and signals
// step 1) for each individual, run 'LW', mask as desired, run 'LW-MASK'


//
// code to populate the 'lw/ folder 
//   lw/samples       first two columns of the sample list
//   lw/signals       gives channel order and colors
//   lw/clocs         [optional] gives channel locations (for topo-plots)    
//   lw/inds/{indiv-id}/   

// in each individual folder, we scan for the following:
//   staging.db              epoch-level; just 'STAGE' variable
//   stage-summary.db        from 'HYPNO' ; baseline and cycle level only 
//   psd-epoch.db            assumes 0.5 hz bins (0.5 .. 20Hz)
//   psd-band-epoch.db       log(abs) band power per epoch
//   mask.db

//   Annotations
//     annot-{annot}.db

//   Spindles

//   coh.db                  slow, so base this on a random sampling of epochs


//  We assume the following variables are set:
//    eeg   which are the EEG channels


lw_prep_t::lw_prep_t( edf_t & edf , param_t & param )
{

  //
  // Root folder, for 'lw/'
  //

  std::string folder = "";

  if ( param.has( "dir" ) ) 
    {
      folder = param.value( "dir" ) ;
      if ( folder.size() > 0 ) 
	if ( folder[ folder.size() -1 ] != globals::folder_delimiter ) 
	  folder += globals::folder_delimiter;
    }


  //
  // Optional smoothing of per-epoch power
  //

  denoise = param.has( "lambda" );
  
  if ( denoise )
    lambda = param.requires_dbl( "lambda" );

  //
  // two modes -- generated all per-epoch/HYPNO measures, OR dump the mask (i.e. 'mask' option set)
  // i.e. assume we will run first to generate measures for all epochs, then run masks, then set MASK
  // and dump with LW mask
  //

  bool dump_mask = param.has( "mask" );
  
  
  //
  // the newly created luna-web folder will be 'lw' and placed (we
  // assume) in the same directory as the sample list; importantly, if
  // using docker this must also be the same folder for which relative
  // EDFs paths will work (i.e. EDFs cannot be outside the mounted
  // script)
  //

  folder += "lw"; 
  folder += globals::folder_delimiter;
  folder += "inds"; 
  folder += globals::folder_delimiter;
  folder += edf.id ;
  

  //
  // check folder exists (should have previously been created
  //
  
  std::string syscmd = "mkdir -p " + folder ;
  system( syscmd.c_str() );


  //
  // hijack the output stream, i.e. we want retval_t directly
  //
  
  writer.nodb();
  
  writer.clear();
  
  retval_t ret;
  
  writer.use_retval( &ret );
  
  writer.id( edf.id , edf.filename );
  


  //
  // mask mode?
  //

  if ( dump_mask )
    {
      
      edf.timeline.first_epoch();
        
      logger << " recording " << edf.timeline.num_total_epochs() << " epochs, of which " 
	     << edf.timeline.num_epochs() << " are unmasked, to " << folder << "/mask.db\n";
      
      // any entry means that that epoch (assuming 30-seconds, 1-based encoding) is masked
      sstore_t ss( folder + "/mask.db" );      
      
      while ( 1 ) 
	{
	  
	  int epoch = edf.timeline.next_epoch_ignoring_mask();      
	  
	  if ( epoch == -1 ) break;
	  
	  if ( edf.timeline.masked_epoch( epoch ) )
	    ss.insert_epoch( epoch + 1 , "MASK" , 1 , NULL , NULL );
	  
	}
      
      ss.index();
      
      ss.dettach();
      
      // all done for creating masks.db

      writer.use_retval( NULL );
      
      writer.clear();

      writer.nodb();

      return;
      
    }

 
  //
  // now run the following commands and save as sstore databases within 'lw/inds/' 
  // i.e. to create the set-up that luna-web will be expecting
  //
  
  cmd_t cmd( "HYPNO & ANNOTS & PSD max=20 bin=0.5 epoch-spectrum sig=${eeg}" );
    
    
  //
  // Evaluate on the current EDF
  //

  cmd.eval( edf );



  //
  // Annotations
  //
  
  std::set<std::string> annots = get_annots( ret , edf.id );

  std::set<std::string>::const_iterator aa = annots.begin();
  while ( aa != annots.end() ) 
    {

      logger << " making lw/inds/" << edf.id << "/annot-" <<  Helper::sanitize( *aa ) << ".db"; 
      
      sstore_t ss( folder + "/annot-" + Helper::sanitize( *aa ) + ".db" );
      
      insert_annot2ints( ret , edf.id , *aa , &ss );
      
      ss.index();  
      
      ss.dettach();

      ++aa;
    }


  //
  // staging.db : epoch-level stages
  //

  sstore_t staging( folder + "/staging.db" );

  insert_epoch2stage( ret , edf.id , &staging );
  
  staging.index();  

  staging.dettach();
  
  
  //
  // Stage summaries (by individual, and by cycle)
  //

  sstore_t stage_summary( folder + "/stage-summary.db" );

  insert_stage_summary( ret , edf.id , &stage_summary );
  
  stage_summary.index();  

  stage_summary.dettach();


  //
  // PSD : log(power) for 1) epoch-by-channel-by-band and 2) epoch-by-channel-by-frequency 
  //

  logger << " making " << folder << "/psd-epoch-band.db"; 
 
  sstore_t psd_band( folder + "/psd-epoch-band.db" );

  insert_psd_band( ret , edf.id , &psd_band );
  
  psd_band.index();  

  psd_band.dettach();

  logger << " making " << folder << "/psd-epoch-spec.db";

  sstore_t psd_spec( folder + "/psd-epoch-spec.db" );

  insert_psd_spec( ret , edf.id , &psd_spec );
  
  psd_spec.index();  

  psd_spec.dettach();
  


  //
  // all done, now turn off the retval (and all other DB streams)
  //
  
  writer.clear();
  
  writer.nodb();
  

}



void  lw_prep_t::insert_epoch2stage( retval_t & ret , 
				     const std::string & indiv , 
				     sstore_t * ss )
{

  retval_cmd_t    rv_cmd( "HYPNO" );  
  retval_factor_t rv_fac( "E" ); 
  retval_var_t    rv_var( "STAGE" );
  retval_indiv_t  rv_indiv( indiv );
    
  const std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > > & dat1 = ret.data[ rv_cmd ][ rv_fac ][ rv_var ];
  
  std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > >::const_iterator ii = dat1.begin();
  
  while ( ii != dat1.end() )
    {
      // get epoch
      retval_factor_level_t epoch_lvl = ii->first.find( "E" );
      if ( epoch_lvl.is_int )
	{
	  int e = epoch_lvl.int_level;
	  std::map<retval_indiv_t,retval_value_t>::const_iterator jj = ii->second.find( rv_indiv );
	  if ( jj != ii->second.end() ) 
	    {	      
	      // string value for 'STAGE'
	      ss->insert_epoch( e  , "STAGE" , jj->second.s , NULL , NULL );
	    }
	}      
      ++ii; // next level
    }
  
}


void  lw_prep_t::insert_stage_summary( retval_t & ret , 
				       const std::string & indiv , 
				       sstore_t * ss )
{

  // all baseline variables
  retval_cmd_t    rv_cmd( "HYPNO" );  
  retval_factor_t rv_baseline; 
  retval_strata_t rv_baseline_strata;
  retval_indiv_t  rv_indiv( indiv );

  // baseline
  const std::map<retval_var_t,
    std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > > > & dat1 = ret.data[ rv_cmd ][ rv_baseline ];
  
  std::map<retval_var_t,
    std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > > >::const_iterator ii = dat1.begin();

  // for each variable
  while ( ii != dat1.end() )
    {
      // get variable
      std::string var = ii->first.name;

      // only baseline strata
      std::map<retval_strata_t,std::map<retval_indiv_t,retval_value_t> >::const_iterator jj = ii->second.find( rv_baseline_strata );
      if ( jj != ii->second.end() )
	{
	  std::map<retval_indiv_t,retval_value_t>::const_iterator kk = jj->second.find( rv_indiv );
	  if ( kk != jj->second.end() )
	    {
	      // values could be int or dbl (also allow string)
	      if ( kk->second.is_int )
		ss->insert_base( var , kk->second.i , NULL , NULL );
	      else if ( kk->second.is_dbl )
		ss->insert_base( var , kk->second.d , NULL , NULL );
	      else
		ss->insert_base( var , kk->second.s , NULL , NULL );
	    }
	}
      
      ++ii; // next level
    }


  // by cycle
  retval_factor_t rv_cycle( "C" );

  const std::map<retval_var_t,
    std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > > > & datc = ret.data[ rv_cmd ][ rv_cycle ];
  
  ii = datc.begin();
  
  // for each variable
  while ( ii != datc.end() )
    {
      // get variable
      std::string var = ii->first.name;

      // for each cycle
      std::map<retval_strata_t,std::map<retval_indiv_t,retval_value_t> >::const_iterator jj = ii->second.begin();
      while ( jj != ii->second.end() )
	{

	  std::string cycle_level = jj->first.print();
	  
	  std::map<retval_indiv_t,retval_value_t>::const_iterator kk = jj->second.find( rv_indiv );
	  if ( kk != jj->second.end() )
	    {
	
	      // CH is NULL, but LVL is C=c
	      if ( kk->second.is_int )
		ss->insert_base( var , kk->second.i , NULL , &cycle_level ); 
	      else if ( kk->second.is_dbl )
		ss->insert_base( var , kk->second.d , NULL , &cycle_level ); 
	      else 
		ss->insert_base( var , kk->second.s , NULL , &cycle_level ); 
	    }
	  
	  ++jj;
	}
      
      ++ii; // next level
    }
  
}


std::set<std::string> lw_prep_t::get_annots( retval_t & ret , const std::string & indiv )
{

  retval_cmd_t    rv_cmd( "ANNOTS" );  
  retval_factor_t rv_fac( "ANNOT" ); 
  retval_var_t    rv_var( "COUNT" );
  retval_indiv_t  rv_indiv( indiv );
    
  const std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > > & dat1 = ret.data[ rv_cmd ][ rv_fac ][ rv_var ];
  
  std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > >::const_iterator ii = dat1.begin();
  
  std::set<std::string> rv;

  while ( ii != dat1.end() )
    {
      // get annot
      retval_factor_level_t annot_lvl = ii->first.find( "ANNOT" );
      if ( annot_lvl.is_str )
	rv.insert( annot_lvl.str_level );
      ++ii; 
    }

  return rv;
}


void lw_prep_t::insert_annot2ints( retval_t & ret , const std::string & indiv , const std::string & annot , sstore_t * ss )
{

  retval_cmd_t    rv_cmd( "ANNOTS" );  
  
  std::set<std::string> dummy;
  dummy.insert( "ANNOT" );
  dummy.insert( "INST" );
  dummy.insert( "T1" );  
  dummy.insert( "T2" );
  retval_factor_t rv_fac( dummy ); 
  
  retval_var_t    rv_var1( "START" );
  retval_var_t    rv_var2( "STOP" );

  retval_indiv_t  rv_indiv( indiv );
    
  const std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > > & dat1 = ret.data[ rv_cmd ][ rv_fac ][ rv_var1 ];
  
  std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > >::const_iterator ii = dat1.begin();

  int cnt = 0;

  while ( ii != dat1.end() )
    {
      // get epoch
      retval_factor_level_t annot_lvl = ii->first.find( "ANNOT" );
      retval_factor_level_t inst_lvl = ii->first.find( "INST" );
      
      if ( annot_lvl.is_str )
	{
	  if ( annot_lvl.str_level == annot ) 
	    {

 	      std::map<retval_indiv_t,retval_value_t>::const_iterator jj = ii->second.find( rv_indiv );
 	      if ( jj != ii->second.end() ) 
 		{	      
 		  // double value (SECS) for START
 		  double start = jj->second.d;
		  // assume this pair value will always exist
		  double stop = ret.data[ rv_cmd ][ rv_fac ][ rv_var2 ][ ii->first ][ rv_indiv ].d ;
		  		  
		  // insert interval (no CH, no LVL)
		  // set 'value' just to name (*aa)
		  ss->insert_interval( start , stop , annot , inst_lvl.str_level , NULL , NULL );
		  
		  ++cnt;
 		}

	    }

	}      
      ++ii; // next level
    }
  
  logger << " ... " << cnt << " intervals\n";

}



void  lw_prep_t::insert_psd_band( retval_t & ret , 
				  const std::string & indiv , 
				  sstore_t * ss )
{

  retval_cmd_t    rv_cmd( "PSD" );  

  std::set<std::string> dummy;
  dummy.insert( "E" );
  dummy.insert( "CH" );
  dummy.insert( "B" );
  retval_factor_t rv_fac( dummy ); 

  retval_var_t    rv_var( "PSD" );
  retval_indiv_t  rv_indiv( indiv );
    
  const std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > > & dat1 = ret.data[ rv_cmd ][ rv_fac ][ rv_var ];
  
  std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > >::const_iterator ii = dat1.begin();

  // accumulate a epoch-length vectors 
  // insert as CH x STAGE only

  std::map<std::string,std::map<std::string,std::vector<double> > > pows;
  
  while ( ii != dat1.end() )
    {
      // get epoch
      retval_factor_level_t epoch_lvl = ii->first.find( "E" );
      retval_factor_level_t channel_lvl = ii->first.find( "CH" );
      retval_factor_level_t band_lvl = ii->first.find( "B" );
      
      // only retain certain bands for plotting
      
      bool keeper = band_lvl.str_level == "DELTA" ||
	band_lvl.str_level == "THETA" ||
	band_lvl.str_level == "ALPHA" ||
	band_lvl.str_level == "SIGMA" ||
	band_lvl.str_level == "BETA" ;
      
      if ( keeper ) 
	{
	  if ( epoch_lvl.is_int )
	    {
	      int e = epoch_lvl.int_level;
	      std::map<retval_indiv_t,retval_value_t>::const_iterator jj = ii->second.find( rv_indiv );
	      if ( jj != ii->second.end() ) 
		{	      
		  // take log(power) from 'PSD'
		  pows[ channel_lvl.str_level ][ band_lvl.str_level ].push_back( log( jj->second.d )  );
		}
	    }      
	}

      ++ii; // next level
    }

  
  //
  // insert all (at baseline level)
  //

  int cnt = 0 ;
  std::map<std::string,std::map<std::string,std::vector<double> > >::iterator jj = pows.begin();
  while ( jj != pows.end() )
    {
      std::map<std::string,std::vector<double> > ::iterator kk = jj->second.begin();
      while ( kk != jj->second.end() )
	{
	  //std::cout << "jj->first " << jj->first << " " << kk->first << " " << kk->second.size() << "\n";
	  
	  // smooth?
	  if ( denoise ) 
	    dsptools::TV1D_denoise( kk->second , lambda );

	  ss->insert_base( "PSD" , kk->second , &(jj->first) , &(kk->first) );
	  ++cnt;
	  ++kk;
	}
      ++jj;
    }

  logger << " ... inserted " << cnt << " values\n";
}


void  lw_prep_t::insert_psd_spec( retval_t & ret , 
				  const std::string & indiv , 
				  sstore_t * ss )
{

  retval_cmd_t    rv_cmd( "PSD" );  

  std::set<std::string> dummy;
  dummy.insert( "E" );
  dummy.insert( "CH" );
  dummy.insert( "F" );
  retval_factor_t rv_fac( dummy ); 

  retval_var_t    rv_var( "PSD" );
  retval_indiv_t  rv_indiv( indiv );
    
  const std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > > & dat1 = ret.data[ rv_cmd ][ rv_fac ][ rv_var ];
  
  std::map<retval_strata_t,
    std::map<retval_indiv_t,
    retval_value_t > >::const_iterator ii = dat1.begin();

  // insert as 41-element vectors (0, 0.25, 0.5, ... , 19.75 )
  // i.e. DC and then midpoint of 0.5 Hz bins
  
  // collate here, and then insert
  std::map<int,std::map<std::string,std::vector<double> > > psds;
  
  while ( ii != dat1.end() )
    {
      // get epoch
      retval_factor_level_t epoch_lvl = ii->first.find( "E" );
      retval_factor_level_t channel_lvl = ii->first.find( "CH" );
      retval_factor_level_t freq_lvl = ii->first.find( "F" );
      std::string freq_label = Helper::dbl2str( freq_lvl.dbl_level );

      if ( epoch_lvl.is_int )
	{
	  int e = epoch_lvl.int_level;
	  std::map<retval_indiv_t,retval_value_t>::const_iterator jj = ii->second.find( rv_indiv );
	  if ( jj != ii->second.end() ) 
	    {	      
	      // take log(power) from 'PSD'	      
	      psds[ e ][ channel_lvl.str_level ].push_back( log( jj->second.d ) );
	    }
	}      
      ++ii; // next level
    }
  
  
  int cnt = 0;

  std::map<int,std::map<std::string,std::vector<double> > >::const_iterator ee = psds.begin();
  while ( ee != psds.end() )
    {
      std::map<std::string,std::vector<double> >::const_iterator jj = ee->second.begin();
      while ( jj != ee->second.end() )
	{
	  //	  std::cout << ee->first << "\t" << jj->second.size() << " " << jj->first << "\n";
	  ss->insert_epoch( ee->first  , "PSD" , jj->second , &(jj->first) , NULL );
	  ++jj;
	  ++cnt;
	}
      ++ee;            
    }

  logger << " ... inserted " << cnt << " PSDs\n";
}
