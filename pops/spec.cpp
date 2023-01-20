
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

#ifdef HAS_LGBM

#include "pops.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

std::map<std::string,pops_feature_t> pops_specs_t::lab2ftr;
std::map<pops_feature_t,std::string> pops_specs_t::ftr2lab;
std::set<std::string> pops_specs_t::lvl2;
std::map<std::string,int> pops_specs_t::blocksize;
std::vector<std::string> pops_specs_t::defaults;

extern logger_t logger;
extern writer_t writer;

void pops_specs_t::read( const std::string & f )
{

  // if already done, skipe
  if ( specs.size() ) return;
  
  // ensure maps are initiated
  init();
  init_default();
  
  const bool use_default = f == "." ; 

  // do not allow a default feature file 
  if ( use_default ) Helper::halt( "no feature file specified" );
  
  if ( ( ! use_default ) && ( ! Helper::fileExists( Helper::expand( f ) ) ) ) 
    Helper::halt( "could not open " + f );
  
  // clear any current specifications
  
  specs.clear();

  // track features/channels (each can only be added once)
  std::set<std::string> checker;

  // track block usage
  std::set<std::string> bmap;
  bool level2 = false;

  std::ifstream IN1;

  if ( ! use_default )
    {      
      logger << "  reading feature specification from " << f << "\n";
      IN1.open( Helper::expand( f ).c_str() , std::ios::in );
    }
  else
    logger << "  using the default feature file\n";
  
  int cnt = 0;
  
  while ( 1 ) 
    {
      // all done?
      if ( use_default && cnt == defaults.size() ) break;

      std::string line;

      if ( use_default )
	line = defaults[ cnt++ ];
      else
	Helper::safe_getline( IN1 , line );
      
      if ( ( ! use_default ) && ( IN1.bad() || IN1.eof() ) ) break;
      
      if ( line == "" ) continue;
      if ( line[0] == '%' ) continue;
      
      // format
      // CH <label1> <label2> ... <sample-rate>
      // SELECT <blocks>
      // DROP
      // block: <feature> <channel-label> <key=val>
      
      std::vector<std::string> tok = Helper::parse( line , " \t" );
      if ( tok.size() == 0 ) continue;
      if ( tok.size() < 2 ) Helper::halt( "bad format for line: " + line );
      
      //
      // Channel specifier? find first that matches
      //
      
      if ( Helper::toupper( tok[0] ) == "CH" )
        {	  
	  if ( tok.size() < 4 )
            Helper::halt( "expecing: CH label {label2} {label3} ... SR UNIT" );
	  
	  // last two entries must be sample rate & unit
          int sr ;
          if ( ! Helper::str2int( tok[tok.size()-2] , &sr ) )
            Helper::halt( "bad format: " + line );
	  
	  std::string unit = tok[ tok.size()-1 ];

	  std::string primary_label = tok[1] ;

	  std::set<std::string> aliases;

	  //
	  // replace this label:: note, this ignores any aliases
	  //

	  if ( pops_opt_t::replacements.find( primary_label ) != pops_opt_t::replacements.end() )
	    {
	      primary_label = pops_opt_t::replacements[ primary_label ];
	      // nb. we'll need to replace from the spec file as well when reading
	    }
	  else // ... otherwise, also read in the aliases
	    {
	      // store any aliases (first fetching from the command line)
	      // nb. if the channel has been replaced, then aliases all go to this new channel

	      // from command line
	      aliases = pops_opt_t::aliases[ primary_label ];

	      // from spec file
	      for (int i=2;i<tok.size()-2;i++)
		aliases.insert( tok[i] );
	    }
	  
	  chs[ primary_label ] = pops_channel_t( primary_label , aliases , sr , unit ) ;
	  
          // next line
	  continue;
        }

      //
      // Final SELECT block command(s)
      //
      
      if ( Helper::toupper( tok[0] ) == "SELECT" )
	{
	  // if ( selected.size() > 0 )
	  //   Helper::halt( "can only use SELECT once" );
	  
	  for (int s=1; s<tok.size(); s++)
	    {
	      std::string selected_block = Helper::toupper( tok[s] ) ;
	      if ( bmap.find( selected_block ) == bmap.end() )
		Helper::halt( "could not find SELECT block " + tok[s] );
	      selected.insert( Helper::toupper( tok[s] ) );
	    }	  
	  continue;
	}
            

      //
      // Optional DROP command(s)
      //
      
      if ( Helper::toupper( tok[0] ) == "DROP" )
	{	  
	  // this expects 1+ final  variable names
	  for (int s=1; s<tok.size(); s++)
	    dropped.insert(  Helper::toupper( tok[s] ) );
	  continue;
	}

      
      //
      // Feature specifer: block name
      //
      
      const std::string & block1 = tok[0];
      if ( block1.size() == 1 || block1[ block1.size() - 1 ] != ':' )
	Helper::halt( "expecting colon after block name\n block: <feature> <args>" );      
      const std::string & block = Helper::toupper( tok[0].substr( 0 , block1.size() - 1 ) ) ;

      //
      // Feature type
      //

      const std::string & ftr = tok[1];
      
      if ( lab2ftr.find( Helper::toupper( ftr ) ) == lab2ftr.end() )
	Helper::halt( "feature not recognized: " + ftr );
      
      //
      // Channels, args
      //
      
      std::vector<std::string> tchs;
      std::map<std::string,std::string> targs;
      
      for (int i=2; i<tok.size(); i++)
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[i] , '=' );
	  if ( tok2.size() > 2 ) Helper::halt( "bad format: " + tok[i] );
	  
	  // special case: COVAR lists variable names
	  if ( ftr == "COVAR" )
	    {
	      targs[ tok2[0] ] = ""; // only requires the keys (var names)
	    }
	  // special case: COH has /pairs/ of channels
	  else if ( ftr == "COH" && tok2.size() == 1 )
	    {
	      // assuming CH1,CH2
	      std::vector<std::string> tok3 = Helper::parse( tok2[0] , "," );
	      if ( tok3.size() != 2 ) Helper::halt( "expecting comma-delimited pair of signals 'COH sig1,sig2'" );
	      
	      // check each
	      std::string channel1_label = tok3[0];
	      std::string channel2_label = tok3[1];

	      // is either being replaced? 
	      if ( pops_opt_t::replacements.find( channel1_label ) != pops_opt_t::replacements.end() )
		channel1_label = pops_opt_t::replacements[ channel1_label ];
	      if ( pops_opt_t::replacements.find( channel2_label ) != pops_opt_t::replacements.end() )
		channel2_label = pops_opt_t::replacements[ channel2_label ];
	      
	      // have both channels already been specified via CH?
	      if ( chs.find( channel1_label) == chs.end() )
		Helper::halt( channel1_label + " not specified via 'CH' yet: " + line );	      
	      if ( chs.find( channel2_label) == chs.end() )
		Helper::halt( channel2_label + " not specified via 'CH' yet: " + line );	      
	      
	      if ( channel1_label == channel2_label ) 
		Helper::halt( "cannot set COH channels to be the same: " 
			      + channel1_label + " " + channel2_label );
	      	      
	      // track paired label for mapping	      

	      const std::string paired_label = channel1_label + "," + channel2_label;
	      
	      tchs.push_back( paired_label );
	      
	      
	      // track
	      if ( checker.find( ftr + "::" + paired_label ) != checker.end() )
		Helper::halt( "can only specify a feature/channel pair once" );
	      
	      checker.insert( ftr + "::" + paired_label );
	      	      
	    }
	  // add as channel
	  else if ( tok2.size() == 1 )
	    {

	      std::string channel_label = tok2[0];

	      // is this being replaced? 
	      if ( pops_opt_t::replacements.find( channel_label ) != pops_opt_t::replacements.end() )
		channel_label = pops_opt_t::replacements[ channel_label ];
	      
	      // has the channel already been specified via CH?
	      if ( channel_label != "." && chs.find( channel_label) == chs.end() )
		Helper::halt( channel_label + " not specified via 'CH' yet: " + line );	      
	      
	      tchs.push_back( channel_label );

	      // track
	      if ( checker.find( ftr + "::" + channel_label ) != checker.end() )
		Helper::halt( "can only specify a feature/channel pair once" );

	      checker.insert( ftr + "::" + channel_label );
	      
	    }
	  else // else, as key=val arg
	    {
	      targs[ tok2[0] ] = tok2[1];
	    }
	}

      // if no channels, e.g. could be a time-track or a covariate; denote that it is empty
      if ( tchs.size() == 0 )
	tchs.push_back( "." );
      
      // check blocks: level 1 or 2?
      const bool level1 = pops_specs_t::lvl2.find( ftr ) == pops_specs_t::lvl2.end();

      // enforcee that all level-2 features must come first
      if ( ! level1 )
	level2 = true;
      else if ( level1 && level2 )
	Helper::halt( "cannot specify a level-1 feature after level-2 feature(s)" );
      
      // check block names
      if ( level1 || ftr == "TIME" )
	{
	  // special cases: lvl1 outliers command w/out channels -- set 'channel' as the new block
	  if ( ftr == "OUTLIERS" )
	    {	      	      
	      // requires that we've seen this block already
	      if ( bmap.find( block ) == bmap.end() )
		Helper::halt( "OUTLIERS specified block " + block + " not found" );
	      tchs.clear();
	      tchs.push_back( block );
	    }
	  else if ( ftr == "COVAR" )
	    {
	      tchs.clear();
              tchs.push_back( "." ); // just put empty channel
	    }

	  // mark that we've seen this block
	  bmap.insert( block );
	}
      else
	{
	  
	  // requires 'block' arg
	  if ( targs.find( "block" ) == targs.end() )
	    Helper::halt( "no block argument for " + ftr );
	  
	  // cannot specify an existing block as a target, unless it is the block= arg
	  // (i.e. replace originals)
	  const std::string from_block = Helper::toupper( targs[ "block" ] );
	  
	  // set 'channel' as the prior block
	  tchs.clear();
	  tchs.push_back( from_block );
	  
	  // does this point to an existing block?
	  if ( bmap.find( from_block ) == bmap.end() )
	    Helper::halt( "specified block " + targs[ "block" ] + " not found" );
	  
	  // not self-replacement
	  if ( from_block != block ) 
	    {
	      if ( bmap.find( block ) != bmap.end() )
		Helper::halt( "cannot specify an existing non-self block for a level-2 feature:\n" + line  );
	      // but now insert to track other lvl2 features
	      bmap.insert( block );
	    }	  
	}
      
      
      // add each channel separately (w/ the same args) 
      for (int c=0; c<tchs.size(); c++)
	{
	  pops_spec_t spec;
	  spec.block = block;
	  spec.ftr = lab2ftr[ Helper::toupper( ftr ) ];
	  spec.ch = tchs[c];
	  spec.arg = targs;
	  fcmap[ spec.ftr ][ spec.ch ] = spec;
	  specs.push_back( spec );
	}
    }

  if ( ! use_default ) 
    IN1.close();

  
  // check that at least some features were selected
  if ( selected.size() == 0 )
    Helper::halt( "no features SELECTed in " + f );

  // make sure features have any required args
  check_args();
  
  // track number of channels
  ns = chs.size();

  int nf = total_cols();
  // int nf_selected = select_cols();
  
  // & construct the map of specs/channels to feature columns
  build_colmap();
  
}


void pops_specs_t::init()
{

  // clear/reset
  lvl2.clear();
  lab2ftr.clear();
  ftr2lab.clear();
  defaults.clear();  
  ftr2ch2col.clear();     
  fcmap.clear();
  chs.clear();
  specs.clear();
  blocksize.clear();
  selected.clear();
  dropped.clear();
  col_block.clear();
  col_label.clear();
  col_original_label.clear();
  col_root.clear();
  col_select.clear();
  col_level.clear();
  orig2final.clear();
  final2orig.clear();
  n1 = na = nf = ns = 0; 

  // set
  
  lab2ftr[ "SPEC" ] = POPS_LOGPSD;
  lab2ftr[ "RSPEC" ] = POPS_RELPSD;
  lab2ftr[ "VSPEC" ] = POPS_CVPSD;
  
  lab2ftr[ "BAND" ] = POPS_BANDS;
  lab2ftr[ "RBAND" ] = POPS_RBANDS;
  lab2ftr[ "VBAND" ] = POPS_VBANDS;
  
  lab2ftr[ "COH" ] = POPS_COH;
  
  lab2ftr[ "SLOPE" ] = POPS_SLOPE;
  lab2ftr[ "SKEW" ] = POPS_SKEW;
  lab2ftr[ "KURTOSIS" ] = POPS_KURTOSIS;
  lab2ftr[ "HJORTH" ] = POPS_HJORTH;
  lab2ftr[ "FD" ] = POPS_FD;
  lab2ftr[ "PE" ] = POPS_PE;
  lab2ftr[ "MEAN" ] = POPS_MEAN;
  lab2ftr[ "OUTLIERS" ] = POPS_EPOCH_OUTLIER;
  lab2ftr[ "COVAR" ] = POPS_COVAR;

  lab2ftr[ "TIME" ] = POPS_TIME;
  lab2ftr[ "SMOOTH" ] = POPS_SMOOTH;
  lab2ftr[ "DENOISE" ] = POPS_DENOISE;
  lab2ftr[ "SVD" ] = POPS_SVD;
  lab2ftr[ "NORM" ] = POPS_NORM;
  lab2ftr[ "RESCALE" ] = POPS_RESCALE;
  lab2ftr[ "CUMUL" ] = POPS_CUMUL;
  lab2ftr[ "DERIV" ] = POPS_DERIV;

  ftr2lab[ POPS_LOGPSD ] = "SPEC";
  ftr2lab[ POPS_RELPSD ] = "RSPEC";
  ftr2lab[ POPS_CVPSD ] = "VSPEC";
  
  ftr2lab[ POPS_BANDS ] = "BAND";
  ftr2lab[ POPS_RBANDS ] = "RBAND";
  ftr2lab[ POPS_VBANDS ] = "VBAND";

  ftr2lab[ POPS_COH ] = "COH";

  ftr2lab[ POPS_SLOPE ] = "SLOPE";   
  ftr2lab[ POPS_SKEW ] = "SKEW";
  ftr2lab[ POPS_KURTOSIS ] = "KURTOSIS";
  ftr2lab[ POPS_HJORTH ] = "HJORTH";
  ftr2lab[ POPS_FD ] = "FD";      
  ftr2lab[ POPS_PE ] = "PE";  
  ftr2lab[ POPS_MEAN ] = "MEAN";
  ftr2lab[ POPS_EPOCH_OUTLIER ] = "OUTLIERS";
  ftr2lab[ POPS_COVAR ] = "COVAR";

  ftr2lab[ POPS_TIME ] = "TIME";  
  ftr2lab[ POPS_SMOOTH ] = "SMOOTH";
  ftr2lab[ POPS_DENOISE ] = "DENOISE";
  ftr2lab[ POPS_SVD ] = "SVD";
  ftr2lab[ POPS_NORM ] = "NORM";
  ftr2lab[ POPS_RESCALE ] = "RESCALE";
  ftr2lab[ POPS_CUMUL ] = "CUMUL";
  ftr2lab[ POPS_DERIV ] = "DERIV";


  // track level-2 features
  lvl2.insert( "TIME" );
  lvl2.insert( "SMOOTH" );
  lvl2.insert( "DENOISE" );
  lvl2.insert( "SVD" );
  lvl2.insert( "NORM" );
  lvl2.insert( "RESCALE" );
  lvl2.insert( "CUMUL" );
  lvl2.insert( "DERIV" );

}


int pops_specs_t::total_cols() const
{
  int n = 0;
  // for (int i=0; i<specs.size(); i++)
  //   specs[i].cols(&n);
  return n;
}

int pops_specs_t::select_cols() const
{
  return 0;
}



void pops_specs_t::check_args()
{
  for (int i=0; i<specs.size(); i++)
    {
      pops_spec_t & spec = specs[i];
      
      if ( spec.ftr == pops_feature_t::POPS_LOGPSD ||
	   spec.ftr == pops_feature_t::POPS_RELPSD ||
	   spec.ftr == pops_feature_t::POPS_CVPSD )
	{
	  if ( spec.arg.find( "lwr" ) == spec.arg.end() )
	    Helper::halt( ftr2lab[ pops_feature_t::POPS_LOGPSD ] + " requires 'lwr' arg" );
	  if ( spec.arg.find( "upr" ) == spec.arg.end() )
	    Helper::halt( ftr2lab[ pops_feature_t::POPS_LOGPSD ] + " requires 'upr' arg" );
	  if ( spec.narg( "lwr" ) > spec.narg( "upr" ) )
	    Helper::halt( ftr2lab[ pops_feature_t::POPS_LOGPSD ] + " requires 'lwr' < 'upr' " );
	  if ( spec.narg( "lwr" ) <= 0 || spec.narg( "upr" ) <= 0 )
	    Helper::halt( ftr2lab[ pops_feature_t::POPS_LOGPSD ] + " requires 'lwr' and 'upr' to be > 0 " );
	}

      // the z-lwr/z-upr range does not need to overlap lwr/upr range for RELPSD
      if ( spec.ftr == pops_feature_t::POPS_RELPSD )
	{
	  if ( spec.arg.find( "z-lwr" ) == spec.arg.end() )
            Helper::halt( ftr2lab[ pops_feature_t::POPS_LOGPSD ] + " requires 'z-lwr' arg" );
          if ( spec.arg.find( "z-upr" ) == spec.arg.end() )
            Helper::halt( ftr2lab[ pops_feature_t::POPS_LOGPSD ] + " requires 'z-upr' arg" );
          if ( spec.narg( "z-lwr" ) > spec.narg( "z-upr" ) )
            Helper::halt( ftr2lab[ pops_feature_t::POPS_LOGPSD ] + " requires 'z-lwr' < 'z-upr' " );
          if ( spec.narg( "z-lwr" ) <= 0 || spec.narg( "z-upr" ) <= 0 )
            Helper::halt( ftr2lab[ pops_feature_t::POPS_LOGPSD ] + " requires 'z-lwr' and 'z-upr' to be > 0 " );
	}
      
      // PE
      if ( spec.ftr == pops_feature_t::POPS_PE )
        {
	  if ( spec.arg.find( "from" ) == spec.arg.end() ||
	       spec.arg.find( "to" ) == spec.arg.end() )
	    Helper::halt( "requires from=X to=Y" );
	  int n1 = spec.narg( "from" );
	  int n2 = spec.narg( "to" );
	  if ( n2 < n1 || n1 < 3 || n1 > 7 || n2 < 3 || n2 > 7 )
	    Helper::halt( "from=x and to=y must be between 3 and 7" );	  
        }

      // COVAR (individual-level)
      if ( spec.ftr == pops_feature_t::POPS_COVAR )
	{
	  if ( spec.arg.size() == 0 ) 
	    Helper::halt( "COVAR requires 1+ variable names listed after" );
	}
      
      // time-tracks
      if ( spec.ftr == pops_feature_t::POPS_TIME )
	{
	  if ( spec.arg.find( "order" ) == spec.arg.end() )
	    spec.arg[ "order" ] = "1";
	}

      // smoothing/denoising
      if ( spec.ftr == pops_feature_t::POPS_DENOISE )
	{
	  if ( spec.arg.find( "lambda" ) == spec.arg.end() )
            Helper::halt( ftr2lab[ pops_feature_t::POPS_DENOISE ] + " requires 'lambda' arg" );	  
	}
      
      if ( spec.ftr == pops_feature_t::POPS_SMOOTH )
	{
	  if ( spec.arg.find( "half-window" ) == spec.arg.end() )
            Helper::halt( ftr2lab[ pops_feature_t::POPS_SMOOTH ] + " requires 'half-window' (epochs) arg" );
	  // can also have 'a' argument -- 0 to 1
	  if ( spec.arg.find( "a" ) == spec.arg.end() )
	    {
	      double a = spec.narg( "a" );
	      if ( a < 0 || a > 1 ) Helper::halt( "expecting 'a' arg to be between 0 and 1" );
	    }
	}
      
      // CUMUL
      if ( spec.ftr == pops_feature_t::POPS_CUMUL )
	{	  
	  if ( spec.arg.find( "type" ) == spec.arg.end() )
	    {
	      spec.arg[ "type" ] = "norm" ; 
	    }
	  else if ( ! ( spec.arg[ "type" ] == "pos" 
			|| spec.arg[ "type" ] == "neg" 
			|| spec.arg[ "type" ] == "abs" ) )
	    Helper::halt( ftr2lab[ pops_feature_t::POPS_CUMUL ] + " requires 'type' as pos,neg or abs" );
	}
      
      // DERIV
      if ( spec.ftr == pops_feature_t::POPS_DERIV )
	{
	  if ( spec.arg.find( "half-window" ) == spec.arg.end() )
	    Helper::halt( ftr2lab[ pops_feature_t::POPS_DERIV ] + " requires 'half-window' (epochs) ");	  
	  const int hw = spec.narg( "half-window" );
	  if ( hw <= 0 || hw > 100 ) Helper::halt( "expecting half-window between 1 and 100 for DERIV" );
	}

      
      // SVD
      if ( spec.ftr == pops_feature_t::POPS_SVD )
	{
	  if ( spec.arg.find( "nc" ) == spec.arg.end() )
	    Helper::halt( ftr2lab[ pops_feature_t::POPS_SVD ] + " requires 'nc' arg" );
	  if ( spec.arg.find( "file" ) == spec.arg.end() )
	    Helper::halt( ftr2lab[ pops_feature_t::POPS_SVD ] + " requires 'file' arg" );
	}

      // OUTLIERS
      if ( spec.ftr == pops_feature_t::POPS_EPOCH_OUTLIER )
	{
	  if ( spec.arg.find( "th" ) == spec.arg.end() )
	    Helper::halt( ftr2lab[ pops_feature_t::POPS_EPOCH_OUTLIER ] + " requires 'th' arg" );	  
	}
      
    }
}


void pops_specs_t::build_colmap()
{
  
  // at this point, all specs are loaded; step through to figure out the implied
  // columns, at three levels:
  //   level 1 (things actually extracted from the EDF, and saved to the data files)
  //   level 2 (derived metrics calculated on loading (all) trainers, e.g. smoothing, etc)
  //   selected : final set of selected metrics : may be a subset of the above 
  //   dropped  : drop any indiv variables
  
  //  e.g. eeg: SPEC C3
  //       eog: SPEC LOC
  //       slope: SLOPE C3
  //       
  //       slope: NORM block=slope
  //       eeg2: SMOOTH block=eeg
  //       eog: SMOOTH block=eog

  //       SELECT eeg2 eog slope

  // expanded feature lists

  //std::vector<int> block_columns( const std::string & );  

  col_block.clear();
  col_label.clear();
  col_original_label.clear(); // if replace=X,Y track originals
  col_root.clear();
  col_select.clear();
  col_level.clear();
  
  ftr2ch2col.clear();
  
  int n = 0;
  
  for (int i=0; i<specs.size(); i++)
    {

      int start = n;
      specs[i].cols(&n);
      int end = n; // one past end
            
      const pops_feature_t ftr = specs[i].ftr;
      const std::string & ftrlab = pops_specs_t::ftr2lab[ ftr ];
      const std::string & ch = specs[i].ch;
      const std::string & block = specs[i].block;
      const bool level1 = pops_specs_t::lvl2.find( ftrlab ) == pops_specs_t::lvl2.end();
      
      // track block size
      blocksize[ block ] += specs[i].size;
      
      for (int j=start; j<end; j++)
	{
	  // track that this feature mapped to this channel to this X1 col
	  ftr2ch2col[ ftr ][ ch ].push_back( j );

	  // final label
	  const std::string vlabel = ftrlab + "." + ch + ".V" + Helper::int2str( j - start + 1 );

	  col_label.push_back( vlabel );
	  
	  if ( pops_opt_t::replacements_rmap.find( ch ) != pops_opt_t::replacements_rmap.end() )
	    col_original_label.push_back( ftrlab + "." + pops_opt_t::replacements_rmap[ ch ] + ".V" + Helper::int2str( j - start + 1 ) );
	  else
	    col_original_label.push_back( ftrlab + "." + ch + ".V" + Helper::int2str( j - start + 1 ) );
	  
	  col_root.push_back( ftrlab + "." + ch ); // i.e. channel/block specific
	  
	  col_block.push_back( block );	  

	  col_level.push_back( level1 ? 1 : 2 );

	  // black selected?
	  bool is_selected = selected.find( block ) != selected.end() ;

	  // but this variable dropped?
	  if ( dropped.find( vlabel ) != dropped.end() ) is_selected = false;
	  
	  col_select.push_back( is_selected );
	}
      
    }

  //
  // denote finals
  //

  n = col_block.size();

  // not sure this is needed...
  col_level.resize( n );
  col_select.resize( n );

  
  n1 = 0; // level 1 features
  na = 0; // level 1 + level 2 features (all)
  nf = 0; // final number of selected features
  
  int p = 0;
  for (int f=0; f<col_block.size(); f++)
    {

      if ( col_select[f] )
	{
	  orig2final[ f ] = p;
	  final2orig[ p ] = f;
	  ++p;
	}

      // only level-1 features
      if ( col_level[f] == 1 ) ++n1;

      // all features
      ++na;

      // only seleted (l1+l2) features
      if ( col_select[f] ) ++nf;
      
      //
      // dump to output
      //

      writer.level( f+1 , globals::feature_strat );
      writer.value( "BLOCK" , col_block[f]  );
      writer.value( "INC" , (int)col_select[f] );
      if ( col_select[f] )
	{
	  writer.value( "FINAL" , orig2final[ f ] + 1 );
	  //writer.value( "ORIG" , final2orig[ orig2final[ f ] ] + 1 );
	}      
      writer.value( "LEVEL" , col_level[f] );
      writer.value( "LABEL" , col_label[f] );
      writer.value( "LABEL_ORIG" , col_original_label[f] );
      writer.value( "ROOT" ,  col_root[f] );
      
      writer.unlevel( globals::feature_strat );

    }

  logger << "   " << n1 << " level-1 features, "
	 << na-n1 << " level-2 features\n"
	 << "   " << nf << " of " << na << " features selected in the final feature set\n";
  
}


// give columns for a spec/channel combo
bool pops_specs_t::has( pops_feature_t ftr , const std::string & ch )
{
  std::map<pops_feature_t,std::map<std::string,std::vector<int> > >::const_iterator cc = ftr2ch2col.find( ftr );
  if ( cc == ftr2ch2col.end() ) return false;
  std::map<std::string,std::vector<int> >::const_iterator dd = cc->second.find( ch );
  return dd != cc->second.end(); 
}

std::vector<int> pops_specs_t::cols( pops_feature_t ftr , const std::string & ch )
{
  std::vector<int> dummy;
  std::map<pops_feature_t,std::map<std::string,std::vector<int> > >::const_iterator cc = ftr2ch2col.find( ftr );
  if ( cc == ftr2ch2col.end() ) return dummy;
  std::map<std::string,std::vector<int> >::const_iterator dd = cc->second.find( ch );
  if ( dd == cc->second.end() ) return dummy;
  return dd->second;
}


std::vector<std::string> pops_specs_t::select_labels()
{
  std::vector<std::string> s;
  std::map<int,int>::const_iterator ii = final2orig.begin();
  while ( ii != final2orig.end() )
    {
      s.push_back( col_label[ ii->second ] );
      ++ii;
    }
  return s;
}

std::vector<std::string> pops_specs_t::select_original_labels()
{
  std::vector<std::string> s;
  std::map<int,int>::const_iterator ii = final2orig.begin();
  while ( ii != final2orig.end() )
    {
      s.push_back( col_original_label[ ii->second ] );
      ++ii;
    }
  return s;
}

std::vector<std::string> pops_specs_t::select_roots()
{
  std::vector<std::string> s;
  std::map<int,int>::const_iterator ii = final2orig.begin();
  while ( ii != final2orig.end() )
    {
      s.push_back( col_root[ ii->second ] );
      ++ii;
    }
  return s;
}


std::vector<std::string> pops_specs_t::select_blocks()
{
  std::vector<std::string> s;
  std::map<int,int>::const_iterator ii = final2orig.begin();
  while ( ii != final2orig.end() )
    {
      s.push_back( col_block[ ii->second ] );
      ++ii;
    }
  return s;
  
}



// return implied number of columns

int pops_spec_t::cols( int * t )
{
  
  // PSD is stratified by frequency
  if ( ftr == POPS_LOGPSD
       || ftr == POPS_RELPSD
       || ftr == POPS_CVPSD )
    {
      double lwr = narg( "lwr" );
      double upr = narg( "upr" );
      int n = ( upr - lwr ) / pops_opt_t::spectral_resolution + 1 ;
      size = n;
      *t += n;
      return n;
    }
  
  // 6 fixed bands (per channel, or between /pair/ of channels for COH)
  if ( ftr == POPS_BANDS 
       || ftr == POPS_RBANDS 
       || ftr == POPS_VBANDS 
       || ftr == POPS_COH )
    {
      *t += 6;
      size = 6;
      return size;
    }
  

  // 1 column per channel
  if ( ftr == POPS_SLOPE
       || ftr == POPS_SKEW
       || ftr == POPS_KURTOSIS
       || ftr == POPS_FD
       || ftr == POPS_MEAN )
    {
      *t += 1;
      size = 1;
      return size;
    }
  
  // 2 or 3 values per channel
  // (only include H1 is 'h1=1' option set
  if ( ftr == POPS_HJORTH )
    {
      int n = narg( "h1" );
      size = n > 0.5 ? 3 : 2 ;
      *t += size ;
      return size ;
    }
  
  // PE is 3..7 
  if ( ftr == POPS_PE )
    {
      int n1 = narg( "from" );
      int n2 = narg( "to" );
      size = n2 - n1 + 1 ;
      *t += size ;
      return size ;
    }

  // COVAR
  if ( ftr == POPS_COVAR )
    {
      size = arg.size();
      *t += size;
      return size;
    }
  
  // time-track
  if ( ftr == POPS_TIME )
    {
      int n = narg( "order" );
      if ( n < 1 )
	Helper::halt( "invalid value for TIME order (1-4)" );
      if ( n > 4 )
	Helper::halt( "invalid value for TIME order (1-4)" );
      *t += n;
      return n;
    }
  
  // SVD
  if ( ftr == POPS_SVD )
    {
      size = narg( "nc" );
      *t += size;
      return size;
    }

  // row removal
  if ( ftr == POPS_EPOCH_OUTLIER )
    {
      size = 0;
      return size;
    }
  
  // block transformations (in-place, or copy)
  if ( ftr == POPS_SMOOTH ||
       ftr == POPS_DENOISE ||
       ftr == POPS_RESCALE ||
       ftr == POPS_CUMUL || 
       ftr == POPS_DERIV || 
       ftr == POPS_NORM )
    {

      std::string from_block = Helper::toupper( arg[ "block" ] );

      // either 0 (if an inplace transformation)
      if ( Helper::toupper( block ) == from_block )
	{
	  size = 0;
	  return size;
	}
      
      // or, the 'size' will the original block will be duplicated
      size = pops_specs_t::blocksize[ from_block ];
      *t += size;
      return size;
    }


  
  // unknown/error at this point  
  Helper::halt( "internal error extracting column count for "
		+ pops_specs_t::ftr2lab[ ftr ] );
  
  return 0;
  
}

std::vector<int> pops_specs_t::block_cols( const std::string & b , int n )
{  
  std::vector<int> r;  
  for (int f=0; f<n; f++)
    {
      if ( col_block[f] == b )
	r.push_back( f );
    }  
  return r;
}


void pops_specs_t::init_default()
{
  defaults.clear();
  // single EEG
  defaults.push_back( "CH C4_M1 C4 C4-M1 C4_A1 C4-A1  128 uV" );

  // level 1
  defaults.push_back( "spec1: SPEC C4_M1 lwr=0.5 upr=35" );
  defaults.push_back( "spec2: RSPEC C4_M1 lwr=2 upr=15 z-lwr=30 z-upr=45" );
  defaults.push_back( "misc: SLOPE C4_M1" );
  defaults.push_back( "misc: SKEW C4_M1" );
  defaults.push_back( "misc: KURTOSIS C4_M1" );
  defaults.push_back( "misc: FD C4_M1" );
  defaults.push_back( "misc: PE C4_M1" );
  defaults.push_back( "hjorth: HJORTH C4_M1" );
  
  // lvl1 outlier removal
  defaults.push_back( "hjorth: OUTLIERS th=8" );
  
  // level 2
  defaults.push_back( "svd1: SVD nc=10 block=spec1" );
  defaults.push_back( "svd2: SMOOTH block=svd1 half-window=7" );
  defaults.push_back( "misc2: SMOOTH block=misc half-window=7" );

  // final
  defaults.push_back( "SELECT svd1 svd2 spec2 misc misc2" );
  
}


#endif
