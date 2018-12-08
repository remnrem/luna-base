
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


#include "luna.h"


// ----------------------------------------------------------------------------------------
//
// Process commands from STDIN
//
// ----------------------------------------------------------------------------------------

void cmd_t::replace_wildcards( const std::string & id )
{
  // replace in all 'params' any instances of 'globals::indiv_wildcard' with 'id'
  params = original_params;
  for (int p = 0 ; p < params.size(); p++ ) params[p].update( id , globals::indiv_wildcard );
}


bool cmd_t::read( const std::string * str , bool silent )
{
  
  bool cmdline_mode = str == NULL;   

  if ( std::cin.eof() && cmdline_mode ) return false;
  
  if ( (!cmdline_mode) && str->size() == 0 ) return false;

  reset();
  
  // CMD param=1 p=1,2,3 f=true out=o1 & CMD2 etc ; 
  
  // EITHER read from std::cin, 
  // OR          from -s command line
  // OR          from a string

  std::istringstream allinput;

  if ( ! cmdline_mode ) // read from 'str' 
    {
      // split by lines
      std::vector<std::string> tok = Helper::parse( *str , "\n" );
      
      std::stringstream ss;
      
      for (int l=0;l<tok.size();l++)
	{
	  if ( tok[l] == "" ) continue;
	  if ( l != 0 ) ss << " & ";
	  ss << tok[l];
	}     
      allinput.str( ss.str() );

    }
  else if ( cmd_t::cmdline_cmds == "" )
    {
      std::stringstream ss;
      bool first_cmd = true;
      while ( 1 )
	{
	  std::string s;
	  std::getline( std::cin , s , '\n' );
	  if ( std::cin.eof() ) break;
	  if ( s == "" ) continue;	  
	  if ( ! first_cmd ) ss << " & ";
	  ss << s ;
	  first_cmd = false; 
	}      
      allinput.str( ss.str() );
    }

  else 
    allinput.str( cmd_t::cmdline_cmds );
      
  //std::cout << "allinput [" << allinput.str() << "]\n";

  // read up to the ';'
  //  std::getline( std::cin , line , ';' );
  std::getline( allinput , line , ';' );
  

  // change any '&' to '\n'
  for (int i=0;i<line.size();i++) if ( line[i] == '&' ) line[i] = '\n';

  // skip comments between commands if line starts with '#' (or '\n')
  bool recheck = true;
  while (recheck)
    {     
      if ( line[0] == '#' || line[0] == '\n' ) 
	{
	  line = line.substr( line.find("\n")+1);
	}
      else recheck = false;
    }
  
  // swap in any variables
  Helper::swap_in_variables( &line , vars );

  
  std::vector<std::string> tok = Helper::quoted_parse( line , "\n" );
  if ( tok.size() == 0 ) 
    {
      quit(true);
      return false;
    }

  // command(s)
  for (int c=0;c<tok.size();c++)
    {

      // all done?
      // if ( Helper::iequals( tok[c] , "QUIT" ) ) 
      // 	{
      // 	  quit( true );
      // 	  return false;
      // 	}
      
      std::vector<std::string> ctok = Helper::quoted_parse( tok[c] , "\t " );
      if ( ctok.size() < 1 ) return false;
      cmds.push_back( ctok[0] );
      param_t param;
      for (int j=1;j<ctok.size();j++) param.parse( ctok[j] );
      params.push_back( param );
    }

  // make a copy of the 'orginal' params
  original_params = params;

  // summary

  if ( ! globals::silent ) 
    {
      std::cerr << "-------------------------------------------------------------------\n";
      std::cerr << "EDF source: " << input << "\n";
      std::cerr << "Output DB : " << writer.name() << "\n";
      for (std::set<std::string>::iterator s=signallist.begin();
	   s!=signallist.end();s++) 
	std::cerr << "signals: " << *s << "\n";
      for (int i=0;i<cmds.size();i++)
	{
	  std::cerr << " c" << i+1 
		    << "\t" << cmds[i] << "\t"
		    << params[i].dump("","|") 
		    << "\n";
	}
      std::cerr << "-------------------------------------------------------------------\n";
    }

  return true;
}



//
// Evaluate commands
//

bool cmd_t::eval( edf_t & edf ) 
{

  //
  // Loop over each command
  //
  
  for ( int c = 0 ; c < num_cmds() ; c++ )
    {	        
      
      // was a problem flag raised when loading the EDF?
      
      if ( globals::problem ) return false;
      
      //
      // If this particular command did not explicitly specify
      // signals, then add all
      //
      
      if ( ! param(c).has( "signal" ) )
	param(c).add( "signal" , signal_string() );
      
      //
      // Attach any new annotations (e.g. sleep stages, etc) to
      // the timeline if requested; these will persist as long as
      // the EDF lasts
      //
      
      std::set<std::string> required_annotations = param(c).strset( "annot" , "," );
      std::set<std::string>::const_iterator ai = required_annotations.begin();	  
      
      while ( ai != required_annotations.end() )
	{
	  attach_annot( edf , *ai );
	  ++ai;
	}
      
      
      //
      // Print command
      //

      if ( ! globals::silent ) 
	std::cerr << " CMD #" << c+1 << ": " << cmd(c) << "\n";
      
      writer.cmd( cmd(c) , c+1 , param(c).dump( "" , " " ) );
      
      
      //
      // Now process the command
      //
      
      if      ( is( c, "WRITE" ) )        proc_write( edf, param(c) );
      else if ( is( c, "SUMMARY" ) )      proc_summaries( edf , param(c) );
      else if ( is( c, "DESC" ) )         proc_desc( edf , param(c) );
      else if ( is( c, "STATS" ) )        proc_validate( edf , param(c) );
      
      else if ( is( c, "REFERENCE" ) )    proc_reference( edf , param(c) );
      else if ( is( c, "FLIP" ) )         proc_flip( edf , param(c) );
      else if ( is( c, "uV" ) )           proc_scale( edf , param(c) , "uV" ); 
      else if ( is( c, "mV" ) )           proc_scale( edf , param(c) , "mV" );
      else if ( is( c, "RECORD-SIZE" ) )  proc_rerecord( edf , param(c) );
      
      else if ( is( c, "TIME-TRACK" ) )   proc_timetrack( edf, param(c) );
      
      else if ( is( c, "STAGE" ) || is( c, "HYPNO" ) ) proc_sleep_stage( edf , param(c) );
      
      else if ( is( c, "TSLIB" ) )        pdc_t::construct_tslib( edf , param(c) );
      else if ( is( c, "SSS" ) )          pdc_t::simple_sleep_scorer( edf , param(c) );
      else if ( is( c, "EXE" ) )          pdc_t::similarity_matrix( edf , param(c) );
      
      else if ( is( c, "DUMP" ) )         proc_dump( edf, param(c) );
      else if ( is( c, "DUMP-RECORDS" ) ) proc_record_dump( edf , param(c) );
      else if ( is( c, "DUMP-EPOCHS" ) )  proc_epoch_dump( edf, param(c) );
      else if ( is( c, "ANNOTS" ) )       proc_list_all_annots( edf, param(c) );
      else if ( is( c, "MATRIX" ) )       proc_epoch_matrix( edf , param(c) );
      else if ( is( c, "RESTRUCTURE" ) )  proc_restructure( edf , param(c) );
      else if ( is( c, "SIGNALS" ) )      proc_drop_signals( edf , param(c) );
      else if ( is( c, "RMS" ) || is( c, "SIGSTATS" ) ) proc_rms( edf, param(c) );
      else if ( is( c, "MSE" ) )          proc_mse( edf, param(c) );
      else if ( is( c, "LZW" ) )          proc_lzw( edf, param(c) );
      else if ( is( c, "ZR" ) )           proc_zratio( edf , param(c) );
      else if ( is( c, "ANON" ) )         proc_anon( edf , param(c) );
      else if ( is( c ,"EPOCH" ) )        proc_epoch( edf, param(c) );
      else if ( is( c ,"SLICE" ) )        proc_slice( edf , param(c) , 1 );
      
      else if ( is( c, "MASK" ) )         proc_mask( edf, param(c) );
      else if ( is( c, "FILE-MASK" ) )    proc_file_mask( edf , param(c) );
      else if ( is( c, "DUMP-MASK" ) )    proc_dump_mask( edf, param(c) );
      
      else if ( is( c, "EPOCH-ANNOT" ) )  proc_file_annot( edf , param(c) );
      else if ( is( c, "EPOCH-MASK" ) )   proc_epoch_mask( edf, param(c) );
      
      
      else if ( is( c, "FILTER" ) )       proc_filter( edf, param(c) );
      else if ( is( c, "FILTER-DESIGN" )) proc_filter_design( edf, param(c) );
      //	  else if ( is( c, "LEGACY-FILTER" )) proc_filter_legacy( edf, param(c) );
      else if ( is( c, "TV" ) )           proc_tv_denoise( edf , param(c) );
      
      else if ( is( c, "COVAR" ) )        proc_covar( edf, param(c) );
      else if ( is( c, "PSD" ) )          proc_psd( edf, param(c) );	  
      else if ( is( c, "1FNORM" ) )       proc_1overf_norm( edf, param(c) );
      
      else if ( is( c, "FIP" ) )          proc_fiplot( edf , param(c) );
      
      else if ( is( c, "COH" ) )          proc_coh( edf , param(c) );
      else if ( is( c, "LECGACY-COH" ) )  proc_coh_legacy( edf , param(c) );
      else if ( is( c, "CORREL" ) )       proc_correl( edf , param(c) );
      else if ( is( c, "ED" ) )           proc_elec_distance( edf , param(c) );
      else if ( is( c, "ICA" ) )          proc_ica( edf, param(c) );
      else if ( is( c, "L1OUT" ) )        proc_leave_one_out( edf , param(c) );
      
      else if ( is( c, "EMD" ) )          proc_emd( edf , param(c) );
  
      else if ( is( c, "MI" ) )           proc_mi( edf, param(c) );
      else if ( is( c, "HR" ) )           proc_bpm( edf , param(c) );
      else if ( is( c, "SUPPRESS-ECG" ) ) proc_ecgsuppression( edf , param(c) );
      else if ( is( c, "PAC" ) )          proc_pac( edf , param(c) );
      else if ( is( c, "CFC" ) )          proc_cfc( edf , param(c) );
      else if ( is( c, "TAG" ) )          proc_tag( param(c) );
      else if ( is( c, "RESAMPLE" ) )     proc_resample( edf, param(c) );
      else if ( is( c, "SPINDLES" ) )     proc_spindles( edf, param(c) );	  
      else if ( is( c, "POL" ) )          proc_polarity( edf, param(c) );	  
      
      else if ( is( c, "SW" ) || is( c, "SLOW-WAVES" ) ) proc_slowwaves( edf, param(c) );
      else if ( is( c, "ARTIFACTS" ) )    proc_artifacts( edf, param(c) );
      else if ( is( c, "SPIKE" ) )        proc_spike( edf , param(c) );
      
      else 
	{
	  Helper::halt( "did not recognize command: " + cmd(c) );
	  return false; 
	}
      
      
      //
      // Was a problem flag set?
      //
      
      if ( globals::problem ) 
	{
	  if ( ! globals::silent ) 
	    std::cerr << "**warning: the PROBLEM flag was set, skipping to next EDF...\n";

	  return false;
	}
         
      
    } // next command
  


  return true;
}



//
// ----------------------------------------------------------------------------------------
//
// Wrapper (proc_*) functions below that drive specific commands
//
// ----------------------------------------------------------------------------------------
//

//
// HEADERS  (vmode 1) : print headers
// VALIDATE (vmode 2) : validate entire EDF
//


void proc_summaries( const std::string & edffile , const std::string & rootname , 
		     int v_mode , cmd_t & cmd , param_t & param )
{

  const std::set<std::string> * inp_signals = NULL;      
  if ( cmd.signals().size() > 0 ) inp_signals = &cmd.signals();
  edf_t edf;  
  
  bool okay = edf.attach( edffile , rootname , inp_signals );

  if ( ! okay ) 
    if ( ! globals::silent ) 
      std::cerr << "Problem: " << rootname << " reading EDF " << edffile << "\n";    
  
  // terse HEADER summary
  if ( v_mode == 1 ) edf.terse_summary() ;
  
  // or a full VALIDATE query?
  if ( v_mode == 2 ) 
    {
      bool result = edf.validate( param );
      // std::cout << rootname << "\t"
      // 		<< "VALIDATION\t"
      // 		<< ( result ? "OKAY" : "PROBLEM" ) << "\n";
            
    }
}


// SUMMARY : summarize EDF files (verbose, human-readable)  

void proc_summaries( edf_t & edf , param_t & param )
{
  if ( ! globals::silent ) 
    std::cout << "EDF filename   : " << edf.filename << "\n" 
	      << edf.header.summary() << "\n"
	      << "----------------------------------------------------------------\n\n";
}


// DESC : very brief summary of contents 

void proc_desc( edf_t & edf , param_t & param )
{
  edf.description();
}

// STATS : get basic stats for an EDF

void proc_validate( edf_t & edf , param_t & param )
{
  edf.validate( param );
}

// RMS/SIGSTATS : calculate root-mean square for signals 

void proc_rms( edf_t & edf , param_t & param )
{    
  rms_per_epoch( edf , param );
}

// MSE : calculate multi-scale entropy per epoch

void proc_mse( edf_t & edf , param_t & param )
{    
  mse_per_epoch( edf , param );
}

// LZW : compression index per epoch/signal
void proc_lzw( edf_t & edf , param_t & param )
{    
  lzw_per_epoch( edf , param );
}


// ZR : Z-ratio 
void proc_zratio( edf_t & edf , param_t & param )
{    
  std::string signal = param.requires( "signal" );

  staging_t staging;
  staging.zratio.calc( edf , signal );
}


// ARTIFACTS : artifact rejection using Buckelmueller et al. 

void proc_artifacts( edf_t & edf , param_t & param )	  
{
  std::string signal = param.requires( "signal" );
  annot_t * a = buckelmuller_artifact_detection( edf , param , signal );  
}


// LEGACY-FILTER : band-pass filter, band-pass filter only

// void proc_filter_legacy( edf_t & edf , param_t & param )	  
// {
//   //band_pass_filter( edf , param );
// }

// FILTER : general FIR

void proc_filter( edf_t & edf , param_t & param )	  
{
  dsptools::apply_fir( edf , param );
}


// -fir  from the command line
void proc_filter_design_cmdline()
{
  
  // expect parameters on STDIN
  
  param_t param;
  while ( ! std::cin.eof() )
    {
      std::string x;
      std::cin >> x;      
      if ( std::cin.eof() ) break;
      if ( x == "" ) continue;
      param.parse( x ); 
    }

  dsptools::design_fir( param );
  
}

// TV   total variation 1D denoising
void proc_tv_denoise( edf_t & edf , param_t & param )
{
  dsptools::tv( edf , param );
}


// -cwt  from the command line
void proc_cwt_design_cmdline()
{
  
  // expect parameters on STDIN
  
  param_t param;
  while ( ! std::cin.eof() )
    {
      std::string x;
      std::cin >> x;      
      if ( std::cin.eof() ) break;
      if ( x == "" ) continue;
      param.parse( x ); 
    }
  
  dsptools::design_cwt( param );
  
}


// FILTER-DESIGN : general FIR design

void proc_filter_design( edf_t & edf , param_t & param )	  
{
  dsptools::design_fir( param );
}


// RESAMPLE : band-pass filter

void proc_resample( edf_t & edf , param_t & param ) 
{
  dsptools::resample_channel( edf, param );
}

// PSD : calculate PSD 

void proc_psd( edf_t & edf , param_t & param )	  
{  
  std::string signal = param.requires( "signal" );
  annot_t * power = spectral_power( edf , signal , param );  
}

// 1FNORM : normalization of signals for the 1/f trend

void proc_1overf_norm( edf_t & edf , param_t & param )	  
{  
  dsptools::norm_1overf( edf,  param ) ;
}

// FI-plot : frequency/interval plot

void proc_fiplot( edf_t & edf , param_t & param )	  
{  
  fiplot_wrapper( edf , param );
}

// TAG : analysis tag

void proc_tag( param_t & param )
{
  set_tag( param.requires( "tag" ) );
}

void set_tag( const std::string & t ) 
{
  globals::current_tag = t ; 

  if ( ! globals::silent ) 
    if ( t != "." ) 
      if ( ! globals::silent ) 
	std::cerr << " setting analysis tag to [" << globals::current_tag << "]\n";

  if ( t == "." ) writer.tag( "." , "." );
  else
    {
      std::vector<std::string> tok = Helper::parse( globals::current_tag , "/" );
      if ( tok.size() != 2 ) Helper::halt( "TAG format should be level/factor" );
      writer.tag( tok[0] , tok[1] );
    }
}

// ANON : anonymize EDF 

void proc_anon( edf_t & edf , param_t & param )
{
  if ( ! globals::silent ) 
    std::cerr << " setting subject ID, start time/date to null ('.') for EDF " 
	      << edf.filename << "\n";
  
  edf.header.patient_id = ".";
  edf.header.starttime = ".";
  edf.header.startdate = ".";
}


// DUMP : dump all data

void proc_dump( edf_t & edf , param_t & param )	  
{
  std::string signal = param.requires( "signal" );  
  edf.data_dumper( signal , param );	  
}
      

// EPOCH DUMP 

void proc_epoch_dump( edf_t & edf , param_t & param )
{

  std::set<std::string> * annots = NULL;
  if ( param.has( "annot" ) )
    {
      annots = new std::set<std::string>;
      *annots = param.strset( "annot" ); // parse comma-delim list
    }

  edf.data_epoch_dumper( param , annots );
}


// MATRIX 

void proc_epoch_matrix( edf_t & edf , param_t & param )
{
  std::set<std::string> * annots = NULL;
  if ( param.has( "annot" ) )
    {
      annots = new std::set<std::string>;
      *annots = param.strset( "annot" ); // parse comma-delim list
    }
  if ( param.has("all-annots") )
    {
      std::vector<std::string> Xannots = edf.timeline.annotations.names();
      if ( Xannots.size() > 0 ) 
	{
	  annots = new std::set<std::string>;
	  for (int i=0;i<Xannots.size();i++) annots->insert( Xannots[i] ) ;	  
	}

      // TODO: clarify epoch annotations versus generic annotations 
      //  std::set<std::string> epoch_annotations = edf.timeline.epoch_annotations();
      //  std::set<std::string>::const_iterator ii = epoch_annotations.begin();
      //  while ( ii != epoch_annotations.end() ) { std::cout << "S " << *ii << "\n"; ++ii; } 

    }
  
  edf.epoch_matrix_dumper( param , annots );
  
  if ( annots != NULL) delete annots;

}


// INTERVALS : raw signal data from an interval list

void proc_intervals( param_t & param , const std::string & data )	  
{  

  std::string ints = param.requires( "intervals" );
  // e.g.: INTERVAL edfs=../scratch/edf.list 
  // intervals from STDIN
  dump_intervals( ints , data );
}



// COVAR : covariance between two signals (not implemented)

void proc_covar( edf_t & edf , param_t & param )
{  
  std::string signals1 = param.requires( "signal1" );
  std::string signals2 = param.requires( "signal2" );
  edf.covar(signals1,signals2);
}


// SPINDLES : spindle detection using CWT or bandpass/RMS

void proc_spindles( edf_t & edf , param_t & param )
{	

  // default to wavelet analysis
  std::string method = param.has( "method" ) ? param.value( "method" ) : "wavelet" ; 
  
  annot_t * a = NULL;

  if      ( method == "bandpass" ) a = spindle_bandpass( edf , param );
  else if ( method == "wavelet" ) a = spindle_wavelet( edf , param );
  else Helper::halt( "SPINDLE method not recognized; should be 'bandpass' or 'wavelet'" );

}

// POL : polarity check for EEG N2/N3 

void proc_polarity( edf_t & edf , param_t & param )
{	
  dsptools::polarity( edf , param );
}

// SW || SLOW-WAVES : detect slow waves, do time-locked FFT on rest of signal

void proc_slowwaves( edf_t & edf , param_t & param )
{	

  // find slow-waves
  slow_waves_t sw( edf , param );
  
}


// WRITE : write a new EDF to disk (but not annotation files, they 
// must always be anchored to the 'original' EDF

void proc_write( edf_t & edf , param_t & param )
{
    
  // add 'tag' to new EDF
  std::string filename = edf.filename;
  if ( Helper::file_extension( filename, "edf" ) || 
       Helper::file_extension( filename, "EDF" ) ) 
    filename = filename.substr(0 , filename.size() - 4 );
  
  filename += "-" + param.requires( "edf-tag" ) + ".edf";

  // optionally, allow directory change
  if ( param.has( "edf-dir" ) )
    {
      const std::string outdir = param.value("edf-dir");
      
      if ( outdir[ outdir.size() - 1 ] != globals::folder_delimiter ) 
	Helper::halt("edf-dir value must end in '/' to specify a folder" );

      int p=filename.size()-1;
      int v = 0;
      for (int j=p;j>=0;j--)
	{
	  if ( filename[j] == globals::folder_delimiter ) { v=j+1; break; }
	}
      filename = outdir + filename.substr( v );            

      // create folder if it does not exist
      std::string syscmd = "mkdir -p " + param.value( "edf-dir" );
      int retval = system( syscmd.c_str() );

      if ( param.has("sample-list") )
	{	  
	  std::string file = param.value("sample-list");

	  // open/append
	  if ( ! globals::silent ) 
	    std::cerr << " appending " << filename << " to sample-list " << file << "\n";
	  
	  std::ofstream FL( file.c_str() , std::ios_base::app );
	  FL << edf.id << "\t"
	     << filename << "\n";
	  
	  FL.close();
	}
    }

  //
  // prep EDF for writing then write to disk
  //

  // arbitrary, but need epochs if not set it seems
  if ( !edf.timeline.epoched() ) 
    edf.timeline.set_epoch( 30 , 30 );

  // if a mask has been set, this will restructure the mask
  edf.restructure(); 

  bool saved = edf.write( filename );

  if ( saved ) 
    if ( ! globals::silent ) 
      std::cerr << " saved new EDF, " << filename << "\n";

}


// EPOCH : set epoch 

void proc_epoch( edf_t & edf , param_t & param )
{
  double dur = 0 , inc = 0;

  // default = 30 seconds, non-overlapping
  if ( ! param.has( "epoch" ) )
    {
      dur = 30; inc = 30;
    }
  else
    {
      std::string p = param.requires( "epoch" );
      std::vector<std::string> tok = Helper::parse( p , "," );
      if ( tok.size() > 2 || tok.size() < 1 ) Helper::halt( "expcting epoch=length{,increment}" );
      
      if ( ! Helper::str2dbl( tok[0] , &dur ) ) Helper::halt( "invalid epoch length" );
      if ( tok.size() == 2 ) 
	{
	  if ( ! Helper::str2dbl( tok[1] , &inc ) ) 
	    Helper::halt( "invalid epoch increment" );
	}
      else inc = dur;
    }

  int ne = edf.timeline.set_epoch( dur , inc );  

  if ( ! globals::silent ) 
    std::cerr << " set epochs, length " << dur << " (step " << inc << "), " << ne << " epochs\n";
  
  if ( param.has("require") )
    {
      int r = param.requires_int( "require" );
      if ( ne < r ) 
	{
	  if ( ! globals::silent ) 
	    std::cout << "EPOCH-PROBLEM\t"
		      << edf.id << "\t"
		      << "[" << globals::current_tag << "]\t"
		      << "required=" << r << "\t"
		      << "observed=" << ne << "\n";
	  globals::problem = true;
	}
    }
}


// FILE-MASK : apply a mask from a file

void proc_file_mask( edf_t & edf , param_t & param )
{ 
  std::string f = "";
  bool exclude = true;

  if      ( param.has("include") ) { f = param.requires("include"); exclude = false; }
  else if ( param.has("exclude") ) f = param.requires("exclude"); 
  else Helper::halt( "need either include or exclude for MASK-FILE");
  
  if ( param.has( "intervals" ) )
    edf.timeline.load_interval_list_mask( f , exclude );
  else
    edf.timeline.load_mask( f , exclude );
}

// EPOCH-MASK  : based on epoch-annotations, apply mask

void proc_epoch_mask( edf_t & edf , param_t & param )
{
  std::set<std::string> vars;
  std::string onelabel;
  
  if ( param.has( "if" ) ) 
    {    
      if ( param.has( "ifnot" ) ) Helper::halt( "both if & ifnot specified" );
      vars = param.strset( "if" );
      onelabel = param.value("if");
      if ( ! globals::silent ) 
	std::cerr << " masking epochs that match " << onelabel << "\n";
    }
  else if ( param.has( "ifnot" ) ) 
    {
      vars = param.strset( "ifnot" );
      onelabel = param.value("ifnot");
      if ( ! globals::silent ) 
	std::cerr << " masking epochs that do not match " << onelabel << "\n";
    }
  else
    Helper::halt( "no if/ifnot specified" );
 
  edf.timeline.apply_simple_epoch_mask( vars , onelabel , param.has("if") );  

}

// EPOCH-ANNOT : directly apply epoch-level annotations

void proc_file_annot( edf_t & edf , param_t & param )
{ 
  std::string f = param.requires( "file" );
  std::vector<std::string> a;
  
  std::map<std::string,std::string> recodes;
  if ( param.has( "recode" ) )
    {
      std::vector<std::string> tok = Helper::quoted_parse( param.value( "recode" ) , "," );
      for (int i=0;i<tok.size();i++)
	{
	  std::vector<std::string> tok2 = Helper::quoted_parse( tok[i] , "=" );
	  if ( tok2.size() == 2 ) 
	    recodes[ Helper::unquote( tok2[1] ) ] = Helper::unquote( tok2[0] );
	  else
	    Helper::halt( "bad format for " + tok[i] );
	}
    }
  
  if ( ! Helper::fileExists( f ) ) Helper::halt( "could not find " + f );

  std::ifstream IN1( f.c_str() , std::ios::in );
  while ( ! IN1.eof() )
    {
      std::string x;
      std::getline( IN1 , x );
      if ( IN1.eof() ) break;
      if ( recodes.find(x) != recodes.end() ) 
	{
	  x = recodes[x];
	}
      a.push_back( x );      
    }
  IN1.close();
  if ( ! globals::silent ) 
    std::cerr << " read " << a.size() << " epoch/annotations from " << f << "\n";
  edf.timeline.annotate_epochs( a );
  
}


// DUMP-MASK : output the current mask as an .annot file

void proc_dump_mask( edf_t & edf , param_t & param )
{
  if ( ! param.has("tag") )
    {
      edf.timeline.dumpmask();
      return;
    }

  // otherwise, create an ANNOT file from the MASK, i.e. for viewing
  std::string tag = param.requires( "tag" );
  std::string path = param.has( "path" ) ? param.value("path") : ".";
  edf.timeline.mask2annot( path, tag );
}


// COUNT/LIST-ANNOTS : show all annotations for the EDF

void proc_list_annots( const std::string & edffile , 
		       const std::string & rootname,
		       const std::vector<std::string> & tok )
{

  const std::set<std::string> inp_signals;
  
  edf_t edf;  

  bool okay = edf.attach( edffile , rootname , &inp_signals );

  if ( ! okay ) 
    if ( ! globals::silent ) 
      std::cerr << "Problem: " << rootname << " reading EDF " << edffile << "\n";

  for (int i=2;i<tok.size();i++) 
    edf.populate_alist( tok[i] );	 

  writer.var( "ANNOT_N" , "Number of occurrences of an annotation" );
  
  std::set<std::string> anames = edf.available_annotations();
  std::set<std::string>::const_iterator ii = anames.begin();
  while ( ii != anames.end() ) 
    {
      // annot as 'level'
      writer.level( *ii , globals::annot_strat );
      writer.value( "ANNOT_N" , edf.aoccur[*ii] );
      ++ii;
    }
}


// LIST-ANNOTS : list all annotations

void proc_list_all_annots( edf_t & edf , param_t & param )
{
  edf.timeline.list_all_annotations( param );
}

// TIME-TRACK : make EDF+

void proc_timetrack( edf_t & edf , param_t & param )
{
  edf.add_continuous_time_track();
}

// RESTRUCTURE : flush masked records

void proc_restructure( edf_t & edf , param_t & param )
{
  // just drop MASK'ed records, then reset mask
  edf.restructure( );
}


// DUMP-RECORDS : show all records (i.e. raw data)

void proc_record_dump( edf_t & edf , param_t & param )
{
  edf.add_continuous_time_track();  
  edf.record_dumper();
}

// STAGE : set sleep stage labels

void proc_sleep_stage( edf_t & edf , param_t & param )
{

//   std::string wake   = "W";
//   std::string nrem1  = "N1";
//   std::string nrem2  = "N2";
//   std::string nrem3  = "N3";
//   std::string nrem4  = "N4";
//   std::string rem    = "R";
//   std::string misc   = "?";
  
  std::string wake   = "";
  std::string nrem1  = "";
  std::string nrem2  = "";
  std::string nrem3  = "";
  std::string nrem4  = "";
  std::string rem    = "";
  std::string misc   = "";
  
  if ( param.has( "W" ) ) wake = param.value("W");
  if ( param.has( "N1" ) ) nrem1 = param.value("N1");
  if ( param.has( "N2" ) ) nrem2 = param.value("N2");
  if ( param.has( "N3" ) ) nrem3 = param.value("N3");
  if ( param.has( "N4" ) ) nrem4 = param.value("N4");
  if ( param.has( "R" ) ) rem = param.value("R");
  if ( param.has( "?" ) ) misc = param.value("?");

  
  if ( param.has( "file" ) )
    {
      std::vector<std::string> ss = Helper::file2strvector( param.value( "file" ) );
      edf.timeline.hypnogram.construct( &edf.timeline , ss );
    }
  else
    {      
      edf.timeline.annotations.make_sleep_stage( wake , nrem1 , nrem2 , nrem3 , nrem4 , rem , misc );
      edf.timeline.hypnogram.construct( &edf.timeline );
    }

  // and output...
  edf.timeline.hypnogram.output();

}


// ED : compute 'electrical distance' measure of bridging

void proc_elec_distance( edf_t & edf , param_t & param )
{
  dsptools::elec_distance( edf , param );
}


// L1OUT : leave-one-out validation via interpolation of all signals
void proc_leave_one_out( edf_t & edf , param_t & param )
{
  dsptools::leave_one_out( edf , param );
}

// COH : calculate cross spectral coherence, using legacy code

void proc_coh_legacy( edf_t & edf , param_t & param )
{
  dsptools::coherence( edf , param , true );
}

// EMD : Empirical Mode Decomposition 
void proc_emd( edf_t & edf , param_t & param )
{
  dsptools::emd_wrapper( edf , param );
}

// ICA : fastICA on sample by channel matrix (whole trace)

void proc_ica( edf_t & edf , param_t & param )
{
  ica_wrapper( edf , param );
}

// COH : calculate cross spectral coherence, using new/default code

void proc_coh( edf_t & edf , param_t & param )
{
  dsptools::coherence( edf , param );
}

// CORREL : correlation

void proc_correl( edf_t & edf , param_t & param )
{
  dsptools::correlate_channels( edf , param );
}

// MI : mutual information

void proc_mi( edf_t & edf , param_t & param )
{
  dsptools::compute_mi( edf , param );
}

// SPIKE : spike in a new bit of signal

void proc_spike( edf_t & edf , param_t & param )
{

  // create a new signal?
  std::string ns = "";
  
  if ( param.has( "new" ) ) ns = param.value( "new" );
  
  signal_list_t from_signal = edf.header.signal_list( param.requires( "from" ) );  
  signal_list_t to_signal   = edf.header.signal_list( param.requires( "to" ) );  
  
  if ( from_signal.size() != 1 ) Helper::halt( "no from={signal}" );
  if ( to_signal.size() != 1 ) Helper::halt( "no to={signal}" );
  
  int s1 = to_signal(0);
  int s2 = from_signal(0);
  
  double wgt = param.requires_dbl( "wgt" );

  spike_signal( edf , s1 , s2 , wgt , ns );
}

// PAC : phase amplitude coupling

void proc_pac( edf_t & edf , param_t & param )
{
  dsptools::pac( edf , param );
}

// CFC : generic cross-frequency coupling methods (other than PAC as above)

void proc_cfc( edf_t & edf , param_t & param )
{
  dsptools::cfc( edf , param );
}

// SUPPRESS-ECG : ECG supression

void proc_ecgsuppression( edf_t & edf , param_t & param )
{
  dsptools::ecgsuppression( edf , param );
}

void proc_bpm( edf_t & edf , param_t & param )
{
  dsptools::bpm( edf , param );
}

// DROP : drop a signal

void proc_drop_signals( edf_t & edf , param_t & param )
{
  
  std::set<std::string> keeps, drops;
  if ( param.has( "keep" ) ) keeps = param.strset( "keep" );
  if ( param.has( "drop" ) ) drops = param.strset( "drop" );
  
  if ( param.has( "keep" ) && param.has( "drop" ) )
    Helper::halt( "can only specify keep or drop with SIGNALS" );
  
  if ( ! ( param.has( "keep" ) || param.has( "drop" ) ) ) 
    Helper::halt( "need to specify keep or drop with SIGNALS" );

  // if a keep list is specified, means we keep 
  if ( keeps.size() > 0 )
    {

      const int ns = edf.header.ns;

      for (int s = 0 ; s < ns ; s++ )
	{
	  std::string label = edf.header.label[s];

	  // is this signal on the keep list?
	  if ( keeps.find( label ) == keeps.end() )
	    {
	      
	      // or, does this signal have an alias that is on the keep list?
	      if ( cmd_t::label_aliases.find( label ) != cmd_t::label_aliases.end() )
		{
		  //std::cout << " has alias " << cmd_t::label_aliases[ label ]  << "\n";
		  if ( keeps.find( cmd_t::label_aliases[ label ] ) == keeps.end() )
		    {
		      //std::cout << "drps " << label << "\n";
		      drops.insert( label );
		      // OR ?? drops.insert( cmd_t::label_aliases[ label ] ) ; should be equiv.
		    }
		}
	      else
		drops.insert( label );
	    }
	  
	} // next signal
    }

  
  std::set<std::string>::const_iterator dd = drops.begin();
  while ( dd != drops.end() )
    {
      if ( edf.header.has_signal( *dd ) )
	{	  	  
	  int s = edf.header.signal( *dd );
	  //std::cerr << "  dropping " << *dd << "\n";
	  edf.drop_signal( s );	  
	}
	++dd;
    }
  
}


// SLICE : pull out slices, based on 'file'

void proc_slice( edf_t & edf , param_t & param , int extract )
{
  
  // x = 1 (extract)
  // x = 0 (exclude)
  
  std::string filename = Helper::expand( param.requires( "file" ) );
  
  std::set<interval_t> intervals;

  if ( ! Helper::fileExists( filename ) ) Helper::halt( "could not find " + filename );


  std::ifstream IN1( filename.c_str() , std::ios::in );  
  while ( ! IN1.eof() )
    {
      interval_t interval;
      IN1 >> interval.start >> interval.stop;
      if ( IN1.eof() ) break;
      if ( interval.stop <= interval.start ) Helper::halt( "problem with interval line" );
      intervals.insert(interval);
    }
  IN1.close();

  if ( ! globals::silent ) 
    std::cerr << " read " << intervals.size() << " from " << filename << "\n";
  
  edf.slicer( intervals , param , extract );
  
}


// Reference tracks

void proc_reference( edf_t & edf , param_t & param )
{
  std::string refstr = param.requires( "ref" );
  signal_list_t references = edf.header.signal_list( refstr );
  
  std::string sigstr = param.requires( "signal" );
  signal_list_t signals = edf.header.signal_list( sigstr );

  edf.reference( signals , references );
  
}

// change record size for one or more signals

void proc_rerecord( edf_t & edf , param_t & param )
{
  double rs = param.requires_dbl( "dur" ); 

  if ( ! globals::silent ) 
    std::cerr << " altering record size from " << edf.header.record_duration << " to " <<  rs << " seconds\n";
  
  edf.reset_record_size( rs );

  if ( ! globals::silent ) 
    std::cerr << " now WRITE'ing EDF to disk, and will set 'problem' flag to skip to next EDF\n";

  proc_write( edf , param );
  globals::problem = true;
}

// uV or mV : set units for tracks

void proc_scale( edf_t & edf , param_t & param , const std::string & sc )
{
  std::string sigstr = param.requires( "signal" );
  signal_list_t signals = edf.header.signal_list( sigstr );
  const int ns = signals.size();  
  for (int s=0;s<ns;s++) 
    edf.rescale( signals(s) , sc );
}


// FLIP : change polarity of signal

void proc_flip( edf_t & edf , param_t & param  )
{
  std::string sigstr = param.requires( "signal" );
  signal_list_t signals = edf.header.signal_list( sigstr );
  const int ns = signals.size();  
  for (int s=0;s<ns;s++) 
    edf.flip( signals(s) );
}





//
// Helper functions
//
	      
void attach_annot( edf_t & edf , const std::string & astr )
{
  
  // annotation is either 
  
  // 1) in memory (and may or may not be in a file,
  // i.e. just created)
  
  annot_t * a = edf.timeline.annotations( astr );
	      
  // 2) or, not in memory but can be retrieved from a file
  
  if ( a == NULL )
    {
      
      // search for file for any 'required' annotation
      // (e.g. annot=stages,spindles)
      
      std::string annot_file = edf.annotation_file( astr );
      
      if ( annot_file == "" ) 
	{
	  if ( ! globals::silent ) 
	    std::cerr << " no instances of annotation [" 
		      << astr << "] for " << edf.id << "\n";		  
	}
      else
	{
	  // add to list and load in data (note, if XML, load() does nothing
	  // as all XML annotations are always added earlier
	  
	  a = edf.timeline.annotations.add( astr );
	  a->load( annot_file );
	}
    }
}

