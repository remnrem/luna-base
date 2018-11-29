
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

#include "stats/cluster.h"

#include "dsp/fir.h"
#include "defs/defs.h"
#include "edf/edf.h"
#include "tinyxml/tinyxml.h"
#include "helper/xml-parser.h"
#include "timeline/timeline.h"
#include "fftw/fftwrap.h"
#include "stats/glm.h"
#include "dsp/dsp.h"
#include "ica/ica.h"
#include "spindles/spindles.h"
#include "spindles/spectral.h"
#include "annot/annot.h"
#include "miscmath/miscmath.h"
#include "artifacts/artifacts.h"
#include "cwt/cwt.h"
#include "pdc/pdc.h"
#include "staging/staging.h"
#include "tinyxml/xmlreader.h"
#include "clocs/topo.h"

#include "main.h"

#include "db/db.h"
#include "dsp/emd.h"
#include "dsp/cfc.h"
#include "dsp/mi.h"

#include "graphics/graphics.h"

#include <thread>
#include <iostream>
#include <unistd.h>
#include <dirent.h>


globals global;

writer_t writer;

std::set<std::string>              cmd_t::commands;
std::string                        cmd_t::input = "";
std::string                        cmd_t::cmdline_cmds = "";
std::string                        cmd_t::stout_file = "";
bool                               cmd_t::append_stout_file = false;
std::map<std::string,std::string>  cmd_t::vars;
std::set<std::string>              cmd_t::signallist;
std::map<std::string,std::string>  cmd_t::label_aliases;
std::map<std::string,std::vector<std::string> >  cmd_t::primary_alias;

#include "dsp/reduce.h"

int main(int argc , char ** argv )
{
  
  //
  // initiate global defintions
  //

  global.init_defs();
  

  //
  // special command-line driven functions that do not involve iterating through a sample list
  //

  bool cmdline_proc_fir_design = false;
  bool cmdline_proc_cwt_design = false;
  bool cmdline_proc_pdlib      = false;

  //
  // parse command line
  //
  
  if ( argc == 2 && strcmp( argv[1] ,"-d" ) == 0 )
    { 
      proc_dummy(); 
      exit(0); 
    } 
  else if ( argc == 3 && strcmp( argv[1] , "--xml" ) == 0 )
    {
      annot_t::dumpxml( argv[2] , false );
      std::exit(0);
    }
  else if ( argc == 3 && strcmp( argv[1] , "--xml-dump" ) == 0 )
    {
      annot_t::dumpxml( argv[2] , true );
      std::exit(0);
    }  
  else if ( argc == 2 && strcmp( argv[1] , "--pdlib" ) == 0 ) 
    {
      param_t param = build_param_from_cmdline();
      writer.nodb();
      writer.begin();      
      writer.id( "." , "." );
      pdc_t::construct_pdlib( param );
      writer.commit();
      std::exit(0);
    }
  else if ( argc >= 2 )
    {  
      
      // pick off any special functions here:
      // i.e. the first element will not be interpreted 
      // as a file list, as we will run a cmdline proc and 
      // then quit

      if ( argc >= 2 )
	{
	  if ( strcmp( argv[1] , "--fir-design" ) == 0 ||  
	       strcmp( argv[1] , "--fir" ) == 0 ) 	   
	    cmdline_proc_fir_design = true;
	  else 
	    if ( strcmp( argv[1] , "--cwt-design" ) == 0 ||  
		 strcmp( argv[1] , "--cwt" ) == 0 ) 	   
	      cmdline_proc_cwt_design = true;
	}

      // otherwise, first element will be treated as a file list
      
      cmd_t::input = argv[1];

      // commands should be in form
      
      //   @variable-file 
      //   key=value
      //   signals=S1,S2,...
      //   --flag
      //   exclude={file}
      //   -o  output db
      //   -a  output db  { as above, except append } 
      //   -s  { rest of line is script }


      // comma-separated strings (-->signals)
      // a single number or range (n, n-m)
      // var=value

      int specified = 0;

      for (int i=2;i<argc;i++)
	{
	  std::vector<std::string> tok = 
	    Helper::quoted_parse( argv[i] , "=" );
	  
	  if ( tok.size() == 2 ) 
	    {
	      // add signals?
	      if ( Helper::iequals( tok[0] , "signal" ) 
		   || Helper::iequals( tok[0] , "signals" ) )
		{		  
		  
		  std::vector<std::string> tok2 
		    = Helper::quoted_parse( tok[1] , "," );		  
		  
		  for (int s=0;s<tok2.size();s++) 
		    cmd_t::signallist.insert(Helper::unquote(tok2[s]));		  
		}

	      // do not read EDF annotations
	      else if ( Helper::iequals( tok[0] , "force-edf" ) )
		{
		  if ( tok[1] == "1" || tok[1] == "Y" || tok[1] == "y" )
		    globals::skip_edf_annots = true;
		}
	      
	      // signal alias?
	      else if ( Helper::iequals( tok[0] , "alias" ) )
		{
		  cmd_t::signal_alias( tok[1] );
		}
	      
	      // exclude individuals?
	      else if ( Helper::iequals( tok[0] , "exclude" ) )
		{

		  std::string xfile = Helper::expand( tok[1] );

		  if ( Helper::fileExists( xfile ) ) 
		    {
		      std::ifstream XIN( xfile.c_str() , std::ios::in );
		      while ( !XIN.eof() ) 
			{
			  // format: ID {white-space} any notes (ignored)
			  std::string line2;
			  std::getline( XIN , line2);
			  if ( XIN.eof() || line2 == "" ) continue;
			  std::vector<std::string> tok = Helper::parse( line2 , "\t " );
			  if ( tok.size() == 0 ) continue;			      
			  std::string xid = tok[0];
			  globals::excludes.insert( xid );
			}
		      std::cerr << "excluding " << globals::excludes.size() 
				<< " individuals from " << xfile << "\n";
		      XIN.close();
		    }
		  else std::cerr << "**warning: exclude file " << xfile << " does not exist\n";
		}
	      else if ( tok[0][0] == '-' )
		{
		  globals::param.add( tok[0].substr(1) , tok[1] );
		}
	      else // a standard variable
		{		  
		  cmd_t::vars[ tok[0] ] = tok[1];
		}
	    }
	  else if ( Helper::iequals( tok[0] , "-o" ) || Helper::iequals( tok[0] , "-a" ) )
	    {
	      // next arg will be DB
	      if ( i + 1 >= argc ) Helper::halt( "expecting database name after -o" );
	      cmd_t::stout_file = argv[ ++i ];
	      if ( Helper::iequals( tok[0] , "-a" ) ) cmd_t::append_stout_file = true;
	    }
	  else if ( Helper::iequals( tok[0] , "-s" ) )
	    {
	      // rest of args will be cmd script
	      for (int j=i+1;j<argc;j++) cmd_t::add_cmdline_cmd( argv[j] );
	      break;		
	    }
	  else if ( argv[i][0] == '-' )
	    {
	      std::string f = argv[i];	      
	      globals::param.add( f.substr(1) );
	    }
	  else if ( argv[i][0] == '@' )
	    {
	      // an 'include' 
	      std::string filename = argv[i];
	      // expand() expands out any ~/ notation to full path
	      filename = Helper::expand( filename.substr(1).c_str() );
	      if ( ! Helper::fileExists( filename ) ) Helper::halt( "could not open " + filename );
	      std::ifstream INC( filename , std::ios::in );
	      if ( INC.bad() ) Helper::halt("could not open file: " + filename );
	      while ( ! INC.eof() )
		{
		  std::string line;
		  std::getline( INC , line);
		  if ( INC.eof() || line == "" ) continue;
		  std::vector<std::string> tok = Helper::quoted_parse( line , "\t" );
		  if ( tok.size() != 2 ) Helper::halt("badly formatted line in " + filename );
		  
		  // a signal?
		  if ( Helper::iequals( tok[0] , "signal" ) || Helper::iequals( tok[0] , "signals" ) )
		    {		  
		      std::vector<std::string> tok2 = Helper::quoted_parse( tok[1] , "," );		  
		      for (int s=0;s<tok2.size();s++) 
			cmd_t::signallist.insert(Helper::unquote(tok2[s]));					  
		    }
		  
		  // signal alias?
		  else if ( Helper::iequals( tok[0] , "alias" ) )
		    {
		      cmd_t::signal_alias( tok[1] );
		    }

		  // skip EDF Annotations channels
		  else if ( Helper::iequals( tok[0] , "force-edf" ) )
		    {
		      if ( tok[1] == "1" || tok[1] == "Y" || tok[1] == "y" ) 
			globals::skip_edf_annots = true;
		    }
		  
		  // default annot folder?
		  else if ( Helper::iequals( tok[0] , "annots" ) ) 
		    {
		      if ( tok[1][ tok[1].size() - 1 ] != globals::folder_delimiter )
			globals::annot_folder = tok[1] + "/";
		      else
			globals::annot_folder = tok[1];		      
		    }
		  else if ( Helper::iequals( tok[0] , "exclude" ) )
		    {
		      std::string xfile = Helper::expand( tok[1] );
		      if ( Helper::fileExists( xfile ) ) 
			{
			  std::ifstream XIN( xfile.c_str() , std::ios::in );
			  while ( !XIN.eof() ) 
			    {
			      std::string line2;
			      std::getline( XIN , line2);
			      if ( XIN.eof() || line2 == "" ) continue;
			      std::vector<std::string> tok = Helper::parse( line2 , "\t " );
			      if ( tok.size() == 0 ) continue;			      
			      std::string xid = tok[0];
			      globals::excludes.insert( xid );
			    }
			  std::cerr << "excluding " << globals::excludes.size() 
				    << " individuals from " << xfile << "\n";
			  XIN.close();
			}
		      else std::cerr << "** warning: exclude file " << xfile << " does not exist\n";
		    }
		  else // ... or just a variable?
		    {		  
		      cmd_t::vars[ tok[0] ] = tok[1];	  
		    }
		}
	      INC.close();
	    }
	  else
	    {
	      int x;
	      if ( ! Helper::str2int( argv[i] , &x ) )
		{
		  // assume this is an ID (i.e. must be a string)
		  globals::sample_list_id = argv[i];		  
		  specified = 2;// i.e. done selecting
		}
	      else if ( specified == 0 )
		{
		  globals::sample_list_min = x;
		  ++specified;
		}
	      else if ( specified == 1 )
		{
		  globals::sample_list_max = x;
		  ++specified;
		}
	      else
		Helper::halt( "cannot parse command line: more than two sample lines specified" );
	      
	    }
	  
	  // only a single sample specified?
	  if ( specified == 1 ) 
	    {
	      globals::sample_list_max = globals::sample_list_min;
	    }
	  else if ( globals::sample_list_max < globals::sample_list_min )
	    {
	      int x = globals::sample_list_max;
	      globals::sample_list_max = globals::sample_list_min;
	      globals::sample_list_min = x;
	    }
	  
	  if ( globals::sample_list_min < 0 ) globals::sample_list_min = -1;
	  if ( globals::sample_list_max < 0 ) globals::sample_list_max = -1;
	    
	}
            
    }
  
  else if ( argc < 2 || ( isatty(STDIN_FILENO) || argc != 1 ) )  
    {
      std::cerr << "luna version " << globals::version << "; " << globals::date 
		<< "\nusage: luna [sample-list|EDF] [n1] [n2] [signal=s1,s2] [v1=val1] [@parameter-file] < command-file\n" ;
      std::exit(1);
    }

  if ( std::cin.eof() || ! std::cin.good() ) 
    Helper::halt( "no input, quitting" );



  //
  // initialize output to a STOUT db or not?
  //

  if ( ! cmd_t::append_stout_file ) 
    Helper::deleteFile( cmd_t::stout_file );
  
  if ( cmd_t::stout_file != "" )
    writer.attach( cmd_t::stout_file );  
  else 
    writer.nodb();
  

  //
  // specify types for common stratifiers
  //

  writer.numeric_factor( globals::epoch_strat );
  writer.numeric_factor( globals::freq_strat );
  writer.numeric_factor( globals::cycle_strat );
  writer.string_factor( globals::band_strat );
  writer.string_factor( globals::signal_strat );
  writer.string_factor( globals::stage_strat );
  writer.numeric_factor( globals::count_strat );
  writer.numeric_factor( globals::time_strat );

  
  //
  // branch off to run any cmdline driven special functions, then quit
  //

  if ( cmdline_proc_fir_design )
    {

      writer.begin();      
      writer.id( "." , "." );
       
      // expects input from std::cin
      proc_filter_design_cmdline();

      writer.commit();

      std::exit(0);
    }


  if ( cmdline_proc_cwt_design )
    {
      writer.begin();      
      writer.id( "." , "." );
       
      // expects input from std::cin
      proc_cwt_design_cmdline();

      writer.commit();

      std::exit(0);

    }


  //
  // iterate through the primary sample-list
  //
  
  int processed = 0, failed = 0;
  int actually_processed = 0;

  while ( ! std::cin.eof() )
    {
      cmd_t cmd; // scans STDIN for next command
      
      if ( cmd.empty() ) break; 
      
      ++processed;

      if ( ! cmd.valid() )
	++failed;
      else
	{
	  
	  
	  // process command ( most will iterate over 1 or more EDFs)
	  if ( cmd.process_edfs() ) 
	    {
	      process_edfs(cmd);
	    }
	  else // handle any exceptions 
	    {

	      // i.e. commands where we do not simply iterate over the
	      // EDF filelist; note, these currently are only
	      // applicable in the single command mode;

	      if ( cmd.is( 0 , "INTERVALS" ) ) 
		proc_intervals( cmd.param(0) , cmd.data() );
	      
	    }
	  
	}


      // if received command from the -s option, we are all done
      if ( cmd_t::cmdline_cmds != "" ) break;
      
    }
  
  //
  // wrap up
  //

  std::cerr << "...processed " << processed << " command(s), ";
  if ( failed == 0 ) std::cerr << " all of which passed\n";
  else std::cerr << failed << " of which failed\n";

  exit(0);
  
}




void process_edfs( cmd_t & cmd )
{
  
  //
  // Iterate over some or all of a list of EDFs and annotations,
  // performing one or more commands
  //

  if ( cmd.num_cmds() == 0 ) return;  // nothing to do
  
  if ( ! Helper::fileExists( cmd.data() ) ) 
    Helper::halt( "could not find file list, " + cmd.data() );

  
  //
  // Open sample-list
  //

  std::string f = cmd.data();
  f = f.substr( f.size() - 4 );
  bool single_edf = f == ".edf" || f == ".EDF";

  std::ifstream EDFLIST;
  
  if ( ! single_edf ) 
    EDFLIST.open( cmd.data().c_str() , std::ios::in );


  //
  // Start iterating through it
  //
  
  int processed = 0;
  int actual = 0;

  while ( single_edf || ! EDFLIST.eof() )
    {
            
      // each line should contain (tab-delimited)  
      //  1  ID
      //  2  EDF file
      //  3+ other optional ANNOT files for that EDF
      
      std::string rootname;
      std::string edffile;
      std::vector<std::string> tok;
      
      if ( ! single_edf )
	{
	  std::string line;
	  std::getline( EDFLIST , line);
	  
	  if ( line == "" ) continue;

	  //
	  // If we are only looking at a subset of the sample list, 
	  // might skip here	  
	  //
	  
	  if ( globals::sample_list_min != -1 || globals::sample_list_max != -1 )
	    {
	      int line_n = processed + 1;
	      if ( line_n < globals::sample_list_min || 
		   line_n > globals::sample_list_max ) 
		{
		  ++processed;
		  continue;
		}
	    }
	  
	  // parse by tabs
	  
	  tok = Helper::parse( line , "\t" );      
	  if ( tok.size() < 2 ) 
	    Helper::halt( "requires (ID) | EDF file | (optional ANNOT files)" );
	  
	  rootname = tok[0];
	  edffile  = tok[1];
	  
	  // else, do we have an 'ID' check?
	  if ( globals::sample_list_id != "" )
	    {	      
	      if ( rootname != globals::sample_list_id )
		{
		  ++processed;
		  continue;
		}
	    }
	}
      else 
	{
	  edffile = cmd.data();
	  rootname = edffile;
	}

      //
      // File in exclude list?
      //
          
      if ( globals::excludes.find( rootname ) != globals::excludes.end() )
	{
	  std::cerr << "\nskipping EDF " << rootname << "\n";
	  ++processed;
	  continue; // to the next EDF in the list
	}
      

      //
      // Begin running through the series of commands
      //

      
      std::cerr << "\nProcessing EDF " << rootname 
		<< " [ #" << processed+1 << " ]\n";
      
      
      //
      // Begin transaction
      //

      writer.begin();

      
      writer.id( rootname , edffile );
      

      //
      // Update any indiv-wildcards in the command list
      //

      cmd.replace_wildcards( rootname );

      
      //
      // Handle 'EDF-validation' commands separately
      //
      // That is for some commands, we do not want to attach the EDF
      // yet; i.e. note, this means that any other commands on the
      // same line will be ignored; these only operate on the EDFs, so
      // and on all signals, so no need to attach annotations, or
      // process signals to extract, etc
      //
      
      bool validation_only = false;
      
      for (int c=0;c<cmd.num_cmds();c++)
	{
	  
	  if ( cmd.is( c,"HEADERS" ) )   
	    { 
	      writer.cmd( cmd.cmd(c) , c+1 , cmd.param(c).dump( "" , " " ) );
	      validation_only = true; 
	      proc_summaries( edffile , rootname , 1 , cmd , cmd.param(c) ); 
	    } 
	  else if ( cmd.is( c, "COUNT-ANNOTS" ) )
	    {
	      writer.cmd( cmd.cmd(c) , c+1 , cmd.param(c).dump( "" , " " ) );
	      validation_only = true; 
	      proc_list_annots( edffile , rootname , tok );
	    }
	}
      
      // skip to next EDF if any commands were 'validation' commands
      
      if ( validation_only ) 
	{ 
	  ++processed; 
	  ++actual; 	  
	  writer.commit();
	  if ( single_edf ) break;
	  std::cerr << "  skipping any other commands\n";	  
	  continue; 
	}
      
      

      
      //
      // Otherwise, assume one or more full commands
      //

      //
      // Unset 'problem' flag (i.e. for bailing for this individual)
      //
      
      globals::problem = false;


      //
      // Limited to specific signals to load in?
      //
      
      const std::set<std::string> * inp_signals = NULL;
      
      if ( cmd.signals().size() > 0 ) inp_signals = &cmd.signals();
      
      //
      // load EDF
      //

      edf_t edf;

      bool okay = edf.attach( edffile , rootname , inp_signals );
      
      if ( ! okay ) 
	{
	  globals::problem = true;
	  std::cerr << "**warning: problem loading " 
		    << edffile << ", skipping...\n";

	  writer.commit();
	  continue;
	}

      
    
      //
      // Attach annotations (not loaded, just store list of available ones)
      //
      
      if ( ! single_edf ) 
	{
	  for (int i=2;i<tok.size();i++) 
	    {
	      if ( tok[i][ tok[i].size() - 1 ] == '/' ) 
		{
		  // this means we are specifying a folder, in which case search for all files that 
		  // start id_<ID>_* and attach thoses
		  DIR * dir;		  
		  struct dirent *ent;
		  if ( (dir = opendir ( tok[i].c_str() ) ) != NULL )
		    {
		      /* print all the files and directories within directory */
		      while ((ent = readdir (dir)) != NULL)
			{
			  std::string fname = ent->d_name;
			  
			  // only feature lists
 			  if ( Helper::file_extension( fname , "ftr" ) ||
 			       Helper::file_extension( fname , "xml" ) )
			    edf.populate_alist( tok[i] + fname );	 			   
			}
		      closedir (dir);
		    }
		  else 
		    Helper::halt( "could not open folder " + tok[i] );
		}
	      else
		{
		  edf.populate_alist( tok[i] );	 
		}
	      
	    }
	}
      
	    
      //
      // Loop over each command
      //

      for ( int c = 0 ; c < cmd.num_cmds() ; c++ )
	{	        
	  
	  // was a problem flag raised when loading the EDF?
	  
	  if ( globals::problem ) break;
	  
	  //
	  // If this particular command did not explicitly specify
	  // signals, then add all
	  //
	  
	  if ( ! cmd.param(c).has( "signal" ) )
	    cmd.param(c).add( "signal" , cmd.signal_string() );
	  
	  //
	  // Attach any new annotations (e.g. sleep stages, etc) to
	  // the timeline if requested; these will persist as long as
	  // the EDF lasts
	  //
	  
	  std::set<std::string> required_annotations = cmd.param(c).strset( "annot" , "," );
	  std::set<std::string>::const_iterator ai = required_annotations.begin();	  

	  while ( ai != required_annotations.end() )
	    {
	      attach_annot( edf , *ai );
	      ++ai;
	    }
	  
	  
	  //
	  // Print command
	  //
	  
	  std::cerr << " CMD #" << c+1 << ": " 
		    << cmd.cmd(c) << "\n";
	  
	  
	  writer.cmd( cmd.cmd(c) , c+1 , cmd.param(c).dump( "" , " " ) );
	  

	  //
	  // Now process the command
	  //
	  
	  if      ( cmd.is( c, "WRITE" ) )        proc_write( edf, cmd.param(c) );
	  else if ( cmd.is( c, "SUMMARY" ) )      proc_summaries( edf , cmd.param(c) );
	  else if ( cmd.is( c, "DESC" ) )         proc_desc( edf , cmd.param(c) );
	  else if ( cmd.is( c, "STATS" ) )        proc_validate( edf , cmd.param(c) );
	  
	  else if ( cmd.is( c, "REFERENCE" ) )    proc_reference( edf , cmd.param(c) );
	  else if ( cmd.is( c, "FLIP" ) )         proc_flip( edf , cmd.param(c) );
	  else if ( cmd.is( c, "uV" ) )           proc_scale( edf , cmd.param(c) , "uV" ); 
	  else if ( cmd.is( c, "mV" ) )           proc_scale( edf , cmd.param(c) , "mV" );
	  else if ( cmd.is( c, "RECORD-SIZE" ) )  proc_rerecord( edf , cmd.param(c) );

	  else if ( cmd.is( c, "TIME-TRACK" ) )   proc_timetrack( edf, cmd.param(c) );
	  
	  else if ( cmd.is( c, "STAGE" ) || cmd.is( c, "HYPNO" ) ) proc_sleep_stage( edf , cmd.param(c) );
	
	  else if ( cmd.is( c, "TSLIB" ) )        pdc_t::construct_tslib( edf , cmd.param(c) );
	  else if ( cmd.is( c, "SSS" ) )          pdc_t::simple_sleep_scorer( edf , cmd.param(c) );
	  else if ( cmd.is( c, "EXE" ) )          pdc_t::similarity_matrix( edf , cmd.param(c) );
	  
	  else if ( cmd.is( c, "DUMP" ) )         proc_dump( edf, cmd.param(c) );
	  else if ( cmd.is( c, "DUMP-RECORDS" ) ) proc_record_dump( edf , cmd.param(c) );
	  else if ( cmd.is( c, "DUMP-EPOCHS" ) )  proc_epoch_dump( edf, cmd.param(c) );
	  else if ( cmd.is( c, "ANNOTS" ) )       proc_list_all_annots( edf, cmd.param(c) );
	  else if ( cmd.is( c, "MATRIX" ) )       proc_epoch_matrix( edf , cmd.param(c) );
	  else if ( cmd.is( c, "RESTRUCTURE" ) )  proc_restructure( edf , cmd.param(c) );
	  else if ( cmd.is( c, "SIGNALS" ) )      proc_drop_signals( edf , cmd.param(c) );
	  else if ( cmd.is( c, "RMS" ) || cmd.is( c, "SIGSTATS" ) ) proc_rms( edf, cmd.param(c) );
	  else if ( cmd.is( c, "MSE" ) )          proc_mse( edf, cmd.param(c) );
	  else if ( cmd.is( c, "LZW" ) )          proc_lzw( edf, cmd.param(c) );
	  else if ( cmd.is( c, "ZR" ) )           proc_zratio( edf , cmd.param(c) );
	  else if ( cmd.is( c, "ANON" ) )         proc_anon( edf , cmd.param(c) );
	  else if ( cmd.is( c ,"EPOCH" ) )        proc_epoch( edf, cmd.param(c) );
	  else if ( cmd.is( c ,"SLICE" ) )        proc_slice( edf , cmd.param(c) , 1 );

	  else if ( cmd.is( c, "MASK" ) )         proc_mask( edf, cmd.param(c) );
	  else if ( cmd.is( c, "FILE-MASK" ) )    proc_file_mask( edf , cmd.param(c) );
	  else if ( cmd.is( c, "DUMP-MASK" ) )    proc_dump_mask( edf, cmd.param(c) );
	  
	  else if ( cmd.is( c, "EPOCH-ANNOT" ) )  proc_file_annot( edf , cmd.param(c) );
	  else if ( cmd.is( c, "EPOCH-MASK" ) )   proc_epoch_mask( edf, cmd.param(c) );
	  
	  
	  else if ( cmd.is( c, "FILTER" ) )       proc_filter( edf, cmd.param(c) );
	  else if ( cmd.is( c, "FILTER-DESIGN" )) proc_filter_design( edf, cmd.param(c) );
	  //	  else if ( cmd.is( c, "LEGACY-FILTER" )) proc_filter_legacy( edf, cmd.param(c) );
	  
	  else if ( cmd.is( c, "COVAR" ) )        proc_covar( edf, cmd.param(c) );
	  else if ( cmd.is( c, "PSD" ) )          proc_psd( edf, cmd.param(c) );	  
	  else if ( cmd.is( c, "1FNORM" ) )       proc_1overf_norm( edf, cmd.param(c) );

	  else if ( cmd.is( c, "FIP" ) )          proc_fiplot( edf , cmd.param(c) );

	  else if ( cmd.is( c, "COH" ) )          proc_coh( edf , cmd.param(c) );
	  else if ( cmd.is( c, "LECGACY-COH" ) )  proc_coh_legacy( edf , cmd.param(c) );
	  else if ( cmd.is( c, "CORREL" ) )       proc_correl( edf , cmd.param(c) );
	  else if ( cmd.is( c, "ED" ) )           proc_elec_distance( edf , cmd.param(c) );
	  else if ( cmd.is( c, "ICA" ) )          proc_ica( edf, cmd.param(c) );
	  else if ( cmd.is( c, "L1OUT" ) )        proc_leave_one_out( edf , cmd.param(c) );

	  else if ( cmd.is( c, "EMD" ) )          proc_emd( edf , cmd.param(c) );

	  else if ( cmd.is( c, "MI" ) )           proc_mi( edf, cmd.param(c) );
	  else if ( cmd.is( c, "HR" ) )           proc_bpm( edf , cmd.param(c) );
	  else if ( cmd.is( c, "SUPPRESS-ECG" ) ) proc_ecgsuppression( edf , cmd.param(c) );
	  else if ( cmd.is( c, "PAC" ) )          proc_pac( edf , cmd.param(c) );
	  else if ( cmd.is( c, "CFC" ) )          proc_cfc( edf , cmd.param(c) );
	  else if ( cmd.is( c, "TAG" ) )          proc_tag( cmd.param(c) );
	  else if ( cmd.is( c, "RESAMPLE" ) )     proc_resample( edf, cmd.param(c) );
	  else if ( cmd.is( c, "SPINDLES" ) )     proc_spindles( edf, cmd.param(c) );	  
	  else if ( cmd.is( c, "POL" ) )          proc_polarity( edf, cmd.param(c) );	  
	  
	  else if ( cmd.is( c, "SW" ) || cmd.is( c, "SLOW-WAVES" ) ) proc_slowwaves( edf, cmd.param(c) );
	  else if ( cmd.is( c, "ARTIFACTS" ) )    proc_artifacts( edf, cmd.param(c) );
	  else if ( cmd.is( c, "SPIKE" ) )        proc_spike( edf , cmd.param(c) );
	  
    	  else Helper::halt( "did not recognize command: " + cmd.cmd(c) );

	  //
	  // Was a problem flag set?
	  //
	  
	  if ( globals::problem ) 
	    {
	      std::cerr << "**warning: the PROBLEM flag was set, skipping to next EDF...\n";
	      break;
	    }

	} // next command
          
      ++processed;
      ++actual;


      //
      // commit output to DB
      //
     
      writer.commit();


      //
      // all done / next EDF
      //

      if ( single_edf ) break;

    }


  //
  // Close sample-list if open
  //
  
  if ( ! single_edf ) 
    EDFLIST.close();
  

  //
  // All done
  //

  std::cerr << "\n...processed " << actual << " EDFs, done.\n";

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
  if ( ! okay ) std::cerr << "Problem: " << rootname << " reading EDF " << edffile << "\n";    
  
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
  if ( t != "." ) std::cerr << " setting analysis tag to [" << globals::current_tag << "]\n";
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
  std::cerr << " set epochs, length " << dur << " (step " << inc << "), " << ne << " epochs\n";
  if ( param.has("require") )
    {
      int r = param.requires_int( "require" );
      if ( ne < r ) 
	{
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
      std::cerr << " masking epochs that match " << onelabel << "\n";
    }
  else if ( param.has( "ifnot" ) ) 
    {
      vars = param.strset( "ifnot" );
      onelabel = param.value("ifnot");
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
  if ( ! okay ) std::cerr << "Problem: " << rootname << " reading EDF " << edffile << "\n";

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
  std::cerr << " altering record size from " << edf.header.record_duration << " to " <<  rs << " seconds\n";
  edf.reset_record_size( rs );
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

bool cmd_t::read()
{
  
  if ( std::cin.eof() ) return false;
  
  reset();
  
  // CMD param=1 p=1,2,3 f=true out=o1 & CMD2 etc ; 
  
  // EITHER read from std::cin, OR from -s command line
  
  std::istringstream allinput;
  
  if ( cmd_t::cmdline_cmds == "" )
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

  return true;
}




// DUMMY : a generic placeholder/scratchpad for templating new things

void proc_dummy()
{

  Data::Matrix<double> D( 10 , 10 );
  for (int i=0;i<10;i++) 
    for (int j=0;j<10;j++) 
      std::cin >> D(i,j) ;
  
  cluster_t cl;
  cl.build( D );

  
  std::exit(1);

  if ( 1 ) 
    {

      topo_t topo;
      
      int ch = topo.load("/Users/shaun/dropbox/ebooks/eegbook/clocs.theta.rad.txt");

      topo.max_radius( 0.55 );
      
      topo.grid( 67 , 67 );
      
      std::map<std::string, double> data;  
      while ( ! std::cin.eof() ) 
	{
	  std::string l;
	  double x,y,z;
	  std::cin >> l;
	  if ( std::cin.eof() ) continue;
	  if ( l == "" ) continue;
	  std::cin >> z ;
	  data[l] = z;
	  
	}
      
      Data::Matrix<double> I = topo.interpolate( data );
      
      std::cout << I.print() << "\n";
      
      std::exit(0);
    }


  if ( 0 ) 
    {
      

      std::vector<double> x,y;
      std::vector<double> z;
      double xmin = 99, xmax = -99;
      double ymin = 99, ymax = -99;
      
      while ( ! std::cin.eof() ) 
	{
	  std::string ch;
	  double x1,y1,z1;
	  std::cin >> ch ;
	  if ( ch == "" ) continue;
	  std::cin >> x1 >> y1 >> z1;
	  x.push_back( x1  );
	  y.push_back( y1 );
	  z.push_back( z.size() );
	  if ( x1 < xmin ) xmin = x1;
	  if ( y1 < ymin ) ymin = y1;
	  
	  if ( x1 > xmax ) xmax = x1;
	  if ( y1 > ymax ) ymax = y1;

	}

      std::cout << "read " << z.size() << "\n";
      
      for (int i=0;i<x.size();i++) 
	{ 
	  x[i] = ( x[i] - xmin ) / ( xmax - xmin ) ; 
	  y[i] = ( y[i] - ymin ) / ( ymax - ymin ) ; 
	}
      

      Data::Matrix<double> results =   
	dsptools::interpolate2D( x , y , z , 0 , 1 , 67 , 0 , 1 , 67 );
      
      std::cout << "\n" << results.print() << "\n";
      
    }


  std::exit(1);
  
  writer.nodb();

  writer.begin();      
  writer.id( "." , "." );
  
  pdc_t pdc;
  pdc.test();

  writer.commit();

  std::exit(1);
  
  std::vector<double> x;
  //  std::vector<double> h;
  std::ifstream X1( "eeg.txt", std::ios::in );
  //std::ifstream H1( "h.txt", std::ios::in );

  while ( ! X1.eof() )
    {
      double xx;
      X1 >> xx;
      if ( X1.eof() ) break;
      x.push_back( xx ) ;
    }
  X1.close();

  // while ( ! H1.eof() )
  //   {
  //     double xx;
  //     H1 >> xx;
  //     if ( H1.eof() ) break;
  //     h.push_back( xx ) ;
  //   }
  // H1.close();

  std::cerr << "read " << x.size() << " x values\n";
  //  std::cerr << "read " << h.size() << " h values\n";

  
  int kaiserWindowLength;
  double beta;

  double transitionWidthHz = 0.5;
  double ripple = 0.001;
  int sampFreq = 128;
  double trans1Freq = 10;
  double trans2Freq = 16;

  fir_t fir;
  fir.calculateKaiserParams( ripple , transitionWidthHz, sampFreq, &kaiserWindowLength, &beta);
  
  std::cout << "KWL, beta = " << kaiserWindowLength << "\t" << beta << "\n";
  if ( kaiserWindowLength % 2 == 0 ) ++kaiserWindowLength;
  std::vector<double> bpf = fir.create2TransSinc(kaiserWindowLength, trans1Freq, trans2Freq,  sampFreq, fir_t::BAND_PASS);
  std::vector<double> fc = fir.createKaiserWindow(&bpf, beta);

  //
  // Convolve
  //
  
  for (int tt=0;tt<10;tt++)
    {
      std::cerr << "tt + " << tt << "\n";

      fir_impl_t ft( fc );      
      std::vector<double> ff1 = ft.filter( &x );
      //std::vector<double> ff2 = ft.fft_filter( &x );      
      //for (int j=0;j<ff1.size();j++) std::cout  << "conv\t" << j << "\t" << ff1[j] << "\n";     
      // for (int j=0;j<ff2.size();j++) std::cout  << "fftconv\t" << j << "\t" << ff2[j] << "\n";
    }
  
  std::exit(1);
    
  //
  // Method 2
  //


  // fir_impl_t ft( fc );
  // std::vector<double> ff2 = ft.filter( &x );

  // for ( int  i = 0 ; i < ff2.size() ; i ++ ) 
  //   std::cout << "ff2_" << i << "\t" << x[i] << "\t" << ff2[i] << "\n";

  // for ( int  i = 0 ; i < c.size() ; i ++ ) 
  //   std::cout << "cc1_" << i << "\t" << c[i] << "\n";
  
  std::exit(1);
  

  std::vector<double> signal( 5 , 1 );
  std::vector<double> kernel( 2 , 1 );
  signal.push_back( 2222 );
  signal.push_back( 2 );
  signal.push_back( 2 );
  signal.push_back( 2 );
  signal.push_back( 2 );

  std::vector<double> conv = dsptools::convolve( signal , kernel );
  
  for (int i = 0 ; i < conv.size() ; i++ ) 
    {
      std::cout << "\t" << conv[i] << "\n";
    }
  
  std::exit(1);
  
  
  
  edf_t edf;
  timeline_t t(&edf);
  t.load_interval_list_mask( "n2.int" );
  std::exit(0);


  //
  // Load signal from STDIN
  //

  //   std::vector<double> x;

//   while ( ! std::cin.eof() ) 
//     {
//       double t;
//       std::cin >> t;
//       if ( std::cin.eof() ) break;
//       x.push_back(t);
//     }
//   std::cerr << "read " << x.size() << " values\n";

   while ( ! std::cin.eof() ) 
     {
       double t;
       std::cin >> t;
       if ( std::cin.eof() ) break;
       x.push_back(t);
     }
   std::cerr << "read " << x.size() << " values\n";


//   std::exit(1);
  
   //
   // CFC/GLM
   //

   cfc_t cfc( x , 0.5 , 1.5 , 9 , 11 , 200 );

   bool okay = cfc.glm();

   if ( ! okay ) Helper::halt( "problem" );

   std::cout << "r_PAC\t" << cfc.r_PAC << "\n"
	     << "C_AMP\t" << cfc.c_AMP << "\n"
	     << "Z_AMP\t" << cfc.z_AMP << "\n"
	     << "R2_TOT\t" << cfc.r2_TOT <<"\n";
   
   std::exit(1);


//   // Test
  
//   int index_length = x.size();
//   int my_Fs = 1000;
//   int index_start = 0;


//   FFT fftseg( index_length , my_Fs , FFT_FORWARD , WINDOW_NONE );
  
//   fftseg.apply( &(x[index_start]) , index_length );
  
//   int my_N = fftseg.cutoff;
//   for (int f=0;f<my_N;f++)
//     {
//       std::cout << f << "\t" << fftseg.frq[f] << "\t" << fftseg.X[f] << "\n";
//     }

//   std::exit(1);

//   //
//   // EMD
//   //
  // 
  // fiplot_t fp( x , 200 , 0 , 5 , 0.1 , 
  // 	       15 , 15 , 1 ); // frequencies


//   double x2 = 0.9243;

//   for (int j=1;j<=10;j++)
//     {
//       std::vector<double> lp = legendre( j , x2 );
//       for (int k=0;k<j;k++) std::cout << " " << j << " " << lp.size() << " " << lp[k] << "\n";
//     }
  
//   std::exit(1);

  //
  //
  
  clocs_t clocs;
  clocs.load_cart( "ex.clocs" );
  ///Users/shaun/dropbox/sleep/projects/grins-test/clocs.txt" );

  std::exit(1);

  // read data : 64 by 
  int ns = 64;
  int np = 63360;
  
  Data::Matrix<double> X( ns , np );

  for (int c=0;c<ns;c++) // channel by data points
    for (int r=0;r<np;r++)
      {
 	double x1;
 	std::cin >> x1;
 	if ( std::cin.eof() ) Helper::halt("prob");
 	X(r,c) = x1;
       }
  
  signal_list_t good_signals;
  signal_list_t bad_signals;

  int si = 0;
  std::ifstream IN1("good.sigs",std::ios::in);
  while ( !IN1.eof() ) 
    {
      std::string l;
      IN1 >> l;
      if ( IN1.eof() ) break;
      good_signals.add( si++ , l );
    }
  IN1.close();

//   std::map<double,double> tvals;
//   double thf;
//   double empirical_threshold = MiscMath::threshold( x , 0.1 , 20 , 0.1 , &thf, &tvals );

//   std::cerr << "et = " << empirical_threshold << "\n";

//   std::map<double,double>::const_iterator tt = tvals.begin();
//   while ( tt != tvals.end() ) 
//     {
//       std::cerr << tt->first << "\t" 
// 		<< tt->second << "\n";
//       ++tt;
//     }

//   std::exit(1);

//   emd_t emd( x , 1000 );

//   std::exit(1);

  // bad
  si = 0;
  std::ifstream IN2("bad.sigs",std::ios::in);
  while ( !IN2.eof() ) 
    {
      std::string l;
      IN2 >> l;
      if ( IN2.eof() ) break;
      bad_signals.add( si++ , l );
    }
  IN2.close();


  Data::Matrix<double> invG;
  Data::Matrix<double> Gi;
  clocs.make_interpolation_matrices( good_signals , bad_signals , &invG , &Gi );
  std::vector<int> gi;
  for (int i=11;i<=64;i++) gi.push_back(i-1);

  Data::Matrix<double> interp = clocs.interpolate( X , gi , invG , Gi );
  
  std::exit(1);

  //
  // Load signal from STDIN
  //

  if ( 0 ) 
    {
      while ( ! std::cin.eof() ) 
	{
	  double t;
	  std::cin >> t;
	  if ( std::cin.eof() ) break;
	  x.push_back(t);
	}
      std::cerr << "read " << x.size() << " values\n";
    }

  if ( 0 ) 
    {
      // 
//       fiplot_t fp( x , 1000 , 
// 		   0   , 5    , 0.05 ,   // min/max for display, time interval 
// 		   0.5 , 25   , 0.5 );   // frequencies
      
      std::exit(1);
    }


  if ( 0 ) 
    {
      int rows = 1000;
      int cols = 2;
      int compc = 2;
      
      mat X = mat_create( rows , cols );
      
      for (int i=0;i<rows;i++)
	for (int j=0;j<cols;j++)
	  {
	    double x;
	    std::cin >> x;
	    X[i][j] = x;
	  }
      
      //    double ** pW = mat_create(compc, compc);
      //    double ** pA = mat_create(compc, compc);
      //    double ** pK = mat_create(cols, compc);
      //    double ** pS = mat_create(rows, cols);     
      
      //      ica_t ica( X , rows , cols , compc );
      
      //     fastICA(X, rows, cols, compc, pK, pW, pA, pS);
      
      //
      // EMD
      //
      
      ica_t ica( X , rows , cols , compc );
      
//     fastICA(X, rows, cols, compc, pK, pW, pA, pS);
      

   std::cout << "W <comp x comp>\n";

      std::cout << "W <comp x comp>\n";
      
      for (int i=0;i<ica.W.size();i++)
	{
	  for (int j=0;j<ica.W[i].size();j++)
	    std::cout << ica.W[i][j] << " ";
	  std::cout << "\n";
	}
      
      std::cout << "A <comp x comp>\n";
      for (int i=0;i<ica.A.size();i++)
	{
	  for (int j=0;j<ica.A[i].size();j++)
	    std::cout << ica.A[i][j] << " ";
	  std::cout << "\n";
	}
      
      std::cout << "K <cols x comp>\n";
      for (int i=0;i<ica.K.size();i++)
	{
	  for (int j=0;j<ica.K[i].size();j++)
	    std::cout << ica.K[i][j] << " ";
	  std::cout << "\n";
	}
      
      std::cout << "S <rows x cols>\n";
      for (int i=0;i<ica.S.size();i++)
	{
	  for (int j=0;j<ica.S[i].size();j++)
	    std::cout << ica.S[i][j] << " ";
	  std::cout << "\n";
	}
      
    }


  if ( 0 ) 
    {
      hilbert_t hilbert( x );
      std::vector<double> ph = * hilbert.phase();
      for (int k=0;k<ph.size();k++) std::cout << x[k] << " " << ph[k] << "\n";
      std::cout << "\n";
      std::vector<double> f = hilbert.instantaneous_frequency( 1000 );
      for (int k=0;k<f.size();k++) std::cout << k << " " << f[k] << "\n";
      
      std::exit(1);
    }

}


struct cmdsyn_t
{

  bool spacer;

  std::string name;
  std::string desc;
  
  // option-2-description
  std::map<std::string,std::string> req;
  std::map<std::string,std::string> opt;

  cmdsyn_t( const std::string & n , const std::string & d )
    : name(n) , desc(d) , spacer(false) { } 

  cmdsyn_t() : spacer(true) { } 
  
  void requires( const std::string & o , const std::string & d )
  { req[o] = d; }

  void optional( const std::string & o , const std::string & d )
  { opt[o] = d; }
  
  std::string display( const bool verbose = true )
  {

    std::stringstream ss;
    if ( spacer ) { ss << "\n"; return ss.str(); } 

    ss << name << "\t" << desc << "\n";

    if ( req.size() > 0 ) std::cerr << "  required: \n";
    std::map<std::string,std::string>::const_iterator ii = req.begin();
    while( ii != req.end() )
      {
	if ( ii != req.begin() ) ss << "            "; 
	ss << ii->first << "\t" << ii->second << "\n";
	++ii;
      }

    if ( opt.size() > 0 ) std::cerr << "  optional: \n";
    ii = opt.begin();
    while( ii != opt.end() )
      {
	if ( ii != opt.begin() ) ss << "            "; 
	ss << ii->first << "\t" << ii->second << "\n";
	++ii;
      }
    return ss.str();
  }
  
};

void list_cmds()
{
  
  std::vector<cmdsyn_t> cmds;
  
  cmdsyn_t c_write( "WRITE" , "Write a new EDF file" );
  c_write.requires( "tag" , "New tag to add to EDF filename: oldname-tag.edf" );
  c_write.optional( "outdir" , "Set a new output directory for EDFs, must end in '/'") ;
  c_write.optional( "sample-list" , "Append to a sample list for the new EDFs" );
  
  cmdsyn_t c_summary( "SUMMARY" , "Display EDF header information" );
  
  cmdsyn_t c_stats( "STATS" , "Summary statistics for an EDF" );
  
  cmdsyn_t c_uv( "uV" , "Change scale from mV or V to uV" );
  cmdsyn_t c_mv( "mV" , "Change scale from uV or V to mV" );
  
  cmdsyn_t c_timetrack( "TIME-TRACK" , "Add a continuous time-track to an EDF" ); 

  cmdsyn_t c_stage( "STAGE" , "Specify sleep stage labels and generate hypnogram metrics" );
  c_stage.optional( "W" , "WAKE label (default 'W')" );
  c_stage.optional( "N1", "N1 label (default 'N1')" );
  c_stage.optional( "N2", "N2 label (default 'N2')" );
  c_stage.optional( "N3", "N3 label (default 'N3')" );
  c_stage.optional( "R" , "REM label (default 'N3')" );
  c_stage.optional( "?" , "Unscored/unknown label (default '?')" );

  cmdsyn_t c_dump( "DUMP" , " " );
  c_dump.optional( "signal" , "Specify signals" );

  cmdsyn_t c_dump_records( "DUMP-RECORDS" , "" );

  cmdsyn_t c_dump_epochs( "DUMP-EPOCHS" , "" );

  cmdsyn_t c_restructure( "RESTRUCTURE" , "Restructure an EDF (drop masked epochs/channels)" );

  cmdsyn_t c_signals( "SIGNALS" , "Drop/retain specified channels" );
  c_signals.optional( "keep" , "Keep these signals" );
  c_signals.optional( "drop" , "Drop these signals" );
  
  
  cmdsyn_t c_sigstats( "SIGSTATS" , "Signal statistics and epoch-filtering" );
  c_sigstats.optional( "mask" , "");
  c_sigstats.optional( "threshold" , "SD unit outlier removal, can be iterative, e.g. threshold=2,2,2");
  c_sigstats.optional( "lzw" , "Lempel-Ziv-Welch compression index, lzw=nbins,nsmooth default 20,1" );
    
  cmdsyn_t c_mse( "MSE" , "Per-epoch multiscale entropy" );  
  c_mse.optional( "m" , "default 2" );
  c_mse.optional( "r" , "default 0.15" );
  c_mse.optional( "s" , "scale lower/upper and increment; default {lwr,upr,inc} 1,10,2" );

  cmdsyn_t c_zr( "ZR" , "Z-ratio" );

  cmdsyn_t c_anon( "ANON" , "Strip identifiers from EDF headers" );
  
  cmdsyn_t c_epoch( "EPOCH" , "Set epoch duration (sec)" );
  c_epoch.requires( "epoch" , "Epoch duration in seconds, default is 30" );
  
  cmdsyn_t c_slice( "SLICE" , "" );

  cmdsyn_t c_mask( "MASK" , "Apply a mask to hide parts of the data (applied to all signals)" );
  c_mask.optional( "force" , "" );
  c_mask.optional( "unmask" , "" );
  c_mask.optional( "mask" , "" );
  c_mask.optional( "clear" , " (also 'include-all' or 'none')" );
  c_mask.optional( "total" , " (also 'exclude-all' or 'all')" );
  c_mask.optional( "random" , "random=n where n is number of epochs" );
  c_mask.optional( "first" , "first=n where n is number of epochs" );
  c_mask.optional( "leading" , "leading={annot}, e.g. leading=W select leading wake" );
  c_mask.optional( "flanked" , "flanked={annot},n where n is number of epochs either side" );
  c_mask.optional( "include" , "" );
  c_mask.optional( "excldue" , "" );
  c_mask.optional( "label" , "?? still used?" );
  c_mask.optional( "flag" , "?? still used?" );
  
  cmdsyn_t c_epochmask( "EPOCH-MASK" , "based on epoch annotations; ?? difference from MASK??" );
  c_epochmask.optional( "include" , "" );
  c_epochmask.optional( "exclude" , "" );
  
  cmdsyn_t c_filemask ( "FILE-MASK" , "mask from file" );
  c_filemask.optional( "include" , "include=filename" );
  c_filemask.optional( "exclude" , "exclude=filename" );
  
  cmdsyn_t c_dumpmask( "DUMP-MASK" , "write current epoch mask to a file" );
  c_dumpmask.optional( "tag" , "create an .annot file from the mask, rather than standard output" );
  c_dumpmask.optional( "path" , "specify path for this file" );
  
  cmdsyn_t c_epochannot( "EPOCH-ANNOT" , "" );
  c_epochannot.optional( "file" , "" );
  c_epochannot.optional( "recode" , "x=y" );
  
  cmdsyn_t c_filter( "FILTER" , "Apply FIR filter" );
  c_filter.optional( "lower" , "lower HZ" );
  c_filter.optional( "upper" , "upper HZ" );
  c_filter.optional( "num_taps" , "filter order" );
  c_filter.optional( "signal" , "" );
  
  cmdsyn_t c_psd( "PSD" , "Spectral density and band power" );
  c_psd.optional( "spectrum" , "");
  c_psd.optional( "epoch" , "");
  c_psd.optional( "epoch-spectrum" , "" );
  c_psd.optional( "mse" , "" );
  c_psd.optional( "fast-slow-sigma" , "" );
  
  c_psd.optional( "segment-sec" , "Welch algorithm window size, default 4" );
  c_psd.optional( "segment-overlap" , "Window overlap, default 2" );

  c_psd.optional( "ranges" , "ranges=lwr,upr,inc in Hz"  );
  c_psd.optional( "epoch-ranges" , "boolean"  );

  cmdsyn_t c_covar( "COVAR" , "signal covariance" );
  
  cmdsyn_t c_coh ("COH" , "Spectral coherence" );
  c_coh.optional( "sr" , "Sample rate" );
  c_coh.optional( "epoch" , "Output per-epoch band-coherence measures" );

  cmdsyn_t c_bpm( "HR" , "Find R peaks and estimate BPM from an ECG channel" );
  c_bpm.optional( "ecg" , "ECG channel" );
  
  cmdsyn_t c_suppress_ecg( "SUPPRESS-ECG" , "Detect/correct for ECG contamination in signals" );
  c_suppress_ecg.requires( "ecg" , "" );
  c_suppress_ecg.optional( "no-suppress" , "do not update signal" );
  c_suppress_ecg.optional( "sr" , "" );
  
  cmdsyn_t c_pac( "PAC" , "" );
  
  cmdsyn_t c_cfc( "CFC" , "" );
  
  cmdsyn_t c_tag( "TAG" , "" );
  c_tag.requires( "tag" , "" );

  cmdsyn_t c_resample( "RESAMPLE" , "" );
  c_resample.requires( "" , "" );
  
  cmdsyn_t c_spindles( "SPINDLES" , "Detect spindles" );
  c_spindles.optional( "fc" , "" );

  cmdsyn_t c_sw( "SW" , "Detect slow waves" );
  cmdsyn_t c_artifacts( "ARTIFACTS" , "Detect EEG artifacts" );  
  cmdsyn_t c_spike( "SPIKE" , "" );
  
}

param_t build_param_from_cmdline()
{

  param_t param;
  
  while ( ! std::cin.eof() )
    {
      std::string x;
      std::cin >> x;      
      if ( std::cin.eof() ) break;
      if ( x == "" ) continue;
      param.parse( x ); 
    }
  return param;
}

