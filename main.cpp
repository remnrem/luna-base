
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
#include "main.h"
#include "param.h"

//
// global resources
//

extern globals global;

extern writer_t writer;

extern logger_t logger;

extern freezer_t freezer;


int main(int argc , char ** argv )
{
   
  //
  // initial check for display of all commands 
  // and save, if writing to log=luna.log
  //
  
  std::string cmd_dump = log_commands( argc , argv );
    
  //
  // initiate global defintions
  //
  
  std::set_new_handler(NoMem);

  global.init_defs();


  //
  // display version info?
  //
  
  bool show_version = argc >= 2 
    && ( strcmp( argv[1] ,"-v" ) == 0
	 || strcmp( argv[1] ,"--version" ) == 0 );
  
  if ( show_version )  
    {
      global.api();
      std::cerr << luna_base_version() ;
      std::cerr << "Eigen library v"
		<< EIGEN_WORLD_VERSION << "."
		<< EIGEN_MAJOR_VERSION << "."
		<< EIGEN_MINOR_VERSION << "\n";
      std::cerr << "sqlite v"
		<< sqlite3_libversion() << "\n";
      std::exit( globals::retcode );
    }
  
  
  //
  // primary usage
  //
  
  std::string usage_msg = luna_base_version() +
    "url: http://zzz.bwh.harvard.edu/luna/\n"
    "primary usage: luna [sample-list|EDF] [n1] [n2] [id=ID] [@param-file] \n"
    "                    [sig=s1,s2] [var1=val1] [-o out.db] [-s COMMANDS] [< command-file]\n";
  
  //
  // degenerate command line?
  //
  
  if ( argc == 1 && isatty(STDIN_FILENO) )  
    {      
      logger << usage_msg << "\n";
      logger.off();
      std::exit(1);
    }
  
  if ( std::cin.eof() || ! std::cin.good() ) 
    Helper::halt( "no input, quitting" );
  

  //
  // debug/dummy options (internal testing)
  //
  
  if ( argc >= 2 && strcmp( argv[1] ,"-d" ) == 0 )
    { 
      std::string p = argc >= 3 ? argv[2] : "";
      std::string p2 = argc >= 4 ? argv[3] : "";
      global.api();
      proc_dummy( p , p2 ); 
      std::exit( globals::retcode ); 
    } 
  
  
  //
  // help mode
  //

  if ( argc >= 2 && ( strcmp( argv[1] , "-h" ) == 0 || strcmp( argv[1] , "-H" ) == 0 ) ) 
    {
      
      global.api();

      const bool primary = strcmp( argv[1] , "-h" ) == 0;
      
      if ( argc == 2 )  // -h 
	{

	  std::cerr << "\n" << usage_msg << "\n";

	  std::cerr << "List of domains\n"
		    << "---------------\n";
	  
	  std::cerr << globals::cmddefs().help_domains()
		    << "\n";

	  std::cerr << "For commands within a domain, add the domain label after -h, e.g.\n"
		    << "  luna -h annot\n\n";

	  std::cerr << "For options and output for a given command, add the (upper-case) command after -h, e.g.\n"
		    << "  luna -h SIGSTATS\n\n";

	}
      else // --h all
	{
	  std::string p = argv[2] ;

	  // 'all'  list all commands for all domains
	  if ( p == "all" ) 
	    {
	      std::cerr << globals::cmddefs().help_commands() << "\n";

	    }
	  // -h {domain}  list all commands (non-verbose)
	  else if ( globals::cmddefs().is_domain(p) ) 
	    {
	      std::string str = p + " : " + globals::cmddefs().help_domain(p);
	      std::cerr << "\n"
			<< str << "\n"
			<< std::string(str.size(),'-') 
			<< "\n\n"
			<< globals::cmddefs().help_commands( p , primary )
			<< "\n";
	    }
	  
	  // -h {cmd}  list all options/tables (verbose)
	  else if ( globals::cmddefs().is_cmd(p) ) 
	    {
	      std::cerr << globals::cmddefs().help( p , true , true , primary ) << "\n";
	    }
	  
	  // otherwise, complain
	  else  
	    std::cerr << "option [" << p << "] not recognized as a domain or command\n";
	}      
      std::exit(0);
    }



  
  //
  // Parse the command line
  //  --> this may execute special commands (e.g. --build)
  //      and quit, i.e. as we don't want to invoke the
  //      full output mechansism/banner etc
  //
  
  int param_from_command_line;
  
  cmdline_proc_t cmdline = parse_cmdline( argc , argv , &param_from_command_line );

  //
  // mirror log to file?
  //
  
  if ( globals::write_log )
    {
      logger.write_log( Helper::expand( globals::log_file ) );
      
      // if we previously did a command-line dump (w/ --log) then
      // write that out now - it will already have been sent to std::cerr
      // before the logger was started properly

      if ( cmd_dump != "" ) logger.print_to_file( cmd_dump );

    }
  
  //
  // banner
  //
  
  logger.banner( globals::version , globals::date );


  //
  // initialize output to a STOUT db or not?
  //
  
  if ( cmd_t::stout_file.find( globals::indiv_wildcard ) != std::string::npos )
    {
      cmd_t::has_indiv_wildcard = true;
      cmd_t::stout_template = cmd_t::stout_file;
    }
  
  //
  // text-table mode?
  //
  
  if ( cmd_t::plaintext_mode )
    {
      writer.use_plaintext( cmd_t::plaintext_root );
    }
  // else was an output db specified?
  else if ( cmd_t::stout_file != "" )
    {
      
      // if using indiv-specific output databases, postpone this...
      if ( ! cmd_t::has_indiv_wildcard ) 
	{
	  // if not append-mode, first wipe it
	  if ( ! cmd_t::append_stout_file ) 
	    Helper::deleteFile( cmd_t::stout_file );
	  
	  writer.attach( cmd_t::stout_file );  
	}
    }
  // otherwise, just send to std out
  else 
    writer.nodb();

  
  //
  // branch off to run any cmdline driven special functions not already done
  // then quit
  //
  
  if ( cmdline != NOPROC )
    {
      exec_cmdline_procs( cmdline , argc, argv, param_from_command_line ) ;
      std::exit(0);
    }
  
  
  //
  // if here, we're expecting to apply luna commands to one or more data inputs
  //   --> iterate through the primary sample-list
  //
  
  int processed = 0, failed = 0;
  int actually_processed = 0;
  
  while ( ! std::cin.eof()  )
    {
      cmd_t cmd ; // scans STDIN for next command
	    
      if ( cmd.empty() ) break; 

      ++processed;

      if ( ! cmd.valid() )
	++failed;
      else
	{
	  
	  // process command ( most will iterate over 1 or more EDFs)
	  if ( cmd.process_edfs() ) 
	    process_edfs(cmd);
	  else // handle any exceptions 
	    {

	      // [ currently, this single exception not used ...
	      //   --> if needed, better to hook in as a command-line
	      //       function
	      
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

  logger << "...processed " << processed << " command set(s), ";

  if ( failed == 0 ) logger << " all of which passed" << "\n";
  else logger << failed << " of which failed\n";
  
  std::exit( globals::retcode );
  
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
  // Open sample-list, or working with single EDF (or ASCII file)?
  //

  std::string f = cmd.data();

  // use .edf (or .EDF extension, or .rec or .edfz or .edf.gz ) to
  // indicate 'single EDF' mode, '.rec'

  f = f.substr( (int)f.size() - 4 >= 0 ? (int)f.size() - 4 : 0 );

  bool single_edf = Helper::iequals( f , ".edf" ) || Helper::iequals( f , ".rec" ) ;

  if ( ! single_edf ) 
    {

      // also test .edfz or .edf.gz
      f = cmd.data();
      f = f.substr( (int)f.size() - 5 >= 0 ? (int)f.size() - 5 : 0 );
      
      if ( Helper::iequals( f , ".edfz" ) ) 
	single_edf = true;
      else 
	{
	  f = cmd.data();
	  f = f.substr( (int)f.size() - 7 >= 0 ? (int)f.size() - 7 : 0 );
	  if ( Helper::iequals( f , ".edf.gz" ) ) 
	    single_edf = true;	      
	}      
    }
  

  // also allow .sedf for Luna summary EDF

  if ( ! single_edf )
    {
      std::string f2 = cmd.data();
      f2 = f2.substr( (int)f2.size() - 5 >= 0 ? (int)f2.size() - 5 : 0 );
      if ( Helper::iequals( f2 , ".sedf" ) ) single_edf = true;
    }
  

  // use presence of --fs command-line option to indicate 'single ASCII file' mode

  bool single_txt = globals::param.has( "-fs" );
  if ( single_txt ) single_edf = true;
  

  // use presence of '.' name to indicate an empty EDF
  // also (below) allow '.' from SLIST, i.e. incase we have annot-only datasets

  bool empty_edf = f == "." ;
  if ( empty_edf ) single_edf = true;
    
  std::ifstream EDFLIST;
  
  if ( ! single_edf ) 
    EDFLIST.open( cmd.data().c_str() , std::ios::in );
  
  //
  // Do we have a search path for EDFs and ANNOTs?
  //

  bool has_project_path = globals::param.has( "path" );

  if ( has_project_path )
    {

      if ( single_edf ) Helper::halt( "cannot specify project path in single EDF mode" );
      
      globals::project_path = globals::param.value( "path" );

      // does this folder exist?
      if ( ! Helper::fileExists( globals::project_path ) )
	Helper::halt( "could not find project path , " + globals::project_path );
      
      if ( globals::project_path[ globals::project_path.size() - 1 ] != globals::folder_delimiter )
	globals::project_path = globals::project_path + globals::folder_delimiter ; 
      
      logger << "path    : " << globals::project_path << "\n";
                 
    }
  

  //
  // Report on any slices
  //

  if ( ! single_edf )
    {
      if ( globals::sample_list_min != -1 || globals::sample_list_max != -1 )
	{
	  if ( globals::sample_list_min == globals::sample_list_max ) 
	    logger << "row(s)  : #" << globals::sample_list_min << " ( n = 1 )\n";
	  else
	    logger << "row(s)  : #" << globals::sample_list_min << " -> #"
		   << globals::sample_list_max
		   << " ( n = " << globals::sample_list_max - globals::sample_list_min + 1 << " )\n";
	}
      else if ( globals::sample_list_ids.size() )
	logger << "rows(s) : #'s from id arg ( n = "
	       << globals::sample_list_ids.size() << " )\n";
      else	
	logger << "row(s)  : all\n";
    }
  

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
	  
	  Helper::safe_getline( EDFLIST , line);
	  
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

	  // allow '.' missing value for annots?

	  if ( tok.size() == 3 && tok[2] == "." ) tok.resize(2);

	  // ignore SL annots?

	  if ( globals::skip_sl_annots ) tok.resize(2);
	  
	  // allow annot field to be comma delimited? expand out here
	  
	  if ( tok.size() == 3 )
	    {
	      std::vector<std::string> annot_fields = Helper::parse( tok[2] , globals::file_list_delimiter );
	      if ( annot_fields.size() > 1 )
		{
		  tok.resize( 2 + annot_fields.size() );
		  for (int a=0; a<annot_fields.size();a++)
		    tok[a+2] = annot_fields[a];
		}
	    }

	  // remove any quotes (e.g. we're okay to have spaces in names now, as tab delim

	  for (int t=1;t<tok.size();t++)
	    tok[t] = Helper::unquote( tok[t] );
	  
	  // add in project path to relative paths?
	  // (but keep absolute paths as they are)

	  if ( has_project_path )
	    {
	      for (int t=1;t<tok.size();t++)
		{
		  if ( tok[t][0] != globals::folder_delimiter )
		    tok[t] = globals::project_path + tok[t];
		}
	    }
	  
	  // extract main items (ID, signal EDF)

	  rootname = tok[0];
	  edffile  = tok[1];
	  
	  // swap in new ID?
	  rootname = cmd_t::remap_id( rootname );
	  
	  // else, do we have an 'ID' check? (id=ID does not match so skip)
	  if ( globals::sample_list_ids.size() )
	    {	      
	      if ( globals::sample_list_ids.find( rootname )
		   == globals::sample_list_ids.end() )
		{
		  ++processed;
		  continue;
		}
	    }

	  // skip=ID matches
	  if ( globals::sample_list_ids_skips.size() )
	    {	      
	      if ( globals::sample_list_ids_skips.find( rootname )
		   != globals::sample_list_ids_skips.end() )
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
	  
	  // remove .edf from ID, making file name ==> ID 
	  if ( Helper::file_extension( rootname , "edf" ) )
	    rootname = rootname.substr( 0 , rootname.size() - 4 );
	  else if ( Helper::file_extension( rootname , "edfz" ) )
	    rootname = rootname.substr( 0 , rootname.size() - 5 );
	  else if (  Helper::file_extension( rootname , "edf.gz" ) )
	    rootname = rootname.substr( 0 , rootname.size() - 7 );
	  
	  tok.resize(2);
	  tok[0] = rootname;
	  tok[1] = edffile;
	}


      //
      // File in exclude list? (or not in an include list?)
      //

      bool include = true;

      if ( globals::id_excludes.find( rootname ) != globals::id_excludes.end() )
	include = false;

      if ( globals::id_includes.size() != 0
	   && globals::id_includes.find( rootname ) == globals::id_includes.end() )
	include = false;

      if ( ! include )
	{
	  logger << "\n"
		 << "___________________________________________________________________\n"
		 << "  **********************************\n"
		 << "  * Skipping EDF " << rootname << "\n"
		 << "  **********************************\n"
		 << "\n";

	  ++processed;
	  continue; // to the next EDF in the list
	}
      

      //
      // Begin running through the series of commands
      //
      
      logger  << "\n___________________________________________________________________\n"
	      << "Processing: " << rootname 
	      << " [ #" << processed+1 << " ]" << "\n";
      

      //
      // Do we need to open an individual-specific out-db?
      //
      
      if ( cmd_t::has_indiv_wildcard ) 
	{

	  if ( cmd_t::plaintext_mode )
	    Helper::halt( "cannot specify -t and have ^ wild card" );

	  cmd_t::stout_file = cmd_t::resolved_outdb( rootname , cmd_t::stout_template );
	  
	  // if not append-mode, first wipe it
	  if ( ! cmd_t::append_stout_file ) 
	    Helper::deleteFile( cmd_t::stout_file );
	  
	  writer.attach( cmd_t::stout_file );  
	  
	  logger << " writing to " << writer.name() << "\n";

	}
      

      //
      // Begin transaction
      //

      writer.begin();

      writer.clear_tags();

      writer.id( rootname , edffile );
          
      
      //
      // Unset 'problem' and 'empty' flags (i.e. for bailing for this individual)
      //
      
      globals::problem = false;

      globals::empty = false; 

      
      //
      // Limited to specific signals to load in?
      //
      
      const std::set<std::string> * inp_signals = NULL;
      
      if ( cmd.signals().size() > 0 ) inp_signals = &cmd.signals();

      
      //
      // load EDF
      //

      annotation_set_t annotations;
      
      edf_t edf( &annotations );      
      
      bool okay = true; 

      //
      // Read a single input from ASCII
      //

      if ( single_txt )
	{
	  
	  int fs = globals::param.requires_int( "-fs" );
	  std::string startdate = globals::param.has("-date") ? globals::param.value( "-date" ) : "01.01.00" ; 
	  std::string starttime = globals::param.has("-time") ? globals::param.value( "-time" ) : "00.00.00" ; 
	  std::string id = globals::param.has("-id") ? globals::param.value( "-id" ) : rootname ; 

	  std::vector<std::string> labels;
	  if ( globals::param.has("-chs" ) ) labels = globals::param.strvector( "-chs" ); 
	  
	  okay = edf.read_from_ascii( edffile , 
				      id , 
				      fs , 
				      labels , 
				      startdate , starttime );

	}
      else if ( empty_edf || edffile == "." ) // i.e. sample-list empty EDF called for
	{	  

	  //
	  // Generate an empty EDF
	  //

	  const int nr = globals::param.has( "-nr" ) 
	    ? globals::param.requires_int( "-nr" ) 
	    : 60 * 60 * 6 ; // 6 hr default

	  const int rs = globals::param.has( "-rs" )
	    ? globals::param.requires_int( "-rs" )
	    : 1 ; // in full seconds (integer)

	  const std::string startdate = globals::param.has("-date")
	    ? globals::param.value( "-date" )
	    : "01.01.00" ;

	  const std::string starttime = globals::param.has("-time")
	    ? globals::param.value( "-time" )
	    : "00.00.00" ;

	  const std::string id = globals::param.has("-id")
	    ? globals::param.value( "-id" )
	    : rootname ;

	  okay = edf.init_empty( id , nr , rs , startdate , starttime );
	}
      else // attach an EDF from disk
	okay = edf.attach( edffile , rootname , inp_signals ); 



      //
      // Check status of attached data, whatever the source
      //

      if ( ! okay ) 
	{
	  
	  globals::problem = true;

	  logger << "**warning: problem loading " 
		 << edffile << ", skipping..." << "\n";
	  
	  if ( globals::write_naughty_list )
	    {
	      logger << "**writing ID " << edf.id << " to " <<  globals::naughty_list << "\n";
	      std::ofstream PROBLEMS( globals::naughty_list.c_str() , std::ios_base::app );
	      PROBLEMS << edf.id << "\n";
	      PROBLEMS.close();
	    }
	 	  
	  writer.commit();

	  continue;
	}
      

      //
      // ensure we've cleaned the freezer (except this new self)
      //
      
      freezer.clean( &edf );
      
      //
      // Check labels are still unique given aliases
      //

      edf.header.check_channels();

      //
      // Give annotations some basic details about the EDF
      //

      edf.annotations->set( &edf );
      
      //
      // Add additional annotations? 
      //

      for (int i=0;i<globals::annot_files.size();i++)
	{
	  // if absolute path given, add in as in  /home/joe/etc

	  if ( globals::annot_files[i][0] == globals::folder_delimiter ) 
	    tok.push_back( globals::annot_files[i] );
	  else  // project path may be "" if not set; but if set, will end in /
	    tok.push_back( globals::project_path + globals::annot_files[i] );
	}
    
      //
      // Attach annotations
      //
      
      if ( ! globals::skip_nonedf_annots ) 
	{

	  for (int i=2;i<tok.size();i++) 
	    {
	      
	      std::string fname = Helper::expand( tok[i] );
	      
	      if ( fname[ fname.size() - 1 ] == globals::folder_delimiter ) 
		{
		  // this means we are specifying a folder, in which case search for all files that 
		  // start id_<ID>_* and attach thoses
		  DIR * dir;		  
		  struct dirent *ent;
		  if ( (dir = opendir ( fname.c_str() ) ) != NULL )
		    {
		      /* print all the files and directories within directory */
		      while ((ent = readdir (dir)) != NULL)
			{
			  std::string fname2 = ent->d_name;
			  // only annot files (.xml, .ftr, .annot, .eannot)
			  if ( Helper::file_extension( fname2 , "annot" ) ||
			       Helper::file_extension( fname2 , "txt" ) ||
			       Helper::file_extension( fname2 , "tsv" ) ||
			       Helper::file_extension( fname2 , "xml" ) ||
			       Helper::file_extension( fname2 , "ameta" ) ||
			       Helper::file_extension( fname2 , "stages" ) ||
			       Helper::file_extension( fname2 , "eannot" ) )   
			    {
			      edf.load_annotations( fname + fname2 );	 			   
			    }			 
			}
		      closedir (dir);
		    }
		  else 
		    Helper::halt( "could not open folder " + fname );
		}
	      else
		{
		  
		  // only annot files (.xml, .ftr, .annot, .eannot)                                            
		  // i.e. skip .sedf files that might also be specified as 
		  // attached to this EDF
		  if ( Helper::file_extension( fname , "annot" ) ||
		       Helper::file_extension( fname , "txt" ) ||
		       Helper::file_extension( fname , "tsv" ) ||
		       Helper::file_extension( fname , "xml" ) ||
		       Helper::file_extension( fname , "ameta" ) ||
		       Helper::file_extension( fname , "stages" ) ||
		       Helper::file_extension( fname , "eannot" ) )
		    {
		      edf.load_annotations( fname );	 
		    }
		  else
		    Helper::halt( "did not recognize annotation file extension: " + fname );
		} 
	      
	    }
	}

      
      //
      // Attach EDF Annotations, potentially
      //

      if ( edf.header.edfplus )
	{
	  // must read if EDF+D (but only the time-track will be taken in)
	  // if EDF+C, then look at 'skip-edf-annots' flag
	    
	  if ( edf.header.continuous && ! globals::skip_edf_annots )
	    edf.annotations->from_EDF( edf , edf.edfz_ptr() );
	  else if ( ! edf.header.continuous )
	    edf.annotations->from_EDF( edf , edf.edfz_ptr() );

	}
      
      
      //
      // Now, all annotations (except EPOCH-ANNOT) are attached and can be reported on
      //
      
      std::vector<std::string> names = edf.annotations->names();
      
      if ( names.size() > 0 ) logger << "\n annotations:\n";
      
      for (int a = 0 ; a < names.size() ; a++ )
	{
	  
	  annot_t * annot = edf.annotations->find( names[a] );
	  
	  if ( annot == NULL ) Helper::halt( "internal problem in list_all_annotations()" );

	  // do not show special annots [ duration_hms, duration_sec, epoch_sec, start_hms ]
	  if ( annot->special() ) continue;
	  
	  const int num_events = annot->num_interval_events();
	  const int nf = annot->types.size();


	  // verbose mode

	  if ( globals::verbose )
	    {
	      
	      logger << "  [" << names[a] << "] " 
		     << num_events << " instance(s)"
		     << " (from " << annot->file << ")\n";
	      
	      // list instance IDs (up to 4) if multiple or differnt from annot name
	      // but only if there are >1 unique value, *and* the number of unique values
	      // does not equal the total instance count (i.e. do not print if just time-stamp
	      // or count for each ID, only if some coding
	      
	      std::set<std::string> instance_ids = annot->instance_ids();
	      
	      if ( instance_ids.size() > 0 && instance_ids.size() != num_events ) 
		{
		  
		  if ( ! ( instance_ids.size() == 1 
			   && ( *instance_ids.begin()  == names[a] || *instance_ids.begin() == "." ) ) )
		    {
		      logger << "   " << instance_ids.size() << " instance IDs: ";
		      std::set<std::string>::const_iterator ii = instance_ids.begin();
		      int icnt = 0 ; 
		      while ( ii != instance_ids.end() )
			{
			  logger << " " << *ii ;
			  ++icnt;
			  if ( icnt > 4 ) { logger << " ..." ; break;  }
			  ++ii;		  
			}
		      logger << "\n";
		    }
		}
	  	
	      // lists meta-data
	      
	      if ( nf > 1 )
		{
		  logger << "   w/ " << nf << " field(s):";
		  std::map<std::string,globals::atype_t>::const_iterator aa = annot->types.begin();
		  while ( aa != annot->types.end() )
		    {
		      logger << " " << aa->first << "[" << globals::type_name[ aa->second ] << "]";
		      ++aa;
		    }
		  logger << "\n";
		}

	    }

	  //
	  // else non-verbose annotaiton listing
	  //

	  else
	    {
	      if ( a != 0 && a % 4 == 0 ) logger << "\n ";
	      else if ( a != 0 ) logger << " |";
	      else logger << " ";
	      logger << " " << names[a] << " (x" << num_events << ")";
	    }
	  
	}

      logger << "\n";
      
      
      //
      // Automatically generate channel type variables based on attached EDF
      //

      cmd.define_channel_type_variables( edf );


      //
      // Set special 'id' variable to EDF ID
      //

      cmd_t::ivars[ edf.id ][ "id" ] = edf.id;

      
      //
      // List any individual level variables (including new channel type variables, and ${id})
      //

      if ( cmd_t::ivars.find( edf.id ) != cmd_t::ivars.end() )
	{
	  logger << "\n variables:\n";
	  
	  const std::map<std::string,std::string> & newvars = cmd_t::ivars.find( edf.id )->second;
	  std::map<std::string,std::string>::const_iterator vv = newvars.begin();
	  int icnt = 0 ;
	  while ( vv != newvars.end() )
	    {

	      if ( vv->second != "" )
		{
		  if ( globals::verbose )
		    {
		      logger << "  " << vv->first << "=" << vv->second << "\n";
		    }
		  else		    
		    {
		      if ( icnt % 5 == 4 ) logger << "\n ";
		      else if ( icnt != 0 ) logger << " |";
		      else logger << " ";
		      logger << " " << vv->first << "=" << Helper::brief( vv->second , 13 ) ;		    
		      ++icnt;
		    }
		}	      
	      ++vv;
	    }
	  logger << "\n";
	}

      
      //
      // Swap in (indiv-level) variables into the command file
      //  - update any indiv-wildcards in the command list
      //  - include any @includes
      //

      cmd.replace_wildcards( rootname );

      
      //
      // Evaluate all commands
      //

      bool cmd_okay = cmd.eval( edf );
      

      //
      // done 
      //
      
      ++processed;
      ++actual;

      //
      // commit output to DB
      //

      writer.commit();


      //
      // clean up if using individual-specific outdb
      // or if in plaintext mode (i.e. as one folder per individual)
      //
      
      if ( cmd_t::has_indiv_wildcard || cmd_t::plaintext_mode ) 
	writer.close();


      //
      // clean the freezer
      //

      freezer.clean( &edf );


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
  // All done, goodnight.
  //

  logger << "\n"
	 << "___________________________________________________________________"
	 << "\n"
	 << "...processed " << actual << " EDFs, done."
	 << "\n";

}



// ------------------------------------------------------------
// ------------------------------------------------------------
//                 Misc helper functions 
// ------------------------------------------------------------
// ------------------------------------------------------------


// populate paramters from stdin instead of the command line

void build_param_from_stdin( param_t * param )
{
  
  while ( ! std::cin.eof() )
    {
      std::string x;
      std::cin >> x;      
      if ( std::cin.eof() ) break;
      if ( x == "" ) continue;
      param->parse( x ); 
    }

  // swap in wildcards: here, means @{includes} 
  param->update( "." , globals::indiv_wildcard  );

}


//
// construct parameters for command-line functions (e.g. build) 
// from the command line / stdin

void build_param( param_t * param , int argc , char** argv , int start )
{

  // get arguments from stdin (rather than the command line options)?
  
  if ( start == 0 )
    {
      build_param_from_stdin( param );
      return;
    }

  //
  // this only triggered w/ command-line commands, e.g. --fft, etc
  // where we have an --options argument;  start equals the arg after that
  //
  
  for (int i=start; i<argc; i++)
    {
      std::string x = argv[i];
      if ( x == "" ) continue;
      param->parse( x ); 
    }
  
  // swap in wildcards: here, means @{includes} 
  param->update( "." , globals::indiv_wildcard  );

}


//
// report Luna version
//

std::string luna_base_version() 
{
  std::stringstream ss;
  ss << "luna-base version " << globals::version << " (release date " << globals::date << ")\n";
  ss << "luna-base build date/time " << __DATE__ << " " << __TIME__ << "\n";
  return ss.str();
}


//
// "handle" out-of-memory conditions
//

void NoMem()
{
  std::cerr << "*****************************************************\n"
	    << "* FATAL ERROR    Exhausted system memory            *\n"
	    << "*                                                   *\n"
	    << "* You need a smaller dataset or a bigger computer...*\n"
	    << "*                                                   *\n"
	    << "* Forced exit now...                                *\n"
	    << "*****************************************************\n\n";
  std::exit(1);
}


//
// write log output to file
//

std::string log_commands( int argc , char ** argv )
{

  std::stringstream ss;
  
  bool dump = false; 
  
  for (int i=0; i<argc; i++)
    {
      if ( strcmp( argv[i] ,"--log" ) == 0 )
	{
	  dump = true;
	  break;
	}
    }
  
  if ( ! dump ) return "";

  bool has_s = false;
  
  // note - anything in ' ' quotes is read as a single item; need to parse out first;
  std::string str ;
  for (int i=0; i<argc; i++)
    {
      str.append( argv[i] );
      str.push_back( ' ' );
    }

  // remove any tabs and newlines
  str = Helper::search_replace( str , '\n' , ' ');
  str = Helper::search_replace( str , '\t' , ' ');
  
  std::vector<std::string> tok = Helper::parse( str , ' ' );

  const int n = tok.size();
  
  ss << "\n"
     << "# " << std::string( 78 , '=' ) << "\n"
     << "\n";
  
  for (int i=0; i<n; i++)
    {
      std::string s = tok[i] ;
      if ( s == "--log" ) continue;
      
      std::string spc = i == 0 ? "" : " ";
      if ( ( ! has_s )  && ( s[0] == '-' || s.find( "=" ) != std::string::npos ) ) spc = " \\\n     ";
      else if ( s == "&" ) spc = " \n        ";  // no newline '\' char

      // now in command string?
      if ( s == "-s" ) { s = "-s '"; has_s = true; }
      
      ss << spc << s ;
      
    }
  
  if ( has_s ) ss << " '";

  ss << "\n\n"
     << "# " << std::string( 78 , '-' ) << "\n"
     << "\n";

  // send to std::cerr
  std::cerr << ss.str();
  
  // & return string (if want to save to file later, via log=file.log)
  return ss.str();
  
}


// 
// parse the command line
//

cmdline_proc_t parse_cmdline( int argc , char ** argv , int * param_from_command_line )
{

  // no inputs?
  if ( argc == 1 ) return NOPROC;

  // implement '--options' feature (i.e. args following
  //  go into a param_t that is passed to the special
  //  cmdline-function
    
  // use standard input versus command line for
  // command-line options (e.g. --massoc, --psc, etc)
  //
  // (default) get options from stdin: (param_from_command_line == 0)
  //
  //   echo "load=file.dat rows" | luna --massoc -o out.db @param.txt vars=phe.txt
  //
  // else, if --options appears as a command line option, then take
  //   everything after that as options for --command i.e. to build
  //   the param_t object
  //
  //   luna --massoc -o out.db @param.txt --options vars=phe.txt load=file.data rows
  //
  
  *param_from_command_line = 0 ; 
  
  for (int i=1; i<argc; i++)
    {
      if ( strcmp( argv[i] , "--options" ) == 0 ||
	   strcmp( argv[i] , "--opt" ) == 0 )
	{
	  *param_from_command_line = i+1;
	  break;
	}
    }	
  
  
  //
  // pick off any special functions here:
  // if the first arg is recognized as a function
  // rather than some form of data input
  //
  
  std::string arg1( argv[1] );
  
  // map of command line options

  std::map<std::string,cmdline_proc_t> clmap;
  
  // primary
  clmap[ "--validate" ]     = PROC_VALIDATE ;
  clmap[ "--build" ]        = PROC_BUILD;
  clmap[ "--repath" ]       = PROC_REPATH;
  clmap[ "--merge" ]        = PROC_MERGE;
  clmap[ "--bind" ]         = PROC_BIND;
  clmap[ "--xml" ]          = PROC_XML1;
  clmap[ "--xml2" ]         = PROC_XML2;
  clmap[ "--eval" ]         = PROC_EVAL;
  clmap[ "--eval-verbose" ] = PROC_EVAL_VERBOSE;  
 
  // signal processing 
  clmap[ "--fft" ]          = PROC_FFT ;
  clmap[ "--cwt-design" ]   = PROC_CWT_DESIGN ;
  clmap[ "--cwt" ]          = PROC_CWT_DESIGN ;
  clmap[ "--fir-design" ]   = PROC_FIR_DESIGN ;

  // linear models / analysis
  clmap[ "--gpa-prep" ]     = PROC_GPA_PREP ;
  clmap[ "--gpa" ]          = PROC_GPA_RUN ;
  clmap[ "--cpt" ]          = PROC_CPERM_TEST ;
  clmap[ "--overlap" ]      = PROC_OVERLAP ;
  
  // dimension reduction (SVD/NMF)
  clmap[ "--psc" ]          = PROC_PSC;
  clmap[ "--nmf" ]          = PROC_NMF;

  // microstate analysis helper functions
  clmap[ "--kmer" ]         = PROC_MS_KMER;
  clmap[ "--cmp-maps" ]     = PROC_MS_CMP_MAPS;
  clmap[ "--correl-maps" ]  = PROC_MS_CORR_MAPS;
  clmap[ "--label-maps" ]   = PROC_MS_LABEL_MAPS;

  // POPS staging
  clmap[ "--pops" ]         = PROC_POPS;
  clmap[ "--eval-stages" ]  = PROC_EVAL_STAGES;
  clmap[ "--priors" ]       = PROC_POPS_ESPRIORS;
  
  // some further redundant / obscure options 
  clmap[ "--otsu" ]         = PROC_OTSU;
  clmap[ "--copy-suds" ]    = PROC_COPY_SUDS;
  clmap[ "--combine-suds" ] = PROC_COMBINE_SUDS;
  clmap[ "--lgbm" ]         = PROC_LGBM;
  clmap[ "--assoc" ]        = PROC_ASSOC;
  clmap[ "--massoc" ]       = PROC_MASSOC;
  clmap[ "--pdc" ]          = PROC_PDC;
  clmap[ "--pdlib" ]        = PROC_PDLIB;
  clmap[ "--mapper" ]       = PROC_MAPPER;
  clmap[ "--mapper-html" ]  = PROC_MAPPER_HTML;

  // do we recognize the first arg as a special cmd?
  // if so, stop now but pass the cmdline function
  // to main();  we still want to process other command line
  // options (before --options) as they may contain e.g.
  // info on the -o out.db etc
  
  cmdline_proc_t cmdline = NOPROC;
  
  if ( clmap.find( arg1 ) != clmap.end() )
    cmdline = clmap[ arg1 ];

  
  
  //
  // otherwise, first element will be treated as the data source
  // (sample list, EDF or text file) 
  //

  if ( cmdline == NOPROC )
    cmd_t::input = argv[1];

  //
  // other args should be in form
  //
  //   @param-file
  //   key=value
  //   a single number or range (n, n-m) or slice (n/m)
  //   sig=S1,S2,...  (e.g. special variable)
  //   v1=N2          (e.g. standard project-level variable)
  //   --flag  ( added minus fist '-' to globals::param ), e.g. --fs, --labels
  
  //   -o  output db
  //   -a  output db  { as above, except append } 
  //   -t  output dir { text-table outputs }
  //   -s  { rest of line is script }


  //
  // process each arg
  //

  int specified = 0;
  
  for (int i=2;i<argc;i++)
    {
      
      // if we've had an --options command (which sets 'param_from_command_line')
      // then ignore any options at or past this value, i.e. as they are specific
      // for the command line tool
      
      if ( *param_from_command_line != 0 && i >= *param_from_command_line ) 
	continue;
      
      // parse for a key=value form
      
      std::vector<std::string> tok = 
	Helper::quoted_parse( argv[i] , "=" );
	  
      const bool is_key_value = tok.size() == 2;
      
      // if so, parse to see if it is a special variale
      // and make any settings changes;  this also
      // saves this as a standard project-level variable
      
      if ( is_key_value )
	{
	  cmd_t::parse_special( tok[0] , tok[1] ); 	      
	  continue;
	}
      
	  
      // otherwise, see if this is a cmdline flag (e.g. -s, -o etc or
      // an @include file)

      // -o : specify database for output
	  
      if ( Helper::iequals( tok[0] , "-o" ) || Helper::iequals( tok[0] , "-a" ) )
	{
	  // next arg will be DB
	  if ( i + 1 >= argc )
	    Helper::halt( "expecting database name after -o/a" );
	  cmd_t::stout_file = argv[ ++i ];
	  if ( Helper::iequals( tok[0] , "-a" ) ) cmd_t::append_stout_file = true;
	  continue;
	}
      

      // -t : specify database for output
      
      if ( Helper::iequals( tok[0] , "-t" ) )
	{
	  // next arg will be root (folder) for plain-text
	  if ( i + 1 >= argc )
	    Helper::halt( "expecting database name after -t" );
	  cmd_t::plaintext_root = argv[ ++i ];
	  cmd_t::plaintext_mode = true;
	  continue;
	}
      
      // -s : luna-script from command line
      
      if ( Helper::iequals( tok[0] , "-s" ) )
	{
	  if ( cmdline != NOPROC )
	    Helper::halt( "cannot specify -s with a special command-line function" );
	  
	  // rest of args will be cmd script: add it here, then we're done
	  for (int j=i+1;j<argc;j++)
	    cmd_t::add_cmdline_cmd( argv[j] );

	  return NOPROC;
	}
      
      // other special variable? e.g. --fs=200
      
      if ( argv[i][0] == '-' )
	{
	  std::string f = argv[i];	      
	  globals::param.add( f.substr(1) );
	  continue;
	}
      
      
      // @param file
      
      if ( argv[i][0] == '@' )
	{
	  include_param_file( argv[i] );
	  continue;
	}
      

      //
      // by this point, we're now expecting either an ID, a
      // numeric code (or two) or an n/m slice, i.e. args to
      // slice a sample list;  therefore, if we have a special
      // command, we can skip this and continue to the next args
      //

      if ( cmdline != NOPROC ) continue;

      //
      // otherwise, determine if we're subsetting the sample-list
      //
      
      int x;
      std::string term( argv[i] );
      const bool is_number = Helper::str2int( term , &x );
      
      // an n/m 'slice' of the sample list? 
      const bool nm_slice = term.find( "/" ) != std::string::npos; 
      
      // an string ID? 
      if ( ! ( nm_slice || is_number ) )
	{
	  // add to the ID list
	  globals::sample_list_ids.insert( term );
	  
	  // either way, we're now done selecting n m values
	  specified = 2;
	  
	  continue;
	}
      	  
      // n/m slice? 
      
      if ( nm_slice ) 
	{		  
	  
	  // set slice mode to fix globals::sample_list_min
	  // and globals::sample_list_max		  
	  
	  std::vector<std::string> tok = Helper::parse( term , "/" );
	  if ( tok.size() != 2 )
	    Helper::halt( "expecting n/m format ID specification: use id= if the ID contains '/'" );
	  
	  int n = 0 , m = 0;
	  if ( ! Helper::str2int( tok[0] , &n ) )
	    Helper::halt( "expecting integer n/m format for sample specification" );
	  if ( ! Helper::str2int( tok[1] , &m ) )
	    Helper::halt( "expecting integer n/m format for sample specification" );
	  if ( n < 1 || m < 1 || n > m )
	    Helper::halt( "expecting integer n/m format for sample specification, n & m >= 1 and n <= m" );
	  
	  int s1 = 0, s2 = 0;
	  bool okay = Helper::sl_slicer( cmd_t::input , n , m , &s1, &s2 );
	  if ( ! okay ) Helper::halt( "problem setting n/m sample slice" );
	  globals::sample_list_min = s1;
	  globals::sample_list_max = s2;
	  
	  // we're now done selecting
	  specified = 2;
	  
	  continue;
	}
      

      // at this point, should only be left w/ a numeric (int) value
      // we'll test nonetheless

      if ( is_number ) 
	{
	  
	  // first or second value?
	  if ( specified == 0 )
	    {
	      globals::sample_list_min = globals::sample_list_max = x;
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

      
    } // process next arg

  
  //
  // check any final sample list selections
  //

  if ( globals::sample_list_max < globals::sample_list_min )
    {
      int x = globals::sample_list_max;
      globals::sample_list_max = globals::sample_list_min;
      globals::sample_list_min = x;
    }
  
  if ( globals::sample_list_min < 0 )
    globals::sample_list_min = -1;
  
  if ( globals::sample_list_max < 0 )
    globals::sample_list_max = -1;

  // let caller know if we're expecting a special command-line call
  
  return cmdline;
      
}


//
// include a @param file
//

void include_param_file( const std::string & paramfile )
{

  // otherwise parse as a normal line: i.e. two tab-delim cols
  // alternatively, allow space/equals delimiters
  std::string delim = "\t";
  if ( globals::allow_space_param ) delim += " ";
  if ( globals::allow_equals_param ) delim += "=";
  
  bool parse_line = true;

  std::string last_grp = "";
  
  // allow missing parameter file "." i.e. to make scripting 
  // easier for LSF submission script that need to pass this 
  
  if ( paramfile.size() > 1 && paramfile != "@." ) 
    {
		  
      // expand() expands out any ~/ notation to full path
      const std::string filename = Helper::expand( paramfile.substr(1).c_str() );
      if ( ! Helper::fileExists( filename ) ) Helper::halt( "could not open " + filename );
		  
      std::ifstream INC( filename.c_str() , std::ios::in );
      if ( INC.bad() ) Helper::halt("could not open file: " + filename );
      while ( ! INC.eof() )
	{
	  
	  std::string line;
	  
	  Helper::safe_getline( INC , line );	  
	  if ( INC.eof() || line == "" ) continue;
		      
	  // skip % comments
	  if ( line[0] == '%' ) continue;
	  
	  // is this an include/exclude section
	  // +group  include only if matches group, otherwise skip
	  // -group  exclude if matches group, otherwise parse
	  
	  if ( line[0] == '+' || line[0] == '-' )
	    {
	      const std::string grp = line.substr(1);
	      
	      if ( grp == "" ) continue;
	      
	      if ( last_grp == "" ) last_grp = line;
	      else if ( last_grp != line )
		Helper::halt( "cannot nest +group/-group lines" );
	      else last_grp = "";
	      
	      bool has_grp =
		cmd_t::vars.find( grp ) != cmd_t::vars.end() ?
		Helper::yesno( cmd_t::vars[ grp ] ) : false ;
	      
	      if ( line[0] == '-' &&   has_grp ) parse_line = ! parse_line;
	      if ( line[0] == '+' && ! has_grp ) parse_line = ! parse_line;
	      
	      // skip to next line now
	      continue;
	    }
	  else
	    {
	      // if not a control line +grp or -grp, and if we are not parsing, then skip
	      if ( ! parse_line ) continue;
	    }
	  
		      
	  
	  std::vector<std::string> tok = Helper::quoted_parse( line , delim );

	  if ( tok.size() != 2 )
	    Helper::halt("badly formatted line parameter file line ( # cols != 2 ) in " + filename
			 + "\n" + line
			 + "\n - see param-spaces and param-equals option possibly, which are set to T by default" 
			 );
	  
	  cmd_t::parse_special( tok[0] , tok[1] );
	  
	}
      
      INC.close();
    }
  
}


//
// execute special command line functions, then quit
//


void exec_cmdline_procs( cmdline_proc_t & cmdline , int argc , char ** argv, int param_from_command_line )
{

  if ( cmdline == NOPROC ) return;
  
  //
  // PSC (or NMF)
  //
  
  if ( cmdline == PROC_PSC || cmdline == PROC_NMF )
    {
      if ( cmdline == PROC_PSC && cmdline == PROC_NMF )
	Helper::halt( "cannot specify both --psc and --nmf" );
      
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );            
      writer.begin();      
      // writer.id( "." , "." );
      const std::string clab = cmdline == PROC_PSC ? "PSC" : "NMF" ;       
      writer.cmd( clab , 1 , "" );
      writer.level( clab, "_" + clab );

      // PSC , or NMF mode:
      psc_t psc;
      psc.construct( param , cmdline == PROC_NMF );
    
      writer.unlevel( "_" + clab );
      writer.commit();
      std::exit(0);
    }


  //
  // GPA
  //

  if ( cmdline == PROC_GPA_PREP || cmdline == PROC_GPA_RUN )
    {
      if ( cmdline == PROC_GPA_PREP && cmdline == PROC_GPA_RUN )
	Helper::halt( "cannot specify both --prep-gpa and --gpa" );
      
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );            
      writer.begin();
      writer.id( "." , "." );
      writer.cmd( "GPA" , 1 , "" );
      writer.level( "GPA", "_GPA" );
      
      // run GPA (prep or run)
      gpa_t gpa( param , cmdline == PROC_GPA_PREP ); 
      
      writer.unlevel( "_GPA" );
      writer.commit();
      std::exit(0);
    }

  //
  // POPS
  //
  
  if ( cmdline == PROC_POPS )
    {
#ifdef HAS_LGBM
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );
      writer.begin();
      writer.id( "." , "." );
      writer.cmd( "POPS" , 1 , "" );
      writer.level( "POPS", "_POPS" );

      pops_t pops( param );
      pops.make_level2_library( param );      

      writer.unlevel( "_POPS" );
      writer.commit();
#else
      Helper::halt( "LGBM support not compiled in" );
#endif
      std::exit(0);
    }

  //
  // POPS EVAL STAGES
  //

  if ( cmdline == PROC_EVAL_STAGES ) 
    {
#ifdef HAS_LGBM
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );
      writer.begin();
      writer.id( "." , "." );
      writer.cmd( "EVAL-STAGES" , 1 , "" );
      writer.level( "EVAL-STAGES", "_EVAL-STAGES" );

      pops_indiv_t indiv( param , param.requires( "file" ) , param.value( "file2" ) );
      
      writer.unlevel( "_EVAL-STAGES" );
      writer.commit();
#else
      Helper::halt( "LGBM support not compiled in" );
#endif

      std::exit(0);
    }


  //
  // Make es-priors (as a standalone function)
  //
  
  if ( cmdline == PROC_POPS_ESPRIORS )
    {
#ifdef HAS_LGBM
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );
      writer.begin();
      writer.id( "." , "." );
      writer.cmd( "POPS" , 1 , "" );
      writer.level( "POPS", "_POPS" );

      pops_t pops( param );
      pops.make_espriors( param );      

      writer.unlevel( "_POPS" );
      writer.commit();
#else
      Helper::halt( "LGBM support not compiled in" );
#endif
      std::exit(0);
    }

  
  //
  // Simple LGBM wrapper
  //

  if ( cmdline == PROC_LGBM )
    {
#ifdef HAS_LGBM
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );
      writer.begin();
      writer.id( "." , "." );
      writer.cmd( "LGBM" , 1 , "" );
      writer.level( "LGBM", "_LGBM" );

      lgbm_cli_wrapper( param );      

      writer.unlevel( "_LGBM" );
      writer.commit();
#else
      Helper::halt( "LGBM support not compiled in" );
#endif
      std::exit(0);
    }
  

  //
  // LGBM-based ASSOC
  //

  if ( cmdline == PROC_ASSOC )
    {
#ifdef HAS_LGBM
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );
            
      writer.begin();
      writer.id( "." , "." );
      writer.cmd( "ASSOC" , 1 , "" );
      writer.level( "ASSOC", "_ASSOC" );

      assoc_t assoc( param );

      writer.unlevel( "_ASSOC" );
      writer.commit();
#else
      Helper::halt( "LGBM support not compiled in" );
#endif
      std::exit(0);
    }


  //
  // LGBM-based MASSOC (i.e. TLOCK-based)
  //
  
  if ( cmdline == PROC_MASSOC )
    {
#ifdef HAS_LGBM
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );
      writer.begin();
      writer.id( "." , "." );
      writer.cmd( "MASSOC" , 1 , "" );
      writer.level( "MASSOC", "_MASSOC" );

      massoc_t massoc( param );

      writer.unlevel( "_MASSOC" );
      writer.commit();
#else
      Helper::halt( "LGBM support not compiled in" );
#endif
      std::exit(0);
    }

  //
  // Basic FFT on stdin signal
  //

  if ( cmdline == PROC_FFT )
    {
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );
      writer.begin();      
      writer.id( "." , "." );
      writer.cmd( "FFT" , 1 , "" );
      writer.level( "FFT", "_FFT" );      

      dsptools::cmdline_fft( param );

      writer.unlevel( "_FFT" );
      writer.commit();
      std::exit(0);      
    }


  //
  // --validate sample list
  //

  if ( cmdline == PROC_VALIDATE )
    {
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );      
      writer.begin();
      writer.id( "." , "." );
      writer.cmd( "VALIDATE" , 1 , "" );
      writer.level( "VALIDATE", "_VALIDATE" );

      Helper::validate_slist( param );

      writer.unlevel( "_VALIDATE" );
      writer.commit();
      std::exit(0);
    }

  
  //
  // OVERLAP enrichment (multi-sample case)
  //

  if ( cmdline == PROC_OVERLAP )
    {
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );      
      writer.begin();      
      writer.id( "." , "." );
      writer.cmd( "OVERLAP" , 1 , "" );
      writer.level( "OVERLAP", "_OVERLAP" );      

      annotate_t annotate( param );

      writer.unlevel( "_OVERLAP" );
      writer.commit();
      std::exit(0);      
    }

  //
  // Otsu thresholding
  //

  if ( cmdline == PROC_OTSU )
    {
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );      
      writer.begin();      
      writer.id( "." , "." );
      writer.cmd( "OTSU" , 1 , "" );
      writer.level( "OTSU", "_OTSU" );      

      dsptools::cmdline_otsu( param );

      writer.unlevel( "_OTSU" );
      writer.commit();
      std::exit(0);      
    }

  //
  // Cluster permitation test (CPT)
  //

  if ( cmdline == PROC_CPERM_TEST )
    {      
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );
      writer.begin();      
      writer.id( "." , "." );
      writer.cmd( "CPT" , 1 , "" );
      writer.level( "CPT", "_CPT" );      

      cpt_wrapper( param );

      writer.unlevel( "_CPT" );
      writer.commit();
      std::exit(0);
    }



  //
  // KMER analysis
  //

  if ( cmdline == PROC_MS_KMER )
    {
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );            
      writer.begin();     
      writer.id( "." , "." );      
      writer.cmd( "KMER" , 1 , "" );
      writer.level( "KMER", "_KMER" );

      dsptools::ms_kmer_wrapper( param );
      
      writer.unlevel( "_KMER" );
      writer.commit();
      std::exit(0);
    }
  

  //
  // MS map label
  //
  
  if ( cmdline == PROC_MS_LABEL_MAPS )
    {
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );
      writer.begin();      
      writer.id( "." , "." );      
      writer.cmd( "LABEL-MAPS" , 1 , "" );
      writer.level( "LABEL-MAPS", "_LABEL-MAPS" );

      dsptools::ms_label_maps( param );
      
      writer.unlevel( "_LABEL-MAPS" );
      writer.commit();
      std::exit(0);
    }



  //
  // MS correl maps
  //

  if ( cmdline == PROC_MS_CORR_MAPS )
    {
      param_t param;      
      build_param( &param, argc, argv, param_from_command_line );      
      writer.begin();      
      writer.id( "." , "." );      
      writer.cmd( "CORREL-MAPS" , 1 , "" );
      writer.level( "CORREL-MAPS", "_CORREL-MAPS" );

      dsptools::ms_correl_maps( param );
                  
      writer.unlevel( "_CORREL-MAPS" );
      writer.commit();
      std::exit(0);
    }


  //
  // MS compare maps
  //

  if ( cmdline == PROC_MS_CMP_MAPS )
    {
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );      
      writer.begin();      
      writer.id( "." , "." );      
      writer.cmd( "CMP-MAPS" , 1 , "" );
      writer.level( "CMP-MAPS", "_CMP-MAPS" );

      dsptools::ms_cmp_maps( param );
	      
      writer.unlevel( "_CMP-MAPS" );
      writer.commit();
      std::exit(0);
    }



  //
  // FIR design
  //

  if ( cmdline == PROC_FIR_DESIGN )
    {

      writer.begin();      
      writer.id( "." , "." );      
      writer.cmd( "FIR-DESIGN" , 1 , "" );
      writer.level( "FIR-DESIGN", "_FIR-DESIGN" );
      
      // expects input from std::cin
      proc_filter_design_cmdline();
      
      writer.unlevel( "_FIR-DESIGN" );
      writer.commit();
      std::exit(0);
    }


  if ( cmdline == PROC_CWT_DESIGN )
    {
      writer.begin();      
      writer.id( "." , "." );
      writer.cmd( "CWT-DESIGN" , 1 , "" );
      writer.level( "CWT-DESIGN", "_CWT-DESIGN" );

      // expects input from std::cin
      proc_cwt_design_cmdline();

      writer.commit();
      std::exit(0);

    }

  
  if ( cmdline == PROC_COPY_SUDS )
    {
      // expects input from std::cin
      proc_copy_suds_cmdline();
      std::exit(0);
    }


  if ( cmdline == PROC_COMBINE_SUDS )
    {
      // expects input from std::cin
      proc_combine_suds_cmdline();
      std::exit(0);
    }
  
  
  // --eval
  
  if ( cmdline == PROC_EVAL )
    {
      global.api();
      proc_eval_tester( false ); 
      std::exit(0);       
    }
  

  // --eval-verbose
  
  if ( cmdline == PROC_EVAL_VERBOSE )
    {
      // as above, but w/ verbose output
      global.api();
      proc_eval_tester( true ); 
      std::exit(0);       
    }


  // --pdc (PDC helper)
  
  if ( cmdline == PROC_PDC )
    {
      param_t param;      
      build_param_from_stdin( &param );      
      writer.nodb();
      writer.begin();      
      writer.id( "." , "." );
      pdc_t pdc;
      pdc.external( param );
      writer.commit();
      std::exit(0);
    }
    

  // --xml / --xml2 (dump an XML file)
  
  if ( argc == 3 && ( cmdline == PROC_XML1 || cmdline == PROC_XML2 ) )
    {
      global.api();
      annot_t::dumpxml( argv[2] , cmdline == PROC_XML2 );
      std::exit(0);
    }

  
  // --build: construct a project list
  
  if ( argc >=2 && cmdline == PROC_BUILD )
    {
      global.api();
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++)
	tok.push_back( argv[i] );
      Helper::build_sample_list( tok );
      std::exit(0);
    }

   
  // --repath: change paths for a sample list
  
  if ( argc >=2 && cmdline == PROC_REPATH )
    {
      global.api();
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++)
	tok.push_back( argv[i] );
      Helper::repath_SL( tok );
      std::exit(0);
    }


  // --merge: merge EDFs (rbind) 

  if ( argc >=2 && cmdline == PROC_MERGE )
    {      
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++)
	tok.push_back( argv[i] );
      Helper::merge_EDFs( tok );
      std::exit(0);
    }

  
  //
  // --bind: merge EDFs (cbind)
  //

  if ( argc >=2 && cmdline == PROC_BIND )
    {      
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++)
	tok.push_back( argv[i] );
      Helper::bind_EDFs( tok );
      std::exit(0);
    }

  // --pdlib

  if ( argc == 2 && cmdline == PROC_PDLIB )
    {
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );      
      writer.nodb();
      writer.begin();      
      writer.id( "." , "." );
      pdc_t::construct_pdlib( param );
      writer.commit();
      std::exit(0);
    }

  //
  // map channels/ annots
  //

  if ( argc >=2 && cmdline == PROC_MAPPER )
    {
      global.api();      
      // expecting form: cmap=xxx amap=xxx c=xxx a=yyy 
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++)
	tok.push_back( argv[i] );
      Helper::channel_annot_mapper( tok , false ) ;
      std::exit(0);
    }
  
  //
  // map channels/ annots, HTML style output
  //
  
  if ( argc >=2 && cmdline == PROC_MAPPER_HTML )
    {
      global.api();
      
      // expecting form: cmap=xxx amap=xxx c=xxx a=yyy 
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++)
	tok.push_back( argv[i] );
      Helper::channel_annot_mapper( tok , true ) ;
      std::exit(0);
    }
  
}


  
//
// misc helper: EVAL expresions
//

void proc_eval_tester( const bool verbose )
{

  // read a single line
  std::string expr;
  Helper::safe_getline( std::cin , expr );

  std::map<std::string,annot_map_t> inputs;

  instance_t out;
  
  Eval tok( expr );

  tok.bind( inputs , &out );
  
  bool is_valid = tok.evaluate( verbose );
  
  bool retval;
  
  if ( ! tok.value( retval ) ) is_valid = false;

  std::cout << "parsed as a valid expression : " << ( is_valid ? "yes" : "no" ) << "\n";
  std::cout << "return value                 : " << tok.result() << "\n";
  std::cout << "return value (as T/F)        : " << ( retval ? "true" : "false" ) << "\n";
  std::cout << "assigned meta-data           : " << out.print() << "\n";  
  std::exit(1);

}
