
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

#include "helper/token.h"

#include "helper/token-eval.h"

#include "luna.h"

#include "main.h"


extern globals global;

extern writer_t writer;

extern logger_t logger;

int main(int argc , char ** argv )
{

  bool show_version = argc >= 2 && ( strcmp( argv[1] ,"-v" ) == 0 || strcmp( argv[1] ,"--version" ) == 0 );
  
  //
  // initiate global defintions
  //

  global.init_defs();

  if ( show_version )  
    {
      global.api();
      std::cerr << luna_base_version() ;
      std::exit(0);
    }

  
  //
  // Some initial options (called prior to the main banner, etc)
  //

  if ( argc >= 2 && strcmp( argv[1] ,"-d" ) == 0 )
    { 
      std::string p = argc >= 3 ? argv[2] : "";
      global.api();
      proc_dummy( p ); 
      exit(0); 
    } 
  
  //
  // help mode
  //

  else if ( argc >= 2 && strcmp( argv[1] , "-h" ) == 0 )
    {

      global.api();

      
      if ( argc == 2 )  // -h 
	{

	  std::cerr << "\nusage: luna [sample-list|EDF] [n1] [n2] [@parameter-file] [sig=s1,s2] [v1=val1] < command-file\n\n";

	  std::cerr << "List of domains\n"
		    << "---------------\n\n";
	  
	  std::cerr << globals::cmddefs.help_domains()
		    << "\n";

	  std::cerr << "for commands within a domain, add the domain label after -h, e.g.\n"
		    << "\n  luna -h annot\n\n";

	  std::cerr << "for options and output for a given command, add the (upper-case) command after -h, e.g.\n"
		    << "\n  luna -h SIGSTATS\n\n";

	}
      else // --h all
	{
	  std::string p = argv[2] ;

	  // 'all'  list all commands for all domains
	  if ( p == "all" ) 
	    {
	      std::cerr << globals::cmddefs.help_commands() << "\n";

	    }
	  // -h {domain}  list all commands (non-verbose)
	  else if ( globals::cmddefs.is_domain(p) ) 
	    {
	      std::cerr << "\n" << globals::cmddefs.help_commands( p ) << "\n";
	    }
	  
	  // -h {cmd}  list all options/tables (verbose)
	  else if ( globals::cmddefs.is_cmd(p) ) 
	    {
	      std::cerr << globals::cmddefs.help( p , true , true ) << "\n";
	    }
	  
	  // otherwise, complain
	  else  
	    std::cerr << "option [" << p << "] not recognized as a domain or command\n";
	}      
      std::exit(0);
    }

  //
  // EVAL from the command line 
  //
  else if ( argc == 2 && strcmp( argv[1] , "--eval" ) == 0 ) 	   
    {
      global.api();
      proc_eval_tester( false ); 
      exit(0);       
    }

  //
  // Verbose EVAL
  //

  else if ( argc == 2 && strcmp( argv[1] , "--eval-verbose" ) == 0 ) 
    {
      // as above, but w/ verbose output
      global.api();
      proc_eval_tester( true ); 
      exit(0);       
    }


  //
  // DUMP an XML file
  //

  else if ( argc == 3 && strcmp( argv[1] , "--xml" ) == 0 )
    {
      global.api();
      annot_t::dumpxml( argv[2] , false );
      std::exit(0);
    }


  //
  // build a project list
  //

  else if ( argc >=2 && strcmp( argv[1] , "--build" ) == 0 )
    {
      global.api();
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++) tok.push_back( argv[i] );
      Helper::build_sample_list( tok );
      std::exit(0);
    }


  //
  // banner
  //

  logger.banner( globals::version , globals::date );


  //
  // special command-line driven functions that do not involve
  // iterating through a sample list
  //

  bool cmdline_proc_fir_design = false;
  bool cmdline_proc_cwt_design = false;
  bool cmdline_proc_pdlib      = false;

  //
  // parse command line
  //

  if ( argc == 2 && strcmp( argv[1] , "--pdlib" ) == 0 ) 
    {
      param_t param;
      build_param_from_cmdline( &param );
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
      //   --flag  ( added minus fist '-' to globals::param ), e.g. --fs, --labels
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

	  // is this a special variable, e.g. path=, alias=, signal=, etc
	  
	  if ( tok.size() == 2 ) 
	    {
	      cmd_t::parse_special( tok[0] , tok[1] ); 	      
	    }
	  
	  // specify database for output
	  
	  else if ( Helper::iequals( tok[0] , "-o" ) || Helper::iequals( tok[0] , "-a" ) )
	    {
	      // next arg will be DB
	      if ( i + 1 >= argc ) Helper::halt( "expecting database name after -o/a" );
	      cmd_t::stout_file = argv[ ++i ];
	      if ( Helper::iequals( tok[0] , "-a" ) ) cmd_t::append_stout_file = true;
	    }
	  
	  // specify database for output
	  
	  else if ( Helper::iequals( tok[0] , "-t" ) )
	    {
	      
	      // next arg will be root (folder) for plain-text
	      if ( i + 1 >= argc ) Helper::halt( "expecting database name after -t" );
	      cmd_t::plaintext_root = argv[ ++i ];
	      cmd_t::plaintext_mode = true;
	    }
	  
	  // luna-script from command line
	  
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

	  
	  // param file
	  
	  else if ( argv[i][0] == '@' )
	    {

	      // an 'include' 
	      std::string filename = argv[i];
	      
	      // expand() expands out any ~/ notation to full path
	      filename = Helper::expand( filename.substr(1).c_str() );
	      if ( ! Helper::fileExists( filename ) ) Helper::halt( "could not open " + filename );
	      
	      std::ifstream INC( filename.c_str() , std::ios::in );
	      if ( INC.bad() ) Helper::halt("could not open file: " + filename );
	      while ( ! INC.eof() )
		{

		  std::string line;
		  
		  //std::getline( INC , line);		  
		  Helper::safe_getline( INC , line );
		  
		  if ( INC.eof() || line == "" ) continue;
		  
		  // skip % comments
		  if ( line[0] == '%' ) continue;

		  std::vector<std::string> tok = Helper::quoted_parse( line , "\t" );
		  if ( tok.size() != 2 ) Helper::halt("badly formatted line in " + filename );
		  
		  cmd_t::parse_special( tok[0] , tok[1] );		  

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

      logger << "usage: luna [sample-list|EDF] [n1] [n2] [@parameter-file] [sig=s1,s2] [v1=val1] < command-file"
	     << "\n";
      logger.off();
      std::exit(1);
    }

  if ( std::cin.eof() || ! std::cin.good() ) 
    Helper::halt( "no input, quitting" );



  //
  // initialize output to a STOUT db or not?
  //
  
  if ( cmd_t::stout_file.find( globals::indiv_wildcard ) != std::string::npos )
    {
      cmd_t::has_indiv_wildcard = true;
      cmd_t::stout_template = cmd_t::stout_file;
    }
  
  // text-table mode?
  if ( cmd_t::plaintext_mode )
    {
      logger << "in text-table mode, writing to " << cmd_t::plaintext_root  << "\n";
      writer.use_plaintext( cmd_t::plaintext_root );
    }
  // was an output db specified?
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
  // branch off to run any cmdline driven special functions, then quit
  //

  if ( cmdline_proc_fir_design )
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
  
  while ( ! std::cin.eof()  )
    {
      
      cmd_t cmd ; // else scans STDIN for next command
      
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

  logger << "...processed " << processed << " command set(s), ";
  if ( failed == 0 ) logger << " all of which passed" << "\n";
  else logger << failed << " of which failed\n";

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
  // Open sample-list, or working with single EDF (or ASCII file)?
  //

  std::string f = cmd.data();

  // use .edf or .EDF extension to indicate 'single EDF' mode
  f = f.substr( (int)f.size() - 4 >= 0 ? (int)f.size() - 4 : 0 );
  bool single_edf = f == ".edf" || f == ".EDF";
  
  // use presence of --fs command-line option to indicate 'single ASCII file' mode
  bool single_txt = globals::param.has( "-fs" );
  if ( single_txt ) single_edf = true;

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

	  tok.resize(2);
	  tok[0] = rootname;
	  tok[1] = edffile;
	}


      //
      // File in exclude list?
      //

      if ( globals::excludes.find( rootname ) != globals::excludes.end() )
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

	  if ( cmd_t::plaintext_mode ) Helper::halt( "cannot specify -t and have ^ wild card" );

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
      // Update any indiv-wildcards in the command list
      //


      cmd.replace_wildcards( rootname );

      
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
      
      bool okay = true; 
      
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
      else
	okay = edf.attach( edffile , rootname , inp_signals ); // read EDF
      
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
      // Check labels are still unique given aliases
      //

      edf.header.check_channels();

      
      //
      // Add additional annotations? 
      //

      for (int i=0;i<globals::annot_files.size();i++)
	{
	  // if absolute path given, add in as in  /home/joe/etc

	  if ( globals::annot_files[i][0] == globals::folder_delimiter ) 
	    tok.push_back( globals::annot_files[i] );
	  else  // projecyt path may be "" if not set; but if set, will end in /
	    tok.push_back( globals::project_path + globals::annot_files[i] );
	}
    
          
      //
      // Attach annotations
      //
      
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
		      if ( Helper::file_extension( fname2 , "ftr" ) ||
			   Helper::file_extension( fname2 , "xml" ) ||
			   Helper::file_extension( fname2 , "eannot" ) ||
			   Helper::file_extension( fname2 , "annot" ) )
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
	      edf.load_annotations( fname );	 
	    }
	  
	}

	  
      //
      // Now, all annotations (except EPOCH-ANNOT) are attached and can be reported on
      //
      
      std::vector<std::string> names = edf.timeline.annotations.names();
      
      if ( names.size() > 0 ) logger << "\n annotations:\n";
      
      for (int a = 0 ; a < names.size() ; a++ )
	{
	  
	  annot_t * annot = edf.timeline.annotations.find( names[a] );
	  
	  if ( annot == NULL ) Helper::halt( "internal problem in list_all_annotations()" );
	  
	  const int num_events = annot->num_interval_events();
	  const int nf = annot->types.size();
	  
	  logger << "  [" << names[a] << "] " 
		 << num_events << " instance(s)"
		 << " (from " << annot->file << ")\n";
	  
	  // list instance IDs (up to 8) if multiple or differnt from annot name
	  
	  std::set<std::string> instance_ids = annot->instance_ids();
	  
	  if ( instance_ids.size() > 0 ) 
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
		      if ( icnt > 8 ) { logger << " ..." ; break;  }
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

  logger << "\n"
	 << "___________________________________________________________________"
	 << "\n"
	 << "...processed " << actual << " EDFs, done."
	 << "\n";

}



// EVAL expresions

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



// DUMMY : a generic placeholder/scratchpad for templating new things

void proc_dummy( const std::string & p )
{
  
  if ( p == "cmddefs" ) 
    {
      
      globals::cmddefs.add_domain( "misc" , "misc" ,  "Misc" );
      
      globals::cmddefs.add_cmd( "misc" , "comm1" , "this is a dummy command" );
      globals::cmddefs.add_table( "comm1" , "XX" , "A simple table" , false );
      globals::cmddefs.add_var( "comm1" , "XX" , "X", "Var X" );
      globals::cmddefs.add_var( "comm1" , "XX" , "Y" , "Var Y" );
      
      globals::cmddefs.add_cmd( "misc" , "comm2" , "this is a dummy command" );
      globals::cmddefs.add_table( "comm2" , "CH,B" , "A nice table" , true );
      globals::cmddefs.add_var( "comm2" , "CH,B" , "V1" , "Variable 1" );
      globals::cmddefs.add_var( "comm2" , "CH,B" , "V2" , "Variable 2" );
      
      //   std::cout << globals::cmddefs.help( "comm1" , true )  << "\n\n\n";
      std::cout << globals::cmddefs.help( "comm2" , true )  << "\n";
      
      // add a dummy tag
      globals::cmddefs.add_tag( "Z" );
      
      zfiles_t files( "folder1" , "indiv1" ); 
      
      zfile_t * z1 = files.file( "comm1" , NULL , "XX" ) ; 
      
      param_t param2;
      
      zfile_t * z2 = files.file( "comm2" , &param2 , "CH,B,Z" ) ; 
      
      z1->write_header();
      
      z2->write_header();
      
      //z1->display() ;
      //z2->display() ;
      
      z1->set_stratum( "XX" , "L1" );
      z1->set_value( "X" , 22 );
      z1->set_value( "Y" , 23 );
      z1->write_buffer();
      z1->set_stratum( "XX" , "L2" );
      z1->set_value( "X" , 24 );
      z1->set_value( "Y" , 25 );
      z1->write_buffer();
      
      z2->set_stratum( "CH" , "C3" );
      z2->set_stratum( "B" , "ALPHA" );
      z2->set_stratum( "Z" , "R1" );
      z2->set_value( "V1" , 22 );
      z2->set_value( "V2" , 23 );
      z2->write_buffer();
      
      files.close();
      
      std::exit(1);
      
    }
  
  //
  // Test cmddefs_t
  //
  
  if ( p == "cmddefs" )
    {
      
      globals::cmddefs.add_cmd( "misc"   , "NEWONE" , "A test command" );
      globals::cmddefs.add_table( "NEWONE" , "" , "Table 0, baseline" );
      globals::cmddefs.add_table( "NEWONE" , "CH" , "Table 1, by channel" );
      globals::cmddefs.add_table( "NEWONE" , "CH,X" , "Table 2, by channel and X" );
      globals::cmddefs.add_table( "NEWONE" , "CH,X,Y" , "Table 2a, by channel and X/Y" , true );
      globals::cmddefs.add_table( "NEWONE" , "CH,X,Z" , "Table 2b, by channel and X/Z"  );

      globals::cmddefs.add_var( "NEWONE" , "" , "V1" , "some var1" );
      globals::cmddefs.add_var( "NEWONE" , "" , "V2" , "some var2" );

      globals::cmddefs.add_var( "NEWONE" , "CH" , "V1" , "some var1" );
      globals::cmddefs.add_var( "NEWONE" , "CH" , "V2" , "some var2" );
      globals::cmddefs.add_var( "NEWONE" , "CH" , "V3" , "some var3" );

      globals::cmddefs.add_var( "NEWONE" , "CH,X" , "V2a" , "some var2" );
      globals::cmddefs.add_var( "NEWONE" , "CH,X" , "V3a" , "some var3" );

      globals::cmddefs.add_var( "NEWONE" , "CH,X,Y" , "V2a" , "some var2" );
      globals::cmddefs.add_var( "NEWONE" , "CH,X,Y" , "V3a" , "some var3" );

      globals::cmddefs.add_var( "NEWONE" , "CH,X,Z" , "V2a" , "some var2" );
      globals::cmddefs.add_var( "NEWONE" , "CH,X,Z" , "V3a" , "some var3" );

      //  std::cout << globals::cmddefs.help_domains() << "\n";
      //std::cout << globals::cmddefs.help_commands() << "\n";
      
      globals::cmddefs.set_compressed( "NEWONE" , tfac_t( "CH,X,Z" ) );
      globals::cmddefs.set_compressed( "NEWONE" , tfac_t( "CH,X,Y" ) , false );
      
      std::cout << globals::cmddefs.help( "NEWONE" , true ) << "\n";
      
      
      std::exit(0);

      param_t param;
      param.add( "epoch" );
      param.add( "ep" );
      std::set<std::string> unk;
      std::cout << "st = " << globals::cmddefs.check( "ANNOTS" , param.keys() , &unk ) << "\n";
      
      std::set<std::string>::const_iterator uu = unk.begin();
      while ( uu != unk.end() ) { std::cout << " bad param " << *uu << "\n"; ++uu; } 



      
      

    }

  //
  // Straight FFT of stdin
  //

  std::vector<double> x;
  
  if ( p == "fir" || p == "fft" || p == "mtm" || p == "tv" 
       || p == "dynam" || p == "ica" || p == "fip" || p == "sl" ) 
    {

      int cnt= 0;
      while ( ! std::cin.eof() )
	{
	  double xx;
	  std::cin >> xx;
	  if ( std::cin.bad() ) { std::cerr << "bad input\n"; std::exit(1);  } 
	  if ( std::cin.eof() ) break;	  
	  x.push_back( xx );	  
	  if ( ++cnt % 100000  == 0 ) std::cerr << cnt << "\n";
	}
      std::cerr << x.size() << " values read\n";

    }

  if ( p == "fip" )
    {
      const int sr = 256;
      const uint64_t fs = globals::tp_1sec / sr;
      std::vector<uint64_t> tp( x.size() );
      for (int i=0;i<tp.size();i++) tp[i] = i * fs;
      double th = 0;
      bool norm = false;
      bool logit = false;
      double t_lwr = 0.1; double t_upr = 5;  double t_inc = 0.1;
      double f_lwr = 1;   double f_upr = 20; double f_inc = 0.5;
      bool logspace = false;
      bool cycles = false;

      int num_cyc = 7;

      fiplot_t fp( x , &tp , sr , 
		   th , norm , logit , 
		   t_lwr, t_upr, t_inc , cycles , 
		   f_lwr, f_upr, f_inc , num_cyc , logspace );

      std::exit(0);
    }


  if ( p == "fir" ) 
    {

      //
      // FIR
      //
      
      double ripple = 0.1;
      double tw = 3;
      double f1 = 2;
      double f2 = 15;
      double fs = 1000;

      std::vector<double> fc = dsptools::design_bandpass_fir( ripple , tw , fs , f1, f2 );
      
      //       std::vector<double> fc = dsptools::design_bandpass_fir( 0.001 , 0.1 , 1000 , 2 , 15 );
      std::cerr << "bandpass FIR order " << fc.size() << "\n";
      fir_impl_t fir_impl ( fc );
      x = fir_impl.filter( &x );
      for (int i=0;i<x.size();i++) std::cout << x[i] << "\n";
      //      x = MiscMath::Z( x ) ;
      
      std::exit(1);
    }

  if ( p == "fft" )
    {
      
      
      int index_length = x.size();
      int my_Fs = 1000; // arbitrary
      int index_start = 0;

      FFT fftseg( index_length , my_Fs , FFT_FORWARD , WINDOW_NONE );
      
      fftseg.apply( &(x[index_start]) , index_length );

      // Extract the raw transform
      std::vector<std::complex<double> > t = fftseg.transform();

      // Extract the raw transform scaled by 1/n
      std::vector<std::complex<double> > t2 = fftseg.scaled_transform();
      
      int my_N = fftseg.cutoff;

      //std::cout << t.size() << "\t" << t2.size() << "\t" << my_N << "\n";
      
      for (int f=0;f<my_N;f++)
	{
	  std::cout << "M" << f << "\t" << fftseg.frq[f] << "\t" << fftseg.X[f] << "\n";
	}
 
//       for (int f=0;f<t.size();f++)
// 	std::cout << "V2\t" << t[f] << "\t" << t2[f] << "\n";
      
      std::exit(1);
    }
  
  if ( p == "sl" ) 
    {

      // CLOCLS
      
      clocs_t clocs;
      
      clocs.load_cart( "/Users/shaun/dropbox/projects/ltest/clocs.eegbook" );

      int i = 0 ; 
      signal_list_t signals;
      
      std::ifstream II( "/Users/shaun/dropbox/projects/ltest/clocs.eegbook" );
      while ( ! II.eof() ) 
	{
	  std::string l;
	  double x, y, z;
	  II >> l >> x >> y >> z;
	  if ( II.eof() ) break;
	  signals.add( i++ , l );
	  
	}
      II.close();

      sl_t sl( clocs , signals );
      
      
      // assume 64 channels; rows = channels; cols = time-points
      const int ns = 64;
      const int np = x.size() / ns;
      Data::Matrix<double> X( np , ns );

      i = 0;
      for (int c=0;c<ns;c++)
	for (int t=0;t<np;t++) 
	  X[t][c] = x[i++];
      
      Data::Matrix<double> O;

      sl.apply( X , O );
      
      
    }

  if ( p == "dynam" )
    {
      
      dynam_t dynam( x );

      double beta, rsq = 0;

      dynam.linear_trend( &beta , &rsq );
      
      std::cout << "beta = " << beta << "\n";
      std::cout << "rsq = " << rsq  << "\n";
      
      std::exit(0);
    }

  if ( p == "mse" )
    {

      while ( ! std::cin.eof() )
	{
	  double xx;
	  std::cin >> xx;
	  if ( std::cin.eof() ) break;	  
	  x.push_back( xx );	  
	}
      std::cerr << x.size() << " values read\n";
      
      mse_t mse( 1,20,1,2,0.15 );
      
      std::map<int,double> mses = mse.calc( x );

      std::map<int,double>::const_iterator ii = mses.begin();
      while ( ii != mses.end() )
	{
	  std::cout << ii->first << "\t" << ii->second << "\n";
	  ++ii;
	}

      std::exit(1);
    }

  //
  // db -> retval
  //

  if ( p == "db" )
    {
      //      retval_t retval = writer_t::dump_to_retval( "/Users/shaun/dropbox/projects/nsrr/tutorial/luna/out.db" , "nsrr01" );
      retval_t retval = writer_t::dump_to_retval( "/Users/shaun/dropbox/projects/nsrr/tutorial/luna/out.db" );

      retval.dump(); 
      std::cout << "\n";
      std::exit(1);
    }

  
  //
  // ICA
  //

  if ( p == "ica" ) 
    {


      // asssume two signals for now
      const int n2 = x.size();
      const int ns = 61;

      int rows = x.size() / ns;
      int cols = ns;
      Data::Matrix<double> X( rows , cols );

      int p = 0;
      
      // nb ex.dat is col major
      for (int j=0;j<ns;j++)
	for (int i=0;i<rows;i++) 	
	  X[i][j] = x[p++];
      
      int compc = 2;
      
      std::cout << "want to perform ICA on " << rows << " x " << cols << " matrix\n";
         
      ica_t ica( X , compc );

      std::cout << "done HERE\n";
      
      for (int i=0;i<rows;i++)
	{
	  std::cout << "ICA " << i ;
	  
// 	  for (int j=0;j<cols;j++) 
// 	    std::cout << "\t" << pX[i][j] << "\t";
	  
  	  for (int j=0;j<compc;j++) 
  	    std::cout << "\t" << ica.S[i][j] ;  // component then time-points (ie. col-major for S )       
	 
	  std::cout << "\n";
 
	}
      
      
      // other matrices
      
      // K : cols x compc
      // A : compc x compc
      // W : compc x compc
      // S : as original data
      
      std::cout << "K\n";
      for (int i=0;i<cols;i++)
	{
	  for (int j=0;j<compc;j++) std::cout << "\t" << ica.K[i][j];
	  std::cout << "\n";
	}
       

      std::cout << "W\n";
      for (int i=0;i<compc;i++)
 	{
 	  for (int j=0;j<compc;j++) std::cout << "\t" << ica.W[i][j];
 	  std::cout << "\n\n";
 	}
      
      std::cout << "A\n";
      for (int i=0;i<compc;i++)
 	{
 	  for (int j=0;j<compc;j++) std::cout << "\t" << ica.A[i][j];
 	  std::cout << "\n\n";
 	}      
      
    }
  

  //
  // retval test
  //

  if ( p == "retval" )
    {
      
      edf_t edf;
      edf.attach( "/Users/shaun/my-dropbox/my-sleep/Purcell06072016.edf" , "smp" );

      // mimic R leval() behavior
      retval_t retval ;
      
      writer.use_retval( &retval );

      // set command string                                                                                                                                                                                                                                                               
      cmd_t cmd( "PSD epoch sig=EEG1 & SPINDLES fc=11,15 sig=EEG1" );

      cmd.eval( edf );

      writer.use_retval( NULL );

      retval.dump();
	
      std::exit(1);
    }
  
  //
  // Windows
  //

  if ( p == "windows" )
    {
      const int N = 100;
      std::vector<double> W1(N), W2(N), W3(N);
      W1 = MiscMath::tukey_window(N,0.5);      
      W2 = MiscMath::hann_window(N);    
      W3 = MiscMath::hamming_window(N);
      
      for (int i=0;i<N;i++)
	std::cout << W1[i] << "\t"
		  << W2[i] << "\t"
		  << W3[i] << "\n";


      
      std::exit(1);
    }
  
  
  //
  // MTM
  //
  
  if ( p == "mtm" )
    {
      const int npi = 3;
      
      const int nwin = 5;
      
      mtm_t mtm( npi , nwin );
      
      mtm.display_tapers = false;
      
      mtm.apply( &x , 20 );
      
      std::cout << mtm.f.size() << "\t" << mtm.spec.size() << "\n";
      
      for (int f=0;f<mtm.f.size();f++)
	std::cout << "MTM\t" << f << "\t" << mtm.f[f] << "\t" << mtm.spec[f] << "\n";
      
      std::exit(0);
    }
  
  
  //
  // TV
  //

  if ( p == "tv" )
    {
      
      double lambda = 10;

      std::vector<double> y = dsptools::TV1D_denoise_copy( x , lambda );
      
      for (int i=0;i<x.size();i++)
	std::cout << x[i] << "\t" << y[i] << "\n";
      
      std::exit(1);
    }


  //
  // topop
  //

  if ( p == "topo" )
    {

      // read map from 'example.topo'   CH  THETA  RAD 
      // read data from std::cin,       CH  VALUE
      
      topo_t topo;
      
      int ch = topo.load( "example.topo" );
      
      topo.max_radius( 0.55 );
      
      topo.grid( 67 , 67 );
      
      std::map<std::string, double> data;  
      
      while ( ! std::cin.eof() ) 
	{
	  std::string l;
	  double z;
	  std::cin >> l;
	  if ( std::cin.eof() ) continue;
	  if ( l == "" ) continue;
	  std::cin >> z ;
	  data[l] = z;	  
	}

      std::cerr << "read topo for " << ch << " channels\n";

      std::cerr << "read data for " <<  data.size() << " channels\n";

      Data::Matrix<double> I = topo.interpolate( data );
      
      std::cout << I.dump() << "\n";
      
      std::exit(0);
    }



  if ( p == "clocs" )
    {

  
      clocs_t clocs;
      clocs.load_cart( "ex.clocs" );
      ///Users/shaun/dropbox/sleep/projects/grins-test/clocs.txt" );
      

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
      
      //   logger << "et = " << empirical_threshold << "\n";
      
      //   std::map<double,double>::const_iterator tt = tvals.begin();
      //   while ( tt != tvals.end() ) 
      //     {
      //       logger << tt->first << "\t" 
      // 		<< tt->second << "\n";
      //       ++tt;
      //     }
      
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
      
    }

  
  //
  // end of dummy()
  //
      
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

    if ( req.size() > 0 ) logger << "  required: \n";
    std::map<std::string,std::string>::const_iterator ii = req.begin();
    while( ii != req.end() )
      {
	if ( ii != req.begin() ) ss << "            "; 
	ss << ii->first << "\t" << ii->second << "\n";
	++ii;
      }

    if ( opt.size() > 0 ) logger << "  optional: \n";
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

void build_param_from_cmdline( param_t * param )
{
  
  while ( ! std::cin.eof() )
    {
      std::string x;
      std::cin >> x;      
      if ( std::cin.eof() ) break;
      if ( x == "" ) continue;
      param->parse( x ); 
    }
}


std::string luna_base_version() 
{
  std::stringstream ss;
  ss << "luna-base version " << globals::version << " (release date " << globals::date << ")\n";
  ss << "luna-base build date/time " << __DATE__ << " " << __TIME__ << "\n";
  return ss.str();
}

