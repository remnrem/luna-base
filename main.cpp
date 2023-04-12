
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
#include "miscmath/crandom.h"
#include <fstream>

#include "utils/cgi-utils.h"

extern globals global;

extern writer_t writer;

extern logger_t logger;

void log_commands( int argc , char ** argv );

int main(int argc , char ** argv )
{

  //
  // initial check for display of all commands
  //
  
  log_commands( argc , argv );
  
  //
  // display version info?
  //
  
  bool show_version = argc >= 2 
    && ( strcmp( argv[1] ,"-v" ) == 0 || strcmp( argv[1] ,"--version" ) == 0 );
  
  //
  // initiate global defintions
  //
  
  std::set_new_handler(NoMem);

  global.init_defs();
  
  if ( show_version )  
    {
      global.api();
      std::cerr << luna_base_version() ;
      std::exit( globals::retcode );
    }

  
  //
  // usage
  //
  
  std::string usage_msg = luna_base_version() +
    "url: http://zzz.bwh.harvard.edu/luna/\n"
    "primary usage: luna [sample-list|EDF] [n1] [n2] [id=ID] [@param-file] \n"
    "                    [sig=s1,s2] [var1=val1] [-o out.db] [-s COMMANDS] [< command-file]\n";
  
  
  //
  // Some initial options (called prior to the main banner, etc)
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

  else if ( argc >= 2 && strcmp( argv[1] , "-h" ) == 0 )
    {

      global.api();

      
      if ( argc == 2 )  // -h 
	{

	  std::cerr << "\n" << usage_msg << "\n\n";

	  std::cerr << "List of domains\n"
		    << "---------------\n\n";
	  
	  std::cerr << globals::cmddefs().help_domains()
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
	      std::cerr << globals::cmddefs().help_commands() << "\n";

	    }
	  // -h {domain}  list all commands (non-verbose)
	  else if ( globals::cmddefs().is_domain(p) ) 
	    {
	      std::cerr << "\n" << globals::cmddefs().help_commands( p ) << "\n";
	    }
	  
	  // -h {cmd}  list all options/tables (verbose)
	  else if ( globals::cmddefs().is_cmd(p) ) 
	    {
	      std::cerr << globals::cmddefs().help( p , true , true ) << "\n";
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
  // PDC helper
  //
  
  else if ( argc == 2 && strcmp( argv[1] , "--pdc" ) == 0 )
    {
      param_t param;      
      build_param_from_stdin( &param );
      //build_param_from_args( &param , argc , argv );
      writer.nodb();
      writer.begin();      
      writer.id( "." , "." );
      pdc_t pdc;
      pdc.external( param );
      writer.commit();
      std::exit(0);
    }
    
  //
  // DUMP an XML file
  //

  else if ( argc == 3 && ( strcmp( argv[1] , "--xml" ) == 0 || strcmp( argv[1] , "--xml2" ) == 0 ) )
    {
      bool raw_format =  strcmp( argv[1] , "--xml2" ) == 0 ;
      global.api();
      annot_t::dumpxml( argv[2] , raw_format );
      std::exit(0);
    }

    
  //
  // build a project list
  //

  if ( argc >=2 && strcmp( argv[1] , "--build" ) == 0 )
    {
      global.api();
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++) tok.push_back( argv[i] );
      Helper::build_sample_list( tok );
      std::exit(0);
    }

  
  //
  // change paths
  //

  if ( argc >=2 && strcmp( argv[1] , "--repath" ) == 0 )
    {
      global.api();
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++) tok.push_back( argv[i] );
      Helper::repath_SL( tok );
      std::exit(0);
    }


  //
  // merge EDfs  
  //

  if ( argc >=2 && strcmp( argv[1] , "--merge" ) == 0 )
    {
      //global.api();
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++) tok.push_back( argv[i] );
      Helper::merge_EDFs( tok );
      std::exit(0);
    }

  
  //
  // map channels/ annots
  //

  if ( argc >=2 && strcmp( argv[1] , "--mapper" ) == 0 )
    {
      global.api();
      
      // expecting form: cmap=xxx amap=xxx c=xxx a=yyy 
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++) tok.push_back( argv[i] );
      Helper::channel_annot_mapper( tok , false ) ;
      std::exit(0);
    }
  
  //
  // map channels/ annots, HTML style output
  //
  
  if ( argc >=2 && strcmp( argv[1] , "--mapper-html" ) == 0 )
    {
      global.api();
      
      // expecting form: cmap=xxx amap=xxx c=xxx a=yyy 
      std::vector<std::string> tok;
      for (int i=2;i<argc;i++) tok.push_back( argv[i] );
      Helper::channel_annot_mapper( tok , true ) ;
       std::exit(0);
    }
  
  
  //
  // special command-line driven functions that do not involve
  // iterating through a sample list
  //
   
  bool cmdline_proc_fir_design    = false;
  bool cmdline_proc_cwt_design    = false;
  bool cmdline_proc_pdlib         = false;
  bool cmdline_proc_psc           = false;
  bool cmdline_proc_nmf           = false;
  bool cmdline_proc_ms_kmer       = false;
  bool cmdline_proc_ms_cmp_maps   = false;
  bool cmdline_proc_ms_corr_maps  = false;
  bool cmdline_proc_ms_label_maps = false;
  bool cmdline_proc_copy_suds     = false;
  bool cmdline_proc_combine_suds  = false;
  bool cmdline_proc_cperm_test    = false;
  bool cmdline_proc_lgbm          = false;
  bool cmdline_proc_assoc         = false;
  bool cmdline_proc_massoc        = false;
  bool cmdline_proc_pops          = false;
  bool cmdline_proc_pops_espriors = false;
  bool cmdline_proc_eval_stages   = false;
  bool cmdline_proc_otsu          = false;
  bool cmdline_proc_fft           = false;
  bool cmdline_proc_overlap       = false;

  //
  // use standard input versus command line for
  // command-line options (e.g. --massoc, --psc, etc)
  //
  
  // e.g. default : get options from stdin: (param_from_command_line == 0)
  //
  //   echo "load=file.dat rows" | luna --massoc -o out.db @param.txt vars=phe.txt
  //
  // else, if --options appears as a command line option, then take everything after that as options for --command
  //   i.e. to build the param_t object 
  //
  //   luna --massoc -o out.db @param.txt --options vars=phe.txt load=file.data rows
  //
  
  int param_from_command_line = 0 ; 
  
  for (int i=0; i<argc; i++)
    {
      if ( strcmp( argv[i] , "--options" ) == 0 ||
	   strcmp( argv[i] , "--opt" ) == 0 )
	{
	  param_from_command_line = i+1;
	  break;
	}
    }	
  
  
  //
  // parse command line
  //

  if ( argc == 2 && strcmp( argv[1] , "--pdlib" ) == 0 ) 
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
	  else if ( strcmp( argv[1] , "--cwt-design" ) == 0 ||  
		    strcmp( argv[1] , "--cwt" ) == 0 ) 	   
	    cmdline_proc_cwt_design = true;
	  else if ( strcmp( argv[1] , "--psc" ) == 0 )
	    cmdline_proc_psc = true;
	  else if ( strcmp( argv[1] , "--nmf" ) == 0 )
	    cmdline_proc_nmf = true;
	  else if ( strcmp( argv[1] , "--cpt" ) == 0 )
	    cmdline_proc_cperm_test = true;
	  else if ( strcmp( argv[1] , "--kmer" ) == 0 )
	    cmdline_proc_ms_kmer = true;
	  else if ( strcmp( argv[1] , "--cmp-maps" ) == 0 )
	    cmdline_proc_ms_cmp_maps = true;
	  else if ( strcmp( argv[1] , "--correl-maps" ) == 0 )
	    cmdline_proc_ms_corr_maps = true;
	  else if ( strcmp( argv[1] , "--label-maps" ) == 0 )
	    cmdline_proc_ms_label_maps = true;
	  else if ( strcmp( argv[1] , "--copy-suds" ) == 0 ) 
	    cmdline_proc_copy_suds = true;
	  else if ( strcmp( argv[1] , "--combine-suds" ) == 0 ) 
	    cmdline_proc_combine_suds = true;
	  else if ( strcmp( argv[1] , "--lgbm" ) == 0 )
	    cmdline_proc_lgbm = true;
	  else if ( strcmp( argv[1] , "--assoc" ) == 0 )
	    cmdline_proc_assoc = true;	  
	  else if ( strcmp( argv[1] , "--massoc" ) == 0 )
	    cmdline_proc_massoc = true;	  
	  else if ( strcmp( argv[1] , "--pops" ) == 0 )
	    cmdline_proc_pops = true;
	  else if ( strcmp( argv[1] , "--eval-stages" ) == 0 )
	    cmdline_proc_eval_stages = true;
	  else if ( strcmp( argv[1] , "--priors" ) == 0 )
	    cmdline_proc_pops_espriors = true;
	  else if ( strcmp( argv[1] , "--otsu" ) == 0 )
	    cmdline_proc_otsu = true;
	  else if ( strcmp( argv[1] , "--fft" ) == 0 )
	    cmdline_proc_fft = true;
	  else if ( strcmp( argv[1] , "--overlap" ) == 0 )
	    cmdline_proc_overlap = true;
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

	  // if we've had an --options command (which sets 'param_from_command_line')
	  // then ignore any options at or past this value, i.e. as they are specific
	  // for the command line tool

	  if ( param_from_command_line != 0 && i >= param_from_command_line ) 
	    continue;

	  // parse for a key=value form: if we do not have this, we assume
	  // it is a special variable (e.g. @include) or a numeric SL row-range
	  
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
	      
	      bool parse_line = true;
	      std::string last_grp = "";
	      
	      // allow missing parameter file "." i.e. to make scripting 
	      // easier for LSF submission script that need to pass this 

	      if ( filename.size() > 1 && filename != "@." ) 
		{
		  
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
		      
		      
		      // otherwise parse as a normal line: i.e. two tab-delim cols		      
		      std::vector<std::string> tok = Helper::quoted_parse( line , "\t" );
		      if ( tok.size() != 2 )
			Helper::halt("badly formatted line ( # tabs != 2 ) in " + filename + "\n" + line );
		      
		      cmd_t::parse_special( tok[0] , tok[1] );
		      
		    }
		  
		  INC.close();
		}
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

      logger << usage_msg
	     << "\n";
      logger.off();
      std::exit(1);
    }

  if ( std::cin.eof() || ! std::cin.good() ) 
    Helper::halt( "no input, quitting" );


  
  //
  // -------- done parsing command args --------
  //


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
  
  // text-table mode?
  if ( cmd_t::plaintext_mode )
    {
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


  //
  // PSC, or NMF
  //

  if ( cmdline_proc_psc || cmdline_proc_nmf )
    {
      if ( cmdline_proc_psc && cmdline_proc_nmf )
	Helper::halt( "cannot specify both --psc and --nmf" );
      
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );
            
      writer.begin();
      
      // writer.id( "." , "." );

      const std::string clab = cmdline_proc_psc ? "PSC" : "NMF" ; 
      
      writer.cmd( clab , 1 , "" );
      writer.level( clab, "_" + clab );

      // PSC , or NMF mode:
      psc_t psc;
      psc.construct( param , cmdline_proc_nmf );
    
      writer.unlevel( "_" + clab );

      writer.commit();
      std::exit(0);
    }


  //
  // POPS
  //
  
  if ( cmdline_proc_pops )
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

  if ( cmdline_proc_eval_stages ) 
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
  
  if ( cmdline_proc_pops_espriors )
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

  if ( cmdline_proc_lgbm )
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

  if ( cmdline_proc_assoc )
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
  
  if ( cmdline_proc_massoc )
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

  if ( cmdline_proc_fft )
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
  // OVERLAP enrichment (multi-sample case)
  //

  if ( cmdline_proc_overlap )
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

  if ( cmdline_proc_otsu )
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

  if ( cmdline_proc_cperm_test )
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

  if ( cmdline_proc_ms_kmer )
    {
      param_t param;
      build_param( &param, argc, argv, param_from_command_line );
            
      writer.begin();
      
      writer.id( "." , "." );      
      writer.cmd( "KMER" , 1 , "" );
      writer.level( "KMER", "_KMER" );

      std::string infile = Helper::expand( param.requires( "file" ) );
      int nreps = param.has( "nreps" ) ? param.requires_int( "nreps" ) : 1000;
      int k1 = param.has( "k1" ) ? param.requires_int( "k1" ) : 2;
      int k2 = param.has( "k2" ) ? param.requires_int( "k2" ) : 6;
      if ( param.has( "k" ) ) k1 = k2 = param.requires_int( "k" );
      // global versus local picks
      int w = param.has( "w" ) ? param.requires_int( "w" ) : 0 ; 
      
      // require at least L sequences; only take the first L
      const int req_len = param.has( "req-len" ) ? param.requires_int( "req-len" ) : 0 ; 

      
      // load from STDIN
      std::map<std::string,std::string> data0;
      std::vector<std::string> ids;
      if ( ! Helper::fileExists( infile ) ) Helper::halt( "could not open " + infile );
      std::ifstream IN1( infile.c_str() , std::ios::in );
      int rejected = 0; 
      while ( ! IN1.eof() )
	{
	  std::string id, s;
	  IN1 >> id >> s;	  
	  if ( IN1.eof() || id == "" || s.size() < 2 ) continue;

	  bool okay = req_len == 0 || s.size() >= req_len ; 	  
	  if ( ! okay ) { ++rejected; continue;  } 

	  // add, either whole sequence, or subset (1..s)
	  data0[ id ] = req_len ? s.substr( 0 , req_len ) : s;
	  ids.push_back( id );
	}
      IN1.close();

      if ( req_len ) 
	logger << "  " << data0.size() << " of " 
	       << data0.size() + rejected 
	       << " individuals included (analysis of first " << req_len << " states only)\n";

      //
      // Splice out '?' and ensure no similar sequences
      //
      
      std::map<std::string,std::string> data;
      std::map<std::string,std::string>::const_iterator ss = data0.begin();
      while ( ss != data0.end() )
	{
	  const std::string & s0 = ss->second;
	  std::vector<char> c;
	  const int n = s0.size();
	  char last = '?';
	  for (int i=0; i<n; i++)
	    {
	      if ( s0[i] == '?' ) continue;
	      if ( s0[i] == last ) continue;
	      c.push_back( s0[i] );
	      last = s0[i];
	    }
	  data[ ss->first ] = std::string( c.begin() , c.end() );

	  // std::cout << "orig = " << s0 << "\n"
	  // 	    << "new = " << data[ ss->first ] << "\n\n";
	  // next indiv.
	  ++ss;
	}
      
      //
      // report indiv-level enrichment? (versus group?)
      //

      const bool indiv_enrichment = param.has( "indiv-enrichment" );
      
      //
      // phenotypes?      
      //

      if ( param.has( "vars" ) )
	cmd_t::attach_ivars( param.value( "vars" ) );
      
      const std::string phe_label = param.has( "phe" ) ? param.value( "phe" ) : "" ;
      const bool grp = phe_label != "";

      if ( grp && indiv_enrichment ) 
	Helper::halt( "cannot specify both indiv-enrichment and phe" );

      std::map<std::string,int> phe;

      if ( grp )
	{
	  phe = cmd_t::pull_ivar( ids , phe_label );
	  int cases = 0 , controls = 0 , missing = 0;
	  
	  std::map<std::string,int>::const_iterator ii = phe.begin();
	  while ( ii != phe.end() )
	    {
	      if ( ii->second == 0 ) ++controls;
	      else if ( ii->second == 1 ) ++cases;
	      else ++missing;
	      ++ii;
	    }

	  logger << "  of " << ids.size() << " total individuals, for "
		 << phe_label << " "
		 << cases << " cases, "
		 << controls << " controls and "
		 << missing << " unknown\n"; 	    

	  if ( cases == 0 || controls == 0 )
	    Helper::halt( "did not observe both cases and controls: cannot run a phenotype-based analysis" );
	  
	}
      

      //
      // show within equivalence-group stats? (W_)
      //
      
      const bool wstats = param.has( "w-stats" );
      
      

      //
      // only show verbose mode for group level analysis 
      //

      const bool verbose_output = ! indiv_enrichment ; 

      //
      // Run analyses (either group level, in which case do once)
      //  or individual-level, in which case, we will iterate
      //  other all groups
      //

      std::map<std::string,std::string>::const_iterator ii = data.begin();

      while ( 1 ) 
	{
	  //
	  // copy over data for this analysis
	  //
	  
	  std::map<std::string,std::string> data1;
	  
	  if ( indiv_enrichment ) 
	    data1[ ii->first ] = ii->second;
	  else
	    data1 = data;

	  //
	  // do kmer enrichment: indiv, or group (w/ or w/out phenotype)
	  //

	  ms_kmer_t kmers( data1 , k1 , k2 , nreps , w, grp ? &phe : NULL , verbose_output );
      
	  
	  //
	  // individual level output?
	  //
	  
	  if ( indiv_enrichment )
	    {
	      logger << "  processing " << ii->first << ", L=" << ii->second.size() << " sequence length\n";
	      writer.id( ii->first , "." ); // ID, EDF
	      // track sequence length for each indiv
	      writer.value( "N" , (int)ii->second.size() );
	    }

	  //
	  // report output: OBS and WITHIN-GROUP
	  //
	  
	  std::map<std::string,double>::const_iterator pp = kmers.basic.pval.begin();
	  while ( pp != kmers.basic.pval.end() )
	    {
	      
	      writer.level( (int)pp->first.size() , "L" );
	      writer.level( pp->first , "S" );
	      
	      bool valid_equiv = kmers.equiv_set_size[ pp->first ] > 1;
	      
	      writer.value( "NG" , kmers.equiv_set_size[ pp->first ] );
	      writer.value( "SG" ,  kmers.obs2equiv[ pp->first ] );
	      
	      writer.value( "OBS" , kmers.basic.obs[ pp->first ] );
	      writer.value( "EXP" , kmers.basic.exp[ pp->first ] );
	      writer.value( "RAT" , kmers.basic.obs[ pp->first ] / (double) kmers.basic.exp[ pp->first ] );	      
	      writer.value( "P" , pp->second );
	      writer.value( "Z" , kmers.basic.zscr[ pp->first ] );	  
	      
	      if ( valid_equiv && wstats )
		{
		  writer.value( "W_OBS" , kmers.equiv.obs[ pp->first ] );
		  writer.value( "W_EXP" , kmers.equiv.exp[ pp->first ] );
		  writer.value( "W_RAT" , kmers.equiv.obs[ pp->first ] / (double)kmers.equiv.exp[ pp->first ]);
		  writer.value( "W_P" , kmers.equiv.pval[ pp->first ] );
		  writer.value( "W_Z" , kmers.equiv.zscr[ pp->first ] );
		}
	      
	      // C/C contrasts?
	      
	      if ( grp )
		{
		  
		  writer.level( "CASE" , "PHE" );
		  
		  writer.value( "NG" , kmers.equiv_set_size[ pp->first ] );
		  writer.value( "SG" ,  kmers.obs2equiv[ pp->first ] );
		  writer.value( "OBS" , kmers.basic_cases.obs[ pp->first ] );
		  writer.value( "EXP" , kmers.basic_cases.exp[ pp->first ] );
		  writer.value( "RAT" , kmers.basic_cases.obs[ pp->first ] /(double)kmers.basic_cases.exp[ pp->first ] );
		  writer.value( "P" , kmers.basic_cases.pval[ pp->first ] );
		  writer.value( "Z" , kmers.basic_cases.zscr[ pp->first ] );
		  
		  if ( valid_equiv && wstats )
		    {
		      writer.value( "W_OBS" , kmers.equiv_cases.obs[ pp->first ] );
		      writer.value( "W_EXP" , kmers.equiv_cases.exp[ pp->first ] );
		      writer.value( "W_RAT" , kmers.equiv_cases.obs[ pp->first ] /(double)kmers.equiv_cases.exp[ pp->first ] );
		      writer.value( "W_P" , kmers.equiv_cases.pval[ pp->first ] );
		      writer.value( "W_Z" , kmers.equiv_cases.zscr[ pp->first ] );
		    }
		  
		  writer.level( "CONTROL" , "PHE" );
		  
		  writer.value( "NG" , kmers.equiv_set_size[ pp->first ] );
		  writer.value( "SG" ,  kmers.obs2equiv[ pp->first ] );
		  writer.value( "OBS" , kmers.basic_controls.obs[ pp->first ] );
		  writer.value( "EXP" , kmers.basic_controls.exp[ pp->first ] );
		  writer.value( "RAT" , kmers.basic_controls.obs[ pp->first ] /(double)kmers.basic_controls.exp[ pp->first ] );
		  writer.value( "P" , kmers.basic_controls.pval[ pp->first ] );
		  writer.value( "Z" , kmers.basic_controls.zscr[ pp->first ] );
		  
		  if ( valid_equiv && wstats )
		    {
		      writer.value( "W_OBS" , kmers.equiv_controls.obs[ pp->first ] );
		      writer.value( "W_EXP" , kmers.equiv_controls.exp[ pp->first ] );
		      writer.value( "W_RAT" , kmers.equiv_controls.obs[ pp->first ] /(double)kmers.equiv_controls.exp[ pp->first ] );
		      writer.value( "W_P" , kmers.equiv_controls.pval[ pp->first ] );
		      writer.value( "W_Z" , kmers.equiv_controls.zscr[ pp->first ] );
		    }
		  
		  writer.level( "DIFF" , "PHE" );
		  
		  writer.value( "NG" , kmers.equiv_set_size[ pp->first ] );
		  writer.value( "SG" ,  kmers.obs2equiv[ pp->first ] );
		  //writer.value( "OBS" , kmers.basic_diffs.obs[ pp->first ] );
		  //writer.value( "EXP" , kmers.basic_diffs.exp[ pp->first ] );
		  writer.value( "Z" , kmers.basic_diffs.zscr[ pp->first ] );
		  
		  if ( valid_equiv && wstats )
		    {
		      // writer.value( "W_OBS" , kmers.equiv_diffs.obs[ pp->first ] );
		      // writer.value( "W_EXP" , kmers.equiv_diffs.exp[ pp->first ] );
		      writer.value( "W_Z" , kmers.equiv_diffs.zscr[ pp->first ] );
		    }
		  
		  writer.unlevel( "PHE" );
		}
	      
	      ++pp;
	    }  
	  
	  writer.unlevel( "S" );
	  writer.unlevel( "L" );
	  
	  //
	  // repeat for EQ groups
	  //
	  
	  pp = kmers.group.pval.begin();
	  while ( pp != kmers.group.pval.end() )
	    {
	      writer.level( (int)pp->first.size() , "L" );
	      writer.level( pp->first , "SG" );
	      
	      writer.value( "NG" , kmers.equiv_set_size[ pp->first ] );
	      
	      writer.value( "OBS" , kmers.group.obs[ pp->first ] );
	      writer.value( "EXP" , kmers.group.exp[ pp->first ] );
	      writer.value( "RAT" , kmers.group.obs[ pp->first ] / (double)kmers.group.exp[ pp->first ] );
	      writer.value( "P" , pp->second );
	      writer.value( "Z" , kmers.group.zscr[ pp->first ] );	  
	      
	      // C/C contrasts?
	      
	      if ( grp )
		{
		  
		  writer.level( "CASE" , "PHE" );
		  
		  writer.value( "NG" , kmers.equiv_set_size[ pp->first ] );	      
		  writer.value( "OBS" , kmers.group_cases.obs[ pp->first ] );
		  writer.value( "EXP" , kmers.group_cases.exp[ pp->first ] );
		  writer.value( "RAT" , kmers.group_cases.obs[ pp->first ] / (double) kmers.group_cases.exp[ pp->first ] );
		  writer.value( "P" , kmers.group_cases.pval[ pp->first ] );
		  writer.value( "Z" , kmers.group_cases.zscr[ pp->first ] );
		  
		  writer.level( "CONTROL" , "PHE" );
		  
		  writer.value( "NG" , kmers.equiv_set_size[ pp->first ] );
		  writer.value( "OBS" , kmers.group_controls.obs[ pp->first ] );
		  writer.value( "EXP" , kmers.group_controls.exp[ pp->first ] );
		  writer.value( "RAT" , kmers.group_controls.obs[ pp->first ] / (double) kmers.group_controls.exp[ pp->first ] );
		  writer.value( "P" , kmers.group_controls.pval[ pp->first ] );
		  writer.value( "Z" , kmers.group_controls.zscr[ pp->first ] );
		  
		  writer.level( "DIFF" , "PHE" );
		  
		  writer.value( "NG" , kmers.equiv_set_size[ pp->first ] );
		  // writer.value( "OBS" , kmers.group_diffs.obs[ pp->first ] );
		  // writer.value( "EXP" , kmers.group_diffs.exp[ pp->first ] );
		  writer.value( "Z" , kmers.group_diffs.zscr[ pp->first ] );
		  
		  writer.unlevel( "PHE" );
		}
	      
	      ++pp;
	    }        
	  
	  writer.unlevel( "SG" );
	  writer.unlevel( "L" );
	  
	  //
	  // if processing indiv-by-indiv, loop back
	  //
	  
	  if ( indiv_enrichment ) 
	    {
	      ++ii;
	      if ( ii == data.end() ) break;
	    }
	  else // if group mode, all done 
	    break;
	  
	}
      
      //
      // all done
      //
      
      writer.unlevel( "_KMER" );
      writer.commit();
      std::exit(0);
    }
  

  //
  // MS map label
  //
  
  if ( cmdline_proc_ms_label_maps )
    {

      param_t param;
      
      build_param( &param, argc, argv, param_from_command_line );
      
      writer.begin();
      
      writer.id( "." , "." );      
      writer.cmd( "LABEL-MAPS" , 1 , "" );
      writer.level( "LABEL-MAPS", "_LABEL-MAPS" );

      logger << " running LABEL-MAPS\n";

      //
      // options
      //

      bool verbose = param.has( "verbose" );

      // minimize sum(1-r)^p

      double p = param.has( "p" ) ? param.requires_dbl( "p" ) : 2 ;

      logger << "  minimizing sum_k (1-r)^" << p << "\n";
      
      //
      // Threshold of min spatial correl? 
      //

      const double th = param.has( "th" ) ? param.requires_dbl( "th" ) : 0 ; 

      if ( th < 0 || th > 1 ) Helper::halt( "invalid 'th' value - expecting 0 -- 1" );

      if ( th > 0 )
	logger << "  only assigning maps with spatial r >= " << th << " to matched template\n";
      else
	logger << "  no spatial correlation threshold ('th') set\n";
      
      //
      // Get template (with labels) and copy the (static) labels
      //
      
      ms_prototypes_t map_template;
      
      std::string template_map = Helper::expand( param.requires( "template" ) );	  
      
      map_template.read( template_map );
      
      std::vector<char> template_labels = map_template.ms_labels;

      //
      // Get to-be-labelled maps (updates ms_labels?)
      //
      
      ms_prototypes_t sol1;
      
      std::string sol1_file = Helper::expand( param.requires( "sol" ) );	  
      
      sol1.read( sol1_file );
      
      std::vector<char> sol1_labels = sol1.ms_labels;

      
      //
      // do mapping (based on maximal spatial correlation), updating static map labels [ to match sol1 ] 
      //  this will also edit 'sol1' to match polarity to closest to the template (for viz)
      //
      
      ms_prototypes_t::ms_labels = ms_cmp_maps_t::label_maps( map_template , template_labels , &sol1 , sol1_labels ,
							      th , p , verbose ); 
							      
      
      
      //
      // Re-write 'sol' (and updated labels will be includede)
      //

      std::string sol1_newfile = Helper::expand( param.requires( "new" ) );
      
      sol1.write( sol1_newfile );
      
      //
      // all done
      //
      
      writer.unlevel( "_LABEL-MAPS" );
      writer.commit();
      std::exit(0);
    }



  //
  // MS correl maps
  //

  if ( cmdline_proc_ms_corr_maps )
    {

      param_t param;
      
      build_param( &param, argc, argv, param_from_command_line );
      
      writer.begin();
      
      writer.id( "." , "." );      
      writer.cmd( "CORREL-MAPS" , 1 , "" );
      writer.level( "CORREL-MAPS", "_CORREL-MAPS" );

      logger << " running CORREL-MAPS\n";
      
      
      //
      // Get to-be-labelled maps (updates ms_labels?)
      //
      
      ms_prototypes_t A;
      
      std::string sol_file = Helper::expand( param.requires( "sol" ) );	  
      
      A.read( sol_file );
      
      //
      // Spatial correlations
      //
      
      const int nk = A.K;
      
      Eigen::MatrixXd R = Eigen::MatrixXd::Zero( nk , nk );
      for (int i=0; i<nk; i++)
	for (int j=0; j<nk; j++)
	  R(i,j) = ms_prototypes_t::spatial_correlation( A.A.col(i) , A.A.col(j) );
      
      for (int i=0; i<nk; i++)
	std::cout << "\t" << A.ms_labels[i];
      std::cout << "\n";
      
      for (int i=0; i<nk; i++)
	{
	  std::cout << A.ms_labels[i];
	  for (int j=0; j<nk; j++) 
	    std::cout << "\t" << R(i,j);
	  std::cout << "\n";	  
	}
      
      //
      // all done
      //
      
      writer.unlevel( "_CORREL-MAPS" );
      writer.commit();
      std::exit(0);
    }


  //
  // MS compare maps
  //

  if ( cmdline_proc_ms_cmp_maps )
    {

      param_t param;

      build_param( &param, argc, argv, param_from_command_line );
      
      writer.begin();
      
      writer.id( "." , "." );      
      writer.cmd( "CMP-MAPS" , 1 , "" );
      writer.level( "CMP-MAPS", "_CMP-MAPS" );

      logger << " running CMP-MAPS\n";
      
      //
      // number of permutations to perform 
      //

      const int nreps = param.has( "nreps" ) ? param.requires_int( "nreps" ) : 1000;

      //
      // to define global similarity: greedy or brute-force (default) enumeration of all possibilities?
      //
      
      // minimize sum(1-r)^p                                                                                                                                                                               
      const double p = param.has( "p" ) ? param.requires_dbl( "p" ) : 2 ;
      logger << "  matching based on minimizing sum_k (1-r)^" << p << "\n";
      
      //
      // Either all-case compared to all-controls : stat = d( concordant pairs ) / d( discordant pairs )
      // OR given a fixed map=M: stat =  ( d( case - X ) - d( control - X )^2 
      // 

      const bool use_fixed = param.has( "template" );
      
      ms_prototypes_t fixed;
      
      if ( use_fixed )
	{
	  // read a standard prototype map file (sol format, i.e. no ID)
	  std::string fixed_map = Helper::expand( param.value( "template" ) );	  
	  fixed.read( fixed_map );
	}
      
      //
      // Load maps
      //

      // expect a file as output from MS A matrix in long format
      // ID	CH	K	A
      
      std::string infile = Helper::expand( param.requires( "file" ) );

      // ID -> K -> CH -> 'A'

      std::map<std::string,std::map<std::string,std::map<std::string,double> > > data;
      if ( ! Helper::fileExists( infile ) ) Helper::halt( "could not open " + infile );
      std::ifstream IN1( infile.c_str() , std::ios::in );
      // header...
      std::string id, ch, k, dummy;
      IN1 >> id >> ch >> k >> dummy;
      if ( id != "ID" || ch != "CH" || k != "K" || dummy != "A" )
	Helper::halt( "bad format" );
      
      while ( ! IN1.eof() )
	{
	  std::string id, ch, k;
	  double a;
	  IN1 >> id >> ch >> k >> a;
	  if ( IN1.eof() || id == "" ) continue;
	  data[ id ][ k ][ ch ] = a;	  
	}
      IN1.close();

      
      //
      // phenotypes?      
      //

      if ( param.has( "vars" ) )
	cmd_t::attach_ivars( param.value( "vars" ) );
      
      const std::string phe_label = param.has( "phe" ) ? param.value( "phe" ) : "" ;
      const bool grp = phe_label != "";
      
      std::map<std::string,int> phe;
      
      if ( grp )
	{
	  std::vector<std::string> ids;
	  std::map<std::string,std::map<std::string,std::map<std::string,double> > >::const_iterator qq =  data.begin();
	  while ( qq != data.end() ) { ids.push_back( qq->first ); ++qq; } 
		
	  phe = cmd_t::pull_ivar( ids , phe_label );
	  int cases = 0 , controls = 0 , missing = 0;
	  
	  std::map<std::string,int>::const_iterator ii = phe.begin();
	  while ( ii != phe.end() )
	    {
	      if ( ii->second == 0 ) ++controls;
	      else if ( ii->second == 1 ) ++cases;
	      else ++missing;
	      ++ii;
	    }

	  logger << "  of " << data.size() << " total individuals, for "
		 << phe_label << " "
		 << cases << " cases, "
		 << controls << " controls and "
		 << missing << " unknown\n";

	  if ( cases == 0 || controls == 0 )
	    Helper::halt( "did not observe both cases and controls: cannot run a phenotype-based analysis" );
	      
	}

      //
      // do analysis (& writes output too)
      //
      
      ms_cmp_maps_t cmp_maps( data ,
			      use_fixed ? &(fixed.A) : NULL ,
			      use_fixed ? &(fixed.chs) : NULL ,
			      phe ,
			      nreps ,
			      p );
      
      //
      // all done
      //
      
      writer.unlevel( "_CMP-MAPS" );
      writer.commit();
      std::exit(0);
    }



  //
  // FIR design
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

      writer.cmd( "CWT-DESIGN" , 1 , "" );
      writer.level( "CWT-DESIGN", "_CWT-DESIGN" );

      // expects input from std::cin
      proc_cwt_design_cmdline();

      writer.commit();

      std::exit(0);

    }

  
  if ( cmdline_proc_copy_suds )
    {
      // expects input from std::cin                                                                                                                                                               
      proc_copy_suds_cmdline();
      std::exit(0);
    }


  if ( cmdline_proc_combine_suds )
    {
      // expects input from std::cin                                                                                                                                                              
      proc_combine_suds_cmdline();
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

  // use .edf (or .EDF extension) to indicate 'single EDF' mode, '.rec'
  f = f.substr( (int)f.size() - 4 >= 0 ? (int)f.size() - 4 : 0 );
  bool single_edf = Helper::iequals( f , ".edf" ) || Helper::iequals( f , ".rec" ) ;
  
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

	  // remove .edf from ID, making file name ==> ID 
	  if ( Helper::file_extension( rootname , "edf" ) )
	    rootname = rootname.substr( 0 , rootname.size() - 4 );

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
      else if ( empty_edf )
	{	  
	  const int nr = globals::param.requires_int( "-nr" );
	  const int rs = globals::param.requires_int( "-rs" ); // in full seconds (integer)
	  const std::string startdate = globals::param.has("-date") ? globals::param.value( "-date" ) : "01.01.00" ;
          const std::string starttime = globals::param.has("-time") ? globals::param.value( "-time" ) : "00.00.00" ;
          const std::string id = globals::param.has("-id") ? globals::param.value( "-id" ) : rootname ;
	  okay = edf.init_empty( id , nr , rs , startdate , starttime );
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
      // Give annotations some basic details about the EDF
      //

      edf.timeline.annotations.set( &edf );
      
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
	    edf.timeline.annotations.from_EDF( edf , edf.edfz_ptr() );
	  else if ( ! edf.header.continuous )
	    edf.timeline.annotations.from_EDF( edf , edf.edfz_ptr() );
	  
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

void proc_dummy( const std::string & p , const std::string & p2 )
{

  if ( p == "randomize-kmer" )
    {

      // e.g.
      // awk ' { print $2 } ' seq.1 | sed -e 's/\(.\)/\1\'$'\n/g' | awk ' NF>0 ' | luna -d randomize-kmer > seq.2
      // from an existing sequence.

      // echo "file=seq.1 k=4 nreps=100 w=5" | luna --kmer -o out-local500-v1.db
      // echo "file=seq.2 k=4 nreps=100 w=5" | luna --kmer -o out-local500-v2.db 
      
      std::vector<char> s;
      std::map<char,int> u;
      
      while ( 1 )
	{
	  std::string c;
	  std::cin >> c;
	  if ( c == "" || std::cin.eof() ) break;
	  if ( c.size() != 1 ) break;
	  s.push_back( c[0] );
	  ++u[c[0]];
	}

      const int n = s.size();
      
      std::cerr << " read " << n << " elements\n";
      std::string s1( n , '.' );
      for (int i=0; i<n; i++)
	s1[i] = s[i];
	  
      std::map<char,int>::const_iterator uu = u.begin();
      while ( uu != u.end() )
	{
	  std::cerr << " " << uu->first << " = " << uu->second << "\n";
	  ++uu;
	}

      ms_kmer_t ms1;

      int w = 0;
      if ( p2 != "" )
        if ( ! Helper::str2int( p2 , &w ) )
	  Helper::halt( "expecting integer w as second parameter" );

      std::cerr << " w = " << w << "\n";

      std::string s2 = ms1.modified_random_draw( s1 , w );
      
      std::cout << "ID1\t" << s2 << "\n";
      
      std::exit(0);
    }

  if ( p == "cgi" )
    {
      std::string res = exec_system( "ls -l" );

      std::cout << " my result\n"
		<< "----------\n"
		<< res
		<< "\n-----------\n";
      
      std::exit(0);
    }
  
  if ( p == "peaks" )
    {
      std::vector<double> x(100);
      for (int i=0;i<100;i++) x[i] = i*i;
      x[9] += 2000;
      x[20] += 2000;
      double m1, m2;
      std::vector<double> s1, s2, s3;
      psd_shape_metrics( x , x , 5 , &m1, &m2, &s1, &s2, &s3 );
      for (int i=0;i<s1.size();i++)
	std::cout << x[i] << "\t"
		  << s1[i] << "\t"
		  << s2[i] << "\t"
		  << s3[i] << "\n";

      std::cout << "m1\t" << m1 << "\n"
		<< "m2\t" << m2 << "\n";
     
      std::exit(0);
    }
  
  if ( p == "kmeans" )
    {
      const int nc = 4;
      const int nr = 150;
      const int nk = 3;
      Data::Matrix<double> X( nr , nc );
      for (int r=0;r<nr; r++)
	for (int c=0;c<nc; c++)
	  std::cin >> X(r,c);
      
      std::cout << "X. " << X.print() << "\n";

      kmeans_t kmeans;
      
      std::vector<int> sol;
      kmeans.kmeans( X , 3 , &sol );
      std::cout << "SOL\n";
      for (int i=0;i<150;i++)
	std::cout << sol[i] << "\n";

      std::exit(1);
    }
  
  if ( p == "json" )
    {

      // store a string in a JSON value
      nlohmann::json j_string = "this is a string";
      
      // retrieve the string value
      auto cpp_string = j_string.get<std::string>();
      // retrieve the string value (alternative when an variable already exists)
      std::string cpp_string2;
      j_string.get_to(cpp_string2);
      
      // retrieve the serialized value (explicit JSON serialization)
      std::string serialized_string = j_string.dump();
      
      // output of original string
      std::cout << cpp_string << " == " << cpp_string2 << " == " << j_string.get<std::string>() << '\n';
      // output of serialized value
      std::cout << j_string << " == " << serialized_string << std::endl;
      
      
      std::exit(1);
      
    }
  
  if ( p == "runs" )
    {
      std::vector<std::string> d = { "S", "S", "S", "F", "S", "F", "F", "F", "F" };


      // runs_test({'M':{1:29, 2:10, 3:8, 4:3, 5:1, 6:1, 7:1, 8:1, 12:2}, 'D':{1:33, 2:17, 3:6, 5:1}}, path=False)

      
      //std::vector<std::string> d = { "S", "S", "F", "S", "F", "F", "F", "F", "S", "S", "S", "F", "F" };
      //std::vector<std::string> d = { "S", "S", "S", "F", "S", "F", "S", "S", "S", "F", "F", "S", "S", "S" };
      std::cout << "runs p = " << Statistics::runs_test( d ) << "\n";
      std::exit(1);
// #x = ["S", "S", "S", "F", "S", "F", "S", "S", "S", "F", "F", "S", "S", "S"]
// #x = ["B","B","A","C","C","A","C","C","C","A","B","A","A","A","B","A","A","B","B","A","B","A","A","B","A","B","B"]
// #x = ["S","S","S","S","S","F", "F", "F", "F", "F","S","S","S","S","S","F", "F", "F", "F", "F", "S", "F"]
// #recalculation of results in O'Brien
// runs_test({'M':{1:29, 2:10, 3:8, 4:3, 5:1, 6:1, 7:1, 8:1, 12:2}, 'D':{1:33, 2:17, 3:6, 5:1}}, path=False)
// #this will produce an exception
// #runs_test(["S", "S", "S", "F", "S", "F", "S", "S", "S", "F", "F", "S", "S", "S"])
// runs_test(["B","B","A","C","C","A","C","C","C","A","B","A","A","A","B","A","A","B","B","A","B","A","A","B","A","B","B"])
//runs_test(["S", "S", "S", "F", "S", "F", "F", "F", "F"])
    }
     

  if ( p == "cwt" )
    {
      CWT cwt;
      cwt.set_sampling_rate( 400 );
      cwt.add_wavelets( 0.5 , 5 , 30 , 0.25 , 35 , 20 ) ;


      std::vector<dcomp> w1 = cwt.alt_wavelet( 0 );
      for (int i=0;i<w1.size();i++)
	std::cout << i << "\t" << w1[i] << "\n";
      
      std::exit(1);
    }

  if ( p == "cancor" )
    {
      const int nrows = 100;
      const int nvars = 10;

      Eigen::MatrixXd X( nrows , nvars );
      Eigen::MatrixXd Y( nrows , nvars );

      int i = 0 , j = 0;
      std::ifstream INX( Helper::expand( "~/x.txt" ).c_str() , std::ios::in );
      while ( ! INX.eof() )
        {
          double d;
          INX >> d;
          if ( INX.eof() ) break;
          X(i,j) = d;
	  ++j;
	  if ( j == nvars ) { ++i; j=0; }
	}
      INX.close();
      
      i = j = 0;
      std::ifstream INY( Helper::expand( "~/y.txt" ).c_str() , std::ios::in );
      while ( ! INY.eof() )
        {
          double d;
          INY >> d;
          if ( INY.eof() ) break;
          Y(i,j) = d;
	  ++j;
	  if ( j == nvars ) { ++i; j=0; }
	}
      INY.close();

      // std::cout << "X\n" << X << "\n\n";
      // std::cout << "Y\n" << Y << "\n\n";
      
      Eigen::VectorXd CCA = eigen_ops::canonical_correlation( X , Y );
      
      std::cout << " CCA \n"
       		<< CCA << "\n";
      
      std::exit(1);
    }
  
  if ( p == "qda" )
    {
      std::vector<std::string> y;
      const int nrows = 1257;
      const int nvars = 18;
      
      Eigen::MatrixXd X( nrows , nvars );
      
      std::ifstream INY( Helper::expand("~/y.txt" ).c_str() , std::ios::in );
      int k = 0;
      while ( ! INY.eof() )
	{
	  std::string s;
	  INY >> s;
	  if ( INY.eof() ) break;
	  y.push_back(s);
	}

      INY.close();
      int i = 0 , j = 0;
      std::ifstream INX( Helper::expand( "~/x.txt" ).c_str() , std::ios::in );
      while ( ! INX.eof() )
        {
          double d;
          INX >> d;
          if ( INX.eof() ) break;
          X(i,j) = d;
	  ++j;
	  if ( j == nvars ) { ++i; j=0; }
	}
      INX.close();
      std::cerr << "done reading\n";

      std::cout << " set up ..\n";
      qda_t qda( y , X );

      std::cout << " fitting...\n";
      
      qda_model_t fit = qda.fit();

      std::cout << " predictiong...\n";

      qda_posteriors_t pp = qda.predict( fit , X );

      for (int i=0;i<pp.pp.rows() ;i++)
	{
	  for (int j=0;j<pp.pp.cols();j++)
	    std::cout << " " << pp.pp(i,j);
	  std::cout << "\t" << pp.cl[i] << "\n";

	}      
      
      std::exit(1);
      
    }

  
  if ( p == "lda" )
    {
      std::vector<std::string> y;

      Eigen::MatrixXd X(500,10);
      
      std::ifstream INY( Helper::expand("~/y.txt" ).c_str() , std::ios::in );
      int k = 0;
      while ( ! INY.eof() )
	{
	  std::string s;
	  INY >> s;
	  if ( INY.eof() ) break;
	  y.push_back(s);
	}

      INY.close();
      int i = 0 , j = 0;
      std::ifstream INX( Helper::expand( "~/x.txt" ).c_str() , std::ios::in );
      while ( ! INX.eof() )
        {
          double d;
          INX >> d;
          if ( INX.eof() ) break;
          X(i,j) = d;
	  ++j;
	  if ( j == 10 ) { ++i; j=0; }
	}
      INX.close();
      std::cerr << "done reading\n";
      
      lda_t lda( y , X );

      lda_model_t fit = lda.fit();

      lda_posteriors_t pp = lda.predict( fit , X );

      for (int i=0;i<pp.pp.rows() ;i++)
	{
	  for (int j=0;j<pp.pp.cols();j++)
	    std::cout << " " << pp.pp(i,j);
	  std::cout << "\t" << pp.cl[i] << "\n";

	}      
      
      std::exit(1);
      
    }

  if ( p == "cache" )
    {
      ctest();
      std::exit(0);
    }

  if ( p == "kmer" )
    {
      std::vector<int> x;
      while ( ! std::cin.eof() )
	{
	  int i;
	  std::cin >> i;
	  if ( std::cin.eof() ) break;
	  x.push_back(i);
	}

      ms_kmer_t kmers( x , 2 , 6 , 1000 , 0 );
      
      std::map<std::string,double>::const_iterator pp = kmers.basic.pval.begin();
      while ( pp != kmers.basic.pval.end() )
	{
	  std::cout << pp->first << "\t"
		    << pp->first.size() << "\t"
		    << kmers.basic.obs[ pp->first ] << "\t"
		    << pp->second << "\n";
	  ++pp;
	}
      
      std::exit(1);
    }
  
  if ( p == "cmddefs" ) 
    {
      
      globals::cmddefs().add_domain( "misc" , "misc" ,  "Misc" );
      
      globals::cmddefs().add_cmd( "misc" , "comm1" , "this is a dummy command" );
      globals::cmddefs().add_table( "comm1" , "XX" , "A simple table" , false );
      globals::cmddefs().add_var( "comm1" , "XX" , "X", "Var X" );
      globals::cmddefs().add_var( "comm1" , "XX" , "Y" , "Var Y" );
      
      globals::cmddefs().add_cmd( "misc" , "comm2" , "this is a dummy command" );
      globals::cmddefs().add_table( "comm2" , "CH,B" , "A nice table" , true );
      globals::cmddefs().add_var( "comm2" , "CH,B" , "V1" , "Variable 1" );
      globals::cmddefs().add_var( "comm2" , "CH,B" , "V2" , "Variable 2" );
      
      //   std::cout << globals::cmddefs().help( "comm1" , true )  << "\n\n\n";
      std::cout << globals::cmddefs().help( "comm2" , true )  << "\n";
      
      // add a dummy tag
      globals::cmddefs().add_tag( "Z" );
      
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
      
      globals::cmddefs().add_cmd( "misc"   , "NEWONE" , "A test command" );
      globals::cmddefs().add_table( "NEWONE" , "" , "Table 0, baseline" );
      globals::cmddefs().add_table( "NEWONE" , "CH" , "Table 1, by channel" );
      globals::cmddefs().add_table( "NEWONE" , "CH,X" , "Table 2, by channel and X" );
      globals::cmddefs().add_table( "NEWONE" , "CH,X,Y" , "Table 2a, by channel and X/Y" , true );
      globals::cmddefs().add_table( "NEWONE" , "CH,X,Z" , "Table 2b, by channel and X/Z"  );

      globals::cmddefs().add_var( "NEWONE" , "" , "V1" , "some var1" );
      globals::cmddefs().add_var( "NEWONE" , "" , "V2" , "some var2" );

      globals::cmddefs().add_var( "NEWONE" , "CH" , "V1" , "some var1" );
      globals::cmddefs().add_var( "NEWONE" , "CH" , "V2" , "some var2" );
      globals::cmddefs().add_var( "NEWONE" , "CH" , "V3" , "some var3" );

      globals::cmddefs().add_var( "NEWONE" , "CH,X" , "V2a" , "some var2" );
      globals::cmddefs().add_var( "NEWONE" , "CH,X" , "V3a" , "some var3" );

      globals::cmddefs().add_var( "NEWONE" , "CH,X,Y" , "V2a" , "some var2" );
      globals::cmddefs().add_var( "NEWONE" , "CH,X,Y" , "V3a" , "some var3" );

      globals::cmddefs().add_var( "NEWONE" , "CH,X,Z" , "V2a" , "some var2" );
      globals::cmddefs().add_var( "NEWONE" , "CH,X,Z" , "V3a" , "some var3" );

      //  std::cout << globals::cmddefs().help_domains() << "\n";
      //std::cout << globals::cmddefs().help_commands() << "\n";
      
      globals::cmddefs().set_compressed( "NEWONE" , tfac_t( "CH,X,Z" ) );
      globals::cmddefs().set_compressed( "NEWONE" , tfac_t( "CH,X,Y" ) , false );
      
      std::cout << globals::cmddefs().help( "NEWONE" , true ) << "\n";
      
      
      std::exit(0);

      param_t param;
      param.add( "epoch" );
      param.add( "ep" );
      std::set<std::string> unk;
      std::cout << "st = " << globals::cmddefs().check( "ANNOTS" , param.keys() , &unk ) << "\n";
      
      std::set<std::string>::const_iterator uu = unk.begin();
      while ( uu != unk.end() ) { std::cout << " bad param " << *uu << "\n"; ++uu; } 



      
      

    }


  //
  // LightGBM
  //

  if ( p == "lgbm" )
    {
#ifdef HAS_LGBM
      
      lgbm_t lgbm( "train.conf" ) ;
      
      lgbm.load_training_data( "binary.train" );

      const int n1 = lgbm_t::rows( lgbm.training );
      const int n2 = lgbm_t::cols( lgbm.training );

      lgbm_label_t labels( "luna.wgt" );
      std::cout << " from luna.wgt " << labels.n << "\n";

      lgbm.add_label_weights( lgbm.training , &lgbm.training_weights, labels );

      lgbm.apply_weights( lgbm.training , &lgbm.training_weights );
      
      //lgbm.load_weights( lgbm.training , "RENAMED_binary.train.weight" );

      std::vector<int> l = lgbm_t::labels( lgbm.training );
      std::vector<double> w = lgbm_t::weights( lgbm.training );
      
      std::cout << " l = " << l.size() << " ... \n";
      for (int i=0; i<30; i++) std::cout << l[i] << "\t" << w[i] << "\n";
                 
      lgbm.load_validation_data( "binary.test" );
      
      lgbm.create_booster();
      
      lgbm.save_model( "my-model.1" );      
      
#endif
      std::exit(0);
    }



  if ( p == "lgbm2" )
    {
#ifdef HAS_LGBM

      Eigen::MatrixXd X = eigen_ops::load_mat( "binary.test" );

      // remove first col
      X = X.rightCols( X.cols() - 1 );

      lgbm_t lgbm;

      lgbm.load_model( "my-model.1" );
      //lgbm.load_model( "LightGBM_model.txt" );
      //lgbm.load_model( "R.bst" );

      lgbm.predict( X );
      
#endif
      std::exit(0);
    }


  //
  // test date/time functions
  //

  if ( p == "datetime" )
    {
      
      if ( 0 )
	{
	  for (int c=0; c<10000; c++)
	    {
	      std::string ds = date_t::datestring( c );
	      date_t dt( ds );
	      int c2 = date_t::count( dt );
	      std::cout << ( c != c2 ? "*****" : "" ) 
			<< c << "\t"
			<< c2 << "\t"
			<< ds << "\n";
	    }
	}

      if ( 1 )
	{
	  std::string inp1, inp2;
	  std::cin >> inp1 >> inp2;
	  clocktime_t t1( inp1 );
	  clocktime_t t2( inp2 );

	  std::cout << "t1: "
		    << t1.valid << "\t"
		    << t1.as_string( ) << "\t"
		    << t1.as_datetime_string( ) << "\n";
	  //<< t1.d << " - " << t1.h << ":" << t1.m << ":" << t1.s << "\n";
	  
	  std::cout << "t2: "
		    << t2.valid << "\t"
		    << t2.as_string( ) << "\t"
		    << t2.as_datetime_string( ) << "\n";
	  //<< t2.d << " - " << t2.h << ":" << t2.m << ":" << t2.s << "\n";

	  const int earlier = clocktime_t::earlier( t1 , t2 ) ; 
	  
	  const double difh = earlier == 0 ? 0 : ( earlier == 1 ? clocktime_t::difference_hours( t1 , t2 ) : clocktime_t::difference_hours( t2 , t1 )   );
	  const double difs = earlier == 0 ? 0 : ( earlier == 1 ? clocktime_t::difference_seconds( t1 , t2 ) : clocktime_t::difference_seconds( t2 , t1 ) );

	  clocktime_t midpoint;
	  midpoint.midpoint( t1 , t2 );
	  
 	  std::cout << " earlier = " << earlier << "\n";
	  std::cout << " t1 - t2 (hours) = " << difh << "\n";
	  std::cout << " t1 - t2 (secs) = " << difs << "\t" << difs / 3600.0 << "\n";
	  std::cout << " midpoint = " << midpoint.as_datetime_string() << "\n";
	  std::cout << "\n";
	  clocktime_t nt = t1;
	  for (int i=0;i<48;i++)
	    {
	      //nt.advance_hrs( 1.222 );
	      nt.advance( clocktime_t( "+1:30" ) );
	      std::cout << "  --> " << nt.as_string() << "\t" << nt.as_datetime_string() << "\n";
	    }
	}
      
      
      
      std::exit(0);
    }
  

  //
  // TRANS
  //

  if ( p == "trans" )
    {      
      std::string line;
      Helper::safe_getline( std::cin , line );
      std::vector<std::string> hdr = Helper::parse( line );
      const int k = hdr.size();
      
      std::cerr << " expr [" << p2 << "]\n";
      
      std::map<std::string,std::vector<double> > inputs;

      int rows = 0;
      
      while ( 1 )
	{
	  std::string line;
	  Helper::safe_getline( std::cin , line );
	  if ( line == "" ) break;
	  if ( std::cin.eof() || std::cin.bad() ) break;
	  std::vector<std::string> tok = Helper::parse( line );
	  if ( tok.size() != k ) Helper::halt( "wrong numbr of columns" );
	  for (int i=0; i<k; i++)
	    {
	      double d;
	      if ( ! Helper::str2dbl( tok[i] , &d ) ) Helper::halt( "bad numeric value" );
	      inputs[hdr[i]].push_back( d );
	    }
	  ++rows;
	}
      
      std::cerr << "read " << rows << "\n";
      
      // output
      instance_t out;
      
      // expression
      Eval expr( p2 );
      
      // bind input/output data to token evaluator
      expr.bind( inputs , &out );

      // evaluate
      bool is_valid = expr.evaluate( );

      // returned a valid bool? (single value)
      bool retval;      
      bool is_valid_retval = true;
      if ( ! expr.value( retval ) ) is_valid_retval = false;

      std::cerr << "parsed as a valid expression : " << ( is_valid ? "yes" : "no" ) << "\n";
      if ( is_valid_retval ) 
	std::cerr << "boolean return value         : " << ( retval ? "true" : "false" ) << "\n";
      std::cerr << "assigned meta-data           : " << out.print() << "\n";  

      //
      // actual output
      //
      
      std::vector<double> rr = expr.value().as_float_vector();
      for (int i=0; i<rr.size();i++)
	std::cout << rr[i] << "\n";

      std::exit(1);
    }
  
  //
  // Straight FFT of stdin
  //

  std::vector<double> x;
  

  if ( p == "fir" || p == "fft" || p == "dfa" || p == "fft-test" || p == "mtm" || p == "tv" || p == "psi" 
       || p == "dynam" || p == "ica" || p == "robust" || p == "fip" || p == "sl" || p == "acf" || p == "otsu"
       || p == "desats" || p == "zpks" || p == "gc" || p == "detrend" || p == "emd" || p == "tri" ) 
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

  if ( p == "desats" )
    {
      hb_find_desats_t r = hb_t::find_desats( eigen_ops::copy_array( x ) , 32 , 1.5 );

      std::cout << r.MagDown << "\n\n";

      std::cout << r.dsatStEnd << "\n\n";

      std::exit(1);
    }

  if ( p == "emd" )
    {

      emd_t emd;
      
      const int nk = emd.proc( &x );
      const int nr = emd.residual.size();
      
      for (int i=0; i<nr; i++)
	{
	  std::cout << x[i] ;
	  for (int j=0; j<nk; j++)
	    std::cout << "\t" << emd.imf[j][i] ;
	  std::cout << "\t" << emd.residual[i] << "\n";
	}
      
      std::exit(1);
    }

  if ( p == "detrend" )
    {
      const int n = x.size();
      
      std::vector<double> x2 = x;
      double beta , intercept;
      MiscMath::detrend( &x2 , &intercept , &beta );
      
      std::cout << "m, b = " << intercept << " " << beta << "\n";
      for (int i=0; i<n; i++)
	std::cout << x[i] << "\t" << x2[i] << "\n";

      
      Eigen::VectorXd T( n );
      for (int i=0; i<n; i++) T[i] = x[i];
      
      std::cout << "EIGEN\n orig T\n" << T << "\n";

      std::cout << "DT\n";
      eigen_ops::detrend( T ) ;
      std::cout << T << "\n";
      
    }

  if ( p == "gc" )
    {
      
      int order = 3;

      if ( p2 != "" )
        if ( ! Helper::str2int( p2 , &order ) ) Helper::halt( "expecting integer model order as second parameter" );

      Eigen::MatrixXd X( x.size() / 2 , 2 );
      int cnt = 0;
      for (int r=0; r< x.size() / 2 ; r++)
	{
	  X(r,0) = x[cnt++] ;
	  X(r,1) = x[cnt++] ; 
	}

      std::cerr << "read " << x.size() / 2 << " observations (pairs)\n";

      // fix
      const int Nr = 99;
      const int Nl = 51;
      
      order = 15; 

      // armorf_t ar1( X.col(0) , Nr , Nl , order );
      
      // std::cout << "ar1.coeff\n" << ar1.coeff << "\n"
      // 		<< "ar1.E\n" << ar1.E << "\n";

      // armorf_t ar2( X.col(1) , Nr , Nl , order );
      
      // std::cout << "ar2.coeff\n" << ar2.coeff << "\n"
      // 		<< "ar2.E\n" << ar2.E << "\n";

      // armorf_t ar12( X , Nr , Nl , order );

      // std::cout << "ar12.coeff\n" << ar12.coeff << "\n"
      // 		<< "ar12.E\n" << ar12.E << "\n";

      signal_list_t signals;
      signals.add( 0, "S1" );
      signals.add( 1, "S2" );
      const int sr = 256;
      // std::vector<double> frqs(3);
      // frqs[0] = 1; frqs[30] ; frqs[2] = 10;

      std::vector<double> frqs = MiscMath::logspace( 10 , 40 , 15 );

      gc_t gc( X , signals , sr , 200 , 60 , &frqs );

      gc.report( signals );
      std::exit(0);

    }


  if ( p == "zpks" ) 
    {
      std::vector<interval_t> ints;
      //std::vector<int> s = MiscMath::smoothedZ( x , 30*256 , 3 , 0 , 128 , 20 , 2 , 256 , true , &ints) ;
      std::vector<int> s = MiscMath::smoothedZ( x , 400 , 3 , 0 , 96 , 0 , 0 , 0 , true , &ints , true ) ;
      
      for (int i=0; i<ints.size(); i++)
	std::cout << i << "\t" << ints[i].start << " -- " << ints[i].stop << "  " << ints[i].stop - ints[i].start << " " << ( ints[i].stop - ints[i].start ) / 256.0 << "\n";

      // std::vector<int> smoothedZ( const std::vector<double> & x , int lag , double threshold , double influence = 0 , 
      // 				  int mindur = 0 , double max = 0 , 
      // 				  double threshold2 = 0 , int mindur2 = 0 , 
      // 				  bool noneg = false , 
      // 				  std::vector<interval_t> * regions = NULL );
      

      //std::vector<int> s = MiscMath::smoothedZ( x , 30 , 5 ) ;
      // for (int i=0; i<s.size(); i++)
      // 	std::cout << x[i] << "\t" << s[i] << "\n";
       std::exit(1);
    }
 
  if ( p == "psi" )
    {
      const int n = x.size() / 2 ;
      Data::Matrix<double> data( n , 2 );
      int r = 0;
      for (int i=0;i<n;i++)
	{
	  data(i,0) = x[ r++ ];
	  data(i,1) = x[ r++ ];
	}

      psi_t psi( &data , 100 , 200 , 200 );

      psi.calc();

      signal_list_t signals;
      signals.add( 0, "S1" );
      signals.add( 1, "S2" );
      
      psi.report( signals );
	    
      std::exit(0);
    }
  
  if ( p == "robust" )
    {
      
      const int n = x.size();
      Eigen::MatrixXd m( n , 1 );
      for (int i=0;i<n;i++) m(i,0) = x[i];
  
      eigen_ops::robust_scale( m , true , true , 0.05 );
      std::cout << "\n" << m  << "\n";
      std::exit(0);
    }

  if ( p == "otsu" )
    {
      std::map<double,double> tvals, fvals;      
      double f;
      double th = MiscMath::threshold2( x , &f, 0 , &fvals , &tvals );
      
      std::cout << "best th = " << th << "\n";

      std::map<double,double>::const_iterator ii = tvals.begin();
      while ( ii != tvals.end() )
       	{
       	  std::cout << "th = " << ii->first << "\t varB = " << ii->second << "\t F = " << fvals[ ii->first ] << "\n";
       	  ++ii;
       	}
      
      std::exit(0);
    }


  if ( p == "acf" )
    {
      acf_t acf( x );
      std::vector<double> rr = acf.acf();
      for (int i=0;i<rr.size();i++)
	std::cout << "lag = " << i << "\t" << rr[i] << "\n";
      std::exit(0);
    }


  if ( p == "anova" )
    {
      std::vector<std::string> group;
      Data::Vector<double> x;
      while ( ! std::cin.eof() ) 
	{
	  std::string g;
	  double t;
	  std::cin >> g >> t;
	  if ( std::cin.eof() ) break;
	  group.push_back( g );
	  x.push_back( t );
	}

      std::cout << Statistics::anova( group , x );
      std::exit(0);
      
    }
  
  if ( p == "fip" )
    {


      int sr = 256;

      if ( p2 != "" )
        if ( ! Helper::str2int( p2 , &sr ) ) Helper::halt( "expecting integer sample rate as second parameter" );

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


  if ( p == "fft-test" )
    {

      // test 1 : equivalence w/ mt_get_spec() and real_FFT
      
      //mtm::mt_get_spec ( b, npoints, klen, amp);  
      int fs = 256;
      double dt = 1/(double)fs;
      double * series = &(x)[0];
      int inum = x.size();
      int npoints = inum;
      int klen = MiscMath::nextpow2( inum );
      int num_freqs = 1 + klen/2;
      std::vector<double> amp( klen , 0 );
      
      //void  mtm::mt_get_spec(double *series, int inum, int klength, double *amp)
      // series = input time series
      // inum   = length of time series
      // klength = number of elements in power spectrum (a power of 2)
      // amp = returned power spectrum

      int             i, j, isign = 1;  
      unsigned long   nn;
      double          tsv;    
      nn = klen;
      
      /* copy amp onto series and apply zero padding to  klength */
      amp = x;
      amp.resize( klen );
      // for (i = 0; i < inum; i++) { amp[i] = series[i];  }
      // zero_pad(amp, inum, klength);
  
      
      /*  Fast Fourier Transform Routine:  here we are using the Numerical Recipes
	  routine jrealft which returns the fft in the 1-D input array
	  packed as pairs of real numbers.
	  The jrealft routine requires the input array to start at index=1
	  so we must decrement the index of amp
      */

      void jrealft(double data[], unsigned long n, int isign);
      double * pamp = &(amp)[0];
      jrealft(pamp-1, nn, isign);
      double anrm = sqrt(npoints/dt);  
      double sum = 0.0;

      // get spectrum from real fourier transform      
      double norm = 1.0/(anrm*anrm);
      
      for(int i=1; i<num_freqs-1; i++)
	{
	  if( 2*i+1 > klen) Helper::halt( "mtm_t error in index");  
	  double sqramp = pow( amp[2*i+1] , 2 ) + pow( amp[2*i] , 2 );
	  std::cout << 2 * norm * (sqramp) << "\n";
	  
	  // sqr_spec[i+kf] = norm * (sqramp);	   == real_FFT()
	  // sum += sqr_spec[i+kf];
	}

      std::cout << "DC " << norm * pow(fabs(amp[0]),2) << "\n"
		<< "NQ " << norm * pow(fabs(amp[1]),2) << "\n";
      
      // sqr_spec[0+kf] = norm * pow(fabs(amp[0]),2);
      // sqr_spec[num_freqs-1+kf] = norm * pow(fabs(amp[1]),2);
  
      std::exit(0);

      // test 2 : real_FFT()

      
      int index_length = x.size();	    
      int my_Fs = 256; // arbitrary

      std::cout << index_length << " is size\n";

      int index_start = 0;
      
      //FFT fftseg( index_length , index_length , my_Fs , FFT_FORWARD , WINDOW_NONE );
      real_FFT fftseg( index_length , index_length , my_Fs , WINDOW_NONE );
      
      const int reps = 5000;

      for ( int i=0; i<reps; i++)
	{
	  std::cout << "i\t" << i << "\n";
	  fftseg.apply( &(x[index_start]) , index_length );
	  
	  int my_N = fftseg.cutoff;
	  
	  // for (int f=0;f<my_N;f++)
	  //   {	  
	  //     // t : raw transform
	  //     // t2 : scaled by 1/N
	  //     // 
	  //     std::cout << f << "\t"
	  // 		<< fftseg.frq[f] << "\t"		    
	  // 		<< fftseg.X[f] << "\n";

	  //   }
	}
      std::exit(1);
    }

  if ( p == "dfa" )
    {
      std::vector<double> w(0);
      dfa_t dfa;
      int nw = 100;
      dfa.set_windows( 200 );
      
      dfa.proc( &x );
      
      for (int i=0; i<nw; i++)
	std::cout << dfa.w[i] << "\t"
		  << dfa.fluctuations[i] << "\t"
		  << dfa.slopes[i] << "\n";
      
      std::exit(1);
    }
    

  if ( p == "tri" )
    {
      int n = x.size();
      int h = 7 ;
      double w = 0.05;

      Eigen::VectorXd Y = eigen_ops::copy_array( x );

      Eigen::VectorXd Y2 = eigen_ops::tri_moving_average( Y , h , w );
      Eigen::VectorXd Y3 = eigen_ops::moving_average( Y , h );

      for (int i=0; i<n; i++)
	std::cout << Y[i] << "\t" << Y2[i] << "\t" << Y3[i] << "\n";
      
    }
  

  
  if ( p == "fft" )
    {
      
      int index_length = x.size();

      int my_Fs = 256; // arbitrary

      if ( p2 != "" )
	if ( ! Helper::str2int( p2 , &my_Fs ) ) Helper::halt( "expecting integer sample rate as second parameter" );

      int index_start = 0;

      FFT fftseg( index_length , index_length , my_Fs , FFT_FORWARD , WINDOW_NONE );
      
      fftseg.apply( &(x[index_start]) , index_length );

      // Extract the raw transform
      std::vector<std::complex<double> > t = fftseg.transform();

      // Extract the raw transform scaled by 1/n
      std::vector<std::complex<double> > t2 = fftseg.scaled_transform();
      
      int my_N = fftseg.cutoff;      
      
      //std::cout << "t sz = " << t.size() << " " << my_N << "\n";
      //for ( int i=0;i<t.size(); i++) std::cout << i << "\t" << std::real( t[i] ) << "\t" << std::imag( t[i] ) << "\n";
      //std::exit(0);
	    
      std::cout << "N" << "\t"
		<< "F" << "\t"
		<< "RE" << "\t"
		<< "IM" << "\t"
		<< "UNNORM_AMP" << "\t"
		<< "NORM_AMP" << "\t"
		<< "PSD" << "\t"
		<< "log10(PSD)" << "\n";

      for (int f=0;f<my_N;f++)
	{	  
	  // t : raw transform
	  // t2 : scaled by 1/N
	  // 
	  std::cout << f << "\t"
		    << fftseg.frq[f] << "\t"		    
		    << std::real( t[f] ) << "\t"
		    << std::imag( t[f] ) << "\t"
		    << fftseg.mag[f] << "\t"
		    << ( f == 0 ? 1 : 2 ) * fftseg.mag[f] / (double)index_length << "\t"
		    << fftseg.X[f] << "\t"
		    << log10( fftseg.X[f] ) << "\n";
	} 
      
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
      const std::string db = p2;

      //      retval_t retval = writer_t::dump_to_retval( "/Users/shaun/dropbox/projects/nsrr/tutorial/luna/out.db" , "nsrr01" );

      retval_t retval = writer_t::dump_to_retval( db );
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
      const int ns = 2;
      
      int rows = x.size() / ns;
      int cols = ns;

      Eigen::MatrixXd X( rows , cols );

      int p = 0;
      for (int i=0;i<rows;i++) 	
	for (int j=0;j<ns;j++)
	  X(i,j) = x[p++];
      
      int compc = 2;
      
      std::cerr << "performing ICA on " << rows << " x " << cols << " matrix\n";
      
      eigen_ica_t ica( X , compc );

      std::cerr << "K\n" << ica.K << "\n";
      std::cerr << "W\n" << ica.W << "\n";
      std::cerr << "A\n" << ica.A << "\n";
      std::cout <<  ica.S << "\n";
      
      std::exit(1);
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
      const int npi = 5;
      
      const int nwin = 9;

      const double segment_sec = 5;
      const double segment_step = 1;
      
      mtm_t mtm( npi , nwin );
      
      mtm.apply( &x , 256 , 256 * segment_sec , 256 * segment_step , true );
      
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

void build_param( param_t * param , int argc , char** argv , int start )
{

  //
  // get arguments from stdin (rather than the command line options?
  //
  
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



std::string luna_base_version() 
{
  std::stringstream ss;
  ss << "luna-base version " << globals::version << " (release date " << globals::date << ")\n";
  ss << "luna-base build date/time " << __DATE__ << " " << __TIME__ << "\n";
  return ss.str();
}


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


void log_commands( int argc , char ** argv )
{
  bool dump = false; 

  for (int i=0; i<argc; i++)
    {
      if ( strcmp( argv[i] ,"--log" ) == 0 )
	{
	  dump = true;
	  break;
	}
    }
  
  if ( ! dump ) return;

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
  
  std::cerr << "\n"
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
      
      std::cerr << spc << s ;
      
    }

  if ( has_s ) std::cerr << "'";

  std::cerr << "\n\n"
	    << "# " << std::string( 78 , '-' ) << "\n"
	    << "\n";

  
}

