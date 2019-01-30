
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
  
  //
  // initiate global defintions
  //

  global.init_defs();

  
  //
  // Some initial options (called prior to the main banner, etc)
  //
  
  if ( argc <= 3 && strcmp( argv[1] ,"-d" ) == 0 )
    { 
      std::string p = argc == 3 ? argv[2] : "";
      global.api();
      proc_dummy( p ); 
      exit(0); 
    } 
  else if ( argc == 2 && strcmp( argv[1] , "--eval" ) == 0 ) 
    {
      // simple test evaluation of expressions from the command line
      global.api();
      proc_eval_tester(); 
      exit(0);       
    }
  else if ( argc == 3 && strcmp( argv[1] , "--validate" ) == 0 ) 
    {
      // 1: sample-list or EDF
      // 2: --validate 
      Helper::halt( "--validate not implemented... ");
    }
  else if ( argc == 3 && strcmp( argv[1] , "--xml" ) == 0 )
    {
      global.api();
      annot_t::dumpxml( argv[2] , false );
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
      //   --flag  ( added minus fist '-' to globals::param
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
	      if ( i + 1 >= argc ) Helper::halt( "expecting database name after -o" );
	      cmd_t::stout_file = argv[ ++i ];
	      if ( Helper::iequals( tok[0] , "-a" ) ) cmd_t::append_stout_file = true;
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
	      
	      std::ifstream INC( filename , std::ios::in );
	      if ( INC.bad() ) Helper::halt("could not open file: " + filename );
	      while ( ! INC.eof() )
		{

		  std::string line;
		  std::getline( INC , line);
		  if ( INC.eof() || line == "" ) continue;
		  
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

      logger << "usage: luna [sample-list|EDF] [n1] [n2] [signal=s1,s2] [v1=val1] [@parameter-file] < command-file"
	     << std::endl;
      logger.off();
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
  writer.string_factor( globals::annot_strat );
  writer.string_factor( globals::annot_instance_strat );
  writer.string_factor( globals::annot_meta_strat );
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
  if ( failed == 0 ) logger << " all of which passed" << std::endl;
  else logger << failed << " of which failed" << std::endl;

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
      
      logger << "path    : " << globals::project_path << std::endl;
                 
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
	  logger << "\nskipping EDF " << rootname << std::endl;
	  ++processed;
	  continue; // to the next EDF in the list
	}
      

      //
      // Begin running through the series of commands
      //

      
      logger  << "\n___________________________________________________________________\n"
	      << "Processing: " << rootname 
	      << " [ #" << processed+1 << " ]" << std::endl;
      
      
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
      // Handle  luna --validate {edf|sample-list} separaetely (NOT IMPLEMENTED YET)
      //

//       if ( globals::validation_mode )
// 	{
	  
// 	  writer.cmd( "VALIDATE" , 1 , "" );

// 	  // test: can we attach this EDF?
// 	  edf_t edf;  
	  
// 	  bool okay = edf.attach( edffile , rootname );

// 	  if ( ! okay ) logger.warning( "could not attach " + edffile );
	  
// 	  writer.value( "VALID" , (int)okay );
	  
// 	  // note that we've done this
// 	  ++processed;
// 	  ++actual;

// 	  writer.commit();

// 	  if ( single_edf ) break;
	  
// 	  // go to next EDF
// 	  continue;
	  
// 	}

	  
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

	  logger << "**warning: problem loading " 
		 << edffile << ", skipping..." << std::endl;
	  
	  writer.commit();

	  continue;
	}
      
      
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
      // Attach annotations.  (run for single_edf mode, as we may have added an annotation file as above)
      //
      
//       if ( ! single_edf ) 
// 	{
	  
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
			  //std::cerr << " fname [" << fname << "]\n";
			  // only annot files (.xml, .ftr, .annot)
 			  if ( Helper::file_extension( fname , "ftr" ) ||
 			       Helper::file_extension( fname , "xml" ) ||
			       Helper::file_extension( fname , "eannot" ) ||
			       Helper::file_extension( fname , "annot" ) )
			    {
			      edf.load_annotations( tok[i] + fname );	 			   
			    }
			}
		      closedir (dir);
		    }
		  else 
		    Helper::halt( "could not open folder " + tok[i] );
		}
	      else
		{
		  edf.load_annotations( tok[i] );	 
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
		  if ( ! ( instance_ids.size() == 1 && ( *instance_ids.begin()  == names[a] || *instance_ids.begin() == "." ) ) )
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


	  //	} // end of single-EDF check 

	    
      //
      // Evaluate all commands
      //

      bool cmd_okay = cmd.eval( edf );
      
      // TODO: we need to track what passed/failed here?? 


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

  logger << std::endl
	 << "___________________________________________________________________"
	 << std::endl
	 << "...processed " << actual << " EDFs, done."
	 << std::endl;

}



// EVAL expresions

void proc_eval_tester()
{

  // read a single line
  std::string expr;
  std::getline( std::cin , expr );

  std::map<std::string,annot_map_t> inputs;

  instance_t out;
  
  Eval tok( expr );

  tok.bind( inputs , &out );

  bool is_valid = tok.evaluate();
  
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
 
  //
  // Straight FFT of stdin
  //

  std::vector<double> x;
  
  if ( p == "fft" || p == "mtm" ) 
    {

      
      while ( ! std::cin.eof() )
	{
	  double xx;
	  std::cin >> xx;
	  if ( std::cin.eof() ) break;	  
	  x.push_back( xx );	  
	}
      std::cerr << x.size() << " values read\n";

    }

  if ( p == "fft" )
    {
      int index_length = x.size();
      int my_Fs = 20; // arbitrary
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
  
  
  std::string expr;
  
  std::getline( std::cin , expr );
  
  std::map<std::string,annot_map_t> inputs;
  
  annot_map_t amap1, amap2;
  
  std::vector<int> ints(3);
  ints[0] = 637; ints[1] = 87; ints[3] = -23;

  instance_t in1;
  in1.set( "v1" , 10.0 );
  in1.set( "v2" , std::string("A") );
  in1.set( "v3" , true );

  instance_t in2;
  in2.set( "v1" , 92.1 );
  in2.set( "v2" , std::string( "B" ) );
  in2.set( "v3" , true );

  instance_t in3;
  in3.set( "v1" , 108.5 );
  in3.set( "v2" , std::string( "C" ) );
  in3.set( "v3" , false );

  annot_t annot1( "a1" );
  
  amap1[ instance_idx_t( &annot1 , interval_t(0,100) , "inst1" ) ] = &in1;
  amap1[ instance_idx_t( &annot1 , interval_t(200,300) , "inst2" ) ] = &in2;
  amap1[ instance_idx_t( &annot1 , interval_t(200,300) , "inst3" ) ] = &in3;

//   amap2[ instance_idx_t( &annot2 , interval_t(0,100) , "." ) ] = &in3;
//   amap2[ instance_idx_t( &annot2 , interval_t(200,300) , "." ) ] = &in4;

  // collate all annotations
  inputs[ "a1" ] = amap1;
  //  inputs[ "annot2" ] = amap2;
  
  instance_t out;
  
  Eval tok( expr );

  tok.bind( inputs , &out );


  bool is_valid = tok.evaluate();
  
  bool retval; 
  
  bool is_valid2 = tok.value( retval );
  
  std::cout << "\n\n\n";
  std::cout << "valid " << is_valid << "\t" << is_valid2 << "\n";
  std::cout << "value  = " << retval << "\n";
  std::cout << "result = " << tok.result() << "\n";

  std::cout << "\n\n\n";
  
  std::cout << "out : " << out.print() << "\n";
  
  std::exit(1);
  
  
  
  //std::cout << "annot " << nsrr_t::remap( "Stage 4 sleep|4" ) << "\n";

  std::exit(0);

  
  std::vector<double> dx;
  while ( ! std::cin.eof() ) 
    {
      double dxx;
      std::cin >> dxx;
      if ( std::cin.eof() ) break;
      dx.push_back( dxx );
    }
  std::cout << "read " << dx.size() << " points\n";

  // nw = 4, i.e. 2*4-1 = 7 tapers
  // 1 window ??

//   mtm_t::mtm_t( const int npi , const int nwin ) : npi(npi) , nwin(nwin)
//   {

//     kind = 1;
//     inorm = 1;
//   }

  // P number of points 
  // W bandwidth 
  // W = P/N
  // P = NW  , ... 4 
 
  // N = 512, i.e. 
  // or 320?  i.e. length of sequence
  // P = NW
  // 

  //         npi , nwin 
  //          P     # tapers (should be max 2P-1)
  //          P=NW
  mtm_t mtm(   4 , 7 );

  //  mtm.inorm = 0; // normalization factor
  //   mtm.kind = 1; // 1 or 2 (weights)

  // adaptive weights and 1/N normalization 
  // seems to give similar results to ML pmtm
  
  // weights
  // 1 hi-res, 2 adaptive, 3 none

  mtm.kind = 2 ; 

  //normalization 1/N
  mtm.inorm = 3 ; 

  // inorm               norm
  // 0  1 (none)      -> 1
  // 1  N             -> 1/N^2
  // 2  1/dt          -> dt^2
  // 3  sqrt(N)       -> 1/N
  
  // norm = 1/(W^2)

  mtm.apply( &dx , 1024 );

  

  //
  // standard
  //

//   int index_length = dx.size();
//   int my_Fs = 1024;
//   int index_start = 0;
  

//   FFT fftseg( index_length , my_Fs , FFT_FORWARD , WINDOW_NONE );

//   fftseg.apply( &(dx[index_start]) , index_length );

//   int my_N = fftseg.cutoff;
//   for (int f=0;f<my_N;f++)
//     {
//       std::cout << f << "\t" << fftseg.frq[f] << "\t" << fftseg.X[f] << "\n";
//    }
  
  
    
  std::exit(1);


  std::string edf_file2 = "/home/shaun/Dropbox/my-sleep/Purcell06072016.edf";  
  edf_t edf2; 
  edf2.attach( edf_file2 , "id-007" );
  
    
  retval_t ret;
  writer.nodb();
  writer.use_retval( &ret );

  // evaluate a command
  cmd_t cmd( "EPOCH & SPINDLES method=wavelet fc=13" );

  cmd.eval( edf2 ); 

  ret.dump();

  std::exit(0);
  
  // retval_strata_t s;
  // s.add( retval_factor_level_t( "F" , 11 ) ) ;
  // s.add( retval_factor_level_t( "SS" , "N2" ) ) ;
  // s.add( retval_factor_level_t( "B" , 12.4 ) ) ; 

  // ret.add( retval_cmd_t("c1") , retval_var_t("v1") , s , 22 );
  // ret.add( retval_cmd_t("c2") , retval_var_t("v1") , s , 22.2 );
  // ret.add( retval_cmd_t("c2") , retval_var_t("v2") , s , 22.5 );
  // ret.dump();
  
  std::exit(1);
  
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

  logger << "read " << x.size() << " x values" << std::endl;
  //  logger << "read " << h.size() << " h values\n";
  
  
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
      logger << "tt + " << tt << "\n";

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
//   logger << "read " << x.size() << " values\n";

   while ( ! std::cin.eof() ) 
     {
       double t;
       std::cin >> t;
       if ( std::cin.eof() ) break;
       x.push_back(t);
     }
   logger << "read " << x.size() << " values\n";


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

//   logger << "et = " << empirical_threshold << "\n";

//   std::map<double,double>::const_iterator tt = tvals.begin();
//   while ( tt != tvals.end() ) 
//     {
//       logger << tt->first << "\t" 
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
      logger << "read " << x.size() << " values\n";
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

