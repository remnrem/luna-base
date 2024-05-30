
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


// 1) upload signal
// 2) db -> rrtables 
// 3) error/handling/checkpoints

#include "lunapi/lunapi.h"
#include <stdexcept>
#include <memory>

extern logger_t logger;

extern globals global;

extern cmd_t cmd;

// --------------------------------------------------------------------------------
// global singleton class functions

void lunapi_bail_function( const std::string & msg )
{
  throw( std::runtime_error( msg ) );
}
  
void lunapi_msg_function( const std::string & msg )
{
  std::cerr << " [lunapi] :: " << msg << "\n";
  return;
}

void lunapi_t::init()
{

  global.init_defs();

  globals::bail_function = &lunapi_bail_function;
  
  globals::bail_on_fail = false;
  
  global.R( 1 ); // 1 means to cache the log
  
  writer.nodb();
  
  logger << "** luna " 
   	 << globals::version 
   	 << " " << globals::date
   	 << "\n";

  logger.print_buffer();

}

void lunapi_t::silence( const bool b )
{
  globals::silent = b;
}

bool lunapi_t::is_silenced() const
{
  return globals::silent;
}

void lunapi_t::var( const std::string & key , const std::string & value )
{
  
  // special treatment for `sig`.   This will just append
  // to an existing signallist;  as we don't always want to have
  // to lreset(), make this one case so that sig clears signlist prior 
  // to setting if the signal list is "." 
  
  if ( key == "sig" && value == "." )
    cmd_t::signallist.clear();
  else
    {
      //std::cout << " parsing [" << key << "] -> [" << value << "]\n";
      cmd_t::parse_special( key , value );
    }
  //logger << "setting " << key << " = " << value << "\n";
}


void lunapi_t::dropallvars()
{
  cmd_t::vars.clear();
}

void lunapi_t::dropvars( const std::vector<std::string> & keys )
{
  for (int i=0; i<keys.size(); i++) dropvar( keys[i] );
}

void lunapi_t::dropvar( const std::string & key )
{
  std::map<std::string,std::string>::iterator ii = cmd_t::vars.find( key );
  if ( ii != cmd_t::vars.end() )
    cmd_t::vars.erase( ii );
  return;  
}

std::map<std::string,std::string> lunapi_t::vars() const
{
  return cmd_t::vars;
}

std::map<std::string,std::variant<std::monostate,std::string> > lunapi_t::vars( const std::vector<std::string> & keys ) const
{
  std::map<std::string,std::variant<std::monostate,std::string> > r;
  for (int i=0; i<keys.size(); i++) 
    r[ keys[i] ] = var( keys[i] );
  return r;
}

std::variant<std::monostate,std::string> lunapi_t::var( const std::string & key ) const
{  
  if ( cmd_t::vars.find( key ) == cmd_t::vars.end() ) return std::monostate{};
  return cmd_t::vars[ key ];
}


void lunapi_inst_t::ivar( const std::string & key , const std::string & value )
{
  cmd_t::ivars[ id ][ key ] = value;
}

std::variant<std::monostate,std::string> lunapi_inst_t::ivar( const std::string & key ) const
{
  if ( cmd_t::ivars[ id ].find( key ) == cmd_t::ivars[ id ].end() ) return std::monostate{};
  return cmd_t::ivars[ id ][ key ];
}

std::map<std::string,std::variant<std::monostate,std::string> > lunapi_inst_t::ivars() const
{
  std::map<std::string,std::variant<std::monostate,std::string> > r;
  std::map<std::string,std::map<std::string,std::string> >::const_iterator ii = cmd_t::ivars.find( id );
  if ( ii == cmd_t::ivars.end() ) return r;
  const std::map<std::string,std::string> & v = ii->second;
  std::map<std::string,std::string>::const_iterator vv = v.begin();
  while ( vv != v.end() )
    {
      r[ vv->first ] = vv->second;
      ++vv;
    }
  return r;
}


void lunapi_inst_t::clear_selected_ivar( const std::set<std::string> & keys )
{

  // only for this indiv
  std::map<std::string,std::map<std::string,std::string> >::const_iterator ii = cmd_t::ivars.find( id );
  if ( ii == cmd_t::ivars.end() ) return;

  // copy over non-cleared parts
  std::map<std::string,std::string> cp ;
  
  // iterate over values
  const std::map<std::string,std::string> & v = ii->second;
  std::map<std::string,std::string>::const_iterator vv = v.begin();
  while ( vv != v.end() )
    {
      if ( keys.find( vv->first ) == keys.end() )
	cp[ vv->first ] = vv->second;
      ++vv;
    }

  // update
  cmd_t::ivars[ id ] = cp;
  
}

void lunapi_inst_t::clear_ivar()
{
  if ( cmd_t::ivars.find( id ) == cmd_t::ivars.end() ) return;
  cmd_t::ivars[ id ].clear();
}

// global clear all ivars
void lunapi_t::clear_ivars()
{
  cmd_t::ivars.clear();
}

void lunapi_t::reset() const
{
  // status flags: note - a fundamental problem that we now allow multiple EDFs
  // to be attached but are still working with global status flags
  
  // -- for typical workflows this should not be a problem... but we should
  //    fix at some point (i.e. make empty a property of edf_t at least
  
  globals::problem = false;  
  globals::empty = false;
}

void lunapi_t::flush()
{
  logger.flush_cache();
}

void lunapi_t::re_init()
{
  // clear all variables (both user-defined and special)
  // but do not alter the EDF attachment

  // clear all user-defined variables, signal lists and aliases  
  cmd_t::clear_static_members();
  
  // also reset global variables that may have been changed since
  global.init_defs();
  
  // but need to re-indicate that we are running inside API
  global.R( 1 ); // 1 means to cache

  reset();
}

// read a Luna @include file and set variables
int lunapi_t::includefile( const std::string & f )
{
  const std::string filename = Helper::expand( f );
  
  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "cannot open " + filename );
  
  // nb. - should make this a single function to share w/ main()

  int tokens = 0;

  bool parse_line = true;
  std::string last_grp = "";

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
	{
	  Helper::halt("badly formatted line ( # tabs != 2 ) in " + filename + "\n" + line );
	  return tokens;
	}
      
      ++tokens;

      logger << "  setting " << tok[0] << " = " << tok[1] << "\n";

      cmd_t::parse_special( tok[0] , tok[1] );
      
    }
  
  INC.close();

  return tokens;
  
}


// read a Luna command file from text and parse as a string
std::string lunapi_t::cmdfile( const std::string & f )
{

  // 1) remove % comments
  // 2) append lines starting with space 

  const std::string filename = Helper::expand( f );

  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "cannot open " + filename );

  bool first = true; 

  std::string cmdstr; 

  std::ifstream IN1( filename.c_str() , std::ios::in );
  while ( ! IN1.eof() )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );

      if ( IN1.eof() || line == "" ) continue;
      
      // skip % comments                                                                                                      
      if ( line[0] == '%' ) continue;
      if ( line.find( "%" ) != std::string::npos )
	line = line.substr( 0 , line.find( "%" ) );
      if ( line.size() == 0 ) continue;
      
      // append same command (if space ident, else new command)
      if ( line[0] != ' ' )
	{
	  if ( ! first ) cmdstr += " & ";
	  else first = false;
	}
  
      // append
      cmdstr += line;
      
    }

  IN1.close();

  return cmdstr;
}


// aliases/remap table

std::vector<std::vector<std::string> > lunapi_t::aliases() const
{
  std::vector<std::vector<std::string> > t;

  // channels
  // mapped --> alias
  std::map<std::string,std::string>::const_iterator ss = cmd_t::label_aliases.begin();
  while ( ss != cmd_t::label_aliases.end() )
    {
      std::vector<std::string> line = { "CH" , ss->second , ss->first };
      t.push_back(line);
      ++ss;
    }

  // annotations
  // alias --> orig
  std::map<std::string,std::string>::const_iterator aa = nsrr_t::amap.begin();
  while ( aa != nsrr_t::amap.end() )
    {
      std::vector<std::string> line = { "ANNOT" , aa->second , aa->first };
      t.push_back(line);
      ++aa;
    }

  return t;
}


// import a generic Luna db

std::vector<std::string> lunapi_t::import_db( const std::string & dbfile )
{
  std::set<std::string> ids0;
  return import_db( dbfile , ids0 );
}

// import subset of indivs from a generic luna db
std::vector<std::string> lunapi_t::import_db( const std::string & dbfile , const std::set<std::string> & ids )
{

  // this gets populated by the IDs actually read
  std::vector<std::string> obs_ids;

  if ( ! Helper::fileExists( dbfile ) ) return obs_ids;
    
  retval_t ret = writer_t::dump_to_retval( dbfile , &ids , &obs_ids );
  
  logger << "  read data on " << obs_ids.size() << " individuals from " << dbfile << "\n" ;

  // store in the internal rtables cache
  rtables = rtables_t( ret );
  
  return obs_ids;

}


// --------------------------------------------------------------------------------
// sample list functions

// --------------------------------------------------------------------------------
// sample list functionality

int lunapi_t::build_sample_list( const std::vector<std::string> & toks )
{
  // clear any existing sample list
  clear();

  // build up the SL, saving to 'sl'
  slist_t sl;
  Helper::build_sample_list( toks , &sl );
  
  // populate this class
  for (int i=0;i<sl.size();i++)
    {
      if ( std::get<0>(sl[i]) != "" && std::get<1>(sl[i]) != "" )
	insert_inst( std::get<0>(sl[i]) , std::get<1>(sl[i]) , std::get<2>(sl[i]) );
    }
  return nobs();
}


int lunapi_t::read_sample_list( const std::string & file )
{
  const std::string filename = Helper::expand( file );
  if  ( ! Helper::fileExists( filename ) ) Helper::halt( "could not open sample list " + filename );

  const bool has_project_path = globals::param.has( "path" );
  
  if ( has_project_path ) 
    globals::project_path = globals::param.value( "path" ) ;
  
  std::ifstream IN1( filename.c_str() , std::ios::in );
  while ( ! IN1.eof() )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      
      if ( IN1.eof() || line == "" ) continue;

      std::vector<std::string> tok = Helper::parse( line , "\t" );
      if ( tok.size() == 0 ) continue;
      if ( tok.size() < 2 || tok.size() > 3 )
	continue;

      // splice in project path?
      if ( has_project_path )
	{
	  // EDF
	  if ( tok[1][0] != globals::folder_delimiter )
	    tok[1] = globals::project_path + tok[1];	  
	}

      // annotations

      std::set<std::string> aset;
      if ( tok.size() == 3 )
	{
	  std::vector<std::string> toka = Helper::parse( tok[2] , "," );
	  for (int a=0;a<toka.size();a++)
	    {
	      if ( has_project_path )
		if ( toka[a][0] != globals::folder_delimiter )
		  toka[a] = globals::project_path + toka[a];	      
	      aset.insert( toka[a] );
	    }
	} 


      // insert
      insert_inst( tok[0] , tok[1] , aset );
            
    }
  
  return nobs();
}

std::vector<std::tuple<std::string,std::string,std::set<std::string> > > lunapi_t::sample_list() const
{
  std::vector<std::tuple<std::string,std::string,std::set<std::string> > > r;
  std::map<std::string,std::string>::const_iterator ee = edfs.begin();
  while ( ee != edfs.end() )
    {
      std::optional<int> n = get_n( ee->first );
      if ( n ) 
	r.push_back( std::make_tuple( ee->first , ee->second , get_annot( *n ) ) );
      ++ee;
    }
  return r;
}

void lunapi_t::insert_inst( const std::string & id ,
			    const std::string & edf ,
			    const std::set<std::string> & annot )
{
  int curr = edfs.size();
  edfs[ id ] = edf;
  annots[ id ] = annot;
  id2n[ id ] = curr;
  n2id[ curr ] = id;
}

int lunapi_t::nobs() const
{
  return edfs.size();
}

void lunapi_t::clear()
{
  edfs.clear();
  annots.clear();
  id2n.clear();
  n2id.clear();
}
  
std::optional<std::string> lunapi_t::get_id( const int i ) const
{
  std::map<int,std::string>::const_iterator ii = n2id.find( i ) ;
  if ( ii == n2id.end() ) return std::nullopt;
  return ii->second;
}

std::string lunapi_t::get_edf( const int i ) const
{
  std::map<int,std::string>::const_iterator ii = n2id.find( i ) ;
  if ( ii == n2id.end() ) return "";
  return edfs.find( ii->second )->second;
}

std::set<std::string> lunapi_t::get_annot( const int i ) const
{
  std::map<int,std::string>::const_iterator ii = n2id.find( i ) ;
  if ( ii == n2id.end() ) return std::set<std::string>();
  return annots.find( ii->second )->second;
}

std::optional<int> lunapi_t::get_n( const std::string & id ) const
{
  std::map<std::string,int>::const_iterator ii = id2n.find( id ) ;
  if ( ii == id2n.end() ) return std::nullopt;
  return ii->second;
}


lunapi_inst_ptr lunapi_t::inst( const std::string & id ) const
{
  reset();
  lunapi_inst_ptr p( new lunapi_inst_t( id ) );
  return p;
}
  
lunapi_inst_ptr lunapi_t::inst( const std::string & id , const std::string & edf ) const
{
  reset();
  lunapi_inst_ptr p( new lunapi_inst_t( id ) );
  p->attach_edf( edf );
  return p;
}

lunapi_inst_ptr lunapi_t::inst( const std::string & id , const std::string & edf , const std::string & annot ) const
{
  reset();
  lunapi_inst_ptr p( new lunapi_inst_t( id ) );
  p->attach_edf( edf );
  p->attach_annot( annot );
  return p;
}

lunapi_inst_ptr lunapi_t::inst( const std::string & id , const std::string & edf , const std::set<std::string> & annots ) const
{
  reset();
  lunapi_inst_ptr p( new lunapi_inst_t( id ) );
  p->attach_edf( edf );
  std::set<std::string>::const_iterator aa = annots.begin();
  while ( aa != annots.end() )
    {
      p->attach_annot( *aa );
      ++aa;
    }
  return p;
}

std::optional<lunapi_inst_ptr> lunapi_t::inst( const int i ) const
{
  reset();
  std::optional<std::string> id = get_id( i );
  if ( ! id ) return std::nullopt;
  
  lunapi_inst_ptr p( new lunapi_inst_t( *id ) );

  // edf
  p->attach_edf( get_edf( i ) );

  // annots
  std::set<std::string> a = get_annot( i );
  std::set<std::string>::const_iterator aa = a.begin();
  while ( aa != a.end() )
    {
      p->attach_annot( *aa );
      ++aa;
    }
  
  return p;
}


//
// --------------------------------------------------------------------------------
// project level desc() convenience function


std::vector<std::vector<std::string> > lunapi_t::desc()
{
  std::vector<std::vector<std::string> > r;
  for (int i=0; i<nobs(); i++)
    {
      std::optional<lunapi_inst_ptr> l1 = inst( i );
      if ( l1 )
        {
          lunapi_inst_ptr p1 = *l1;
	  r.push_back( p1->desc() );
        }
    }  
  return r;
}


//
// --------------------------------------------------------------------------------
// evaluate Luna commands across multiple individuals
//   

rtables_return_t lunapi_t::eval( const std::string & cmdstr )
{
  
  retval_t accumulator;
  
  writer.clear();
  writer.set_types(); // likely not needed, but harmless to keep
  writer.use_retval( &accumulator );
    
  for (int i=0; i<nobs(); i++)
    {
      std::optional<lunapi_inst_ptr> l1 = inst( i );
      if ( l1 ) 
	{
	  // clear any problem flags
	  reset();
	  lunapi_inst_ptr p1 = *l1;
	  p1->eval_project( cmdstr , &accumulator );
	}
    }

  //  accumulator.dump();
    
  // get all results
  rtables = rtables_t( accumulator );

  writer.use_retval( NULL );
  writer.clear();
  writer.set_types();
      
  //  rtables.dump();
  
  return rtables.data();
  
}


rtable_t lunapi_t::table( const std::string & cmd , const std::string & faclvl ) const
{
  return rtables.table( cmd , faclvl );
}

std::vector<std::string> lunapi_t::variables( const std::string & cmd , const std::string & faclvl ) const
{
  return rtables.table( cmd , faclvl ).cols;
}

rtable_return_t lunapi_t::results( const std::string & cmd , const std::string & faclvl ) const
{
  return rtables.data( cmd , faclvl );
}

rtables_return_t lunapi_t::results() const
{
  return rtables.data() ;
}



//
//
// --------------------------------------------------------------------------------
// instance functions
//
//


std::string lunapi_inst_t::get_id() const
{
  return id;
}

int lunapi_inst_t::get_state() const
{
  // 0 empty; +1 attached okay, -1 problem
  return state;
}

double lunapi_inst_t::last_sec() const
{
  // last addressable timepoint (in current internal EDF)
  return ( edf.timeline.last_time_point_tp + 1LLU ) * globals::tp_duration;
}

double lunapi_inst_t::last_sec_original() const
{
  // last addressable timepoint (from original EDF)
  return ( edf.header.last_time_point_tp_orig + 1LLU ) * globals::tp_duration;
}

std::string lunapi_inst_t::get_edf_file() const
{
  return edf_filename;
}

std::string lunapi_inst_t::get_annot_files() const
{
  return Helper::stringize( annot_filenames );
}

void lunapi_inst_t::refresh()
{
  if ( state != 1 ) 
    {
      Helper::halt( "lunapi_inst_t::refresh(): no attached EDF" );
      return;      
    }
  
  // drop edf_t  
  edf.init();
  
  // reattach EDF (and this will remake the timeline too)
  attach_edf( edf_filename );
  
  if ( state != 1 ) 
    {
      Helper::halt( "lunapi_t::refresh(): problem reattaching EDF" );
      return;
    }
  
  // reload annotations
  std::set<std::string>::const_iterator aa = annot_filenames.begin();
  while ( aa != annot_filenames.end() )
    {
      edf.load_annotations( *aa );
      ++aa;
    }
  
}

void lunapi_inst_t::drop()
{
  // clear out
  edf.init();
  // set to empty
  edf_t empty;
  edf = empty;

  // track meta-data
  state = 0;
  id = "";
  edf_filename = "";
  annot_filenames.clear();
}

std::vector<std::string> lunapi_inst_t::desc()
{
  std::vector<std::string> ret;
  param_t p0;
  p0.add( "sig" , "*" );
  edf.description( p0, &ret ); 
  return ret;  
}


std::map<std::string,datum_t> lunapi_inst_t::status() const 
{
  
  std::map<std::string,datum_t> r;

  r[ "state" ]  = state;
  
  if ( state != 1 ) return r;
  
  r[ "edf_file" ] = edf_filename;
  
  r[ "annotation_files" ] = Helper::stringize( annot_filenames );

  int n_data_channels = 0;
  int n_annot_channels = 0;
  for (int i=0;i<edf.header.ns;i++)
    {
      if ( edf.header.is_data_channel( i ) ) ++n_data_channels;
      else ++n_annot_channels;
    }

  r[ "id" ] = edf.id;
  r[ "ns" ] = n_data_channels;
  r[ "nt" ] = edf.header.ns_all;
  r[ "na" ] = (int)edf.timeline.annotations.names().size();
  
  // Record duration, as hh:mm:ss string

  uint64_t duration_tp = globals::tp_1sec
    * (uint64_t)edf.header.nr
    * edf.header.record_duration;  
  std::string total_duration_hms = Helper::timestring( duration_tp );
  
  r [ "duration" ] = total_duration_hms;

  // epoch/mask info
  if ( edf.timeline.epoched() )
    {
      r[ "ne" ] = edf.timeline.num_epochs();
      r[ "elen" ] = edf.timeline.epoch_length();
      r[ "nem" ] = edf.timeline.num_total_epochs() - edf.timeline.num_epochs();
    }

  return r;
  
}



bool lunapi_inst_t::attach_edf( const std::string & _filename )
{
  
  const std::string filename = Helper::expand( _filename );
  
  if ( ! Helper::fileExists( filename ) ) 
    Helper::halt( "cannot find " + filename );

  // restrict to limited set of input signals?
  const std::set<std::string> * inp_signals = cmd.signals().size() > 0 ? &cmd.signals() : NULL;
  
  // load EDF
  bool okay = edf.attach( filename , id , inp_signals );

  if ( ! okay )
    {
      state = -1;
      return false;
    }

  edf_filename = filename;

  // EDF+ annotations?
  if ( edf.header.edfplus )
    {
      // must read if EDF+D (but only the time-track will be taken in)                                          
      // if EDF+C, then look at 'skip-edf-annots' flag                                                          
      if ( edf.header.continuous && ! globals::skip_edf_annots )
        edf.timeline.annotations.from_EDF( edf , edf.edfz_ptr() );
      else if ( ! edf.header.continuous )
        edf.timeline.annotations.from_EDF( edf , edf.edfz_ptr() );
    }

  cmd_t::define_channel_type_variables( edf );

  state = 1; 
    
  return true;
  
}


bool lunapi_inst_t::attach_annot( const std::string & annotfile )
{
  
  if ( annotfile.size() == 0 ) return false;

  // is 'annotfile' in fact a folder (i.e. ending in '/') ?
  
  if ( annotfile[ annotfile.size() - 1 ] == globals::folder_delimiter )
    {
      
      // this means we are specifying a folder, in which case search for all files that
      // start id_<ID>_* and attach those

      DIR * dir;
      struct dirent *ent;
      if ( (dir = opendir ( annotfile.c_str() ) ) != NULL )
        {
          /* print all the files and directories within directory */
          while ((ent = readdir (dir)) != NULL)
            {
              std::string fname = ent->d_name;
	      
              if ( Helper::file_extension( fname , "ftr" ) ||
                   Helper::file_extension( fname , "xml" ) ||
                   Helper::file_extension( fname , "eannot" ) ||
                   Helper::file_extension( fname , "annot" ) )
                {
                  edf.load_annotations( annotfile + fname );
		  annot_filenames.insert( annotfile + fname );
                }
            }
          closedir (dir);
        }
      else
        {
          Helper::halt( "could not open folder " + annotfile );
	  return false;
        }
    }
  
  //
  // else a single file, load it
  //
  else
    {
      edf.load_annotations( Helper::expand( annotfile ) );
      annot_filenames.insert( annotfile );        
    }
  
  return true;

}



// #1 eval returning all output to caller 
std::tuple<std::string,rtables_return_t> lunapi_inst_t::eval_return_data( const std::string & cmdstr )
{  
  const std::string s = eval( cmdstr );
  return std::make_tuple( s , rtables.data() );
}

std::string lunapi_inst_t::eval( const std::string & cmdstr )
{
  return eval1( cmdstr , NULL );
}

std::string lunapi_inst_t::eval_project( const std::string & cmdstr , retval_t * accumulator )
{
  return eval1( cmdstr , accumulator );
}

// #2: eval, but not returning outputs to caller (stores in lunapi_t::rtables)
std::string lunapi_inst_t::eval1( const std::string & cmdstr , retval_t * accumulator )
{

  //
  // set up retval_t mechanism to catch outputs
  //

  retval_t ret;

  if ( ! accumulator )
    {
      writer.clear();
      writer.set_types(); // not sure this is needed now...
      writer.use_retval( &ret );
    }
  
  //
  // set ID 
  //

  writer.id( id , edf_filename );

  //
  // set command string
  //
  
  cmd_t cmd( cmdstr );
  
  //
  // replace any variables (or @includes, conditionals,etc) into command
  //

  cmd.replace_wildcards( id );

  
  //
  // eval on the current EDF
  //

  cmd.eval( edf );

  
  //
  // switch off the retval stream (which is local to this function and
  // so will be deleted when leaving this scope) and clear the writer
  // (ensures prior strata not applied to next run)
  //

  if ( ! accumulator )
    {
      writer.use_retval( NULL );
      writer.clear();
      writer.set_types();
    }

  //
  // get any results
  //

  if ( ! accumulator ) 
    rtables = rtables_t( ret );

  //
  // was a problem flag set?
  //
  
  if ( globals::problem ) 
    Helper::halt( "problem flag set: likely no unmasked records left?" );
  
  //
  // all done
  //

  if ( accumulator ) return "";
  
  return logger.print_buffer();
  
}


rtable_t lunapi_inst_t::table( const std::string & cmd , const std::string & faclvl ) const
{
  return rtables.table( cmd , faclvl );
}

std::vector<std::string> lunapi_inst_t::variables( const std::string & cmd , const std::string & faclvl ) const
{
  return rtables.table( cmd , faclvl ).cols;
}

rtable_return_t lunapi_inst_t::results( const std::string & cmd , const std::string & faclvl ) const
{
  return rtables.data( cmd , faclvl );
}

rtables_return_t lunapi_inst_t::results() const
{
  return rtables.data() ;
}

//
// fetch signal data :
//   - given either list of epochs or intervals
//   - either combining all data into a single frame, or keeping separate


lint_t lunapi_inst_t::epochs2intervals( const std::vector<int> & epochs )
{

  lint_t r;

  if ( state != 1 )
    return r;
  
  edf.timeline.ensure_epoched();  
  
  int total_epochs = edf.timeline.num_total_epochs();
  
  for (int epoch=0;epoch<epochs.size();epoch++)
    {
      // passed 1-based e-counts
      if ( epochs[epoch] < 1 || epochs[epoch] > total_epochs ) continue;

      // internally, 0-based
      int epoch0 = epochs[epoch] - 1;
      
      interval_t interval = edf.timeline.epoch( epoch0 );

      r.push_back( std::make_tuple( interval.start , interval.stop ) );
    }
  
  return r;
}
  
lint_t lunapi_inst_t::seconds2intervals( const std::vector<std::tuple<double,double> > & s )
{
  lint_t r;
  
  for (int i=0;i<s.size();i++)
    r.push_back( std::make_tuple( std::get<0>(s[i]) * globals::tp_1sec ,
				  std::get<1>(s[i]) * globals::tp_1sec ) );
  
  return r;
}
  


bool lunapi_inst_t::proc_channots( const std::string & chstr ,
				   const std::string & anstr ,
				   std::vector<std::string> * columns,
				   signal_list_t * signals , 
				   std::map<std::string,int> * atype )
{
  
  if ( state != 1 )
    return false;

  // Annotations: 0 not found, 1 interval, 2 epoch
  //  do not support epoch-annots right now
  std::vector<std::string> ans = Helper::parse( anstr , "," );
     
  for (int i=0;i<ans.size();i++)
    {
      if ( edf.timeline.annotations( ans[i] ) != NULL ) // is this an interval annotation? 
	(*atype)[ ans[i] ] = 1;
      else
	(*atype)[ ans[i] ] = 0;	  
    }
  
  // alphabetical order of annots:
  std::map<std::string,int>::const_iterator aa = atype->begin();
  while ( aa != atype->end() ) { columns->push_back( aa->first ) ; ++aa; } 
  
  // get signals  
  *signals = edf.header.signal_list( chstr );
  
  // check similar SRs  
  int fs = -1; 
  for (int s=0; s< signals->size(); s++) 
    {      
      if ( edf.header.is_data_channel( (*signals)(s) ) )
	{
	  columns->push_back( signals->label(s) );
	  if ( fs < 0 ) fs = edf.header.sampling_freq( (*signals)(s) );
	  else if ( edf.header.sampling_freq( (*signals)(s) ) != fs ) 
	    Helper::halt( "requires uniform sampling rate across signals" );	
	}
    }
  return true;
}


ldat_t lunapi_inst_t::data( const std::vector<std::string> & chs ,			
			    const std::vector<std::string> & anns ,
			    const bool time_track )
{
  const interval_t whole = edf.timeline.wholetrace();
  lint_t w;
  w.push_back( std::make_tuple( whole.start , whole.stop ) );  
  return slice( w , chs , anns , time_track ); 
}


ldat_t lunapi_inst_t::slice( const lint_t & intervals , 
			     const std::vector<std::string> & chs ,			
			     const std::vector<std::string> & anns ,
			     const bool time_track )
{
  
  if ( state != 1 ) 
    return std::make_tuple( std::vector<std::string>(0), Eigen::MatrixXd::Zero(0,0) );

  const std::string chstr = Helper::stringize( chs );
  const std::string anstr = Helper::stringize( anns );
  
  // labels
  std::vector<std::string> columns;
  if ( time_track ) columns.push_back( "T" );  
  std::map<std::string,int> atype;     
  signal_list_t signals;
  
  // proc channels/annots
  if ( ! proc_channots( chstr , anstr , &columns , &signals, &atype ) )
    return std::make_tuple( std::vector<std::string>(0), Eigen::MatrixXd::Zero(0,0) );
  
  // pull data
  return std::make_tuple( columns, matrix_internal( intervals, signals, atype , time_track ) );
  
}


ldats_t lunapi_inst_t::slices( const lint_t & intervals , 
			       const std::vector<std::string> & chs ,
			       const std::vector<std::string> & anns ,
			       const bool time_track )
{
  
  
  if ( state != 1 ) 
    return std::make_tuple( std::vector<std::string>(0), std::vector<Eigen::MatrixXd>(0) );

  const std::string chstr = Helper::stringize( chs );
  const std::string anstr = Helper::stringize( anns );

  std::vector<std::string> columns;
  if ( time_track ) columns.push_back( "T" );
  std::map<std::string,int> atype;
  signal_list_t signals;
  
  // get/check channel labels etc
  if ( ! proc_channots( chstr , anstr , &columns , &signals, &atype ) )
    return std::make_tuple( std::vector<std::string>(0), std::vector<Eigen::MatrixXd>(0) );
  
  // iterate over each interval
  std::vector<Eigen::MatrixXd> data;
  for ( int i=0; i<intervals.size(); i++)
    {
      lint_t i1( 1, intervals[i] );
      data.push_back( matrix_internal( i1 , signals, atype , time_track ) );
    }

  // return all 
  return std::make_tuple( columns , data );
  
}



Eigen::MatrixXd lunapi_inst_t::matrix_internal( const lint_t & intervals , 
						const signal_list_t & signals , 
						const std::map<std::string,int> & atype ,
						const bool time_track )
  
  
{
  
  const int ni = intervals.size();

  const int na = atype.size();

  // count signals
  int ns = 0;
  for (int s = 0 ; s < signals.size() ; s++ )
    if ( edf.header.is_data_channel( signals(s) ) ) ++ns;  
  if ( ns == 0 ) 
    Helper::halt( "requires at least one channel/data signal" );
    
  // # of columns: (T) + NS + NA   
  const int ncols = (int)time_track + ns + na;
  
  // # of rows... pull records 
  int nrows = 0;    
  for (int i=0;i<ni;i++)
    {      
      // Interval 
      const interval_t interval( std::get<0>(intervals[i]) , std::get<1>(intervals[i]) );
      
      // arbitary: first signal
      slice_t slice( edf , signals(0) , interval );      
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      nrows += tp->size();      
    }

  // allocate matrix
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero( nrows , ncols );  

  // first signal starts (after T and) annotations
  int s_col = (int)time_track + na;
  bool first = true;
  
  // Iterate over signals
  for (int s=0; s<ns; s++) 
    {
      
      uint64_t row = 0;
      
      // Consider each interval
      
      for (int i=0;i<ni;i++)
	{
	  
	  const interval_t interval( std::get<0>(intervals[i]) , std::get<1>(intervals[i]) );
	  
	  // Get data
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  const std::vector<double> * data = slice.pdata();
	  
	  const std::vector<uint64_t> * tp = slice.ptimepoints();
	  
	  int nrows_per_interval = tp->size();
	  
	  // Populate signals
	  
	  for (int r=0;r<nrows_per_interval;r++)
	    {
	      
	      // Only add T and annotations once
	      if ( first )
		{		    
		  
		  // elapsed time in seconds
		  if (time_track) 
		    X(row,0) = (*tp)[r] * globals::tp_duration ;
		  
		  // Annotations (0/1) E,S 		  
		  int a_col = (int)time_track;	// start at 1 or 0
		  
		  std::map<std::string,int>::const_iterator aa = atype.begin();
		  while ( aa != atype.end() )
		    {		      
		      if ( aa->second == 0 )
			X(row,a_col) = std::numeric_limits<double>::quiet_NaN();
		      else if ( aa->second == 1 )
			{
			  // get exact point      
			  interval_t interval2 = interval_t( (*tp)[r] , (*tp)[r] + 1LLU );
			  annot_t * annot = edf.timeline.annotations( aa->first );
			  annot_map_t events = annot->extract( interval2 );
			  bool has_annot = events.size() ;
			  X(row,a_col) = (int)(has_annot ? 1 : 0 );
			}
		      // else if ( aa->second == 2 )
		      // 	X(row,a_col) = (int)( edf.timeline.epoch_annotation( aa->first , epoch0 ) ? 1 : 0 ) ;
		      
		      // next annotation
		      ++a_col;
		      ++aa;
		    }
		
		} // end of special case (T/ANNOTS)
	      	      
	      // Signal data	      
	      X(row,s_col) = (*data)[r];
	      
	      // next row 	      
	      ++row;
	    }

	  // Next interval
	}

      // no need to add T/ANNOTS again if looping back
      first = false;

      // advance to next signal
      ++s_col;
      
      // use Helper::sanitize() to make col names
    }
  
  return X;
  
}




std::vector<std::string> lunapi_inst_t::channels()
{
  std::vector<std::string> chs;
  if ( state != 1 ) return chs;
  signal_list_t signals = edf.header.signal_list( "*" );
  const int ns = signals.size();
  for (int s=0;s<ns;s++) 
    if ( edf.header.is_data_channel( signals(s) ) ) 
      chs.push_back( signals.label(s) );
  return chs; 
}  


std::vector<bool> lunapi_inst_t::has_channels( const std::vector<std::string> & chs )
{
  std::vector<bool> res;
  if ( state != 1 ) return res;
  res.resize( chs.size() );  
  const int ns = chs.size();
  for (int s=0;s<ns;s++)
    res[s] = edf.header.has_signal( chs[s] );
  return res;
}


std::vector<bool> lunapi_inst_t::has_annots( const std::vector<std::string> & anns )
{
  std::vector<bool> res;
  if ( state != 1 ) return res;
  res.resize( anns.size() );
  const int ns = anns.size();
  for (int s=0;s<ns;s++)
    res[s] = edf.timeline.annotations.find( anns[s] ) != NULL;
  return res;  
}

bool lunapi_inst_t::has_staging() 
{
  // get staging                                                                                                             
  edf.timeline.annotations.make_sleep_stage( edf.timeline );

  // valid?                                                                                                                  
  param_t empty_param;
  bool has_staging = edf.timeline.hypnogram.construct( &(edf.timeline) , empty_param , false );
  
  // valid, but empty?                                                                                                       
  if ( has_staging && edf.timeline.hypnogram.empty() )
    has_staging = false;
  
  return has_staging;
}

std::vector<std::string> lunapi_inst_t::annots() const
{
  if ( state != 1 ) return std::vector<std::string>(0);
  return edf.timeline.annotations.names();
}


//
// Inserters 
//

bool lunapi_inst_t::update_signal( const std::string & label , const std::vector<double> & x )
{
  if ( state != 1 ) return false;
  if ( ! edf.header.has_signal( label ) ) return false;
  const int slot = edf.header.signal( label );
  
  // void update_signal_retain_range( int s , const std::vector<double> * );

  // void update_signal( int s , const std::vector<double> * , int16_t * dmin = NULL , int16_t * dmax = NULL , 
  // 		      const double * pmin = NULL , const double * pmax = NULL );
  
  edf.update_signal( slot , &x ); 
  
  return true;
}


bool lunapi_inst_t::insert_signal( const std::string & label ,
				   const std::vector<double> & x ,
				   const int sr )
{
  if ( state != 1 ) return false;

  // edf.add_signal( const std::string & label , const int n_samples , const std::vector<double> & data ,
  //                  double pmin = 0 , double pmax = 0 ,
  //                  int16_t dmin = 0 , int16_t dmax = 0 );
  
  edf.add_signal( label , sr , x );
  
  return true;
}



lannot_t lunapi_inst_t::fetch_annots( const std::vector<std::string> & anns , const double interp ) const
{

  // interp is a special fix for scope - to chop stage annotations into units of no greater than 30
  // seconds, to make the plotting of hypnograms easier...  if -1 then ignore
  const bool do_interp = interp > 0 ;
  
  // typedef std::vector<std::tuple<std::string,double,double> > lannot_t;
  lannot_t r;
  if ( state != 1 ) return r;

  const int na = anns.size();
  for (int a=0; a<na; a++)
    {
      
      annot_t * annot = edf.timeline.annotations.find( anns[a] );
      if ( annot == NULL ) continue;      
      if ( annot->interval_events.size() == 0 ) continue;

      annot_map_t::const_iterator ii = annot->interval_events.begin();   
      while ( ii != annot->interval_events.end() )
	{
	  const instance_idx_t & instance_idx = ii->first;

	  if ( do_interp )
	    {
	      
	      uint64_t s = instance_idx.interval.start;
	      uint64_t w = interp * globals::tp_1sec;
	      
	      while ( 1 )
		{
		  // all done
		  if ( s >= instance_idx.interval.stop ) break;

		  // add right length (not past end)
		  const uint64_t s2 = s + w > instance_idx.interval.stop ? instance_idx.interval.stop : s + w ; 

		  // add
		  r.push_back( std::make_tuple( anns[a] , s * globals::tp_duration , s2 * globals::tp_duration ) );
		  
		  // next chunk
		  s += w;
		}
	      
	    }
	  else // add as a simple one
	    r.push_back( std::make_tuple( anns[a] ,
					  instance_idx.interval.start * globals::tp_duration,
					  instance_idx.interval.stop * globals::tp_duration ) );
	  ++ii;
	}
    }  
  
  return r;
}

lannot_full_t lunapi_inst_t::fetch_full_annots( const std::vector<std::string> & anns ) const
{
  // typedef std::vector<std::tuple<std::string,std::string,std::string,std::string,double,double> > lannot_full_t;
  lannot_full_t r;
  if ( state != 1 ) return r;
  
  const int na = anns.size();
  for (int a=0; a<na; a++)
    {
      
      annot_t * annot = edf.timeline.annotations.find( anns[a] );    
      if ( annot == NULL ) continue;      
      if ( annot->interval_events.size() == 0 ) continue;

      annot_map_t::const_iterator ii = annot->interval_events.begin();   
      while ( ii != annot->interval_events.end() )
	{
	  const instance_idx_t & instance_idx = ii->first;

	  const instance_t * inst = ii->second;
	  
	  std::string meta_data;
	  if ( inst->data.size() == 0 ) meta_data = ".";
	  else
	    {
	      std::stringstream ss;
	      std::map<std::string,avar_t*>::const_iterator dd = inst->data.begin();
	      while ( dd != inst->data.end() )
		{
		  if ( dd != inst->data.begin() ) ss << "|";
		  ss << *dd->second; // meta-data value
		  ++dd;		  
		}
	      meta_data = ss.str();
	    }

	  // add to return list
	  r.push_back( std::make_tuple( anns[a] ,
					( instance_idx.id == "" ? "." : instance_idx.id ) ,
					( instance_idx.ch_str == "" ? "." : instance_idx.ch_str ) ,
					meta_data , 
					instance_idx.interval.start * globals::tp_duration,
					instance_idx.interval.stop * globals::tp_duration ) );
	  ++ii;
	}
    }  
  
  return r;
}
				

bool lunapi_inst_t::insert_annotation( const std::string & class_label ,
				       const std::vector<std::tuple<double,double > > & x ,
				       const bool durcol2 )
{

  if ( state != 1 ) return false;
  if ( x.size() == 0 ) return false;
  if ( class_label == "" ) return false;
  
  const int n = x.size();

  // okay if class_label already exists, this will append new intervals
  annot_t * annot = edf.timeline.annotations.add( class_label );
  
  for (int i=0; i<n; i++ )
    {
      // skip bad elements
      if ( std::get<0>( x[i] ) < 0 || std::get<1>(x[i]) < 0 ) continue;
      
      const uint64_t start = std::get<0>(x[i]) * globals::tp_1sec ;
      const uint64_t stop  = std::get<1>(x[i]) * globals::tp_1sec + ( durcol2 ? start : 0LLU ) ; 
      
      annot->add( "." , // dummy instance ID
                  interval_t( start, stop ) ,
		  "." );  // channel ID dummy
    }

  return true;
}




