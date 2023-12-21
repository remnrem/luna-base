
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

#include "lunapi/lunapi.h"

extern logger_t logger;

extern globals global;

void lunapi_bail_function( const std::string & msg )
{
  logger << "*** [lunapi] error: " << msg << "\n";
  return;
}
  

void lunapi_t::init()
{

  global.init_defs();

  globals::bail_function = &lunapi_bail_function;

  //globals::logger_function = &lunapi_message_function;

  // turn off external output
  global.api();
  
  logger << "** lunapi " 
	 << globals::version 
	 << " " << globals::date
	 << "\n";
  
}


void lunapi_t::reset()
{
  globals::problem = false;
  globals::empty = false;
}



bool lunapi_t::attach_edf( const std::string & filename )
{
  if ( ! Helper::fileExists( filename ) ) 
    Helper::halt( "cannot find " + filename );

  bool okay = edf.attach( filename , id , NULL );

  if ( ! okay )
    {
      state = -1;
      return false;
    }

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

  return true;
  
  return okay;
  
}


bool lunapi_t::attach_annot( const std::string & annotfile )
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
                }
            }
          closedir (dir);
        }
      else
        {
          Helper::halt( "could not open folder " + annotfile );
        }
    }
  
  //
  // else a single file, load it
  //
  else
    {
      edf.load_annotations( annotfile );
    }
  
  return true;

}

retval_t lunapi_t::eval( const std::string & cmdstr )
{

  //
  // set up retval_t mechanism to catch outputs
  //

  retval_t ret;
  
  writer.clear();
  writer.set_types(); // not sure this is needed now...
  writer.use_retval( &ret );

  //
  // set ID 
  //

  writer.id( id , filename );

  logger << "evaluating...\n";

  // track errors? (for moonlight) 
  // R_last_eval_failed = false;
  // R_last_eval_errmsg = "";

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

  writer.use_retval( NULL );


  writer.clear();
  writer.set_types();

  //
  // was a problem flag set?
  //

  if ( globals::problem ) 
    Helper::halt( "problem flag set: likely no unmasked records left?" );

  
  //
  // all done
  //

  
  return ret;
  
}
