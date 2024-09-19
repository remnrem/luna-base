
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

#include "helper/validate.h"

#include "defs/defs.h"
#include "edf/edf.h"
#include "db/db.h"
#include "helper/logger.h"
#include "helper/helper.h"

#include <fstream>
#include <cstdio>
#include <dirent.h>

extern writer_t writer;
extern logger_t logger;

void Helper::validate_slist( param_t & param )
{

  
  // take an existing sample list
  // check each EDF and all annotations
  // output if we cannot read
  // generate an exclude list

  globals::validation_mode = true;
    
  const std::string slist = Helper::expand( param.requires( "slist" ) );

  logger << "  validating files in sample list " << slist << "\n\n";  

  if ( ! Helper::fileExists( slist ) ) Helper::halt( "could not open sample-list " + slist );

  std::ifstream IN1( slist.c_str() , std::ios::in );


  //
  // trackers
  //

  std::set<std::string> exclude_edf;
  std::map<std::string,std::set<std::string>> exclude_annots;
  std::set<std::string> exclude;  
  
  //
  // set up 
  //
  
  const bool has_project_path = globals::param.has( "path" );
  
  int goodn = 0 , badn = 0;
  
  while ( 1 )
    {

      //
      // clear any flags      
      //

      globals::problem = false;

      //
      // read SL
      //

      std::string line;
      Helper::safe_getline( IN1 , line );

      // all done?
      if ( IN1.eof() || IN1.bad() ) break;

      if ( line == "" ) continue;


      //
      // get SL line
      //
      
      std::vector<std::string> tok = Helper::parse( line , "\t" );

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

      //
      // extract main items (ID, signal EDF)
      //
      
      std::string rootname = tok[0];
      std::string edffile  = tok[1];

      //
      // swap in new ID?
      //

      rootname = cmd_t::remap_id( rootname );
      
      
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
	  
	  continue; // to the next EDF in the list
	}


      //
      // set up writer      
      //
      
      writer.id( rootname , "." );

      //
      // try to load 
      //

      edf_t edf;

      bool edf_okay = edf.attach( edffile , rootname , NULL , true );

      if ( ! edf_okay )
	exclude_edf.insert( rootname );
	    
      
      // ------------------------------------------------------------
      //
      // try to load annotations, but only for good annots
      //
      // ------------------------------------------------------------

      //
      // init an empty EDF in case the above was left in a weird state
      //

      // if attaching an .eannot, we need the right number of records, etc
      // to get impllied epoch length, however

      // if bad EDF< default = 24 hr empty EDF, although this should not matter
      // (other than for .eannots)
      
      const int nr = edf_okay ? edf.header.nr : 24 * 60 ; 
      const int rs = edf_okay ? edf.header.record_duration : 60 ;
      const std::string startdate = edf_okay ? edf.header.startdate : "01.01.00" ;
      const std::string starttime = edf_okay ? edf.header.starttime : "00.00.00" ;
      const std::string id = edf_okay ? rootname : "__bad_EDF__";
      
      edf_t dummy;

      bool empty_okay = dummy.init_empty( id , nr , rs , startdate , starttime );

      if ( ! empty_okay )
	Helper::halt( "internal error constructing an empty EDF to evaluate annotations" );
      
     
      // some basic set-up

      dummy.timeline.annotations.set( &dummy );
      

      // Add additional annotations? (outside of slist)
      // not sure this is needed here

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
			      
			      bool okay = dummy.load_annotations( fname + fname2 );

			      if ( ! okay )
				{
				  exclude_annots[ rootname ].insert( fname + fname2 );
				  exclude.insert( rootname );
				}
			    }			 
			}
		      closedir (dir);
		    }
		  else
		    {
		      Helper::vmode_halt( "could not open folder " + fname );
		      exclude_annots[ rootname ].insert( fname );		      
		      exclude.insert( rootname );
		    }
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
		      bool okay = dummy.load_annotations( fname );
		      if ( ! okay )
			{
			  exclude_annots[ rootname ].insert( fname );
			  exclude.insert( rootname );
			}
		    }
		  else
		    {
		      Helper::vmode_halt( "did not recognize annotation file extension: " + fname );
		      exclude_annots[ rootname ].insert( fname );
		      exclude.insert( rootname );
		    } 
		}
	    }
	}
      
      const bool annots_okay = exclude_annots.find( rootname ) == exclude_annots.end();
      
      writer.value( "ANNOTS" , (int)annots_okay );


      //
      // track
      //

      if ( edf_okay && annots_okay ) ++goodn;
      else // if either EDF or annots bad:
	{
	  ++badn;	  
	  exclude.insert( rootname );	  
	}
      writer.value( "EDF" , edf_okay );



      //
      // Final outputs
      //

      
      std::map<std::string,std::set<std::string>>::const_iterator ii = exclude_annots.find( rootname );
      if ( ii != exclude_annots.end() )
	{
	  const std::set<std::string> & f = ii->second;
	  std::set<std::string>::const_iterator ff = f.begin();
	  while ( ff != f.end() )
	    {
	      writer.level( *ff , "FILE" );
	      writer.value( "EXC" , 1 );
	      ++ff;
	    }
	  writer.unlevel( "FILE" );
	}

      std::set<std::string>::const_iterator jj = exclude_edf.find( rootname );
      if ( jj != exclude_edf.end() )
	{
	  writer.level( edffile , "FILE" );
	  writer.value( "EXC" , 1 );
	  writer.unlevel( "FILE" );
	  ++jj;
	}
      

      
      //
      // Next individual
      //
      
    }

  if ( badn ) 
    logger << "\n  " << badn << " of " <<  goodn + badn << " observations scanned had corrupt/missing EDF/annotation files\n";
  else
    logger << "  all good, no problems detected in " << goodn << " observations scanned\n";
    
  IN1.close();
  

  
  
  //
  // write exclude lists?
  //

  // any?
  if ( param.has( "exclude-list" ) && exclude.size() != 0 )
    {
      logger << "  writing exclude list (based on either EDF or annotation issues) to " << param.value( "exclude-list" ) << "\n";
      std::ofstream O1( Helper::expand( param.requires( "exclude-list" ) ).c_str() , std::ios::out );
      std::set<std::string>::const_iterator ii = exclude.begin();
      while ( ii != exclude.end() )
	{
	  O1 << *ii << "\n";
	  ++ii;
	}      
      O1.close();
    }

  // EDF only
  if ( param.has( "edf-exclude-list" ) && exclude_edf.size() != 0 )
    {
      logger << "  writing exclude list (based on EDF issues only) to " << param.value( "edf-exclude-list" ) << "\n";
      std::ofstream O1( Helper::expand( param.requires( "edf-exclude-list" ) ).c_str() , std::ios::out );      
      std::set<std::string>::const_iterator ii = exclude_edf.begin();
      while ( ii != exclude_edf.end() )
	{
	  O1 << *ii << "\n";
	  ++ii;
	}      
      O1.close();
    }

  // annots only
  if ( param.has( "annot-exclude-list" ) && exclude_annots.size() != 0 )
    {
      logger << "  writing exclude list (based on annotation EDF issues only) to " << param.value( "edf-exclude-list" ) << "\n";
      std::ofstream O1( Helper::expand( param.requires( "annot-exclude-list" ) ).c_str() , std::ios::out );
      std::map<std::string,std::set<std::string>>::const_iterator ii = exclude_annots.begin();
      while ( ii != exclude_annots.end() )
        {
          O1 << ii->first << "\n";
          ++ii;
        }
      O1.close();
    }

  
  //
  // all done
  //
  
  globals::validation_mode = false;

}


//
// because of the different interfaces, probably cleaner to split off this functionality (but keep here rather than lunapi/)
//


std::vector<std::tuple<std::string,std::string,bool> > Helper::validate_slist_lunapi_mode( const std::vector<std::tuple<std::string,std::string,std::set<std::string> > > & sl )
{
  
  
  std::vector<std::tuple<std::string,std::string,bool> > r ;

  globals::validation_mode = true;
  
  std::vector<std::tuple<std::string,std::string,std::set<std::string> > >::const_iterator ss = sl.begin();
  while ( ss != sl.end() )
    {
	
      //
      // get ID (possibly remapping) and file names
      //

      //const std::string rootname = cmd_t::remap_id( std::get<0>(*ss) );
      const std::string rootname = std::get<0>(*ss) ;
      const std::string edffile = std::get<1>(*ss);
      const std::set<std::string> annots = std::get<2>(*ss);

      
      //
      // clear any problem flags
      //

      globals::problem = false;

      //
      // include/exclude?
      //
      
      bool include = true;
      
      if ( globals::id_excludes.find( rootname ) != globals::id_excludes.end() )
        include = false;
      
      if ( globals::id_includes.size() != 0
           && globals::id_includes.find( rootname ) == globals::id_includes.end() )
        include = false;

      if ( ! include ) { ++ss; continue; }
      
      // else, do we have an 'ID' check? (id=ID does not match so skip)                                                                                                                             
      if ( globals::sample_list_ids.size() )
	{
	  if ( globals::sample_list_ids.find( rootname ) == globals::sample_list_ids.end() )
	    {
	      ++ss;
	      continue;
	    }
	}

      // skip=ID matches                                                                                                                                                                            
      if ( globals::sample_list_ids_skips.size() )
	{
	  if ( globals::sample_list_ids_skips.find( rootname ) != globals::sample_list_ids_skips.end() )
	    {
	      ++ss;
	      continue;
	    }
	}
      
      
      //
      // try EDF
      //

      edf_t edf;
      
      bool edf_okay = edf.attach( edffile , rootname , NULL , true );

      r.push_back( std::make_tuple( rootname , edffile , edf_okay ) );
      

      //
      // try to load annotations, but only for good annots
      //   -- see above for logic
      //
      
      const int nr = edf_okay ? edf.header.nr : 24 * 60 ; 
      const int rs = edf_okay ? edf.header.record_duration : 60 ;
      const std::string startdate = edf_okay ? edf.header.startdate : "01.01.00" ;
      const std::string starttime = edf_okay ? edf.header.starttime : "00.00.00" ;
      const std::string id = edf_okay ? rootname : "__bad_EDF__";
      
      edf_t dummy;
      
      bool empty_okay = dummy.init_empty( id , nr , rs , startdate , starttime );
      
      if ( ! empty_okay )
	Helper::halt( "internal error constructing an empty EDF to evaluate annotations" );
      
     
      // some basic set-up

      dummy.timeline.annotations.set( &dummy );

     
      //
      // try anntatioons
      //
      
      if ( ! globals::skip_nonedf_annots ) 
	{

	  std::set<std::string>::const_iterator aa = annots.begin();

	  while ( aa != annots.end() )
	    {
	      
	      std::string fname = Helper::expand( *aa );
	      
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
			      
			      bool okay = dummy.load_annotations( fname + fname2 );

			      // track
			      r.push_back( std::make_tuple( rootname , fname + fname2 , okay ) );

			    }			 
			}
		      closedir (dir);
		    }
		  else
		    {
		      Helper::vmode_halt( "could not open folder " + fname );

		      // track                                                                                                                                                                  
		      r.push_back( std::make_tuple( rootname , fname , false ) );
		    }
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
		      bool okay = dummy.load_annotations( fname );

		      // track
		      r.push_back( std::make_tuple( rootname , fname , okay ) );

		    }
		  else
		    {
		      Helper::vmode_halt( "did not recognize annotation file extension: " + fname );
		      // track 
                      r.push_back( std::make_tuple( rootname , fname , false ) );
		    } 
		}

	      // next annotation
	      ++aa;
	    }
	}
      
      //
      // Next individual
      //

      ++ss;
    }

  globals::validation_mode = false;
  
  return r;
}
