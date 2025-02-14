
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


#include "helper/helper.h"

#include "defs/defs.h"
#include "edf/edf.h"
#include <iostream>


int fn_luna_slbuilder( const std::string & filename );


//
// POSIX directory traversal function (called from ftw() below)
//

#ifndef WINDOWS

#include <ftw.h>

int fn_luna_slbuilder_ftw(const char * fpath, const struct stat *ptr, int type )
{
  if (type == FTW_F) 
    {
      std::string filename( fpath );
      return fn_luna_slbuilder( filename );
    }  
  return 0;
}


#else

//
// Windows API folder traveral 
//

#include <windows.h>

void fn_luna_slbuilder_win( const std::string & fpath )
{

  //  std::cout << "search [" << fpath << "]\n" ;
  
  WIN32_FIND_DATA FindFileData;
  
  HANDLE hFind = FindFirstFile( (fpath+"\\*").c_str() ,  &FindFileData);

  if ( hFind == INVALID_HANDLE_VALUE )
    {
      Helper::halt( "search failed on " + fpath ) ;
    }
  
  //
  // look at all results
  //

  bool cont = true;

  while ( cont )
    {
      
      std::string fname =  FindFileData.cFileName ;
      
      //std::cout << "do... [" <<  fname << "]\n";

      if ( fname != "." && fname != ".." )
	{
	  
	  // directory? if so, search (recursively)
	  
	  if ( FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
	    {	      
	      //std::cout << "folder... [" <<  fname << "]\n";
	      fn_luna_slbuilder_win( fpath + "\\" + fname );
	    }
      
	  // else, process as a file
	  
	  else
	    {
	      //std::cout << "testing " << fname << "\n";
	      fn_luna_slbuilder( fpath + "\\" + fname );
	    }
	  
	  
	}

      cont = FindNextFile( hFind , &FindFileData );
      
    }


  if ( GetLastError() != ERROR_NO_MORE_FILES )
    Helper::halt( "FindNextFile failed" );
  
  
  if ( ! FindClose(hFind) ) Helper::halt( "FindClose failed" );
  

}




#endif



//
// Generic code
//

int fn_luna_slbuilder( const std::string & filename )
{
 
  // need to parse foldes?
  std::vector<std::string> tok;
  if ( globals::sl_link_across_folders || ( ! globals::sl_visit_edf ) )
    {
      std::string dl( 1 , globals::folder_delimiter );	      
      tok = Helper::parse( filename , dl );
    }
  
  const bool has_edf_ext = Helper::file_extension( filename , "edf" ) ; // 4 chars
  const bool has_edfz_ext = Helper::file_extension( filename , "edfz" ) ; // 5 
  const bool has_edfgz_ext = Helper::file_extension( filename , "edf.gz" ) ; // 7  
  
  if ( has_edf_ext || has_edfz_ext || has_edfgz_ext ) 
    {
      
      const int xchar = has_edf_ext ? 4 : ( has_edfz_ext ? 5 : 7 ) ;

      std::string id;
      
      // get ID from EDF? 
      if ( globals::sl_visit_edf ) 
	{
	  edf_t edf;
	  edf.attach( filename , "." );
	  id = edf.header.patient_id;
	  if ( id == "" ) 
	    {
	      std::cerr << " *** empty Patient ID header for " << filename << ", so going to set ID to filename\n";
	      id = tok[tok.size()-1].substr( 0 , tok[ tok.size() - 1 ].size() - xchar );
	    }
	}
      else
	{	      
	  id = tok[tok.size()-1].substr( 0 , tok[ tok.size() - 1 ].size() - xchar );
	}
      
      //
      // store EDF/ID pair (locked on either full of local )
      //   for EDF:   /path/to/file.edf
      
      //   either link EDF to ANNOTs based on ; 
      //       /path/to/file    ( if -nospan ) 
      //   or:
      //       file             ( the default) 
      //   i.e. the latter will match annots and EDFs in different folders
      //        but it will also clash if different EDFs have same names but are in different folders
      
      std::string key = globals::sl_link_across_folders ? 
	tok[tok.size()-1].substr( 0 , tok[ tok.size() - 1 ].size() - xchar ) : 
	filename.substr( 0 , filename.size() - xchar );
      
      globals::sl_data[ key ].id = id;
      globals::sl_data[ key ].edf = filename;
      
    }
  else 
    {
      std::set<std::string>::const_iterator ii = globals::sl_annot_extensions.begin();
      while ( ii != globals::sl_annot_extensions.end() )
	{
	  
	  bool match_with_period = Helper::file_extension( filename , *ii , true );
	  bool match_no_period = Helper::file_extension( filename , *ii , false );
	  
	  if ( match_with_period || match_no_period )
	    {

	      int len = ii->size() + ( match_with_period ? 1 : 0 ); 
	      
	      std::string key = globals::sl_link_across_folders ? 
		tok[tok.size()-1].substr( 0 , tok[ tok.size() - 1 ].size() - len ) :
		filename.substr( 0 , filename.size() - len );

	      globals::sl_data[ key ].annots.insert( filename );
	      
	    }
	  
	  ++ii;
	}
    }
  
  return 0;

}


void Helper::build_sample_list( const std::vector<std::string> & tok0 , slist_t * slist )
{
  // write a sample list to stdout
  // use ftw() to recursively scan folders for .edf and related files
  // optionally, look into EDFs for IDs
  // link same-named 

  globals::sl_visit_edf = false;
  globals::sl_link_across_folders = true;  // index just on file file, not full path

  // if there are distinct EDFs/ANNOTs with the same name, but in
  // different paths, this will fail; but if EDFs and ANNOTs are in
  // different folders, then this needs to be set (the default)

  // optionally save to slist is slist is non-NULL (and do not write to std::cout)
  // (this is when called internally from lunapi_t)

  const bool write_cout = slist == NULL;
  
  std::vector<std::string> tok;
 
  bool specified_extensions = false;

  bool show_path = true;

  bool allow_numeric_ids = false;

  // added in ID is a number (i.e. better to have string IDs in Luna)
  const std::string id_prefix = "id_";
  
  for (int t=0;t<tok0.size();t++)
    {
      if ( tok0[t][0] == '-' ) 
	{

	  // use ID from EDF header, rather than file name
	  if ( tok0[t].substr(1) == "edfid" ) globals::sl_visit_edf = true;

	  // only link EDF and annotation files within the same folder
	  if ( tok0[t].substr(1) == "nospan" ) globals::sl_link_across_folders = false;

	  // allow numeric IDs?
	  if ( tok0[t].substr(1) == "allow-numeric-ids" ) allow_numeric_ids = true;
	  
	  // special case for -nsrr.xml extensions
	  if ( tok0[t].substr(1) == "nsrr" ) {	    
	    globals::sl_annot_extensions.insert( "-nsrr.xml" ); 
	    specified_extensions = true;
	  }

	  // do not show path in s.lst
	  if (  tok0[t].substr(1) == "rel" ) {
	    show_path = false;
	  }
	  
	  // other user-defined extensions
	  if ( tok0[t].substr(1,4) == "ext=" ) 
	    {
	      // -ext=txt,eannot,xls.eannot
	      std::string s = tok0[t].substr(5);
	      std::vector<std::string> tok2 = Helper::parse( s , "," );
	      for (int i=0;i<tok2.size();i++) globals::sl_annot_extensions.insert( tok2[i] );
	      specified_extensions = true;
	    }

	}
      else tok.push_back( tok0[t] );
    }
  
  //
  // standard annotation types, if none given above
  //

  if ( ! specified_extensions )
    {
      globals::sl_annot_extensions.insert( "xml" );
      globals::sl_annot_extensions.insert( "annot" );
      globals::sl_annot_extensions.insert( "eannot" );
      globals::sl_annot_extensions.insert( "txt" );
      globals::sl_annot_extensions.insert( "tsv" );
    }

  //
  // traverse folders
  //

  for (int t=0;t<tok.size();t++)
    {      
      
#ifndef WINDOWS
      // expands any ~ home folder
      if ( ftw( Helper::expand( tok[t] ).c_str() , fn_luna_slbuilder_ftw, 10 ) != 0 ) 
	Helper::halt( "problem traversing folder " + tok[t] );      
#else
      fn_luna_slbuilder_win( tok[t] );
#endif      
    }
  

  //
  // Check 
  //
  
  std::vector<std::string> annot_wout_edf; // annots w/out EDFs
  std::set<std::string> dumped_annot; // annots we wrote
  std::map<int,int> edf_annot_count;
  int edfs = 0;
  
  std::set<std::string> ids;
  
  //
  // Write sample list
  //
  
  bool dupes = false;

  std::map<std::string,sample_list_t>::const_iterator ii = globals::sl_data.begin();
  while ( ii != globals::sl_data.end() )
    {
      // does this have an EDF? 
      
      bool has_edf = ii->second.edf != "";
      int  num_annots = ii->second.annots.size();
      
      if ( has_edf && ids.find( ii->second.id ) != ids.end() ) 
	{
	  dupes = true;
	  logger << "*** warning *** ID " << ii->second.id << " is duplicated\n";
	}

      if ( has_edf ) 
	{
	  ids.insert( ii->second.id );
	  edf_annot_count[ num_annots ]++;
	  ++edfs; 
	}

      //
      // output
      //

      std::string exp_id = "";
      std::string exp_edf = "";
      std::set<std::string> exp_annots;
      
      if ( has_edf )
	{

	  //
	  // ID
	  //
	  
	  if ( allow_numeric_ids ) 
	    {
	      if ( write_cout )
		std::cout << ii->second.id << "\t";
	      else
		exp_id = ii->second.id;
	    }
	  else
	    {
	      double num_test; // if this cases to a valiud number

	      bool is_numeric = Helper::str2dbl( ii->second.id , &num_test ) ;

	      // but allow 123-567, etc ; so test for - _ as common delimiters 
	      if      ( ii->second.id.find( "-" ) != std::string::npos ) is_numeric = false;
	      else if ( ii->second.id.find( "_" ) != std::string::npos ) is_numeric = false;
	      else if ( ii->second.id.find( ":" ) != std::string::npos ) is_numeric = false;
	      else if ( ii->second.id.find( "+" ) != std::string::npos ) is_numeric = false;

	      if ( is_numeric ) 
		{
		  if ( write_cout )
		    std::cout << id_prefix << ii->second.id << "\t";
		  else
		    exp_id = id_prefix + ii->second.id;
		}
	      else
		{
		  if ( write_cout )
		    std::cout << ii->second.id << "\t";
		  else
		    exp_id = ii->second.id;
		}
	    }

	  //
	  // EDF
	  //

	  if ( write_cout )
	    std::cout << ii->second.edf   // was Helper::quote_spaced( ii->second.edf )
		      << "\t";
	  else
	    exp_edf = ii->second.edf;
	}

      //
      // Annots
      //

      // new sl-format: can add '.' in third slot if no annots
      if ( ii->second.annots.size() == 0 )
	{
	  if ( write_cout )
	    std::cout << ".";
	  // else just return empty exp_annots set
	}

      // new sl-format: can add comma (or globals::file_list_delimiter) separated values in one tab slot
      bool first = true;
      std::set<std::string>::const_iterator jj = ii->second.annots.begin();
      while ( jj != ii->second.annots.end() )
	{	  
	  if ( has_edf ) 
	    {
	      if ( write_cout )
		{
		  if ( ! first ) std::cout << globals::file_list_delimiter;
		  std::cout << *jj ; // was  Helper::quote_spaced( *jj ) ;
		}
	      else
		{
		  exp_annots.insert( *jj );
		}
	      
	      first = false;
	      dumped_annot.insert( *jj );
	    }
	  else
	    annot_wout_edf.push_back( *jj );	  
	  ++jj;
	}
      
      if ( has_edf ) 
	if ( write_cout )
	  std::cout << "\n";
      
      
      //
      // save to slist?
      //
      
      if ( ! write_cout ) 
	{
	  slist->push_back( std::make_tuple( exp_id , exp_edf , exp_annots ) );
	}
      

      //
      // Next 
      //

      ++ii;
    }



  //
  // extra outputs?
  //
  
  if ( write_cout )
    {
      std::cerr << "\nwrote " << edfs << " EDFs to the sample list\n";
      
      std::map<int,int>::const_iterator kk = edf_annot_count.begin();
      while ( kk != edf_annot_count.end() )
	{
	  std::cerr << "  " << kk->second << " of which had " << kk->first << " linked annotation files\n";
	  ++kk;
	}
      
      if ( annot_wout_edf.size() > 0 ) 
	{
	  // need to exclude e.g.
	  //     file.edf    file.eannot
	  //  --> this will think there is a   'file.e.edf'  w/out ?
	  
	  std::vector<std::string> annot_wout_edf2;
	  for (int i=0;i<annot_wout_edf.size() ; i++)
	    if ( dumped_annot.find( annot_wout_edf[i] ) == dumped_annot.end() )
	      annot_wout_edf2.push_back( annot_wout_edf[i] );
	  
	  if ( annot_wout_edf2.size() > 0 )
	    {
	      std::cerr << "\nWarning: also found " 
			<< annot_wout_edf2.size() 
			<< " annotation files without a matching EDF:\n";
	      
	      for (int i=0;i<annot_wout_edf2.size() ; i++) 
		std::cerr << annot_wout_edf2[i] << "\n";
	    }
	}
      
      if ( dupes ) 
	std::cerr << "\nWarning: duplicate IDs encountered\n";
    }
  
  
}

