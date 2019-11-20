
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
#include <ftw.h>


int fn_luna_slbuilder(const char * fpath, const struct stat *ptr, int type )
{
  if (type == FTW_F)
    {
      std::string filename( fpath );
      
      // need to parse foldes?
      std::vector<std::string> tok;
      if ( globals::sl_link_across_folders || ( ! globals::sl_visit_edf ) )
	{
	  std::string dl( 1 , globals::folder_delimiter );	      
	  tok = Helper::parse( filename , dl );
	}
      
      if ( Helper::file_extension( filename , "edf" ) )
	{
	  
	  std::string id;
	  
	  // get ID from EDF? 
	  if ( globals::sl_visit_edf ) 
	    {
	      edf_t edf;
	      edf.attach( filename , "." );
	      id = edf.header.patient_id;
	    }
	  else
	    {	      
	      id = tok[tok.size()-1].substr( 0 , tok[ tok.size() - 1 ].size() - 4 );
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
	    tok[tok.size()-1].substr( 0 , tok[ tok.size() - 1 ].size() - 4 ) : 
	    filename.substr( 0 , filename.size() - 4 );
	  
	  globals::sl_data[ key ].id = id;
	  globals::sl_data[ key ].edf = filename;

	}
      else 
	{
	  std::set<std::string>::const_iterator ii = globals::sl_annot_extensions.begin();
	  while ( ii != globals::sl_annot_extensions.end() )
	    {

	      if ( Helper::file_extension( filename , *ii ) )
		{
		  
		  std::string key = globals::sl_link_across_folders ? 
		    tok[tok.size()-1].substr( 0 , tok[ tok.size() - 1 ].size() - 4 ) :
		    filename.substr( 0 , filename.size() - 4 );
		  
		  globals::sl_data[ key ].annots.insert( filename );
		  
		}
	      
	      ++ii;
	    }
	}
      
      
    }


  return 0;
}


void Helper::build_sample_list( const std::vector<std::string> & tok0 )
{
  // write a sample list to stdout
  // use ftw() to recursively scan folders for .edf and related files
  // optionally, look into EDFs for IDs
  // link same-named 

  globals::sl_visit_edf = true;
  globals::sl_link_across_folders = true;  // index just on file file, not full path

  // if there are distinct EDFs/ANNOTs with the same name, but in
  // different paths, this will fail; but if EDFs and ANNOTs are in
  // different folders, then this needs to be set (the default)
  
  std::vector<std::string> tok;
  
  for (int t=0;t<tok0.size();t++)
    {
      if ( tok0[t][0] == '-' ) 
	{
	  if ( tok0[t].substr(1) == "usefilename" ) globals::sl_visit_edf = false;
	  if ( tok0[t].substr(1) == "nospan" ) globals::sl_link_across_folders = false;
	  if ( tok0[t].substr(1,4) == "ext=" ) 
	    {
	      // -ext=txt,eannot,xls.eannot
	      std::string s = tok0[t].substr(5);
	      std::cerr << "s [" << s << "]\n";
	      std::vector<std::string> tok2 = Helper::parse( s , "," );
	      for (int i=0;i<tok2.size();i++) globals::sl_annot_extensions.insert( tok2[i] );
	    }
	}
      else tok.push_back( tok0[t] );
    }


  
  globals::sl_annot_extensions.insert( "xml" );
  
  for (int t=0;t<tok.size();t++)
    {      
      // expands any ~ home folder
      if ( ftw( Helper::expand( tok[t] ).c_str() , fn_luna_slbuilder, 10 ) != 0 ) 
	Helper::halt( "problem traversing folder " + tok[t] );      
    }
  

  //
  // Check 
  //
  
  std::vector<std::string> annot_wout_edf;
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
	  std::cerr << "*** warning *** ID " << ii->second.id << " is duplicated\n";
	}

      if ( has_edf ) 
	{
	  ids.insert( ii->second.id );
	  edf_annot_count[ num_annots ]++;
	  ++edfs; 
	}

      if ( has_edf ) 
	std::cout //<< ii->first << "\t" 
		  << ii->second.id << "\t"
		  << ii->second.edf ;
      
      
      std::set<std::string>::const_iterator jj = ii->second.annots.begin();
      while ( jj != ii->second.annots.end() )
	{
	  
	  if ( has_edf ) 
	    std::cout << "\t" << *jj;
	  else
	    annot_wout_edf.push_back( *jj );
	  
	  ++jj;
	}

      if ( has_edf ) 
	std::cout << "\n";

      ++ii;
    }  
  
  std::cerr << "\nwrote " << edfs << " EDFs to the sample list\n";

  std::map<int,int>::const_iterator kk = edf_annot_count.begin();
  while ( kk != edf_annot_count.end() )
    {
      std::cerr << "  " << kk->second << " of which had " << kk->first << " linked annotation files\n";
      ++kk;
    }
  
  if ( annot_wout_edf.size() > 0 ) 
    {
      std::cerr << "\nWarning: also found " 
		<< annot_wout_edf.size() 
		<< " annotation files without a matching EDF:\n";
      for (int i=0;i<annot_wout_edf.size() ; i++) 
	std::cerr << annot_wout_edf[i] << "\n";
    }

  if ( dupes ) 
    std::cerr << "\nWarning: duplicate IDs encountered\n";
}

