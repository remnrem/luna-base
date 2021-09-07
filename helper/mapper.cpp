
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

#include "helper/mapper.h"
#include "annot/nsrr-remap.h"

extern globals global;


void Helper::channel_annot_mapper( const std::vector<std::string> & tok , bool html )
{
  
  const int n = tok.size(); 

  // remapping files
  std::string cmap = "" , amap = "";
  
  // canonical channel file
  std::vector<std::string> csfiles;
  
  // channels/annots to be remapped
  std::map<std::string,std::string> anns;
  std::vector<std::string> nonuniq_chs, chs;
  std::vector<bool> anchan;
  
  for (int i=0; i<n; i++)
    {
      std::vector<std::string> tok2 = Helper::quoted_parse( tok[i] , "=" );
      
      // ignore bad stuff
      if ( tok2.size() != 2 ) continue;
      
      if ( tok2[0] == "cmap" ) cmap = Helper::expand( tok2[1] );
      else if ( tok2[0] == "amap" ) amap = Helper::expand( tok2[1] );
      else if ( tok2[0] == "csfile" ) csfiles.push_back( Helper::expand( tok2[1] ) );
      else if ( tok2[0] == "c" )
	{
	  std::vector<std::string> tok3 = Helper::quoted_parse( tok2[1] , "," );
	  for (int j=0; j<tok3.size(); j++)
	    nonuniq_chs.push_back( tok3[j] );
	}
      else if ( tok2[0] == "a" )
	{
	  std::vector<std::string> tok3 = Helper::quoted_parse( tok2[1] , "," );
	  for (int j=0; j<tok3.size(); j++)
	    anns[ tok3[j] ] = "";
	}
            
    }

  const bool do_amap = Helper::fileExists( amap );
  const bool do_cmap = Helper::fileExists( cmap );
  const bool do_cansigs = csfiles.size() == 0 || Helper::fileExists( csfiles[0] );
  
  
  //
  // expecting canonical cmap and amap to be in simple form, i.e.  
  //

  if ( cmap != "" &&  ! Helper::fileExists( cmap ) ) cmap = "";
  if ( amap != "" &&  ! Helper::fileExists( amap ) ) amap = "";
  
  if ( cmap == "" && amap == "" ) 
    {
      std::cout << "no mappings given, quitting\n";
      return;
    }



  //
  // Annotations + nsrr-remap --> F 
  //

  if ( do_amap && anns.size() > 0 ) 
    {
      std::ifstream INC( amap.c_str() , std::ios::in );
      if ( INC.bad() ) Helper::halt( "could not open file: " + amap );
      
      while ( ! INC.eof() )
        {
          std::string line;
          Helper::safe_getline( INC , line );
          if ( INC.eof() || line == "" ) continue;
	  
          // skip % comments and, here, conditionals                                                                                        
          if ( line[0] == '%' ) continue;
          if ( line[0] == '+' || line[0] == '-' ) continue;
	  
          // otherwise parse as a normal line: i.e. two tab-delim cols                                                                      
          std::vector<std::string> tok = Helper::quoted_parse( line , "\t" );
          if ( tok.size() != 2 ) continue;
	  
	  // allow nsrr-remap command to be executed - i.e. to clear any internal map first
	  if ( tok[0] != "remap" && tok[0] != "nsrr-remap" ) continue;
	  
	  cmd_t::parse_special( tok[0] , tok[1] );
        }

      INC.close();
    }


  //
  // Channels
  //
  
  if ( do_cmap && nonuniq_chs.size() > 0 )
    {
      std::ifstream INC( cmap.c_str() , std::ios::in );
      if ( INC.bad() ) Helper::halt( "could not open file: " + cmap );
      
      while ( ! INC.eof() )
        {                     
          std::string line;       
          Helper::safe_getline( INC , line );                 
          if ( INC.eof() || line == "" ) continue;
                      
          // skip % comments and, here, conditionals
          if ( line[0] == '%' ) continue;
          if ( line[0] == '+' || line[0] == '-' ) continue;
          
          // otherwise parse as a normal line: i.e. two tab-delim cols
          std::vector<std::string> tok = Helper::quoted_parse( line , "\t" );
          if ( tok.size() != 2 ) continue;      
          if ( tok[0] != "alias" ) continue;
          cmd_t::parse_special( tok[0] , tok[1] );
        }
      
      INC.close();
    }


  

  //
  // Map annotations
  //
  
  // ensure we only report things that are mapped:
  nsrr_t::whitelist = true;
  
  std::map<std::string,std::string>::iterator aa = anns.begin();
  while ( aa != anns.end() )
    {
      
      // remap? (and if so, track)
      std::string y = nsrr_t::remap( aa->first ) ;
      
      // not white-listed?
      if ( y != "" ) aa->second = y;
      
      ++aa;
    }


  
  //
  // Map channels --> aliases
  //
  
  std::set<std::string> slabels;

  std::vector<std::string>::const_iterator cc = nonuniq_chs.begin();
  while ( cc != nonuniq_chs.end() )
    {
      
      // signal label, trim leading/trailing spaces
      std::string l = Helper::trim( *cc );
      
      // annotation?
      if ( Helper::imatch( l , "EDF Annotation" , 14 ) ) {
	chs.push_back( *cc );
	anchan.push_back( true );
	++cc;
	continue;
      } 

      // swap internal spaces?
      if ( globals::replace_channel_spaces )
   	l = Helper::search_replace( l , ' ' , globals::space_replacement );
      
      // make all data-channels upper case?
      if ( globals::uppercase_channels )
   	l = Helper::toupper( l );

      // key on UC version
      std::string uc_l = Helper::toupper( l );

      // find any alias?
      std::map<std::string,std::string>::const_iterator ff = cmd_t::label_aliases.find( uc_l );
      
      // either store alias, or "cleaned" version (spaces, uniqify)
      if ( ff != cmd_t::label_aliases.end() )
	{
	  l = ff->second;
	  uc_l = Helper::toupper( l );
	}

      
      // does this exist already? if so, uniqify 
      if ( slabels.find( uc_l ) != slabels.end() )
   	{
   	  int inc = 1;
   	  while ( 1 ) 
   	    {
   	      // new unique label?
   	      if ( slabels.find( uc_l + "." + Helper::int2str( inc )  ) == slabels.end() )
   		{
   		  l = l + "." + Helper::int2str( inc );
   		  uc_l = uc_l + "." + Helper::int2str( inc );
   		  break;
   		}
   	      else // keep trying
   		++inc;
   	    }
	}
      else
	slabels.insert( uc_l );
      
      //
      // store
      //

      chs.push_back( l );
      anchan.push_back( false );
      // next
      ++cc;
    }

  

  //
  // Construct a templete EDF header, and do any canonical mappings
  //

  edf_t edf;
  
  edf.header.nr = 0 ;
  edf.header.nr_all = 0;

  // set number of channels
  edf.header.ns = edf.header.ns_all = chs.size();
  edf.header.label.resize( edf.header.ns_all );
  edf.header.annotation_channel.resize( edf.header.ns_all );

  for (int c=0; c<chs.size(); c++)
    {
      edf.header.label2header[ Helper::toupper( chs[c] ) ] = c;
      edf.header.label[c] = chs[c];

      // kludge: signal_list() in make_canonical() checks
      // whether the signal is an annotation or not, so we
      // need to specify that
      edf.header.annotation_channel[c] = anchan[c];
    }


  //
  // For this dummy EDF, see which canonical signals can be made
  //

  cansigs_t cs0;
  
  if ( do_cansigs )
    {
      
      // csfiles --> definitions (files)
      // '.'     --> no special group; can add this ability later
      // false   --> do not make new signals      
      // false   --> do not drop originals
      // NULL    --> attempt all CS in the file
      // true    --> only check labels (i.e. do not even look at SR etc) 
      
      cs0 = edf.make_canonicals( csfiles , "." , false , false , "" , NULL , true ); 
      
    }
  
  //
  // Final report
  // 
  
  aa = anns.begin();
  while ( aa != anns.end() )
    {
      
      std::cout << "annot" << "\t"
		<< aa->first << "\t"
		<< ( aa->second == "" ? "-unmapped-" : aa->second ) << "\n";
      ++aa;
    }

  const int nc = chs.size();
  for ( int c=0; c<nc; c++)
    {
      const bool used = cs0.used.find( chs[c] ) != cs0.used.end();
      
      std::cout << "ch" << "\t"
		<< nonuniq_chs[c] << "\t"
		<< chs[c] << "\t"
		<<  ( used ? "Y" : "N" ) 
		<< "\n";
    }

  //
  // mapped canonical signals:
  //

    
  std::map<std::string,bool>::const_iterator ss = cs0.okay.begin();
  while ( ss != cs0.okay.end() )
    {
      
      if ( ss->second )
	{
	  std::cout << ss->first << "\t"
		    << cs0.sig[ ss->first ] ;
	  
	  if ( cs0.ref[ ss->first ] != "." ) 
	    std::cout << " - " << cs0.ref[ ss->first ] ;
	  
	  std::cout << "\n";
	}
      
      ++ss;
    }
  
  
}
