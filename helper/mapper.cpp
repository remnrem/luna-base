
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
#include "edf/canonical.h"

extern globals global;

void Helper::channel_annot_mapper( const std::vector<std::string> & tok , bool html )
{

  // input
  const int n = tok.size(); 

  // for (int i=0; i<n; i++)
  //   std::cout << " tok[" << i << "] = [" << tok[i] << "]<br>";
  // std::cout << "</p>";
  
  // remapping files
  std::string cmap = "" , amap = "";

  // Or, can also taken alias/remap instructions on-the-fly remppings
  //  tok =  alias="primary"|"secondary"
  //  tok =  remap="primary"|"secondary"
  // canonical channel file: harmonized EDF, then base EDF

  std::vector<std::string> csfiles_harm, csfiles_base;
  
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
      else if ( tok2[0] == "cs-harm" ) csfiles_harm.push_back( Helper::expand( tok2[1] ) );
      else if ( tok2[0] == "cs-base" ) csfiles_base.push_back( Helper::expand( tok2[1] ) );      
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
  const bool do_cansigs_harm = csfiles_harm.size() == 0 || Helper::fileExists( csfiles_harm[0] );
  const bool do_cansigs_base = csfiles_base.size() == 0 || Helper::fileExists( csfiles_base[0] );
  
  
  //
  // expecting canonical cmap and amap to be in simple form, i.e.  
  //

  if ( cmap != "" &&  ! Helper::fileExists( cmap ) ) cmap = "";
  if ( amap != "" &&  ! Helper::fileExists( amap ) ) amap = "";
  
  // if ( cmap == "" && amap == "" ) 
  //   {
  //     std::cout << "no mappings given, quitting\n";
  //     return;
  //   }



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
	  
	  // add an annotation mapping?
	  if ( tok[0] != "remap" ) continue;
	  
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
  // Any additional, user-specified mappings? (only add if valid)
  //

  for (int i=0; i<n; i++)
    {
      std::vector<std::string> tok2 = Helper::quoted_parse( tok[i] , "=" );      
      if ( tok2.size() != 2 ) continue;      

      if ( tok2[0] == "alias" )
	{
	  cmd_t::parse_special( tok2[0] , tok2[1] );      
	}
      else if ( tok2[0] == "remap" )
	{
	  cmd_t::parse_special( tok2[0] , tok2[1] );

	  // std::vector<std::string> tt = Helper::quoted_parse( tok2[1] , "|" );

	  // if ( tt.size() == 2 )
	  //   {
	  //     // std::string from = Helper::unquote( Helper::toupper( tt[1] ) );
	  //     // std::string to =  Helper::unquote( Helper::toupper( tt[0] ) );

	  //     // bool okay = true;
	  
 	  //     // std::cout << " testing " << from
	  //     // 		<< " " << to << "  sz = " 
	  //     // 		<< nsrr_t::pmap.size() << "</p>";
	      
	  //     // if ( nsrr_t::pmap.find( from ) != nsrr_t::pmap.end() )
	  //     // 	{
	  //     // 	  okay = false;
	  //     // 	  std::cout << "already  found " << from << " as a primary (to-var)</p>"; 
	  //     // 	}

	  //     // std::cout << " HERE HERE " << from << "</p>";
	      
	  //     //if ( okay )
	      
	  //     // else
	  //     // 	std::cout << " not okay " << from << " " << to << "\n";
	  //   }
	}
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
  // For this dummy EDF, see which canonical signals ( orig --> HARM ) 
  //
  
  cansigs_t cs0;
  
  if ( do_cansigs_harm )
    {
      
      // csfiles --> definitions (files)
      // '.'     --> no special group; can add this ability later
      // false   --> do not make new signals      
      // false   --> do not drop originals
      // NULL    --> attempt all CS in the file
      // true    --> only check labels (i.e. do not even look at SR etc) 
      
      cs0 = edf.make_canonicals( csfiles_harm , "." , false , false , "" , NULL , true ); 
      
    }

  int ns_harm = 0;
  std::map<std::string,bool>::const_iterator ss = cs0.okay.begin();
  while ( ss != cs0.okay.end() )
    {
      if ( ss->second ) ++ns_harm;
      ++ss;
    }

  
  //
  // Make a second EDF (harmonized) --> base EDF
  //
  
  edf_t edf1;
  edf1.header.nr = 0 ;
  edf1.header.nr_all = 0;

  // set number of channels
  edf1.header.ns = edf1.header.ns_all = ns_harm;
  edf1.header.label.resize( edf1.header.ns_all );
  edf1.header.annotation_channel.resize( edf1.header.ns_all );

  ss = cs0.okay.begin();
  int c = 0;
  while ( ss != cs0.okay.end() )
    {
      if ( ss->second )
	{
	  edf1.header.label2header[ Helper::toupper( ss->first ) ] = c;
	  edf1.header.label[c] = ss->first;
	  edf1.header.annotation_channel[c] = false ;
	  ++c;
	}
      ++ss;
    }
  

  //
  // For this dummy EDF, see which canonical signals ( orig --> HARM ) 
  //
  
  cansigs_t cs1;
  
  if ( do_cansigs_base )
    cs1 = edf1.make_canonicals( csfiles_base , "." , false , false , "" , NULL , true ); 


  //
  // Counts
  //

  const int nc = chs.size();
  const int na = anns.size();

  int a_mapped = 0;
  int a_aliased = 0;
    
  int c_aliased = 0;
  int c_harm = 0;
  int c_base = 0;
  int c_used = 0;
 
  aa = anns.begin();
  while ( aa != anns.end() )
    {
      if (  aa->second != "" )
	{
	  ++a_mapped;
	  if ( aa->second !=  aa->first)
	    ++a_aliased;
	}      
      ++aa;
    }

  for ( int c=0; c<nc; c++)
    {
      if ( nonuniq_chs[c] != chs[c] ) ++c_aliased;
      if ( cs0.used.find( Helper::toupper( chs[c] ) ) != cs0.used.end() ) ++c_used;
    }
  
  ss = cs0.okay.begin();
  while ( ss != cs0.okay.end() )
    {
      if ( ss->second ) ++c_harm;
      ++ss;
    }
  
  ss = cs1.okay.begin();
  while ( ss != cs1.okay.end() )
    {
      if ( ss->second ) ++c_base;
      ++ss;
    }
  
  
  //
  // Final report
  // 

  std::cout << "<table width=100%><tr><td width=50%>";
  
  std::cout << "<li>" << c_used << " of " << nc << " channels mapped, aliasing " << c_aliased ;
  std::cout << "<li>" << c_harm << " harmonized EDF channels constructed";
  std::cout << "<li>" << c_base << " base EDF channels constructed";

  if ( c_used < nc )
    {
      std::cout << "<li><b><em> " << nc-c_used << " unmapped channels:</em></b>";
      for ( int c=0; c<nc; c++)
        {
          const bool used = cs0.used.find( Helper::toupper( chs[c] )  ) != cs0.used.end();
	  if ( ! used ) std::cout << " " << chs[c] ;
	}
      std::cout << "</li>";
    }

  std::cout << "</ul>";
    
  std::cout << "</td><td width=50%>";
  std::cout << "<ul>";
  
  std::cout << "<li>" << a_mapped << " of " << na << " annotations mapped, aliasing " << a_aliased;

  if ( a_mapped < na )
    {
      std::cout << "<li><b><em> " << na - a_mapped << " unmapped annotations:</em></b>";
      aa = anns.begin();
      while ( aa != anns.end() )
	{
	  if (  aa->second == "" ) std::cout << " " << aa->first ;
	  ++aa;
	}
      std::cout << "</li>";
    }

  std::cout << "</ul>";

  std::cout << "</td></tr></table>";
  

  //
  // First level channel aliases
  //

  if ( nc > 0 )
    {
      std::cout << "<hr><h4>Channels</h4>";
      std::cout << "<table border=0 width=100%><tr valign=\"top\"><td width=30%>";
     
      std::cout << "<em>Channel aliases</em><br>"
		<< "<table width=100%>"
		<< "<tr><th style=\"border: 1px solid\" > &nbsp; Original &nbsp; </th>"
		<< "<th style=\"border: 1px solid\" > &nbsp; Mapped &nbsp; </th>"
		<< "<th style=\"border: 1px solid\" > &nbsp; Used? &nbsp; </th>"
		<< "</tr>";

      //
      // list of original channels and any aliases
      //
      
      for ( int c=0; c<nc; c++)
	{

	  const bool used = cs0.used.find( Helper::toupper( chs[c] )  ) != cs0.used.end();
	  const bool mapped = nonuniq_chs[c] != chs[c];
	    
	  std::cout << "<tr>";
	  std::cout << "<td style=\"text-align: center; background: " << (!mapped ? "#eeeeee" : "#ffffff" ) << "\">" << nonuniq_chs[c] << "</td>"
		    << "<td style=\"text-align: center; background: #eeeeee; "  << (used ? "font-weight: bold" : "color: orange" ) << "\">" << chs[c] << "</td>"
		    << "<td style=\"text-align: center; background: #eeeeee\">" << (used ? "Y" : "N" ) << "</td>";
	  
	  std::cout << "</tr>";
	    
	}
      
      std::cout << "</table>";
      
      std::cout << "</td> <td>&nbsp;</td> <td width=30%>";
      
      //
      // Harmonized EDF channels
      //
      
      std::cout << "<em>Harmonized EDF</em><br>"
		<< "<table width=100%>"
                << "<tr><th style=\"border: 1px solid\" > &nbsp; Harmonized &nbsp; </th>"
                << "<th style=\"border: 1px solid\" > &nbsp; Pri. &nbsp; </th>"
                << "<th style=\"border: 1px solid\" > &nbsp; Ref. &nbsp; </th>"
		<< "</tr>";
      
      ss = cs0.okay.begin();
      while ( ss != cs0.okay.end() )
	{
	  if ( ss->second )
	    {
	      std::cout << "<tr>"
			<< "<td style=\"text-align: center; background: #eeeeee\"><b>" << ss->first << "</b></td>"
			<< "<td style=\"text-align: center; background: #eeeeee\">" << cs0.sig[ ss->first ] << "</td>"
			<< "<td style=\"text-align: center; background: #eeeeee\">" << cs0.ref[ ss->first ] << "</td>"
			<< "</tr>";
	    }
	  
	  ++ss;
	}

      std::cout << "</table>";
      std::cout << "<td> &nbsp; </td>";
      std::cout << "</td><td width=30%>";
	    
      //
      // Base EDF channels
      //
      
      std::cout << "<em>Base EDF</em><br>"
		<< "<table width=100%>"
                << "<tr><th style=\"border: 1px solid\" > &nbsp; Base &nbsp; </th>"
                << "<th style=\"border: 1px solid\" > &nbsp; Pri. &nbsp; </th>"
                << "<th style=\"border: 1px solid\" > &nbsp; Ref. &nbsp; </th>"
		<< "</tr>";
      
      ss = cs1.okay.begin();
      while ( ss != cs1.okay.end() )
        {
	  
          if ( ss->second )
            {
		
	      std::cout << "<tr>"
                        << "<td style=\"text-align: center; background: #eeeeee\"><b>" << ss->first << "</b></td>"
                        << "<td style=\"text-align: center; background: #eeeeee\">" << cs1.sig[ ss->first ] << "</td>"
                        << "<td style=\"text-align: center; background: #eeeeee\">" << cs1.ref[ ss->first ] << "</td>"
                        << "</tr>";
            }
	  
	  ++ss;
	}
      
      std::cout << "</table>";

      std::cout << "</td></tr>";
      std::cout << "</table>";
      
    }
        

  //
  // Anntotation mappings
  //

  if ( anns.size() > 0 )
    {

      std::cout << "<hr><h4>Annotations</h4>";
      
      std::cout << "<table border=0 width=100%><tr valign=\"top\"><td width=50%>";

      std::cout << "<em>Mapped annotations</em><br>"
		<< "<table width=100%>"
		<< "<tr><th style=\"border: 1px solid\" > &nbsp; Original &nbsp; </th>"
		<< "<th style=\"border: 1px solid\" > &nbsp; Mapped  &nbsp; </th></tr>";
      
      aa = anns.begin();
      while ( aa != anns.end() )
	{
	  if ( aa->second != "" )
	    {
	      std::cout << "<tr>";
	      
	      bool mapped = aa->second !=  aa->first;
	      
	      if ( mapped )
		std::cout << "</td><td style=\"text-align: center; background: #ffffff\">";
	      else
		std::cout << "</td><td style=\"text-align: center; background: #eeeeee\">";

	      std::cout << aa->first 
			<< "</td>";

	      std::cout << "<td style=\"text-align: center; background: #eeeeee\">"
			<< aa->second ;

	      std::cout << "</tr>";
	      
	    }
	  ++aa;
	}

      std::cout << "</table>";
      
      std::cout << "</td> <td>&nbsp;</td> <td>";

      std::cout << "<em>Unmapped annotations</em><br>"
		<< "<table width=100%>"
		<< "<tr><th style=\"border: 1px solid\" > &nbsp; Label &nbsp; </th></tr>";

      aa = anns.begin();
      while ( aa != anns.end() )
        {
          if ( aa->second == "" )
            std::cout << "<tr><td style=\"text-align: center; background: #eeeeee; color: orange\">"
		      << aa->first
	              << "</td></tr>";
          ++aa;
	}

      std::cout << "</table>";
      std::cout << "</td></tr>";
      std::cout << "</table>";
      
    }
  


  //
  // All done
  //
  
}
