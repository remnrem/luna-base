
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


#ifndef __NSRR_REMAP_H__
#define __NSRR_REMAP_H__

#include <map>
#include <set>
#include <string>
#include <vector>

struct nsrr_t { 

  // set up stage mappings 
  static void init();
  
  // set up all NSRR annot mappings
  static void init_nsrr_mappings();
  
  // indicates that no remapping is to be done)
  //static bool do_remap;
    
  // do mapping
  static std::string remap( const std::string & ); 
  
  // annotation map data (one-to-one)
  static std::map<std::string,std::string> amap; // alias --> orig

  // UC-primary->preferred primary case mapping
  static std::map<std::string,std::string> pmap; // ALIAS -> Alias

  // from primary --> multiple (can be one to many) 
  static std::map<std::string,std::vector<std::string > > bmap; // orig --> alias(es)

  // add a new annotation remap (in 'alias/remap' format  canonical|alias1|"alias2 |2"
  static void annot_remapping( const std::string & s );
  
  // add a new annotation remap
  static void add( const std::string & primary , const std::string & alias );

  // clear all existing 
  static void clear();

  // only return annots that are white-listed
  static bool whitelist;

  // only return annots that are non white-listed
  static bool unmapped;

  // for EDF+ annotations only: make as class instead of edt_annot_t --> inst.
  static void edf_annot_class( const std::string & s );
  static std::set<std::string> edf_class;
  static bool all_edf_class;  // make /all/ EDF annotations class level (i.e. no edf_annot )
  static bool as_edf_class( const std::string & s );
  
};

#endif
