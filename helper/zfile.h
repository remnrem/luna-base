
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


#ifndef __LUNA_ZFILE_H__
#define __LUNA_ZFILE_H__

#include "helper/zfstream.h"
#include <fstream>
#include "cmddefs.h"
#include "defs/defs.h"
#include "eval.h"

#include <string>
#include <map>

extern globals global;

struct zfiles_t;

struct zfile_t { 

 public:

 zfile_t( zfiles_t * p , 
	  const std::string & n , 
	  const std::string & indiv , 
	  const std::string & cmd , 
	  const std::string & table ,
	  const param_t * param = NULL , 
	  bool compressed = true ) 
 : parent(p) , indiv(indiv) , cmd(cmd) , table(table), compressed(compressed) 
  { 
    
    if ( compressed ) 
       zout.open( n.c_str() , std::ios_base::out );
     else
       out.open( n.c_str() );
     
     // this (via the tfac_t() conversion) purposefully ignores any
     // TAG-defined factors... i.e. get straight to core variables
     // as represented in the cmddefs_t
     
     vars = globals::cmddefs().variables( cmd , param , tfac_t( table , "_" ) );
     
     // however, in generating the output, we'll need to respect the TAG-defined
     // factors... so track these separately
     
     facs = str2set( table , "_" );
     
     // and write the header
     
     write_header();

   }

  void display() const;
  
  template <class T>
  void print ( const T & rhs ) {
    if ( compressed ) zout << rhs;
    else out << rhs; 
  }

  void print ( const std::string & rhs ) {
    if ( compressed ) zout << rhs;
    else out << rhs; 
  }

  friend zfile_t & operator<< ( zfile_t & zf , const std::string & rhs );
  
  void close()
  { 
    write_buffer();
    if ( out.is_open() ) out.close();
    if ( zout.is_open() ) zout.close();
  } 
  
 ~zfile_t() 
   { 
     close();
   } 

  
  bool set_stratum( const std::string &  , const std::string &  );
  bool set_stratum( const std::map<std::string,std::string> &  );
  bool set_value( const std::string &  , const std::string &  );
  bool set_value( const std::string & , int  );
  bool set_value( const std::string &  , double  );
  void write_header() ;
  void write_buffer();
  
private:
  
  zfiles_t * parent;
  
  gzofstream zout;
  
  std::ofstream out;
  
  std::string indiv, cmd, table;
  
  bool compressed;
  
  // factors (defined in cmddefs_t)
  std::set<std::string> facs; 
  
  // expceted variable names
  std::set<std::string> vars;
  
  //
  // data
  //

  // levels (i.e. defining current write-strata)
  std::map<std::string,std::string> stratum;
  
  // data buffer (for this strata)
  std::map<std::string,std::string> buf;
 
  
  // helper function
  std::set<std::string> str2set( const std::string & str , const std::string & delim = "," ) 
    {
      std::vector<std::string> tok = Helper::parse( str , delim );
      std::set<std::string> r;
      for (int i=0; i<tok.size(); i++) r.insert( tok[i] );
      return r;
    }

};


// need track of 'N' output files, but here 
// only for one individual at a time

struct zfiles_t { 
  
  // root/indiv/{prepend}command-table{append}.txt.gz
  // root/indiv/{prepend}command-table{append}.txt
  
  zfiles_t( const std::string & fileroot , 
	    const std::string & _indiv ) 
  { 
    
    indiv = _indiv;
    
    // create this folder if it does not exist    
    folder = fileroot + globals::folder_delimiter + indiv + globals::folder_delimiter ; 
    
    // create folder if it does not exist (will need to change for windows...)
    std::string syscmd = "mkdir -p " + folder ; 
    
    int dummy = system( syscmd.c_str() );
    
    //std::cerr << "Writing to " << folder << "\n";
    
    // by default, show ID col and header row
    mode(true,true);


  } 

  void mode( const bool indiv_col = true , const bool header_row = true )
  {
    show_indiv_col = indiv_col;
    show_header_row = header_row;
  }

  // check whether file exists 
  zfile_t * exists( const std::string & cmd , const std::string & table ) 
  {
    std::map<std::string,std::map<std::string,zfile_t*> >::const_iterator cc = files.begin();
    if ( cc == files.end() ) return NULL;
    const std::map<std::string,zfile_t*> & t = cc->second;
    std::map<std::string,zfile_t*>::const_iterator tt = t.find( table );
    if ( tt == t.end() ) return NULL;
    return tt->second;
  }
  
  // access file by command-name/table ID
  zfile_t * file( const std::string & cmd , const param_t * param , const std::string & table  ) 
  {
    
    std::map<std::string,std::map<std::string,zfile_t*> >::const_iterator cc = files.find( cmd );
    
    if ( cc == files.end() )
	return new_file( cmd , param , table );

    const std::map<std::string,zfile_t*> & t = cc->second;

    std::map<std::string,zfile_t*>::const_iterator tt = t.find( table );

    if ( tt == t.end() ) return new_file( cmd , param , table );

    return tt->second;
  }

  void close()
  {
    std::map<std::string,std::map<std::string,zfile_t*> >::iterator cc = files.begin();
    while ( cc != files.end() ) {
      std::map<std::string,zfile_t*> & t = cc->second;
      std::map<std::string,zfile_t*>::iterator tt = t.begin();
      while ( tt != t.end() ) { 
	if ( tt->second != NULL ) 
	  {
	    tt->second->close();
	    delete tt->second;
	    tt->second = NULL;
	  }
	++tt;
      }
      ++cc;
    }
    files.clear();
  }

  
  // ensure all files are closed
  ~zfiles_t()
  {
    close();
  }

  bool show_header_row;

  bool show_indiv_col;

  std::string folder;
  
  std::string indiv;

private:

  // cmd -> table -> file
  std::map<std::string,std::map<std::string,zfile_t*> > files;
  
  
  zfile_t * new_file( const std::string & cmd , 		      
		      const param_t * param , 
		      const std::string & table ) 
  {
    
    // figure out whether this table should be compressed, according to cmddefs

    tfac_t tfac( table , "_" );
    
    // is this a valid table?  
	 
    if ( ! globals::cmddefs().exists( cmd , tfac ) )
      return NULL;
      	 
    // should this be compressed by default?
	 bool compressed = globals::cmddefs().out_compressed( cmd , tfac );

    // create a new file, store a pointer to it, and return that pointer
    //  prepend-COMMAND-F1_F2_F3{_append}.txt{.gz}
 
    // COMM.txt
    // COMM_F1_F2.txt
    // COMM_F1-V1.txt  // sets F1 to "V1"
    // {prepend}COMM{_FAC_FAC}{append}.txt

    std::string filename = folder 
      + globals::txt_table_prepend 
      + cmd 
      + ( table == "" ? "" : "_" + table ) 
      + globals::txt_table_append 
      + ( compressed ? ".txt.gz" : ".txt" ) ;

    zfile_t * p = new zfile_t( this , filename , indiv , cmd , table , param , compressed );
    files[ cmd ][ table ] = p;
    return p;
  }

  
}; 



#endif
