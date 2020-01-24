
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

#include "helper/zfile.h"

extern globals global;

void zfile_t::display() const { 
  
  std::set<std::string>::const_iterator ii = vars.begin();
  while ( ii != vars.end() )
    {
      std::cout << " v = " << *ii << "\n";
      ++ii;
    }
  
  ii = facs.begin();
  while ( ii != facs.end() )
    {
      std::cout << " f = " << *ii << "\n";
      ++ii;
    }

}




bool zfile_t::set_stratum( const std::string & f , const std::string & l )
{
    
  // this is called if a 'level()' was called by writer_t
  
  // we are changing strata, so write out buffer if it has been written to 
  // if empty, nothing will happen (i.e. as will be the case when first setting 
  // a level

  write_buffer();
  
  // 'facs' should only contain normal factors, and tags; i.e. no
  // commands, but

  if ( facs.find( f ) == facs.end() ) Helper::halt( "factor " + f + " not specified" );
 
  stratum[ f ] = l;
  
  // other levels are kept as is, which should be the right thing to do

  return true;
}


bool zfile_t::set_stratum( const std::map<std::string,std::string> & fl )
{
    
  // this is called if a 'level()' was called by writer_t
  
  // we are changing strata, so write out buffer if it has been written to 
  // if empty, nothing will happen (i.e. as will be the case when first setting 
  // a level
  //std::cout << "in set al stratum\n";

  write_buffer();
  
  // 'facs' should only contain normal factors, and tags; i.e. no
  // commands, but
  
  std::map<std::string,std::string>::const_iterator ii = fl.begin();
  while ( ii != fl.end() )
    {

      if ( facs.find( ii->first ) == facs.end() ) 
	Helper::halt( "factor " + ii->first + " not specified" );
      ++ii;
    }
  
  stratum = fl;
  
  return true;
}


bool zfile_t::set_value( const std::string & k , const std::string & v )
{
  buf[ k ] = v;
  return true;
}

bool zfile_t::set_value( const std::string & k , int v )
{
  buf[ k ] = Helper::int2str( v );
  return true;
}

bool zfile_t::set_value( const std::string & k , double v )
{
  buf[ k ] = Helper::dbl2str( v );
  return true;
}

void zfile_t::write_header() 
{
  
  bool first_col = true;
  if ( parent->show_indiv_col ) { print( "ID" ); first_col = false; } 

  std::set<std::string>::const_iterator ff = facs.begin();
  while ( ff != facs.end() )
    {
      if ( ! first_col ) print( "\t" );      
      print( *ff );
      first_col = false;
      ++ff;
    }

  std::set<std::string>::const_iterator vv = vars.begin();
  while ( vv != vars.end() )
    {
      if ( ! first_col ) print( "\t" );      
      print( *vv );
      first_col = false;
      ++vv;
    }

  print( "\n" );

}

void zfile_t::write_buffer()
{
  
  //  std::cout << "in write buffer " << buf.size() << "\n";
  
  // if no variables specified, nothing to do [ i.e. do not write a blank row ] 
  // as this will be called every time a level is set, this means it will write out
  // completed buffers, but should not give an empty row

  if ( buf.size() == 0 ) return;
  
  // all factors must be specified 
  //  std::cout << "facs st = " << facs.size() << " " << stratum.size() << " " << vars.size() << "\n";

  if ( facs.size() != stratum.size() ) Helper::halt( "not all levels specified" );

  if ( parent->show_indiv_col ) print ( indiv );

  std::map<std::string,std::string>::const_iterator ss = stratum.begin();
  while ( ss != stratum.end() ) 
    {
      print( "\t" );
      print( ss->second );
      //std::cout << "st = " << ss->second << "\n";
      ++ss;
    }
  
  // okay if some variables are missing.. we'll just write "NA" instead
  
  std::set<std::string>::const_iterator vv = vars.begin();
  while ( vv != vars.end() )
    {
      //std::cout << "v = " << *vv << "\n";
      print( "\t" );
      bool missing = buf.find( *vv ) == buf.end() ;
      if ( missing ) print( "NA" );
      else print( buf.find( * vv )->second );
      
      //std::cout << *vv << " = " << ( missing ? "NA" : buf.find( * vv )->second ) << "\n";
      ++vv;
    }
  
  print( "\n" );
  
  // clear buffers
  stratum.clear(); 
  buf.clear();
  
}
