
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


#include "lunapi/rtables.h"

rtable_t::rtable_t()
{
  nrows = -1;
}

std::string rtable_t::dump()
{

  if ( nrows == -1 ) return "<empty>";
  
  std::stringstream ss;

  if ( ! check() ) return "";
  
  const int ncols = cols.size();

  // 012 = str dbl int
  std::map<std::string,int> ctype;
  for (int j=0; j<ncols; j++)
    {
      if ( strcols.find( cols[j] ) != strcols.end() ) ctype[ cols[j] ] = 0;
      else if ( dblcols.find( cols[j] ) != dblcols.end() ) ctype[ cols[j] ] = 1;
      else ctype[ cols[j] ] = 2;
    }
  
  // header
  for (int j=0; j<ncols; j++)
    {
      if ( j ) ss << "\t";
      ss << cols[j];
    }
  ss << "\n";
  
  // data
  for (int i=0;i<nrows;i++)
    {
      for (int j=0; j<ncols; j++)
	{
	  if ( j ) ss << "\t";

	  const std::string & col = cols[j];
	  
	  int ct = ctype[ col ];

	  if ( ct == 0 )
	    {
	      if ( strmiss[ col ][ i ] ) ss << ".";
	      else ss << strcols[ col ][ i ] ;
	    }
	  else if ( ct == 1 )
	    {
	      if ( dblmiss[ col ][ i ] ) ss << "NA";
              else ss << dblcols[ col ][ i ] ;	      
	    }
	  else
	    {
	      if ( intmiss[ col ][ i ] ) ss << "NA";
              else ss << intcols[ col ][ i ] ;	      
	    }
	}
      ss << "\n";
    }  
  
  return ss.str();
  
}



bool rtable_t::check() const
{
  std::set<std::string> ucols;
  for (int i=0;i<cols.size();i++) ucols.insert( cols[i] );
  return ucols.size() == cols.size();
}


void rtable_t::checkrows( int n )
{
  if ( nrows == -1 )
    nrows = n;
  else if ( nrows != n )
    Helper::halt( "internal problem building an rtable_t" );    
}

void rtable_t::add( const std::string & v , const std::vector<std::string> & x )
{
  checkrows( x.size() ); 
  std::vector<bool> missing( nrows , false );
  add( v, x , missing );
}

void rtable_t::add( const std::string & v , const std::vector<std::string> & x , const std::vector<bool> & m )
{
  cols.push_back(v);
  checkrows( x.size() );
  checkrows( m.size() );
  strcols[ v ] = x;
  strmiss[ v ] = m;
}

// doubles

void rtable_t::add( const std::string & v , const std::vector<double> & x )
{
  checkrows( x.size() ); 
  std::vector<bool> missing( nrows , false );
  add( v, x , missing );
}

void rtable_t::add( const std::string & v , const std::vector<double> & x , const std::vector<bool> & m )
{
  cols.push_back(v);
  checkrows( x.size() );
  checkrows( m.size() );
  dblcols[ v ] = x;
  dblmiss[ v ] = m;
}

// doubles

void rtable_t::add( const std::string & v , const std::vector<int> & x )
{
  checkrows( x.size() ); 
  std::vector<bool> missing( nrows , false );
  add( v, x , missing );
}

void rtable_t::add( const std::string & v , const std::vector<int> & x , const std::vector<bool> & m )
{
  cols.push_back(v);
  checkrows( x.size() );
  checkrows( m.size() );
  intcols[ v ] = x;
  intmiss[ v ] = m;
}

rtables_t::rtables_t( const retval_t & retval )
{
  tables = retval.make_tables();
}

std::vector<std::string> rtables_t::commands() const
{
  std::vector<std::string> r;
  std::map<std::string,std::map<std::string,rtable_t> >::const_iterator tt = tables.begin();
  while ( tt != tables.end() )
    {
      r.push_back( tt->first );
      ++tt;
    }
  return r;
}

std::vector<std::vector<std::string> > rtables_t::list() const
{
  std::vector<std::vector<std::string> > r(2);
  std::map<std::string,std::map<std::string,rtable_t> >::const_iterator tt = tables.begin();
  while ( tt != tables.end() )
    {
      std::map<std::string,rtable_t>::const_iterator ss = tt->second.begin();
      while ( ss != tt->second.end() )
	{
	  r[0].push_back( tt->first );
	  r[1].push_back( ss->first );
	  ++ss;
	}
      ++tt;
    }
  return r;
}

  
rtable_t rtables_t::table( const std::string & cmd , const std::string & strata )
{
  std::map<std::string,std::map<std::string,rtable_t> >::const_iterator tt = tables.find( cmd );
  if ( tt == tables.end() ) return rtable_t();
  std::map<std::string,rtable_t>::const_iterator ss = tt->second.find( strata );
  if ( ss == tt->second.end() ) return rtable_t();
  return ss->second;    
}

