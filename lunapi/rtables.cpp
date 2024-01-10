
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

std::string rtable_t::dump() const
{

  if ( nrows == -1 ) return "<empty>";
  
  std::stringstream ss;

  const int ncols = cols.size();
  
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

	  const rtable_elem_t & e = data[ j ][ i ];

	  if ( std::holds_alternative<double>( e ) )
	    ss << std::get<double>( e ) ;
	  else if ( std::holds_alternative<double>( e ) )
	    ss << std::get<int>( e ) ;
	  else if ( std::holds_alternative<std::string>( e ) )
	    ss << std::get<std::string>( e ) ;
	  else
	    ss << ".";
	  
	}
      ss << "\n";
    } 
  
  return ss.str();
  
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
  std::vector<rtable_elem_t> d( nrows , std::monostate{} );
  for (int i=0;i<nrows;i++)
    if ( ! m[i] ) d[i] = x[i] ;
  data.push_back( d );
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
  std::vector<rtable_elem_t> d( nrows , std::monostate{} );
  for (int i=0;i<nrows;i++)
    if ( ! m[i] ) d[i] = x[i] ;
  data.push_back( d );
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
  std::vector<rtable_elem_t> d( nrows , std::monostate{} );
  for (int i=0;i<nrows;i++)
    if ( ! m[i] ) d[i] = x[i] ;
  data.push_back( d );
}


//
// rtables
//

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

std::vector<std::pair<std::string,std::string> > rtables_t::list() const
{
  std::vector<std::pair<std::string,std::string> > r;
  std::map<std::string,std::map<std::string,rtable_t> >::const_iterator tt = tables.begin();
  while ( tt != tables.end() )
    {
      std::map<std::string,rtable_t>::const_iterator ss = tt->second.begin();
      while ( ss != tt->second.end() )
	{
	  r.push_back( std::make_pair( tt->first , ss->first ) );
	  ++ss;
	}
      ++tt;
    }
  return r;
}

  
rtable_t rtables_t::table( const std::string & cmd , const std::string & strata ) const
{
  std::map<std::string,std::map<std::string,rtable_t> >::const_iterator tt = tables.find( cmd );
  if ( tt == tables.end() ) return rtable_t();
  std::map<std::string,rtable_t>::const_iterator ss = tt->second.find( strata );
  if ( ss == tt->second.end() ) return rtable_t();
  return ss->second;    
}

rtable_return_t rtables_t::data( const std::string & cmd , const std::string & strata ) const
{
  std::map<std::string,std::map<std::string,rtable_t> >::const_iterator tt = tables.find( cmd );
  if ( tt == tables.end() ) return rtable_return_t();
  std::map<std::string,rtable_t>::const_iterator ss = tt->second.find( strata );
  if ( ss == tt->second.end() ) return rtable_return_t();
  return std::make_tuple( table( cmd , strata ).cols , ss->second.data );    
}


rtables_return_t rtables_t::data() const
{

  std::map<std::string,std::map<std::string,rtable_return_t> > r;

  std::map<std::string,std::map<std::string,rtable_t> >::const_iterator tt = tables.begin();
  while ( tt != tables.end() )
    {
      std::map<std::string,rtable_t>::const_iterator ss = tt->second.begin();
      while ( ss != tt->second.end() )
	{
	  r[ tt->first ][ ss->first ] = std::make_tuple( table( tt->first , ss->first ).cols , ss->second.data );
	  ++ss;
	}
      ++tt;
    }  
  return r;
}


void rtables_t::dump() const
{
  std::map<std::string,std::map<std::string,rtable_t> >::const_iterator tt = tables.begin();
  while ( tt != tables.end() )
    {
      std::map<std::string,rtable_t>::const_iterator ss = tt->second.begin();
      while ( ss != tt->second.end() )
	{
	  std::cout << tt->first << "\t" << ss->first << "\n"
		    << ss->second.dump() << "\n"
		    << std::string( 80, '-' ) << "\n";
	  ++ss;
	}
      ++tt;
    }
}
