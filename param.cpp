
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


#include "eval.h"

#include "luna.h"

extern logger_t logger;

extern writer_t writer;

extern freezer_t freezer;


//
// param_t 
//

void param_t::add( const std::string & option , const std::string & value ) 
{
  // set key=value pairs to opt[]

  if ( option == "" ) return;
  
  // we're assuming this is on the same luna command line, which makes
  // no sense to have multiple versions

  // two special cases:
  //  1) if key+=value, then ","-append to any exists list
  //  2) in API mode, allow key to already exist

  const bool append_mode = option[ option.size() - 1 ] == '+';

  if ( append_mode )
    {
      const std::string option1 = option.substr( 0 , option.size() - 1 );
      if ( option1 == "" ) return;
      if ( opt.find( option1 ) == opt.end() )
	opt[ option1 ] = value;
      else
	opt[ option1 ] = opt[ option1 ] + "," + value;
      return;
    }

  // else check no doubles unless in API mode
  if ( ! globals::api_mode ) 
    if ( opt.find( option ) != opt.end() ) 
      Helper::halt( option + " parameter specified twice, only one value would be retained" );

  opt[ option ] = value; 
  
}  


void param_t::add_hidden( const std::string & option , const std::string & value ) 
{
  add( option , value );
  hidden.insert( option );
}


int param_t::size() const 
{ 
  // handle hidden things...
  return opt.size() - hidden.size();
}


void param_t::parse( const std::string & s )
{
  std::vector<std::string> tok = Helper::quoted_parse( s , "=" );
  if ( tok.size() == 2 )     add( tok[0] , tok[1] );
  else if ( tok.size() == 1 ) add( tok[0] , "__null__" );
  else // ignore subsequent '=' signs in 'value'  (i.e. key=value=2  is 'okay', means "value=2" is set to 'key')
    {
      std::string v = tok[1];
      for (int i=2;i<tok.size();i++) v += "=" + tok[i];
      add( tok[0] , v );
    }
}


void param_t::update( const std::string & id , const std::string & wc )
{

  // replace all instances of 'globals::indiv_wildcard' with 'id'
  // for all values

  // TODO: amend this so that we can edit keys (or generic text)
  // with variables and includes;  currently, only the values of keys
  //  e.g.  for x=y, only y can be a variable/include  

  std::map<std::string,std::string>::iterator ii = opt.begin();
  while ( ii != opt.end() ) 
    {
      
      std::string v = ii->second;
      bool changed = false;

      // 1. replace indiv wildcard (e.g. ^) with this person's ID
      //    nb.  this will also happen via ${id} which is a special, automatic individual-level
      //         variable
      
      while ( v.find( wc ) != std::string::npos )
	{
	  int p = v.find( wc );
	  v = v.substr( 0 , p ) + id + v.substr(p+1);
	  changed = true;
	}
      
      // 2. for any @{includes}, insert contents of file (comma-delimited)
      
      if ( Helper::swap_in_includes( &v ) )
	changed = true;

      //
      // needs updating?
      //

      if ( changed || v != ii->second ) 
	ii->second = v;
      
      ++ii;
    }
  
}

void param_t::clear() 
{ 
  opt.clear(); 
  hidden.clear(); 
} 

bool param_t::has(const std::string & s ) const 
{
  return opt.find(s) != opt.end(); 
} 

bool param_t::empty(const std::string & s ) const
{
  if ( ! has( s ) ) return true; // no key
  return opt.find( s )->second == "__null__";
}

bool param_t::yesno(const std::string & s ) const
{
  if ( ! has( s ) ) return false;
  return Helper::yesno( opt.find( s )->second ) ; 
}

std::string param_t::value( const std::string & s , const bool uppercase ) const 
{ 
  if ( has( s ) )
    return uppercase ?
      Helper::remove_all_quotes( Helper::toupper( opt.find( s )->second ) )
      : Helper::remove_all_quotes( opt.find( s )->second );
  else
    return "";
}

bool param_t::single() const 
{ 
  return size() == 1; 
}

std::string param_t::single_value() const 
{ 
  if ( ! single() ) Helper::halt( "no single value" ); 
  
  std::map<std::string,std::string>::const_iterator ii = opt.begin();
      
  while ( ii != opt.end() ) 
    {
      if ( hidden.find( ii->first ) == hidden.end() )
	return Helper::remove_all_quotes( ii->first );
      ++ii;
    }
  return ""; // should not happen
}

std::string param_t::single_pair( std::string * value ) const 
{ 
  if ( ! single() ) Helper::halt( "no single value/pair" ); 
  
  std::map<std::string,std::string>::const_iterator ii = opt.begin();
  
  while ( ii != opt.end() ) 
    {
      if ( hidden.find( ii->first ) == hidden.end() )
	{
	  *value = Helper::remove_all_quotes( ii->second );
	  return Helper::remove_all_quotes( ii->first );
	}
      ++ii;
    }
  // should not happen
  *value = "";
  return ""; 
}

std::string param_t::requires( const std::string & s , const bool uppercase ) const
{
  if ( ! has(s) ) Helper::halt( "command requires parameter " + s );
  return value(s, uppercase );
}

int param_t::requires_int( const std::string & s ) const
{
  if ( ! has(s) ) Helper::halt( "command requires parameter " + s );
  int r;
  if ( ! Helper::str2int( value(s) , &r ) ) 
    Helper::halt( "command requires parameter " + s + " to have an integer value" );
  return r;
}

double param_t::requires_dbl( const std::string & s ) const
{
  if ( ! has(s) ) Helper::halt( "command requires parameter " + s );
  double r;
  if ( ! Helper::str2dbl( value(s) , &r ) ) 
    Helper::halt( "command requires parameter " + s + " to have a numeric value" );
  return r;
}

std::string param_t::dump( const std::string & indent , const std::string & delim ) const
{
  std::map<std::string,std::string>::const_iterator ii = opt.begin();
  int sz = opt.size();
  int cnt = 1;
  std::stringstream ss;
  while ( ii != opt.end() ) 
    {

      if ( ii->second != "__null__" )
	ss << indent << ii->first << "=" << ii->second; 
      else
	ss << indent << ii->first ;

      if ( cnt != sz )
	ss << delim; 
      
      ++cnt;
      ++ii;
    }
  return ss.str();
}

std::set<std::string> param_t::strset( const std::string & k , const std::string delim , const bool uppercase ) const
{
  std::set<std::string> s;
  if ( ! has(k) ) return s;
  std::vector<std::string> tok = Helper::quoted_parse( value(k , uppercase ) , delim );
  for (int i=0;i<tok.size();i++)
    s.insert( Helper::unquote( tok[i]) );
  return s;
}

std::set<std::string> param_t::strset_xsigs( const std::string & k , const std::string delim , const bool uppercase ) const
{
  std::set<std::string> s;
  if ( ! has(k) ) return s;
  const std::string t = Helper::incexc( Helper::xsigs( value(k,uppercase) ) );
  std::vector<std::string> tok = Helper::quoted_parse( t , delim );
  for (int i=0;i<tok.size();i++)
    s.insert( Helper::unquote( tok[i]) );
  return s;
}

std::vector<std::string> param_t::strvector( const std::string & k , const std::string delim , const bool uppercase ) const
{
  std::vector<std::string> s;
  if ( ! has(k) ) return s;
  std::vector<std::string> tok = Helper::quoted_parse( value(k,uppercase) , delim );
  for (int i=0;i<tok.size();i++)
    s.push_back( Helper::unquote( tok[i]) );
  return s;
}

// embed [x][y] expansion in strvector, as well as [x] & [-x] inc/exc
std::vector<std::string> param_t::strvector_xsigs( const std::string & k , const std::string delim , const bool uppercase ) const
{  
  std::vector<std::string> s;
  if ( ! has(k) ) return s;

  // first get string and process for xsigs, then tokenize
  const std::string t = Helper::incexc( Helper::xsigs( value(k,uppercase) ) );  
  //  std::cout << "\n\nt [" << t << "]\n";
  std::vector<std::string> tok = Helper::quoted_parse( t , delim );
  for (int i=0;i<tok.size();i++)
    s.push_back( Helper::unquote( tok[i]) );
  return s;
}

std::vector<double> param_t::dblvector( const std::string & k , const std::string delim ) const
{
  std::vector<double> s;
  if ( ! has(k) ) return s;
  std::vector<std::string> tok = Helper::quoted_parse( value(k) , delim );
  for (int i=0;i<tok.size();i++) 
    {
      std::string str = Helper::unquote( tok[i]);
      double d = 0;
      if ( ! Helper::str2dbl( str , &d ) ) Helper::halt( "Option " + k + " requires a double value(s)" );
      s.push_back(d); 
    }
  return s;
}

std::vector<int> param_t::intvector( const std::string & k , const std::string delim ) const
{
  std::vector<int> s;
  if ( ! has(k) ) return s;
  std::vector<std::string> tok = Helper::quoted_parse( value(k) , delim );
  for (int i=0;i<tok.size();i++) 
    {
      std::string str = Helper::unquote( tok[i]);
      int d = 0;
      if ( ! Helper::str2int( str , &d ) ) Helper::halt( "Option " + k + " requires an integer value(s)" );
      s.push_back(d);
    }
  return s;
}


std::set<std::string> param_t::keys() const
{
  std::set<std::string> s;
  std::map<std::string,std::string>::const_iterator ii = opt.begin();
  while ( ii != opt.end() )
    {
      s.insert( ii->first );
      ++ii;
    }
  return s;
}

