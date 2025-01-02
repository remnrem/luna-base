
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


#include "token.h"
#include "token-eval.h"

#include "helper.h"
#include "miscmath/crandom.h"

#include "annot/annot.h"

#include <sstream>
#include <cmath>
#include <algorithm>

#include <iostream>


// TODO
// for now, vector x scalar comparisons do not allow for type conversion


void TokenFunctions::attach( instance_t * m )
{ 
  meta = m; 
  accumulator = NULL;
  global_vars = NULL;
}

void TokenFunctions::attach( instance_t * m , instance_t * a , const std::set<std::string> * gv )
{ 
  meta = m;   
  accumulator = a;
  global_vars = gv;
}

void Token::set()
{
    ttype = UNDEF;
}

void Token::set( const std::string & s )
{
    ttype = STRING;
    sval = s;
}

void Token::set( const double d )
{
    ttype = FLOAT;
    fval = d;
}

void Token::set( const int i )
{
    ttype = INT;
    ival = i;
}

void Token::set( const bool b )
{
    ttype = BOOL;
    bval = b;
}

// vectors

void Token::set( const std::vector<std::string> & s )
{
  if ( s.size() == 1 ) 
    set( s[0] );
  else
    {
      ttype = STRING_VECTOR;
      svec = s;
      unmask();
    }
}

void Token::set( const std::vector<double> & d )
{
  if ( d.size() == 1 ) 
    set( d[0] );
  else 
    {
      ttype = FLOAT_VECTOR;
      fvec = d;
      unmask();
    }
}

void Token::set( const std::vector<int> & i )
{
  if ( i.size() == 1 ) 
    set( i[0] ); 
  else 
    {
      ttype = INT_VECTOR;
      ivec = i;
      unmask();
    }
}

void Token::set( const std::vector<bool> & b )
{
  if (b.size()==1) 
    set(b[0]);
  else 
    { 
      ttype = BOOL_VECTOR;
      bvec = b;
      unmask();
    }
}

// vectors + masks

void Token::set( const std::vector<std::string> & s , const std::vector<int> & m )
{
  if ( s.size() == 1 ) 
    set( s[0] );
  else
    {
      ttype = STRING_VECTOR;
      svec = s;
      subset( m );
    }
}

void Token::set( const std::vector<double> & d , const std::vector<int> & m )
{
  if ( d.size() == 1 ) 
    set( d[0] );
  else 
    {
      ttype = FLOAT_VECTOR;
      fvec = d;
      subset( m );
    }
}

void Token::set( const std::vector<int> & i , const std::vector<int> & m )
{
  if ( i.size() == 1 ) 
    set( i[0] ); 
  else 
    {
      ttype = INT_VECTOR;
      ivec = i;
      subset( m );
    }
}

void Token::set( const std::vector<bool> & b , const std::vector<int> & m )
{
  if (b.size()==1) 
    set(b[0]);
  else 
    { 
      ttype = BOOL_VECTOR;
      bvec = b;
      subset( m );
    }
}

// updates

void Token::update( const std::vector<std::string> & s )
{
  if ( ttype != STRING_VECTOR ) Helper::halt( "type conflict" );
  if ( ve.size() != s.size() ) Helper::halt( "size conflict in vector subset update" );
  // update
  for (int j=0; j<ve.size(); j++)
    svec[ ve[j] ] = s[j];
  // after update, unmask
  unmask();
}

void Token::update( const std::vector<double> & d )
{
  if ( ttype != FLOAT_VECTOR ) Helper::halt( "type conflict" );
  if ( ve.size() != d.size() ) Helper::halt( "size conflict in vector subset update" );
  // update
  for (int j=0; j<ve.size(); j++)
    fvec[ ve[j] ] = d[j];
  // after update, unmask
  unmask();
}

void Token::update( const std::vector<int> & i )
{

  if ( ve.size() != i.size() )
    Helper::halt( "size conflict in vector subset update" );
  
  if ( ttype == INT_VECTOR )
    {
      for (int j=0; j<ve.size(); j++)
	ivec[ ve[j] ] = i[j];
    }
  else if (  ttype == FLOAT_VECTOR )
    {
      for (int j=0; j<ve.size(); j++)
	{
	  //	  std::cout << " updaing " << ve[j] << " with " << i[j] << "\n";
	  fvec[ ve[j] ] = i[j];
	}
    }
  else
    Helper::halt( "type conflict" );

  // after update, unmask
  unmask();


}

void Token::update( const std::vector<bool> & b )
{
  if ( ttype != BOOL_VECTOR ) Helper::halt( "type conflict" );
  if ( ve.size() != b.size() ) Helper::halt( "size conflict in vector subset update" );
  // update
  for (int j=0; j<ve.size(); j++)
    bvec[ ve[j] ] = b[j];
  // after update, unmask
  unmask();
}


// functions, operators

void Token::function( const std::string & fn )
{
  ttype = FUNCTION;
  tname = fn;
}

void Token::oper( Token::tok_type t )
{
  ttype = t;    
}

void Token::variable( const std::string & mf )
{
  ttype = VARIABLE;
  tname = mf;
}

void Token::init()
{
  
  tok_map[ "*" ]  = MULTIPLY_OPERATOR;
  //  tok_map[ "^" ]  = POWER_OPERATOR;
  tok_map[ "/" ]  = DIVIDE_OPERATOR;
  tok_map[ "%" ]  = MOD_OPERATOR;
  tok_map[ "%%" ] = MOD_OPERATOR;
  tok_map[ "+" ]  = ADD_OPERATOR;
  tok_map[ "-" ]  = SUBTRACT_OPERATOR;
  tok_map[ "&&" ] = AND_OPERATOR;
  tok_map[ "&" ]  = AND_OPERATOR;
  tok_map[ "||" ] = OR_OPERATOR;
  tok_map[ "|" ]  = OR_OPERATOR;
  tok_map[ "=" ]  = ASSIGNMENT_OPERATOR;
  tok_map[ "==" ] = EQUAL_OPERATOR;
  tok_map[ "=~" ] = HAS_OPERATOR; // like equal but vector =~ scalar returns scalar, any( lhs == rhs )
  tok_map[ "!=" ] = UNEQUAL_OPERATOR;
  tok_map[ "!" ]  = NOT_OPERATOR;
  tok_map[ "~" ]  = NOT_OPERATOR;
  tok_map[ ">" ]  = GREATER_THAN_OPERATOR;
  tok_map[ ">=" ] = GREATER_THAN_OR_EQUAL_OPERATOR;
  tok_map[ "<" ]  = LESS_THAN_OPERATOR;
  tok_map[ "<=" ] = LESS_THAN_OR_EQUAL_OPERATOR;
  
  
  //
  // Reverse mapping
  //
  
  std::map<std::string,Token::tok_type>::iterator i = tok_map.begin();
  while ( i != tok_map.end() ) 
    {
      tok_unmap[ i->second ] = i->first;
      ++i;
    }
  
  
  //
  // Token function map
  //
  
  fn_map[ "if" ]     = 1;  // number of args
  fn_map[ "ifnot" ]  = 1;  // complement of if()
  fn_map[ "sqrt" ]   = 1;  // square-root
  fn_map[ "rand" ]   = 1;  // rand(10)
  fn_map[ "rnd" ]    = 0;  // rnd() 0..1

  fn_map[ "round" ]  = 1;
  fn_map[ "floor" ]  = 1;  
  fn_map[ "abs" ]    = 1;
  fn_map[ "sqr"  ]   = 1;  // X^2
  fn_map[ "log"  ]   = 1;
  fn_map[ "log10"]   = 1;
  fn_map[ "exp"  ]   = 1;
  fn_map[ "pow"  ]   = 2;  // X^N    
  fn_map[ "ifelse" ] = 3;  // ifelse( cond , T , F )
  
  // vector functions
  
  fn_map[ "element" ] = 2;  // element(Y,i)      extract element 'i' from vector 'Y'
  fn_map[ "length" ]  = 1;  // length(Y)   length of Y
  fn_map[ "size" ]    = 1;  // size(Y)   length of Y
  fn_map[ "min" ]     = 1;  // min(Y)      minimum element in Y
  fn_map[ "max" ]     = 1;  // max(Y)      max elenebt in Y
  fn_map[ "sum" ]     = 1;  // sum(Y)      sum of elements in Y
  fn_map[ "mean" ]    = 1;  // mean(Y)     mean of elements in Y
  fn_map[ "sd" ]      = 1;  // sd(Y)       SD of elements in Y (N-1 denom)
  fn_map[ "sort" ]    = 1;  // sort(Y)     returns sorted vector Y (asc.)
  
  // vector creation 
  
  fn_map[ "num_func" ]    = -1;  // num( 1,0,1 ) -- floating point vector --> num(3,1,0,1)
  fn_map[ "int_func" ]    = -1;  // int( 1,0,1 )  ints
  fn_map[ "txt_func" ]    = -1;  // txt( '1','0','1' )  strings
  fn_map[ "bool_func" ]   = -1;  // bool( 1,0,1 )  bools
  fn_map[ "c_func" ]      = -1;  // c( expr1, expr2 , expr3 ... )         concatenate similar types   

  // misc

  fn_map[ "any" ]      = 1;   // any( expr1 )              returns BOOL , countif(x,T)>0
  fn_map[ "all" ]      = 1;   // all( expr1 )              returns BOOL , countif(x,T) == size(x)
  fn_map[ "contains" ] = 2;   // contains( expr1 , expr )  returns BOOL , countif(x,y)>0
  fn_map[ "countif" ]  = 2;   // countif( expr1, expr2 )   returns INT,  for # of elements in expr1 that match expr2


  
}



bool Token::is_bool(bool * b ) const
{ 
  if ( ttype == BOOL ) 
    {
      if ( b ) *b = bval;
      return true;
    }
  return false;
}


bool Token::is_string( std::string * s ) const 
{ 
  if ( ttype == STRING ) 
    {
      if ( s ) *s = sval;
      return true;
    }
  return false;
}


bool Token::is_float( double * f ) const
{ 
  if ( ttype == FLOAT ) 
    {
      if ( f ) *f = fval;
      return true;
    }
  return false;
}


bool Token::is_int( int * i ) const
{ 
  if ( ttype == INT ) 
    {
      if ( i ) *i = ival;
      return true;
    }
  return false;
}



bool Token::is_bool_vector( std::vector<bool> * b ) const
{ 
  if ( ttype == BOOL_VECTOR ) 
    {
      if ( b ) *b = bvec;
      return true;
    }
  return false;
}

bool Token::is_string_vector( std::vector<std::string> * s ) const 
{ 
  if ( ttype == STRING_VECTOR ) 
    {
      if ( s ) *s = svec;
      return true;
    }
  return false;
}


bool Token::is_float_vector( std::vector<double> * f ) const
{ 
  if ( ttype == FLOAT_VECTOR ) 
    {
      if ( f ) *f = fvec;
      return true;
    }
  return false;
}


bool Token::is_int_vector( std::vector<int> * i ) const
{ 
  if ( ttype == INT_VECTOR ) 
    {
      if ( i ) *i = ivec;
      return true;
    }
  return false;
}



bool Token::is_operator() const
{
  // note -- treat parenthesis separately

  return ttype == EQUAL_OPERATOR ||
    ttype == HAS_OPERATOR ||
    ttype == UNEQUAL_OPERATOR ||
    ttype == ASSIGNMENT_OPERATOR ||
    ttype == NOT_OPERATOR ||
    ttype == AND_OPERATOR ||
    ttype == OR_OPERATOR ||     
    ttype == GREATER_THAN_OPERATOR ||
    ttype == GREATER_THAN_OR_EQUAL_OPERATOR ||
    ttype == LESS_THAN_OPERATOR ||
    ttype == LESS_THAN_OR_EQUAL_OPERATOR ||
    ttype == MOD_OPERATOR ||
    ttype == MULTIPLY_OPERATOR || 
    ttype == DIVIDE_OPERATOR ||
    ttype == ADD_OPERATOR ||
    ttype == SUBTRACT_OPERATOR;
}


bool Token::is_scalar() const
{
  return ttype == INT || ttype == FLOAT || ttype == STRING || ttype == BOOL;
}

bool Token::is_vector() const
{
  return ttype == INT_VECTOR || ttype == FLOAT_VECTOR || ttype == STRING_VECTOR || ttype == BOOL_VECTOR;
}


bool Token::is_function() const
{
  return ttype == FUNCTION;
}

bool Token::is_ident() const
{
  return ! ( is_operator() || is_function() || is_left_paren() || is_right_paren() || is_separator() );
}

bool Token::is_variable() const
{
  return ttype == VARIABLE;
}

Token::Token( const std::string & s )
{
  ttype = STRING;
  sval = s;
}

Token::Token( const double d )
{
  ttype = FLOAT;
  fval = d;
}

Token::Token( const int i )
{
  ttype = INT;
  ival = i;
}

Token::Token( const bool b )
{
  ttype = BOOL;
  bval = b;
}


Token::Token( const std::vector<std::string> & s )
{
  ttype = STRING_VECTOR;
  svec = s;
  unmask();
}

Token::Token( const std::vector<double> & d )
{
  ttype = FLOAT_VECTOR;
  fvec = d;
  unmask();
}

Token::Token( const std::vector<int> & i )
{
  ttype = INT_VECTOR;
  ivec = i;
  unmask();
}

Token::Token( const std::vector<bool> & b )
{
  ttype = BOOL_VECTOR;
  bvec = b;
  unmask();
}



Token::Token( const Token & rhs )
{
    *this = rhs;
}


void Token::unmask()
{
  if ( ! is_masked() ) return;
  ve.resize( fullsize() );
  for (int i=0; i<fullsize(); i++) ve[i] = i;
}


void Token::prune()
{
  if ( ! is_vector() ) return;
  if ( ! is_masked() ) return;

  if ( ttype == BOOL_VECTOR )
    {
      std::vector<bool> bvec2;
      for (int i=0; i<ve.size(); i++) bvec2.push_back( bvec[ve[i]] );
      bvec = bvec2;
      unmask();
      return;
    }

    if ( ttype == INT_VECTOR )
    {
      std::vector<int> ivec2;
      for (int i=0; i<ve.size(); i++) ivec2.push_back( ivec[ve[i]] );
      ivec = ivec2;
      unmask();
      return;
    }

    if ( ttype == FLOAT_VECTOR )
    {
      std::vector<double> fvec2;
      for (int i=0; i<ve.size(); i++) fvec2.push_back( fvec[ve[i]] );
      fvec = fvec2;
      unmask();
      return;
    }
    
    if ( ttype == STRING_VECTOR )
    {
      std::vector<std::string> svec2;
      for (int i=0; i<ve.size(); i++) svec2.push_back( svec[ve[i]] );
      svec = svec2;
      unmask();
      return;
    }
   
}


void Token::subset( const std::vector<int> & s )
{

  // sets the ve[] vector element mask
  
  if ( ! is_vector() ) return;

  // size of current subset ( ss <= fs )
  const int ss = size();
  
  // size of full vector
  const int fs = fullsize();
  
  // we may want to relax this. e.g. if we can repeat elements
  //  but for now, just think of the mask as unique

  if ( s.size() > fs )
    Helper::halt( "subset length > full vector length" );
  
  // update index: which may be subsetting an existing ve[] index. e.g. x[c(1,2,3)][2]
  // Although.. not sure if Eval will respect that syntax yet, so for now this is a moot
  // point, i.e. subsetting will always be from an unmasked vector, but set here so that
  // it will work in the future if/when we allow X[][] indexing

  std::vector<int> ve2 = ve;

  ve.clear();

  // check all index values are valid and unique
  std::set<int> vi;

  for (int i=0; i<s.size(); i++)
    {
      if ( s[i] < 0 || s[i] >= ss )
	Helper::halt( "bad index value for vector subsetting" );
      
      vi.insert( ve2[s[i]] );      
      ve.push_back( ve2[s[i]] );
      
    }
  
  if ( vi.size() != ve.size() )
    Helper::halt( "cannot have repeated vector element index values currently" );
}



Token Token::operator!() const
{
  // handles only bools and ints
  // both scalars and vectors
  
  if ( is_bool() ) return Token( ! bval ); 
  else if ( is_int() ) return Token( ival == 0 );
  else if ( is_bool_vector() )
    {
      std::vector<bool> ans( ve.size() );
      for ( int i=0; i<ve.size(); i++) ans[i] = ! bvec[ve[i]];
      return Token( ans );
    }  
  else if ( is_int_vector() )
    {
      std::vector<bool> ans( ve.size() );
      for ( int i=0; i<ve.size(); i++) ans[i] = ! ivec[ve[i]];
      return Token( ans );
    }  
  else return Token();
}


Token & Token::operator=(const Token & rhs)
{
  Token ret;
  
  ttype = rhs.ttype;
  tname = rhs.tname;
  
  ival = rhs.ival;
  sval = rhs.sval;
  fval = rhs.fval;
  bval = rhs.bval;
  
  ivec = rhs.ivec;
  svec = rhs.svec;
  fvec = rhs.fvec;
  bvec = rhs.bvec;
  ve = rhs.ve;
  
  return *this;
}

// if a masked vector, return masked size
int Token::size() const
{
  // i.e. for 'unmasked' vector, ve.size() == storage size anyway
  if      ( is_scalar() ) return 1;
  else if ( ! is_vector() ) return 0;
  else return ve.size(); 
  return 0;
}

// full size
int Token::fullsize() const
{
  if      ( is_scalar() ) return 1;
  else if ( ! is_vector() ) return 0;
  else if ( ttype == INT_VECTOR ) return ivec.size();
  else if ( ttype == FLOAT_VECTOR ) return fvec.size();
  else if ( ttype == STRING_VECTOR ) return svec.size();
  else if ( ttype == BOOL_VECTOR ) return bvec.size();
  return 0;
}


Token Token::operator!=(const Token & rhs ) const
{

  if ( is_vector() && rhs.is_vector() ) 
    {
      
      if ( size() != rhs.size() ) return Token();
      const int sz = size();
      std::vector<bool> ans( sz );      

      if ( rhs.is_int_vector() ) 
	{	  
	  if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] != rhs.ivec[rhs.ve[i]];
	  else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] != rhs.ivec[rhs.ve[i]]; 
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] != rhs.ivec[rhs.ve[i]]; 
	  return Token( ans );
	}
      else if ( rhs.is_float_vector() )
	{
	  if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] != rhs.fvec[rhs.ve[i]];
	  else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] != rhs.fvec[rhs.ve[i]]; 
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] != rhs.fvec[rhs.ve[i]]; 
	  return Token( ans );	 
	}
      else if ( rhs.is_bool_vector() )
	{
	  if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] != rhs.bvec[rhs.ve[i]];
	  else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] != rhs.bvec[rhs.ve[i]]; 
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] != rhs.bvec[rhs.ve[i]]; 
	  return Token( ans );
	}	  
      else if ( rhs.is_string_vector() )
	{
	  if ( is_string_vector() )  for (int i=0; i<sz; i++) ans[i] = svec[ve[i]] != rhs.svec[ve[i]]; 
	  else return Token();
	  return Token( ans );
	}	  
      else
	return Token();
    }


  // vector != scalar 

  if ( is_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();

      std::vector<bool> ans( sz );      

      if ( rhs.is_int() ) 
	{
	  if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] != rhs.ival;
	  else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] != rhs.ival; 
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] != rhs.ival; 	  
	  return Token( ans );
	}
      else if ( rhs.is_float() )
	{
	  if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] != rhs.fval;
	  else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] != rhs.fval; 
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] != rhs.fval; 	  
	  return Token( ans );
	}
      else if ( rhs.is_bool() ) 
	{
	  if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] != rhs.bval;
	  else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] != rhs.bval; 
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] != rhs.bval; 	  
	  return Token( ans );
	}
      else if ( rhs.is_string() ) 
	{
	  if ( is_string_vector() ) for (int i=0; i<sz; i++) ans[i] = svec[ve[i]] != rhs.sval;  
	  else return Token(); 
	  return Token( ans );
	}
      else return Token();
    }
  

  // scalar != vector

  if ( rhs.is_vector() )
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      
      if ( is_int() )
	{
	  if      ( rhs.is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = rhs.ivec[rhs.ve[i]] != ival;
	  else if ( rhs.is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = rhs.fvec[rhs.ve[i]] != ival;
	  else if ( rhs.is_string_vector() ) return Token();
	  else if ( rhs.is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = rhs.bvec[rhs.ve[i]] != ival;
	  return Token( ans );      
	}
      else if ( is_float() )
	{
	  if      ( rhs.is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = rhs.ivec[rhs.ve[i]] != fval;
	  else if ( rhs.is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = rhs.fvec[rhs.ve[i]] != fval;
	  else if ( rhs.is_string_vector() ) return Token();
	  else if ( rhs.is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = rhs.bvec[rhs.ve[i]] != fval;
	  return Token( ans );     
	}
      else if ( is_bool() )
	{
	  if      ( rhs.is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = rhs.ivec[rhs.ve[i]] != bval;
	  else if ( rhs.is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = rhs.fvec[rhs.ve[i]] != bval;
	  else if ( rhs.is_string_vector() ) return Token();
	  else if ( rhs.is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = rhs.bvec[rhs.ve[i]] != bval;
	  return Token( ans );     
	}
      else if ( is_string() )
	{
	  if ( rhs.is_string_vector() ) for (int i=0; i<sz; i++) ans[i] = rhs.svec[rhs.ve[i]] != sval;
	  else return Token();
	  return Token( ans );     
	}
      else return Token();
    }

      
  // scalar x scalar

  if ( is_bool()   && rhs.is_bool()   )  return Token( bval != rhs.bval ); 
  if ( is_int()    && rhs.is_int()    )  return Token( ival != rhs.ival ); 
  if ( is_float()  && rhs.is_float()  )  return Token( fval != rhs.fval ); 
  if ( is_string() && rhs.is_string() )  return Token( sval != rhs.sval ); 

  // also allow int / bool and  int / float comparisons
  if ( is_int()   && rhs.is_bool()    )  return Token( ival != rhs.bval );
  if ( is_bool()  && rhs.is_int()     )  return Token( bval != rhs.ival );

  if ( is_float() && rhs.is_bool()    )  return Token( fval != rhs.bval );
  if ( is_bool()  && rhs.is_float()   )  return Token( bval != rhs.fval );

  if ( is_float() && rhs.is_int()     )  return Token( fval != rhs.ival );
  if ( is_int()   && rhs.is_float()   )  return Token( ival != rhs.fval );

  return Token();
}


Token Token::contains( const Token & rhs) const 
{

  // first test for two vectors (of any length)
  // otherwise, we can use existing operator==() code to test for scalar/vector and scalar/scalar
  
  if ( is_vector() && rhs.is_vector() )
    {
	  
      const int sz1 = size();
      const int sz2 = rhs.size();
      
      if ( rhs.is_int_vector() ) 
	{	  
	  if      ( is_int_vector() )    
	    {
	      for (int i=0; i<sz1; i++) for (int j=0; j<sz2; j++) if ( ivec[ve[i]] == rhs.ivec[rhs.ve[j]] ) return true;
	      return false;
	    }
	  
	  else if ( is_float_vector() )
	    {
	      for (int i=0; i<sz1; i++) for (int j=0; j<sz2; j++) if ( fvec[ve[i]] == rhs.ivec[rhs.ve[j]] ) return true; 
	      return false;
	    }
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )   
	    {
	      for (int i=0; i<sz1; i++) for (int j=0; j<sz2; j++) if ( bvec[ve[i]] == rhs.ivec[rhs.ve[j]] ) return true;
	      return false;
	    }
	  return Token();
	}
      else if ( rhs.is_float_vector() )
	{
	  if      ( is_int_vector() )    
	    {
	      for (int i=0; i<sz1; i++) for (int j=0; j<sz2; j++) if ( ivec[ve[i]] == rhs.fvec[rhs.ve[j]] ) return true;
	      return false;
	    }
	  else if ( is_float_vector() )  
	    {
	      for (int i=0; i<sz1; i++) for (int j=0; j<sz2; j++) if ( fvec[ve[i]] == rhs.fvec[rhs.ve[j]] ) return true;
	      return false;
	    }
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )  
	    {
	      for (int i=0; i<sz1; i++) for (int j=0; j< sz2; j++) if ( bvec[ve[i]] == rhs.fvec[rhs.ve[j]] ) return true;
	      return false;
	    }
	  return Token( );	 
	}
      else if ( rhs.is_bool_vector() )
	{
	  if      ( is_int_vector() )    
	    {
	      for (int i=0; i<sz1; i++) for (int j=0; j<sz2; j++) if ( ivec[ve[i]] == rhs.bvec[rhs.ve[j]] ) return true;
	      return false;
	    }
	  else if ( is_float_vector() )  
	    {
	      for (int i=0; i<sz1; i++)  for (int j=0; j<sz2; j++) if ( fvec[ve[i]] == rhs.bvec[rhs.ve[j]] ) return true;
	      return false;
	    }
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )   
	    {
	      for (int i=0; i<sz1; i++)  for (int j=0; j<sz2; j++) if ( bvec[ve[i]] == rhs.bvec[rhs.ve[j]] ) return true;
	      return false;
	    }
	  return Token();
	}	  
      else if ( rhs.is_string_vector() )
	{
	  if ( is_string_vector() )  
	    {
	      for (int i=0; i<sz1; i++) for (int j=0; j< sz2; j++) if ( svec[ve[i]] == rhs.svec[rhs.ve[j]] ) return true;
	      return false;
	    }
	  else return Token();
	  return Token();
	}	  	  
      else
	return Token();
    }
  
  //
  // Otherwise, handle cases involving at least one scalar
  //

  
  if ( ! is_set() ) return false;

  if ( ! rhs.is_set() ) return false;
  
  Token res = *this == rhs;

  if ( ! res.is_set() ) return false;

  if ( res.is_scalar() ) return res;

  if ( ! res.is_bool_vector() ) 
    Helper::halt("internal error");
  
  for (int i=0;i<res.size();i++) 
    if ( res.bvec[res.ve[i]] ) return Token(true);
  
  return Token(false);
  
}

Token Token::operator==(const Token & rhs) const
{

  // vector x vector comparison defined for same-length vectors (element-wise comparison)

  if ( is_vector() && rhs.is_vector() ) 
    {

      //
      // element-wise comparisons, if sizes match
      //
      
      if ( size() == rhs.size() ) 
	{

	  const int sz = size();
	  std::vector<bool> ans( sz );      
	  
	  if ( rhs.is_int_vector() ) 
	    {	  
	      if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] == rhs.ivec[rhs.ve[i]];
	      else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] == rhs.ivec[rhs.ve[i]]; 
	      else if ( is_string_vector() ) return Token();
	      else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] == rhs.ivec[rhs.ve[i]]; 
	      return Token( ans );
	    }
	  else if ( rhs.is_float_vector() )
	    {
	      if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] == rhs.fvec[rhs.ve[i]];
	      else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] == rhs.fvec[rhs.ve[i]]; 
	      else if ( is_string_vector() ) return Token();
	      else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] == rhs.fvec[rhs.ve[i]]; 
	      return Token( ans );	 
	    }
	  else if ( rhs.is_bool_vector() )
	    {
	      if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] == rhs.bvec[rhs.ve[i]];
	      else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] == rhs.bvec[rhs.ve[i]]; 
	      else if ( is_string_vector() ) return Token();
	      else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] == rhs.bvec[rhs.ve[i]]; 
	      return Token( ans );
	    }	  
	  else if ( rhs.is_string_vector() )
	    {
	      if ( is_string_vector() )  for (int i=0; i<sz; i++) ans[i] = svec[ve[i]] == rhs.svec[rhs.ve[i]]; 
	      else return Token();
	      return Token( ans );
	    }	  	  
	  else
	    return Token();
	}
      else 
	return Token();  // comparisons of different length vectors not allowed 

    }


  // vector == scalar 
  
  if ( is_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();

      std::vector<bool> ans( sz );      

      if ( rhs.is_int() ) 
	{
	  if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] == rhs.ival;
	  else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] == rhs.ival; 
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] == rhs.ival; 	  
	  return Token( ans );
	}
      else if ( rhs.is_float() )
	{
	  if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] == rhs.fval;
	  else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] == rhs.fval; 
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] == rhs.fval; 	  
	  return Token( ans );
	}
      else if ( rhs.is_bool() ) 
	{
	  if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] == rhs.bval;
	  else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] == rhs.bval; 
	  else if ( is_string_vector() ) return Token();
	  else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] == rhs.bval; 	  
	  return Token( ans );
	}
      else if ( rhs.is_string() ) 
	{
	  if ( is_string_vector() ) for (int i=0; i<sz; i++) ans[i] = svec[ve[i]] == rhs.sval;  
	  else return Token(); 
	  return Token( ans );
	}
      else return Token();
    }
  

  // scalar == vector

  if ( rhs.is_vector() )
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      
      if ( is_int() )
	{
	  if      ( rhs.is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = rhs.ivec[rhs.ve[i]] == ival;
	  else if ( rhs.is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = rhs.fvec[rhs.ve[i]] == ival;
	  else if ( rhs.is_string_vector() ) return Token();
	  else if ( rhs.is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = rhs.bvec[rhs.ve[i]] == ival;
	  return Token( ans );      
	}
      else if ( is_float() )
	{
	  if      ( rhs.is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = rhs.ivec[rhs.ve[i]] == fval;
	  else if ( rhs.is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = rhs.fvec[rhs.ve[i]] == fval;
	  else if ( rhs.is_string_vector() ) return Token();
	  else if ( rhs.is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = rhs.bvec[rhs.ve[i]] == fval;
	  return Token( ans );     
	}
      else if ( is_bool() )
	{
	  if      ( rhs.is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = rhs.ivec[rhs.ve[i]] == bval;
	  else if ( rhs.is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = rhs.fvec[rhs.ve[i]] == bval;
	  else if ( rhs.is_string_vector() ) return Token();
	  else if ( rhs.is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = rhs.bvec[rhs.ve[i]] == bval;
	  return Token( ans );     
	}
      else if ( is_string() )
	{
	  if ( rhs.is_string_vector() ) for (int i=0; i<sz; i++) ans[i] = rhs.svec[rhs.ve[i]] == sval;
	  else return Token();
	  return Token( ans );     
	}
      else return Token();
    }

            
  // scalar x scalar
  if ( is_bool()   && rhs.is_bool()   )  return Token( bval == rhs.bval ); 
  if ( is_int()    && rhs.is_int()    )  return Token( ival == rhs.ival ); 
  if ( is_float()  && rhs.is_float()  )  return Token( fval == rhs.fval ); 
  if ( is_string() && rhs.is_string() )  return Token( sval == rhs.sval ); 

  // also allow int / bool and  int / float comparisons
  if ( is_int()   && rhs.is_bool()    )  return Token( ival == (int)rhs.bval );
  if ( is_bool()  && rhs.is_int()     )  return Token( (int)bval == rhs.ival );

  if ( is_float() && rhs.is_bool()    )  return Token( (bool)fval == rhs.bval );
  if ( is_bool()  && rhs.is_float()   )  return Token( bval == (bool)rhs.fval );

  if ( is_float() && rhs.is_int()     )  return Token( fval == rhs.ival );
  if ( is_int()   && rhs.is_float()   )  return Token( ival == rhs.fval );

  return Token();

}



Token Token::operator+(const Token & rhs) const
{

  // vector x vector comparison defined for same-length vectors

  if ( is_vector() && rhs.is_vector() ) 
    {

      if ( size() != rhs.size() ) return Token(); 
      const int sz = size();      
      
      // concatenate strings
      if ( is_string_vector() && rhs.is_string_vector() )
	{
	  std::vector<std::string> ans( sz );      
	  for (int i=0; i<sz; i++) ans[i] = svec[ve[i]] + rhs.svec[rhs.ve[i]];
	  return Token( ans );            
	}
      
      if ( is_int_vector() ) 
	{
	  std::vector<int> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] + rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] + (int)rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] + (int)rhs.bvec[rhs.ve[i]];
	  else return Token();
	  return Token( ans );            
	}
      
      if ( is_float_vector() )
	{
	  std::vector<double> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] + rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] + rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] + (double)rhs.bvec[rhs.ve[i]];
	  return Token( ans );                  
	}
      
      else if ( is_bool_vector() )
	{
	  std::vector<double> ans( sz );      
	  if ( rhs.is_int_vector() )        for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] + rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] + rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] + (double)rhs.bvec[rhs.ve[i]];
	  return Token( ans );                  
	}

      return Token();
    }

  
  // concatenate strings ( in which case, A+B != B+A, so do before next step)
  if ( is_string_vector() && rhs.is_string() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<std::string> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = svec[ve[i]] + rhs.sval;
      return Token( ans );            
    }
  
  if ( is_string() && rhs.is_string_vector() )
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<std::string> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = sval + rhs.svec[rhs.ve[i]];
      return Token( ans );            
    }
  
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] + rhs.ival;      
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] + (int)rhs.bval;
      else if ( rhs.is_float() ) 
	{
	  std::vector<double> ans( sz );   	
	  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] + rhs.fval;
	  return Token( ans );            
	}
      else return Token();
      return Token( ans );            
    }
  
  if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival      + rhs.ivec[rhs.ve[i]];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (int)bval + rhs.ivec[rhs.ve[i]];
      else if ( is_float() ) 
	{
	  std::vector<double> ans( sz );      
	  for (int i=0; i<sz; i++) ans[i] = fval + rhs.ivec[rhs.ve[i]];
	  return Token( ans );            
	}
      else return Token();
      return Token( ans );            
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] + rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] + rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] + (double)rhs.bval;
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         + rhs.fvec[rhs.ve[i]];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         + rhs.fvec[rhs.ve[i]];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval + rhs.fvec[rhs.ve[i]];
      return Token( ans );            
    }

  else if ( is_bool_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] + rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] + rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] + (double)rhs.bval;
      return Token( ans );                  
    }

  else if ( rhs.is_bool_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         + rhs.bvec[rhs.ve[i]];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         + rhs.bvec[rhs.ve[i]];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval + rhs.bvec[rhs.ve[i]];
      return Token( ans );            
    }


  // scalar + scalar

  if ( is_int() ) 
    {
      if ( rhs.is_int() ) return Token( ival + rhs.ival );
      if ( rhs.is_bool() ) return Token( ival + rhs.bval );
      if ( rhs.is_float() ) return Token( ival + rhs.fval );
    }

  if ( is_float() ) 
    {
      if ( rhs.is_int() ) return Token( fval + rhs.ival );
      if ( rhs.is_bool() ) return Token( fval + rhs.bval );
      if ( rhs.is_float() ) return Token( fval + rhs.fval );
    }

  if ( is_bool() ) 
    {
      if ( rhs.is_int() ) return Token( bval + rhs.ival );
      if ( rhs.is_bool() ) return Token( (int)bval + (int)rhs.bval );
      if ( rhs.is_float() ) return Token( bval + rhs.fval );
    }

  if ( is_string() ) // concatenate
    {
      if ( rhs.is_string() ) return Token( sval + rhs.sval );
    }

  return Token();
}

Token Token::operator-(const Token & rhs) const
{

  // vector x vector comparison defined for same-length vectors
  if ( is_vector() && rhs.is_vector() ) 
    {
      
      if ( size() != rhs.size() ) return Token();
      const int sz = size();      
      
      if ( is_int_vector() ) 
	{
	  std::vector<int> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] - rhs.ivec[rhs.ve[i]];	  
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] - (int)rhs.bvec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) 
	    {
	      std::vector<double> ans( sz );      
	      for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] - rhs.fvec[rhs.ve[i]];
	      return Token( ans );            
	    }
	  else return Token();
	  return Token( ans );            
	}
      
      if ( is_float_vector() )
	{
	  std::vector<double> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] - rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] - rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] - (double)rhs.bvec[rhs.ve[i]];
	  return Token( ans );                  
	}
      
      else if ( is_bool_vector() )
	{
	  std::vector<double> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] - rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] - rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] - (double)rhs.bvec[rhs.ve[i]];
	  return Token( ans );                  
	}
      
      return Token();
 
    }


  // vector x scalar ops
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] - rhs.ival;     
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] - (int)rhs.bval;
      else if ( rhs.is_float() ) 
	{
	  std::vector<double> ans( sz );      
	  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] - rhs.fval;
	  return Token( ans );            
	}
      else return Token();
      return Token( ans );            
    }
  
  if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival      - rhs.ivec[rhs.ve[i]];      
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (int)bval - rhs.ivec[rhs.ve[i]];
      else if ( is_float() ) 
	{
	  std::vector<double> ans( sz );      
	  for (int i=0; i<sz; i++) ans[i] = fval - rhs.ivec[rhs.ve[i]];
	  return Token( ans );            
	}
      else return Token();
      return Token( ans );            
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] - rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] - rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] - (double)rhs.bval;
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         - rhs.fvec[rhs.ve[i]];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         - rhs.fvec[rhs.ve[i]];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval - rhs.fvec[rhs.ve[i]];
      return Token( ans );            
    }


  // scalar + scalar

  if ( is_int() ) 
    {
      if ( rhs.is_int()   ) return Token( ival - rhs.ival );
      if ( rhs.is_bool()  ) return Token( ival - rhs.bval );
      if ( rhs.is_float() ) return Token( ival - rhs.fval );
    }

  if ( is_float() ) 
    {
      if ( rhs.is_int()   ) return Token( fval - rhs.ival );
      if ( rhs.is_bool()  ) return Token( fval - rhs.bval );
      if ( rhs.is_float() ) return Token( fval - rhs.fval );
    }

  if ( is_bool() ) 
    {
      if ( rhs.is_int()   ) return Token( bval - rhs.ival );
      if ( rhs.is_bool()  ) return Token( (int)bval - (int)rhs.bval );
      if ( rhs.is_float() ) return Token( bval - rhs.fval );
    }

  return Token();

}

Token Token::operator*(const Token & rhs) const
{

  // vector x vector comparison defined for same-length vectors

  if ( is_vector() && rhs.is_vector() ) 
    {

      if ( size() != rhs.size() ) return Token();
      const int sz = size();      
      
      if ( is_int_vector() ) 
	{
	  std::vector<int> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] * rhs.ivec[rhs.ve[i]];	  
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] * (int)rhs.bvec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) 
	    {
	      std::vector<double> ans( sz );      	      
	      for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] * rhs.fvec[rhs.ve[i]];
	      return Token( ans );
	    }	  
	  else return Token();
	  return Token( ans );            
	}
      
      if ( is_float_vector() )
	{
	  std::vector<double> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] * rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] * (double)rhs.bvec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] * rhs.fvec[rhs.ve[i]];
	  return Token( ans );                  
	}
      
      else if ( is_bool_vector() )
	{
	  std::vector<double> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] * rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] * rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] * (double)rhs.bvec[rhs.ve[i]];
	  return Token( ans );                  
	}

      return Token();

    }


  // vector x scalar ops
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] * rhs.ival;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] * (int)rhs.bval;
      else if ( rhs.is_float() ) 
	{
	  std::vector<double> ans( sz );      
	  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] * rhs.fval;
	  return Token( ans );            
	}
      else return Token();
      return Token( ans );            
    }
  
  if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival      * rhs.ivec[rhs.ve[i]];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (int)bval * rhs.ivec[rhs.ve[i]];
      else if ( is_float() ) 
	{
	  std::vector<double> ans( sz );      
	  for (int i=0; i<sz; i++) ans[i] = fval * rhs.ivec[rhs.ve[i]];
	  return Token( ans );            
	}
      else return Token();
      return Token( ans );            
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] * rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] * rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] * (double)rhs.bval;
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         * rhs.fvec[rhs.ve[i]];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         * rhs.fvec[rhs.ve[i]];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval * rhs.fvec[rhs.ve[i]];
      return Token( ans );            
    }

  else if ( is_bool_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] * rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] * rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] * (double)rhs.bval;
      return Token( ans );                  
    }

  else if ( rhs.is_bool_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         * rhs.bvec[rhs.ve[i]];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         * rhs.bvec[rhs.ve[i]];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval * rhs.bvec[rhs.ve[i]];
      return Token( ans );            
    }

  // scalar * scalar

  if ( is_int() ) 
    {
      if ( rhs.is_int()   ) return Token( ival * rhs.ival );
      if ( rhs.is_float() ) return Token( ival * rhs.fval );
      if ( rhs.is_bool()  ) return Token( ival * rhs.bval );
    }

  if ( is_float() ) 
    {
      if ( rhs.is_int()   ) return Token( fval * rhs.ival );
      if ( rhs.is_float() ) return Token( fval * rhs.fval );
      if ( rhs.is_bool()  ) return Token( fval * rhs.bval );
    }

  if ( is_bool() ) 
    {
      if ( rhs.is_int()   ) return Token( bval * rhs.ival );
      if ( rhs.is_float() ) return Token( bval * rhs.fval );
      if ( rhs.is_bool()  ) return Token( bval * rhs.bval );
    }

  return Token();
}

Token Token::operator^(const Token & rhs) const
{
  Helper::halt("^ operator not supported, use pow() or sqr()" );

  // if ( rhs.is_vector() ) Helper::halt( "not allowed vector expression 'x' ^ vector" );

  // vector ^ scalar
  // scalar ^ scalar

  // if ( is_int_vector() ) 
  //   {
  //     const int sz = size();
  //     if ( sz == 0 ) return Token();
  //     std::vector<double> ans( sz );      
  //     if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = pow( ivec[i] , rhs.ival );
  //     else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = pow( ivec[i] , rhs.fval );
  //     return Token( ans );
  //   }

  // if ( is_float_vector() ) 
  //   {
  //     const int sz = size();
  //     if ( sz == 0 ) return Token();
  //     std::vector<double> ans( sz );      
  //     if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = pow( fvec[i] , rhs.ival );
  //     else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = pow( fvec[i] , rhs.fval );
  //     return Token( ans );
  //   }

  // if ( is_int() ) 
  //   {
  //     if ( rhs.is_int() ) return Token( pow( ival , rhs.ival ) );
  //     if ( rhs.is_float() ) return Token( pow( ival , rhs.fval ) );
  //   }
  // if ( is_float() ) 
  //   {
  //     if ( rhs.is_int() ) return Token( pow( fval , rhs.ival ) );
  //     if ( rhs.is_float() ) return Token( pow( fval , rhs.fval ) );
  //   }
  return Token();
}

Token Token::operator/(const Token & rhs) const
{
  
  // always return a float
  // allow bool as numerator, but not demoninator
  // no strings

  // vector x vector comparison defined for same length vector
  if ( is_vector() && rhs.is_vector() ) 
    {

      if ( size() != rhs.size() ) return Token();
      const int sz = size();      
      
      if ( is_int_vector() ) 
	{
	  std::vector<double> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] / rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] / (int)rhs.fvec[rhs.ve[i]];
	  else return Token();
	  return Token( ans );            
	}
      
      if ( is_float_vector() )
	{
	  std::vector<double> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] / (double)rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] / rhs.fvec[rhs.ve[i]];	  
	  return Token( ans );                  
	}
      
      else if ( is_bool_vector() )
	{
	  std::vector<double> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = (int)bvec[ve[i]] / (double)rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = (int)bvec[ve[i]] / rhs.fvec[rhs.ve[i]];	  
	  return Token( ans );                  
	}

      return Token();

    }

  // vector / scalar ops
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] / (double)rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] / rhs.fval;
      return Token( ans );            
    }
  
  if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival / (double)rhs.ivec[rhs.ve[i]];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval / (double)rhs.ivec[rhs.ve[i]];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (int)bval / (double)rhs.ivec[rhs.ve[i]];
      return Token( ans );
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] / (double)rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] / rhs.fval;      
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         / rhs.fvec[rhs.ve[i]];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         / rhs.fvec[rhs.ve[i]];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval / rhs.fvec[rhs.ve[i]];
      return Token( ans );            
    }
  
  else if ( is_bool_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] / (double)rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] / rhs.fval;
      
      return Token( ans );                  
    }
  

  // scalar / scalar  (always return a double)

  if ( is_int() ) 
    {
      if ( rhs.is_int()   ) return Token( ival / (double)rhs.ival );
      if ( rhs.is_float() ) return Token( ival / rhs.fval );
    }

  if ( is_float() ) 
    {
      if ( rhs.is_int()   ) return Token( fval / (double)rhs.ival );
      if ( rhs.is_float() ) return Token( fval / rhs.fval );
    }

   if ( is_bool() ) 
     {
       if ( rhs.is_int()   ) return Token( bval / (double)rhs.ival );
       if ( rhs.is_float() ) return Token( bval / rhs.fval );
     }

  return Token();

}

Token Token::operator%(const Token & rhs) const
{

  if ( rhs.is_vector() ) Helper::halt( "not allowed vector expression 'x' % vector" );
  if ( ! rhs.is_int() ) return Token();

  // scalar % scalar
  // vector % scalar
  
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = (int)(ivec[ve[i]] % rhs.ival );
      return Token(ans);
    }
  
  if ( is_int() ) return Token( (int)(ival % rhs.ival) );

  return Token();
}


Token Token::operator<(const Token & rhs) const
{

  // vector x vector comparison defined for same length entities
  if ( is_vector() && rhs.is_vector() ) 
    {

      if ( size() != rhs.size() ) return Token();
      const int sz = size();      
      
      // allow string-to-string comparison
      if ( is_string_vector() && rhs.is_string_vector() )
	{
	  std::vector<bool> ans( sz );      
	  for (int i=0; i<sz; i++) ans[i] = svec[ve[i]] < rhs.svec[rhs.ve[i]];
	  return Token( ans );            
	}
      
      if ( is_int_vector() ) 
	{
	  std::vector<bool> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] < rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] < rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] < (int)rhs.bvec[rhs.ve[i]];
	  else return Token();
	  return Token( ans );            
	}
      
      if ( is_float_vector() )
	{
	  std::vector<bool> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] < rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] < rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] < (double)rhs.bvec[rhs.ve[i]];
	  return Token( ans );                  
	}
      
      else if ( is_bool_vector() )
	{
	  std::vector<bool> ans( sz );      
	  if ( rhs.is_int_vector() )        for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] < rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] < rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] < (double)rhs.bvec[rhs.ve[i]];
	  return Token( ans );
	}

      return Token();

    }
 
 
  // vector x scalar ops
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] < rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] < rhs.fval;
      return Token( ans );            
    }
  
  else if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival < rhs.ivec[rhs.ve[i]];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval < rhs.ivec[rhs.ve[i]];
      return Token( ans ); 
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] < rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] < rhs.fval;
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         < rhs.fvec[rhs.ve[i]];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         < rhs.fvec[rhs.ve[i]];
      return Token( ans );            
    }


  else if ( is_string_vector() && rhs.is_string() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = svec[ve[i]] < rhs.sval;
      return Token( ans );            
    }
  
  else if ( is_string() && rhs.is_string_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = sval < rhs.svec[ rhs.ve[i] ];;
      return Token( ans );            
    }

  
  // scalar < scalar  
  else if ( is_int() ) 
    {
      if ( rhs.is_int() )   return Token( ival < rhs.ival );
      if ( rhs.is_float() ) return Token( ival < rhs.fval );
      if ( rhs.is_bool() )  return Token( ival < rhs.bval );
    }
  
  else if ( is_float() ) 
    {
      if ( rhs.is_int() )   return Token( fval < rhs.ival );
      if ( rhs.is_float() ) return Token( fval < rhs.fval );
      if ( rhs.is_bool() )  return Token( fval < rhs.bval );
    }
  
  else if ( is_bool() )
    {
      if ( rhs.is_bool() )   return Token( bval < rhs.bval );
      if ( rhs.is_int()   )  return Token( bval < rhs.ival );
      if ( rhs.is_float() )  return Token( bval < rhs.fval );
    }

  else if ( is_string() ) 
    {
      if ( rhs.is_string() )  return Token( sval < rhs.sval );
    }
  
  return Token();
}


Token Token::operator>(const Token & rhs) const
{

  // vector x vector comparison not defined
  if ( is_vector() && rhs.is_vector() ) 
    {
      if ( size() != rhs.size() ) return Token();
      const int sz = size();      
      
      // allow string-to-string comparison
      if ( is_string_vector() && rhs.is_string_vector() )
	{
	  std::vector<bool> ans( sz );      
	  for (int i=0; i<sz; i++) ans[i] = svec[ve[i]] > rhs.svec[rhs.ve[i]];
	  return Token( ans );            
	}
      
      if ( is_int_vector() ) 
	{
	  std::vector<bool> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] > rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] > rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] > (int)rhs.bvec[rhs.ve[i]];
	  else return Token();
	  return Token( ans );            
	}
      
      if ( is_float_vector() )
	{
	  std::vector<bool> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] > rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] > rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] > (double)rhs.bvec[rhs.ve[i]];
	  return Token( ans );                  
	}
      
      else if ( is_bool_vector() )
	{
	  std::vector<bool> ans( sz );      
	  if ( rhs.is_int_vector() )        for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] > rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_float_vector() ) for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] > rhs.fvec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] > (double)rhs.bvec[rhs.ve[i]];
	  return Token( ans );                  
	}

      return Token();

    }
  
  // vector x scalar ops
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] > rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] > rhs.fval;
      return Token( ans );            
    }
  
  else if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival      > rhs.ivec[rhs.ve[i]];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] =      fval > rhs.ivec[rhs.ve[i]];
      return Token( ans );            
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] > rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[ve[i]] > rhs.fval;
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         > rhs.fvec[rhs.ve[i]];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         > rhs.fvec[rhs.ve[i]];
      return Token( ans );            
    }

  else if ( is_string_vector() && rhs.is_string() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = svec[ve[i]] > rhs.sval;
      return Token( ans );            
    }
  
  else if ( is_string() && rhs.is_string_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = sval > rhs.svec[rhs.ve[i]];
      return Token( ans );            
    }



  if ( is_int() ) 
    {
      if ( rhs.is_int() )  return Token( ival > rhs.ival );
      if ( rhs.is_float() ) return Token( ival > rhs.fval );
      if ( rhs.is_bool() ) return Token( ival > rhs.bval );
    }
  
  if ( is_float() ) 
    {
      if ( rhs.is_int() )  return Token( fval > rhs.ival );
      if ( rhs.is_float() ) return Token( fval > rhs.fval );
      if ( rhs.is_bool() ) return Token( fval > rhs.bval );
    }
  
  else if ( is_bool() )
    {
      if ( rhs.is_bool() )  return Token( bval > rhs.bval );
      if ( rhs.is_int() )   return Token( bval > rhs.ival );
      if ( rhs.is_float() ) return Token( bval > rhs.fval );

    }

  if ( is_string() ) 
    {
      if ( rhs.is_string() )  return Token( sval > rhs.sval );
    }
  

  return Token();
}

Token Token::operator>=(const Token & rhs) const
{
  return ! ( *this < rhs );
}

Token Token::operator<=(const Token & rhs) const
{
  return ! ( *this > rhs );
}


Token Token::operator&&(const Token & rhs) const
{

  // if *either* side is NULL, return NULL
  if ( ! ( is_set() && rhs.is_set() ) ) return Token();

  // only defined for int and bool

  // vector x vector comparison defined for same-length vectors
  if ( is_vector() && rhs.is_vector() ) 
    {
      if ( size() != rhs.size() ) return Token();
      const int sz = size();      
      
      if ( is_int_vector() ) 
	{
	  std::vector<bool> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] && rhs.ivec[rhs.ve[i]];	  
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] && rhs.bvec[rhs.ve[i]];
	  else return Token();
	  return Token( ans );            
	}
      
      else if ( is_bool_vector() )
	{
	  std::vector<bool> ans( sz );      
	  if ( rhs.is_int_vector() )        for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] && rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] && rhs.bvec[rhs.ve[i]];
	  return Token( ans );                  
	}

      return Token();

    }
  
  // TODO: Note -- currently no vector x scalar implementataion....


  // scalar AND scalar

  if ( is_bool() && rhs.is_bool() ) return Token( bval && rhs.bval );
  if ( is_bool() && rhs.is_int()  ) return Token( bval && rhs.ival );
  if ( is_int( ) && rhs.is_bool() ) return Token( ival && rhs.bval );
  if ( is_int( ) && rhs.is_int()  ) return Token( ival && rhs.ival );

  // undefined for other types
  return Token();
  
}


Token Token::operator||(const Token & rhs) const
{
  
  // if *both* sides NULL, return NULL
  if ( ! ( is_set() || rhs.is_set() ) ) return Token();
  
  // vector x vector comparison not defined
  if ( is_vector() && rhs.is_vector() ) 
    {
      if ( size() != rhs.size() ) return Token();
      const int sz = size();
      
      if ( is_int_vector() ) 
	{
	  std::vector<bool> ans( sz );      
	  if      ( rhs.is_int_vector() )   for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] || rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = ivec[ve[i]] || rhs.bvec[rhs.ve[i]];
	  else return Token();
	  return Token( ans );            
	}
      
      else if ( is_bool_vector() )
	{
	  std::vector<bool> ans( sz );      
	  if ( rhs.is_int_vector() )        for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] || rhs.ivec[rhs.ve[i]];
	  else if ( rhs.is_bool_vector() )  for (int i=0; i<sz; i++) ans[i] = bvec[ve[i]] || rhs.bvec[rhs.ve[i]];
	  return Token( ans );                  
	}

      return Token();
    }

  // TODO: no vector x scalar ops

  
  //
  // scalar OR scalar
  //
  
  bool left_valid = is_bool() || is_int() ;
  bool right_valid = rhs.is_bool() || rhs.is_int();
  
  // only NULL if BOTH sides are invalid (not bool or int)
  if ( ( ! left_valid ) && ( ! right_valid ) ) return Token();
  
  bool left_val  = left_valid ? ( is_bool() ? bval : ival ) : false ; 
  bool right_val = right_valid ? ( rhs.is_bool() ? rhs.bval : rhs.ival ) : false ;
  
  return Token( left_val || right_val );
  
}


Token Token::operands( Token & t)
{
  if ( ttype == NOT_OPERATOR ) return !t;
  else return Token();
}

Token Token::operands( Token & right, Token & left )
{
  switch( ttype )
    {
    case ASSIGNMENT_OPERATOR : return right;
    case ADD_OPERATOR : return left + right ;      
    case SUBTRACT_OPERATOR : return left - right ;
    case MULTIPLY_OPERATOR : return left * right ;
    case DIVIDE_OPERATOR : return left / right ;
    case MOD_OPERATOR : return left % right ;
    case AND_OPERATOR : return left && right ;
    case OR_OPERATOR : return left || right ;
    case LESS_THAN_OPERATOR : return left < right ;
    case LESS_THAN_OR_EQUAL_OPERATOR : return left <= right ;
    case GREATER_THAN_OPERATOR : return left > right ;
    case GREATER_THAN_OR_EQUAL_OPERATOR : return left >= right ;
    case EQUAL_OPERATOR : return left == right ;
    case HAS_OPERATOR : return left.contains( right );
    case UNEQUAL_OPERATOR : return left != right ;
    default : return Token();
}
  return Token();
}



int Token::as_int() const
{

  switch ( ttype ) 
    {
    case INT : return ival;
    case FLOAT : return (int)fval;
    case BOOL : return bval ? 1 : 0;
    case STRING :
      {
	int i;
	if ( Helper::from_string<int>( i , sval, std::dec ) ) 
	  return i;
      }
    default: return 0; // UNDEF
    }
  
  return 0; // UNDEF
}


double Token::as_float() const
{
 
  switch ( ttype ) 
    {
    case INT : return (double)ival;
    case FLOAT : return fval;
    case BOOL : return bval ? 1.0 : 0.0;
    case STRING :
      {
	double d;
	if ( Helper::from_string<double>( d , sval, std::dec ) ) 
	  return d;
      }
    default : return 0.0;
   }
  return 0.0; // UNDEF   
}


std::string Token::as_string() const
{
  // unlike the above scalar as_type() functions, 
  // this automatically converts vectors to strings, comma-delimited

  if ( ttype == STRING ) return sval;

  std::stringstream ss;

  if ( ttype == INT ) ss << ival;
  else if ( ttype == FLOAT ) ss << fval;
  else if ( ttype == BOOL ) ss << ( bval ? "true" : "false" );
  else if ( ttype == STRING_VECTOR )
    {
      for (int i=0;i<ve.size();i++) 
	ss << ( i ? "," : "" ) << svec[ve[i]]; 	
    }
  else if ( ttype == INT_VECTOR ) 
    {
      for (int i=0;i<ve.size();i++) 
	ss << ( i ? "," : "" ) << ivec[ve[i]]; 	
    }
  else if ( ttype == FLOAT_VECTOR ) 
    {
      for (int i=0;i<ve.size();i++) 
	ss << ( i ? "," : "" ) << fvec[ve[i]]; 	
    }
  else if ( ttype == BOOL_VECTOR )
    {
      for (int i=0;i<ve.size();i++) 
	ss << ( i ? "," : "" ) << ( bvec[ve[i]] ? "true" : "false" );
    }
  else 
    ss << ".";

  return ss.str();
}




bool Token::as_bool() const
{

  // unlike as_int and as_float, as_bool automatically converts boolean vectors
  // to a scalar vector, based on any(T)

  if      ( ttype == BOOL )   return bval;    
  else if ( ttype == INT )    return ival;
  else if ( ttype == FLOAT )  return fval;
  else if ( ttype == STRING ) return !( sval == "" || sval == "." || sval == "0" || sval == "false" || sval == "FALSE" ) ;

  else if ( ttype == BOOL_VECTOR )  for (int i=0;i<ve.size(); i++) { if ( bvec[ve[i]] ) return true; } 
  else if ( ttype == INT_VECTOR )   for (int i=0;i<ve.size(); i++) { if ( ivec[ve[i]] ) return true; }
  else if ( ttype == FLOAT_VECTOR ) for (int i=0;i<ve.size(); i++) { if ( fvec[ve[i]] ) return true; }
  else if ( ttype == STRING_VECTOR ) 
    for (int i=0;i<ve.size(); i++)
      {
	if ( ! ( svec[ve[i]] == "." || svec[ve[i]] == ""
		 || svec[ve[i]] == "0" || svec[ve[i]] == "false"
		 || svec[ve[i]] == "FALSE") )
	  return true;
      }

  return false;
}



int Token::int_element(const int i) const
{
  if ( i < 0 || i >= size() ) 
    Helper::halt( "out of range for " + name() + " (" + Helper::int2str(i+1) + " of " + Helper::int2str(size() ) +")" );
  if ( ttype == INT_VECTOR ) return ivec[ve[i]];
  if ( ttype == INT ) return ival;
  return 0;      
}

double Token::float_element(const int i) const
{
  if ( i < 0 || i >= size() ) 
    Helper::halt( "out of range for " + name() + " (" + Helper::int2str(i+1) + " of " + Helper::int2str(size() ) +")" );
  if ( ttype == FLOAT_VECTOR ) return fvec[ve[i]];
  if ( ttype == FLOAT ) return fval;
  return 0;      
}


std::string Token::string_element(const int i) const
{
  if ( i < 0 || i >= size() ) 
    Helper::halt( "out of range for " + name() + " (" + Helper::int2str(i+1) + " of " + Helper::int2str(size() ) +")" );
  if ( ttype == STRING_VECTOR ) return svec[ve[i]];
  if ( ttype == STRING ) return sval;
  return ".";
}

bool Token::bool_element(const int i) const
{
  if ( i < 0 || i >= size() ) 
    Helper::halt( "out of range for " + name() + " (" + Helper::int2str(i+1) + " of " + Helper::int2str(size() ) +")" );
  if ( ttype == BOOL_VECTOR ) return bvec[ve[i]];
  if ( ttype == BOOL ) return bval;
  return false;  
}



int Token::as_int_element(const int i) const
{
  if ( i < 0 || i >= size() ) 
    Helper::halt( "out of range for " + name() + " (" + Helper::int2str(i+1) + " of " + Helper::int2str(size() ) +")" );  
  if ( ttype == INT_VECTOR ) return ivec[ve[i]];
  if ( ttype == INT ) return ival;
  if ( ttype == FLOAT_VECTOR ) return (int)fvec[ve[i]];
  if ( ttype == FLOAT ) return (int)fval;  
  if ( ttype == BOOL_VECTOR ) return bvec[ve[i]];
  if ( ttype == BOOL ) return bval;  
  return 0;      
}

double Token::as_float_element(const int i) const
{
  if ( i < 0 || i >= size() ) 
    Helper::halt( "out of range for " + name() + " (" + Helper::int2str(i+1) + " of " + Helper::int2str(size() ) +")" );  
  if ( ttype == FLOAT_VECTOR ) return fvec[ve[i]];
  if ( ttype == FLOAT ) return fval;
  if ( ttype == INT_VECTOR ) return ivec[ve[i]];
  if ( ttype == INT ) return ival;
  if ( ttype == BOOL_VECTOR ) return bvec[ve[i]];
  if ( ttype == BOOL ) return bval;  
  return 0;      
}


std::string Token::as_string_element(const int i) const
{

  if ( i < 0 || i >= size() ) 
    Helper::halt( "out of range for " + name() + " (" + Helper::int2str(i+1) + " of " + Helper::int2str(size() ) +")" );  

  if ( ttype == STRING_VECTOR ) return svec[ve[i]];
  if ( ttype == STRING ) return sval;

  if ( ttype == INT_VECTOR ) return Helper::int2str( ivec[ve[i]] );
  if ( ttype == INT ) return Helper::int2str( ival );

  if ( ttype == FLOAT_VECTOR ) return Helper::dbl2str( fvec[ve[i]] );
  if ( ttype == FLOAT ) return Helper::dbl2str( fval );

  if ( ttype == BOOL_VECTOR ) return bvec[ve[i]] ? "true" : "false" ; 
  if ( ttype == BOOL ) return bval ? "true" : "false" ; 

  return ".";      
}

bool Token::string2bool( const std::string & sval ) const
{ 
  return !( sval == "" || sval == "." || sval == "0" || sval == "false" || sval == "FALSE" ) ; 
}
  
bool Token::as_bool_element(const int i) const
{
  if ( i < 0 || i >= size() ) 
    Helper::halt( "out of range for " + name() + " (" + Helper::int2str(i+1) + " of " + Helper::int2str(size() ) +")" );
  if ( ttype == BOOL_VECTOR ) return bvec[ve[i]];
  if ( ttype == BOOL ) return bval;
  if ( ttype == INT_VECTOR ) return ivec[ve[i]];
  if ( ttype == INT ) return ival;
  if ( ttype == FLOAT_VECTOR ) return fvec[ve[i]];
  if ( ttype == FLOAT ) return fval;
  if ( ttype == STRING_VECTOR ) return string2bool( svec[ve[i]] ); 
  if ( ttype == STRING ) return string2bool( sval ); 
  return false;  
}



std::vector<int> Token::as_int_vector() const
{

  if ( ttype == INT_VECTOR && ! is_masked() ) return ivec;

  std::vector<int> ans( size() );
  
  if ( ttype == FLOAT_VECTOR ) 
    {
      for (int i=0; i<ve.size(); i++) ans[i] = (int)fvec[ve[i]];
      return ans;
    }

    if ( ttype == INT_VECTOR ) 
    {
      for (int i=0; i<ve.size(); i++) ans[i] = (int)ivec[ve[i]];
      return ans;
    }

  if ( ttype == BOOL_VECTOR ) 
    {
      for (int i=0; i<ve.size(); i++) ans[i] = bvec[ve[i]] ? 1 : 0 ;
      return ans;
    }
  
  if ( ttype == STRING_VECTOR ) 
    {
      for (int i=0; i<ve.size(); i++) 
	if ( ! Helper::from_string<int>( ans[i] , svec[ve[i]], std::dec ) ) ans[i] = 0;	
      return ans;
    }

  switch ( ttype ) 
    {
    case INT   : ans[0] = ival; return ans;
    case FLOAT : ans[0] = (int)fval; return ans;
    case BOOL  : ans[0] = bval ? 1 : 0; return ans;
    case STRING : 
      {      
      if ( ! Helper::from_string<int>( ans[0] , sval, std::dec ) ) ans[0] = 0;
      return ans;
      }
    default : return ans; // 0-size vector
    }
  
  // should be size 0 vector
  return ans; // UNDEF
  
}

std::vector<double> Token::as_float_vector() const
{
  if ( ttype == FLOAT_VECTOR && ! is_masked() ) return fvec;

  std::vector<double> ans( size() );
  
  if ( ttype == INT_VECTOR ) 
    {
      for (int i=0; i<ve.size(); i++) ans[i] = ivec[ve[i]];
      return ans;
    }
  
  if ( ttype == FLOAT_VECTOR )
    {
      for (int i=0; i<ve.size(); i++) ans[i] = fvec[ve[i]];
      return ans;
    }

  if ( ttype == BOOL_VECTOR ) 
    {
      for (int i=0; i<ve.size(); i++) ans[i] = bvec[ve[i]] ? 1 : 0 ;
      return ans;
    }
  
  if ( ttype == STRING_VECTOR ) 
    {
      for (int i=0; i<ve.size(); i++) 
	if ( ! Helper::from_string<double>( ans[i] , svec[ve[i]], std::dec ) ) ans[i] = 0;
      return ans;
    }
  
  switch ( ttype ) 
    {
    case INT   : ans[0] = ival; return ans;
    case FLOAT : ans[0] = fval; return ans;
    case BOOL  : ans[0] = bval ? 1 : 0; return ans;
    case STRING: 
    {      
      if ( ! Helper::from_string<double>( ans[0] , sval, std::dec ) ) ans[0] = 0;
      return ans;
    }
    default : return ans; 
    }
  
  // should be size 0 vector
  return ans; // UNDEF

}


std::vector<std::string> Token::as_string_vector() const
{
  if ( ttype == STRING_VECTOR && ! is_masked() ) return svec;

  std::vector<std::string> ans;
  
  if ( is_scalar() ) { ans.push_back( as_string() ); return ans; }
  
  ans.resize( size() );
  for (int i=0; i<size(); i++) 
    ans[i] = as_string_element( i );
  
  return ans;
}


std::vector<bool> Token::as_bool_vector() const
{
  if ( ttype == BOOL_VECTOR && ! is_masked() ) return bvec;

  std::vector<bool> ans;
  if ( is_scalar() ) { ans.push_back( as_bool() ); return ans; }
  
  ans.resize( size() );
  if      ( ttype == INT_VECTOR )   { for (int i=0; i<ve.size(); i++) ans[i] = ivec[ve[i]]; }  
  else if ( ttype == FLOAT_VECTOR ) { for (int i=0; i<ve.size(); i++) ans[i] = fvec[ve[i]]; }
  else if ( ttype == BOOL_VECTOR ) { for (int i=0; i<ve.size(); i++) ans[i] = bvec[ve[i]]; }
  else if ( ttype == STRING_VECTOR ) { for (int i=0; i<ve.size(); i++) ans[i] = Helper::yesno( svec[ve[i]] ); }
  
  return ans;
}



//
// Token functions
//


Token TokenFunctions::fn_set( const Token & tok ) const
{
    return tok.is_set();
} 

Token TokenFunctions::fn_notset( const Token & tok ) const
{
    return ! tok.is_set();
} 

Token TokenFunctions::fn_sqrt( const Token & tok ) const
{
  if ( tok.is_int() ) return Token( sqrt( tok.as_int() ) );
  if ( tok.is_float() ) return Token( sqrt( tok.as_float() ) );
  if ( tok.is_int_vector() || tok.is_float_vector() ) 
    {
      std::vector<double> ans = tok.as_float_vector();
      for (int i=0; i<ans.size(); i++) ans[i] = sqrt( ans[i] );
      return Token( ans );
    }
  return Token();
}


Token TokenFunctions::fn_rnd() const
{
  return Token( CRandom::rand() );
}


Token TokenFunctions::fn_rnd( const Token & tok ) const
{
  
  // make int return value, base-1
  if ( tok.is_int() || tok.is_float() ) 
    return Token( 1 + CRandom::rand( tok.as_int() ) );
  
  // oherwise 0..1
  return fn_rnd();
}

Token TokenFunctions::fn_floor( const Token & tok ) const
{
  // only applies to numeric data
  if ( tok.is_float() )
    return Token( floor( tok.as_float() ) );
  
  if ( tok.is_float_vector() ) 
    {
      std::vector<double> ans = tok.as_float_vector();
      for (int i=0; i<ans.size(); i++) ans[i] = floor( ans[i] );
      return Token( ans );
    }
  return Token();
}


Token TokenFunctions::fn_abs( const Token & tok ) const
{
  // only applies to numeric data
  if ( tok.is_float() )
    return Token( fabs( tok.as_float() ) );
  
  if ( tok.is_int() )
    return Token( abs( tok.as_int() ) );

  if ( tok.is_float_vector() ) 
    {
      std::vector<double> ans = tok.as_float_vector();
      for (int i=0; i<ans.size(); i++) ans[i] = fabs( ans[i] );
      return Token( ans );
    }

  if ( tok.is_int_vector() ) 
    {
      std::vector<int> ans = tok.as_int_vector();
      for (int i=0; i<ans.size(); i++) ans[i] = abs( ans[i] );
      return Token( ans );
    }

  return Token();
}

Token TokenFunctions::fn_round( const Token & tok ) const
{
  // only applies to numeric data
  if ( tok.is_float() )
    return Token( round( tok.as_float() ) );
  
  if ( tok.is_float_vector() ) 
    {
      std::vector<double> ans = tok.as_float_vector();
      for (int i=0; i<ans.size(); i++) ans[i] = round( ans[i] );
      return Token( ans );
    }
  return Token();
}



Token TokenFunctions::fn_log( const Token & tok ) const
{
  if ( tok.is_int() ) return Token( log( tok.as_int() ) );
  if ( tok.is_float() ) return Token( log( tok.as_float() ) );
  if ( tok.is_int_vector() || tok.is_float_vector() ) 
    {
      std::vector<double> ans = tok.as_float_vector();
      for (int i=0; i<ans.size(); i++) ans[i] = log( ans[i] );
      return Token( ans );
    }
  return Token();
}


Token TokenFunctions::fn_log10( const Token & tok ) const
{
  if ( tok.is_int() ) return Token( log10( tok.as_int() ) );
  if ( tok.is_float() ) return Token( log10( tok.as_float() ) );
  if ( tok.is_int_vector() || tok.is_float_vector() ) 
    {
      std::vector<double> ans = tok.as_float_vector();
      for (int i=0; i<ans.size(); i++) ans[i] = log10( ans[i] );
      return Token( ans );
    }
  return Token();
}

Token TokenFunctions::fn_exp( const Token & tok ) const
{
  if ( tok.is_int() ) return Token( exp( tok.as_int() ) );
  if ( tok.is_float() ) return Token( exp( tok.as_float() ) );
  if ( tok.is_int_vector() || tok.is_float_vector() ) 
    {
      std::vector<double> ans = tok.as_float_vector();
      for (int i=0; i<ans.size(); i++) ans[i] = exp( ans[i] );
      return Token( ans );
    }
  return Token();
}


Token TokenFunctions::fn_pow( const Token & tok , const Token & tok2 ) const
{
  bool is_int_type = tok.is_int() || tok.is_int_vector();
  bool is_float_type = tok.is_float() || tok.is_float_vector();
  if ( ! ( is_int_type || is_float_type ) ) return Token();
  if ( ! ( tok2.is_int() || tok2.is_float() ) ) return Token();
  
  if ( is_int_type && tok2.is_int() ) 
    {
      if ( tok.is_scalar() ) return Token( (int)(pow( tok.as_int() , tok2.as_int() )) );
      std::vector<int> ans = tok.as_int_vector();
      const int ev = tok2.as_int();
      for (int i=0; i<ans.size(); i++) ans[i] = (int)(pow( ans[i] , ev ) );
      return Token( ans );
    }
  
  double ev = tok2.as_float();

  if ( tok.is_int() || tok.is_float() ) return Token( pow( tok.as_float() , ev ) );
  else if ( tok.is_int_vector() || tok.is_float_vector() )
    {
      std::vector<double> ans = tok.as_float_vector();
      for (int i=0; i<ans.size(); i++) ans[i] = pow( ans[i] , ev ) ;
      return Token( ans );
    }

  return Token();
}





Token TokenFunctions::fn_ifelse( const Token & cond , const Token & opt1 , const Token & opt2 ) const
{
  
  // cond ? opt1 : opt2

  // cond == T --> opt1
  // cond == F --> opt2

  // condition must evaluate to a scalar boolean, or be easily converted
  // no automatic conversion of vector conditions
  
  bool b;
  
  if ( ! cond.is_bool(&b) ) 
    {
      // allow implicit conversion of int -> bool
      if ( cond.is_int() ) b = cond.as_bool();
      else return Token();
    }  

  // opt1 and opt2 must have the same type, or be easily converted

  if ( opt1.type() == opt2.type() ) 
    {
      return b ? opt1 : opt2 ;
    }    
  
  // If opt1 and opt2 are not the same, attempt to convert in some cases

  // Safe automatic upcasts
  // F I --> F F
  // F B --> F F
  // I B --> I I (0/1)
  // do not allow any STRINGs to be converted
  
  Token tmp1 = opt1;
  Token tmp2 = opt2;

  Token::tok_type t1 = tmp1.type();
  Token::tok_type t2 = tmp2.type();

  if ( t1 == Token::UNDEF || t2 == Token::UNDEF ) return Token();
  
  if ( t1 == Token::STRING || t2 == Token::STRING ) 
    Helper::halt("ifelse(?,T,F) cannot specify incompatible return types");
  
  if      ( t1 == Token::FLOAT )  tmp2 = Token( tmp2.as_float() );
  else if ( t2 == Token::FLOAT )  tmp1 = Token( tmp1.as_float() );
  else if ( t1 == Token::INT   )  tmp2 = Token( tmp2.as_int()   ); 
  else if ( t2 == Token::INT   )  tmp1 = Token( tmp1.as_int()   );  
  else return Token();  // all other comparisons not defined

  // Now evaluate the original expression
  return b ? tmp1 : tmp2 ;
  
}


Token TokenFunctions::fn_vec_length( const Token & tok) const 
{
  return Token( tok.size() ) ;
}


Token TokenFunctions::fn_vec_any( const Token & tok1 , const Token & tok2 ) const
{
  return Token( fn_vec_count( tok1 , tok2 ) > 0 );  
}

Token TokenFunctions::fn_vec_any( const Token & tok1 ) const
{
  return Token( fn_vec_count( tok1 , Token( true ) ) > 0 );  
}

Token TokenFunctions::fn_vec_all( const Token & tok1 ) const
{
  return Token( fn_vec_count( tok1 , Token( true ) ) == tok1.size() );  
}

Token TokenFunctions::fn_vec_count( const Token & tok1 , const Token & tok2 ) const
{
  return fn_vec_sum( tok1 == tok2 );
  //  return Token();  
}

Token TokenFunctions::fn_vec_sum( const Token & tok ) const
{
  
  Token::tok_type ttype = tok.type();
  
  if ( tok.is_scalar() ) return tok;
  
  if ( ttype == Token::INT_VECTOR ) 
    {
      int sm = 0;
      std::vector<int> tmp = tok.as_int_vector();
      for (int i=0;i<tmp.size(); i++) sm += tmp[i];      
      return Token( sm );
    }
  
  if ( ttype == Token::FLOAT_VECTOR ) 
    {
      double sm = 0;
      std::vector<double> tmp = tok.as_float_vector();
      for (int i=0;i<tmp.size(); i++) sm += tmp[i];
      return Token( sm );
    }

  // convert to INT
  if ( ttype == Token::BOOL_VECTOR )
    {
      int sm = 0;
      std::vector<bool> tmp = tok.as_bool_vector();
      for (int i=0;i<tmp.size(); i++) sm += tmp[i];
      return Token( sm );
    }

  // undefined
  return Token();
       
}

  
Token TokenFunctions::fn_vec_cat( const std::vector<Token> & tok ) const
{
  if ( tok.size() == 0 ) return Token();
  if ( tok.size() == 1 ) return tok[0];
  // reverse order
  const int ns = tok.size() ;
  Token ans = tok[ns-1];
  for (int i=ns-2;i>=0;i--)
    ans = fn_vec_cat( ans , tok[i] );

  return ans;
}

Token TokenFunctions::fn_vec_cat( const Token & tok1 , const Token & tok2 ) const
{

  Token::tok_type mode ; 
  
  if ( ( tok1.type() == Token::INT || tok1.type() == Token::INT_VECTOR ) )
    {
      if ( ! ( tok2.type() == Token::INT || tok2.type() == Token::INT_VECTOR ) ) 
	Helper::halt( "can only concatenate similar types" );
      else 
	mode = Token::INT_VECTOR;
    }
  
  if ( ( tok1.type() == Token::FLOAT || tok1.type() == Token::FLOAT_VECTOR ) ) 
    {
      if ( ! ( tok2.type() == Token::FLOAT || tok2.type() == Token::FLOAT_VECTOR ) ) 
	Helper::halt( "can only concatenate similar types" );
      else 
	mode = Token::FLOAT_VECTOR;
    }

  if ( ( tok1.type() == Token::STRING || tok1.type() == Token::STRING_VECTOR ) )
    {
      if ( ! ( tok2.type() == Token::STRING || tok2.type() == Token::STRING_VECTOR ) ) 
	Helper::halt( "can only concatenate similar types" );
      else 
	mode = Token::STRING_VECTOR;
    }

  if ( ( tok1.type() == Token::BOOL || tok1.type() == Token::BOOL_VECTOR ) )
    {
      if ( ! ( tok2.type() == Token::BOOL || tok2.type() == Token::BOOL_VECTOR ) ) 
	Helper::halt( "can only concatenate similar types" );
      else
	mode = Token::BOOL_VECTOR;
    }

  if ( mode == Token::INT_VECTOR )
    {
      std::vector<int> res1 = tok1.as_int_vector();
      std::vector<int> res2 = tok2.as_int_vector();
      for (int i=0;i<res2.size();i++) res1.push_back( res2[i] );
      Token tok( res1 );
      return tok;
    }

  if ( mode == Token::FLOAT_VECTOR )
    {
      std::vector<double> res1 = tok1.as_float_vector();
      std::vector<double> res2 = tok2.as_float_vector();
      for (int i=0;i<res2.size();i++) res1.push_back( res2[i] );
      Token tok( res1 );
      return tok;
    }

  if ( mode == Token::STRING_VECTOR )
    {
      std::vector<std::string> res1 = tok1.as_string_vector();
      std::vector<std::string> res2 = tok2.as_string_vector();
      for (int i=0;i<res2.size();i++) res1.push_back( res2[i] );
      Token tok( res1 );
      return tok;
    }

  if ( mode == Token::BOOL_VECTOR )
    {
      std::vector<bool> res1 = tok1.as_bool_vector();
      std::vector<bool> res2 = tok2.as_bool_vector();
      for (int i=0;i<res2.size();i++) res1.push_back( res2[i] );
      Token tok( res1 );
      return tok;
    }

  // undefined
  return Token();

}


Token TokenFunctions::fn_vec_mean( const Token & tok1 ) const
{
  return fn_vec_sum( tok1 ) / fn_vec_length( tok1 ); 
}

Token TokenFunctions::fn_vec_sd( const Token & tok1 ) const
{
  if ( ! ( tok1.is_float_vector() || tok1.is_int_vector() || tok1.is_bool_vector() ) ) return Token();
  std::vector<double> d = tok1.as_float_vector();
  return Token( MiscMath::sdev( d ) );
}


Token TokenFunctions::fn_vec_sort( const Token & tok ) const
{

  if ( ! tok.is_vector() ) return tok;  

  Token::tok_type ttype = tok.type();        

  if ( ttype == Token::INT_VECTOR ) 
    {
      std::vector<int> t = tok.as_int_vector();
      std::sort( t.begin() , t.end() );
      return Token( t );
    }
  else if ( ttype == Token::FLOAT_VECTOR ) 
    {
      std::vector<double> t = tok.as_float_vector();
      std::sort( t.begin() , t.end() );
      return Token( t );
    }
  else if  ( ttype == Token::STRING_VECTOR ) 
    {
      std::vector<std::string> t = tok.as_string_vector();
      std::sort( t.begin() , t.end() );
      return Token( t );  
    }
  else if ( ttype == Token::BOOL_VECTOR ) 
    {
      std::vector<bool> t = tok.as_bool_vector();
      std::sort( t.begin() , t.end() );
      return Token( t );
    }
  return Token();
}


Token TokenFunctions::fn_vec_extract( const Token & tok , const Token & idx ) const
{
  
  if ( ! ( idx.is_int() || idx.is_int_vector() || idx.is_bool_vector() ) ) 
    Helper::halt( "index for vector subscripting is not an integer value, integer vector or boolean vector" );
  
  //
  // NEW, for subsetting, only change the ve[] mask
  //

  // if scalar, return whole thing (i.e. could only have been x[0] so pointless)                                           
  if ( ! tok.is_vector() ) return tok;
  
  // nb, the vector may already be subsetted; size() will return current ve[] if so
  //  i.e. X = int(10,11,12) ;   X[ c(2,3) ][ 1 ] = 11 
  //   nb - not X[][] not allowed yet
  
  const int s = tok.size(); 
  //  std::cout << "fn_vec_extract size = " << s << " " << tok.fullsize() << "\n";
  
  std::vector<int> tokve = tok.get_subset_index();
  std::vector<int> ve;
  
  // subset() does range checking, so no need to repeat it here
  // just build the new subset (w/ 0-based counting)

  if ( idx.is_int() )
    {
      ve.push_back( tokve[ idx.as_int() - 1 ] );
    }
  else if ( idx.is_int_vector() )
    {
      std::vector<int> ii = idx.as_int_vector();

      for (int i=0; i<ii.size(); i++)
	{
	  int ix = ii[i] - 1;
	  if ( ix < 0 || ix >= tokve.size() )
	    Helper::halt( "bad index" );
	  ve.push_back( tokve[ ix ] );
	}
    }
  else if ( idx.is_bool_vector() )
    {
      // here require that length of index matches (no cycling, as per R)
      if ( idx.size() != s )
        Helper::halt( "boolean index vector should be of similar size to matching vector "
                      + tok.name() + " " + Helper::int2str( idx.size() ) );
      
      std::vector<bool> bi = idx.as_bool_vector();
      for (int i=0; i<s; i++)
	if ( bi[i] ) ve.push_back( tokve[i] );
      
    }


  Token res = tok;

  res.subset( ve );

  // std::cout <<  " new tok " << tok.size() << " " << tok.fullsize() << "\n";
  // std::cout <<  " new res " << res.size() << " " << res.fullsize() << "\n";
  
  // all done

  return res;
  
  
  //
  // OLD ... can delete 
  //

  // // subscript out of range? (1-based)
  // if ( i < 1 || i > s ) 
  // 	Helper::halt( "out of range for "
  //                 + tok.name()
  //                 + " (" + Helper::int2str(i)
  //                 + " of "
  //                 + Helper::int2str( s ) +")" );
  
	
  if ( 0 )
    {
      // extract out a single element and return as a scalar
  
  if ( idx.is_int() )
    {
      
      int i = idx.as_int();
      
      // subscript out of range? (1-based)
      if ( i < 1 || i > tok.size() ) 
	Helper::halt( "out of range for "
		      + tok.name()
		      + " (" + Helper::int2str(i)
		      + " of "
		      + Helper::int2str(tok.size() ) +")" );
      //return Token();
      
      // if scalar, return whole thing (i.e. could only have been x[0] so pointless)
      if ( ! tok.is_vector() ) return tok;
      
      Token::tok_type ttype = tok.type();      
      
      if      ( ttype == Token::INT_VECTOR )    return Token( tok.int_element(i-1) ); 
      else if ( ttype == Token::FLOAT_VECTOR )  return Token( tok.float_element(i-1) );
      else if ( ttype == Token::STRING_VECTOR ) return Token( tok.string_element(i-1) );
      else if ( ttype == Token::BOOL_VECTOR )   return Token( tok.bool_element(i-1) );
      
      return Token();

    }

  // extract out a list 
  
  if ( idx.is_int_vector() )
    {
      
      Token::tok_type ttype = tok.type();
      
      if ( ttype == Token::INT_VECTOR ) 
	{
	  std::vector<int> ans;
	  for (int i=0;i<idx.size();i++) ans.push_back( tok.int_element( idx.int_element(i) - 1 ) );
	  return Token(ans);
	}
      else if ( ttype == Token::FLOAT_VECTOR )
	{
	  std::vector<double> ans;
	  for (int i=0;i<idx.size();i++) ans.push_back( tok.float_element( idx.int_element(i) - 1 ) );
	  return Token(ans);
	}
      else if ( ttype == Token::STRING_VECTOR )
	{
	  std::vector<std::string> ans;
	  for (int i=0;i<idx.size();i++) ans.push_back( tok.string_element( idx.int_element(i) - 1 ) );
	  return Token(ans);
	}
      else if ( ttype == Token::BOOL_VECTOR )
	{
	  std::vector<bool> ans;
	  for (int i=0;i<idx.size();i++) ans.push_back( tok.bool_element( idx.int_element(i) - 1 ) );
	  return Token(ans);
	}
      return Token();
    }

  
  // extract from a boolean vector
  
  if ( idx.is_bool_vector() )
    {
  
      // here require that length of index matches (no cycling, as per R)
      if ( idx.size() != tok.size() ) 
	Helper::halt( "boolean index vector should be of similar size to matching vector " 
		      + tok.name() + " " + Helper::int2str( idx.size() ) );
      
      Token::tok_type ttype = tok.type();
      
      if ( ttype == Token::INT_VECTOR ) 
	{
	  std::vector<int> ans;
	  for (int i=0;i<idx.size();i++) if ( idx.bool_element(i) ) ans.push_back( tok.int_element(i) );
	  return Token(ans);
	}
      else if ( ttype == Token::FLOAT_VECTOR )
	{
	  std::vector<double> ans;
	  for (int i=0;i<idx.size();i++) if ( idx.bool_element(i) ) ans.push_back( tok.float_element(i) );
	  return Token(ans);
	}
      else if ( ttype == Token::STRING_VECTOR )
	{
	  std::vector<std::string> ans;
	  for (int i=0;i<idx.size();i++) if ( idx.bool_element(i) ) ans.push_back( tok.string_element(i) );
	  return Token(ans);
	}
      else if ( ttype == Token::BOOL_VECTOR )
	{
	  
	  std::vector<bool> ans;
	  for (int i=0;i<idx.size();i++) if ( idx.bool_element(i) ) ans.push_back( tok.bool_element(i) );
	  return Token(ans);
	}
      return Token();
    
    }

    }

  
  return Token();
}


Token TokenFunctions::fn_vec_min( const Token & tok ) const
{
  if ( ! tok.is_vector() ) return tok;
  Token s = fn_vec_sort(tok);
  Token::tok_type ttype = tok.type();  
  if      ( ttype == Token::INT_VECTOR ) return Token( s.int_element(0) );
  else if ( ttype == Token::FLOAT_VECTOR ) return Token( s.float_element(0) );
  else if ( ttype == Token::BOOL_VECTOR ) return Token( s.bool_element(0) );
  else if ( ttype == Token::STRING_VECTOR ) return Token( s.string_element(0) );
  return Token();
}


Token TokenFunctions::fn_vec_maj( const Token & tok ) const
{
  if ( ! tok.is_vector() ) return tok;
  Token s = fn_vec_sort(tok);
  Token::tok_type ttype = tok.type();  
  const int sz = tok.size() - 1;
  if      ( ttype == Token::INT_VECTOR ) return Token( s.int_element(sz) );
  else if ( ttype == Token::FLOAT_VECTOR ) return Token( s.float_element(sz) );
  else if ( ttype == Token::BOOL_VECTOR ) return Token( s.bool_element(sz) );
  else if ( ttype == Token::STRING_VECTOR ) return Token( s.string_element(sz) );
  return Token();
}


Token TokenFunctions::fn_vec_new_float( const std::vector<Token> & tok ) const
{
  // read off in reverse (RPN) 
  // concatenates scalars and vectors
  if ( tok.size() == 0 ) return Token();
  std::vector<double> d;
  for (int i=tok.size()-1; i >= 0 ; i-- ) 
    for (int j=0;j<tok[i].size();j++) 
      d.push_back( tok[i].as_float_element(j) );
  return Token( d );
}


Token TokenFunctions::fn_vec_new_int( const std::vector<Token> & tok ) const
{
  if ( tok.size() == 0 ) return Token();
  std::vector<int> d;
  for (int i=tok.size()-1; i >= 0 ; i-- ) 
    for (int j=0;j<tok[i].size();j++) 
      d.push_back( tok[i].as_int_element(j) );      
  return Token( d );
}

Token TokenFunctions::fn_vec_new_str( const std::vector<Token> & tok ) const
{
  if ( tok.size() == 0 ) return Token();
  std::vector<std::string> d;
  for (int i=tok.size()-1; i >= 0 ; i-- ) 
    for (int j=0;j<tok[i].size();j++) 
      d.push_back( tok[i].as_string_element(j) );
  return Token( d );
}


Token TokenFunctions::fn_vec_new_bool( const std::vector<Token> & tok ) const
{
  if ( tok.size() == 0 ) return Token();
  std::vector<bool> d;
  for (int i=tok.size()-1; i >= 0 ; i-- ) 
    for (int j=0;j<tok[i].size();j++) 
      d.push_back( tok[i].as_bool_element(j) );
  return Token( d );
}


Token TokenFunctions::fn_assign( Token & lhs , const Token & rhs0 )
{
  
  // std::cout << "orig value val fs " << lhs.fullsize() << "\n";
  // std::cout << "orig value val " << lhs << "\n";

  
  // assignment is either to the 'local' or 'global' (e.g. accumulator) instance
  
  // global variables, if being used, need to start with an underscore to distinguish them
  // also, the variable must already exist (i.e. have been initialized with EVAL globals=X,Y,Z for example

  // enfore rule that you cannot assign to existing variables; this means if we have a '.' in the name
  // i.e. it is a member of an instance that already exists, then we cannot assign
  // give an error

  if ( lhs.name().find( "." ) != std::string::npos ) 
    Helper::halt( "cannot assign new values to existing instance meta-data:" + lhs.name() );

  bool global = accumulator != NULL && global_vars->find( lhs.name() ) != global_vars->end() ;
  
  instance_t * m = global ? accumulator : meta ; 
  
  if ( ! m ) return Token();
  
  // still check that we already found a global var, i.e. is initialized
  if ( global && m->find( lhs.name() ) == NULL ) 
    Helper::halt( "internal error: did not initialize global variable " + lhs.name() ); 

  //
  // prune the RHS if needed
  //
  
  Token rhs = rhs0;
  
  if ( rhs.is_masked() )
    rhs.prune();


  //
  // Permissable assignments
  //
  
  // scalar <- scalar
  // scalar <- vector  :: will change scalar to vector
  // unmasked vector <- vector   : must be same size, but can change type of LHS (i.e. to RHS type)
  
  // LHS type and size stays the same
  // masked vector <- vector     : fills in masked values with vector of matching size; will cast to LHS type 
  // (masked) vector <- scalar   : sets all values to scalar: will attempt type conversion to LHS
  
  const bool lhs_is_vec = lhs.is_vector();
  const bool lhs_is_masked = lhs.is_masked();
  const bool rhs_is_vec = rhs.is_vector();  
  
  if ( lhs_is_vec && rhs_is_vec && lhs.size() != rhs.size() )
    Helper::halt( "assigning vector of different length not valid" );

  //
  // scalar <- scalar 
  //
  
  if ( ! ( lhs_is_vec || rhs_is_vec ) )
    {

      // changes type of LHS to match RHS

      bool b;
      int i;
      double f;
      std::string s;
      
      if ( rhs.is_float(&f) ) 
	{
	  m->set( lhs.name() , f );
	  lhs.set( f );
	}
      else if ( rhs.is_int(&i) )
	{      
	  m->set( lhs.name() , i );
	  lhs.set( i );	  
	}
      else if ( rhs.is_bool(&b) )
	{      
	  lhs.set( b );
	  m->set( lhs.name() , b );	  
	}
      else if ( rhs.is_string(&s) )
	{
	  m->set( lhs.name() , s );
	  lhs.set( s );
	}
      else return Token() ;
      
      return Token( true );
      
    }

  //
  // unmasked vector <- vector (may change type of LHS)
  // or
  // scalar/undefined <- vector (will make LHS a vector)
  //
  
  if ( rhs_is_vec && ( ( lhs_is_vec && ! lhs_is_masked ) || ! lhs_is_vec ) )
    {
      
      // changes type of LHS to match RHS
      // but the size must be the same, i.e not allowing size to change
      
      std::vector<double> fv;
      std::vector<bool> bv;
      std::vector<int> iv;
      std::vector<std::string> sv;
 
      if ( rhs.is_float_vector(&fv) ) 
	{
	  m->set( lhs.name() , fv );
	  lhs.set( fv );	  
	}
      else if ( rhs.is_int_vector(&iv) )
	{                  
	  m->set( lhs.name() , iv );
	  lhs.set( iv );
	}
      else if ( rhs.is_bool_vector(&bv) )
	{      
	  m->set( lhs.name() , bv );
	  lhs.set( bv );	  
	}
      else if ( rhs.is_string_vector(&sv) )
	{
	  m->set( lhs.name() , sv );
	  lhs.set( sv );	  
	}
      else return Token();

      return Token( true );
    }

  //
  // (masked) vector <- scalar 
  //

  if ( ! rhs_is_vec )
    {

      // type and size of LHS will not change

      if ( lhs.is_float_vector() )
	{
	  std::vector<double> ff( lhs.size() , rhs.as_float() );
	  // do update 
	  lhs.update( ff );
	  lhs.unmask();
	  m->set( lhs.name() , lhs.as_float_vector() );
	  return Token( true );
	}
      
      else if ( lhs.is_int_vector() )
	{
	  std::vector<int> ii( lhs.size() , rhs.as_int() );      
	  // do update 
	  lhs.update( ii );
	  lhs.unmask();
	  m->set( lhs.name() , lhs.as_int_vector() );
	  return Token( true );
	}
      
      else if ( lhs.is_bool_vector() )
        {
          std::vector<bool> bb( lhs.size() , rhs.as_bool() );
          // do update   
          lhs.update( bb );
          lhs.unmask();
          m->set( lhs.name() , lhs.as_bool_vector() );
          return Token( true );
        }
      
      else if ( lhs.is_string_vector() )
        {
          std::vector<std::string> ss( lhs.size() , rhs.as_string() );
          // do update
          lhs.update( ss );
          lhs.unmask();
          m->set( lhs.name() , lhs.as_string_vector() );
          return Token( true );
	}
      
      return Token( true );  
    }

  
  //
  // masked vector <- vector 
  //
  
  if ( lhs.is_float_vector() )
    {      
      std::vector<double> ff = rhs.as_float_vector();
      // do update
      lhs.update( ff );
      lhs.unmask();
      m->set( lhs.name() , lhs.as_float_vector() );
      return Token( true );
    }
  
  else if ( lhs.is_int_vector() )
    {
      std::vector<int> ii = rhs.as_int_vector();
      // do update
      lhs.update( ii );
      lhs.unmask();
      m->set( lhs.name() , lhs.as_int_vector() );
      return Token( true );
    }

    else if ( lhs.is_bool_vector() )
    {
      std::vector<bool> bb = rhs.as_bool_vector();
      // do update
      lhs.update( bb );
      lhs.unmask();
      m->set( lhs.name() , lhs.as_bool_vector() );
      return Token( true );
    }

    else if ( lhs.is_string_vector() )
    {
      std::vector<std::string> ss = rhs.as_string_vector();
      // do update
      lhs.update( ss );
      lhs.unmask();
      m->set( lhs.name() , lhs.as_string_vector() );
      return Token( true );
    }
  
  return Token();
}

  
