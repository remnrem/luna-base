
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


#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "token-eval.h"
#include "annot/annot.h"

// initialise static members
std::map<std::string,Token::tok_type> Token::tok_map; 
std::map<Token::tok_type,std::string> Token::tok_unmap; 
std::map<std::string,int> Token::fn_map; 


void Eval::init( bool na )
{
  is_valid = false;
  
  no_assignments = na;

  errs = "";
}

bool Eval::get_token( std::string & input ,  Token & tok )
{

  // no more input 

  if ( input.size() == 0 ) return false;

  std::string c = input.substr(0,1);
  
  // To interpret whether a + or - sign, or '.' is the start of a number 
  // rather than a binary operator
  // e.g 
  // (C - 2)  vs ( C < -2 )
  
  // eat leading space
  
  while ( 1 ) 
    {
      if ( c == " " ) input = input.substr(1);
      else break;
      if ( input.size() == 0 ) return false;	    
      c = input.substr(0,1);
    }
  


  // token is either func()
  //                 number
  //                 variable
  //                 operator 
  

  std::map<std::string,Token::tok_type>::iterator i = 
    Token::tok_map.find( c ) ;

  // Numeric, allow to start with + or - or .
  
  if ( ( c >= "0" && c <= "9" ) || c == "." || ( (!previous_value) && ( c=="-" || c=="+" ) ) )
    {

      unsigned int p = 1;
      
      bool need_leading_zero = c == ".";
      
      while ( 1 ) 
	{
	  
	  if ( p == input.size() ) break;
	  char d = input[p];
	  if ( d >= '0' && d <= '9' ) c += input.substr(p,1);
	  else if ( d == '.' ) c += input.substr(p,1); 
	  else if ( d == 'e' || d == 'E' ) c += input.substr(p,1);
	  else if ( d == '+' || d == '-' ) {
	    // + and - okay, if previous char was 'e' or 'E'
	    if ( input.substr(p-1,1) == "E" || input.substr(p-1,1) == "e" )
	      c += input.substr(p,1);
	    else break;
	  }
	  else break;
	  ++p;
	}

      // a valid int or float? 
      int i;
      double d;
      bool err = false;

      if ( ! Helper::from_string<int>( i , need_leading_zero ? "0"+c : c , std::dec ) ) err = true;
      if ( ! Helper::from_string<double>( d , need_leading_zero ? "0"+c : c , std::dec ) ) err = true;

      if ( ! err ) 
	{
	  if ( i == d ) tok.set(i);
	  else tok.set(d);
	  previous_value = true;
	}
      
    }
    
  else if ( i != Token::tok_map.end() )
    {

      // looks like an operator
      // if potentially a two-char operator, look for 2nd char
      
      if      ( c == "%" && input.substr(1,1) == "%" ) c = "%%";
      else if ( c == "<" && input.substr(1,1) == "=" ) c = "<=";
      else if ( c == ">" && input.substr(1,1) == "=" ) c = ">=";
      else if ( c == "&" && input.substr(1,1) == "&" ) c = "&&";
      else if ( c == "|" && input.substr(1,1) == "|" ) c = "||";
      else if ( c == "!" && input.substr(1,1) == "=" ) c = "!=";
      else if ( c == "=" && input.substr(1,1) == "=" ) c = "==";
      else if ( c == "=" && input.substr(1,1) == "~" ) c = "=~";
      
      tok.oper( Token::tok_map[ c ] ); 
      
      previous_value = false;

      // eat input name 
      input = input.substr( c.size() );
      return true;
      
    }
  

  // If not an operator or a function, must be 
  // a variable (meta-field of a literal)
  // literal strings should be in 'single quotes'
  
  
  // Handle some other special case syntax first: ( ) and ,
  
  if ( c == "(" ) 
    {
      tok.oper( Token::LEFT_PARENTHESIS );
      previous_value = false;
    }
  
  else if ( c == ")" )
    {
      tok.oper( Token::RIGHT_PARENTHESIS );
      previous_value = true;
    }
  
  else if ( c == "," )
    {
      tok.oper( Token::ARG_SEPARATOR );
      previous_value = false;      
    }

  // a literal string
  else if ( c == "'" )
    {
      unsigned int p = 1;
      bool found = false;
      while ( 1 ) 
	{
	  if ( p == input.size() ) break;
	  if ( input.substr(p,1) == "'" ) { c+= "'"; found = true; break; }
	  c += input.substr(p,1);
	  ++p;
	}
      
      tok.set( c.substr(1,c.size()-2) ); 
      previous_value = true;
    }
  
  // a literal string with specific start/stop quotes, allowing nesting
  else if ( c == "{" )
    {
      int cnt = 1;
      unsigned int p = 1;
      bool found = false;
      while ( 1 ) 
	{
	  
	  if      ( p == input.size() ) break;
	  else if ( input.substr(p,1) == "{" ) ++cnt;
	  else if ( input.substr(p,1) == "}" ) 
	    { 
	      --cnt;
	      if ( cnt == 0 ) 
		{
		  c+= "}"; 
		  found = true; break; 
		}
	    }
	  
	  c += input.substr(p,1);
	  ++p;
	}      

      tok.set( c.substr(1,c.size()-2) ); 

      previous_value = true;
    }

  // otherwise, assume a variable that start: a-z, A-Z or _
  // or a function name if ends in "("
  
  else if ( ( c >= "a" && c <= "z" ) || 
	    ( c >= "A" && c <= "Z" ) ||
	    c == "_" ) 
    {
            
      // Is this a function?
      bool isfn = false;

      // read until space or next operator char, or comma

      unsigned int p = 1;
      while ( 1 ) 
	{
	  if ( p == input.size() ) break;
	  
	  if ( input.substr(p,1) == "(" ) 
	    { isfn=true; break; }
	  
	  // e.g. end of set(AB) for 'AB'	  
	  if ( input.substr(p,1) == ")" ) break; 
	  
	  if ( input.substr(p,1) == "," ) break;

	  if ( input.substr(p,1) == " " ) break;

	  if ( Token::tok_map.find( input.substr(p,1) ) 
	       != Token::tok_map.end() ) 
	    break;
	  
	  c += input.substr(p,1);
	  ++p;
	}
      
      if ( isfn )
	{
	
	  // does this look like a valid function name?
	  std::map<std::string,int>::iterator f = 
	    Token::fn_map.find( c );
	  
	  if ( f == Token::fn_map.end() ) 
	    Helper::halt( "did not recognize function " + c + "()" );
	  
	  // store function token
	  tok.function( c );	    
	}
      else if ( c == "true" )
	{
	  tok.set( true );
	}
      else if ( c == "false" )
	{
	  tok.set( false );
	}
      else // add token as variable
	{	  
	  tok.variable( c );
	}
      previous_value = true;
    }
  
  // remove this input

  input = input.substr( c.size() );

  return true;
  
}


int Eval::op_preced( const Token & c )
{
  
  Token::tok_type t = c.type();
  
  switch(t) 
    {
    case Token::NOT_OPERATOR :  return 9;
      
    case Token::MULTIPLY_OPERATOR :  
    case Token::DIVIDE_OPERATOR :  
    case Token::MOD_OPERATOR :  return 8;
      
    case Token::ADD_OPERATOR :
    case Token::SUBTRACT_OPERATOR : return 7;
      
    case Token::LESS_THAN_OPERATOR : 
    case Token::LESS_THAN_OR_EQUAL_OPERATOR : 
    case Token::GREATER_THAN_OPERATOR : 
    case Token::GREATER_THAN_OR_EQUAL_OPERATOR : return 6;
      
    case Token::EQUAL_OPERATOR : 
    case Token::HAS_OPERATOR :
    case Token::UNEQUAL_OPERATOR : return 5;

    case Token::AND_OPERATOR : return 4;

    case Token::OR_OPERATOR : return 3;
      
    case Token::LEFT_PARENTHESIS : 
    case Token::RIGHT_PARENTHESIS : return 2;
      
    case Token::ASSIGNMENT_OPERATOR : return 1;
      
    case Token::ARG_SEPARATOR : return 0;

    default:
      break;
      
    }
  return 0;
}

bool Eval::op_left_assoc(const Token & tok )
{
  
  Token::tok_type t = tok.type();
  
  switch(t)    
    {
    case Token::MULTIPLY_OPERATOR :
    case Token::DIVIDE_OPERATOR :
    case Token::MOD_OPERATOR :
    case Token::ADD_OPERATOR :
    case Token::EQUAL_OPERATOR :
    case Token::HAS_OPERATOR :
    case Token::AND_OPERATOR :
    case Token::OR_OPERATOR :
    case Token::ARG_SEPARATOR : 
    case Token::SUBTRACT_OPERATOR : return true;
      
    case Token::ASSIGNMENT_OPERATOR :
    case Token::NOT_OPERATOR : return false;

    default:
      break;
    }
  return false;
}

unsigned int Eval::op_arg_count( const Token & tok )
{

  Token::tok_type t = tok.type();

  //  std::cout << "tok type " << tok.name() << " " << t << "\n";

  switch( t )  
    {
      
    case Token::NOT_OPERATOR :  return 1;
      
    case Token::ASSIGNMENT_OPERATOR :       
    case Token::MULTIPLY_OPERATOR :  
    case Token::DIVIDE_OPERATOR :  
    case Token::MOD_OPERATOR : 
    case Token::ADD_OPERATOR :
    case Token::SUBTRACT_OPERATOR : 
    case Token::LESS_THAN_OPERATOR : 
    case Token::LESS_THAN_OR_EQUAL_OPERATOR : 
    case Token::GREATER_THAN_OPERATOR : 
    case Token::GREATER_THAN_OR_EQUAL_OPERATOR : 
    case Token::EQUAL_OPERATOR : 
    case Token::HAS_OPERATOR : 
    case Token::UNEQUAL_OPERATOR : 
    case Token::AND_OPERATOR : 
    case Token::OR_OPERATOR : return 2;
            
    case Token::FUNCTION : 
      if ( Token::fn_map.find( tok.name() ) == Token::fn_map.end() )
	Helper::halt( "did not recognize function " + tok.name() );
      return Token::fn_map[ tok.name() ];
      
    default:
      Helper::halt( "did not recognize operator " + tok.name() );
      break;
      
    }

  return 0; 
}


bool Eval::shunting_yard( const std::string & oinput, 
			  std::vector<Token> & output )
{
  
  std::string input = oinput;

  // collect output tokens here
  output.resize( input.size() );
  output.clear();
  
  Token sc;
  
  std::vector<Token> stack;
  
  // redundant, but use 'sl' to keep track of stack size also
  
  unsigned int sl = 0;
  
  previous_value = false;

  // Parse input string token by token
  
  while( 1 )
    {
      Token c;
      
      // get next token -- if end of input break
      
      if ( ! get_token( input , c ) ) break;
      
      // If the token is a number (identifier), then add it to output queue.
      
      if( c.is_ident() )  
	{
	  output.push_back(c);
	}
      
      // If the token is a function token, then push it onto the stack.
      
      else if ( c.is_function() )   
	{
	  stack.push_back( c );
	  ++sl;
	}
	
      // If the token is a function argument separator (e.g., a comma):
      
      else if ( c.is_separator() )
	{
	  bool pe = false;
	  
	  while( sl > 0 )   	    
	    {	
	      sc = stack.back();

	      if( sc.is_left_paren() )
		{		    
		  pe = true;
		  break;
		}
	      else  
		{
		  // Until the token at the top of the stack is a left
		  // parenthesis, pop operators off the stack onto the
		  // output queue.
		  
		  output.push_back( sc );
		  stack.pop_back();
		  sl--;
		}
	    }
	  
	  // If no left parentheses are encountered, either the
	  // separator was misplaced
	  // or parentheses were mismatched.
	  
	  if( !pe )   
	    {
	      Helper::halt( "separator or parentheses mismatched" );
	      return false;
	    }
	}
      

      // If the token is an operator, op1, then:
      
      else if ( c.is_operator() )  
	{
	  
	  while ( sl > 0 )    
	    {
	      
	      sc = stack.back();
	      
	      // While there is an operator token, o2, at the top of
	      // the stack op1 is left-associative and its precedence
	      // is less than or equal to that of op2, or op1 is
	      // right-associative and its precedence is less than
	      // that of op2,
	      
	      if( sc.is_operator() &&
		  ( ( op_left_assoc(c) && ( op_preced(c) <= op_preced(sc) ) ) ||
		    ( !op_left_assoc(c) && ( op_preced(c) < op_preced(sc) ) )))   {
		
		// Pop o2 off the stack, onto the output queue;
		
		output.push_back( sc );
		stack.pop_back();
		sl--;
	      }
	      else
		{
		  break;
		}
	    }
	  
	  // push op1 onto the stack.
	  
	  stack.push_back( c );
	  ++sl;
	}
      
      
      // If the token is a left parenthesis, then push it onto the stack.
      
      else if ( c.is_left_paren() ) 
	{	    
	  stack.push_back( c );
	  ++sl;
	}
      
      
      // If the token is a right parenthesis:
      
	else if ( c.is_right_paren() )
	  {
	    
	    bool pe = false;
	    
	    // Until the token at the top of the stack is a left
	    // parenthesis, pop operators off the stack onto the
	    // output queue
	    
	    while( sl > 0 )     
	      {
		sc = stack.back();
		
		if( sc.is_left_paren() )
		  {
		    pe = true;
		    break;
		  }
		else  
		  {
		    output.push_back( sc );
		    stack.pop_back();
		    sl--;
		  }
	      }
	    
	    // If the stack runs out without finding a left
	    // parenthesis, then there are mismatched parentheses.
	    
	    if( ! pe )  
	      {
		Helper::halt( "parentheses mismatched" );
		return false;
	      }
	    
	    
	    // Pop the left parenthesis from the stack, but not onto
	    // the output queue.
	    
	    stack.pop_back();
	    sl--;
	    

	    // If the token at the top of the stack is a function
	    // token, pop it onto the output queue.

	    if( sl > 0 )   
	      {
		sc = stack.back();
		
		if( sc.is_function() )
		  {
		    output.push_back( sc );
		    stack.pop_back();
		    sl--;
		  }
	    }
	  }
	else  
	  {	    
	    Helper::halt( "unknown token" ) ;
	    return false; 
	  }
      
    }

  
  
  // When there are no more tokens to read: While there are still
  // operator tokens in the stack:
  
  while ( sl > 0 )  
    {
      sc = stack.back();
      
      if ( sc.is_left_paren() || sc.is_right_paren() )
	{
	  Helper::halt( "parentheses mismatched" );
	  return false;
        }
      
      output.push_back( sc );
      stack.pop_back();
      --sl;
    }
  
  return true;
}


bool Eval::execute( const std::vector<Token> & input )
{

  Token sc;
  
  std::vector<Token> stack;
  
  // redundant, but keep track of stack size here also
  unsigned int sl = 0;
  
  if ( verbose ) 
    {
      std::cout << "-----------------------------------------------------------\n";
      std::cout << "evaluating " << input.size() << " tokens\n";  
      for (int i=0;i<input.size();i++)
	std::cout << " token " << i << "\t" << input[i] << "\n";
      std::cout << "\n";
    }

  // While there are input tokens left
  
  for (unsigned int i = 0 ; i < input.size() ; i++ )
    {
      
      // Read the next token from input.
      
      Token c = input[i];

      if ( verbose ) 
	std::cout << " considering token " << c << "\n";
      
      // If the token is a value or identifier
      
      if ( c.is_ident() )
	{
	  if ( verbose ) 
	    std::cout << "  pushing onto stack\n";
	  
	  // Push it onto the stack.	  
	  stack.push_back( c );
	  ++sl;
        }
      
      // Otherwise, the token is an operator
      // (operator here includes both operators, and functions).
      
      else if ( c.is_operator() || c.is_function() )
	{
	  
	  // It is known a priori that the operator takes n arguments.
	  
	  int nargs = op_arg_count(c);
	  
	  // If there are fewer than n values on the stack
	  
	  if ( (int)sl < nargs && nargs != -1 ) 
	    {
	      Helper::halt( "not enough arguments for " + c.name() ) ;
	      return false;
	    }
	  
	  //
	  // Else, Pop the top n values from the stack.
	  // Evaluate the operator, with the values as arguments.
	  //
	  
	  Token res;

	  if ( c.is_function() ) 
	    {
	  
	      std::vector<Token> args;
	      
	      // for fixed 'nargs' functions, just count them off 
	      
	      if ( nargs == -1 ) 
		{
		  nargs = stack.back().as_int();
		  stack.pop_back(); 
		  sl--;
		}
	      
	      while ( nargs > 0 ) 
		{		  
		  sc = stack.back();
		  stack.pop_back(); 
		  sl--;
		  args.push_back(sc);
		  --nargs;

		  if ( verbose )
		    std::cout << "  popping argument off stack: " << sc << "\n";
		}
	      

	      //
	      // Perform function calls:
	      //
	      
	      // correct # of arguments? 
	      // -1 means variable length
	      
	      const int nargs = Token::fn_map[ c.name() ]; 
	      if ( args.size() != nargs && nargs != -1 ) 
		{
		  Helper::halt( "wrong number of arguments for " + c.name() );
		  return false;
		}


	      // note: args are in reverse order here
	      
	      if      ( c.name() == "if" )     res = func.fn_set( args[0] );
	      else if ( c.name() == "ifnot" )  res = func.fn_notset( args[0] );
	      
	      else if ( c.name() == "sqrt" )   res = func.fn_sqrt( args[0] );
	      else if ( c.name() == "sqr" )    res = func.fn_sqr( args[0] );
	      else if ( c.name() == "pow" )    res = func.fn_pow( args[1] , args[0] );

	      else if ( c.name() == "rnd" )    res = func.fn_rnd();
	      else if ( c.name() == "rand" )   res = func.fn_rnd( args[0] );

	      else if ( c.name() == "exp" )    res = func.fn_exp( args[0] );
	      else if ( c.name() == "log" )    res = func.fn_log( args[0] );
	      else if ( c.name() == "log10" )  res = func.fn_log10( args[0] );
	      
	      else if ( c.name() == "floor" )  res = func.fn_floor( args[0] );
	      else if ( c.name() == "round" )  res = func.fn_round( args[0] );
				      
	      else if ( c.name() == "ifelse" ) res = func.fn_ifelse( args[2], args[1], args[0] );
	      
	      // vector functions
	
	      else if ( c.name() == "element" )  res = func.fn_vec_extract( args[1] , args[0] );
	      
	      else if ( c.name() == "length" )   res = func.fn_vec_length( args[0] );	      
	      else if ( c.name() == "size" )     res = func.fn_vec_length( args[0] );	      
	      	      
	      else if ( c.name() == "min" )      res = func.fn_vec_min( args[0] );	      
	      else if ( c.name() == "max" )      res = func.fn_vec_maj( args[0] );
		
	      else if ( c.name() == "sum" )      res = func.fn_vec_sum( args[0] );
	      else if ( c.name() == "mean" )     res = func.fn_vec_mean( args[0] );
	      else if ( c.name() == "sd" )       res = func.fn_vec_sd( args[0] );
	      else if ( c.name() == "sort" )     res = func.fn_vec_sort( args[0] );

	      else if ( c.name() == "num_func" )  res = func.fn_vec_new_float( args );
	      else if ( c.name() == "int_func" )  res = func.fn_vec_new_int( args );
	      else if ( c.name() == "txt_func" )  res = func.fn_vec_new_str( args );     
	      else if ( c.name() == "bool_func" ) res = func.fn_vec_new_bool( args ); 
	      else if ( c.name() == "c_func" )    res = func.fn_vec_cat( args );
	      
	      else if ( c.name() == "any" )       res = func.fn_vec_any( args[0] );	      
	      else if ( c.name() == "all" )       res = func.fn_vec_all( args[0] );	      
	      else if ( c.name() == "contains" )  res = func.fn_vec_any( args[1] , args[0] );
	      else if ( c.name() == "countif" )   res = func.fn_vec_count( args[1] , args[0] );
	      	    	      
	      else Helper::halt( "did not recognize function " + c.name() );

	      
	    }
	  else
	    {
	      
	      if ( nargs == 1 )  // unary operator
		{

		  sc = stack.back();
		  stack.pop_back(); sl--; 
		  res = c.operands( sc );

		  if ( verbose )
		    {
		      std::cout << "  popping 1 value off stack, " << sc << "\n";		      
		    }
		  
		}
	      else // binary operator
		{
		  
		  Token t0 = stack.back();
		  stack.pop_back(); sl--;
		  
		  sc = stack.back();
		  stack.pop_back(); sl--;		    
		  
		  if ( verbose )
		    {
		      std::cout << "  popping 2 values off stack, " << t0 << " and " << sc << "\n";
		    }

		  // For assignment (that impacts meta-information)
		  // need to call a function also

		  if ( c.is_assignment() ) 
		    {
		      
		      if ( no_assignments ) 
			Helper::halt( "no A = B assigments allowed in this expression" );
		      
		      res = func.fn_assign( sc, t0 );  // res==T

		      // and bind new value if needed
		      bind( &sc );

		    }		  
		  else
		    res = c.operands( t0 , sc );
		  		  
		}
	    }
	  
	  if ( verbose ) 
	    std::cout << "  pushing result on stack, " << res << "\n";

	  // Push the returned results, if any, back onto the stack.
	  stack.push_back( res );
	  ++sl;
        }


      if ( verbose )
	{
	  std::cout << " current stack size n=" << stack.size() << "\n";
	  for (int ss = 0 ; ss < stack.size() ; ss++ )
	    std::cout << "  " << stack[ss] << "\n";
	  std::cout << "\n\n";
	}

    }


  // If there is only one value in the stack
  // That value is the result of the calculation.
  
  if ( sl != 1 || stack.size() != 1 ) 
    {
      Helper::halt( "badly formed eval expression" );
      return false;
    }
  
  if ( sl == 1 ) 
    {

      sc = stack.back(); 
      stack.pop_back();
      sl--;

      if ( verbose )
	{
	  std::cout << "final value " << sc << "\n";
	  std::cout << "ALL DONE.\n\n";
	}
      
      // store result in primary slot, e
      e = sc;
      
      return true;
    }
  
  
  // If there are more values in the stack
  // (Error) The user input has too many values.
  
  Helper::halt( "badly formed expression: too many values" );
  return false;
}

bool Eval::parse( const std::string & input )
{

  delete_symbols();    
  
  std::string input2 = input;

  if ( ! expand_indices( &input2 ) ) return false;

  if ( ! expand_vargs( &input2 ) ) return false;
  
  // this may contain several statements, delimited by ";"
  // evaluate each sequential, to perform any assignments into
  // meta data, but only the last will be reflected in the 'e' 
  // endpoint

  std::vector<std::string> etok0 = Helper::parse( input2, ";" );
  
  // strip out any empty expressions (e.g. if there was a final ; but no expression after)
  
   std::vector<std::string> etok;
   for (int i=0;i<etok0.size();i++)
     {
       std::string j = Helper::ltrim( etok0[i] );
       j = Helper::rtrim( j );
       if ( j.size() > 0 ) etok.push_back( j );
     }
  
  // set number of evals we need to do

  neval = etok.size();
  
  output.resize( neval );
  
  is_valid = true;


  for (int i=0; i<etok.size(); i++)
    {      
      output[i].clear();    

      errs = "";

      if ( ! shunting_yard( etok[i], output[i] ) )
	is_valid = false;
      
    }

  // set pointers to all variables now construction of tokens is complete
  for (int i=0; i<etok.size(); i++)  
    locate_symbols( output[i] ); 
  
  return is_valid;

  
}

void Eval::assign_to( instance_t & m )
{
  func.attach( &m );
}

void Eval::bind( const std::map<std::string,annot_map_t> & inputs , 
		 instance_t * outputs , 
		 instance_t * accumulator , 
		 const std::set<std::string> * global_vars , 
		 bool reset )   
{

  
  if ( reset ) reset_symbols();  
  
  //
  // Input: bind information from inputs and also from accumulator 
  //
  
  //
  // Output: either to outputs (i.e. which is assumed to be local) OR to the 'global' accumulator (i.e. 
  // which is assumed to persist across evaluations
  //
  

  //
  // create a single instance that organizes all the information in
  // the annot_map_t above in a sane variable naming scheme
  // 
  // nb. some redundancy here, as we copy stuff over that we might not need in the 
  // expression;  can edit here to only copy variables that are in the vartb (i.e. 
  // those that will actually be used in evaluating the expression)
  //
  
  std::map<std::string,std::vector<std::string> > accum_txt;
  std::map<std::string,std::vector<int> > accum_int;
  std::map<std::string,std::vector<double> > accum_dbl;
  std::map<std::string,std::vector<bool> > accum_bool;

  std::map<std::string,annot_map_t>::const_iterator ii = inputs.begin();
  
  while ( ii != inputs.end() )
    {
      const std::string & annot_name = ii->first;
      const annot_map_t & annot_map  = ii->second;
      
      // consider every instance for this annotation; make a txtvector that 
      // lists all instance IDs 
      std::vector<std::string> instance_labels;
      
      annot_map_t::const_iterator mm = annot_map.begin();
      while ( mm != annot_map.end() )
	{
	  
	  const instance_idx_t & instance_idx = mm->first;
	  const std::string & instance_id = instance_idx.id;
	  const instance_t * instance = mm->second;
	  
	  // store instance IDs as variable name == 'annot_name'
	  accum_txt[ annot_name ].push_back( instance_id );
	  
	  // store interval duration (in seconds) as variable name == 'annot_name
	  accum_dbl[ annot_name + "_sec" ].push_back( instance_idx.interval.duration_sec() );

	  // store arbitrary meta-data
	  std::map<std::string,avar_t*>::const_iterator kk = instance->data.begin();
	  while ( kk != instance->data.end() )
	    {
	      
	      // variable name is annot.var
	      
	      const std::string meta_name = annot_name + "." + kk->first;
	      
	      const avar_t * value = kk->second;
	      
	      // for now... we cannot read IN vectors
	      // is okay, as currently no way to specify those from files in any case!
	      // although note.. could be generated via a prior EVAL statement...
	      
	      globals::atype_t type = value->atype();
	      
	      if      ( type == globals::A_TXT_T ) accum_txt[ meta_name ].push_back( value->text_value() ); 
	      else if ( type == globals::A_DBL_T ) accum_dbl[ meta_name ].push_back( value->double_value() );
	      else if ( type == globals::A_INT_T ) accum_int[ meta_name ].push_back( value->int_value() );
	      else if ( type == globals::A_BOOL_T ) accum_bool[ meta_name ].push_back( value->bool_value() );
	      
	      // next meta-data
	      ++kk;
	    } 
	  
	  // next instance
	  ++mm;
	} 
      
      // next annotation
      ++ii;
    }
  

  //
  // add all to a single instance_t to hand to the token parser
  //

  instance_t m;
  
  std::map<std::string,std::vector<std::string> >::const_iterator ii1 = accum_txt.begin();
  while ( ii1 != accum_txt.end() )
    {
      if ( ii1->second.size() == 1 ) m.set( ii1->first , ii1->second[0] );
      else m.set( ii1->first , ii1->second );
      ++ii1;
    }

  std::map<std::string,std::vector<double> >::const_iterator ii2 = accum_dbl.begin();
  while ( ii2 != accum_dbl.end() )
    {
      if ( ii2->second.size() == 1 ) m.set( ii2->first , ii2->second[0] );
      else m.set( ii2->first , ii2->second );
      ++ii2;
    }

  std::map<std::string,std::vector<int> >::const_iterator ii3 = accum_int.begin();
  while ( ii3 != accum_int.end() )
    {
      if ( ii3->second.size() == 1 ) m.set( ii3->first , ii3->second[0] );
      else m.set( ii3->first , ii3->second );
      ++ii3;
    }
  
  std::map<std::string,std::vector<bool> >::const_iterator ii4 = accum_bool.begin();
  while ( ii4 != accum_bool.end() )
    {
      if ( ii4->second.size() == 1 ) m.set( ii4->first , ii4->second[0] );
      else m.set( ii4->first , ii4->second );
      ++ii4;
    }


  //
  // now add any meta-data from the global input/output class, directly to 'm'
  //
    
  if ( accumulator != NULL )
    {
      std::map<std::string,avar_t*>::const_iterator kk = accumulator->data.begin();
      while ( kk != accumulator->data.end() )
	{
	  
	  // variable name (which /should/ always start with _ if it has been written to a global store)
	  
	  const std::string meta_name =  kk->first;
	  
	  const avar_t * value = kk->second;
	  
	  // for now... we cannot read IN vectors
	  // is okay, as currently no way to specify those from files in any case!
	  
	  globals::atype_t type = value->atype();
	  
	  if      ( type == globals::A_TXT_T ) m.set( meta_name , value->text_value() ); 
	  else if ( type == globals::A_DBL_T ) m.set( meta_name , value->double_value() );
	  else if ( type == globals::A_INT_T ) m.set( meta_name , value->int_value() );
	  else if ( type == globals::A_BOOL_T ) m.set( meta_name , value->bool_value() );
	  else if ( type == globals::A_TXTVEC_T ) m.set( meta_name , value->text_vector() ); 
	  else if ( type == globals::A_DBLVEC_T ) m.set( meta_name , value->double_vector() );
	  else if ( type == globals::A_INTVEC_T ) m.set( meta_name , value->int_vector() );
	  else if ( type == globals::A_BOOLVEC_T ) m.set( meta_name , value->bool_vector() );
	  
	  // next meta-data
	  ++kk;
	} 
      
    }
  
  // Check...
  //    std::cout << m.print() << "\n";

//    if ( accumulator != NULL ) 
//      std::cout << "acc\n" << accumulator->print() << "\n";

  
  //
  // Two scenarios: 
  //  1) a value might already exist (in which case we use that type and value to specify the token )
  //  2) this could be used as part of an assignment to an instance; here we would not know the 
  //     type of the to-be-assigned field yet [ IF we are allowing assignments ] 
  //
  
  std::map<std::string,std::set<Token*> >::iterator i = vartb.begin();
  while ( i != vartb.end() )
    { 
      
      const std::string & var_name = i->first;

      std::set<Token*>::iterator tok = i->second.begin();
      
      while ( tok != i->second.end() )
	{
	  
	  globals::atype_t mt = m.type( var_name );

	  avar_t * a = m.find( var_name );
	  
	  // numeric, int, text and bool scalars
	  // also allow vectors downstream
	  
	  if ( a == NULL ) (*tok)->set(); // UNDEFINED

	  else if ( a->atype() == globals::A_INT_T  )  { (*tok)->set( a->int_value() ) ;  }
	  else if ( a->atype() == globals::A_DBL_T  )  { (*tok)->set( a->double_value() ); }
	  else if ( a->atype() == globals::A_TXT_T  )  { (*tok)->set( a->text_value() ); }
	  else if ( a->atype() == globals::A_BOOL_T )  { (*tok)->set( a->bool_value() );   }	      

	  else if ( a->atype() == globals::A_INTVEC_T  )  { (*tok)->set( a->int_vector() ) ;  }
	  else if ( a->atype() == globals::A_DBLVEC_T  )  { (*tok)->set( a->double_vector() ); }
	  else if ( a->atype() == globals::A_TXTVEC_T  )  { (*tok)->set( a->text_vector() ); }
	  else if ( a->atype() == globals::A_BOOLVEC_T )  { (*tok)->set( a->bool_vector() );   }	      
	  
	  else (*tok)->set(); // UNDEFINED
	    	  
	  ++tok;
	}
      ++i;
    }    

  //
  // For assignment
  //
  
  func.attach( outputs , accumulator , global_vars );
  
}





void Eval::bind( const std::map<std::string,std::vector<double> > & inputs , 
		 instance_t * outputs )
{
  
  // assuming input from signal data, i.e. numeric vector(s) 
  reset_symbols();  
  
  //
  // create a single instance that organizes all the information in
  // the annot_map_t above in a sane variable naming scheme
  // 
  // nb. some redundancy here, as we copy stuff over that we might not need in the 
  // expression;  can edit here to only copy variables that are in the vartb (i.e. 
  // those that will actually be used in evaluating the expression)
  //
  
  std::map<std::string,std::vector<double> > accum_dbl;
  

  //
  // add 'inputs' to a single instance_t to hand to the token parser
  //

  instance_t m;
  
  std::map<std::string,std::vector<double> >::const_iterator ii = inputs.begin();
  while ( ii != inputs.end() )
    {
      m.set( ii->first , ii->second );
      ++ii;
    }
  
  
  //
  // Two scenarios: 
  //  1) a value might already exist (in which case we use that type and value to specify the token )
  //  2) this could be used as part of an assignment to an instance; here we would not know the 
  //     type of the to-be-assigned field yet [ IF we are allowing assignments ] 
  //
  
  std::map<std::string,std::set<Token*> >::iterator i = vartb.begin();
  while ( i != vartb.end() )
    { 
      
      const std::string & var_name = i->first;

      std::set<Token*>::iterator tok = i->second.begin();
      
      while ( tok != i->second.end() )
	{
	  
	  globals::atype_t mt = m.type( var_name );
	  
	  avar_t * a = m.find( var_name );
	  
	  // numeric, int, text and bool scalars
	  // also allow vectors downstream
	  
	  if ( a == NULL ) (*tok)->set(); // UNDEFINED
	  
	  else if ( a->atype() == globals::A_INT_T  )  { (*tok)->set( a->int_value() ) ;  }
	  else if ( a->atype() == globals::A_DBL_T  )  { (*tok)->set( a->double_value() ); }
	  else if ( a->atype() == globals::A_TXT_T  )  { (*tok)->set( a->text_value() ); }
	  else if ( a->atype() == globals::A_BOOL_T )  { (*tok)->set( a->bool_value() );   }	      

	  else if ( a->atype() == globals::A_INTVEC_T  )  { (*tok)->set( a->int_vector() ) ;  }
	  else if ( a->atype() == globals::A_DBLVEC_T  )  { (*tok)->set( a->double_vector() ); }
	  else if ( a->atype() == globals::A_TXTVEC_T  )  { (*tok)->set( a->text_vector() ); }
	  else if ( a->atype() == globals::A_BOOLVEC_T )  { (*tok)->set( a->bool_vector() );   }	      
	  
	  else (*tok)->set(); // UNDEFINED
	    	  
	  ++tok;
	}
      ++i;
    }    

  //
  // For assignment
  //
  
  func.attach( outputs );
  
}




void Eval::bind( const Token * ntok )
{

  // we have have just created a meta-field from an assignment --
  // check whether this value now needs to be bound to any other
  // variable tokens (i.e. in next expression). ie. this function is
  // only called after an assignment is made and a meta-variable newly
  // created
  
  std::map<std::string,std::set<Token*> >::iterator i = vartb.find( ntok->name() );
  if ( i != vartb.end() )
    { 
      std::set<Token*>::iterator tok = i->second.begin();
      
      while ( tok != i->second.end() )
	{	  

	  if ( ntok == *tok ) { ++tok; continue; } 
	  
	  // we assume this new token has a well-defined type
	  // copy it over
	  
	  *(*tok) = *ntok;
	  
	  ++tok;	  
	}
    }
}

bool Eval::evaluate( const bool v )
{
  verbose = v;
  for (int i=0; i<neval; i++)
    {
      if ( verbose )
	std::cerr << " Prior to expression " << i+1 << " status = " << ( is_valid ? "VALID" : "INVALID" ) << "\n";

      if ( is_valid ) 
	is_valid = execute( output[i] );

      if ( verbose )
	std::cerr << " Post to expression " << i+1 << " status = " << ( is_valid ? "VALID" : "INVALID" ) << "\n";

    }
  
  if ( verbose ) std::cerr << " returning " << ( is_valid ? "VALID" : "INVALID" ) << " token\n";
  return is_valid;
}

bool Eval::valid() const 
{
  return is_valid;
}

std::string Eval::errmsg() const
{
  return errs;
}

void Eval::errmsg( const std::string & e )
{
  errs += e + "\n";
}


Token::tok_type Eval::rtype() const
{
  return e.type();
}

bool Eval::value(bool & b)
{

  if ( e.is_bool(&b) ) return true;

  int i;
  if ( e.is_int(&i) ) { b=i; return true; }

  // for vectors, evaluate to T is at least one T

  std::vector<bool> bv;
  if ( e.is_bool_vector(&bv) ) 
    {
      b = false;
      for (int i=0;i<bv.size();i++) 
	if ( bv[i] ) { b=true; break; }
      return true;
    }
  
  std::vector<int> iv;
  if ( e.is_int_vector(&iv) ) 
    {
      b = false;
      for (int i=0;i<iv.size();i++) 
	if ( iv[i] ) { b=true; break; }
      return true;
    }

  return false;
}

Token Eval::value() const
{
  return e;
}

bool Eval::value(int & i)
{
  if ( e.is_int(&i) ) return true;
  bool b;
  if ( e.is_bool(&b) ) { i=b; return true; }
  return false;
}

bool Eval::value(double & d)
{
  if ( e.is_float(&d) ) return true;
  int i;
  if ( e.is_int(&i) ) { d=i; return true; }
  bool b;
  if ( e.is_bool(&b) ) { d=b; return true; }
  return false;
}

bool Eval::value(std::string & s)
{
  if ( e.is_string(&s) ) return true;
  return false;
}


std::string Eval::result() const
{
  std::stringstream ss;
  ss << e;
  return ss.str();
}


void Eval::reset_symbols()
{
  // clear all variables -- put Tokens in the UNDEFINED state
  std::map<std::string,std::set<Token*> >::iterator i = vartb.begin();
  while ( i != vartb.end() )
    {
      std::set<Token*>::iterator k = i->second.begin();
      while ( k != i->second.end() )
	{
	  (*k)->set();
	  ++k;
	}
      ++i;
    }
  e.set(); // also clear end expression
}


void Eval::delete_symbols() 
{
  vartb.clear();
}


bool Eval::expand_vargs( std::string * s )
{
  // change vec(22,23,24) to vec(22,23,24,3)
  //        vec()         to vec(0)
  //        vec(22)       to vec(22,1)
  // etc i.e. put VA at end (as will come off stack first)


  std::vector<std::string> fname;
  fname.push_back("num(");
  fname.push_back("int(");
  fname.push_back("txt(");
  fname.push_back("bool(");
  fname.push_back("c(");

  for (int f = 0 ; f < fname.size() ; f++)
    while ( 1 ) 
      {
	
	size_t p = s->find( fname[f] );
	
	if ( p != std::string::npos ) 
	  {	      
	    
	    // this will also match 
	    //     blahint(X)
	    
	    if ( p > 1 )
	      {
		if ( ( (*s)[p-1] >= 'A' && (*s)[p-1] <= 'Z' ) || 
		     ( (*s)[p-1] >= 'a' && (*s)[p-1] <= 'z' ) || 
		     ( (*s)[p-1] >= '0' && (*s)[p-1] <= '9' ) || 
		     (*s)[p-1] >= '_' ) { p = std::string::npos; } 
	      }

	  }
	
	bool found = false;

	if ( p != std::string::npos )
	  {
	    
	    found = true;
	    
	    // count # of upper level ARGs, by counting commas until the next matching ")"
	    
	    //  vec( 1 , 2 , 3 )     -->  2 commas , 3 args
	    //  int( pow(2,3) , 7 )  ---> 1 comma  , 2 args
	    	    
	    // look for closing brace
	    
	    int bc = 0;
	    int q = p;
	    int pp = p;
	    int commacnt = 0;

	    while ( ++q )
	      {
		
		// gone past end of string?
		if ( q == s->size() ) return false;
		
		char c = s->substr(q,1)[0];
		
		if ( c == '(' ) 
		  {
		    ++bc;  //includes first paran		    
		  }
		else if ( c == ')' )
		  {
		    --bc;
		    
		    // end of expression
		    if ( bc == 0 ) 
		      {
			break;
		      }
		  }
		else if ( bc == 1 && c == ',' ) ++commacnt;
	      }
	    
	    int len = q-p+1;	    
	    std::string orig = s->substr(p,q-p+1);    
	    orig = fname[f].substr( 0 , fname[f].size() - 1 ) + "_func("  + orig.substr( fname[f].size() ) ;
	    orig = orig.substr(0,orig.size()-1);
	    orig += "," + Helper::int2str( commacnt+1 ) + ")";
	    
	    s->replace( p , len , orig );
	    	    
	  }
	
	// are we sure there are no more gfunc()s in input?
	if ( ! found ) break;
	
      }	  

  return true;
}


bool Eval::expand_indices( std::string * s )
{

  // simply convert  X[i]  to  xfn(X,i)
  
  while ( 1 ) 
    {
      
      // search for opening index
      // assume it comes after a variable
      //  x[y]
      //  backtrack to get 'x'
      //  then forwards to get 'y'
      //  then make extract(x,y)
      // valid delimiters backwards are ' ' \n \t , % < > & | ! ~ = * + - / ; { } (

      // but if we hit ")" before anything else (except whitespace) then we need to jump to the next 
      // "(" and then look for an identifier 
      
      
      int p = s->find( "[" );
      if ( p == std::string::npos ) return true;
      int q = p;
      bool anything = false;

      while ( --q )
	{
	  if ( q < 0 ) return false;
	  if ( q == 0 ) { ++q; break; }
	  char c = s->substr(q,1)[0];
	  
	  if ( c == ')' )
	    {
	      // jump to next '(', allowing nesting
	      int nest = 1;
	      while (1) 
		{
		  --q;
		  if ( s->substr(q,1) == ")" ) ++nest;
		  else if ( s->substr(q,1) == "(" ) --nest;
		  if ( nest == 0 ) break; 
		}
	    }

	  if ( c == ',' || c == '&' || c == '%' || c == '>' || c == '<' || c == '|' || c == '(' ||
	       c == '!' || c == '~' || c == '^' || c == '=' || c == '*' || c == '+' || c == '-' || 
	       c == '/' || c == ';' || c == ':'  )
	    { ++q; break; }	  
	  
	  // whitespace is a delimiter (outside of parens) 
	  
	  if ( c == ' ' || c == '\n' || c == '\t' )
	    {
	      if ( anything ) { ++q; break; } 
	    }
	  else 
	    anything = true; // i.e. we must have encountered an identifier

	}
      
      std::string vec_idx = s->substr(q,p-q) ;
  
      std::string  arg_idx;
      int r = p;      
      while ( ++p )
	{  
	  if ( p == s->size() ) return false;
	  char c = s->substr(p,1)[0];
	  // do not allow nested indexing
	  if ( c == '[' ) return false;	  
	  if ( c == ']' )
	    {
	      arg_idx = s->substr(r+1,p-r-1);
	      break;
	    }  
	}
      
      std::string label = "element(" + vec_idx + "," + arg_idx + ")";
      s->replace( q , (p-q+1) , label );      

    } // search for next []  

  return true;    
}



void Eval::locate_symbols( std::vector<Token> & tok ) 
{
  
  // set pointers to output table for the updating, as this 
  // does not shift between runs
  
  for (int i=0; i<tok.size(); i++)
    if ( tok[i].is_variable() )
      {
	vartb[ tok[i].name() ].insert( &(tok[i]) );
      }
  reset_symbols(); // symbol table
}

