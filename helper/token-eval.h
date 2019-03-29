
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


#ifndef __LUNA_TOKEN_EVAL_H__
#define __LUNA_TOKEN_EVAL_H__

#include "token.h"

#include <cstring>
#include <ostream>
#include <map>
#include <set>
#include <vector>

struct instance_t;
struct instance_idx_t;

typedef std::map<instance_idx_t,instance_t*> annot_map_t;

struct Eval {
    
  friend std::ostream & operator<<( std::ostream & out , Eval & rhs )
  {
    out << rhs.e;
    return out;
  }

  
 public:    
  
  // constructors

  Eval( const bool na = false ) 
    { init( na ); }

  Eval( const std::string & input , const bool na = false ) 
    { init( na ); parse(input); }
  
  void init( const bool ); 

  // primary functions
  
  bool parse( const std::string & input );

  void assign_to( instance_t & m );

  void bind( const std::map<std::string,annot_map_t> & inputs , 
	     instance_t * outputs , 
	     instance_t * accumulator = NULL , 
	     const std::set<std::string> * global_vars = NULL , 
	     bool reset = true );
  
  void bind( const Token * );

  bool evaluate( const bool v = false );
    
  bool valid() const;
  
  std::string errmsg() const;

  void errmsg(const std::string & );

  std::string result() const;

  // queries into value of expression
  bool value(bool & b);
  bool value(int & );
  bool value(double &);
  bool value(std::string &);
  
  // what does this return?
  Token::tok_type rtype() const;
  Token value() const;
  
 private:
  
  //work horses
  bool get_token( std::string & input ,  Token & );
  bool previous_value;
  bool execute( const std::vector<Token> & );
  bool shunting_yard( const std::string & input, std::vector<Token> & );  
  
  // helpers
  int op_preced(const Token & tok );
  bool op_left_assoc(const Token & tok );
  unsigned int op_arg_count(const Token & tok );

  // expression in RPN notation (post parse) (for each eval)
  std::vector< std::vector<Token> >output;   
  
  // keep track of state (errors?)
  bool is_valid;
  std::string errs;
  
  // slot for final expression value
  Token e;
  
  // Symbol table for variables, and related functions
  bool expand_indices( std::string * s );
  bool expand_vargs( std::string * s );
  std::map<std::string,std::set<Token*> > vartb;
  void delete_symbols();
  void reset_symbols();
  void locate_symbols( std::vector<Token> & );

  // Keep token-functions in one place
  TokenFunctions func;

  // the number of evals
  int neval;
  
  // run in no-assignment mode (i.e. MASK)
  bool no_assignments;

  // verbose mode?
  bool verbose;

};


#endif

