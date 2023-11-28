
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

#include "models/model.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "helper/helper.h"

#include "helper/logger.h"

extern logger_t logger;

// special terms required/expected in a model
// title     <- "Title of the model"
// reference <- "PMID, URL, or other citation"
// outcome   <- "label for the predicted measure"
// type      <- "linear or logistic"
// training  <- "brief description of training population (N=XXXX)" 
// data      <- "filepath for training data"

void prediction_model_t::read( const std::string & f , const std::string & id , const bool cacheless )
{
  
  // if cacheless, allow variables in the file to not be specified (i.e. puts in a blank)
  
  const std::string filename = Helper::expand( f );
  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not open " + f ) ;

  terms.clear();
  specials.clear();
  
  // process any variable substitutions?
  std::map<std::string,std::string>  allvars = cmd_t::indiv_var_map( id );
  
  // build current term

  model_term_t term;

  bool in_term = false;

  // load

  std::ifstream IN1( filename.c_str() , std::ios::in );
    
  while ( ! IN1.eof() ) 
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() ) break;
      
      // blank line - skip, but also ends any term being
      // defined

      if ( line == "" )
	{
	  if ( in_term )
	    {
	      terms.insert( term );
	      term.clear();
	      in_term = false;
	    }
	  
	  continue; // skip blanks
	}

      // comments
      if ( line[0] == '%' ) continue; 

      
      // variable substitution : if cacheless == T, will allow variables to
      // be missing (e.g. ${ch} , as we don't need this when reading from 
      // the vars table (instead of the cache) 
      Helper::swap_in_variables( &line , &allvars , cacheless );

      // parse on whitespace
      std::vector<std::string> tok = Helper::quoted_parse( line , "\t " );

      // same as blank (end term) 
      if ( tok.size() == 0 )
	{
	  if ( in_term )
            {
              terms.insert( term );
              term.clear();
              in_term = false;
            }
	  
	  continue;
	}
		    
      
      // special assignments are always on a single line
      //    key <- value
      if ( tok.size() == 3 && tok[1] == "<-" )
	{
	  if ( in_term )
	    Helper::halt( "bad syntax: cannot have a special assignment mid-term: " + line );

	  // string?
	  if ( tok[2][0] == '"' )
	    {
	      specials_str[ tok[0] ] = Helper::unquote( tok[2] );
	    }
	  else // numeric
	    {
	      // allows missing value via '.' 
	      if ( tok[2] != "." )
		{
		  double x;
		  if ( ! Helper::str2dbl( tok[2] , &x ) )
		    Helper::halt( "could not convert to a numeric value (use period for missing value) : " + line );
		  specials[ tok[0] ] = x;
		}
	    }
	  
	  // next line
	  continue;
	}

      //
      // if not in_term, then next term must be a term label (no =)
      //

      int start = 0;
      
      if ( ! in_term )
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[0] , '=' );
	  if ( tok2.size() != 1 ) Helper::halt( "expecting new label, no = assignments" );
	  term.label = tok[0];
	  in_term = true;
	  start = 1;
	}
      
      //
      // parse remainders
      //

      if ( in_term )
	{
	  for (int i=start; i<tok.size(); i++)
	    {

	      // is this a new term label (which ends the old one?) 
	      std::vector<std::string> tok2 = Helper::parse( tok[i] , '=' );
	      if ( tok2.size() == 1 )
		{
		  // save this term
		  terms.insert( term );
		  // start new term
		  term.clear();
		  term.label = tok2[0];
		}
	      else
		{
		  // we have some key=value pairing
		  if ( tok2.size() != 2 ) Helper::halt( "bad key=value syntax: " + line );

		  const std::string key = Helper::toupper( tok2[0] );
		  const std::string value = tok2[1];

		  //std::cout << "key/value=[" << key << "][" << value << "]\n";
		  
		  if      ( key == "CMD" )
		    term.cmd = value;
		  else if ( key == "VAR" )
		    term.var = value;
		  else if ( key == "VALUE" )
		    {
		      term.value = value; // may be missing at this point
		      term.has_value = true;
		    }
		  else if ( key == "CH" )
		    term.chs = Helper::parse( value , ',' );		  
		  else if ( key == "CHS" ) // expecting C1+C2,A1+A2,<etc>
		    term.pairs = Helper::parse( value , ',' ); 
		  else if ( key == "STRATA" && value != "." ) // allow baseline strata = empty
		    term.strata = Helper::mapize( value  , ',' , '/' );
		  else if ( key == "B" )
		    {
		      double x;
		      if ( ! Helper::str2dbl( value , &x ) ) Helper::halt( "bad numeric value: " + line );
		      term.coef = x;
		    }
		  else if ( key == "M" )
		    {
		      double x;
		      if ( ! Helper::str2dbl( value , &x ) ) Helper::halt( "bad numeric value: " + line );
		      term.mean = x;
		    }
		  else if ( key == "SD" )
		    {
		      double x;
                      if ( ! Helper::str2dbl( value , &x ) ) Helper::halt( "bad numeric value: " + line );
		      term.sd = x;
		    }
		  else if ( key == "REQ" )
		    {
		      term.required = Helper::yesno( value );
		    }
		  else if ( key == "LOG" )
		    {
		      term.log_transform = Helper::yesno( value );
		    }		 
		  else if ( key == "DIR" )
		    {
		      term.directed = Helper::yesno( value );
		    } 
		  else
		    Helper::halt( "unrecognized key term: " + key );
		}
	      
	    } // next token

	} // end of processing in-term line 
      
    } // next line

  // final term?
  if ( in_term )
    terms.insert( term );
  
  IN1.close();

  logger << "  read " << terms.size()
	 << " terms and " << specials.size()
	 << " special variables from " << filename << "\n";
  
  // complain if additional variables have not been specified
  
  if ( specials_str.find( "title" ) == specials_str.end() )
    logger << "  *** no 'title' specified ***\n";
  if ( specials_str.find( "outcome" ) == specials_str.end() )
    logger << "  *** no 'outcome' specified ***\n";
  if ( specials_str.find( "reference" ) == specials_str.end() )
    logger << "  *** no 'reference' specified ***\n";
  if ( specials_str.find( "training" ) == specials_str.end() )
    logger << "  *** no 'training' information specified ***\n";
  if ( specials_str.find( "type" ) == specials_str.end() )
    logger << "  *** no 'type' information (linear/logistic) specified ***\n";  
  
}
      
void prediction_model_t::populate()
{
  int nt = size();
  coef = Eigen::VectorXd::Zero( nt );
  mean = Eigen::VectorXd::Zero( nt );
  sd = Eigen::VectorXd::Zero( nt );
  int i = 0;  
  std::set<model_term_t>::const_iterator tt = terms.begin();
  while ( tt != terms.end() )
    {
      coef[i] = tt->coef;
      mean[i] = tt->mean;
      sd[i] = tt->sd;
      ++i;
      ++tt;
    }
}

// std::set<std::string> prediction_model_t::channels() const
// {
//   std::set<std::string> chs;
//   std::set<model_term_t>::const_iterator tt = terms.begin();
//   while ( tt != terms.end() )
//     {
//       // any 'channel' labels?
//       std::map<std::string,std::string>::const_iterator cc = tt->strata.find( "CH" );
//       ++tt;
//     }

//   return chs;
// }


void prediction_model_t::dump() const
{

  std::cout << "% dumping current parsed model\n\n";
  
  if ( specials_str.size() > 0 )
    {      
      std::map<std::string,std::string>::const_iterator qq = specials_str.begin();
      while ( qq != specials_str.end() )
	{
	  std::cout << "  " << qq->first << " <- \"" << qq->second << "\"\n";
	  ++qq;
	}
      std::cout << "\n";
    }
  

  if ( specials.size() > 0 )
    {
      std::map<std::string,double>::const_iterator ss = specials.begin();
      while ( ss != specials.end() )
	{
	  std::cout << "  " << ss->first << " <- " << ss->second << "\n";
	  ++ss;
	}
      std::cout << "\n";
    }
  
  std::set<model_term_t>::const_iterator tt = terms.begin();
  while ( tt != terms.end() )
    {
      if ( tt->has_value )
	std::cout << tt->label << "\n"
		  << "  value=" << tt->value << " "
		  << "req=" << tt->required << " "
		  << "log=" << tt->log_transform << "\n"
		  << "  b=" << tt->coef << " "
		  << "m=" << tt->mean << " "
		  << "sd=" << tt->sd << "\n\n";
      else
	{
	  std::cout << tt->label << "\n"
		    << "  cmd=" << tt->cmd << " "
		    << "var=" << tt->var << " "
		    << "req=" << tt->required << " "
		    << "log=" << tt->log_transform << " ";
	  if ( tt->chs.size() ) 
	    std::cout << "ch=" << Helper::stringize( tt->chs ) << " ";
	  if ( tt->pairs.size() )
	    std::cout << "chs=" << Helper::stringize( tt->pairs ) << " ";

	  std::cout << "strata=" << Helper::ezipam( tt->strata , ',', '/' ) << "\n"
		    << "  b=" << tt->coef << " "
		    << "m=" << tt->mean << " "
		    << "sd=" << tt->sd << "\n\n";
	}
      ++tt;
    }
  
  
}

std::vector<std::string> prediction_model_t::header() const
{
  std::vector<std::string> h;
  std::set<model_term_t>::const_iterator tt = terms.begin();
  while ( tt != terms.end() )
    {
      h.push_back( tt->label );
      ++tt;
    }
  return h;
}
