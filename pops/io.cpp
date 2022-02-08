
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

#ifdef HAS_LGBM

#include "pops/pops.h"
#include "helper/helper.h"
#include "helper/logger.h"

extern logger_t logger;

void pops_indiv_t::save1( const std::string & id , const std::string & f )
{
  // ID, ne, nf { features , row-major }   
  std::ofstream OUT1( Helper::expand( f ).c_str() , std::ios::binary | std::ios::out );
  bwrite( OUT1, id ) ;
  bwrite( OUT1, ne ) ;
  bwrite( OUT1, pops_t::specs.n1 );  
  for (int i=0; i<ne; i++)
    {
      // epoch number
      bwrite( OUT1, E[i] );
      // sleep stage
      bwrite( OUT1, S[i] );
      // l1-features
      for (int j=0; j<pops_t::specs.n1; j++)
	bwrite( OUT1, X1(i,j) );
    }
  OUT1.close();
}


void pops_t::load1( const std::string & f )
{
  std::ifstream IN1( Helper::expand( f ).c_str() , std::ios::binary | std::ios::in );

  int total_epochs = 0;
  int n_indiv = 0;
  
  // get size of data first 
  while ( 1 )
    {
      std::string id = pops_indiv_t::bread_str( IN1 );
      if ( IN1.eof() || IN1.bad() ) break;
      ++n_indiv;
      
      int ne1 = pops_indiv_t::bread_int( IN1 );
      total_epochs += ne1;
      
      int nf1 = pops_indiv_t::bread_int( IN1 );
      if ( nf1 != pops_t::specs.n1 )
	Helper::halt( "data in " + f + " does not match feature-specification file" );

      // skip rest (not, order diff. but obvs doesn't matter):
      // epochs, stages
      pops_indiv_t::bskip_int( IN1 , ne1 * 2 );
      // skip features
      pops_indiv_t::bskip_dbl( IN1 , ne1 * nf1 );

    }
  IN1.close();

  logger << "  reading " << total_epochs << " epochs from " << n_indiv << " individuals\n";
  
  X1.resize( total_epochs , pops_t::specs.n1 );
  S.resize( total_epochs );
  E.resize( total_epochs );
  
  // track when indivs start/stop
  Istart.clear();
  Iend.clear();
  
  // re-read
  std::ifstream IN2( Helper::expand( f ).c_str() , std::ios::binary | std::ios::in );
  total_epochs = 0;
  
  while ( 1 )
    {
      std::string id = pops_indiv_t::bread_str( IN2 );
      if ( IN2.eof() || IN2.bad() ) break;
      
      int ne1 = pops_indiv_t::bread_int( IN2 );      
      int nf1 = pops_indiv_t::bread_int( IN2 );
      
      Istart.push_back( total_epochs );
	    
      for (int i=0; i<ne1; i++)
	{

	  // epoch
	  E[ total_epochs ] = pops_indiv_t::bread_int( IN2 );

	  // stage
	  S[ total_epochs ] = pops_indiv_t::bread_int( IN2 );

	  // features
	  for (int j=0; j<pops_t::specs.n1; j++)
	    X1( total_epochs , j ) = pops_indiv_t::bread_dbl( IN2 );	  
	  ++total_epochs;
	}
      
      Iend.push_back( total_epochs - 1 );
      
    }
  IN2.close();
  
}




void pops_indiv_t::bwrite( std::ofstream & O , const std::string & s ) 
{
  uint8_t l = s.size();
  O.write( (char*)( &l ), sizeof(uint8_t) );
  O.write( s.c_str(), l );
}

void pops_indiv_t::bwrite( std::ofstream & O , int i ) 
{
  O.write( (char*)( &i ), sizeof(int) );
}

void pops_indiv_t::bwrite( std::ofstream & O , double d ) 
{
  O.write( (char*)( &d ), sizeof(double) );
}

std::string pops_indiv_t::bread_str( std::ifstream & I )
{
  uint8_t len;
  I.read( (char*)( &len ), sizeof(uint8_t) );
  std::vector<char> b( len );
  I.read( &b[0] , len );
  std::string s( b.begin() , b.end() );
  return s;
}

int pops_indiv_t::bread_int( std::ifstream & I )
{
  int i;
  I.read( (char*)( &i ), sizeof(int) );
  return i;
}

double pops_indiv_t::bread_dbl( std::ifstream & I )
{
  double d;
  I.read( (char*)( &d ), sizeof(double) );
  return d;
}

void pops_indiv_t::bskip_dbl( std::ifstream & I , const int n )
{
  std::vector<double> dummy( n ) ;
  I.read( (char*)( &dummy[0] ), n * sizeof(double) );
}

void pops_indiv_t::bskip_int( std::ifstream & I , const int n )
{
  std::vector<double> dummy( n ) ;
  I.read( (char*)( &dummy[0] ), n * sizeof(int) );
}



#endif
