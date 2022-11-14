
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

  logger << "  writing binary data (" << X1.rows() << " " << ne << " epochs, " 
	 << X1.cols() << " features) to " << f << "\n";
  
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
  
  ni_validation = 0;
  int ne_validation = 0;
  int ne_training = 0;

  // get size of data first 
  while ( 1 )
    {
      std::string id = pops_indiv_t::bread_str( IN1 );
      if ( IN1.eof() || IN1.bad() ) break;
      ++n_indiv;
      
      int ne1 = pops_indiv_t::bread_int( IN1 );
      total_epochs += ne1;
      
      if ( holdouts.find( id ) != holdouts.end() )
	{
	  ne_validation += ne1;
	  ++ni_validation;
	}
      else
	ne_training += ne1;
      
      
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

  logger << "  reading " << total_epochs << " epochs from " << n_indiv << " individuals"
	 << " (" << ni_validation << " of whom held back for model validation)\n";

  // store these for splitting X1 later when passing to LGBM
  nrows_training = ne_training;
  nrows_validation = ne_validation;
  
  // store : note, here set to total (training + validation size) for X1 (only)
  X1.resize( total_epochs , pops_t::specs.n1 );
  E.resize( ne_training ); // ...
  S.resize( ne_training ); // will push_back() below
  
  // track when indivs start/stop, and who they are
  Istart.clear();
  Iend.clear();
  I.clear();

  // validation people (will be added to end of X1, S and E, etc
  Eigen::MatrixXd X2 = Eigen::MatrixXd::Zero( ne_validation , pops_t::specs.n1 );
  std::vector<int> S2( ne_validation );
  std::vector<int> E2( ne_validation );
  std::vector<int> Istart2, Iend2;
  std::vector<std::string> I2;

  // re-read
  std::ifstream IN2( Helper::expand( f ).c_str() , std::ios::binary | std::ios::in );

  int offset = ne_training; // i.e. for ultimate Istart[i] -> for appended validations
  ne_training = 0;
  ne_validation = 0;

  while ( 1 )
    {
      std::string id = pops_indiv_t::bread_str( IN2 );
      if ( IN2.eof() || IN2.bad() ) break;
      
      const bool is_training = holdouts.find( id ) == holdouts.end();
      
      int ne1 = pops_indiv_t::bread_int( IN2 );      
      int nf1 = pops_indiv_t::bread_int( IN2 );
      
      if ( is_training )
	{
	  I.push_back( id );
	  Istart.push_back( ne_training );
	  for (int i=0; i<ne1; i++)
	    {
	      // epoch
	      E[ ne_training ] = pops_indiv_t::bread_int( IN2 );
	      // stage
	      S[ ne_training ] = pops_indiv_t::bread_int( IN2 );
	      // features
	      for (int j=0; j<pops_t::specs.n1; j++)
		X1( ne_training , j ) = pops_indiv_t::bread_dbl( IN2 );
	      ++ne_training;
	    }
	  Iend.push_back( ne_training - 1 );
	}
      else  // validation individul
	{
	  I2.push_back( id );
	  Istart2.push_back( offset + ne_validation );
	  for (int i=0; i<ne1; i++)
	    {
	      // epoch
	      E2[ ne_validation ] = pops_indiv_t::bread_int( IN2 );
	      // stage
	      S2[ ne_validation ] = pops_indiv_t::bread_int( IN2 );
	      // features
	      for (int j=0; j<pops_t::specs.n1; j++)
		X2( ne_validation , j ) = pops_indiv_t::bread_dbl( IN2 );	  
	      ++ne_validation;
	    }
	  Iend2.push_back( offset + ne_validation - 1 );
	}
    }
  IN2.close();
  
  // now concatenate training and validation samples  
  
  // space was already pre-allocated for X1
  X1.bottomRows( ne_validation ) = X2; 
  
  // extend indiv-level stores
  for (int i=0; i<Istart2.size(); i++)
    {
      Istart.push_back( Istart2[i] );
      Iend.push_back( Iend2[i] );      
      I.push_back( I2[i] );
    }
  
  for (int i=0; i<ne_validation; i++)
    {
      S.push_back( S2[i] );
      E.push_back( E2[i] );
    }

}



void pops_t::load1_stages_only( const std::string & f )
{

  // only populate S, E and I -- based on ALL individuals in the file
  // alignment does not need to worry about trainers versus validation, etc
  
  std::ifstream IN1( Helper::expand( f ).c_str() , std::ios::binary | std::ios::in );
  
  int total_epochs = 0;
  int n_indiv = 0;
  int ne_training = 0;

  // track when indivs start/stop, and who they are
  E.clear();
  S.clear();
  Istart.clear();
  Iend.clear();
  I.clear();
  
  // step through file
  while ( 1 )
    {
      std::string id = pops_indiv_t::bread_str( IN1 );
      if ( IN1.eof() || IN1.bad() ) break;
      ++n_indiv;

      // number of epochs (to read)
      int ne1 = pops_indiv_t::bread_int( IN1 );
      total_epochs += ne1;

      // number of features (to skip) 
      int nf1 = pops_indiv_t::bread_int( IN1 );
      
      // read epochs and stages
      // pops_indiv_t::bskip_int( IN1 , ne1 * 2 );
      
      I.push_back( id );
      Istart.push_back( ne_training );

      for (int i=0; i<ne1; i++)
	{
	  // epoch
	  E.push_back( pops_indiv_t::bread_int( IN1 ) );

	  // stage
	  S.push_back(pops_indiv_t::bread_int( IN1 ) );

	  // skip features
	  pops_indiv_t::bskip_dbl( IN1 , nf1 );

	  // track epoch numbering 
	  ++ne_training;
	}
      Iend.push_back( ne_training - 1 );

      // next individual
    }

  IN1.close();
  
  logger << "  read " << total_epochs << " stages from " << n_indiv << " individuals\n";
  
  // now S, E, I* populated, but not X1
  // i.e this function only called by standalone es-priors call
  
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
