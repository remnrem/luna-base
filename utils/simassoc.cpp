
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


// internal tool: to simulate and assess association between spectra and quantitative
// phenotypes using;  not supported or designed for external use otherwise

#include "simassoc.h"

#include "../helper/helper.h"
#include "../helper/logger.h"
#include "../db/db.h"
#include "../stats/eigen_ops.h"
#include "../eval.h"
#include "../main.h" 
#include "../miscmath/crandom.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

extern logger_t logger;
extern writer_t writer;

int main( int argc , char ** argv )
{

  CRandom::srand( time(0) );
  
  param_t param;
  build_param_from_cmdline( &param );

  simassoc_t sim;  

  //
  // Get data
  //

  sim.load( param.requires( "data" ) );
  
  if ( param.has( "covar" ) )
    sim.load_covar( param.value( "covar" ) );
  else
    sim.Z = Eigen::MatrixXd::Zero( 0 , 0 );

  //
  // Get generative model
  //
  
  sim.generative_model();


  //
  // Begin the primary simulation loop
  //
  
  int nrep = param.requires_int( "nrep" );
  
  //  sim.clear();
  
  for (int r=0; r<nrep; r++)
    {
      sim.simulate( );
      sim.assoc();
      // sim.collate();
    }
  
  sim.output();
  
  return 0;
}

void simassoc_t::load( const std::string & f )
{
  const std::string fdata = Helper::expand( f );
  if ( ! Helper::fileExists( fdata ) )
    Helper::halt( "could not open " + fdata );
  X = eigen_ops::load_mat( fdata , &hdr );
  logger << "read " << X.rows() << " indivs, " << X.cols() << " features from " << fdata << "\n";  
}

void simassoc_t::load_covar( const std::string & f )
{
  const std::string fdata = Helper::expand( f );
  if ( ! Helper::fileExists( fdata ) )
    Helper::halt( "could not open " + fdata );
  Z = eigen_ops::load_mat( fdata , &hdr );
  logger << "read " << Z.rows() << " indivs, " << Z.cols() << " covariates from " << fdata << "\n";  
}


void simassoc_t::describe_cols()
{

  // expect always in the form:: VAR_F/B/N_SS

  //  PSD_SIGMA_N2
  //  PSD_12.3_N2

  // Vars:
  //  PSD
  //  PSC
  //  PER (periodic PSD)
  //  APER (aperiod PSD) -->> or just SLOPE/INTERCEPT?
  
  const int nv = hdr.size();

  // for (int v=0; v<nv; v++)
  //   {
  //     std::vector<std::string> tok = Helper::parse( hdr[v] , "_" );
  //     if ( tok.size() != 3 ) Helper::halt( "problem with format : " + hdr[v] );

  //     const bool is_psd = tok[0] == "PSD";
  //     // const bool is_psd = tok[0] == "PSC";
  //     // const bool is_psd = tok[0] == "PSD";

  //     if ( is_psd )
  // 	{	  
  // 	  double f;	  
  // 	  const bool assume_band = Helper::str2dbl( tok[1] , &f );
  // 	  if ( ! assume_band ) col2f[ hdr[v] ] = f ;
  // 	  else col2f[ hdr[v] ] = 0 ; // i.e. for bands, set all frequencies to 0, as we will not allow clustering	  

	 
  // 	}

  //     col2var[ hdr[v] ] = tok[0];
  //     col2ch1[ hdr[v] ] = "S1";
  //     col2ch2[ hdr[v] ] = "S2";
  //     col2t[ hdr[v] ] = 0;
}


void simassoc_t::generative_model()
{
  const int ni = X.rows();
  const int nv = X.cols();
  W = Eigen::VectorXd::Zero( ni, nv );
}

void simassoc_t::simulate()
{
  const int ni = X.rows();
  Y = Eigen::VectorXd::Zero( ni );
  for (int i=0; i<ni; i++) Y[i] = CRandom::rand();
  std::cout << "Y\n"<< Y << "\n";
}
  
void simassoc_t::assoc()
{

  // 
  
  // cpt_t cpt( Y , X , Z ) ;

  // // will only have variation in 'f'
  // //  i.e. do all this single channel
  
  // cpt.calc_adjacencies( hdr ,
  // 			col2var ,
  // 			col2f ,
  // 			col2t,
  // 			col2ch1 ,
  // 			col2ch2 ,
  // 			freq_threshold ,
  // 			0 ,    // time_threshold ,
  // 			NULL , // no clocs
  // 			0 ,    // no spatial_threshold ,
  // 			true ); // verbose

  // //
  // // Run assoc.
  // //
  
  // cpt_results_t results = cpt.run( nreps , cl_threshold , ! one_sided_test , verbose );
  
  //   //
  // // Report results
  // //

  // for (int y=0; y<vname.size(); y++)
  //   {
  //     const std::string & var =  vname[y] ;
  //     writer.level( var , globals::var_strat );
  //     writer.value( "B"  , results.beta[ var ] );
  //     writer.value( "STAT"  , results.t[ var ] );
  //     writer.value( "PU" , results.emp[ var ] );
  //     writer.value( "PC" , results.emp_corrected[ var ] );      
  //     writer.value( "CLST" , results.inclst[ var ] ); // 0 if not in a cluster

  //   }

  // //
  // // cluster-based results
  // //
  
  // int cln = 0;
  // std::map<std::string,double>::const_iterator qq = results.cluster_emp.begin();
  // while ( qq != results.cluster_emp.end() )
  //   {
  //     const std::set<std::string> & members = results.cluster_members.find( qq->first )->second;
      
  //     writer.level( ++cln , globals::cluster_strat );
      
  //     writer.value( "SEED" , qq->first );
  //     writer.value( "P" ,  qq->second );
  //     writer.value( "N" , (int)members.size() );
      
  //     // members
  //     int memn = 0;
  //     std::set<std::string>::const_iterator mm = members.begin();
  //     while ( mm != members.end() )
  // 	{
  // 	  writer.level( ++memn , "M" );
  // 	  writer.value( "VAR" , *mm );
  // 	  ++mm;
  // 	}
  //     writer.unlevel( "M" );
  //     ++qq;
  //   }
  // writer.unlevel( globals::cluster_strat );

  

}

void simassoc_t::output()
{
  
  

}

void build_param_from_cmdline( param_t * param )
{
  
  while ( ! std::cin.eof() )
    {
      std::string x;
      std::cin >> x;      
      if ( std::cin.eof() ) break;
      if ( x == "" ) continue;
      param->parse( x ); 
    }

}
