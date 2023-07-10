
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

#include "dsp/sl.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "clocs/clocs.h"
#include "clocs/legendre_polynomial.h"

void dsptools::surface_laplacian_wrapper( edf_t & edf , param_t & param )
{

  //
  // requires some clocs
  //
  
  if ( ! edf.clocs.attached() ) 
    edf.clocs.set_default();
  
  //
  // parameters
  //

  int m = param.has("m") ? param.requires_int( "m" ) : 4 ;
  int order = param.has("order") ? param.requires_int("order") : 10;
  double lambda = param.has("lambda") ? param.requires_dbl( "lambda" ) : 1e-5;

  //
  // get signals, dropping any non-data channels
  //

  std::string signal_label = param.requires( "sig" );  
  signal_list_t signals = edf.header.signal_list( signal_label );  
  edf.header.drop_annots_from_signal_list( &signals );
  const int ns = signals.size();

  
  //
  // no valid signals? quit
  //

  if ( ns == 0 ) 
    { 
      logger << "  no signals for SL, leaving\n";
      return; 
    }  
   
  //
  // check that all signals have the same sample rate
  //

  int sr = 0;
  for (int s=0;s<signals.size();s++)
    {
      if ( sr == 0 ) sr = edf.header.sampling_freq( signals(s) ) ;
      if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	Helper::halt( "requires all signals to have similar sampling rate, see RESAMPLE" );      
    }


  //
  // create G and H matrices given channel locations
  //
  

  sl_t sl( edf.clocs , signals , m , order , lambda );
  

  //
  // get signals
  //

  interval_t interval = edf.timeline.wholetrace();
  
  matslice_t mslice( edf , signals , interval );
  
  const Data::Matrix<double> & X = mslice.data_ref();


  //
  // apply to all signals
  //
  
  Data::Matrix<double> L;
  
  sl.apply( X , L );


  //
  // Update original signal
  //
  
  logger << "  updating with spatially-filtered signals\n";
  
  for (int s=0; s<signals.size(); s++)	
    edf.update_signal( signals(s) , L.col(s).data_pointer() );      
  

  //
  // all done
  //

}
 


sl_t::sl_t( const clocs_t & orig_clocs , const signal_list_t & signals , int m_ , int order_  , double lambda_ )
{
  
  //
  // default SL parameters
  //

  m = m_; // 4
  
  order = order_; // 20
   
  lambda = lambda_; // 1e-5

  
  
  //
  // copy clocs
  //

  clocs_t clocs = orig_clocs;
  
  const int ns = signals.size();

  
  //
  // scale channels to unit sphere
  //

  clocs.convert_to_unit_sphere();

//    std::map<std::string,cart_t>::const_iterator cc = clocs.cloc.begin();
//    while ( cc != clocs.cloc.end() )
//      {
//        std::cout << cc->first << "\t" << cc->second.x << "\t" << cc->second.y << "\t" <<cc->second.z << "\n";
//        ++cc;
//      }

  //
  // inter-electrode cosdistance matrix
  //
  

  Data::Matrix<double> D = clocs.interelectrode_distance_matrix( signals );


  //
  // Evaluate Legendre polynomials
  //
  
  std::vector<Data::Matrix<double> > L = legendre( order , D );


  //
  // precompute electrode-independent variables
  //

  std::vector<int> twoN1;
  std::vector<double> gdenom;
  std::vector<double> hdenom;

  for (int i=1;i<=order;i++) 
    { 
      twoN1.push_back( ( 2 * i ) + 1 ) ; 
      gdenom.push_back( pow( i*(i+1)  , m ) ) ;
      hdenom.push_back( pow( i*(i+1)  , m-1 ) ) ;

    }
  

  //
  // compute G and H
  //

  G.resize( ns , ns , 0 );
  H.resize( ns , ns , 0 );
  
  for (int i=0;i<ns;i++)
    for (int j=i;j<ns;j++)
      {
	double g = 0 , h = 0;
	for (int n=0;n<order;n++)
	  {
	    g += (twoN1[n] * L[n](i,j) ) / gdenom[n];
	    h -= (twoN1[n] * L[n](i,j) ) / hdenom[n];
	  }
	G(i,j) = g / ( 4.0 * M_PI );
	G(j,i) = G(i,j);

	H(i,j) = -h / ( 4.0 * M_PI );
	H(j,i) = H(i,j);
      }
  
    
  // 
  // Add lambda to each diagonal element
  //

  for (int i=0;i<ns;i++) G(i,i) = G(i,i) + lambda ; 
  
  //
  // Display
  //
  
  //std::cout << G.print() << "\n";

  //   std::cout << H.print( "H" , 8 , 8  ) << "\n";
  
  //
  // Inverse of G
  //

  bool okay = true;

  invG = Statistics::inverse( G , &okay );
  
  if ( ! okay ) Helper::halt( "problem in sl_t::sl_t() inverting G" );
  
  //  std::cout << std::fixed << std::setprecision(12) << invG.print() << "\n";
  // GsinvS = sum(inv(Gs));                                                                                                                                                                                                                                       
  GsinvS.resize( ns , 0 );
  sumGsinvS = 0;
  for (int i=0;i<ns;i++)
    for (int j=0;j<ns;j++)
      {
	GsinvS[j] += invG(i,j);
	sumGsinvS += invG(i,j);
      }
  //std::cout << "sumGsinvS = " << sumGsinvS << "\n";
 
  //  for (int i=0;i<ns;i++) std::cout << "GsinvS[j] " << GsinvS[i] << "\n";
 
}
  
bool sl_t::apply( const Data::Matrix<double> & data , Data::Matrix<double> & output )
{

  const int np = data.dim1();  

  const int ns = data.dim2();  

  logger << "  applying surface Laplacian for " << ns << " signals to " << np << " sample points\n";

  // dataGs = data'/Gs   [ ( np x ns )  =  (np x ns ) * ( ns x ns )
  // -->  data' * inv(Gs)

  Data::Matrix<double> dataGs( np , ns );
  for (int i=0;i<np;i++)
    for (int j=0;j<ns;j++)
      for (int k=0;k<ns;k++)
	dataGs(i,j) += data(i,k) * invG(k,j);

  // C = dataGs - (sum(dataGs,2)/sum(GsinvS))*GsinvS;
  std::vector<double> sumdataGs( np );
  for (int i=0;i<np;i++)
    {
      for(int j=0;j<ns;j++)
	sumdataGs[i] += dataGs(i,j);
      sumdataGs[i] /= sumGsinvS;
    }

  // sum(dataGs,2) is vector length(tp)
  //  np x 1  *  1 x ns   
  Data::Matrix<double> C( np , ns );
  for (int i=0;i<np;i++)
    for(int j=0;j<ns;j++)
      C(i,j) = dataGs(i,j) - sumdataGs[i] * GsinvS[j];

  
  // (C*H')'
  output.resize( np , ns );
  for (int i=0;i<np;i++)
    for(int j=0;j<ns;j++)
      for (int k=0;k<ns;k++)
	output(i,j) += C(i,k) * H(k,j);
    
  return true;
}


