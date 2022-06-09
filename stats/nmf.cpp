
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

#include "stats/nmf.h"

#include "helper/logger.h"

extern logger_t logger;

// naive NMF implementation...

nmf_t::nmf_t( const Eigen::MatrixXd & V_ , const int maxiter , const double EPS )
  : V(V_) , maxiter(maxiter ) , EPS(EPS)
{
  // scale to positivity?
  const double m = V.minCoeff();

  if ( m < 0 )
    {
      V.array() -= m;
    }
      
    
  // remove null rows
  rows.clear();
  const int nr = V.rows();
  included.resize( nr , true );
  
  for (int r=0;r<nr;r++)
    {
      if ( V.row(r).sum() >= EPS )
	rows.push_back( r );
      else
	included[r] = false;	
    }

  // need to splice out? 
  if ( rows.size() < nr )
    {
      Eigen::MatrixXd V2 = V;
      V.resize( rows.size() , Eigen::NoChange );
      for (int r=0; r<rows.size(); r++)
	V.row(r) = V2.row(rows[r]);      
      logger << " spliced out " << V.rows() << " from " << V_.rows() << "\n";
    }
}

void nmf_t::factorize( const int num_sources )
{
  // dimensions
  const int xs = V.cols();
  const int ys = V.rows();
  
  // randomize initial values for H and W
  H = Eigen::MatrixXd::Zero(num_sources , xs );
  H.setRandom();
  H = H.array().abs();
  
  W = Eigen::MatrixXd::Zero( ys, num_sources );
  W.setRandom();
  W = W.array().abs();
  
  Eigen::MatrixXd ones = Eigen::MatrixXd::Constant( xs, ys, 1.0 );
  
  // iterate  
  iter = 0;

  for(int i = 0; i < maxiter ; i++)
    {
      // update W
      Eigen::MatrixXd WH = (W * H).array() + EPS ;      
      Eigen::MatrixXd numerator = V.array() / WH.array() ;      
      numerator = numerator * H.transpose();  
      Eigen::MatrixXd denominator = (H * ones).array() + EPS ;
      W = W.array() * numerator.array();
      W = W.array() / denominator.transpose().array();

      // update H
      WH = (W * H).array() + EPS;
      numerator = V.array() / WH.array();
      numerator = W.transpose() * numerator;
      denominator = (ones * W).array() + EPS;  
      H = H.array() * numerator.array();
      H = H.array() / denominator.transpose().array();       

      ++iter;
    }
  
}
