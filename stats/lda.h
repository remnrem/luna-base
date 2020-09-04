
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

#ifndef __LDA_H__
#define __LDA_H__

#include "matrix.h"
#include <vector>
#include <map>


struct lda_model_t {
  
  lda_model_t()
  {
    valid = false;
    errmsg = "";
  }

   
  bool valid;

  std::string errmsg;
  
  Data::Vector<double> prior;
  std::map<std::string,int> counts;
  Data::Matrix<double> means;
  Data::Matrix<double> scaling;
  int n;
  Data::Vector<double> svd;
  std::vector<std::string> labels;

  Data::Vector<double> prop_trace() const
  {
    Data::Vector<double> t( svd.size() );

    double sum = 0;

    for (int i=0;i<svd.size();i++)
      {
	t[i] = svd[i] *	svd[i];
	sum += svd[i] * svd[i];
      }

    for	(int i=0;i<svd.size();i++)
      t[i] /= sum;

    return t;
  }
  
};


struct lda_posteriors_t {

  // cols = classes
  // rows = observations
  
  Data::Matrix<double> pp;
  
  // most likely class label for each observation
  std::vector<std::string> cl;
  std::vector<int> cli; // as above, but index

  
};


class lda_t {
  
 public:
  
 lda_t( const std::vector<std::string> y ,
	const Data::Matrix<double> & X )
    : y(y) , X(X)
    {
      tol = 1e-4;
    } 
  
  lda_model_t fit();
  
  static lda_posteriors_t predict( const lda_model_t & , const Data::Matrix<double> & X );

 private:
  
  std::vector<std::string> y;

  Data::Matrix<double> X;

  double tol;
  
};

#endif
