
// https://www.kaggle.com/code/jdarcy/introducing-ssa-for-time-series-decomposition
// https://ssa.cf.ac.uk/ssa2010/a_brief_introduction_to_ssa.pdf
// https://la.mathworks.com/matlabcentral/fileexchange/58968-multichannel-singular-spectrum-analysis-beginners-guide


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


#include "dsp/ssa.h"

#include <iostream>
#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "eval.h"
#include "db/db.h"

extern logger_t logger;
extern writer_t writer;

void dsptools::ssa_wrapper( edf_t & edf , param_t & param )
{
  
}


ssa_t::ssa_t( const std::vector<double> * x , const int l )
{
  
  const int n = x->size();
  
  // copy
  Eigen::VectorXd t = Eigen::VectorXd::Zero( n );
  for (int i=0; i<n; i++) t[i] = (*x)[i];

  fit( t , l );
}

ssa_t::ssa_t( const Eigen::VectorXd & t , const int l )
{
  fit( t , l );
}


void ssa_t::fit( const Eigen::VectorXd & t , const int l )
{
  
  // SSA: decomposes time series t using SSA, assuming equal intervals 

  // t : time series (length n)
  // l : window length ( 2 <= l <= n/2)

  const int n = t.size();
  
  if ( l < 2 || l > n/2 )
    Helper::halt( "window length l must be between 2 and n/2" );
  
  // save memory by not storing these
  const bool store_elementary_matrices = false;

  // columns in trajectory matrix
  const int k = n - l + 1;

  // allocate trajectory matrix X
  X = Eigen::MatrixXd::Zero( l , k );  
  for (int i=0; i<k; i++)
    X.col(i) = t.segment(i,l);
  
  // decompose X w/ SVD

  Eigen::BDCSVD<Eigen::MatrixXd> svd( X , Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXd U = svd.matrixU();
  Eigen::MatrixXd V = svd.matrixV();
  Eigen::MatrixXd W = svd.singularValues();

  // rank of X
  int d = svd.rank();

// # Decompose the trajectory matrix
//   self.U, self.Sigma, VT = np.linalg.svd(self.X)
//         self.d = np.linalg.matrix_rank(self.X)

  // components
  Eigen::MatrixXd TS_comps = Eigen::MatrixXd::Zero( n , d );
        
    //     if not save_mem:
    //         # Construct and save all the elementary matrices
    //         self.X_elem = np.array([ self.Sigma[i]*np.outer(self.U[:,i], VT[i,:]) for i in range(self.d) ])

    //         # Diagonally average the elementary matrices, store them as columns in array.           
    //         for i in range(self.d):
    //             X_rev = self.X_elem[i, ::-1]
    //             self.TS_comps[:,i] = [X_rev.diagonal(j).mean() for j in range(-X_rev.shape[0]+1, X_rev.shape[1])]
            
    //         self.V = VT.T
    //     else:
    //         # Reconstruct the elementary matrices without storing them
    //         for i in range(self.d):
    //             X_elem = self.Sigma[i]*np.outer(self.U[:,i], VT[i,:])
    //             X_rev = X_elem[::-1]
    //             self.TS_comps[:,i] = [X_rev.diagonal(j).mean() for j in range(-X_rev.shape[0]+1, X_rev.shape[1])]
            
    //         self.X_elem = "Re-run with save_mem=False to retain the elementary matrices."
            
    //         # The V array may also be very large under these circumstances, so we won't keep it.
    //         self.V = "Re-run with save_mem=False to retain the V matrix."
        
    //     # Calculate the w-correlation matrix.
    //     self.calc_wcorr()
            
    // def components_to_df(self, n=0):
    //     """
    //     Returns all the time series components in a single Pandas DataFrame object.
    //     """
    //     if n > 0:
    //         n = min(n, self.d)
    //     else:
    //         n = self.d
        
    //     # Create list of columns - call them F0, F1, F2, ...
    //     cols = ["F{}".format(i) for i in range(n)]
    //     return pd.DataFrame(self.TS_comps[:, :n], columns=cols, index=self.orig_TS.index)
            
    
    // def reconstruct(self, indices):
    //     """
    //     Reconstructs the time series from its elementary components, using the given indices. Returns a Pandas Series
    //     object with the reconstructed time series.
        
    //     Parameters
    //     ----------
    //     indices: An integer, list of integers or slice(n,m) object, representing the elementary components to sum.
    //     """
    //     if isinstance(indices, int): indices = [indices]
        
    //     ts_vals = self.TS_comps[:,indices].sum(axis=1)
    //     return pd.Series(ts_vals, index=self.orig_TS.index)
    
    // def calc_wcorr(self):
    //     """
    //     Calculates the w-correlation matrix for the time series.
    //     """
             
    //     # Calculate the weights
    //     w = np.array(list(np.arange(self.L)+1) + [self.L]*(self.K-self.L-1) + list(np.arange(self.L)+1)[::-1])
        
    //     def w_inner(F_i, F_j):
    //         return w.dot(F_i*F_j)
        
    //     # Calculated weighted norms, ||F_i||_w, then invert.
    //     F_wnorms = np.array([w_inner(self.TS_comps[:,i], self.TS_comps[:,i]) for i in range(self.d)])
    //     F_wnorms = F_wnorms**-0.5
        
    //     # Calculate Wcorr.
    //     self.Wcorr = np.identity(self.d)
    //     for i in range(self.d):
    //         for j in range(i+1,self.d):
    //             self.Wcorr[i,j] = abs(w_inner(self.TS_comps[:,i], self.TS_comps[:,j]) * F_wnorms[i] * F_wnorms[j])
    //             self.Wcorr[j,i] = self.Wcorr[i,j]
    
    // def plot_wcorr(self, min=None, max=None):
    //     """
    //     Plots the w-correlation matrix for the decomposed time series.
    //     """
    //     if min is None:
    //         min = 0
    //     if max is None:
    //         max = self.d
        
    //     if self.Wcorr is None:
    //         self.calc_wcorr()
        
    //     ax = plt.imshow(self.Wcorr)
    //     plt.xlabel(r"$\tilde{F}_i$")
    //     plt.ylabel(r"$\tilde{F}_j$")
    //     plt.colorbar(ax.colorbar, fraction=0.045)
    //     ax.colorbar.set_label("$W_{i,j}$")
    //     plt.clim(0,1)
        
    //     # For plotting purposes:
    //     if max == self.d:
    //         max_rnge = self.d-1
    //     else:
    //         max_rnge = max
        
    //     plt.xlim(min-0.5, max_rnge+0.5)
    //     plt.ylim(max_rnge+0.5, min-0.5)

	      


  
  std::cout << "X\n" << X << "\n";
  
}
