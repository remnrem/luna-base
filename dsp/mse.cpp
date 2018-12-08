
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


#include "mse.h"

#include "miscmath/miscmath.h"

  
#include <iostream>

/* file: mse.c			M. Costa		1 August 2004
				Last revised:		4 August 2004 (GBM)
-------------------------------------------------------------------------------
mse: calculates multiscale entropy (MSE) of one or multiple data sets
Copyright (C) 2004 Madalena Costa

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place - Suite 330, Boston, MA 02111-1307, USA.

You may contact the author by e-mail (mcosta@fsa.harvard.edu).  For updates to
this software, please visit PhysioNet (http://www.physionet.org/).
_______________________________________________________________________________

Compile this program by
    gcc -o mse -O mse.c -lm

There are two major steps in the calculations performed by mse:
1. Time series are coarse-grained.
2. Sample entropy (SampEn) is calculated for each coarse-grained time series.

Output file:
1st line: shows the r value.
2nd line: shows the m values. 

After the 2nd line there are several columns: the first column (of integers)
is the scale factor. The following columns are SampEn values for coarse-grained
time series calculated for the values of r and m specified. If the option for
calculating MSE for several r values is chosen a new line containing the new r
value and new columns with the corresponding results are written.
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>


std::map<int,double> mse_t::calc( const std::vector<double> & d )
{
  
  std::map<int,double> retval;
  
  // first normalize input
  std::vector<double> zd = MiscMath::Z( d );
  
  const int n = d.size();
  
  // get SD for data
  //double sdev = SD( zd );
    
  // Iterate over each scale j
  for (int j = 1; j <= scale_max; j += scale_step)
    {      
      
      std::vector<double> y = coarse_graining( zd , j ) ;
            
      // faster version (from mse.c)
      //retval[j] = sample_entropy( y , 1.0 );
      
      // old version (slower)
       retval[j] = sampen( y , m , r );
      
    }
  
  return retval;

}


double mse_t::SD(const std::vector<double> & x)
{
  double sum=0.0, sum2=0.0, sd;
  int j;
  
  const int nlin = x.size();
  
  for (int j = 0; j < nlin; j++) 
    {
      sum += x[j];
      sum2 += x[j] * x[j];
    }
  sd = sqrt((sum2 - sum*sum/nlin)/(nlin - 1));
  return sd;
}


std::vector<double> mse_t::coarse_graining( const std::vector<double> & x , int j )
{
  
  const int n  = x.size();
  const int ns = n/j; 
  
  std::vector<double> y( ns , 0 );
  
  for (int i = 0; i < ns; i++) 
    {      
      for (int k = 0; k < j; k++) y[i] += x[i*j+k];
      y[i] /= j; 
    }
  return y;
}


double mse_t::sample_entropy( const std::vector<double> & y , double sd )
{
  
  // sd   std dev
  
  // m    pattern length parameter  (typically 2)
  // r    pattern match (SD units, i.e. typically SD == 1 ) (typically 0.15)
  
  std::vector<int> cont( m+2 , 0 );  
  
  const int nlin_j = y.size() - m; 
   
  double r_new = r * sd;

  for (int i = 0; i < nlin_j; ++i) 
    {
      for (int l = i+1; l < nlin_j; ++l) 
 	{ /*self-matches are not counted*/
 	  int k = 0;
 	  while (k < m && fabs(y[i+k] - y[l+k]) <= r_new)
 	    cont[++k]++;
 	  if (k == m && fabs(y[i+m] - y[l+m]) <= r_new)
 	    cont[m+1]++;
 	} 
    }     
  
  //    for (i = 1; i <= m; i++)
  //      if (cont[i] == 0 || cont[i-1] == 0)
  //        SE[j][i] = -log((double)1/((nlin_j)*(nlin_j-1)));
  //      else
  //        SE[j][i] = -log((double)cont[i+1]/cont[i]);
  
  if ( cont[m+1] == 0 || cont[m] == 0 ) return -1;

  return -log( cont[m+1]/(double)cont[m] );

  //   std::cout << "v2 = " << cont[m+1] << " / " << cont[m] << " " << -log( cont[m+1]/(double)cont[m] ) <<"\n";
}



// sampen() calculates an estimate of sample entropy 

double mse_t::sampen( const std::vector<double> y , int M , double r )
{
  
  const int n = y.size();
    
  ++M;
  
  std::vector<long> run( n );
  std::vector<long> lastrun( n );
  std::vector<double> A(M);
  std::vector<double> B(M);
  std::vector<double> p(M);

  /* start running */
  for (int i = 0; i < n - 1; i++) 
    {
      int nj = n - i - 1;
      double y1 = y[i];
      for (int jj = 0; jj < nj; jj++) 
	{
	  int j = jj + i + 1;
	  if (((y[j] - y1) < r) && ((y1 - y[j]) < r)) 
	    {
	      run[jj] = lastrun[jj] + 1;
	      int M1 = M < run[jj] ? M : run[jj];
	      for (int mi = 0; mi < M1; mi++) 
		{
		  A[mi]++;
		  if (j < n - 1)
		    B[mi]++;
		}
	    }
	  else
	    run[jj] = 0;
	}/* for jj */
      for (int j = 0; j < nj; j++)
	lastrun[j] = run[j];
    }/* for i */
  
  int N = (long) (n * (n - 1) / 2);
  p[0] = A[0] / N;
  //  printf("SampEn(0,%g,%d) = %lf\n", r, n, -log(p[0]));

  if ( false || 1 ) 
    {
      
      for (int mi = 1; mi < M; mi++) 
	{
	  p[mi] = A[mi] / B[mi - 1];
	  if (p[mi] == 0)
	    {
	      //	std::cout << "sampEn " << m << " " << r << " " << n << " " << "NA" << "\n";
	      //std::cerr << "sampEn ... no matches!\n";
	      //printf("No matches! SampEn((%d,%g,%d) = Inf!\n", m, r, n);
	    }
	  else
	    {
	      //std::cout << "sampEn " << mi << " " << r << " " << n << " " << -log(p[mi]) << "\n";
	      //      printf("SampEn(%d,%g,%d) = %lf\n", m, r, n, 
	    }
	}
    }

  const int mi = m;
  p[mi] = A[mi] / B[mi-1];

  //  std::cout << "final A/B " << A[mi] << " " << B[mi-1] << " (" << A[mi-1] << ") " << -log(p[mi]) << "\n";

  // return SE for m=2 


  if (p[mi] == 0) return -1;
  return -log(p[mi]);

}



