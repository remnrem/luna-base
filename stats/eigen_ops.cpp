
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

#include "stats/eigen_ops.h"
#include "miscmath/miscmath.h"
#include "stats/statistics.h"
#include <vector>
#include "miscmath/crandom.h"


// nb. using Eigen:::Ref<>
// for a writable reference:    Eigen::Ref<Eigen::VectorXd> 
// for a const ref:             const Eigen::Ref<const Eigen::VectorXd> & 

std::vector<double> eigen_ops::copy_vector( const Eigen::VectorXd & e ) 
{
  std::vector<double> v( &e[0] , e.data() + e.size() );
  return v;
}

std::vector<double> eigen_ops::copy_array( const Eigen::ArrayXd & e ) 
{
  std::vector<double> v( &e[0] , e.data() + e.size() );
  return v;
}

Eigen::ArrayXd eigen_ops::copy_array( const std::vector<double> & e ) 
{
  Eigen::ArrayXd v = Eigen::ArrayXd::Zero( e.size() );
  for ( int i=0; i<e.size(); i++) v[i] = e[i];
  return v;
}

void eigen_ops::random_normal( Eigen::MatrixXd & M )
{
  const int rows = M.rows();
  const int cols = M.cols();
  for (int r = 0 ; r < rows ; r++ )
    for (int c = 0 ; c < cols ; c++)
      M(r,c) = Statistics::ltqnorm( CRandom::rand() );
}


bool eigen_ops::detrend( Eigen::Ref<Eigen::MatrixXd> M )
{

  // linear regression by hand for univariate case (detrending)
  const int n = M.rows();
  const int c = M.cols();
  
  // sum of 1.. n = n(n+1) / 2
  const double pmean = ( n + 1 ) / 2.0 ;
  Eigen::ArrayXd p( n , 1 );
  for (int i=0; i<n; i++) 
    p[i] = (i+1) - pmean;
  
  double pvar = p.square().sum()/(double)(n-1);

  for (int j=0; j<c; j++)
    {
      
      // mean center DV 
      const double intercept =  M.col(j).mean();
      Eigen::ArrayXd y = M.col(j).array() - intercept;
      
      // just need slope       
      double beta = ( ( y * p ).sum()/(double)(n-1) ) / pvar ; 

      // get residuals
      M.col(j) = y - ( p * beta );
      
    }

  return true;
}


bool eigen_ops::scale( Eigen::Ref<Eigen::MatrixXd> M , const bool center , const bool normalize ,
		       const bool ignore_invariants , std::vector<int> * zeros )
{

  if ( ! ( center || normalize ) ) return true;
  
  const int N = M.rows();
  
  Eigen::Array<double, 1, Eigen::Dynamic> means = M.colwise().mean();

  if ( normalize )
    {
      Eigen::Array<double, 1, Eigen::Dynamic> sds = ((M.array().rowwise() - means ).square().colwise().sum()/(N-1)).sqrt();

      for (int i=0;i<sds.size();i++) 
       	if ( sds[i] == 0 )
	  {
	    if ( ! ignore_invariants ) 
	      return false;
	    if ( zeros != NULL )
	      zeros->push_back( i );
	    sds[i] = 1.0; // make harmless
	  }
      
      if ( center ) 
	M.array().rowwise() -= means;
      M.array().rowwise() /= sds;
    }
  else
    {
      M.array().rowwise() -= means;
    }
  
  return true;
}


bool eigen_ops::robust_scale( Eigen::Ref<Eigen::MatrixXd> m , const bool center , bool normalize , double w , bool second_rescale ,
			      const bool ignore_invariants , std::vector<int> * zeros )
{
  // 1) winsorize at +/- w 

  const int rows = m.rows();
  const int cols = m.cols();
  for (int c=0;c<cols;c++)
    {

      std::vector<double> v = copy_vector( m.col(c) );
      double median = center ? MiscMath::median( v ) : 0 ;
      
      double iqr = normalize ? MiscMath::iqr( v ) : 0 ;

      // if no variation, set SD to one
      if ( normalize && iqr <= 1e-8  )
	{
	  normalize = false;
	  if ( ! ignore_invariants ) return false;
	  if ( zeros != NULL ) zeros->push_back( c );
	}
      double robust_sd = normalize ? 0.7413 * iqr : 1 ; 
      
      // winsorize?

      if ( w > 0 )
	{
	  double lwr = MiscMath::percentile( v , w );
	  double upr = MiscMath::percentile( v , 1-w );
	  if ( lwr >= upr ) 
	    Helper::halt( "cannot robust_scale().. pls fix me" );

	  for (int i=0; i<rows; i++)
	    {	      
	      if      ( m(i,c) < lwr ) m(i,c) = lwr;
	      else if ( m(i,c) > upr ) m(i,c) = upr;
	    }
	}

      // median / IQR normalize
      if ( center && normalize )  
	{
	  for (int i=0; i<rows; i++)	
	    m(i,c) = ( m(i,c) - median ) / robust_sd;
	}
      else if ( normalize )
	{
	  for (int i=0; i<rows; i++)	
	    m(i,c) = m(i,c) / robust_sd;
	}
      else if ( center )
	{
	  for (int i=0; i<rows; i++)	
	    m(i,c) = m(i,c) - median ;
	}
    }

  
  // finally, also scale by mean/variance too, just to ensure correct 
  // overall scale
  // hmm... unnecessary?  
  
  bool okay = true;
  
  if ( second_rescale ) 
    okay = scale( m , center , normalize , ignore_invariants );
  
  return okay;
}



double eigen_ops::sdev( const Eigen::VectorXd & x )
{
  const int N = x.size();
  double mean = x.mean();
  return sqrt( ( x.array() - mean ).square().sum()/(N-1) );
}



// unit_scale, specifying min/max and truncating at 0 and 1
Eigen::VectorXd eigen_ops::unit_scale( const Eigen::VectorXd & x , double xmin , double xmax )
{
  const int n = x.size();
  if ( n == 0 ) return x;
  if ( xmin >= xmax ) return x;
  
  Eigen::VectorXd r( n );
  for (int i=0;i<n;i++) 
    {
      if ( x[i] <= xmin ) r[i] = 0;
      else if ( x[i] >= xmax ) r[i] = 1;
      else r[i] = ( x[i] - xmin ) / ( xmax - xmin );
    }
  return r;

}


Eigen::VectorXd eigen_ops::unit_scale( const Eigen::VectorXd & x )
{

  const int n = x.size();
  if ( n == 0 ) return x;

  double xmin = x[0] , xmax = x[0];
  for (int i=0;i<n;i++)
    {
      if ( x[i] < xmin ) xmin = x[i];
      else if ( x[i] > xmax ) xmax = x[i];
    }  

  if ( xmin == xmax ) return x;

  Eigen::VectorXd r( n );
  for (int i=0;i<n;i++) r[i] = ( x[i] - xmin ) / ( xmax - xmin );
  return r;
}


Eigen::VectorXd eigen_ops::moving_average( const Eigen::VectorXd & x , int s )
{
  
  if ( s == 1 ) return x;

  const int n = x.size();

  if ( n == 0 ) return x;

  if ( s >= n ) 
    {
      std::cerr << "warning: in moving_average(), vector size is less than window size\n";
      s = n-1; 
      if ( s % 2 == 0 ) --s; // check that it remains odd
      if ( s < 2 ) return x; // bail out
    }

  if ( s % 2 == 0 ) Helper::halt( "require an odd-number for moving average" );

  double z = 0;
  
  const int edge = (s-1)/2;  
  const int start = edge;
  const int stop  = n - edge - 1;

  Eigen::VectorXd a = Eigen::VectorXd::Zero( n ) ;
  const double fac = 1.0/(double)s;
  for (int i=0;i<n;i++) a[i] = fac;
  
  // accumulate first sum
  for (int i=0;i<s;i++) z += x[i];

  // the main sets
  for (int i=start; i<=stop; i++)
    {
      a[i] *= z;
      if ( i == stop ) break;
      z -= x[i-edge];
      z += x[i+edge+1];      
    }

  // fill in at ends  
  for (int i=0;i<start;i++) a[i] = a[start];
  for (int i=stop+1;i<n;i++) a[i] = a[stop];
  return a;
  
}


Eigen::VectorXd eigen_ops::median_filter( const Eigen::VectorXd & x , const int n )
{

  bool odd = n % 2 ; 
  
  // For N odd, Y(k) is the median of X( k-(N-1)/2 : k+(N-1)/2 ).
  // For N even, Y(k) is the median of X( k-N/2 : k+N/2-1 ).
  
  const int t = x.size();

  Eigen::VectorXd ret( t );

  int v1 = odd ? (n-1)/2 : n/2;
  int v2 = odd ? (n-1)/2 : n/2-1;
  
  for (int i = 0 ; i < t ; i++ ) 
    {
      std::vector<double> y(n,0);
      int cnt = 0;
      for ( int j = i - v1 ; j <= i + v2 ; j++ )
	if ( j >= 0 && j < t ) y[cnt++] = x[j] ;
      
      // get median
      ret[i] = median_destroy( &y[0] , cnt );
      
    }
  
  return ret;
  
}


std::map<int,std::vector<double> > eigen_ops::group_means( const Eigen::MatrixXd & x , const std::vector<int> & g )
{
  std::map<int,std::vector<double> > m;
  std::map<int,int> c;

  const int n = g.size();
  if ( n != x.rows() )
    Helper::halt( "bad inputs to Statistics::group_means()" );

  if ( n == 0 )
    Helper::halt( "empty Statistics::group_means()" );

  const int p = x.cols();

  // initialize
  std::vector<double> t(p,0);
  for (int i=0; i<n; i++)
    if ( m.find(g[i]) == m.end() )
      m[g[i]] = t;

  // count
  for (int i=0; i<n; i++)
    {
      c[g[i]]++;
      std::vector<double> t(p,0);
      for (int j=0; j<p; j++)
	m[g[i]][j] += x(i,j);      
    }

  std::map<int,std::vector<double> >::iterator gg = m.begin();
  while ( gg != m.end() )
    {
      for (int j=0; j<p; j++)
	gg->second[j] /= (double)c[ gg->first ];
      ++gg;
    }
  
  return m;

}





// // apply function fx() with parameter param, to each matrix element

// void eigen_ops::apply_fx( Eigen::MatrixXd & M, double (*fx)(double,double), double param)
// {
//   const int rows = M.rows();
//   const int cols = M.cols();  
//   for (int i=0; i<rows; i++)
//     for (int j=0; j<cols; j++)
//       M[i][j] = (*fx)(M[i][j], param);
// }


// // row means
// void eigen_ops::mat_mean_rows( Eigen::MatrixXd & M, Eigen::ArrayXd & v)
// {

//   const int rows = M.dim1();
//   const int cols = M.dim2();
  
//   double sum;
  
//   for (int i=0; i< rows; i++) {
//     sum = 0;
//     for (int j=0; j<cols; j++)
//       sum += M[i][j];
//     v[i] = sum / cols;
//   }

// }

//  // * Returns the maximal element on the diagonal
//  // * of the matrix M.

// double eigen_ops::mat_max_diag( Eigen::MatrixXd & M )
// {

//   const int rows = M.dim1();
//   double max = M[0][0];
//   for (int i=1; i<rows; i++)
//     if (M[i][i] > max)
//       max = M[i][i];
//   return max;
// }

//  // Returns the maximal absolute element on the diagonal
//  // of the matrix M.
 
// double eigen_ops::mat_max_abs_diag(Eigen::MatrixXd & M )
// {
//   const int rows = M.dim1();
//   double max = fabs( M[0][0] ) ;
//   for (int i=1; i<rows; i++)
//     if ( fabs( M[i][i] ) > max)
//       max = fabs( M[i][i] );
//   return max;
// }

// // Creates a diagonal matrix from vector v.
 
// void eigen_ops::mat_diag( Eigen::ArrayXd & v , Eigen::MatrixXd & R )
// {
//   const int n = v.size();
//   mat_zeroize( R );
//   for (int i=0; i<n; i++)
//     R[i][i] = v[i];
// }

// // Transponse matrix M.

// void eigen_ops::mat_transpose( Eigen::MatrixXd & M , Eigen::MatrixXd & R )
// {
//   const int rows = M.dim1();
//   const int cols = M.dim2();
//   for (int i=0; i<rows; i++)
//     for(int j=0; j<cols; j++)
//       R[j][i] = M[i][j]; 
// }

// // Centers mat M. (Subtracts the mean from every column)

// void eigen_ops::mat_center( Eigen::MatrixXd & M )
// {
//   //
// }

// void eigen_ops::mat_center( Eigen::MatrixXd & M , Eigen::ArrayXd & means )
// {

//   // int rows = M.dim1();
//   // int cols = M.dim2();
  
//   // vect_zeroize( means, cols );
  
//   // for (int i=0; i<rows; i++)
//   //   for(int j=0; j<cols; j++)
//   //     means[j] += M[i][j];		
//   // for (int i=0; i<cols; i++)
//   //   means[i] /= rows; 
//   // for (int i=0; i<rows; i++)
//   //   for(int j=0; j<cols; j++)
//   //     M[i][j] -= means[j];	
// }


// void eigen_op::mat_decenter( Eigen::MatrixXd & M , Eigen::ArrayXd &  means )
// {

//   const int rows = M.dim1();
//   const int cols = M.dim2();
  
//   for(int i=0; i<rows; i++)
//     for(int j=0; j<cols; j++)
//       M[i][j] += means[j]; 
// }




double eigen_ops::between_within_group_variance( const std::vector<std::string> & g , const Eigen::VectorXd & x )
{
  // this routine is ONLY to be used by SUDS as a special case.
  // note - here we know 'x' is already standardized.   Just get the max within-class variance and return
  // i.e. want to check it is not (e.g. ) more than twece the total variance (implies too)
  // note - function name is now misleading.... change in future.. this is simpoly the max within-group variance
  // (given X is normalized overall)
  
  const int n = x.size();
  const double grand_mean = x.sum() / (double)n;
  const double grand_sumsq = x.array().square().sum();
  
  // group means/ sumsqs  
  std::map<std::string,double> group_s;
  std::map<std::string,int> group_n;
  std::map<std::string,double> group_sumsq;;
  std::map<std::string,double> group_mean;
  for (int i=0; i<n; i++)
    {
      group_s[ g[i] ] += x[i];
      group_sumsq[ g[i] ] += x[i] * x[i];
      group_n[ g[i] ]++;      
    }

  const int ng = group_n.size();
  
  // not enough valid groups
  if ( ng < 2 ) return 0;

  std::map<std::string,double>::iterator gg = group_s.begin();
  while ( gg != group_s.end() )
    {      
      group_mean[ gg->first ] = group_s[ gg->first ] / (double) group_n[ gg->first ] ;      
      ++gg;
    }
        
  //
  // each group variance
  //
  
  double wmax = 0;
  
  gg = group_sumsq.begin();
  while ( gg != group_sumsq.end() )
    {
      if ( group_n[ gg->first ] >= 2 )
	{
	  double within_variance = gg->second - group_n[ gg->first ] * group_mean[ gg->first ] * group_mean[ gg->first ] ;
	  within_variance /= (double) group_n[ gg->first ] - 1 ; 	  
	  if ( within_variance > wmax ) wmax = within_variance ;
	}
      ++gg;
    }
  
  return wmax ;
    
}


