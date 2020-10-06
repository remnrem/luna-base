
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

#include "miscmath.h"
#include "helper/helper.h"
#include "stats/statistics.h"

#include <cmath>
#include <vector>

#include <iostream>

//
// Find next largest power of 2
//

long int MiscMath::nextpow2( const int a )
{
  // for now, just go up to 2^31
  for (int i=1;i<32;i++)
    {
      long int t = pow(2,i);
      if ( a <= t ) return t;
    }
  Helper::halt("value too large in nextpow2()");
  return 0;
}


//
// logspace function
//

std::vector<double> MiscMath::logspace(double a, double b, int n)
{
  if ( n < 2 ) Helper::halt( "logspace requires at least two values" );
  // for a->b, generate a series of logarithmically (log10) equally spaced intervals
  a = log10(a);
  b = log10(b);
  const double st = (b-a)/(double)(n-1);
  std::vector<double> r(n);
  r[0] = pow(10,a);
  r[n-1] = pow(10,b);
  for (int i=1;i<n-1;i++) r[i] = pow( 10, a + i*st );
  return r;
}


//
// linspace function
//

std::vector<double> MiscMath::linspace(double a, double b, int n)
{
  if ( n < 2 ) Helper::halt( "linspace requires at least two values" );
  const double st = (b-a)/(double)(n-1);
  std::vector<double> r(n);
  r[0] = a; r[n-1] = b;
  for (int i=1;i<n-1;i++) r[i] = a + i*st ;
  return r;
}

//
// Differences 
//

std::vector<double> MiscMath::diff( const std::vector<double> & x )
{
  const int n = x.size();
  if ( n < 2 ) 
    Helper::halt( "problem in diff() -- input less than two elements" );  
  std::vector<double> r( n - 1 );
  for (int i=1;i<n;i++)
    r[i-1] = x[i] - x[i-1];
  return r;
}

//
// Proportion of clipped points
//

double MiscMath::clipped( const std::vector<double> & x , double mn , double mx )
{
  
  const double rng = mx - mn;
  
  // Use a fixed tolerance, if range if ~0. Basically just intereted
  // in the case where everything is flat (in which case return 1.0 as 
  // completely 'clipped')
  
  if ( rng < 1e-12 ) return 1.0;
  
  const double tol = rng * 0.0001;
  
  int c = 0;
  
  const int n = x.size();
  for (int i=0;i<n;i++)
    {
      if ( fabs( x[i] - mx ) < tol ) ++c;
      if ( fabs( x[i] - mn ) < tol ) ++c;
    }

  // allow for the two 'index' points (i.e. identified min/max)
  c -= 2;
  if ( c < 0 ) c = 0;
  return c / (double)(n-2);
}

double MiscMath::clipped( const std::vector<double> & x )
{
  const int n = x.size();
  double mx = 0, mn = 0;
  
  for (int i=0;i<n;i++)
    {
      if ( x[i] > mx ) mx = x[i];
      if ( x[i] < mn ) mn = x[i];
    }
  return clipped( x , mn , mx );
}

double MiscMath::flat( const std::vector<double> & x , double EPS )
{
  const int n = x.size();
  int c = 0;
  for (int i=1;i<n;i++) if ( fabs( x[i] - x[i-1] ) < EPS ) ++c;
  return c/(double)(n-1);  
}


double MiscMath::max( const std::vector<double> & x , double th )
{
  const int n = x.size();
  int c = 0;
  for (int i=0;i<n;i++) if ( fabs( x[i] ) > th ) ++c;
  return c/(double)(n);  
}



double MiscMath::sqr( const double a )
{
  return a*a;
}

//
// Mean, variance
//

double MiscMath::sum( const std::vector<double> & x )
{
  const int n = x.size();
  if ( n == 0 ) return 0; // silently fail here
  double s = 0;
  for (int i=0;i<n;i++) s += x[i];
  return s;
}

double MiscMath::mean( const std::vector<double> & x )
{
  const int n = x.size();
  if ( n == 0 ) return 0; // silently fail here
  double s = 0;
  for (int i=0;i<n;i++) s += x[i];
  return s/(double)n;
}

double MiscMath::mean( const std::vector<int> & x )
{
  const int n = x.size();
  if ( n == 0 ) return 0; // silently fail here
  double s = 0;
  for (int i=0;i<n;i++) s += x[i];
  return s/(double)n;
}

std::complex<double> MiscMath::max( const std::vector<std::complex<double> > & x )
{
  const int n = x.size();
  if ( n == 0 ) return std::complex<double>(0,0);
  double mm = 0;
  int mi = 0;
  for (int i=0;i<n;i++)
    {
      double m = std::abs( x[i] );
      if ( m > mm ) { mi = i; mm = m; }
    }
  return x[mi];
}


std::complex<double> MiscMath::mean( const std::vector<std::complex<double> > & x )
{
  const int n = x.size();
  if ( n == 0 ) return std::complex<double>(0,0);

  double rl = 0 , img = 0;
  for (int i=0;i<n;i++)
    {
      rl += std::real(x[i]);
      img += std::imag(x[i]);
    }
  rl /= (double)n;
  img /= (double)n;
  return std::complex<double>( rl , img );
}


double MiscMath::skewness( const std::vector<double> & x )
{
  double m    = MiscMath::mean( x );
  double sd   = MiscMath::sdev( x , m );
  return MiscMath::skewness( x , m , sd );
}

double MiscMath::skewness( const std::vector<double> & x , double m , double sd )
{
  double sum = 0;
  const int n = x.size();
  for (int i=0; i<n; i++)
    sum += ( x[i] - m ) * ( x[i] - m ) * ( x[i] - m ) ;
  return sum / (double)( n * sd * sd * sd );  
} 
  

double MiscMath::median( const std::vector<double> & x )
{
  //  const double * a = &(x[0]);
  const int n = x.size();  
  //  return median_preserve( p,n );
  if ( n == 0 ) Helper::halt( "internal problem, taking median of 0 elements");
  if ( n == 1 ) return x[0];
  return MiscMath::kth_smallest_preserve(x,(((n)&1)?((n)/2):(((n)/2)-1)));
}

double MiscMath::iqr( const std::vector<double> & x )
{
  std::vector<double> q(2);
  q[0] = 0.25; q[1] = 0.75;
  std::vector<double> quartiles = MiscMath::quantile<double>( x , q );
  return quartiles[1] - quartiles[0];
}

double MiscMath::percentile( const std::vector<double> & x , double p )
{
  
  const int n = x.size();  
  if ( n == 0 ) Helper::halt( "internal problem, taking percentile of 0 elements");
  if ( n == 1 ) return x[0];
  if ( p < 0 || p > 1 ) Helper::halt( "internal problem, invalid percentile specified" );
  int pn = n*p;
  return MiscMath::kth_smallest_preserve(x,pn);
}


double MiscMath::meansq( const std::vector<double> & x )
{
  const int n = x.size();
  if ( n == 0 ) return 0;
  double s = 0;
  for (int i=0;i<n;i++) s += x[i]*x[i];
  return s/(double)n;
}

double MiscMath::variance( const std::vector<double> & x )
{
  return variance( x , mean(x) );
}

double MiscMath::variance( const std::vector<int> & x )
{
  return variance( x , mean(x) );
}

double MiscMath::variance( const std::vector<double> & x , double m )
{
  const int n = x.size();
  if ( n == 0 ) return 0;
  double ss = 0;
  for (int i=0;i<n;i++)
    {
      const double t = x[i] - m;
      ss += t*t;      
    }
  return ss/(double)(n-1);
}

double MiscMath::variance( const std::vector<int> & x , double m )
{
  const int n = x.size();
  if ( n == 0 ) return 0;
  double ss = 0;
  for (int i=0;i<n;i++)
    {
      const double t = x[i] - m;
      ss += t*t;      
    }
  return ss/(double)(n-1);
}

double MiscMath::sdev( const std::vector<double> & x )
{
  return sqrt( variance(x) );
}

double MiscMath::sdev( const std::vector<double> & x , double m )
{
  return sqrt( variance( x , m ) );
}


//
// Hjorth parameters
//

void MiscMath::hjorth( const std::vector<double> * data , double * activity , double * mobility , double * complexity )
{
  
  if ( activity == NULL || data == NULL || mobility == NULL || complexity == NULL ) 
    Helper::halt( "NULL given to hjorth()" );
  
  const int n = data->size();
  if ( n == 0 ) 
    {
      *activity   = 0;
      *complexity = 0;
      *mobility   = 0;
      return;
    }
  
  std::vector<double> dxV = diff( *data );
  std::vector<double> ddxV = diff( dxV );  
  
  double mx2 = meansq( *data );
  double mdx2 = meansq( dxV );
  double mddx2 = meansq( ddxV );
  
  *activity   = mx2;
  *mobility   = mdx2 / mx2;
  *complexity = sqrt( mddx2 / mdx2 - *mobility );
  *mobility   = sqrt( *mobility );
  
  if ( ! Helper::realnum( *activity ) ) *activity = 0;
  if ( ! Helper::realnum( *mobility ) ) *mobility = 0;
  if ( ! Helper::realnum( *complexity ) ) *complexity = 0;
  
}

//
// Turning rate
//

double MiscMath::turning_rate( const std::vector<double> * d , int sr , int es , int trd , std::vector<double> * sub )
{

  // Epoch size in seconds
  const int nt = d->size() / sr;

  // number of sub-epochs (may truncate some data)
  const int ne = nt / es;

  // sample points per sub-epoch
  const int le = es * sr;

  int p = 0;
  
  std::vector<double> stored;

  double mean = 0;

  for (int e=0;e<ne;e++)
    {
      
      // extract data points
      std::vector<double> extract;

      for (int j=trd; j<le-trd; j++)
	{
	  int pj = p + j;
	  //std::cerr << " pj = " << p << " " << j << " " << pj << " " << d->size() << "\n";
	  if ( extract.size() == 0 ) extract.push_back( (*d)[p+j] );
	  else if ( extract[j-1] != (*d)[p+j] ) extract.push_back( (*d)[p+j] );
	}

      int turns = 0; 
      const int n = extract.size();
      for (int j=1;j<n-1;j++)
	{
	  if      ( extract[j-1] > extract[j] && extract[j] < extract[j+1] ) ++turns;
	  else if ( extract[j-1] < extract[j] && extract[j] > extract[j+1] ) ++turns;
	}
      

      const double tr = turns / (double)(n-2);
      // store for sub-epoch (e.g. 1 sec)
      stored.push_back( tr);

      // get mean over all sub-epochs (i.e. typically the 30-second value)
      mean += tr;
      
      // skip whole sub-epoch
      p += le;
    }
  
  
  // report back original sub-epoch level results
  if ( sub != NULL ) *sub = stored;
  
  return mean / (double)ne;
  
}



//
// Root mean square function
//

double MiscMath::rms( const std::vector<double> & x)
{
  const int N = x.size();
  double d = 0;  
  for (int i=0;i<N;i++) d += x[i] * x[i];
  return sqrt( d / (double)N );
}

double MiscMath::rms( const double * x , const int N)
{
  double d = 0;  
  for (int i=0;i<N;i++) 
    {
      d += (*x) * (*x);
      ++x;
    }
  return sqrt( d / (double)N );
}


//
// Tukey window
//


std::vector<double> MiscMath::tukey_window( int n , double r )
{
  
  double step = 1.0/(double)(n-1);
  double rhalf = r/2.0;
  std::vector<double> w(n);
  
  for (int i=0;i<n;i++)
    {
      double x = i * step;
      if ( x < rhalf ) 
	{
	  w[i] = 0.5 * ( 1 + cos(  ( (  2 * M_PI ) / r ) * ( x - rhalf )   ) );
	}
      else if ( x >= 1 - rhalf )
	{
	  w[i] = 0.5 * ( 1 + cos(  ( (  2 * M_PI ) / r ) * ( x - 1 + rhalf )  ) );
	}
      else
	{	  
	  w[i] = 1;
	}

    }
  return w;
}

void MiscMath::tukey_window( std::vector<double> * d , double r )
{
  std::vector<double> w = tukey_window( d->size() , r );
  for (int i=0;i<d->size();i++) (*d)[i] = (*d)[i] * w[i];
}


//
// Hann window
//

double MiscMath::hann_window(unsigned int i, unsigned int N)
{ 
  return 0.5*(1 - cos(2.0*M_PI*i/(double)(N - 1)));
}

void MiscMath::hann_window( std::vector<double> * d)
{
  std::vector<double> w = hann_window( d->size() );
  for (int i=0;i<d->size();i++) (*d)[i] = (*d)[i] * w[i];
}

std::vector<double> MiscMath::hann_window(int N )
{
  std::vector<double> w(N);
  for (int i = 0; i < N; i++)
    w[i] = 0.5*(1 - cos(2.0*M_PI*i/(double)(N - 1)));
  return w;
}


//
// Hanning  (same as Hann(n-2) w/ 0 at start/end
// i.e. this matches matlab hanning() function 
//

std::vector<double> MiscMath::hanning_window(int N )
{
  if ( N < 3 ) Helper::halt( "bad hanning window" );
  std::vector<double> r(N,0);
  std::vector<double> w = hann_window( N+2 );
  for (int i=0;i<N;i++)
    r[i] = w[i+1];
  return r;
}


//
// Hamming window function (symmetric)
//

void MiscMath::hamming_window( std::vector<double> * d)
{
  std::vector<double> w = hamming_window( d->size() );
  for (int i=0;i<d->size();i++) (*d)[i] = (*d)[i] * w[i];
}

std::vector<double> MiscMath::hamming_window( int N )
{
  std::vector<double> w(N);
  for (int n=0;n<N;n++) w[n] = hamming_window(n,N);
  return w;
}

double MiscMath::hamming_window(unsigned int n, unsigned int N) 
{ 
  // note, 1-based
  return 0.54f - 0.46f * std::cos( 2.0 * M_PI * ( (n) / (double)(N-1) ) ); 
}


//
// Z-score normalisation
//

std::vector<double> MiscMath::Z( const std::vector<double> & x )
{
  const int n = x.size();
  const double m = mean( x );
  const double sd   = sdev( x );
  if ( sd == 0 ) return x;
  std::vector<double> y( n );
  for (int i=0;i<n;i++) y[i] = ( x[i] - m ) / sd;
  return y;
}

//
// Logs
//

std::vector<double> MiscMath::logvector( const std::vector<double> & x )
{
  const int n = x.size();
  std::vector<double> y( n );
  for (int i=0;i<n;i++) y[i] = log( x[i] );
  return y;
}


//
// Wirth median function
//

MiscMath::elem_type MiscMath::kth_smallest_preserve( const std::vector<MiscMath::elem_type> & a , int k )
{
  std::vector<elem_type> cpy = a;
  return kth_smallest_destroy( &(cpy[0]), cpy.size(), k );
}

MiscMath::elem_type MiscMath::kth_smallest_destroy(MiscMath::elem_type a[], int n, int k)
{
  int i,j,l,m ;
  elem_type x ;

  l=0 ; m=n-1 ;
  while (l<m) {
    x=a[k] ;
    i=l ;
    j=m ;
    do {
      while (a[i]<x) i++ ;
      while (x<a[j]) j-- ;
      if (i<=j) {
	ELEM_SWAP(a[i],a[j]) ;
	i++ ; j-- ;
      }
    } while (i<=j) ;
    if (j<k) l=i ;
    if (k<i) m=j ;
  }
  return a[k] ;
}



//
// Epoch functions
//

int MiscMath::position2leftepoch( uint64_t p , uint64_t e_length , uint64_t e_overlap , int e_total )
{  
  // find the 'first' epoch that overlaps this position
  // assumes 0-based positions and epochs
  uint64_t a = p / e_overlap; // right-most epoch
  const uint64_t b = p % e_overlap; // remainder 
  //std::cout << "b el = " << b << "\n";
  if ( b >= e_length ) return -1; // in-between epochs (i.e. where overlap > length)
  const double actual_overlap = e_length - e_overlap;
  if ( actual_overlap < 0 ) return a; // non-overlapping epochs
  const int adjustment = (e_length-b-1)/e_overlap;
  if ( adjustment > a ) a = 0;  
  else a -= adjustment;
  if ( e_total > 0 && a >= e_total ) return -1; // problem -- gone past end
  return a;
}

int MiscMath::position2rightepoch( uint64_t p , uint64_t e_length , uint64_t e_overlap , int e_total )
{
  // find the 'last' epoch that overlaps this position
  // assumes 0-based positions and epochs
  uint64_t a = p / e_overlap; // right-most epoch
  uint64_t b = p % e_overlap; // remainder
  if ( b >= e_length ) return -1; // in-between epochs (i.e. where overlap > length)
  if ( e_total > 0 && a >= e_total ) return -1; // problem -- gone past end
  return a;
}


//
// Median filter
//

std::vector<double> MiscMath::median_filter( const std::vector<double> & x , int n )
{

  bool odd = n % 2 ; 
  
  // For N odd, Y(k) is the median of X( k-(N-1)/2 : k+(N-1)/2 ).
  // For N even, Y(k) is the median of X( k-N/2 : k+N/2-1 ).

  const int t = x.size();

  std::vector<double> ret( t , 0 );
  
  int v1 = odd ? (n-1)/2 : n/2;
  int v2 = odd ? (n-1)/2 : n/2-1;
  
  for (int i = 0 ; i < t ; i++ ) 
    {
      std::vector<double> y( n , 0 );
      int cnt = 0;
      for ( int j = i - v1 ; j <= i + v2 ; j++ )
	{
	  if ( j >= 0 && j < t ) y[cnt] = x[j];
	  ++cnt;
	}

      // get median
      ret[i] = median_destroy( &y[0] , y.size() );
      
    }
  
  return ret;
  
}



//
// Moving average
//

std::vector<double> MiscMath::moving_average( const std::vector<double> & x , int s )
{
  
  if ( s == 1 ) return x;

  const int n     = x.size();

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

  std::vector<double> a( n , 1.0/(double)s );

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



std::vector<double> MiscMath::moving_average_filter( const std::vector<double> & x , int s )
{
  
  // use same approach as filter( 1/n * (ones(1,n) , 1 , x ) 
  // i.e. sum of the last 's' values (including current)

  if ( s == 1 ) return x;
  
  double c = 1.0/(double)s;
   
  const int n     = x.size();
  if ( s >= n ) Helper::halt( "need s < n for moving average" );
  
  
  std::vector<double> r(n,0);
  
  for (int i=0;i<n;i++)
    {
      int j = i - s + 1;
      double z = 0;
      if ( j < 0 ) j = 0;
      for (int k=j;k<=i;k++) z += x[k];
      r[i] = z * c;
    }
  return r;  
}


std::vector<double> MiscMath::detrend( const std::vector<double> & x , double * pa , double * pb )
{
  std::vector<double> r = x;
  detrend(&r,pa,pb);
  return r;
}

void MiscMath::detrend( std::vector<double> * y , double * pa , double * pb )
{
  const int n = y->size();
  // 'x' is 0:(n-1)  
  double yfirst = (*y)[0];
  double ylast  = (*y)[ y->size() - 1 ];
  double b = ( yfirst - ylast ) / (double)( 0 - ( n - 1 ) );
  double a = yfirst ; // x[0] = 0 is intercept
  for (int i=0;i<n;i++) (*y)[i] = (*y)[i] - ( a + b * i );   
  *pa = a;
  *pb = b;
}


std::vector<double> MiscMath::centre( const std::vector<double> & x )
{
  std::vector<double> r = x;
  centre(&r);
  return r;
}

void MiscMath::centre( std::vector<double> * x )
{
  const int n = x->size();
  double s = 0;
  for (int i=0;i<n;i++) s += (*x)[i];
  double mean = s / (double)n;
  for (int i=0;i<n;i++) (*x)[i] -= mean;
}


double MiscMath::covariance( const std::vector<double> & x ,
			     const std::vector<double> & y ,
			     const int w )
{

  if ( w < 1 ) return 0;
  if ( x.size() != y.size() ) return 0;
  
  if ( w == 1 ) 
    {      



      
      const int n = x.size();
      if ( n < 2 ) return 0;

      double mx = MiscMath::mean( x );
      double my = MiscMath::mean( y );

      // cross product
      double sxy = 0;      
      for (int i=0;i<n;i++)
	sxy += ( x[i] - mx ) * ( y[i] - my );
      
      return sxy/(double)(n-1);
      
    }
  else
    {
      
      std::vector<double> x2 = MiscMath::moving_average( x , w );
      std::vector<double> y2 = MiscMath::moving_average( y , w );
      
      const int n = x2.size();
      if ( n < 2 ) return 0;

      double mx = MiscMath::mean( x2 );
      double my = MiscMath::mean( y2 );

      // cross product
      double sxy = 0;      
      for (int i=0;i<n;i++)
	sxy += ( x2[i] - mx ) * ( y2[i] - my );
      
      return sxy/(double)(n-1);

    }
  
  return 0;
}



double MiscMath::overdispersion( const std::vector<int> & a , double * pv )
{

  const int n = a.size();
  // for now, just add some upper bound, come back and fix this
  int mx = 0;
  for (int i=0;i<n;i++) if ( a[i] > mx ) mx = a[i];
  if ( mx > 100 ) Helper::halt("bailed in overdispersion test..." );
  
  const double m = mean( a );
  const double s2 = variance( a );

  // generate expected values
  std::vector<double> e(mx+1);
  for (int i=0;i<=mx;i++) e[i] = n * poisson( i , m );
  
  // and tabulate corresponding observed table
  std::vector<double> o(mx+1);
  for (int i=0;i<n;i++) ++o[a[i]];

  double pval = chisq( o , e );

  // std::cout << "stat  = " << s2/m << " " << pval << "\n";
  // double stat2 = sqrt( (n-1)/2.0 ) * ( (s2/m) - 1 );
  // double pval2 = Statistics::chi2_prob( stat2 , 1 );
  // std::cout << "stat2 = " << stat2 << " " << pval2 << "\n";
  
  if ( pv != NULL ) *pv = pval;

  // return over-dispersion parameter
  return m == 0 ? 0 : s2 / m;

}

  
double MiscMath::poisson(const double k, const double lambda) 
{
  return exp(k * log(lambda) - lgamma(k + 1.0) - lambda);
}

double MiscMath::chisq( const std::vector<double> & o , const std::vector<double> & e )
{

  const int n = o.size();
  if ( e.size() != n ) Helper::halt( "problem in chisq()" );

  double x = 0;
  int k = 0;

  for (int i=0;i<n;i++) 
    {
      if ( e[i] > 1.0 )
	{
	  x += ( ( o[i] - e[i] ) * ( o[i] - e[i] ) ) / e[i] ;
	  ++k;
	}
    }
  if ( k < 2 ) return 1;
  return Statistics::chi2_prob( x , k-1 );
}


double MiscMath::empirical_pvalue( const double s , const std::vector<double> & x )
{
  double pv = 1;
  for (int i=0;i<x.size();i++) if ( x[i] >= s ) ++pv;
  return pv / (double)( x.size() + 1.0 );
}


void MiscMath::normalize( std::vector<double> * x , double * mn , double * mx)
{
  // get min/max
  minmax( *x , mn , mx );
  double denom = *mx - *mn;
  const int n = x->size();
  for (int i=0;i<n;i++) (*x)[i] = ( (*x)[i] - *mn ) / denom ;   
}

// min/max normalization with mask

void MiscMath::normalize( std::vector<double> * x , const std::vector<bool> & include_mask )
{
  std::vector<double> nx;
  std::vector<int> ox;
  if ( x->size() != include_mask.size() ) Helper::halt( "error in normalize()" );
  
  
  for (int i=0;i<x->size();i++) 
    {
      if ( include_mask[i] ) 
	{
	  nx.push_back( (*x)[i] );
	  ox.push_back( i );
	}
    }

  const int n = nx.size();

  if ( n == 0 ) return;
 
  double mn, mx;

  minmax( nx , &mn , &mx );

  double denom = mx - mn;

  // alter original 
  for (int i=0;i<n;i++)
    {
      (*x)[ ox[i] ] = ( nx[i] - mn ) / denom;
    }
  
}



void MiscMath::minmax( const std::vector<double> & x , double * mn , double * mx)
{

  const int n = x.size();

  if ( n == 0 ) 
    {
      *mn = *mx = 0;
      return;
    }
  
  *mn = *mx = x[0];
  for (int i=1;i<n;i++)
    {
      if      ( x[i] < *mn ) *mn = x[i];
      else if ( x[i] > *mx ) *mx = x[i]; 
    }
  
}


double MiscMath::threshold( const std::vector<double> & x , double lwr, double upr, double inc , 
			    double * w , 
			    std::map<double,double> * tvals )
{

  // Otsu (1979) A Threshold Selection Method from Gray-Level
  // Histograms

  if ( tvals != NULL ) 
    tvals->clear();  
  
  std::map<double,int> l;

  const int n = x.size();

  double grand_mean = 0;

  for (int i=0;i<n;i++) 
    {
      l[x[i]]++;
      grand_mean += x[i];
    }
  
  grand_mean /= (double)n;

  std::cout << "grand mean = " << grand_mean << "\n";
  
  double cum_sum = 0;
  double cum_f = 0;
  double max_sigma_b = 0;
  double max_t = 0;  // here t implies up to and including 't'
  double max_t2 = 0; // binned version
  double best_f = 0;

  double t = lwr;

  // previous 't' (i.e. lowest possible value)
  double last_t = l.begin()->first;
  
  std::cout << "starting t = " << t << "\n"
	    << "last (previous) t = " << last_t << "\n";


  // iterate over all observed values
  
  std::map<double,int>::const_iterator ii = l.begin();
  while ( ii != l.end() )
    {
      
      const double this_t = ii->first;
      
      std::cerr << "observed value " << this_t << "\n";
      
      // check we don't skip a category
      if ( this_t > t + inc )
	{
	  std::cerr << "updating t... from " << t << "\n";
	  while ( 1 ) 
	    {
	      t += inc;
	      if ( this_t <= t ) break;
	      
	    }
	  std::cerr << "t is now " << t << "\n";
	}

      
      std::cout << "test threshold t = " << t << "\n";
      
      cum_f += ii->second;
      cum_sum += this_t * ii->second;

      std::cout << "updating cumulative sum = " << cum_sum << "\n";	    
           
      // a test-point?
      // i.e. if we've gone one past (or equal to) the current threshold
      
      if ( this_t >= t && last_t < t )
	{

	  std::cout << "  -- triggering evaluation\n";
	  
	  const double f = cum_f / (double)n;
	  const double m = cum_sum / cum_f;
	  
 	  std::cerr << "w = " << f << "\n";
 	  // std::cerr << "m = " << m << "\n";
	  
	  if ( f >= 0 || f <= 1 ) 
	    {
	      
	      const double sigma_b = ( ( grand_mean * f - m ) * ( grand_mean * f - m ) ) / ( f * (1-f) );
	      
	      if ( sigma_b > max_sigma_b ) 
		{
		  max_sigma_b = sigma_b;
		  max_t = this_t;
		  max_t2 = t;
		  best_f = f;
		}
	      
	      if ( tvals != NULL ) (*tvals)[ t ] = sigma_b;

	      std::cout << " sigma_B\t" << sigma_b << "\n";
	      
 	      // std::cout << "details " << t << "\t"
	      // 		<< f << "\t"
	      // 		<< ii->first << "\t"
 	      // 		<< sigma_b << "\t"
 	      // 		<< max_sigma_b << "\t"
 	      // 		<< max_t << "\t" 
 	      // 		<< max_t2 << "\n";
	      
	    }
	  
	  // go to the next t
	  t += inc;	  	    
	  
	  if ( t > upr ) break;
	}
      
      last_t = this_t;
      
      ++ii;
    }

  // normalize tvals
  std::map<double,double>::iterator tt = tvals->begin();
  while ( tt != tvals->end() )
    {
      tt->second /= max_sigma_b;
      ++tt;
    }
  
  // i.e. threshold is x > t 
  //      rather than x >= t

    std::cerr << "maximum threshold is " << max_t << " " << max_t2 << "\n";
  if ( w != NULL ) *w = 1 - best_f;
  return max_t2;
   
}


double MiscMath::rad2deg(double radians) 
{ 
  return radians * (180.0 / M_PI); 
}

double MiscMath::deg2rad(double degrees) 
{ 
  return degrees * (M_PI / 180.0); 
}

double MiscMath::shift_degrees( double d , double x ) 
{
  d += x;
  while ( d >= 360 ) { d -= 360; } 
  while ( d < 0  ) { d += 360; } 
  return d;
}


double MiscMath::as_angle_0_pos2neg( const double r )
{
  // also as angle, 0..360 such that 0 if pos-to-neg cross-over 
  // 90 DN, 180 neg-to-pos, 270 UP
  
  // original                      converted
  
  // -pi      DN       -180        pi/2    90     
  // -pi/2    neg2pos   -90        pi      180   
  // 0        UP          0        3/2pi   270
  // +pi/2    pos2neg    90        0 / 2pi 0/360
  // pi       DN        180        pi/2    90
  
  double a = MiscMath::rad2deg( r ) + 270;

  while ( a >= 360 ) a -= 360;

  return a;
}


double MiscMath::accuracy( const std::vector<int> & a , const std::vector<int> & b , std::vector<double> * spacc )
{
  std::vector<std::string> aa( a.size() );
  std::vector<std::string> bb( b.size() );
  for (int i=0;i<a.size();i++) aa[i] = Helper::int2str( a[i] );
  for (int i=0;i<b.size();i++) bb[i] = Helper::int2str( b[i] );
  return accuracy( aa, bb , spacc );
}

double MiscMath::accuracy( const std::vector<std::string> & a , const std::vector<std::string> & b , std::vector<double> * spacc)
{
  const int n = a.size();
  if ( n != b.size() ) Helper::halt( "mismatched vectors in accuracy()" );

  int m = 0;
  //std::map<std::string,int> match;

  for (int i=0;i<n;i++)
    {
      if ( a[i] == b[i] )
	{
	  ++m;
	  //match[ a[i] ] = 1;
	}
    }
  
  // items specific accuracies
  if ( spacc != NULL )
    {
      // leave for now
    }

  // overall accuracy
  return m / (double)n;
  
}


double MiscMath::kappa( const std::vector<int> & a , const std::vector<int> & b )
{
  std::vector<std::string> aa( a.size() );
  std::vector<std::string> bb( b.size() );
  for (int i=0;i<a.size();i++) aa[i] = Helper::int2str( a[i] );
  for (int i=0;i<b.size();i++) bb[i] = Helper::int2str( b[i] );
  return kappa( aa, bb );
}

double MiscMath::kappa( const std::vector<std::string> & a , const std::vector<std::string> & b )
{
  if ( a.size() != b.size() )
    Helper::halt( "unequal input vectors for kappa()" );

  std::map<std::string,int> allcounts;
  std::map<std::string,double> acounts;
  std::map<std::string,double> bcounts;
  std::map<std::string,std::map<std::string,double > > abcounts;

  const int n = a.size();

  if ( n == 0 ) return 0;
  double inc = 1/(double)n;
  
  for ( int i = 0 ; i < a.size() ; i++ )
    {
      //std::cout << "ZZ\t" << a[i] << "\t" << b[i] << "\n";
      
      allcounts[ a[i] ]++ ; allcounts[ b[i] ]++;      
      acounts[ a[i] ] += inc;
      bcounts[ b[i] ] += inc;
      abcounts[ a[i] ][ b[i] ] += inc;
    }
  
  
  double observed_agreement = 0;
  double chance_agreement = 0;
  //  std::cout << "allcounts = " << allcounts.size() << "\n";
  
  std::map<std::string,int>::const_iterator cc = allcounts.begin();
  while ( cc != allcounts.end() )
    {
      // std::cout << " cc-> " << cc->first << "\t"
      // 		<< abcounts[ cc->first ][ cc->first ] << "\t"
      // 		<< acounts[ cc->first ] << " * " <<  bcounts[ cc->first ] << " = " << acounts[ cc->first ] * bcounts[ cc->first ] << "\n";
      observed_agreement += abcounts[ cc->first ][ cc->first ] ;
      chance_agreement += acounts[ cc->first ] * bcounts[ cc->first ];
      ++cc;
    }

  //  std::cout << " observed_agreement = " << observed_agreement << " " << chance_agreement << "\n";
  
  double kappa = ( observed_agreement - chance_agreement ) / ( 1 - chance_agreement );
  return kappa;
  
}


int MiscMath::nearest_idx( const std::vector<double> & x , double value , int lwr  , int upr )
{
  if ( x.size() == 0 ) return -1;
  int start = lwr >= 0 ? lwr : 0 ;
  int stop  = upr >= 0 ? upr : x.size() - 1 ;
  int nidx = -1;
  double diff = 0;

  for (int i=start;i<=stop;i++)
    {
      if ( nidx == -1 )
	{
	  diff = fabs( x[i] - value );
	  nidx = i;
	}
      else
	{
	  double d = fabs( x[i] - value );
	  if ( d < diff )
	    {
	      nidx = i;
	      diff = d;
	    }
	}	
    }
  return nidx;
}



double MiscMath::pF(const double F, const int df1, const int df2)
{
  return betai(0.5*df2,0.5*df1,(double)df2/(double)(df2+df1*F));
}


double MiscMath::betai(const double a, const double b, const double x)
{
  double bt;
  
  if (x < 0.0 || x > 1.0) Helper::halt("Internal error: bad x in routine betai");
  if (x == 0.0 || x == 1.0) bt=0.0;
  else
    bt=exp(Statistics::gammln(a+b)-Statistics::gammln(a)-Statistics::gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a;
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
}


double MiscMath::betacf(const double a, const double b, const double x)
{
  
  const int MAXIT = 100;
  const double EPS = 3e-7;
  const double FPMIN = 1.0e-30;

  int m,m2;
  double aa,c,d,del,h,qab,qam,qap;
  
  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) <= EPS) break;
  }
  if (m > MAXIT) Helper::halt("Internal error in betacf() function (please report)");
  return h;
}



// second-order Hjorth parameters (window, inc)
void MiscMath::hjorth2( const std::vector<double> * x , double * r , int w , int inc )
{

  // expect r[9] to be a 9 element array; (3x3)
  //  0,1,2   H1H1 H1H2 H1H3
  //  3,4,5   H2H1 H2H2 H2H3
  //  6,7,8   H3H1 H3H2 H3H3

  // w   window size (in sample points)
  // inc window increment (in sample points)

  if ( inc == 0 ) inc = w;

  const int nx = x->size();
  const int nw = nx / w;

  std::vector<double> h1, h2, h3;
 
  for (int i=0; i<nx; i += inc )
    {
      std::vector<double> t;
      for (int j=i;j<i+w;j++)
	t.push_back( (*x)[j] );

      // normalize window
      t = MiscMath::Z( t );
      
      double hp1, hp2, hp3;
      hjorth( &t , &hp1, &hp2, &hp3 );
      h1.push_back( hp1 );
      h2.push_back( hp2 );
      h3.push_back( hp3 );      

    }

  if ( h1.size() != nw ) Helper::halt( "internal error in hjorth2()" );
  
  hjorth( &h1 , r   , r+1 , r+2 );
  hjorth( &h2 , r+3 , r+4 , r+5 );
  hjorth( &h3 , r+6 , r+7 , r+8 );
  
}


int MiscMath::outliers( const std::vector<double> * x , double th ,
			std::vector<bool> * inc , 
			const std::vector<bool> * prior )
{  
  // only set 'inc' to missing; but if 'prior' exists, set any missing in 'prior'
  // to missing in 'inc'
  // i.e. never unset (set to include) any in inc.

  int removed = 0;

  if ( prior == NULL )
    {
      std::vector<double> z = MiscMath::Z(*x);
      const int n = z.size();
      for (int i=0;i<n; i++)
	if ( z[i] < -th || z[i] > th )
	  {
	    ++removed;
	    (*inc)[ i ] = false;
	  }
    }
  else
    {
      std::vector<double> xx;
      std::vector<double> xi;

      for (int i=0;i<x->size();i++)
	if ( (*prior)[i] )
	  {
	    xx.push_back( (*x)[i] );
	    xi.push_back( i );
	  }
	else
	  {
	    (*inc)[i] = false;
	  }

      std::vector<double> z = MiscMath::Z(xx);
      const int n = z.size();
      for (int i=0;i<n; i++)
	if ( z[i] < -th || z[i] > th )
	  {
	    (*inc)[ xi[i] ] = false;
	    ++removed;
	  }
    }
  
  return removed;

}

