
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
#include "miscmath/dynam.h"

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

std::vector<double> MiscMath::log2space(double a, double b, int n)
{
  if ( n < 2 ) Helper::halt( "log2space requires at least two values" );
  // for a->b, generate a series of logarithmically (log2) equally spaced intervals
  a = log2(a);
  b = log2(b);
  const double st = (b-a)/(double)(n-1);
  std::vector<double> r(n);
  r[0] = pow(2,a);
  r[n-1] = pow(2,b);
  for (int i=1;i<n-1;i++) r[i] = pow( 2, a + i*st );
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


double MiscMath::petrosian_FD( const std::vector<double> & x )
{
  const int n = x.size();
  // not defined
  if ( n < 3 ) return 0;
  std::vector<bool> b(n-1);  
  for (int i=1;i<n; i++)
    b[i-1] = x[i] - x[i-1] > 0 ;  
  int n_delta = 0;
  for (int i=1;i<n-1; i++)
    if ( b[i] != b[i-1] ) ++n_delta;    
  double fd = log10(n) / ( log10(n) + log10( n / (double)( n + 0.4 * n_delta ) ) ) ;
  return fd;
}

double MiscMath::kurtosis( const std::vector<double> & x )
{
  std::vector<double> d = x;
  const double m = MiscMath::mean( d ) ;
  for (int i=0; i<d.size(); i++) d[i] -= m;
  return kurtosis0( d );
}

double MiscMath::kurtosis0( const std::vector<double> & x )
{  
  // kurtosis ( assume mean = 0 )
  // assumes mean = 0

  const int n = x.size();
  
  double numer = 0 , denom = 0;
  for (int i=0; i<n; i++)
    {
      numer += pow( x[i] , 4 );
      denom += pow( x[i] , 2 );
    }
  numer /= (double)n;
  denom /= (double)n;
  denom *= denom;
  return numer / denom - 3.0;
  
}

double MiscMath::kurtosis( const std::vector<double> & x , double m )
{
  std::vector<double> d = x;
  for (int i=0; i<d.size(); i++) d[i] -= m;
  return kurtosis0( d );
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
  

double MiscMath::median( const std::vector<double> & x , const bool also_upper )
{
  
  const int n = x.size();  

  const bool is_odd = n % 2;
  
  if ( n == 0 ) Helper::halt( "internal problem, taking median of 0 elements");
  if ( n == 1 ) return x[0];
  
  if ( is_odd ) 
    return MiscMath::kth_smallest_preserve( x , ( n - 1 ) / 2 );
  
  const double lower_median = MiscMath::kth_smallest_preserve( x , n / 2 - 1 );
  
  if ( ! also_upper ) return lower_median;
  
  const double upper_median = MiscMath::kth_smallest_preserve( x , n / 2 );
  
  return ( lower_median + upper_median ) / 2.0 ;
      
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

std::vector<double> MiscMath::remove_median_filter( const std::vector<double> & x , int n , std::vector<double> * p )
{
  std::vector<double> f = MiscMath::median_filter( x , n );
  if ( p != NULL ) *p = f; 
  for (int i=0; i<f.size(); i++)
    f[i] = x[i] - f[i];
  return f;
}
  
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
	if ( j >= 0 && j < t ) y[cnt++] = x[j];
      
      // get median
      ret[i] = median_destroy( &y[0] , cnt );
      
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
  
  // assume equal spacing
  std::vector<double> x( n );
  for (int i=0; i<n; i++) x[i] = i;
    
  // fit line  
  dynam_t spec_slope( *y , x );
  double beta, m;
  spec_slope.linear_trend( &beta , NULL , &m);
  
  // adjust 
  for (int i=0; i<n; i++) (*y)[i] -= m + beta * x[i];
  
  if ( pa ) *pa = m;
  if ( pb ) *pb = beta;

  // OLD

  // // 'x' is 0:(n-1)  
  // double yfirst = (*y)[0];
  // double ylast  = (*y)[ y->size() - 1 ];
  // double b = ( yfirst - ylast ) / (double)( 0 - ( n - 1 ) );
  // double a = yfirst ; // x[0] = 0 is intercept
  // // adjust y
  // for (int i=0;i<n;i++) (*y)[i] = (*y)[i] - ( a + b * i );   
  // // return estimates
  // if ( pa ) *pa = a;
  // if ( pb ) *pb = b;

}


std::vector<double> MiscMath::edge_detrend( const std::vector<double> & x , double * pa , double * pb )
{
  std::vector<double> r = x;
  edge_detrend(&r,pa,pb);
  return r;
}

void MiscMath::edge_detrend( std::vector<double> * y , double * pa , double * pb )
{
  const int n = y->size(); 
  // 'x' is 0:(n-1)  
  double yfirst = (*y)[0];
  double ylast  = (*y)[ y->size() - 1 ];
  double b = ( yfirst - ylast ) / (double)( 0 - ( n - 1 ) );
  double a = yfirst ; // x[0] = 0 is intercept
  // adjust y
  for (int i=0;i<n;i++) (*y)[i] = (*y)[i] - ( a + b * i );   
  // return estimates
  if ( pa ) *pa = a;
  if ( pb ) *pb = b;  
}


std::vector<double> MiscMath::centre( const std::vector<double> & x )
{
  std::vector<double> r = x;
  centre(&r);
  return r;
}

double MiscMath::centre( std::vector<double> * x )
{
  const int n = x->size();
  double s = 0;
  for (int i=0;i<n;i++) s += (*x)[i];
  double mean = s / (double)n;
  for (int i=0;i<n;i++) (*x)[i] -= mean;
  return mean;
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


double MiscMath::max(const std::vector<double> & x )
{
  double mn, mx;
  minmax( x , &mn , &mx );
  return mx;
}

double MiscMath::min(const std::vector<double> & x )
{
  double mn, mx;
  minmax( x , &mn , &mx );
  return mn;
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
      
      //std::cerr << "observed value " << this_t << "\n";
      
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
      
      
      //std::cout << "test threshold t = " << t << "\n";
      
      cum_f += ii->second;
      cum_sum += this_t * ii->second;

      //std::cout << "updating cumulative sum = " << cum_sum << "\n";	    
      
      // a test-point?
      // i.e. if we've gone one past (or equal to) the current threshold
      
      if ( this_t >= t && last_t < t )
	{
	  
	  std::cout << "  -- triggering evaluation\n";
	  
	  const double f = cum_f / (double)n;
	  const double m = cum_sum / cum_f;
	  
 	  std::cerr << "w = " << f << "\n";
	  std::cerr << "m = " << m << "\n";
	  
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
	      
	      std::cout << "details " << t << "\t"
	       		<< f << "\t"
	       		<< ii->first << "\t"
 	       		<< sigma_b << "\t"
 	       		<< max_sigma_b << "\t"
 	       		<< max_t << "\t" 
 	       		<< max_t2 << "\n";
	      
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



double MiscMath::threshold2( const std::vector<double> & x ,
			     double * empf , 
			     const int k , 
			     std::map<double,double> * fvals , 
			     std::map<double,double> * tvals )
{
  
  // Otsu (1979) A Threshold Selection Method from Gray-Level
  // Histograms

  if ( empf == NULL )
    Helper::halt( "internal error calling threshold2() " );
  
  if ( fvals != NULL ) fvals->clear();
  if ( tvals != NULL ) tvals->clear();  
    
  if ( k > 0 && ( fvals == NULL || tvals == NULL ) )
    Helper::halt( "internal error calling threshold2() " );
  
  const int n = x.size();
  
  double sum = 0;

  std::map<double,int> l;
  for (int i=0;i<n;i++) 
    {
      l[ x[i] ]++;
      sum += x[i];
    }
  
  double sumB = 0;
  double wB = 0;
  double wF = 0;
  
  double varMax = 0;
  double threshold = 0;
  
  double cnt = 0;
  
  std::map<double,int>::const_iterator ii = l.begin();
  while ( ii != l.end() )
    {

      // track above-threshold proportion
      if ( fvals != NULL )
	{
	  cnt += ii->second ;
	  (*fvals)[ ii->first ] = cnt / (double)n ; 
	}
      
      wB += ii->second;
      if ( wB == 0 ) { ++ii; continue; }

      wF = n - wB;
      if ( wF == 0 ) break;

      sumB += ii->first * ii->second;
      
      double mB = sumB / wB ;
      double mF = ( sum - sumB ) / wF ;
      
      double varBetween = wB * wF * ( mB - mF ) * ( mB - mF );
      
      if ( tvals != NULL )
	(*tvals)[ ii->first ] = varBetween ; 
      
      if ( varBetween > varMax )
	{
	  varMax = varBetween;
	  threshold = ii->first ;
	  //std::cout << " setting " << cnt << " " << n << " " << cnt/(double)n << "\n";
	  *empf = cnt / (double)n ;
	}

      // next possible value 
      ++ii;
      
    }
  
  //
  // splice out a subset 
  //
  
  if ( k > 0 )
    {

      std::map<double,double> f2 = *fvals;
      std::map<double,double> t2 = *tvals;
      
      fvals->clear();
      tvals->clear();

      
      int k2 = t2.size() / k ; 
      
      int p = 0;
      std::map<double,double>::const_iterator ii = t2.begin();
      while ( ii != t2.end() )
	{
	  if ( p % k2 == 0 )
	    {
	      (*tvals)[ ii->first ] = ii->second;
	      (*fvals)[ ii->first ] = f2[ ii->first ];
	    }
	  ++p;
	  ++ii;
	}
      
    }

  //
  // normalize t-vals
  //

  if ( tvals != NULL )
    {
      std::map<double,double>::iterator ii = tvals->begin();
      while ( ii != tvals->end() )
	{
	  (*tvals)[ ii->first ] = ii->second / varMax;
	  ++ii;
	}
    }
  
  return threshold;

}


double MiscMath::rad2deg(double radians) 
{ 
  return radians * (180.0 / M_PI); 
}

double MiscMath::deg2rad(double degrees) 
{ 
  return degrees * (M_PI / 180.0); 
}

double MiscMath::angle_difference( double a , double b )
{
  // get signed distance from a --> b, allowing for wrapping
  // assumes input = 0..360

  if ( a < 0 || a > 360 || b < 0 || b > 360 )
    Helper::halt(" angle_difference expecting 0 - 360 " );

  //   0     45   50    360
  if ( a == b ) return 0;

  // assuming in same region
  double d1 = fabs( b - a );

  // if wraps
  double d2 = a > b ? b + 360 - a : -( a + 360 - b ) ;
  
  // standard metric 
  if ( d1 < fabs( d2 ) ) 
    return b - a ; 
  
  return d2;
  
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


double MiscMath::accuracy( const std::vector<int> & a , const std::vector<int> & b , 
			   const int unknown , 
			   std::vector<int> * labels ,
			   std::vector<double> * precision ,
			   std::vector<double> * recall ,
			   std::vector<double> * f1 ,
			   double * macro_precision ,
			   double * macro_recall ,
			   double * macro_f1 ,
			   double * avg_weighted_precision ,
			   double * avg_weighted_recall ,
			   double * avg_weighted_f1 ,
			   double * mcc_val )
{
  std::vector<std::string> aa( a.size() );
  std::vector<std::string> bb( b.size() );

  for (int i=0;i<a.size();i++) aa[i] = a[i] == unknown ? "?" : Helper::int2str( a[i] );
  for (int i=0;i<b.size();i++) bb[i] = b[i] == unknown ? "?" : Helper::int2str( b[i] );

  std::vector<std::string> ll;
  if ( labels != NULL )
    {
      ll.resize( labels->size() );
      for (int i=0;i<labels->size();i++) 
	{
	  if ( (*labels)[i] == unknown ) 
	    Helper::halt( "internal error in accuracy(): cannot specify unknown value as an explicit label" );
	  ll[i] = Helper::int2str( (*labels)[i] );
	}
    }

  return accuracy( aa, bb , "?" , 
		   labels == NULL ? NULL : &ll ,
		   precision , recall , f1 ,
		   macro_precision , macro_recall , macro_f1 , 
		   avg_weighted_precision , avg_weighted_recall , avg_weighted_f1 , mcc_val );
}

double MiscMath::accuracy( const std::vector<std::string> & a , const std::vector<std::string> & b ,
			   const std::string & unknown , 
			   std::vector<std::string> * labels , 
                           std::vector<double> * precision ,
                           std::vector<double> * recall ,
                           std::vector<double> * f1 ,
			   double * macro_precision ,
			   double * macro_recall ,
			   double * macro_f1 ,
			   double * avg_weighted_precision ,
			   double * avg_weighted_recall ,
			   double * avg_weighted_f1 ,
			   double * mcc_val )
{
  const int n = a.size();

  if ( n != b.size() ) Helper::halt( "mismatched vectors in accuracy()" );

  int m = 0;

  std::map<std::string,std::map<std::string,int> > table;
  std::map<std::string,int> rows, cols;
  std::set<std::string> obs;
  int n_obs = 0;
  for (int i=0;i<n;i++)
    {
      // only consider known cells
      if ( a[i] != unknown && b[i] != unknown )
	{
	  // overall accuracy
	  if ( a[i] == b[i] ) ++m;
	  // build table
	  ++table[ a[i] ][ b[i] ];
	  ++rows[ a[i] ];
	  ++cols[ b[i] ];
	  obs.insert( a[i] );
	  obs.insert( b[i] );      
	  ++n_obs;
	}
    }

  
  // items specific precision/recall/F1
  if ( labels != NULL )
    {
      // assume a [ row ] is the 'truth'
      //        b [ col ] is the prediction

      // so, precision = ii / sum(i*)
      //     recall    = ii / (*i)

      precision->resize( labels->size() );
      recall->resize( labels->size() );
      f1->resize( labels->size() );

      *macro_recall = 0;
      *macro_precision = 0;
      *macro_f1 = 0;

      *avg_weighted_f1 = 0;
      *avg_weighted_precision = 0;
      *avg_weighted_recall = 0;
      
      int ncat = 0;
      int nobs = 0;
            
      // for (int i=0;i<labels->size();i++)
      // 	if ( rows.find( (*labels)[i] ) != rows.end() ) 
      // 	  std::cout << "ROW " << (*labels)[i] << " " << rows.find((*labels)[i] )->second << "\n";

      // for (int i=0;i<labels->size();i++)
      // 	if ( cols.find( (*labels)[i] ) != cols.end() ) 
      // 	  std::cout << "COL " << (*labels)[i] << " " << cols.find((*labels)[i] )->second << "\n";

      for (int i=0;i<labels->size();i++)
	{
	  // only counts 'observed' labels, so no need to handle missing cats '?'
	  // separately here

	  if ( obs.find( (*labels)[i] ) != obs.end() )
	    {

	      ++ncat;
	      nobs += rows[ (*labels)[i] ] ;
	      
	      (*precision)[i] = table[ (*labels)[i] ][ (*labels)[i] ] / (double)cols[ (*labels)[i] ];
	      (*recall)[i] = table[ (*labels)[i] ][ (*labels)[i] ] / (double)rows[ (*labels)[i] ];
	      (*f1)[i] = 2 * ( (*precision)[i] * (*recall)[i] ) / ( (*precision)[i] + (*recall)[i] ); 
	      
	      *macro_f1 += (*f1)[i];
	      *macro_precision += (*precision)[i];
	      *macro_recall += (*recall)[i];
	      
	      *avg_weighted_f1 += rows[ (*labels)[i] ] * (*f1)[i];
	      *avg_weighted_precision += rows[ (*labels)[i] ] * (*precision)[i];
	      *avg_weighted_recall += rows[ (*labels)[i] ] * (*recall)[i];
	    }
	  else
	    {
	      (*precision)[i] = 0;
	      (*recall)[i] = 0;
	      (*f1)[i] = 0;
	    }

	}

      *macro_f1 /= (double)ncat;
      *macro_recall /= (double)ncat;
      *macro_precision /= (double)ncat;

      *avg_weighted_f1 /= (double)nobs;
      *avg_weighted_precision /= (double)nobs;
      *avg_weighted_recall /= (double)nobs;
      
    }

  // MCC?
  if ( mcc_val != NULL )
    {
      std::vector<std::string> obsvec;
      std::set<std::string>::const_iterator ii = obs.begin();
      while ( ii != obs.end() )  { obsvec.push_back( *ii ) ; ++ii ; } 
      *mcc_val = mcc( table , obsvec );
    }

  // overall accuracy (considering only non-missing obs)
  if ( n_obs == 0 ) return 0;
  return m / (double)n_obs;
  
}


double MiscMath::mcc( std::map<std::string,std::map<std::string,int> > table ,
		      const std::vector<std::string> & labels )
{
  // Implements MCC for multiclass data, based on the gawk implementation below:
  // https://rth.dk/resources/rk/software/rkorrC
  // # Computes the K-category correlation coefficient.
  // # Copyright: Jan Gorodkin, gorodkin@bioinf.kvl.dk
  // # This software is distributed with a 
  // # GNU GENERAL PUBLIC LICENSE, see http://www.gnu.org/licenses/gpl.txt
  // # Copyright (C) 2004  Jan Gorodkin
  // # For publication of results, please cite:
  // #  Comparing two K-category assignments by a K-category correlation coefficient.
  // #  J. Gorodkin, Computational Biology and Chemistry, 28:367-374, 2004.

  const int nk = labels.size();
  int N = 0;
  
  Data::Matrix<double> C( nk, nk );
  for (int r=0;r<nk;r++)
    for (int c=0;c<nk;c++)
      {
	C(r,c) = table[ labels[r] ][ labels[c] ] ;
	N += C(r,c);
      }

  //  std::cout << C.print() << "\n";
  
  // trace

  double tr = 0;
  for (int i=0;i<nk;i++) tr += C(i,i);

  // sum row col dot product:
  
  double rowcol_sumprod = 0;
  for(int r=0; r<nk; r++)
    for(int c=0; c<nk; c++)
      for (int i=0; i<nk; i++)
	rowcol_sumprod += C(r,i) * C(i,c);
  
  // sum over row dot products
  double rowrow_sumprod = 0;
  for(int r1=0; r1<nk; r1++)
    for(int r2=0; r2<nk; r2++)
      for (int i=0; i<nk; i++)
	rowrow_sumprod += C(r1,i) * C(r2,i);

  // sum over col dot products
  double colcol_sumprod = 0;
  for(int c1=0; c1<nk; c1++)
    for(int c2=0; c2<nk; c2++)
      for (int i=0; i<nk; i++)
	colcol_sumprod += C(i,c1) * C(i,c2);

  double Qk = tr / (double)N;

  double COV_XY = N*tr - rowcol_sumprod; //  NOTE, COV_XY, COV_XX and COV_YY
  double COV_XX = N*N - rowrow_sumprod;  //# are HERE ALWAYS INTEGERS!
  double COV_YY = N*N - colcol_sumprod; 
  double denominator = sqrt( COV_XX*COV_YY );

  double RK = 0;

  if ( denominator  >0 ) {
    RK= COV_XY / denominator;
    //double b = COV_XY / COV_XX;
    //double bprime=COV_XY/COV_YY;
    // output 
    // printf "%13.10f  %3d  %13.10f  %13.10f  %13d  %13d  %13d  %12.10f ", RK,  K, b, bprime, COV_XY, COV_XX, COV_YY, QK;
  }
  else if (denominator == 0 )
    {      
      RK=1;
        // for(i=1;i<startpos;i++) printf misclabel[i]"  ";
        // printf "%13s  %3d  %13s  %13s  %13d  %13d  %13d  %12.10f ", "Nan         ",  K, "Nan         ", "Nan         ", COV_XY, COV_XX, COV_YY, QK;
        // for(l=1;l<K+1;l++) for(k=1;k<K+1;k++) printf "  "C[k,l];
        // printf "\n";
     }

  return RK;

}



double MiscMath::kappa( const std::vector<int> & a , const std::vector<int> & b , const int unknown )
{
  std::vector<std::string> aa( a.size() );
  std::vector<std::string> bb( b.size() );
  for (int i=0;i<a.size();i++) aa[i] = a[i] == unknown ? "?" : Helper::int2str( a[i] );
  for (int i=0;i<b.size();i++) bb[i] = b[i] == unknown ? "?" : Helper::int2str( b[i] );
  return kappa( aa, bb , "?" );
}

double MiscMath::kappa( const std::vector<std::string> & a , const std::vector<std::string> & b , const std::string & unknown )
{
  if ( a.size() != b.size() )
    Helper::halt( "unequal input vectors for kappa()" );

  std::map<std::string,int> allcounts;
  std::map<std::string,double> acounts;
  std::map<std::string,double> bcounts;
  std::map<std::string,std::map<std::string,double > > abcounts;

  // only consider cells where both are 'known' 
  std::vector<bool> incl( a.size() );
  int n = 0;
  for (int i=0;i<a.size();i++)
    {
      incl[i] = a[i] != unknown && b[i] != unknown ;
      if ( incl[i] ) ++n;
    }
  
  if ( n == 0 ) return 0;
  double inc = 1/(double)n;
  
  for ( int i = 0 ; i < a.size() ; i++ )
    {
      //std::cout << "ZZ\t" << a[i] << "\t" << b[i] << "\n";
      if ( incl[i] )
	{	   
	  allcounts[ a[i] ]++ ; allcounts[ b[i] ]++;      
	  acounts[ a[i] ] += inc;
	  bcounts[ b[i] ] += inc;
	  abcounts[ a[i] ][ b[i] ] += inc;
	}
    }
  
  double observed_agreement = 0;
  double chance_agreement = 0;
  
  std::map<std::string,int>::const_iterator cc = allcounts.begin();
  while ( cc != allcounts.end() )
    {
      observed_agreement += abcounts[ cc->first ][ cc->first ] ;
      chance_agreement += acounts[ cc->first ] * bcounts[ cc->first ];
      ++cc;
    }

  double kappa = ( observed_agreement - chance_agreement ) / ( 1.0 - chance_agreement );
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

std::vector<double> MiscMath::outliers( const std::vector<double> * x , double th )
{

  if ( th <= 0 ) return *x;

  const int n = x->size();

  std::vector<bool> inc( n , true );
  
  int removed = MiscMath::outliers( x , th , &inc );

  if ( removed == 0 ) return *x;

  std::vector<double> y;
  for (int i=0; i<n; i++) 
    if ( inc[i] ) y.push_back( (*x)[i] );  
  return y;

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


std::vector<int> MiscMath::smoothedZ( const std::vector<double> & x , 
				      int lag , double threshold , double influence ,
				      int mindur , double max  ,
				      double threshold2 , int mindur2 , 
				      bool noneg , 
				      std::vector<interval_t> * regions , 
				      bool verbose )
{
  const int n = x.size();

  // populate lag window with lag .. 2lag-1  for evaluation of the first 0 .. lag - 1 elements
  // i.e. start going backwards in time, otherwise we have to skip the first region
  // nb. this will not ignore influence of signals/peaks in the second window, but 
  // neither does the original version...
  
  // get global robust SD (to fill in , in case the window has no / little variation ) 

  double global_median = MiscMath::median( x );
  double global_iqr = MiscMath::iqr( x );
  double global_robust_sd = 0.7413 * global_iqr;
  const double sd_eps = global_robust_sd * 1e-3;

  std::vector<int> s( n , 0 );
  
  if ( n <= 2 * lag + 1 ) return s;
  
  std::vector<double> y = x;
  double sum = 0 , sumsq = 0; 

  // track values? in abs SD units: only used if
  bool rec_values = max > 0 || threshold2 > 0 || mindur2 > 0 ;
  std::vector<double> scaled( rec_values ? n : 0 );

  // use "second" window e.g. 30 - 60 seconds to burn-in window
  
  for (int i=lag; i<2*lag; i++)
    {
      sum  += x[i];
      sumsq += x[i] * x[i];
    }
  
  double avg = sum / (double)lag;
  double sd = sqrt( ( lag * sumsq - sum * sum ) / ( (double)((lag-1)*lag) ) );
  if ( sd < sd_eps ) sd = global_robust_sd;

  // starts at 0 now
  for (int i=0; i<n; i++)
    {

      // value in |SD units|
      const double value = std::abs( x[i] - avg ) / sd; 
      
      if ( value > 10000 ) std::cerr << "  warning: large " << i << "  " << x[i] << " " << avg << " " << sd << " " << sd_eps << " " << global_robust_sd << "\n";

      // save scaled value? (for any second round of max/expanded thresholding)
      if ( rec_values ) scaled[i] = value;
      
      // track direction, ercord in s[] and update filtered measure [y] 
      if ( value > threshold ) 
	{
	  s[i] = x[i] > avg ? +1 : -1; 
	  y[i] = influence * x[i] + (1 - influence) * y[i-1] ;
	}
      
      // update sums (respecting that first window 0..lag-1 will be using lag..2lag-1 rather than i - lag 
      int rem = i < lag ? 2 * lag - i - 1 : i - lag; 
      sum = sum - y[ rem ] + y[i] ; 
      sumsq = sumsq - ( y[ rem ] * y[ rem ] ) + ( y[i] * y[i] ) ; 

      // get mean/SD
      // nb. formulate for sample S = (  n * SUMSQ - SUM^2 ) / n(n-1) 
      avg = sum / (double)lag; 
      sd = sqrt( ( lag * sumsq - sum * sum ) / ( (double)((lag-1)*lag) ) );  

      // if the window has no variation, copy over the previous
      if ( sd < sd_eps ) sd = global_robust_sd;

      if ( verbose ) 
	   std::cout << x[i] << "\t" << value << "\t" << s[i] << "\t" << avg - threshold * sd << "\t" << avg + threshold * sd << "\n";
      
    }

  // 
  // Ignore negative peaks?
  //

  if ( noneg ) 
    {
      for (int i=0; i<n; i++)
	if ( s[i] == -1 ) s[i] = 0; 
    }


  //
  // Post-processing to get a) min duration (in sample points) and b) reject any intervals above max threshold
  //

  if ( mindur > 0 ) 
    {
      bool in = s[0] != 0;
      int start = 0;
      for (int i=1; i<n; i++)
	{
	  
	  if ( s[i] != 0 && ! in )
	    {
	      // start new? 
	      start = i;
	      in = true;
	    }
	  else if ( in && ( s[i] == 0 || i == n-1 ) )
	    {
	      
	      // 1-past end
	      int end = s[i] == 0 ? i : n ; 
	      
	      // ending a region
	      if ( end - start < mindur ) // 1 past end so no +1
		for ( int j=start; j<end; j++) s[j] = 0;

	      in = false;
	    }
	}
      
      // did we have a new interval start at the last point? 
      // cannot have a single point, i.e. if mindur is > 
      if ( in && mindur > 1 ) 
	s[ n-1 ] = 0;

    }


  //
  // Above-max events? Consider SD outlier threshold based on MAX of detected regions 
  // only.  i.e. given we have 100 intervals detected, are any > X SD units above the 
  // mean?   Base MAX on the original values, not the scaled ones
  //

  if ( max > 0 ) 
    {
      
      std::vector<double> mxs;
      std::vector<int> starts, stops;

      // first round to get the means

      bool in = s[0] != 0;
      int start= 0;
      double mx = x[0];
      for (int i=1; i<n; i++)
	{

          if ( s[i] != 0 && ! in )
            {
	      // start new?                                                                                                                                                                     
              start = i;
              in = true;
	      mx = x[i];
            }
          else if ( in && ( s[i] == 0 || i == n-1 ) )
            {
	      // include this point in max? (if last)
	      if ( s[i] != 0 ) 
		if ( x[i] > mx ) mx = x[i];
	      
	      // 1 past end
	      int end = s[i] != 0 ? n : i ; 
	      
	      // save the max/start/stop
	      mxs.push_back( mx );
	      starts.push_back( start );
	      stops.push_back( end ); 
	      
	      in = false;
            }
	  else if ( in ) 
	    {
	      if ( x[i] > mx ) mx = x[i];
	    }
	}
      
      // just started one?
      if ( in )
	{
	  mxs.push_back( x[n-1] );
	  starts.push_back( n-1 );
	  stops.push_back( n ); // one-past-end
	}
            
      
      // 
      // Get distribution of peaks
      //

      if ( mxs.size() > 1 ) 
	{

	  // robust norms:

	  double median = MiscMath::median( mxs );
	  double iqr = MiscMath::iqr( mxs );
	  double robust_sd = 0.7413 * iqr;
	  
	  for (int j=0; j<mxs.size(); j++)
	    {	      
	      // wipe this region if above threshold
	      if ( mxs[j] > median + robust_sd * max ) 
		for (int k=starts[j]; k<stops[j]; k++)
		  s[k] = 0;
	    }
	}
      
    }


  //
  // Expanded regions (core vs flanking)
  //
  
  if ( threshold2 > 0 ) 
    {
      
      std::vector<int> s2 = s;

      bool in = s[0] != 0;

      int start= 0;
      
      for (int i=0; i<n; i++)
	{
	  
	  if ( s[i] != 0 && ! in )
	    {
	      // start new?
	      start = i;
	      in = true;
	    }
	  else if ( in && ( s[i] == 0 || i == n-1 ) )
	    {  
	      // use exact-end encoding here, not 1-past-end
              int end = s[i] != 0 ? n-1 : i-1 ;

	      int start2 = start;
	      int end2 = end;
	      in = false;
	      
	      // for (int ss=start - 10 ; ss < end + 10 ; ss++ ) 
	      //  	std::cout << "ss " << ss << "\t" << scaled[ss] << "\t" << ( ss >= start && ss <= end ) << "\n";
	      
	      while ( 1 )
		{
		  if ( start2 == 0 ) break;	
		  if ( scaled[ start2 - 1 ] < threshold2 ) break;	      		  
		  --start2;	
		  s2[start2] = s[start];
		}

	      while ( 1 )
		{
		  if ( end2 == n-1 ) break;
		  if ( scaled[end2+1] < threshold2 ) break;	      
		  ++end2;		  
		  s2[end2] = s[end];
		}
	      
	      // new putative event is start2 to end2 [inclusive]
	      // is this long enough? 

	      int fill = s[ start ];

	      if ( mindur2 > 0 && end2 - start2 + 1 < mindur2 ) 
		fill = 0;
	      
	      //std::cout << "fill = " << s[ start ] << " --> " << fill << " " << start << " " << end << " --> " << start2 << " " << end2 << "\n"; 
	      	      
	      // either wipes all, or sets all to an event
	      for (int j=start2; j<=end2; j++)
		{
		  //		  std::cout << " sc = " << scaled[j] << "\n";
		  s2[j] = fill;
		}

	    } // end of core interval	 	  
	} // next sample-point
      
      // copy back over
      s = s2;
    }
  
  
  //
  // Record as intervals?
  //

  if ( regions != NULL ) 
    {

      regions->clear();
      
      bool in = s[0] != 0;

      int start= 0;

      for (int i=1; i<n; i++)
        {

          if ( s[i] != 0 && ! in )
            {
              // start new?
	      start = i;
              in = true;
              
            }
          else if ( in && ( s[i] == 0 || i == n-1 ) )
            {
              // 1 past end
              int end = s[i] != 0 ? n : i ;
	      
	      // but record as sample-points (usual start, end+1 encoding)
	      regions->push_back( interval_t( start , end ) );

              in = false;
            }

        }
      
      // just started one?

      if ( in )
	regions->push_back( interval_t( start , n ) );
      
    }

  return s;
  
}



//
// disjoint sets
//


void MiscMath::disjoint_set_t::make_set(std::vector<int> const & universe)
{
  // create n disjoint sets (one for each item)
  for (int i: universe)
    {
      parent[i] = i;
      rank[i] = 0;
    }
}

// find the root of the set in which element `k` belongs
int MiscMath::disjoint_set_t::find(int k)
{
  // if `k` is not the root
  if (parent[k] != k)
    {
      // path compression
      parent[k] = find(parent[k]);
    }
  
  return parent[k];
}

// perform union of two subsets
void MiscMath::disjoint_set_t::make_union(int a, int b)
{
  // find the root of the sets in which elements x and y belongs
  int x = find(a);
  int y = find(b);
  
  // if x and y are present in the same set
  if (x == y) {
    return;
  }
  
  // always attach a smaller depth tree under the root of the deeper tree.
  if ( rank[x] > rank[y] ) {
    parent[y] = x;
  }
  else if (rank[x] < rank[y]) {
    parent[x] = y;
  }
  else {
    parent[x] = y;
    rank[y]++;
  }
}


void MiscMath::print_sets(std::vector<int> const &universe, disjoint_set_t &ds)
{
  for (int i: universe) {
    std::cout << ds.find(i) << ' ';
  }
  std::cout << "\n";
}

std::map<int,std::set<int> > MiscMath::get_sets( std::vector<int> const &universe, disjoint_set_t &ds )
{
  std::map<int,std::set<int> > r;
  for (int i: universe) 
    r[ ds.find(i) ].insert( i );
  return r;
}

    // usage: disjoint–Set data structure (Union–Find algorithm)

    // universe of items
    // vector<int> universe = { 1, 2, 3, 4, 5 };
    
    // initialize `DisjointSet` class
    // disjoint_set_t ds;
    
    // create a singleton set for each element of the universe
    // ds.make_set(universe);
    // print_sets(universe, ds);
 
    // ds.make_union(4, 3);        // 4 and 3 are in the same set
    // print_sets(universe, ds);
 

void MiscMath::winsorize( std::vector<double> * x , double p )
{
  // assume 'p' between 0 and 0.5
  if ( p < 0 || p > 0.5 ) 
    Helper::halt( "MiscMath::winsorize() with invalid p" );

  if ( p == 0 ) return;
  
  double lwr = MiscMath::percentile( *x , p );
  double upr = MiscMath::percentile( *x , 1-p );
  
  if ( lwr >= upr )
    Helper::halt( "should not happen...pls fix me" );
  
  for (int i=0; i<x->size(); i++)
    {
      if      ( (*x)[i] < lwr ) (*x)[i] = lwr;
      else if ( (*x)[i] > upr ) (*x)[i] = upr;
    }

}
