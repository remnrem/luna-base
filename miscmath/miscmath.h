
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

#ifndef __MISCMATH_H__
#define __MISCMATH_H__

#include <vector>
#include <cstddef>
#include <complex>
#include <algorithm>

#include <stdint.h>
#include <map>

namespace MiscMath
{
  
  // next pow2
  long int nextpow2( const int a );
  std::vector<double> logspace(double a, double b, int n);
  std::vector<double> linspace(double a, double b, int n);
  
  double rms( const std::vector<double> & );
  double rms( const double * , int n );
  
  // differences 
  std::vector<double> diff( const std::vector<double> & x );

  // nearest index [within lwr/upr range]
  int nearest_idx( const std::vector<double> & x , double value , int lwr = -1 , int upr = -1 );

  // mean/variance  
  double mean( const std::vector<double> & x );
  double mean( const std::vector<int> & x );
  std::complex<double> mean( const std::vector<std::complex<double> > & x );

  std::complex<double> max( const std::vector<std::complex<double> > & x );
  double empirical_pvalue( const double s , const std::vector<double> & x );
  
  double skewness( const std::vector<double> & x );
  double skewness( const std::vector<double> & x , double m , double sd );
  
  double median( const std::vector<double> & x );
  double iqr( const std::vector<double> & x );
  double percentile( const std::vector<double> & x , double p );
  
  template<typename T> static inline double Lerp(T v0, T v1, T t)
   {
     return (1 - t)*v0 + t*v1;
   }

 template<typename T> std::vector<T> quantile(const std::vector<T>& inData, const std::vector<T>& probs)
   {
     
     if (inData.empty())
       {
	 return std::vector<T>();
       }
     
     if (1 == inData.size())
       {
	 return std::vector<T>(1, inData[0]);
       }
     
     std::vector<T> data = inData;
     std::sort(data.begin(), data.end());
     std::vector<T> quantiles;
     
     for (size_t i = 0; i < probs.size(); ++i)
       {
	 T poi = MiscMath::Lerp<T>(-0.5, data.size() - 0.5, probs[i]);
	 
	 size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
	 size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));
	 
	 T datLeft = data.at(left);
	 T datRight = data.at(right);
	 
	 T quantile = Lerp<T>(datLeft, datRight, poi - left);
	 
	 quantiles.push_back(quantile);
       }
     
     return quantiles;
   }
 
 
 
 double meansq( const std::vector<double> & x );
 double variance( const std::vector<double> & x );  
 double variance( const std::vector<int> & x );  
 double variance( const std::vector<double> & x , double m );
 double variance( const std::vector<int> & x , double m );
 double sdev( const std::vector<double> & x );
 double sdev( const std::vector<double> & x , double m );
 double sum( const std::vector<double> & x );
 double sum( const std::vector<int> & x );
 double sqr( const double a );

 void minmax( const std::vector<double> & x , double * mn , double * mx);
 void normalize( std::vector<double> * x , double * mn , double * mx);
 void normalize( std::vector<double> * x , const std::vector<bool> & include_mask );
  
  std::vector<double> Z( const std::vector<double> & x );
  std::vector<double> logvector( const std::vector<double> & x );
  
  int outliers( const std::vector<double> * x , double th , 
		std::vector<bool> * inc , const std::vector<bool> * prior = NULL );  
  double kappa( const std::vector<int> & a , const std::vector<int> & b );
  double kappa( const std::vector<std::string> & a , const std::vector<std::string> & b );

  double mcc( std::map<std::string,std::map<std::string,int> > table ,
	      const std::vector<std::string> & labels );
  
  double accuracy( const std::vector<int> & a , const std::vector<int> & b ,
		   std::vector<int> * labels = NULL ,
		   std::vector<double> * precision = NULL ,
		   std::vector<double> * recall = NULL ,
		   std::vector<double> * f1 = NULL ,
		   double * macro_precision = NULL ,
		   double * macro_recall = NULL ,
		   double * macro_f1 = NULL ,		  
		   double * avg_weighted_precision = NULL ,
		   double * avg_weighted_recall = NULL ,
		   double * avg_weighted_f1 = NULL ,
		   double * mcc = NULL );

  double accuracy( const std::vector<std::string> & a , const std::vector<std::string> & b ,
		   std::vector<std::string> * labels = NULL , 
		   std::vector<double> * precision = NULL ,
		   std::vector<double> * recall = NULL ,
		   std::vector<double> * f1 = NULL ,
		   double * macro_precision = NULL ,
		   double * macro_recall = NULL ,
		   double * macro_f1 = NULL ,		   
		   double * avg_weighted_precision = NULL ,
		   double * avg_weighted_recall = NULL ,
		   double * avg_weighted_f1 = NULL ,
		   double * mcc = NULL );

  
  // p-values for F-test
  double pF(const double F, const int df1, const int df2);
  double betai(const double a, const double b, const double x);
  double betacf(const double a, const double b, const double x);
  
  // covariance
  double covariance( const std::vector<double> & x ,
		     const std::vector<double> & y ,
		     const int w = 1 );
  
  double clipped( const std::vector<double> & x );
  double clipped( const std::vector<double> & x , double, double);

  double flat( const std::vector<double> & x , double EPS = 1e-6 );
  double max( const std::vector<double> & x , double th );
  
  // Hjorth parameters
  void hjorth( const std::vector<double> * , double * , double * , double * );

  // second-order Hjorth parameters (window, inc)
  void hjorth2( const std::vector<double> * , double * , int w , int inc = 0 );

  // turning rate
  double turning_rate( const std::vector<double> * , int , int , int , std::vector<double> * sub );

  // test for over-dispersion of Poisson 
  double overdispersion( const std::vector<int> & , double * pv = NULL );
  double poisson( double, double );

  // Hann window
  std::vector<double> hann_window( int n );
  void hann_window( std::vector<double> * );
  double hann_window(unsigned int n, unsigned int N);

  // Hanning window (to match Matlab: same as Hann(n-2), but w/ 0 at ends)
  std::vector<double> hanning_window( int n );

  // Hamming window
  std::vector<double> hamming_window( int n );
  void hamming_window( std::vector<double> * );
  double hamming_window(unsigned int n, unsigned int N);

  // Tukey (cosine-taper) window
  std::vector<double> tukey_window( int n , double r = 0.5 );
  void tukey_window( std::vector<double> * , double r = 0.5 );

  // mean-centre
  std::vector<double> centre( const std::vector<double> & x );
  void centre( std::vector<double> * x );

  // detrend
  std::vector<double> detrend( const std::vector<double> & x , double * pa = NULL , double * pb = NULL );
  void detrend( std::vector<double> * x , double * pa = NULL , double * pb = NULL );

  // thresholding
  double threshold( const std::vector<double> & , double,double,double, double *, std::map<double,double> * t = NULL );


  // angles/phases

  double rad2deg(double radians);

  double deg2rad(double degrees);

  double shift_degrees( double degree , double x );

  double as_angle_0_pos2neg( const double r );

  // epoch/windowing

  int position2leftepoch( uint64_t , uint64_t , uint64_t , int e_total );
  int position2rightepoch( uint64_t , uint64_t , uint64_t , int e_total );

  // chi-sq goodness of fit test
  double chisq( const std::vector<double> & obs , const std::vector<double> & expected );

  // Moving average
  
  std::vector<double> moving_average( const std::vector<double> & x , int n );
  std::vector<double> moving_average_filter( const std::vector<double> & x , int n );
  
  // Median filter

  std::vector<double> median_filter( const std::vector<double> & x , int n );
  
  
  
  /*
   * Algorithm from N. Wirth's book, implementation by N. Devillard.
   * This code in public domain.
   */
  
  
  typedef double elem_type ;
  
#define ELEM_SWAP(a,b) { elem_type t=(a);(a)=(b);(b)=t; }


/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median() macro defined below to get the median. 

                Reference:

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

		  ---------------------------------------------------------------------------*/

  
  elem_type kth_smallest_destroy(elem_type a[], int n, int k);

  elem_type kth_smallest_preserve( const std::vector<elem_type> & a , int k );
  
#define median_destroy(a,n) MiscMath::kth_smallest_destroy(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))
#define median_preserve(a,n) MiscMath::kth_smallest_preserve(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))
#define percentile_preserve(a,n,k) MiscMath::kth_smallest_preserve(a,n,k);

  
}

#endif
  
