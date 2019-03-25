
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

  double rms( const std::vector<double> & );
  double rms( const double * , int n );

  // differences 
 std::vector<double> diff( const std::vector<double> & x );

  // mean/variance  
 double mean( const std::vector<double> & x );
 double mean( const std::vector<int> & x );
 std::complex<double> mean( const std::vector<std::complex<double> > & x );
 
 double median( const std::vector<double> & x );
 double iqr( const std::vector<double> & x );

 
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
  
  std::vector<double> Z( const std::vector<double> & x );
  std::vector<double> logvector( const std::vector<double> & x );

	      
  // covariance
  double covariance( const std::vector<double> & x ,
		     const std::vector<double> & y ,
		     const int w = 1 );
  
  double clipped( const std::vector<double> & x );
  double clipped( const std::vector<double> & x , double, double);
  
  // Hjorth parameters
  void hjorth( const std::vector<double> * , double * , double * , double * );

  // turning rate
  double turning_rate( const std::vector<double> * , int , int , int , std::vector<double> * sub );

  // test for over-dispersion of Poisson 
  double overdispersion( const std::vector<int> & , double * pv = NULL );
  double poisson( double, double );

  // Hann window
  std::vector<double> hann_window( int n );
  void hann_window( std::vector<double> * );
  double hann_window(unsigned int n, unsigned int N);
  
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
  
#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }


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
  
}

#endif
  
