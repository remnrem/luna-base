
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


#ifndef __CLOCS_H__
#define __CLOCS_H__

#include <cmath>
#include <map>
#include <vector>

#include "stats/statistics.h" 

struct edf_t;
struct cart_t;

struct signal_list_t;

struct sph_t { 
  sph_t() { } 
  sph_t( double a, double e, double r ) : azimuth(a) , elevation(e) , r(r) { } 
  double azimuth,elevation,r;
  cart_t cart() const;
};

struct cart_t { 
  cart_t() { } 
  cart_t(double x, double y, double z) : x(x) , y(y) , z(z) { } 
  double x, y, z;  
  sph_t sph() const;
  
};

struct polar_t {
  polar_t() { } 
  polar_t( sph_t & sph ); 
  cart_t cart() const { return cart_t( radius * cos( angle ) , radius * sin( angle ) , 0 ); } 
  
  double angle, radius; 

};



struct clocs_t { 
  
  edf_t * edf;

  void add_cart( const std::string & label , double x , double y , double z ) 
  { 
    cloc[ label ] = cart_t( x , y , z ) ;
  } 
  
  std::map<std::string,cart_t> cloc;
  

  void convert_to_unit_sphere();

  Data::Matrix<double> interelectrode_distance_matrix( const signal_list_t & signals ) const; 
						       
  Data::Matrix<double> interelectrode_distance_matrix( const signal_list_t & signals1 , 
						       const signal_list_t & signals2 ) const;

  bool make_interpolation_matrices( const signal_list_t & good_signals , 
				    const signal_list_t & bad_signals , 
				    Data::Matrix<double> * G , 
				    Data::Matrix<double> * Gi );
  
  
  Data::Matrix<double> interpolate( const Data::Matrix<double> & data , 
				    const std::vector<int> & good_channels , 
				    const Data::Matrix<double> & G , 
				    const Data::Matrix<double> & Gi );

    
  void attach( edf_t * edf_p ) 
  {
    edf = edf_p;
  }
  
  

  bool has( const std::string & cl ) const
  {
    return cloc.find( cl ) != cloc.end() ;
  }
  
  cart_t cart( const std::string & cl ) const
  {
    if ( ! has( cl ) ) Helper::halt( "did not have map position for " + cl ) ;
    return cloc.find( cl )->second; 
  }

  sph_t sph( const std::string & cl ) const
  {
    if ( ! has( cl ) ) Helper::halt( "did not have map position for " + cl ) ;
    return cloc.find( cl )->second.sph(); 
  }
  
  int load_cart( const std::string & filename );



  //
  // Helper functions
  //

  static inline double rad2deg(double radians) { return radians * (180.0 / M_PI); }

  static inline double deg2rad(double degrees) { return degrees * (M_PI / 180.0); }

  static cart_t sph2cart( const sph_t & sph )
  {
    cart_t c;
    c.x = sph.r * cos(sph.elevation) * cos(sph.azimuth);
    c.y = sph.r * cos(sph.elevation) * sin(sph.azimuth);
    c.z = sph.r * sin(sph.elevation);
    return c;
  }

  static sph_t cart2sph( const cart_t & cart )
  {
    sph_t sph;
    sph.azimuth = atan2( cart.y , cart.x );
    sph.elevation = atan2( cart.z, sqrt( cart.x*cart.x + cart.y*cart.y) );
    sph.r = sqrt( cart.x*cart.x + cart.y*cart.y + cart.z*cart.z);
    return sph;

    // as per Matlab : 
    // The notation for spherical coordinates is not standard. For the
    // cart2sph function, elevation is measured from the x-y
    // plane. Notice that if elevation = 0, the point is in the x-y
    // plane. If elevation = pi/2, then the point is on the positive z-axis.																			      
  }
  
  
  
  
};


#endif
