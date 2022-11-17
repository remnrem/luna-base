
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

#ifndef __DYNAM_H__
#define __DYNAM_H__

#include <vector>
#include <map>
#include <set>
#include <string>

// wrapper
void dynam_report( const std::vector<double> & y , 
		   const std::vector<double> & t , 
		   const std::vector<std::string> * g = NULL ); 

void dynam_report_with_log( const std::vector<double> & y , 
			    const std::vector<double> & t , 
			    const std::vector<std::string> * g = NULL ); 


struct dynam_t { 

  dynam_t() { }
  dynam_t( const std::vector<double> & y );
  dynam_t( const std::vector<double> & y , const std::vector<double> & t );
  dynam_t( const std::vector<double> & y , const std::vector<int> & t );
  
  int size() const { return y.size(); } 
  
  void clear() 
  {
    y.clear();
    t.clear();
  }

  // for a dynamic time series, calculate the following properties
  
  void denoise( double lambda );
  
  bool mean_variance( double * mean , double * var );
  
  bool linear_trend( double * beta , double * rsq , double * intercept = NULL );
  
  void hjorth( double * h1 , double *h2 , double *h3 );
    
  // data
  
  std::vector<double> y;

  std::vector<double> t;
  
};


// group (i.e. sleep cycle) dynamics 

struct gdynam_t { 
  
  // between and within 'group' dynamics

  gdynam_t( const std::vector<int> & g , const std::vector<double> & y );
  gdynam_t( const std::vector<int> & g , const std::vector<double> & y , const std::vector<double> & t );
  gdynam_t( const std::vector<int> & g , const std::vector<double> & y , const std::vector<int> & t );

  void clear()
  {
    g.clear();
    y.clear();
    t.clear();
    gmap.clear();
    within.clear();
    between.clear();
  }

  int stratify();
  
  // linear trends: within groups
  // differences between groups;
  
  // data
  
  std::vector<int> g;
  
  std::vector<double> y;

  std::vector<double> t;
  
  // group mapping
  std::map<int,std::set<int> > gmap;
  
  // within-group effects
  std::map<int,dynam_t> within;
  dynam_t between;
  
};

struct dissipation_t 
{
  
  dissipation_t( const std::vector<double> & x , 
		 const int mx = 0 ,  
		 const double winsor = 0.05 );
  
  std::vector<double> plife( const std::vector<double> & ps );
  
  std::vector<double> s;

};


#endif
