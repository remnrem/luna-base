
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

// helper function

struct edf_t;
struct param_t;

bool dynam_compile_cycles( edf_t & edf , std::vector<std::string> * , std::vector<int> * );

// wrapper
void dynam_report( param_t & param,
		   const std::vector<double> & y , 
		   const std::vector<double> & t , 
		   const std::vector<std::string> * g = NULL ); 

void dynam_report_with_log( param_t & param,
			    const std::vector<double> & y , 
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


// separate class for quantitative differences

struct qdynam_results_t
{

  qdynam_results_t()
  {
    sd = 0;
    omean = mean = 0;
    cv = 0;
    tstat1 = tstat2 = 0;
    ne = 0;

    tmax = amax = lmax = rmax = 0;
    tmin = amin = lmin = rmin = 0;
    tminmax = aminmax = lminmax = rminmax = 0;
  }
  
  double omean; // original input mean (post winsor/log)
  
  double sd;
  double mean;
  double cv;
  double tstat1; // based on simple epoch count
  double tstat2; // uses 'actual' (not clock) epoch count

  double tmax; // time from start to max (post smoothing) (epochs)
  double amax; // max amplitude (expressed as max - min) 
  double lmax; // amax * tmax (i.e. how sustained and large the increase)
  double rmax; // amax / tmax (i.e. how quick to max) 

  // as above, but for mins
  double tmin; // time from start to max (post smoothing) (epochs)
  double amin; // max amplitude (expressed as max - min) 
  double lmin; // amax * tmax (i.e. how sustained and large the increase)
  double rmin; // amax / tmax (i.e. how quick to max) 

  // min vs max
  double tminmax; // time from start to max (post smoothing) (epochs)
  double aminmax; // max amplitude (expressed as max - min) 
  double lminmax; // amax * tmax (i.e. how sustained and large the increase)
  double rminmax; // amax / tmax (i.e. how quick to max) 

  int ne; // number of epochs included  
  
};

struct qdynam_t
{
  
  qdynam_t( const int ne , const std::vector<std::string> * cycles = NULL ) ;

  void winsorize( const double p );
  void log_transform( const bool b ) { logscale = b; }
  void set_smoothing_median_window( const int w ) { median_window = w; } 
  void set_smoothing_mean_window( const int w ) { mean_window = w; } 
  void include( const std::vector<int> & );
  void include( const std::vector<bool> & );
  void set_epochs(  const std::vector<int> & );
  void set_min_ne( const int x ) { min_ne = x ; } 
  void proc( const std::vector<double> & x );
  void norm( const bool b ) { norm01 = b; }
  void weight_cycles( const bool b ) { wcycles = b; } 
  void set_nq ( const int x ) { nq = x; } 
  std::vector<double> results() const;
  
  std::map<std::string,std::vector<double> > stratified_results() const;
  
  qdynam_results_t r1; // all
  qdynam_results_t rb; // between cycle
  qdynam_results_t rwa; // average of within-cycle results
  std::map<std::string,qdynam_results_t> rw; // within cycle

  // smoothed & normed series
  std::vector<double> r1_smoothed_series; // copy for TOT
  std::map<std::string,std::vector<double> > rw_smoothed_series; // within cycle
  std::map<std::string,std::vector<int> > rw_epochs; // within cycle

  // quantile traces (ss)
  std::vector<double> r1_q10; 
  std::map<std::string,std::vector<double> > rw_q10; // within cycle  

  // quantile traces (os)
  // repeats, for original (smoothed) series
  std::vector<double> r1_os_q10; 
  std::map<std::string,std::vector<double> > rw_os_q10; // within cycle  
  
private:
  
  std::vector<std::string> cycles;  
  bool has_cycles;  
  int ne; 
  int min_ne; // to include an cycle in the within-mean
  int median_window;
  int mean_window;
  bool norm01;
  bool wcycles;
  int nq;
  
  std::vector<double> ss; // smoothed series (normed)
  std::vector<double> os; // original series (unnormed - but smoothed)
  
  std::vector<bool> incl;
  std::vector<int> epochs;

  double winsor;
  bool logscale;
  
  // results
  
  // main calc function
  qdynam_results_t calc( const std::vector<double> & xx ,
			 const std::vector<double> & ox ,
			 const std::vector<int> & ee ,
			 const bool skip_smoothing = false );
public:
  static void output_helper( const qdynam_results_t & res , const bool verbose );

  static std::vector<double> qnt( const std::vector<double> & x , const int nq = 10 );

  static std::vector<double> smooth( const std::vector<double> & x , const int w1, const int w2 );
  
};


#endif
