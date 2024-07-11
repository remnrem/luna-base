
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

    tmax = amax = rmax = 0;
    tmin = amin = rmin = 0;
    tminmax = aminmax = rminmax = 0;
  }
  
  double omean; // original input mean (post winsor/log)
  
  double sd;
  double mean;
  double cv;
  double tstat1; // based on simple epoch count
  double tstat2; // uses 'actual' (not clock) epoch count

  // max stats
  double tmax; // time from start to max (post smoothing) (epochs)
  double amax; // max amplitude (expressed as max - min) 
  double rmax; // amax / tmax (i.e. how steep the slope to max)
  // could also consider the product? (i.e. how sustained and large the increase?)

  // same for local min
  double tmin; 
  double amin; 
  double rmin; 

  // min vs max stats (i.e. two above relative to 0/start)
  double tminmax; 
  double aminmax; 
  double rminmax; 
  
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
  void set_norm_max( const bool b ) { norm01 = b; }
  void set_norm_mean( const bool b ) { norm_mean = b; }
  void set_norm_cycles( const bool b ) { norm_each_section = b; }
  bool norm_cycles() const { return norm_each_section; }
  void weight_cycles( const bool b ) { wcycles = b; } 
  void set_nq ( const int x ) { nq = x; } 
  void set_max_cycles( const int n );
  void set_cycles( const std::vector<int> n );
  
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

  bool norm01; // scale min-max to range 1.0 
  bool norm_mean; // mean-center each trace
  bool norm_each_section; // repeat norm for each cycle

  bool wcycles; // weight each cycle by # epochs for BETWEEN
  
  std::set<std::string> incl_cycles;
  
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
			 const std::vector<int> & ee ,
			 const bool do_smoothing ,  
			 const bool do_norming );
public:
  static void output_helper( const qdynam_results_t & res , const bool verbose , const bool between = false );

  static std::vector<double> qnt( const std::vector<double> & x , const int nq = 10 );

  static std::vector<double> smooth( const std::vector<double> & x , const std::vector<int> & e ,
				     const int w1, const int w2 );

  static void norm( std::vector<double> * x , const bool do_max , const bool do_mean );
  
};


#endif
