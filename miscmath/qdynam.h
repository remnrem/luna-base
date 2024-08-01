
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

#ifndef __QDYNAM_H__
#define __QDYNAM_H__

#include <vector>
#include <map>
#include <set>
#include <string>

// helper function

struct edf_t;
struct param_t;

// separate class for quantitative differences

struct qdynam_results_t
{

  qdynam_results_t()
  {
    init();
  }
  
  void init()
  {
    sd = 0;
    omean = mean = 0;
    cv = 0;
    corr1 = corr2 = 0;

    lma1 = 0; lmb2 = 0; // time-based                                                                                                                        
    r_lma1 = 0; r_lmb2 = 0; // rank-based                                                                                                                    
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

  double corr1; // rank-order corr
  double corr2; // as tstat2

  double lma1, lmb2;
  double r_lma1, r_lmb2; 
  
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

  // usage:

  // ctr, init
  qdynam_t( ) { } ;
  void init( edf_t & , param_t & );
    
  //   loop over channels, freqs, etc,
  //   add stratified epoch-level metrics at time of compute)
  void add( const std::map<std::string,std::string> & faclvl, 
	    const std::string & metric,
	    const int epoch ,
	    const double value );
  
  //  at end of rutine, call proc() which does over all prior strata
  void proc_all();

  
  
  //
  // internal functions
  //

private:
  
  //  once at start run qdynam_t::dynam_compile_cycles()
  bool dynam_compile_cycles( edf_t & edf );

  void reinit(); // called when starting a new proc()
  
  // main calc
  void proc( const std::vector<double> & x ,
	     const std::vector<int> & e,
	     const std::vector<std::string> & c );

  void winsorize( const double p );
  void log_transform( const bool b ) { logscale = b; }
  void set_smoothing_median_window( const int w ) { median_window = w; } 
  void set_smoothing_mean_window( const int w ) { mean_window = w; } 
  void set_min_ne( const int x ) { min_ne = x ; } 
  void set_trim_epochs( const std::vector<int> & x ) { if ( x.size() == 2 ) trim_epochs = x; }
  void set_norm_max( const bool b ) { norm01 = b; }
  void set_norm_mean( const bool b ) { norm_mean = b; }
  void set_norm_cycles( const bool b ) { norm_each_section = b; }
  bool norm_cycles() const { return norm_each_section; }
  void set_weight_cycles( const bool b ) { wcycles = b; } 
  void set_nq ( const int x ) { nq = x; } 
  void set_max_cycles( const int n );
  void set_cycles( const std::vector<int> n );
  

  //
  // cache results
  //

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

  // faclvl -> var -> epoch -> values (allows for different vars to have diff epoch counts,
  //  e.g. due to different epoch outlier procedures of definitions) 
  std::map<std::map<std::string,std::string>,std::map<std::string,std::map<int,double> > > sequences;  
  std::vector<std::map<std::string,std::string> > osequences;
  
  std::map<int,std::string> cycles;   // epoch2cycle
  bool has_cycles;  

  int min_ne; // to include an cycle in the within-mean
  std::vector<int> trim_epochs; // remove X and Y epochs on either side of a cycle (within only)
  int median_window;
  int mean_window;

  bool norm01; // scale min-max to range 1.0 
  bool norm_mean; // mean-center each trace
  bool norm_each_section; // repeat norm for each cycle

  bool wcycles; // weight each cycle by # epochs for BETWEEN
  
  std::set<std::string> incl_cycles; // uniq cycles to include
  
  int nq;
  
  std::vector<double> ss; // smoothed series (normed)
  std::vector<double> os; // original series (unnormed - but smoothed)
  
  std::vector<int> epochs; // built at start; given seq can be a subset, made per proc
  std::set<int> uepochs;   // track unique epochs - added() epochs must in already present
  
  double winsor;
  bool logscale;

  bool verbose;
  bool epoch_output;
  
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
