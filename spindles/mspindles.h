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


#ifndef __MSPINDLES_H__
#define __MSPINDLES_H__


struct spindle_t;
struct signal_list_t;

struct mspindle_t { 

  mspindle_t()
  {
    spindles.clear();
    run.clear();
    frq = lwr_frq = upr_frq = 0;
    stat = 0;    
  }
  
  void add( const spindle_t * s , const int r , const std::string & l ) { spindles.push_back( s ); run.push_back(r); lab.push_back(l);  }
  
  void summarize() 
  { 
    const int ns = spindles.size();
    if ( ns == 0 ) return;
    
    start = spindles[0]->tp.start;
    stop = spindles[0]->tp.stop;
    

    stat = 0;
    std::vector<double> wgt( ns );
    
    for (int i=0;i<ns;i++) 
      stat += spindles[i]->mean_stat;

    for (int i=0;i<ns;i++) wgt[i] = spindles[i]->mean_stat / stat;
    
    stat /= (double)ns;

    // weighted mean for frequency

    lwr_frq = spindles[0]->fft;
    upr_frq = spindles[0]->fft;
    frq = 0;

    for (int i=0;i<ns;i++) 
      {
	frq += wgt[i] * spindles[i]->fft;
	if ( spindles[i]->fft < lwr_frq ) lwr_frq = spindles[i]->fft;
	if ( spindles[i]->fft > upr_frq ) upr_frq = spindles[i]->fft;

	if ( spindles[i]->tp.start < start ) start = spindles[i]->tp.start;
	if ( spindles[i]->tp.stop > stop ) stop = spindles[i]->tp.stop;
      }
    
    
  }

  std::vector<const spindle_t*> spindles;
  std::vector<int> run;
  std::vector<std::string> lab;

  uint64_t start, stop;

  double frq;
  double lwr_frq;
  double upr_frq;
  
  double stat;

  double dur() const { return ( stop - start + 1 ) / (double)globals::tp_1sec ; }

  double n() const { return spindles.size(); } 
  
};


struct sort_t;

struct mspindles_t 
{ 
  
  mspindles_t() { } 

  mspindles_t( edf_t * edf ) : edf(edf) 
  { 
    interval_th = 0; // proportional overlap
    cross_ch_interval_th = 0;
    within_ch_interval_th = 0;
    frq_th = 1; // Hz
  } 
  
  // parameters

  bool hms;

  bool per_spindle_verbosity;

  // generic interval threshold
  double interval_th ;
  
  // interval threshold applied to intervals of same Fc, but across channels
  double cross_ch_interval_th;

  // interval threshold applied to intervals of different Fc, but across channels
  double within_ch_interval_th;

  // interval window (in sec)
  double window;

  // frequency threshold, i.e. merge two spindles if within t Hz
  double frq_th;

  void add( const std::vector<spindle_t> & spindles , 
	    int fs , uint64_t len , 
	    int fc , int ch , 
	    const std::string & label );

  void collate();
  void proc_overlaps( const std::vector<sort_t> & );
  void output( const signal_list_t & );
  void plot( const std::string & );
  void pairwise_statistics( int, int );
  
private:

  edf_t * edf;

  // all spindle intervals (w/ stats populated)
  std::vector<std::vector<spindle_t> > S;
  
  // Mins duration, i.e. denominator for each interval set
  std::vector<double> mins;

  // track channel and Fc for each set
  std::vector<int> ch;
  std::vector<double> frq;
  std::vector<std::string> run_label;

  // collate() creates this set of combined spindles
  std::vector<mspindle_t> M;
  
};



// helper struct, for grouping all spindles

struct sort_t { 

  // for sortign
  interval_t i; 
  double f; 
  int ch;
  int run; // i.e. could be different thresholds for same 'f' and 'ch', so also track unique 'run' number
  std::string label;  //for output only

  // track back to original spdinle
  const spindle_t * spindle;
  
sort_t( const interval_t & i, const double f, const int ch , const int run, const std::string & label , const spindle_t * spindle ) 
: i(i), f(f), ch(ch), spindle(spindle), label(label), run(run)  { } 

  bool operator<( const sort_t & rhs ) const 
  {
    if ( i < rhs.i ) return true;
    if ( i > rhs.i ) return false;
    if ( run < rhs.run ) return true;
    if ( run > rhs.run ) return false;
    if ( ch < rhs.ch ) return true;
    if ( ch > rhs.ch ) return false;
    return f < rhs.f;
  }  
};



void intersect_spindles( edf_t & edf , 
			 const std::string & la, 
			 const std::vector<spindle_t> & a , 
			 const std::string & lb, 
			 const std::vector<spindle_t> & b ,
			 double th , 
			 double win,
			 const std::vector<double> & save_t_minutes ) ;


#endif
