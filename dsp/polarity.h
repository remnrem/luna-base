#ifndef __POLARITY_H__
#define __POLARITY_H__

#include "edf/edf.h"
#include "intervals/intervals.h"
#include <vector>


namespace dsptools 
{ 

  void polarity( edf_t & edf , const param_t & param );
  
  std::vector<bool> make_mask( const std::vector<double> & x , double th );
  
  void polarity_check( const std::vector<double> & x0 , const std::vector<uint64_t> * tp , int fs , 
		       double th ,  // threshold for extracted filtered segments
		       bool zc2zc , // extract around peaks up to ZC's
		       double flim , // calculate PSD up to flim Hz only
		       double f_lwr , // lower BPF transition frequency
		       double f_upr , // upper 
		       bool mirror_mode ,   // mirror odd numbered up/down segments, i.e. make 'wave-like' signal
		       bool double_up , // instead of mirror-mode alternate segments, double-enter each (up and down)
		       bool analyse_bpf_signal , // analysis the BPF signal, not raw data
		       bool dmode // segment 'upward' and 'downward' rather than pos and neg
		       );


  void ht_polarity_check( const std::vector<double> & x0 , 
			  const std::vector<uint64_t> * tp , 
			  int fs , 
			  double f_lwr , // lower BPF transition frequency
			  double f_upr // upper 
			  );
    
}



#endif
