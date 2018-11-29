#ifndef __CONV_H__
#define __CONV_H__

#include <vector>

namespace dsptools 
{ 
  
  std::vector<double> convolve( const std::vector<double> & signal , 
				const std::vector<double> & kernel ); 
  
}

#endif

