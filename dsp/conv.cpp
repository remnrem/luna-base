
#include "conv.h"

#include <cstddef>
#include <cstdio>


std::vector<double> dsptools::convolve( const std::vector<double> & signal , 
					const std::vector<double> & kernel )
{

  const int nsig = signal.size();
  const int nkern = kernel.size();
  const int nconv = nsig + nkern - 1 ; 

  std::vector<double> result( nconv , 0 );

  for (int n=0; n < nconv; n++ )
    {
      
      size_t kmin, kmax;

      kmin = (n >= nkern - 1) ? n - (nkern - 1) :  0;
      kmax = (n < nsig - 1)   ? n               :  nsig - 1;

      for (size_t k = kmin; k <= kmax; k++)
	result[n] += signal[k] * kernel[n - k];       
    }
  
  return result;

}

