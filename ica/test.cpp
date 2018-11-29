#include <iostream>
#include "ica.h"
#include <vector>

int main()
{
  std::vector<std::vector<double> > X ;

  int n = 2 ;

  ica_t ica( X , n );
  
  std::cout << ica.S[0][1] << "\n";

}
