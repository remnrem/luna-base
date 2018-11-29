#include "sstore.h"
#include <iostream>
#include "../helper/helper.h"

//
// loadints, a really simple interval loader for a sstore_t 
// 

int main(int argc , char ** argv )
{

  if ( argc != 3 ) 
    { 
      std::cerr << "usage: ./loadints {filename} {label} < < input\n"
		<< "where input is the interval ranges (uint64_t tp)"
		<< "\n";          
	std::exit(1); 
    } 
  

  // 
  // Input format, tab-delimited
  //
  
  // interval :  START STOP

  std::string filename = argv[1];
  std::string label    = argv[2];
  
  //
  // Open/create sstore_t
  //
  
  sstore_t ss( filename );
  
  ss.drop_index();
  
  while ( ! std::cin.eof() ) 
    {
      uint64_t a, b;
      std::cin >> a >> b;
      if ( std::cin.eof() ) break;
      
      // only label (as name) : no value, or channel/level stratifiers
      ss.insert_interval( a , b , label , "." , NULL , NULL );
      
      // next row of input
    }
  
  ss.index();

  std::exit(0);
}
