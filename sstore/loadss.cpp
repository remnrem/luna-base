#include "sstore.h"
#include <iostream>
#include "../helper/helper.h"

//
// loadss, a simple loader for a sstore_t 
// 

int main(int argc , char ** argv )
{

   if ( 0 ) 
     {
       
       sstore_t ss( "ss1.db" );

       std::map<int,sstore_data_t> epochs = ss.fetch_epochs();
       
       std::cout << "sz = " << epochs.size() << "\n";

       std::map<int,sstore_data_t>::const_iterator ii = epochs.begin();

       while ( ii != epochs.end() )
	 {
	   
	   std::cout << "epoch " << ii->first << "\n";
	   	   	   
	   const std::map<sstore_key_t,sstore_value_t> & datum   = ii->second.data;
	   
	   std::map<sstore_key_t,sstore_value_t>::const_iterator kk = datum.begin();
	   while ( kk != datum.end() )
	     {
	       
	       std::cout << kk->first.id << " (" << kk->first.ch << ") ";
	       if      ( kk->second.is_text ) std::cout << " str = " << kk->second.str_value << "\n";
	       else if ( kk->second.is_double ) std::cout << " dbl = " << kk->second.dbl_value << "\n";
	       else {
		 std::cout << " vec[" << kk->second.vec_value.size() << "]";
		 for (int i=0;i<kk->second.vec_value.size();i++) std::cout << " " << kk->second.vec_value[i] ;
		 std::cout << "\n";
	       }
	       
	       ++kk;
	       
	     }
	   
	   ++ii;
	 }
       
       std::exit(0);
     }

  
  if ( argc != 3 ) 
    { 
      std::cerr << "usage: ./loadss {ss.db} {strata}  < input\n"
		<< "where ss.db  --> sstore_t database file\n"
		<< "      strata --> [-a|-e|-i|index|unindex] to specify all/epoch/interval data\n"	
		<< "\n";          
	std::exit(1); 
    } 
  

  // 
  // Input format, tab-delimited
  //
  
  // all      :   ID LVL CH              N VALUE(S)
  // epoch    :   ID LVL CH  E           N VALUE(S)
  // interval :   ID LVL CH  START STOP  N VALUE(S)


  // LVL and CH are optional ( set to . if missing)

  std::string filename = argv[1];
  std::string mode     = argv[2];
  
  // special case of indexing/dropping index
  
  if ( mode == "index" ) 
    {
      sstore_t ss( filename );
      ss.index();
      std::exit(0);
    }

  if ( mode == "unindex" ) 
    { 
      sstore_t ss( filename );
      ss.drop_index();
      std::exit(0);
    }

  
  bool mode_baseline = mode == "-a";
  bool mode_epoch    = mode == "-e";
  bool mode_interval = mode == "-i";
  
  if ( ! ( mode_baseline || mode_epoch || mode_interval ) ) 
    {
      Helper::halt( "mode argument should be -a, -e or -i" );
    }
  
  //
  // Open/create sstore_t
  //

  sstore_t ss( filename );

  ss.drop_index();
  
  while ( ! std::cin.eof() ) 
    {
      std::string line;
      std::getline( std::cin , line , '\n' );
      if ( std::cin.eof() ) break;
      std::vector<std::string> tok = Helper::parse( line , "\t" );

      const int t = tok.size();
      if ( tok.size() == 0 ) continue;
      
      //
      // Baseline-level inputs
      //
      
      if ( mode_baseline )
	{
	  
	  if ( t < 5 ) Helper::halt( "base: format problem:\n" + line );
	  
	  int n;
	  if ( ! Helper::str2int( tok[3] , &n ) ) 
	    Helper::halt( "format problem" );
	  
	  int expected = 5; 
	  if ( n > 1 ) expected += n - 1; 
	  if ( t != expected ) Helper::halt( "format problem:\n" + line );
	  	  
	  bool has_level = tok[1] != ".";
	  bool has_channel = tok[2] != ".";
	  
	  const std::string * level_ptr = has_level ? &(tok)[1] : NULL ; 
	  const std::string * channel_ptr = has_channel ? &(tok)[2] : NULL ; 
	  
	  if ( n == 0 ) // text
	    {
	      ss.insert_base( tok[0] , tok[4] , channel_ptr , level_ptr );
	    }
	  else if ( n == 1 ) // double 
	    {
	      double d;
	      if ( ! Helper::str2dbl( tok[4] , &d ) ) 
		Helper::halt( "format problem, expecting double:\n" + line );
	      
	      ss.insert_base( tok[0] , d , channel_ptr , level_ptr );
	      
	    }
	  else // array of doubles
	    {
	      std::vector<double> d( n , 0 );
	      for (int i=0;i<n;i++)
		if ( ! Helper::str2dbl( tok[4+i] , &(d)[i] ) ) 
		  Helper::halt( "format problem, expecting double:\n" + line );

		ss.insert_base( tok[0] , d , channel_ptr , level_ptr );
	      
	    }
	}
      
      //
      // Epoch-level inputs 
      //

      if ( mode_epoch ) 
	{
	  
	  if ( t < 6 ) Helper::halt( "format problem:\n" + line );
	  
	  int n;
	  if ( ! Helper::str2int( tok[4] , &n ) ) 
	    Helper::halt( "format problem:\n" + line );
	  
	  int e;
	  if ( ! Helper::str2int( tok[3] , &e ) ) 
	    Helper::halt( "format problem:\n" + line  );

	  int expected = 6; 
	  if ( n > 1 ) expected += n - 1; 
	  if ( t != expected ) Helper::halt( "format problem:\n" + line );
	  
	  bool has_level   = tok[1] != ".";
	  bool has_channel = tok[2] != ".";

	  const std::string * level_ptr = has_level ? &(tok)[1] : NULL ; 
	  const std::string * channel_ptr = has_channel ? &(tok)[2] : NULL ; 
	  
	  if ( n == 0 ) // text
	    {
	      ss.insert_epoch( e , tok[0] , tok[5] , channel_ptr , level_ptr );	      
	    }
	  else if ( n == 1 ) // double 
	    {
	      double d;
	      if ( ! Helper::str2dbl( tok[5] , &d ) ) 
		Helper::halt( "format problem, expecting double:\n" + line );
	      
	      ss.insert_epoch( e, tok[0] , d , channel_ptr , level_ptr );
	      
	    }
	  else // array of doubles
	    {
	      std::vector<double> d( n , 0 );
	      for (int i=0;i<n;i++)
		if ( ! Helper::str2dbl( tok[5+i] , &(d)[i] ) ) 
		  Helper::halt( "format problem, expecting double:\n" + line );

	      ss.insert_epoch( e, tok[0] , d , channel_ptr , level_ptr );
	      
	    }
	  
	}



      //
      // Interval-level inputs
      //

      if ( mode_interval )
	{
	  
	  if ( t < 7 ) Helper::halt( "format problem:\n" + line );
	  
	  int n;
	  if ( ! Helper::str2int( tok[5] , &n ) ) 
	    Helper::halt( "format problem:\n" + line );
	  
	  uint64_t a;
	  if ( ! Helper::str2int64( tok[3] , &a ) ) 
	    Helper::halt( "format problem:\n" + line  );

	  uint64_t b;
	  if ( ! Helper::str2int64( tok[4] , &b ) ) 
	    Helper::halt( "format problem:\n" + line  );

	  int expected = 7; 
	  if ( n > 1 ) expected += n - 1; 
	  if ( t != expected ) Helper::halt( "format problem:\n" + line );
	  
	  bool has_level   = tok[1] != ".";
	  bool has_channel = tok[2] != ".";
	  
	  const std::string * level_ptr = has_level ? &(tok)[1] : NULL ; 
	  const std::string * channel_ptr = has_channel ? &(tok)[2] : NULL ; 
	  
	  if ( n == 0 ) // text
	    {
	      ss.insert_interval( a , b , tok[0] , tok[6] , channel_ptr , level_ptr );
	    }
	  else if ( n == 1 ) // double 
	    {
	      double d;
	      if ( ! Helper::str2dbl( tok[6] , &d ) ) 
		Helper::halt( "format problem, expecting double:\n" + line );
	      
	      ss.insert_interval( a, b , tok[0] , d , channel_ptr , level_ptr );

	    }
	  else // array of doubles
	    {
	      std::vector<double> d( n , 0 );
	      for (int i=0;i<n;i++)
		if ( ! Helper::str2dbl( tok[6+i] , &(d)[i] ) ) 
		  Helper::halt( "format problem, expecting double:\n" + line );
	      
	      ss.insert_interval( a, b , tok[0] , d , channel_ptr , level_ptr );
	      
	    }
	  
	}
      
      // next row of input
    }
    

  ss.index();

  // all      :   ID CH LVL             N VALUE(S)
  // epoch    :   ID CH LVL E           N VALUE(S)
  // interval :   ID CH LVL START STOP  N VALUE(S)


  std::exit(0);
}
