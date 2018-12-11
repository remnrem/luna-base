#ifndef __LUNA_RETVAL_H__
#define __LUNA_RETVAL_H__

#include <map>
#include <set>
#include <string>
#include <vector>

// when working with Luna via an API or through R, this structure
// provides a way to return results from multiple commands for
// a single dataset (i.e. edf/individual)

// it is designed to be a plug-in for writer, 
// i.e. writer.var() and writer.value(), writer.level(), 
// writer.unlevel(), and writer.epoch()

// hopefully this will become a better, more light-weight 
// replacement for destrat

struct retval_t { 
  
  // cmd --> {table/strata} --> { data-table }  
  //        
  // data-table , i.e. long format
  //    STRATUM1 STRATUM2 ...  VAR  VALUE 

  // strata must be unique 

};

struct retval_strata_t { 
  
  std::set<std::string> factors;
  
  bool operator< (const retval_strata_t & rhs ) const { 
    if ( factors.size() < rhs.factors.size() ) return true;
    if ( factors.size() > rhs.factors.size() ) return false;
    std::set<std::string>::const_iterator ii = factors.begin();
    std::set<std::string>::const_iterator jj = rhs.factors.begin();
    while ( ii != factors.end() ) 
      {
	if ( *ii < *jj ) return true;
	if ( *jj < *ii ) return false;
	++ii;
	++jj;
      }

    return false;
  }

};

struct retval_dblvar_t { 
  std::string var;
  double value;
};

struct retval_intvar_t { 
  std::string var;
  int value;
};

struct retval_strvar_t { 
  std::string var;
  std::string value;
};

struct retval_data_t { 
  
  std::map<retval_strata_t,retval_dblvar_t> dbl_data; 
  std::map<retval_strata_t,retval_intvar_t> int_data; 
  std::map<retval_strata_t,retval_strvar_t> str_data; 

};

#endif
