/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2001    Robert Gentleman, Ross Ihaka 
 *                             and the R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 * Memory Allocation (garbage collected) --- INCLUDING S compatibility ---
 */

#include <vector>
#include "helper/helper.h"

#ifndef __FISHER_H__
#define __FISHER_H__

#ifndef R_EXT_MEMORY_H_
#define R_EXT_MEMORY_H_

#ifdef  __cplusplus
extern "C" {
#endif

  char * vmaxget(void);
  void vmaxset(char*);

  void R_gc(void);

  char * R_alloc(long, int);
  char * S_alloc(long, int);
  char * S_realloc(char*, long, long, int);

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_MEMORY_H_ */


#ifndef R_EXT_BOOLEAN_H_
#define R_EXT_BOOLEAN_H_

#undef FALSE
#undef TRUE

#ifdef  __cplusplus
extern "C" {
#endif
  typedef enum { FALSE = 0, TRUE /*, MAYBE */ } Rboolean;

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_BOOLEAN_H_ */

#ifndef R_EXT_CONSTANTS_H_
#define R_EXT_CONSTANTS_H_

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375
#endif

#define PI             M_PI
#define SINGLE_EPS     FLT_EPSILON
#define SINGLE_BASE    FLT_RADIX
#define SINGLE_XMIN    FLT_MIN
#define SINGLE_XMAX    FLT_MAX
#define DOUBLE_DIGITS  DBL_MANT_DIG
#define DOUBLE_EPS     DBL_EPSILON
#define DOUBLE_XMAX    DBL_MAX
#define DOUBLE_XMIN    DBL_MIN

#endif

// Helper table_t class

class Table { 
 
 public:

  // expect column-major specification

  Table() { resize(2,2); }
  
  Table( const int i , const int j , std::vector<double> & nd )
    {
      if ( nd.size() != i*j ) Helper::halt("internal error specifying Table() dimensions");
      resize(i,j);
      d = nd;
    }

  // Special case for 2x2
  Table( const double a, const double b , const double c , const double e ) 
    { 
      resize(2,2); 
      d[0] = a; d[1] = b; d[2] = c; d[3] = e;
    }

  Table( const int a, const int b , const int c , const int e ) 
    {
      resize(2,2);
      d[0] = a; d[1] = b; d[2] = c; d[3] = e;
    }
  
  void resize(const int i, const int j)
  {
    d.clear();
    r = i; c = j;
    d.resize( i * j , 0 );
  }
  
  int nrow() const { return r; }
  int ncol() const { return c; }
  double elem(const int i, const int j) 
  { 
    if ( i < 0 || i >= r ) return false; 
    if ( j < 0 || j >= c ) return false;
    return d[ i + j*r ] ;
  }
  
  double * data() const { return (double*)&(d[0]); }

  bool odds_ratio( double * odds ) const;

 private:
  
  int r, c;

  std::vector<double> d;
    
};

// Fisher's exact test

bool fexact(int *nrow, 
	    int *ncol, 
	    double *table, 
	    int *ldtabl,
	    double *expect, 
	    double *percnt, 
	    double *emin, 
	    double *prt,
	    double *pre, 
	    int *workspace);

bool fisher( const Table &  , double * );

#endif
