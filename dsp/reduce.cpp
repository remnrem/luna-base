
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


#include "reduce.h"
#include "miscmath/miscmath.h"
#include "helper/helper.h"

#include <iostream>
#include <cmath>

#define EPS_reduce 1e-4

reduce_t::reduce_t( const std::vector<double> & x , int np )
{
  reduced = false;
  
  // if xs > np, reduce data to np
  const int xs = x.size();  
  if ( xs <= np ) return;
  
  // reduce down to np
  reduced = true;
  hi.resize( np );
  lo.resize( np );
  mean.resize( np );
  sd.resize( np );  
  
  const double t    = xs/(double)np;
  const int    tu   = ceil(t);
  const int    tl   = floor(t);
  const double td   = t - (double)tl;
  
  int m = 0;
  double last = 0;
  
  for (int i=0; i<xs; i+= tu)
    {
      if ( i + tl > xs ) break;
      int ii = i;
      std::vector<double> w;
      double y = 1.0 - last;
      w.push_back( y );
      bool leftovers = true;
      while ( 1 )
	{
	  const double remaining = t - y ; 
	  if ( remaining >= 1 ) 
	    { y++; w.push_back(1); } 
	  else
	    {
	      if ( remaining > EPS_reduce ) 
		{
		  w.push_back( remaining );
		  last = remaining;
		  if ( 1.0 - last < EPS_reduce ) leftovers = false;
		}
	      else
		{
		  last = 0;
		  leftovers = false;
		}
	      break;
	    }	  
	}
      
      // do we need to adjust the count?
      if ( w.size() == tu && leftovers && tl != tu ) --i;
      
      // get (Weighted) mean/ SD and min/max for this bin
      double s = 0;
      double v = 0, vv=0;
      double xmin = x[ii];
      double xmax = x[ii];
      
      for (int j=0;j<w.size();j++)
	{
	  s += w[j] * x[ii+j];
	  v += w[j];
	  vv += w[j]*w[j];
	  if      ( x[ii+j] < xmin ) xmin = x[ii+j];
	  else if ( x[ii+j] > xmax ) xmax = x[ii+j];	  
	}

      double wm = s/v;
      
      double ssq = 0;      
      for (int j=0;j<w.size();j++)
	ssq += w[j] * ( x[ii+j] - wm )  * ( x[ii+j] - wm );
      double wsd = sqrt( ssq / ( v - ( vv/v ) ) );
      
      hi[m] = xmax;
      lo[m] = xmin;
      mean[m] = wm;
      sd[m] = wsd;
      
      //       std::cout << "bin " << m << " " 
      // 		<< wm << " " << wsd << " " << xmin << "  " << xmax << "\n";

      // next bin
      ++m;
            
      // all done?
      if ( m >= np ) break;
      

    }
  
}




reduce_t::reduce_t( const std::vector<double> * x , 
		    const std::vector<uint64_t> * t , uint64_t a , uint64_t b , int np )
{
  // span is from a..b
  // we need to divide 'x' into 'np' segments
  // some of these may have no data (i.e. EDF+D)
  
  uint64_t span = b - a + 1;
  uint64_t each = span / (uint64_t)np;
  //  std::cerr << "each bin spans " << each << " tps\n";
  
  const int nx = x->size();
  
  mean.resize(np,0);
  sd.resize(np,0);
  n.resize(np,0);
  lo.resize(np,0);
  hi.resize(np,0);
  
  std::vector<double> sx(np,0);
  std::vector<double> sxx(np,0);

  //  std::cout << "aiming " << nx << " to " << np << "\n";


  uint64_t next = a + each;
  int p = 0;
  bool first = true;
  for (int i=0;i<nx;i++)
    {

      // moving into the next bin?
      // may be some rounding issues -- okay to dump into the last bin for now..
      if ( (*t)[i] >= next ) 
	{ 
	  first = true; 
	  ++p; 
	  next += each;
	  if ( p >  np ) { std::cerr << "prblemo!\n"; std::exit(1); } 
	  if ( p == np ) { p--; }
	  
	}
      
//       std::cout << "comp: " << i << " " << p << "/" << np << "\t"
// 		<< (*t)[i] << "\t"
// 		<< next << "\n";
      


      n[p]++;
      sx[p] += (*x)[i];
      sxx[p] += (*x)[i] * (*x)[i];

      if ( first ) 
	{
	  lo[p] = (*x)[i];
	  hi[p] = (*x)[i];
	  first = false;
	}
      else
	{
	  if ( (*x)[i] < lo[p] ) lo[p] = (*x)[i];
	  if ( (*x)[i] > hi[p] ) hi[p] = (*x)[i];
	}
    }
  
  for (int p=0;p<np;p++)
    {
      if ( n[p] > 0 )
	{
	  mean[p] = sx[p] / (double)n[p];
	  if ( n[p] > 2 )
	    sd[p] = sqrt( ( sxx[p] - (sx[p]*sx[p])/(double)n[p] ) / ((double)n[p]-1.0) );      
	  else
	    sd[p] = 0;
	}
      
//       std::cout << "p = " << p << "\t"
// 		<< n[p] << "\t"
// 		<< mean[p] << "\t"
// 		<< sd[p] << "\t"
// 		<< lo[p] << "\t"
// 		<< hi[p] << "\n";
    }

}
