
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


#include "lzw.h"

#include <map>
#include <iostream>

#include "helper/helper.h"
#include "miscmath/miscmath.h"
#include "stats/statistics.h"

coarse_t::coarse_t( const std::vector< std::vector<double> > & d , const int nbins , const int nsmooth )  
{
  
  std::vector< std::vector<double> > c1;

  // smooth time-series first?
  if ( nsmooth > 1 )
    {
      for (int e=0;e<d.size();e++)
	{
	  std::vector<double> c2;
	  for (int i=0; i<d[e].size();i += nsmooth)
	    {
	      double value = d[e][i];
	      int mx = i+nsmooth-1;
	      int mxi = nsmooth;
	      if ( mx >= d[e].size() ) 
		{
		  int dif = mx - d[e].size() + 1 ;
		  mxi -= dif;	      
		}
	      for (int j=1;j<mxi;j++) value += d[e][i+j];
	      value /= (double) mxi;	    
	      
	      //if ( mxi <= 0 ) std::cout << "mxi = " << mxi << "\n";
	      c2.push_back( value );
	    }
	  c1.push_back(c2);
	}
    }
  
  // otherwise, just grab the whole thing
  if ( nsmooth == 1 ) c1 = d;
  
  // normalize each epoch 
  for (int e=0;e<c1.size();e++)
    c1[e] = MiscMath::Z( c1[e] );

  // get bins based on standard normal
  if ( nbins < 2 || nbins > 100 ) Helper::halt( "bad nbins" );
  
  double inc = 1.0 / (double) nbins; 
  std::vector<double> t; // size = nbins + 1;
  t.push_back( -99999 );
  double p = inc;
  for (int i = 0; i < nbins-1 ; i++ )
    {
      t.push_back( Statistics::ltqnorm( p ) );
      p += inc;      
    }
  t.push_back( 99999 );
  

  // coarse-ify    
  recoded.resize( c1.size() );
  
  for (int e=0;e<c1.size();e++)
    {
      int last = -1;
      
      recoded[e] = std::string( c1[e].size() , ' ' );
      
      for (int i=0;i<c1[e].size();i++)
	{

	  int curr = -1;
	  
	  // try previous point
	  
	  if ( last != -1 )
	    {
	      if ( c1[e][i] > t[last] && c1[e][i] <= t[ last+1 ] ) 
		curr = last;	    
	    }

	  if ( curr == -1 )
	    {
	      for (int c = 1; c <= nbins ; c++ )
		{
		  //		  std::cout << "testing " << c1[e][i] << "\t" << t[c-1] << "\t" << t[c] << " " << nbins << "\n";
		  if ( c1[e][i] > t[c-1] && c1[e][i] <= t[c] ) 
		    {
		      curr = c;
		      break;		  
		    }
		}
	    }	  

	  //std::cout << "curr = " << curr << "\n"; 

	  if ( curr == -1 ) Helper::halt( "problem in LZW...." );
	  
	  recoded[e][i] = (char)( curr + 32 ) ;
	  
	  // set last point
	  last = curr;
	  
	}

      //std::cout << "E " << e << " [" << recoded[e] << "]\n"; 
    }

}
	    
std::string coarse_t::epoch(const int e) const
{
  if ( e < 0 || e > recoded.size() ) return "";
  return recoded[e];
}


// Compress a string to a list of output symbols.
// The result will be written to the output iterator
// starting at "result"; the final iterator is returned.

template <typename Iterator> Iterator lzw_t::compress(const std::string &uncompressed, Iterator result) 
{
  // Build the dictionary.
  int dictSize = 256;
  std::map<std::string,int> dictionary;
  for (int i = 0; i < 256; i++)
    dictionary[std::string(1, i)] = i;
  
  std::string w;
  for (std::string::const_iterator it = uncompressed.begin();
       it != uncompressed.end(); ++it) {
    char c = *it;
    std::string wc = w + c;
    if (dictionary.count(wc))
      w = wc;
    else {
      *result++ = dictionary[w];
      // Add wc to the dictionary.
      dictionary[wc] = dictSize++;
      w = std::string(1, c);
    }
  }
  
  // Output the code for w.
  if (!w.empty())
    *result++ = dictionary[w];
  return result;
}
 
// Decompress a list of output ks to a string.
// "begin" and "end" must form a valid range of ints
template <typename Iterator> std::string lzw_t::decompress(Iterator begin, Iterator end) 
{
  // Build the dictionary.
  int dictSize = 256;
  std::map<int,std::string> dictionary;
  for (int i = 0; i < 256; i++)
    dictionary[i] = std::string(1, i);
  
  std::string w(1, *begin++);
  std::string result = w;
  std::string entry;
  for ( ; begin != end; begin++) {
    int k = *begin;
    if (dictionary.count(k))
      entry = dictionary[k];
    else if (k == dictSize)
      entry = w + w[0];
    else
      throw "Bad compressed k";
 
    result += entry;
 
    // Add w+entry[0] to the dictionary.
    dictionary[dictSize++] = w + entry[0];
 
    w = entry;
  }
  return result;
}
 

lzw_t::lzw_t( const coarse_t & x )
{
  sizes.clear();
  const int ne = x.size();
  for (int e = 0; e < ne ; e++ ) 
    {      
      std::vector<int> compressed;
      compress( x.epoch(e) , std::back_inserter(compressed));
      sizes.push_back( compressed.size() );
    }
}




