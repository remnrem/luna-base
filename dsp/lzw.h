
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

#ifndef __LZW_H__
#define __LZW_H__

#include <string>
#include <vector>
#include <iterator>
struct coarse_t 
{

  coarse_t( const std::vector< std::vector<double> > & , const int , const int );

  std::string epoch(const int) const;
  int size() const { return recoded.size(); } 

private:

  std::vector< std::string > recoded;

};

struct lzw_t 
{ 
  
  //  lzw_t( const std::string & x , bool compress = true );
  lzw_t( const coarse_t & x );
  
  std::string decompress() { return decompress(compressed.begin(), compressed.end()); }
  int size() const { return sizes.size(); } 
  int size(const int e ) const { return sizes[e]; } 

private:
  
  template <typename Iterator> Iterator compress(const std::string &uncompressed, Iterator result);
  template <typename Iterator> std::string decompress(Iterator begin, Iterator end);
  
  std::vector<int> sizes;
  std::vector<int> compressed;  


};

#endif
