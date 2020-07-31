
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

#include "matrix.h"
#include "statistics.h"
#include <sstream>
#include "helper/helper.h"


template<class T> void Data::Matrix<T>::cbind( const Data::Matrix<T> & rhs )
{
  if ( nrow != rhs.dim1() ) 
    Helper::halt( "cbind() for matrices with unequal number of rows" );
  for (int c=0; c<rhs.dim2(); c++)
    add_col( rhs.col(c) );
}
    
template<class T> void Data::Matrix<T>::add_row( const Vector<T> & r ) 
{ 
  if ( r.size() != ncol ) 
    { 
      if ( nrow == 0 ) { ncol = r.size(); resize(0,r.size()); }
      else { Helper::warn("bad row addition"); return; }
    }
  
  for( int i=0; i<ncol; i++ ) data[i].push_back( r[i] );
  ++nrow;
}

template<class T> void Data::Matrix<T>::add_row( const std::vector<T> & r ) 
{ 
  if ( r.size() != ncol ) 
    {
      if ( nrow == 0 ) { ncol = r.size(); resize(0,r.size()); }
      else { Helper::warn("bad row addition"); return; }
    }
  
  for( int i=0; i<ncol; i++ ) data[i].push_back( r[i] );
  ++nrow;
}

template<class T> Data::Vector<T> Data::Vector<T>::operator*( const Data::Matrix<T> & rhs ) const
{
  return Statistics::matrix_multiply( *this , rhs );
}

template<class T> Data::Vector<T> Data::Vector<T>::operator-( const Data::Vector<T> & rhs ) const
{
  Data::Vector<T> r( rhs.dim1() );
  for (int i=0; i<rhs.dim1(); i++) r[i] = (*this)[i] - rhs[i];
  return r;
}

template<class T> Data::Vector<T> Data::Vector<T>::operator+( const Data::Vector<T> & rhs ) const
{
  Data::Vector<T> r( rhs.dim1() );
  for (int i=0; i<rhs.dim1(); i++) r[i] = (*this)[i] + rhs[i];
  return r;
}

template<class T> void Data::Vector<T>::inplace_add( const double x )
{
  for (int i=0; i<dim1(); i++) data[i] += x;
}

template<class T> void Data::Vector<T>::inplace_multiply( const double x )
{
  for (int i=0; i<dim1(); i++) data[i] *= x;
}

template<class T> void Data::Matrix<T>::inplace_add( const double x )
{
  for (int i=0; i<dim1(); i++) 
    for (int j=0; j<dim2(); j++) 
      (*this)(i,j) += x;
}

template<class T> void Data::Matrix<T>::inplace_multiply( const double x )
{
  for (int i=0; i<dim1(); i++) 
    for (int j=0; j<dim2(); j++) 
      (*this)(i,j) *= x;
}

template<class T> Data::Matrix<T> Data::Matrix<T>::operator*( const Data::Matrix<T> & rhs ) const
{
  return Statistics::matrix_multiply( *this , rhs );
}

template<class T> Data::Vector<T> Data::Matrix<T>::operator*( const Data::Vector<T> & rhs ) const
{
  return Statistics::matrix_multiply( *this , rhs );
}

template<class T> Data::Matrix<T> Data::Matrix<T>::operator-( const Data::Matrix<T> & rhs ) const
{
  Data::Matrix<T> r( rhs.dim1() , rhs.dim2() );
  for (int i=0; i<rhs.dim1(); i++) 
    for (int j=0; j<rhs.dim2(); j++) 
      r(i,j) = (*this)(i,j) - rhs(i,j);
  return r;
}

template<class T> Data::Matrix<T> Data::Matrix<T>::operator+( const Data::Matrix<T> & rhs ) const
{
  Data::Matrix<T> r( rhs.dim1() , rhs.dim2() );
  for (int i=0; i<rhs.dim1(); i++) 
    for (int j=0; j<rhs.dim2(); j++) 
      r(i,j) = (*this)(i,j) + rhs(i,j);
  return r;
}


// pretty-printers
template<class T> std::string Data::Vector<T>::print( const std::string & label , const int nelem ) const
{
  int aelem =  nelem == 0 || nelem  > size() ? size() : nelem ;

  std::stringstream ss;
  if ( label != "" ) ss << label << "\n";
  for (int r=0;r<aelem;r++)
    {
      ss << " [ " << data[r] << " ]\n";
    }
  return ss.str();
}


template<class T> std::string Data::Matrix<T>::print( const std::string & label , const int nrow , const int ncol) const
{
  int arow =  nrow == 0 || nrow > dim1() ? dim1() : nrow ; 
  int acol =  ncol == 0 || ncol > dim2() ? dim2() : ncol ;

  std::stringstream ss;
  if ( label != "" ) ss << label << "\n";

  for (int r=0;r<arow;r++)
    {
      ss << " [ " ;
      for (int c=0;c<acol;c++)
	ss << " " << (*this)(r,c) ;
      ss << " ]\n";
    }
  return ss.str();
}


template<class T> std::string Data::Matrix<T>::dump() const
{
  int arow =  dim1() ;
  int acol =  dim2() ;

  std::stringstream ss;

  for (int r=0;r<arow;r++)
    {
      for (int c=0;c<acol;c++)
	ss << (c ? "\t" : "" ) << (*this)(r,c) ;
      ss << "\n";
    }
  return ss.str();
}


// added to avoid linker errors, given we've defined these templated functions in the .cpp file.
template class Data::Matrix<double>;
template class Data::Vector<double>;
 
