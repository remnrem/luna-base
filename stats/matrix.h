
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

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>
#include <string>

namespace Data { 

  template<class T> class Matrix;

  template<class T = double> class Vector {

    public:

    Vector() { } 
    Vector(int i) { resize(i); mask.resize(i,false); }
    Vector(const std::vector<T> & x ) { data = x; mask.resize( data.size() , false ); }
    
    void clear() { data.clear(); mask.clear(); }
    void resize(const int i) { data.resize(i); mask.resize(i,false); }
    void resize(const int i, const T & t) { data.resize(i,t); mask.resize(i,false); }
    void push_back( const T & t ) { data.push_back(t); mask.push_back(false); }
    int size() const { return data.size(); }
    int dim1() const { return data.size(); }
    
    T  operator[] ( const unsigned int i) const { return data[i]; } 
    T & operator[] ( const unsigned int i) { return data[i]; } 
    
    T operator() (const unsigned int i ) const { return data[i]; }
    T & operator() (const unsigned int i ) { return data[i]; }
    
    Vector<T> operator*( const Matrix<T> & rhs ) const;
    Vector<T> operator+( const Vector<T> & rhs ) const;
    Vector<T> operator-( const Vector<T> & rhs ) const;

    // some convenience functions
    
    void inplace_add( const double x );
    void inplace_multiply( const double x );

    std::string print( const std::string & label = "" , const int nelem = 0 ) const;
    
    void set_elem_mask( const int r , const bool val = true )
    {
      if ( r < 0 || r >= mask.size() ) return;
      mask[r] = val;
    }

    bool masked(const int r ) const
    {
      if ( r < 0 || r >= data.size() ) return false;
      return mask[r];
    }
    
    Vector purge_rows() 
    {
      int sz = 0;
      for (int i=0; i<mask.size(); i++) if ( !mask[i] ) ++sz;
      Vector<T> v( sz );
      sz = 0;
      for (int i=0; i<mask.size(); i++) if ( !mask[i] ) v[sz++] = data[i]; 
      return v;
    }
    
    const std::vector<T> * data_pointer() const { return data.size() ? &data : NULL ; }
    T * elem_pointer( const int i ) { return data.size() ? &data[i] : NULL ; }
    std::vector<T> extract() const { return data; } // ignores mask

    private:
    
    std::vector<T> data;
    std::vector<bool> mask;
  };

  

  template<class T = double> class Matrix {
    
    public:

    // row access

    struct Row
    { 
      Row( Matrix & m , int i ) : mat(m), row(i) { }
      T & operator[](const int j) { return mat(row,j); }     
      //T operator[](const int j) const { return mat(row,j); }     
      private:
      int row;
      Matrix & mat;      
    };
    
    struct ConstRow
    { 
      ConstRow( const Matrix & m , int i ) : mat(m) , row(i) { } 
      T operator[](const int j) const { return mat(row,j); }           
      private:
      int row;
      const Matrix & mat;      
    };
    
    Matrix() { clear(); } 
    Matrix(const int r, const int c) { resize(r,c); }
    Matrix(const int r, const int c, const T & t) { resize(r,c,t); }
    
    T operator() (const unsigned int i, const unsigned int j ) const { return data[j][i]; }
    T & operator() (const unsigned int i, const unsigned int j ) { return data[j][i]; }
    
    Row operator[] ( const unsigned int i) { return Row(*this,i); }
    ConstRow operator[] ( const unsigned int i) const { return ConstRow(*this,i); }
    
    Vector<T> row( const int r ) 
    { 
      Vector<T> d( ncol );
      for (int c=0; c<ncol; c++) d[c] = (*this)(r,c);
      return d;
    } 

    Vector<T> col( const int c ) const { return data[c]; } 
    Vector<T> & col( const int c ) { return data[c]; } 
    const Vector<T> * col_pointer( const int c ) { return &data[c]; } 

    void add_col( const Vector<T> & r ) 
    { 
      if ( ncol == 0 ) nrow = r.size(); 
      data.push_back(r);  // add data
      ++ncol;             // track increase in col count
      
      // propagate case-wise missingness across columns for each row
      for (int i=0; i<r.size(); i++) 
	if( r.masked(i) ) set_row_mask(i); 
    }
    
    void add_col( const std::vector<T> & r ) 
    { 
      if ( ncol == 0 ) nrow = r.size();
      data.push_back( Vector<T>(r) ); ++ncol; 
    }
    
    void cbind( const Data::Matrix<T> & rhs );
    
    void add_row( const Vector<T> & r ) ;

    void add_row( const std::vector<T> & r );
      
    void set_row_mask( int r , const bool b = true ) 
    { 
      if ( r >= 0 && r < nrow ) row_mask[r] = b;
    }
    
    bool masked( const int r ) const
    {
      if ( ncol == 0 ) return false; // no data
      if ( r >= 0 && r < nrow ) return row_mask[r];
      return true; // mask out-of-range items
    }

    void clear() { data.clear(); row_mask.clear(); nrow = ncol = 0; }

    Matrix<T> purge_rows() 
    {
      int sz = 0; 
      for (int i=0; i<row_mask.size(); i++) if ( ! row_mask[i] ) ++sz;
      Matrix<T> v( sz , ncol );
      for (int c = 0 ; c < ncol ; c++ ) 
	{
	  int sz = 0;
	  for (int r=0; r<nrow; r++) if ( ! row_mask[r] ) v( sz++ , c ) = data[c][r]; 
	}      
      return v;
    }
    
    void resize(const int r, const int c) 
    { 
      nrow = r;
      ncol = c;
      row_mask.resize( nrow , false ); // masked-out
      data.resize(c); 
      for (int j=0; j<c; j++) data[j].resize( nrow );
    }
    
    void resize(const int r, const int c, const T & t ) 
    { 
      nrow = r;
      ncol = c;
      row_mask.resize( nrow , false ); // masked-out
      data.resize(c); 
      for (int j=0; j<c; j++) data[j].resize( nrow , t );
    }

    int dim1() const { return nrow; }
    int dim2() const { return ncol; }
    
    // pretty-print
    
    std::string print( const std::string & label = "" , const int nrow = 0 , const int ncol = 0 ) const;

    // some convenience functions
    
    void inplace_add( const double x );
    void inplace_multiply( const double x );

    // op. overloading for common matrix operations
    
    Matrix<T> operator*( const Data::Matrix<T> & rhs ) const;
    Vector<T> operator*( const Data::Vector<T> & rhs ) const;

    Matrix<T> operator+( const Data::Matrix<T> & rhs ) const;
    Matrix<T> operator-( const Data::Matrix<T> & rhs ) const;    

    private:
    
    std::vector< Vector<T> > data;
    std::vector<bool> row_mask;
    int nrow ;
    int ncol ;
  };


  
  

}


#endif
