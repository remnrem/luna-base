
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


#ifndef __LUNA_GRAPHICS_H__
#define __LUNA_GRAPHICS_H__

#ifndef NO_HPDFLIB

#include "hpdf.h"

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

struct rgb_t
{  
  rgb_t() { r=0;g=0;b=0; }
  rgb_t(double r,double g,double b) : r(r) , g(g) , b(b) { } 
  double r,g,b;  
  
  static rgb_t heatmap(double value)
  {
    const int NUM_COLORS = 4;
    
    static float color[NUM_COLORS][3] = { {0,0,1}, {0,1,0}, {1,1,0}, {1,0,0} };
    // A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.
    
    int idx1;        // |-- Our desired color will be between these two indexes in "color".
    int idx2;        // |
    float fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.
    
    if      (value <= 0)  { idx1 = idx2 = 0;            }    // accounts for an input <=0
    else if (value >= 1)  { idx1 = idx2 = NUM_COLORS-1; }    // accounts for an input >=0
    else
      {
	value = value * (NUM_COLORS-1);        // Will multiply value by 3.
	idx1  = floor(value);                  // Our desired color will be after this index.
	idx2  = idx1+1;                        // ... and before this index (inclusive).
	fractBetween = value - float(idx1);    // Distance between the two indexes (0-1).
      }
    
    double red   = (color[idx2][0] - color[idx1][0])*fractBetween + color[idx1][0];
    double green = (color[idx2][1] - color[idx1][1])*fractBetween + color[idx1][1];
    double blue  = (color[idx2][2] - color[idx1][2])*fractBetween + color[idx1][2];
    
    return rgb_t( red , green , blue );
  }

};


struct pdf_t 
{
  
  pdf_t();
  ~pdf_t();
  
  int add_page(int,int);
  bool write(const std::string &);
  void newdoc();

  void outline();

  void set_fontsize(int);

  void set_line_type_round();
  void set_line_type_butt();
  void set_line_type_square();

  void set_line_join_miter();
  void set_line_join_round();
  void set_line_join_bevel();

  void set_font(const std::string &);
  void set_font_color( const std::string & );
  void set_font_color( const rgb_t & rgb ) { set_font_color(rgb.r,rgb.g,rgb.b); }
  void set_font_color(double,double,double); // 0.0 .. 1.0

  void rectangle( double,double,double,double );
  void circle( double,double,double);
  void text( double,double,const std::string & f );
  void textbox( float px1 , float py1 , float px2 , float py2 , const std::string & t , HPDF_TextAlignment align =  HPDF_TALIGN_LEFT );
  
  void heatmap( float px1 , float py1 , float px2 , float py2 , const std::vector<std::vector<double> > & d , const std::vector<double> & );
  
  void set_line_width( double w );
  void set_line_type_dashed();
  void set_line_type_solid();

  void set_grayscale( double g );
  void set_grayscale_fill( double g );
  

  void set_color( const std::string & );
  void set_color_fill( const std::string & );
  void set_color( const rgb_t & );
  void set_color_fill( const rgb_t & );
  
  void move( double x , double y );
  void line( double x , double y );

  void fill();
  void stroke();
  void stroke_fill();
  
  HPDF_Doc pdf;
  HPDF_Page page;

  HPDF_REAL height; 
  HPDF_REAL width;

  HPDF_Font font;
  float fontsize;
  
  float xborder , yborder;

  int grid_nx, grid_ny;
  int grid_currx1 , grid_curry1;
  int grid_currx2 , grid_curry2;
  float grid_spacer , grid_spacerx, grid_spacery;
  float grid_sx, grid_sy;
  
  float x1,x2,y1,y2;

  // effective range (from 0..1, given current grid box(es)

  void set_grid( int a , int b , int a2 = -1 , int b2 = - 1);
  bool set_grid( int n ); // i.e., by rows, then cols, assume a single only
  
  float x( float px ) const;
  float y( float py ) const;

  // colour palette
  std::map<std::string,rgb_t> palette;
  
  
};


#endif

#endif

