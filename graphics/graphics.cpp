
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


#include "graphics.h"

#ifndef NO_HPDFLIB

#include <cstdlib>
#include <cstdio>
#include <setjmp.h>
#include "helper/helper.h"

// Fonts

//Courier
//Courier-Bold
//Courier-Oblique
//Courier-BoldOblique
//Helvetica
//Helvetica-Bold
//Helvetica-Oblique
//Helvetica-BoldOblique
//Times-Roman
//Times-Bold
//Times-Italic
//Times-BoldItalic

// In the default coordinate system of PDF, shown below, the
// lower-left corner is at coordinates (0, 0), and the upper-right
// corner is at coordinates (width, height). The default resolution is
// 72dpi.

jmp_buf env;

#ifdef HPDF_DLL
void  __stdcall
#else
void
#endif
error_handler  (HPDF_STATUS   error_no,
                HPDF_STATUS   detail_no,
                void         *user_data)
{
  printf ("ERROR: error_no=%04X, detail_no=%u\n", (HPDF_UINT)error_no,
	  (HPDF_UINT)detail_no);
  longjmp(env, 1);
}



pdf_t::pdf_t()
{

  pdf = HPDF_New(error_handler, NULL);
  if ( ! pdf ) { Helper::halt( "cannot open PDF" ); }
  if (setjmp(env)) { HPDF_Free(pdf); Helper::halt( "problem in pdf_t(), bailing..."); }

  /* set compression mode */
  HPDF_SetCompressionMode (pdf, HPDF_COMP_ALL);

  // /* set page mode to use outlines. */
  // HPDF_SetPageMode (pdf, HPDF_PAGE_MODE_USE_OUTLINE);
  
  // /* set password */
  // HPDF_SetPassword (pdf, "owner", "user");

  // palette

  palette[ "white" ] = rgb_t(1,1,1);
  palette[ "black" ] = rgb_t(0,0,0);
  palette[ "gray" ] = rgb_t( 0.5, 0.5, 0.5);
  palette[ "silver" ] = rgb_t( 0.75, 0.75, 0.75);
  palette[ "maroon" ] = rgb_t( .5, 0, 0);
  palette[ "red" ] = rgb_t( 1,0 ,0 );
  palette[ "olive" ] = rgb_t( .5, .5, 0);
  palette[ "yellow" ] = rgb_t( 1, 0, 1);
  palette[ "green" ] = rgb_t( 0, .5, 0);
  palette[ "lime" ] = rgb_t( 0, 1, 0);
  palette[ "teal" ] = rgb_t( 0, .5, .5);
  palette[ "aqua" ] = rgb_t( 0, 1, 1);
  palette[ "navy" ] = rgb_t( 0, 0, .5);
  palette[ "blue" ] = rgb_t( 0, 0, 1);
  palette[ "purple" ] = rgb_t( .5, 0, .5);
  palette[ "fuschia" ] = rgb_t( 1, 0, 1 );  

}



pdf_t::~pdf_t()
{
  HPDF_Free(pdf);
}


int pdf_t::add_page( int nx , int ny )
{
  
  HPDF_Page next_page;
  next_page = HPDF_AddPage(pdf);
  HPDF_Page_SetSize ( next_page,  HPDF_PAGE_SIZE_LETTER , HPDF_PAGE_LANDSCAPE);
  
  // record height, etc
  page = next_page;
  height = HPDF_Page_GetHeight (page);
  width = HPDF_Page_GetWidth (page);

  // Font

  font = HPDF_GetFont (pdf, "Helvetica", NULL);
  fontsize = 8;
  set_fontsize( 8 );

  // border (proportionally)
  xborder = 0.025;
  yborder = 0.025;
  
  // plot grid
  grid_nx = nx;
  grid_ny = ny;
  
  grid_spacer = 0.02;
  
  grid_spacerx = grid_spacer * width;
  grid_spacery = grid_spacer * height;
  
  grid_sx = ( width  - ( grid_nx-1 ) * grid_spacerx - width*2*xborder) / (double)grid_nx;
  grid_sy = ( height - ( grid_ny-1 ) * grid_spacery - height*2*yborder) / (double)grid_ny;
  
  grid_currx1 = grid_currx2 = 0;
  grid_curry1 = grid_curry2 = 0;
  
  set_grid(0,0);
  
  
  return true;
}

bool pdf_t::set_grid( int n )
{
  const int t = grid_nx * grid_ny;
  if ( n < 0 || n >= t ) return false;  
  int yy = n / grid_nx; 
  int xx = n - yy * grid_nx;
  set_grid( xx , yy );
  return true;
}

void pdf_t::set_grid( int a , int b , int a2 , int b2 )
{
  if ( a2 == -1 ) a2 = a;
  if ( b2 == -1 ) b2 = b;
  if ( a < 0 || a >= grid_nx ) Helper::halt("bad grid setting");
  if ( b < 0 || b >= grid_ny ) Helper::halt("bad grid setting");
    
  grid_currx1 = a; 
  grid_currx2 = a2;
  
  grid_curry1 = b;
  grid_curry2 = b2;

  // upper left corner of (a,b) 
  // lower right corner of (a2,b2)

  x1 = width * xborder + a * grid_sx + a * grid_spacerx ;
  x2 = width * xborder + a2 * grid_sx + a2 * grid_spacerx + grid_sx - 1 ;

  y1 = height * yborder + b * grid_sy + b * grid_spacery ;
  y2 = height * yborder + b2 * grid_sy + b2 * grid_spacery + grid_sy - 1 ;

  if ( x1 < 0 ) x1 = 0;
  if ( y1 < 0 ) y1 = 0;
  if ( x2 >= width ) x2 = width-1;
  if ( y2 >= height ) y2 = height-1;
   
  // flip Y axis... (i.e. so that 0 is top of screen)
  y1 = height - y1;
  y2 = height - y2;
}

float pdf_t::x( float px ) const
{
  if ( px < 0 ) px = 0;
  if ( px > 1 ) px = 1;
  float xw = x2 - x1 + 1;
  return x1 + px * xw;
}

float pdf_t::y( float py ) const
{
  if ( py < 0 ) py = 0;
  if ( py > 1 ) py = 1;
  float yw = y2 - y1 + 1;
  return y1 + py * yw;
}


bool pdf_t::write(const std::string & f )
{
  HPDF_SaveToFile (pdf, f.c_str() );
  return true;
}


void pdf_t::newdoc()
{
  HPDF_NewDoc (pdf);
}


void pdf_t::outline()
{
  set_font( "Helvetica" );
  HPDF_UINT x, y;
  
  HPDF_Page_SetFontAndSize (page, font, 5);
  HPDF_Page_SetGrayFill (page, 0.5);
  HPDF_Page_SetGrayStroke (page, 0.8);

}


void pdf_t::set_line_type_butt()
{
  HPDF_STATUS st = HPDF_Page_SetLineCap ( page ,  HPDF_BUTT_END  );
}

void pdf_t::set_line_type_round()
{
  HPDF_STATUS st = HPDF_Page_SetLineCap ( page ,   HPDF_ROUND_END   );
}

void pdf_t::set_line_type_square()
{
  HPDF_STATUS st = HPDF_Page_SetLineCap ( page ,   HPDF_PROJECTING_SCUARE_END   );
}

void pdf_t::set_line_join_miter()
{
  HPDF_STATUS st = HPDF_Page_SetLineJoin (page , HPDF_MITER_JOIN );
}

void pdf_t::set_line_join_round()
{
  HPDF_STATUS st = HPDF_Page_SetLineJoin (page ,  HPDF_ROUND_JOIN  );
}


void pdf_t::set_line_join_bevel()
{
  HPDF_STATUS st = HPDF_Page_SetLineJoin (page ,  HPDF_BEVEL_JOIN  );
}


void pdf_t::set_font( const std::string & f )
{
  font = HPDF_GetFont (pdf, f.c_str() , NULL);
}

void pdf_t::set_fontsize( int f )
{
  fontsize = f;
  HPDF_STATUS font_change = HPDF_Page_SetFontAndSize (page , font , fontsize );
}

void pdf_t::set_font_color( const std::string & c )
{
  if ( palette.find( c ) != palette.end() ) set_font_color( palette[ c ] );
}

void pdf_t::set_font_color(double r , double g , double b)
{
  HPDF_STATUS st1 = HPDF_Page_BeginText( page );  
  HPDF_STATUS st2a = HPDF_Page_SetRGBStroke( page , r , g, b );  
  HPDF_STATUS st2b = HPDF_Page_SetRGBFill( page , r , g, b );  
  HPDF_STATUS st3 = HPDF_Page_EndText( page );
}

void pdf_t::circle( double lx, double ly , double r )
{
  float px = x( lx );
  float py = y( ly );
  float radius = r * grid_sx;
  HPDF_STATUS st = HPDF_Page_Circle  (page , px , py , radius );  
}

void pdf_t::rectangle( double lx ,double ly , double ux, double uy )
{
  float lwr_x = x(lx);
  float lwr_y = y(ly);
  float upr_x = x(ux);
  float upr_y = y(uy);
  float rect_width = upr_x - lwr_x;
  float rect_height = upr_y - lwr_y;

  //  std::cout << "rect " << lwr_x << ","<< lwr_y << "  to  " << upr_x << "," << upr_y << "\n";

  HPDF_STATUS st = HPDF_Page_Rectangle( page , lwr_x , lwr_y , rect_width , rect_height );
  
}


void pdf_t::textbox( float px1 , float py1 , float px2 , float py2 , const std::string & t , HPDF_TextAlignment align )
{
  
  HPDF_STATUS st1 = HPDF_Page_BeginText( page );  
  
  px1 = x(px1);
  py1 = y(py1);
  px2 = x(px2);
  py2 = y(py2);
  
  HPDF_STATUS st2 = HPDF_Page_TextRect ( page , px1 , py1 , px2 , py2 , 
					 t.c_str() , align , NULL );
  
  HPDF_STATUS st3 = HPDF_Page_EndText( page );

}

void pdf_t::text( double lx , double ly , const std::string & t )
{
  
  float px = x(lx);
  float py = y(ly);
  
  HPDF_STATUS st1 = HPDF_Page_BeginText( page );  
  HPDF_STATUS st2 = HPDF_Page_TextOut( page , px , py , t.c_str() );
  HPDF_STATUS st3 = HPDF_Page_EndText( page );
  
}

void pdf_t::set_line_width( double w )
{
  HPDF_STATUS st = HPDF_Page_SetLineWidth( page , w );
}

void pdf_t::set_line_type_dashed()
{
  const HPDF_UINT16 DASH_MODE1[] = {3};
  HPDF_STATUS st = HPDF_Page_SetDash( page , DASH_MODE1 , 1 , 1 );
}

void pdf_t::set_line_type_solid()
{
  HPDF_STATUS st = HPDF_Page_SetDash( page , NULL , 0 , 0 );
}

void pdf_t::set_grayscale( double g )
{
  HPDF_STATUS st = HPDF_Page_SetGrayStroke( page , g ); // 0..1
}

void pdf_t::set_grayscale_fill( double g )
{
  HPDF_STATUS st = HPDF_Page_SetGrayFill( page , g ); // 0..1
}

void pdf_t::set_color( const std::string & c )
{
  if ( palette.find( c ) != palette.end() ) set_color( palette[ c ] );
}

void pdf_t::set_color_fill( const std::string & c )
{
  if ( palette.find( c ) != palette.end() ) set_color_fill( palette[ c ] );
}

void pdf_t::set_color( const rgb_t & rgb )
{
  HPDF_STATUS st = HPDF_Page_SetRGBStroke( page , rgb.r , rgb.g , rgb.b  ); 
}

void pdf_t::set_color_fill( const rgb_t & rgb )
{
  HPDF_STATUS st = HPDF_Page_SetRGBFill( page , rgb.r , rgb.g , rgb.b  ); 
}


void pdf_t::move( double lx , double ly )
{
  HPDF_STATUS st = HPDF_Page_MoveTo(page, x(lx) , y(ly) );
}

void pdf_t::line( double lx , double ly )
{
  HPDF_STATUS st = HPDF_Page_LineTo( page , x(lx) , y(ly) );
}

void pdf_t::stroke()
{
  HPDF_STATUS st = HPDF_Page_Stroke( page );
}

void pdf_t::fill()
{
  HPDF_STATUS st = HPDF_Page_Fill( page );
}


void pdf_t::stroke_fill()
{
  HPDF_STATUS st = HPDF_Page_FillStroke( page );
}


void pdf_t::heatmap( float px1 , float py1 , float px2 , float py2 , const std::vector<std::vector<double> > & d , const std::vector<double> & yaxis )
{
  // convert coordinates
  float lwr_x = x(px1);
  float lwr_y = y(py1);
  float upr_x = x(px2);
  float upr_y = y(py2);

  float rect_width = upr_x - lwr_x;
  float rect_height = upr_y - lwr_y;

  // draw frame
  set_line_width( 0.2 );
  set_grayscale( 1 );  
  rectangle( px1 , py1 , px2 , py2 ); // i.e. original co-ordinates
  stroke();
  
  // labels
  set_fontsize( 6 );
  text( px1 - 0.02 , py1 , Helper::dbl2str( yaxis[0] ) );
  text( px1 - 0.02 , py2 , Helper::dbl2str( yaxis[yaxis.size()-1] ) );
  

  // d rows / cols
  if ( d.size() == 0 ) return;
  const int nrow = d.size();
  const int ncol = d[0].size();
  
  const double yrange = fabs( py2 - py1 );
  const double row_inc = yrange / (double)nrow;
  const double col_inc = ( px2 - px1 ) / (double)ncol; 

  // canonical markers
  for (int f=0; f<nrow; f++)
    if ( yaxis[f] == 13.5 )  
      {
	set_font_color( "red" );
	text( px1 - 0.02 , py1 - f * row_inc , "*" );
      }

  for (int r=0; r<nrow; r++)
    for (int c=0; c<ncol; c++)
      {
	rgb_t z = rgb_t::heatmap( d[r][c] );
	set_color_fill( z );
	set_color( z );
	rectangle( px1 + c * col_inc , 
		   py1 - r * row_inc , 
		   px1 + (c+1) * col_inc , 
		   py1 - (r+1) * row_inc ) ;
	stroke_fill();
      }
  
//  HPDF_STATUS st = HPDF_Page_Rectangle( page , lwr_x , lwr_y , rect_width , rect_height );
  
}


#endif


