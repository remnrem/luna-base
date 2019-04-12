#
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


#include "topo.h"
#include "helper/helper.h"
#include "dsp/interpolate.h"

#include <fstream>
#include <locale>
#include <string>

bool topo_t::scaled_xy( const std::string & ch , double * x, double * y)
{
  std::map<std::string,int>::const_iterator nn = lab2n.find( ch );
  if ( nn == lab2n.end() ) return false;
  int n = nn->second;
  chid_t chid( n );
  std::map<chid_t,topoloc_t>::const_iterator ii = cxy.find( chid );
  if ( ii == cxy.end() ) return false;

  *x = ii->second.x;
  *y = ii->second.y;
  return true;
}


int topo_t::load( const std::string & filename )
{

  if ( ! Helper::fileExists( filename ) ) 
    Helper::halt( "could not find " + filename );

  cxy.clear();
  lab2n.clear();

  // assume LABEL theta radius(deg)
  std::ifstream IN1( filename.c_str() , std::ios::in );
  int p = 0;
  while ( !IN1.eof() )
    {
      chid_t ch;;
      ch.n = p;

      topoloc_t l;
      IN1 >> ch.label >> l.th >> l.r ;
      
      if ( IN1.eof() ) break;

      // skip
      if ( l.r > 0.6 ) continue;
      
      //      ch.label = Helper::toupper( ch.label );
      lab2n[  ch.label ] = p;
      ++p;

      // degress->radians
      l.th = M_PI/180.0 * l.th;
      // cartesian co-ordinates
      l.x = l.r * cos( l.th );
      l.y = l.r * sin( l.th );
      cxy[ ch ] = l;
      
      logger << "  channel location " << ch.label << " [TH,R] = " << l.th << " " << l.r << ", [X,Y] = " << l.x << " " << l.y << "\n";

    }
  IN1.close();

  
  logger << " read " << cxy.size() << " channel locations\n"; 

  return cxy.size();

}

bool topo_t::add( const std::string & label , const topoloc_t & loc )
{
  // check doesn't already exist
  if ( lab2n.find( label ) != lab2n.end() ) return false;
  chid_t ch;
  ch.label = label;
  ch.n = cxy.size();
  lab2n[ label ] = ch.n;
  cxy[ ch ] = loc; 
  return true;
}


void topo_t::max_radius( double f )
{
  // set 0.5 as actual radius
  squeeze( 0.5 / f );
  pos();
}

void topo_t::squeeze( double f )
{
  
  std::map<chid_t,topoloc_t>::iterator ii = cxy.begin();
  while ( ii != cxy.end() )
    {
      topoloc_t t = ii->second;
      t.x *= f;
      t.y *= f;
      t.r *= f;
      ++ii;
    }  
}

void topo_t::pos()
{
  // make all coords positive
  double xmin = 99, xmax = -99;
  double ymin = 99, ymax = -99;
  std::map<chid_t,topoloc_t>::iterator ii = cxy.begin();
  while ( ii != cxy.end() )
    {
      if ( ii->second.x < xmin ) xmin = ii->second.x;
      if ( ii->second.x > xmax ) xmax = ii->second.x;

      if ( ii->second.y < ymin ) ymin = ii->second.y;
      if ( ii->second.y > ymax ) ymax = ii->second.y;

      ++ii;
    }

  ii = cxy.begin();
  while ( ii != cxy.end() )
    {
      ii->second.x = ( ii->second.x - xmin ) / ( xmax - xmin );
      ii->second.y = ( ii->second.y - ymin ) / ( ymax - ymin );
      ++ii;
    }
  
}

void topo_t::grid( int nx , int ny ) 
{
  grid(0,1,nx,0,1,ny);
}


void topo_t::grid( double xmin, double xmax, int _nx, 
		   double ymin, double ymax, int _ny)
{
  // make a grid: in future, may want to place some circle 
  // constraints on what we bother to estimate
  nx = _nx; ny = _ny;
  out_xy.clear();
  out_inc.clear();
  double xinc = ( xmax - xmin ) / (double)(nx-1);
  double yinc = ( ymax - ymin ) / (double)(ny-1);

  double xx = xmin, yy = ymin;
  for (int xi=0;xi<nx;xi++)
    {
      double xx = xmin + xi * xinc;
      double xx2 = xx - 0.5;

      for (int yi=0;yi<ny;yi++)
	{ 
	  
	  double yy = ymin + yi * yinc;
	  double yy2 = yy - 0.5;
	  
	  double r = sqrt( xx2*xx2 + yy2*yy2 );
	  //double th = atan2( yy , xx );
	  
	  if ( r < 0.5 ) 
	    {
	      out_xy.push_back( xx ); 
	      out_xy.push_back( yy ); 
	      out_inc.push_back( true );
	    }
	  else 
	    {
	      out_inc.push_back( false );
	    }

	}
    }
  
  out_n = out_xy.size() / 2;
  
}

int topo_t::label2n( const std::string & s ) 
{
  std::map<std::string,int>::const_iterator ii = lab2n.find( s );
  if ( ii == lab2n.end() ) return -1;
  return ii->second;
}


Data::Matrix<double> topo_t::interpolate( const std::map<std::string, double> & data ) 
{
    
  if ( out_n == 0 ) 
    Helper::halt( "need to set topo_t::grid() prior to interpolate()" );

  // count number of channels in data[] for which we have topographic information
  // we previously expect a 'signals' to have been called, so we expect that all 
  // channels in topo_t must also be present in data, otherwise throw an error
  // however, there can be channels w/out topographically information -- i.e. in data[]
  // but not in topo_t;  here, we just skip those

  // input 
  inp_n = 0;
  inp_xy.clear();
  std::vector<double> inp_z;
  
  std::map<std::string,double>::const_iterator ii = data.begin();
  while ( ii != data.end() )
    {
      
      //std::string lupper = Helper::toupper( ii->first );
      
      int l = label2n(  ii->first );
      //      std::cout << "<p>found " << l << " for " <<   ii->first  << "</p>";

      if ( l == -1 ) 
	{ 
	  logger << " no topographical information for " << ii->first << " found, dropping\n"; 
	  ++ii; 
	  continue; 
	}
      
      std::map<chid_t,topoloc_t>::const_iterator jj = cxy.find( l );
      if ( jj != cxy.end() ) 
	{	      
	  ++inp_n;
	  inp_xy.push_back( jj->second.x );
	  inp_xy.push_back( jj->second.y );
	  inp_z.push_back( ii->second );
	}
      
      ++ii;
    }

  //  std::cout << "<p>" << inp_n << " is inp_n </p>";


  //
  // show input
  //

  if ( 0 ) 
    {
      int cc = 0;
      for (int i=0;i<inp_z.size();i++) 
	{
	  chid_t ch = cxy.find( i )->first;
	  std::cout << "z " << i << " " << ch.label << " " << ch.n << " "  << inp_z[i] << "\t" 
		    << inp_xy[cc++];
	  std::cout << " " << inp_xy[cc++] << "\n";
	}
    }

  if ( inp_n < 8 ) Helper::halt( "requires at least 8 channels with x-y coordinate information for topographical plots" );
  
  //
  // perform actual 2D interpolation
  //


  dsptools::interpolate2D( this  , inp_z );


  //
  // assemble outpuit
  //
  
  if ( out_inc.size() != nx * ny ) 
    Helper::halt( "internal problem in grid structure");
  
  // interpolate2D() will have sized and populated out_z[] 
  // with interpolated values as specified by given out_xy
  // (i.e. which was set from grid() ).  but not all grid points 
  // were interpolated: out_inc[] controls this

  int p = 0 , k = 0;
  Data::Matrix<double> R(nx,ny,-999);
   for (int i=0;i<nx;i++)
     for (int j=0;j<ny;j++)
       if ( out_inc[p++] ) 
	 R(i,j) = out_z[k++];
   
  return R;
}


void topo_t::dump()
{
  std::map<chid_t,topoloc_t>::const_iterator ii = cxy.begin();
  while ( ii != cxy.end() )
    {
      std::cout << ii->first.label << "\t" << ii->second.x << "\t" << ii->second.y << "\n";
      ++ii;
    }
			
  std::cout << "\ngrid\n";
  int p = 0;
  for (int i=0;i<out_xy.size();i+=2)
    {
      std::cout << "out_xy[" << i << "]\t" << out_xy[i] << "\t" << out_xy[i+1] << "\n";
      ++p;
    }
  
}

