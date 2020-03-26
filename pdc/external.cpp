
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


#include "pdc.h"

#include <string>
#include <iostream>
#include <cmath>
#include <set>

#include "eval.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "miscmath/miscmath.h"
#include "dsp/resample.h"

#include "db/db.h"
#include "helper/logger.h"

extern writer_t writer;

extern logger_t logger;

void pdc_t::external( param_t & param )
{

  const std::string input = param.requires( "input" );
  const std::string output = param.requires( "output" );

  // read data from file 'input'
  // calculate distance matrix, and write to 'output'
  
  // input format expected:
  // ID E VAL1 VAL2...

  if ( ! Helper::fileExists( input ) ) Helper::halt( "could not find file " + input );

  std::map<std::string,std::vector<std::vector<double> > > data;
	     
  std::ifstream IN( input.c_str() , std::ios::in );

  int nv = 0;
  
  while ( ! IN.eof() )
    {
      std::string line;
      Helper::safe_getline( IN , line );
      if ( IN.eof() ) break;
      std::vector<std::string> tok = Helper::parse( line );
      if ( tok.size() <= 2 ) Helper::halt( "bad line: " + line );

      if ( nv != 0 )
	{
	  if ( nv !=  tok.size() - 2 ) Helper::halt( "uneven number of items" );
	}
      else nv = tok.size() - 2;
      
      if ( data.find( tok[0] ) == data.end() )
	data[ tok[0] ].resize(nv);
      
      for (int c=0; c<nv; c++)
	{
	  double val;
	  Helper::str2dbl( tok[c+2] , &val );
	  data[ tok[0] ][c].push_back( val );
	}
    }

  IN.close();

  logger << "  read time-series for " << data.size() << " individuals\n";
    
  // for now, only expected a single channel

  for (int j=0;j<nv; j++)
    add_channel( "_T" + Helper::int2str( j ) );
  
  //
  // 'm' and 't' values that will be used
  //

  const int encoding_m = param.requires_int( "m" );
  const int encoding_t = param.requires_int( "t" );

  set_param( encoding_m , encoding_t );
  
  //
  // Add signals to pdc_t 
  //

  // Always map against _T test-channel label (i.e. all pdlib channels)
  
  std::map<std::string,std::vector<std::vector<double> > >::const_iterator ii = data.begin();
  while ( ii != data.end() )
    {
      
      // a record
      pdc_obs_t rec( nv );
      
      // record
      for (int c=0; c<nv; c++) {
	rec.ch[c] = true;
	rec.ts[c] = ii->second[c];
      }

      rec.encode( encoding_m , encoding_t );

      // add to the main object
      add( rec );
      
      ++ii;
    }

  encode_ts();
	
  //
  // Generate all-by-all matric
  //

  Data::Matrix<double> D = all_by_all();

  //
  // Write
  //
  
  std::ofstream OUT( output.c_str() , std::ios::out );
  OUT << D.dump();
  OUT.close();
  
}
