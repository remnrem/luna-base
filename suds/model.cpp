
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

#include "suds.h"

#include <vector>
#include <map>
#include <set>
#include <iomanip>

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

#include "stats/eigen_ops.h"
#include "stats/statistics.h"
#include "miscmath/miscmath.h"

extern logger_t logger;

extern writer_t writer;


//
// Specification of SUDS models
//

void suds_model_t::init()
{

  lab2ftr[ "SPEC" ] = SUDS_LOGPSD;
  lab2ftr[ "RSPEC" ] = SUDS_RELPSD;
  lab2ftr[ "VSPEC" ] = SUDS_CVPSD;
  lab2ftr[ "SLOPE" ] = SUDS_SLOPE;
  lab2ftr[ "SKEW" ] = SUDS_SKEW;
  lab2ftr[ "KURTOSIS" ] = SUDS_KURTOSIS;
  lab2ftr[ "HJORTH" ] = SUDS_HJORTH;
  lab2ftr[ "FD" ] = SUDS_FD;
  lab2ftr[ "PE" ] = SUDS_PE;
  lab2ftr[ "MEAN" ] = SUDS_MEAN;
  lab2ftr[ "TIME" ] = SUDS_TIME;
  lab2ftr[ "SMOOTH" ] = SUDS_SMOOTH;
  lab2ftr[ "DENOISE" ] = SUDS_DENOISE;
  lab2ftr[ "SMOOTH2" ] = SUDS_SMOOTH2;
  lab2ftr[ "DENOISE2" ] = SUDS_DENOISE2;

  ftr2lab[ SUDS_LOGPSD] = "SPEC";
  ftr2lab[ SUDS_RELPSD] = "RSPEC";
  ftr2lab[ SUDS_CVPSD] = "VSPEC";
  ftr2lab[ SUDS_SLOPE] = "SLOPE";   
  ftr2lab[ SUDS_SKEW] = "SKEW";
  ftr2lab[ SUDS_KURTOSIS] = "KURTOSIS";
  ftr2lab[ SUDS_HJORTH] = "HJORTH";
  ftr2lab[ SUDS_FD] = "FD";      
  ftr2lab[ SUDS_PE] = "PE";  
  ftr2lab[ SUDS_MEAN] = "MEAN";
  ftr2lab[ SUDS_TIME] = "TIME";    
  ftr2lab[ SUDS_SMOOTH] = "SMOOTH";
  ftr2lab[ SUDS_DENOISE] = "DENOISE";  
  ftr2lab[ SUDS_SMOOTH2] = "SMOOTH2";
  ftr2lab[ SUDS_DENOISE2] = "DENOISE2";  

  // other clears/resets
  nc = 0;  
  chs.clear();
  specs.clear();
  fcmap.clear();
  
}

 
bool suds_model_t::read( const std::string & modelfile ,
			 const std::string & winfile ,
			 const std::string & woutfile ,
			 const std::string & default_channel )
{


  if ( modelfile == "" ) 
    Helper::halt( "error specifying SOAP model file: empty model file name" );
  
  // this file needs to specify:
  //  channels used, and sample rates  (CH)
  //  the number of SVD components that will be extracted (NC)
  //  the features to construct the raw matrix (SPEC, etc)

  suds_t::nc = 0;
  
  // ensure we have initiated the maps
  init();
  

  // clear any current specifications
  specs.clear();

  std::vector<std::string> lines;

  if ( modelfile[0] != '_' ) 
    {
      if ( ! Helper::fileExists( modelfile ) )
	Helper::halt( "could not open " + modelfile );
      
      // read from a file
      std::ifstream IN1( modelfile.c_str() , std::ios::in );
      while ( ! IN1.eof() )
	{
	  std::string line;
	  Helper::safe_getline( IN1 , line );
	  if ( line == "" ) continue;
	  if ( line[0] == '%' ) continue;
	  lines.push_back( line );
	}
      IN1.close();
    }
  else
    {
      // populate 'default' SOAP model

      lines.push_back( "CH " + default_channel + " 128" );	

      if ( modelfile == "_1" )
	{
	  lines.push_back( "SPEC " + default_channel + " lwr=0.5 upr=25" );
	}
      else if ( modelfile == "_2" )
	{
	  lines.push_back( "SPEC " + default_channel + " lwr=0.5 upr=25" );
	  lines.push_back( "RSPEC " + default_channel + " lwr=5 upr=20 z-lwr=30 z-upr=45" );
	  lines.push_back( "SLOPE " + default_channel );
	  lines.push_back( "SKEW " + default_channel + "" );
	  lines.push_back( "KURTOSIS " + default_channel );
	  lines.push_back( "FD " + default_channel );
	  lines.push_back( "PE " + default_channel );
	  lines.push_back( "DENOISE2 lambda=0.5" );
	  //%SMOOTH2 half-window=15
	  lines.push_back( "TIME order=4" );	  
	}
      
      lines.push_back( "NC 10" );
    }  
  
  for (int l=0; l<lines.size(); l++)
    {

      // get line
      const std::string & line = lines[l];
      
      std::vector<std::string> tok = Helper::parse( line , " \t" );

      // expecting format:
      //  CHANNEL CH SR 
      //  FEATURE  { CH } { KEY=VAL }
      //  i.e. either a channel, or (if has an '=') an argument
      
      if ( tok.size() < 2 ) Helper::halt( "bad format for line: " + line );

      // channel specifier?
      if ( Helper::toupper( tok[0] ) == "CH" )
	{
	  if ( tok.size() != 3 )
	    Helper::halt( "expecing: CH label SR" );
	  
	  int sr ;
	  if ( ! Helper::str2int( tok[2] , &sr ) )
	    Helper::halt( "bad format: " + line );
	  
	  // store
	  chs[ tok[1] ] = suds_channel_t( tok[1] , sr ) ;

	  // next line
	  continue;
	}

      // component (NC) specifier?
      if ( Helper::toupper( tok[0] ) == "NC" )
	{
	  if ( tok.size() != 2 )
            Helper::halt( "expecing: NC value" );

	  int nc ;
	  if ( ! Helper::str2int( tok[1] , &nc ) )
            Helper::halt( "bad format: " + line );
	  if ( nc < 1 )
	    Helper::halt( "bad format: " + line );

	  suds_t::nc = nc;
	  continue;
	}

      
      // feature specifier?

      if ( lab2ftr.find( Helper::toupper( tok[0] ) ) == lab2ftr.end() )
	Helper::halt( "feature not recognized: " + tok[0] );

      // get list of channels, args separately
      std::vector<std::string> tchs;
      std::map<std::string,double> targs;

      for (int i=1; i<tok.size(); i++)
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[i] , '=' );
	  
	  if ( tok2.size() > 2 ) Helper::halt( "bad format: " + tok[i] );
	  
	  // add as a channel
	  if ( tok2.size() == 1 )
	    {
	      // has the channel already been specified via CH?
	      if ( tok2[0] != "." && chs.find( tok2[0] ) == chs.end() )
		Helper::halt( tok2[0] + " not specified via 'CH' yet: " + line );	      
	      tchs.push_back( tok2[0] );
	    }
	  else
	    {
	      double val = 0;
	      if ( ! Helper::str2dbl( tok2[1] , &val ) )
		Helper::halt( "bad numeric input: " + tok[i] );
	      targs[ tok2[0] ] = val;
	    }
	}

      // if no channels, e.g. could be a time-track; denote that it is empty
      if ( tchs.size() == 0 )
	tchs.push_back( "." );
      // add each channel separately (w/ the same args)
      
      for (int c=0; c<tchs.size(); c++)
	{
	  suds_spec_t spec;	  
	  spec.ftr = lab2ftr[ Helper::toupper( tok[0] ) ];
	  spec.ch = tchs[c];
	  spec.arg = targs;

	  // but check we only have each feature/channel pair specified
	  // not more than once
	  if ( fcmap.find( spec.ftr ) != fcmap.end() )
	    {
	      std::map<std::string,suds_spec_t>::const_iterator ff = fcmap.find( spec.ftr )->second.find( spec.ch );
	      if ( ff != fcmap.find( spec.ftr )->second.end() )
		Helper::halt( "cannot specify feature/channel pair more than once" );
	    }

	  // track we've done this feature/channel combo
	  fcmap[ spec.ftr ][ spec.ch ] = spec;
	  
	  // add to the main list (in order)
	  specs.push_back( spec );
	}

      // process next line of model file
    } 

  
  //
  // check that NC was specified
  //
  
  if ( suds_t::nc == 0 )
    Helper::halt( "model file did not specify the number of components (NC)" );

  // make sure commands have the required arguments
  check_args();

  // track the implied total # of features
  suds_t::nf = cols();
  suds_t::ns = chs.size();
    
  logger << "  read " << specs.size() << " feature specifications ("
	 << suds_t::nf << " total features on " << suds_t::ns << " channels) from " << modelfile << "\n";
  
  // & construct the map of specs/channels to feature columns
  build_colmap();

  // set weights (default, or read from file)
  if ( winfile != "" )
    read_weights( winfile );
  else
    set_weights();  
  
  // or write (i.e. template to be editted for next read-weights
  if ( woutfile != "" )
    write_weights( woutfile ) ;

  return true;
}
      
void suds_model_t::default_model()
{
  //  channels used, and sample rates  (CH)
  //  the number of SVD components that will be extracted (NC)
  //  the features to construct the raw matrix (SPEC, etc)

  // suds_t::nc = 0;
  
  // // ensure we have initiated the maps
  // init();
  
  // if ( ! Helper::fileExists( modelfile ) )
  //   Helper::halt( "could not open " + modelfile );

  // // clear any current specifications
  // specs.clear();
  
  // std::ifstream IN1( modelfile.c_str() , std::ios::in );
  // while ( ! IN1.eof() )
  //   {
  //     std::string line;
  //     Helper::safe_getline( IN1 , line );
  //     if ( line == "" ) continue;
  //     if ( line[0] == '%' ) continue;
      
  //     // expecting format:
  //     //  CHANNEL CH SR 
  //     //  FEATURE  { CH } { KEY=VAL }
  //     //  i.e. either a channel, or (if has an '=') an argument

  //     std::vector<std::string> tok = Helper::parse( line , " \t" );
      
  //     if ( tok.size() < 2 ) Helper::halt( "bad format for line: " + line );

  //     // channel specifier?
  //     if ( Helper::toupper( tok[0] ) == "CH" )
  // 	{
  // 	  if ( tok.size() != 3 )
  // 	    Helper::halt( "expecing: CH label SR" );
	  
  // 	  int sr ;
  // 	  if ( ! Helper::str2int( tok[2] , &sr ) )
  // 	    Helper::halt( "bad format: " + line );
	  
  // 	  // store
  // 	  chs[ tok[1] ] = suds_channel_t( tok[1] , sr ) ;

  // 	  // next line
  // 	  continue;
  // 	}

  //     // component (NC) specifier?
  //     if ( Helper::toupper( tok[0] ) == "NC" )
  // 	{
  // 	  if ( tok.size() != 2 )
  //           Helper::halt( "expecing: NC value" );

  // 	  int nc ;
  // 	  if ( ! Helper::str2int( tok[1] , &nc ) )
  //           Helper::halt( "bad format: " + line );
  // 	  if ( nc < 1 )
  // 	    Helper::halt( "bad format: " + line );

  // 	  suds_t::nc = nc;
  // 	  continue;
  // 	}

      
  //     // feature specifier?

  //     if ( lab2ftr.find( Helper::toupper( tok[0] ) ) == lab2ftr.end() )
  // 	Helper::halt( "feature not recognized: " + tok[0] );

  //     // get list of channels, args separately
  //     std::vector<std::string> tchs;
  //     std::map<std::string,double> targs;

  //     for (int i=1; i<tok.size(); i++)
  // 	{
  // 	  std::vector<std::string> tok2 = Helper::parse( tok[i] , '=' );
	  
  // 	  if ( tok2.size() > 2 ) Helper::halt( "bad format: " + tok[i] );
	  
  // 	  // add as a channel
  // 	  if ( tok2.size() == 1 )
  // 	    {
  // 	      // has the channel already been specified via CH?
  // 	      if ( tok2[0] != "." && chs.find( tok2[0] ) == chs.end() )
  // 		Helper::halt( tok2[0] + " not specified via 'CH' yet: " + line );	      
  // 	      tchs.push_back( tok2[0] );
  // 	    }
  // 	  else
  // 	    {
  // 	      double val = 0;
  // 	      if ( ! Helper::str2dbl( tok2[1] , &val ) )
  // 		Helper::halt( "bad numeric input: " + tok[i] );
  // 	      targs[ tok2[0] ] = val;
  // 	    }
  // 	}

  //     // if no channels, e.g. could be a time-track; denote that it is empty
  //     if ( tchs.size() == 0 )
  // 	tchs.push_back( "." );
  //     // add each channel separately (w/ the same args)
      
  //     for (int c=0; c<tchs.size(); c++)
  // 	{
  // 	  suds_spec_t spec;	  
  // 	  spec.ftr = lab2ftr[ Helper::toupper( tok[0] ) ];
  // 	  spec.ch = tchs[c];
  // 	  spec.arg = targs;

  // 	  // but check we only have each feature/channel pair specified
  // 	  // not more than once
  // 	  if ( fcmap.find( spec.ftr ) != fcmap.end() )
  // 	    {
  // 	      std::map<std::string,suds_spec_t>::const_iterator ff = fcmap.find( spec.ftr )->second.find( spec.ch );
  // 	      if ( ff != fcmap.find( spec.ftr )->second.end() )
  // 		Helper::halt( "cannot specify feature/channel pair more than once" );
  // 	    }

  // 	  // track we've done this feature/channel combo
  // 	  fcmap[ spec.ftr ][ spec.ch ] = spec;
	  
  // 	  // add to the main list (in order)
  // 	  specs.push_back( spec );
  // 	}
  //   }

  // // check that NC was specified
  // if ( suds_t::nc == 0 )
  //   Helper::halt( "model file did not specify the number of components (NC)" );

  // // make sure commands have the required arguments
  // check_args();

  // // track the implied total # of features
  // suds_t::nf = cols();
  // suds_t::ns = chs.size();
    
  // logger << "  using default model, with " << specs.size() << " feature specifications ("
  // 	 << suds_t::nf << " total features on " << suds_t::ns << " channels)\n";
  
  // // construct the map of specs/channels to feature columns
  // build_colmap();
  
  // // use default weights
  // set_weights();  
  
}

bool suds_model_t::write( const std::string & modelfile )
{
  std::ofstream OUT1( modelfile.c_str() , std::ios::out );
  for (int i=0; i<specs.size(); i++)
    {
      // class
      OUT1 << ftr2lab[ specs[i].ftr ] ;
      
      // channel
      OUT1 << specs[i].ch ;
      
      // args
      std::map<std::string,double>::const_iterator ii = specs[i].arg.begin();
      while ( ii != specs[i].arg.end() )
	{
	  OUT1 << "\t" << ii->first << "=" << ii->second;
	  ++ii;
	}
      OUT1 << "\n";
    }
  OUT1.close();

  return true;
}


int suds_model_t::cols() const
{
  int n = 0;
  for (int i=0; i<specs.size(); i++)
    specs[i].cols(&n);
  return n;
}

void suds_model_t::check_args()
{
  for (int i=0; i<specs.size(); i++)
    {
      suds_spec_t & spec = specs[i];
      
      if ( spec.ftr == suds_feature_t::SUDS_LOGPSD ||
	   spec.ftr == suds_feature_t::SUDS_RELPSD ||
	   spec.ftr == suds_feature_t::SUDS_CVPSD )
	{
	  if ( spec.arg.find( "lwr" ) == spec.arg.end() )
	    Helper::halt( ftr2lab[ suds_feature_t::SUDS_LOGPSD ] + " requires 'lwr' arg" );
	  if ( spec.arg.find( "upr" ) == spec.arg.end() )
	    Helper::halt( ftr2lab[ suds_feature_t::SUDS_LOGPSD ] + " requires 'upr' arg" );
	  if ( spec.arg[ "lwr" ] > spec.arg[ "upr" ] )
	    Helper::halt( ftr2lab[ suds_feature_t::SUDS_LOGPSD ] + " requires 'lwr' < 'upr' " );
	  if ( spec.arg[ "lwr" ] <= 0 || spec.arg[ "upr" ] <= 0 )
	    Helper::halt( ftr2lab[ suds_feature_t::SUDS_LOGPSD ] + " requires 'lwr' and 'upr' to be > 0 " );
	}

      // the z-lwr/z-upr range does not need to overlap lwr/upr range for RELPSD
      if ( spec.ftr == suds_feature_t::SUDS_RELPSD )
	{
	  if ( spec.arg.find( "z-lwr" ) == spec.arg.end() )
            Helper::halt( ftr2lab[ suds_feature_t::SUDS_LOGPSD ] + " requires 'z-lwr' arg" );
          if ( spec.arg.find( "z-upr" ) == spec.arg.end() )
            Helper::halt( ftr2lab[ suds_feature_t::SUDS_LOGPSD ] + " requires 'z-upr' arg" );
          if ( spec.arg[ "z-lwr" ] > spec.arg[ "z-upr" ] )
            Helper::halt( ftr2lab[ suds_feature_t::SUDS_LOGPSD ] + " requires 'z-lwr' < 'z-upr' " );
          if ( spec.arg[ "z-lwr" ] <= 0 || spec.arg[ "z-upr" ] <= 0 )
            Helper::halt( ftr2lab[ suds_feature_t::SUDS_LOGPSD ] + " requires 'z-lwr' and 'z-upr' to be > 0 " );
	}

      // time-tracks
      if ( spec.ftr == suds_feature_t::SUDS_TIME )
	{
	  if ( spec.arg.find( "order" ) == spec.arg.end() )
            Helper::halt( ftr2lab[ suds_feature_t::SUDS_TIME ] + " requires 'order' arg" );
	}

      // smoothing/denoising
      if ( spec.ftr == suds_feature_t::SUDS_DENOISE )
	{
	  if ( spec.arg.find( "lambda" ) == spec.arg.end() )
            Helper::halt( ftr2lab[ suds_feature_t::SUDS_DENOISE ] + " requires 'lambda' arg" );	  
	}
      
      if ( spec.ftr == suds_feature_t::SUDS_SMOOTH )
	{
	  if ( spec.arg.find( "half-window" ) == spec.arg.end() )
            Helper::halt( ftr2lab[ suds_feature_t::SUDS_SMOOTH ] + " requires 'half-window' (epochs) arg" );
	}

      if ( spec.ftr == suds_feature_t::SUDS_DENOISE2 )
	{
	  if ( spec.arg.find( "lambda" ) == spec.arg.end() )
            Helper::halt( ftr2lab[ suds_feature_t::SUDS_DENOISE2 ] + " requires 'lambda' arg" );	  
	}
      
      if ( spec.ftr == suds_feature_t::SUDS_SMOOTH2 )
	{
	  if ( spec.arg.find( "half-window" ) == spec.arg.end() )
            Helper::halt( ftr2lab[ suds_feature_t::SUDS_SMOOTH2 ] + " requires 'half-window' (epochs) arg" );
	}

      
    }
}
  
void suds_model_t::build_colmap()
{
  ftr2ch2col.clear();
  int n = 0;
  for (int i=0; i<specs.size(); i++)
    {
      int start = n;
      specs[i].cols(&n);
      int end = n; // one past end
      const suds_feature_t ftr = specs[i].ftr;
      const std::string & ch = specs[i].ch;

      // for smooth/denoise, n could be 0;
      // doesn't add new cols, but implies we should point back
      // to the original cols
      if ( start == end ) 
       	{
	  std::vector<int> origs(n);
	  for (int i=0;i<n;i++) origs[i] = i;
       	  ftr2ch2col[ ftr ][ ch ] = origs;
	}
      else
	{
	  for (int j=start; j<end; j++)
	    ftr2ch2col[ ftr ][ ch ].push_back( j );
	}
    }
}


// give columns for a spec/channel combo
bool suds_model_t::has( suds_feature_t ftr , const std::string & ch )
{
  std::map<suds_feature_t,std::map<std::string,std::vector<int> > >::const_iterator cc = ftr2ch2col.find( ftr );
  if ( cc == ftr2ch2col.end() ) return false;
  std::map<std::string,std::vector<int> >::const_iterator dd = cc->second.find( ch );
  return dd != cc->second.end(); 
}

std::vector<int> suds_model_t::cols( suds_feature_t ftr , const std::string & ch )
{
  std::vector<int> dummy;
  std::map<suds_feature_t,std::map<std::string,std::vector<int> > >::const_iterator cc = ftr2ch2col.find( ftr );
  if ( cc == ftr2ch2col.end() ) return dummy;
  std::map<std::string,std::vector<int> >::const_iterator dd = cc->second.find( ch );
  if ( dd == cc->second.end() ) return dummy;
  return dd->second;
}


std::vector<std::string> suds_model_t::labels()
{
  std::vector<std::string> l;
  
  for (int i=0; i<specs.size(); i++)
    {
      
      const suds_feature_t ftr = specs[i].ftr; 
      const std::string & ch = specs[i].ch;      
      const std::string l0 = suds_model_t::ftr2lab[ ftr ];
      
      if ( ftr == SUDS_LOGPSD
	   || ftr == SUDS_RELPSD
	   || ftr == SUDS_CVPSD )
	{
	  double lwr = specs[i].arg[ "lwr" ] ;
	  double upr = specs[i].arg[ "upr" ] ;
	  int n = ( upr - lwr ) / suds_t::spectral_resolution + 1;	  	  
	  for (int i=0; i<n; i++)
	    l.push_back( l0 + "_" + ch + "_" + Helper::dbl2str( lwr + i * suds_t::spectral_resolution ) );
	}
      
      
      else if ( ftr == SUDS_SLOPE
		|| ftr == SUDS_SKEW
		|| ftr == SUDS_KURTOSIS
		|| ftr == SUDS_FD
		|| ftr == SUDS_MEAN )
	{	  
	    l.push_back( l0 + "_" + ch );
	}

      else if ( ftr == SUDS_HJORTH )
	{
	  l.push_back( l0 + "2_" + ch );
	  l.push_back( l0 + "3_" + ch );	  
	}
      
      else if ( ftr == SUDS_PE )
	{
	  l.push_back( l0 + "3_" + ch );
	  l.push_back( l0 + "4_" + ch );
	  l.push_back( l0 + "5_" + ch );
	  l.push_back( l0 + "6_" + ch );
	  l.push_back( l0 + "7_" + ch );	
	}
      

      // duplicate current set :

      if ( ftr == SUDS_SMOOTH2 || ftr == SUDS_DENOISE2 )
	{
	  std::vector<std::string> c = l;
	  for (int i=0; i<c.size(); i++)
	    l.push_back( suds_model_t::ftr2lab[ ftr ] + "_" + c[i] );
	}

      // replace current set 
      if ( ftr == SUDS_SMOOTH || ftr == SUDS_DENOISE )
	{
	  for (int i=0; i<l.size(); i++)
	    l[i] = suds_model_t::ftr2lab[ ftr ] + "_" + l[i] ;
	}
      
      // time-track
      if ( ftr == SUDS_TIME )
	{	  
	  int n = specs[i].arg[ "order" ] ;
	  for (int i=0; i<n; i++)
	    l.push_back( l0 + Helper::int2str( i+1 ) );
	}      

    }
      
   return l;    
}


// dump weights to a file
void suds_model_t::write_weights( const std::string & weightfile )
{

  logger << "  writing feature weights to " << weightfile << "\n";

  const std::vector<std::string> l = labels();

  // check
  if (l.size() != W.size() )
    Helper::halt( "internal error in suds_model_t::write_weights()" );

  std::ofstream W1( weightfile.c_str() , std::ios::out );
  for (int i=0; i<l.size(); i++)
    W1 << l[i] << "\t" << W[i] << "\n";
  W1.close();  
}

// read weights from a file
void suds_model_t::read_weights( const std::string & weightfile )
{
  logger << "  reading feature weights from " << weightfile << "\n";
  
  // here we do not check weights... we just check that the total N is
  // as expected, i.e. given the modelfile.   i.e. it is the user's
  // responsibility to not mess things up

  const int n = cols();

  if ( ! Helper::fileExists( weightfile ) )
    Helper::halt( "could not open " + weightfile );

  std::vector<double> wt;
  std::ifstream W1( weightfile.c_str() , std::ios::in );
  while ( ! W1.eof() )
    {      
      std::string dummyvar;
      double w1;
      W1 >> dummyvar >> w1;
      if ( W1.eof() || W1.bad() ) break;
      wt.push_back( w1 );
    }
  W1.close();

  if ( wt.size() != n )
    Helper::halt( "expecting " + Helper::int2str(n)
		  + " but read " + Helper::int2str( (int)wt.size() )
		  + " weights from " + weightfile );

  // update main weight vector
  W.resize( n );
  for (int i=0; i<n; i++) W[i] = wt[i];
  
}



// return implied number of columns

int suds_spec_t::cols( int * t ) const
{
  
  // PSD is stratified by frequencu
  if ( ftr == SUDS_LOGPSD
       || ftr == SUDS_RELPSD
       || ftr == SUDS_CVPSD )
    {
      double lwr = arg.find( "lwr" )->second ;
      double upr = arg.find( "upr" )->second ;
      int n = ( upr - lwr ) / suds_t::spectral_resolution + 1 ;
      *t += n;
      return n;
    }
  
  // 1 column per channel
  if ( ftr == SUDS_SLOPE
       || ftr == SUDS_SKEW
       || ftr == SUDS_KURTOSIS
       || ftr == SUDS_FD
       || ftr == SUDS_MEAN )
    {
      *t += 1;
      return 1;
    }
  
  // 2 values per channel
  if ( ftr == SUDS_HJORTH )
    {
      *t += 2 ;
      return 2 ;
    }
  
  // PE is 3..7
  if ( ftr == SUDS_PE )
    {
      *t += 5 ;
      return 5 ;
    }
  
  // doubles current set
  if ( ftr == SUDS_SMOOTH2 || ftr == SUDS_DENOISE2 )
    {
      int n = *t;
      *t += n;
      return n;
    }

  // replaces current set
  if ( ftr == SUDS_SMOOTH || ftr == SUDS_DENOISE )
    {
      // i.e. total unchanged
      int n = *t;
      return n;
    }

  
  // time-track
  if ( ftr == SUDS_TIME )
    {
      int n = arg.find( "order" )->second;
      if ( n < 0 )
	Helper::halt( "invalid value for TIME order (0-10)" );
      if ( n > 10 )
	Helper::halt( "invalid value for TIME order (0-10)" );
      *t += n;
      return n;
    }
  
  // unknown/error at this point

  Helper::halt( "could not process model file / extracting implied col count" );
  return 0;
  
}

void suds_model_t::set_weights()
{
  // get proper size for weight matrix
  W.resize( suds_t::nf );
  
  // default is that each domain should sum to 1.0
  // except, if a duplocate (ie. from smoothing, then keep as the original)

  std::vector<std::string> l = labels();
  
  int n = 0;
  int p = 0;
  for (int i=0; i<specs.size(); i++)
    {
      // if replace/smooth, do nothing
      if ( specs[i].ftr == suds_feature_t::SUDS_SMOOTH ) continue;
      if ( specs[i].ftr == suds_feature_t::SUDS_DENOISE ) continue;

      // duplicate weights?
      if (  specs[i].ftr == suds_feature_t::SUDS_SMOOTH2 ||
	    specs[i].ftr == suds_feature_t::SUDS_DENOISE2 )
	{	  
	  // just double whatever we already have done.. 
	  const int q = p;	  
	  for (int i=0; i<q; i++)
	    W[ p++ ] = W[ i ];	  
	}
      else // normal features
	{
	  int n1 = specs[i].cols(&n);
	  for (int j=0; j<n1; j++)
	    W[ p++ ] = 1.0 / (double) n1 ; 
	}
    }
  
}

