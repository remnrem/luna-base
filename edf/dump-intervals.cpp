
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

#include "edf/edf.h"
#include "edf/slice.h"

#include "helper/helper.h"
#include "helper/logger.h"

#include "defs/defs.h"
#include "intervals/intervals.h"

extern logger_t logger;

void dump_intervals( const std::string & ints , 
		     const std::string & edfs )
{
  
  //
  // Map ID to EDF file
  //
  
  std::map<std::string,std::string> id2edf;
  
  std::ifstream EDFLIST( edfs.c_str() );
  while ( ! EDFLIST.eof() )
    {      
      std::string line;
      Helper::safe_getline( EDFLIST , line );
      if ( EDFLIST.eof() || line == "" ) continue;
      
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      
      if ( tok.size() < 2 ) 
	Helper::halt( "requires (ID) | EDF file | (optional ANNOT files)" );
      
      id2edf[ tok[0] ] = tok[1];
    }
  
  EDFLIST.close();
  
  //
  // Get interval input from STDIN
  //

  // each line will correspond to one 'chunk' of output
  
  std::map<std::string,std::vector<feature_t> > id2feature;

  std::ifstream INTLIST( ints.c_str() );
  
  while ( ! INTLIST.eof() ) 
    {
      std::string line;
      Helper::safe_getline( INTLIST , line );
      if ( INTLIST.eof() || line == "" ) continue;
      
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      // expecting format: (tab-delim)
      
      // 0 ID
      // 1 signal/channel
      // 2 label 
      // 3 feature start (MSEC)
      // 4 feature stop  (MSEC)
      // 5 window size (seconds around center point, add w/2)
      
      std::string indiv  = tok[0];

      feature_t f;      
      f.signal = tok[1];
      f.label  = tok[2];	

      uint64_t fstart = 0 , fstop = 0;
      double window_sec = 0;

      if ( tok.size() < 6 ) 
	Helper::halt( "requires ID|signal(s)|label|start|stop|window" );
      
      if ( ! Helper::str2int64( tok[3] , &fstart ) ) 
	Helper::halt( "bad fstart value" );

      if ( ! Helper::str2int64( tok[4] , &fstop ) ) 
	Helper::halt( "bad fstop value" );
     
      if ( ! Helper::str2dbl( tok[5] , &window_sec ) ) 
	Helper::halt( "bad fstop value" );
      
      if ( fstop < fstart ) Helper::halt( "fstop < fstart" );
      
      if ( id2edf.find( indiv ) == id2edf.end() )
	Helper::halt( "could not find individual " + indiv );
     
      
      uint64_t mid = (fstop+fstart)/2; 
      f.feature = interval_t( fstart , fstop );

      uint64_t wstart = 0; 
      uint64_t wstop = 0;
      
      double mid_sec = (double)mid / globals::tp_1sec; 
      double whalf = window_sec/2.0;
      
      if ( whalf > mid_sec ) wstart = 0;
      else wstart = mid - uint64_t(window_sec*5e5);
                  
      wstop = mid + uint64_t(window_sec*5e5);

      f.window  = interval_t( wstart , wstop );

      id2feature[ indiv ].push_back( f );
      
    }
  
  INTLIST.close();

  //
  // Header
  //

  std::cout << "REC" << "\t"
	    << "ID" << "\t"
	    << "N" << "\t"
	    << "F" << "\t"
	    << "T" << "\t"
	    << "Y" << "\n";
  
  //
  // Now iterate over each EDF required
  //
  
  int rec = 0;

  std::map<std::string,std::vector<feature_t> >::const_iterator ii = id2feature.begin();

  while ( ii != id2feature.end() )
    {
      
      std::string indiv     = ii->first;
      std::string edffile   = id2edf[ indiv ];
      
      //
      // get all required signals
      //
      
      std::set<std::string> inp_signals;      
      for (int f=0;f<ii->second.size();f++)
	{
	  std::vector<std::string> ss = Helper::parse( ii->second[f].signal , "|" );
	  for (int k=0;k<ss.size();k++) inp_signals.insert( ss[k] );
	}
      
      
      //
      // Load EDF
      //

      edf_t edf;
      
      bool okay = edf.attach( edffile , indiv , &inp_signals );
      
      if ( ! okay ) 
	{
	  logger << "problem loading " << edffile << ", skipping...\n";
	  continue;
	}
      
      //
      // For each requested feature
      //
      
      for (int f=0;f<ii->second.size();f++)
	{

	  std::vector<std::string> ss = Helper::parse( ii->second[f].signal , "|" );	  

	  const interval_t & interval = ii->second[f].window;
	  const interval_t & feature  = ii->second[f].feature;
	  
	  for (int s=0;s<ss.size();s++) 
	    {
	      
	      int signal = edf.header.signal( ss[s] );
	      
	      if ( signal == -1 ) continue;
	      
	      //
	      // Extract raw signal
	      //
	      
	      slice_t slice( edf , signal , interval );
	      
	      const std::vector<double> * d = slice.pdata();
	      
	      const std::vector<uint64_t> * pt = slice.ptimepoints();
	      
	      //
	      // Display
	      //
	      
	      int np = d->size();
	      
	      for (int o=0;o<np;o++)
		{
		  std::cout << rec << "\t"
			    << indiv << "\t"
			    << f << "\t"
			    << ( feature.contains( (*pt)[o] ) ? 1 : 0 ) << "\t" 
			    << (*pt)[o] << "\t"
			    << (*d)[o] << "\n";
		  
		}
	      
	      // next feature
	      ++rec;
	      
	      // done
	      logger << "processed : " << indiv << ", " << f+1 << " features\n";
	    }
	}

      ++ii;
    }
  
  
}



