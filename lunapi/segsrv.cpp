
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

// TODOs
//   define viz-epochs
//   define viz-clock mappings (on 0..1 scale) 
//   spectral / Hjorth epoch-level summaries
//   decimation on-the-fly and/or summaries at higher time-scales
//   bandpass filtering on-the-fly? (or at creation?)
//   sleep stages

#include "lunapi/segsrv.h"
#include "lunapi/lunapi.h"

segsrv_t::segsrv_t( lunapi_inst_ptr inst ) : p( inst ) 
{
  awin = bwin = 0;
}

void segsrv_t::init()
{
  awin = bwin = 0;
  aidx.clear();
  bidx.clear();
  segments = p->edf.timeline.segments();
  gaps = p->edf.timeline.gaps( segments ); 
  //  std::cout << " found " << segments.size() << " segs and " << gaps.size() << " gaps\n";
}


int segsrv_t::populate( const std::vector<std::string> & chs , const std::vector<std::string> & anns )
{

  // for now, can only init this once
  //  (i.e. do a single data pull - this avoids if the internal object is changed
  //   after calling, but the view is still open)
  //  that is, we only touch the EDF/instance on creation, or after update w/ 
  
  init();
  
  int count = 0;
  for (int i=0; i<chs.size(); i++)
    if ( add_channel( chs[i] ) ) ++count;
   
  return count;
}


std::set<std::pair<double,double> > segsrv_t::get_gaps() const
{

  // current window = awin, bwin
  // segments included

  uint64_t atp = awin * globals::tp_1sec ;
  uint64_t btp = bwin * globals::tp_1sec ;

  std::set<std::pair<double,double> > g;
  
  std::set<interval_t>::const_iterator gg = gaps.begin();
  while ( gg != gaps.end() )
    {
      // std::cout << " considering gap: " << gg->as_string() << "\n";
      // std::cout << " window = " << atp << " " << btp << "\n";
      // does any of the window fall in this gap?
      if ( btp > gg->start && atp < gg->stop )
	{
	  //std::cout << " this gap in window\n";
	  uint64_t start1 = atp > gg->start ? atp : gg->start;
	  uint64_t stop1 = btp < gg->stop ? btp : gg->stop;
	  g.insert( std::pair<double,double>( start1 * globals::tp_duration , stop1 * globals::tp_duration ) );
	  //std::cout << "  --> add " << start1 * globals::tp_duration << " " << stop1 * globals::tp_duration << "\n";
	}
      ++gg;
    }

  return g;
  
}

bool segsrv_t::add_channel( const std::string & ch )
{
  
  const int slot = p->edf.header.signal( ch );
  if ( slot == -1 ) return false;

  // sample rate
  int sr = p->edf.header.sampling_freq( slot );
  //  std::cout <<" ch " << ch << " " << slot << " " << sr << "\n";
  
  // get all data
  slice_t slice( p->edf , slot , p->edf.timeline.wholetrace() );
  const std::vector<double> * data = slice.pdata();
  const int n = data->size();
  Eigen::VectorXf d = Eigen::VectorXf::Zero( n );
  for (int i=0; i<n; i++) d[i] = (*data)[i];
  sigmap[ ch ] = d;

  // store
  srmap[ ch ] = sr;

  // do we already have a time-track?
  if ( tidx.find( sr ) == tidx.end() )
    {      
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      std::map<double,int> tt;
      Eigen::VectorXf ts = Eigen::VectorXf::Zero( n );

      // for quick lookup: map of time-in-sec --> idx-int
      for (int i=0; i<tp->size(); i++)
	{
	  const double sec = (*tp)[i] * globals::tp_duration;
	  tt[ sec ] = i;
	  ts[ i ] = sec;
	}
      
      // store
      tidx[ sr ] = tt;

      // also store original (to return actual times)
      tmap[ sr ] = ts;
      
    }
    
  return true;
}

  
// set window
bool segsrv_t::set_window( double a , double b )
{

  // max time (seconds, 1-tp-unit past end) 
  const double tmax = p->last_sec();

  //  std::cout <<" tmax = " << tmax << "\n";
  
  // store seconds 
  awin = a < 0 ? 0 : ( a > tmax ? tmax : a ) ;
  bwin = b < 0 ? 0 : ( b > tmax ? tmax : b ) ;

  if ( awin > bwin )
    {
      double tmp = bwin;
      bwin = awin;
      awin = tmp;
    }

  
  std::set<int> srs;
  std::map<std::string,int>::const_iterator cc = srmap.begin();
  while ( cc != srmap.end() )
    {
      srs.insert( cc->second );
      ++cc;
    }

  bool all_okay = true;

  std::set<int>::const_iterator ss = srs.begin();
  while ( ss != srs.end() )
    {
      
      int aa = 0, bb = 0;
      const bool okay = get_tidx( awin, bwin , *ss , &aa , &bb );
      
      if ( okay )
	{
	  aidx[ *ss ] = aa;
	  bidx[ *ss ] = bb;
	  //std::cout << " set win " << awin << " " << bwin << " --> " << aa << " " << bb << "\n";
	}
      else
	all_okay = false;
      
      ++ss;
    }
  
  return all_okay;
}

Eigen::VectorXf segsrv_t::get_timetrack( const std::string & ch ) const
{
  std::map<std::string,int>::const_iterator ss = srmap.find( ch );
  if ( ss == srmap.end() )
    {
      Eigen::VectorXf empty;
      return empty;
    }
  int sr = ss->second;

  // return SR-specific time-track
  // (which may contain gaps etc so need whole thing)
  
  const Eigen::VectorXf & tt = tmap.find( sr )->second;
  const int aa = aidx.find(sr)->second;
  const int bb = bidx.find(sr)->second;  
  return tt.segment( aa , bb - aa );

}


Eigen::VectorXf segsrv_t::get_signal( const std::string & ch ) const 
{

  std::map<std::string,int>::const_iterator ss = srmap.find( ch );

  if ( ss == srmap.end() )
    {
      Eigen::VectorXf empty;
      return empty;
    }
  
  int sr = ss->second;

  // return data (signal)
  const Eigen::VectorXf & data = sigmap.find( ch )->second;
  const int aa = aidx.find(sr)->second;
  const int bb = bidx.find(sr)->second;  
  return data.segment( aa , bb - aa );

}


// given two times and a sample rate, get indices
bool segsrv_t::get_tidx( double a, double b , int sr , int * aa, int *bb ) const
{

  if ( tidx.find( sr ) == tidx.end() ) return false;

  const std::map<double,int> & ts = tidx.find( sr )->second ;
    
  // iterator equal/greater than start
  std::map<double,int>::const_iterator abound = ts.lower_bound( a );
  if ( abound == ts.end() ) return false;
  
  // one-past the end
  std::map<double,int>::const_iterator bbound = ts.lower_bound( b );
  if ( bbound == ts.end() ) return false;
  
  // std::cout  << "window " << abound->first << " " << bbound->first
  // 	     << " --> " << abound->second << " " << bbound->second << "\n";
  
  *aa = abound->second;
  *bb = bbound->second;
  return true;
    
}

