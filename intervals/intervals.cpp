
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

#include "intervals.h"

#include "../helper/helper.h"
#include "../edf/edf.h"
#include "../defs/defs.h"

#include <fstream>
#include <string>
#include <map>
#include <set>

#include <iostream>

std::ostream & operator<<( std::ostream & out , const interval_t & rhs ) 
{
  out << rhs.start << "-" << rhs.stop;
  return out;
}




int interval_t::intersect( const std::set<interval_t> & a, 
			   const std::set<interval_t> & b, 
			   std::set<interval_t> * botha ,
			   std::set<interval_t> * bothb ,
			   std::set<interval_t> * cons ,
			   std::set<interval_t> * uns  ,
			   std::set<interval_t> * onlya ,
			   std::set<interval_t> * onlyb , 
			   double th , uint64_t win )
{
  
  
  // although this typically won't be the case, we'll allow for
  // overlapping intervals within list 'a' or 'b'
  
  // note -- however, in the case of within-list overlap, we probably 
  // want to handle this differently (i.e. to decide whether to report
  // merged intervals (based on either consensus/union) within list
  // as it stands, we will get a different number for each pair (i.e. could
  // be a different number of 
  
  if ( a.size() == 0 || b.size() == 0 ) 
    {
      *onlya = a;
      *onlyb = b;
      botha->clear();
      bothb->clear();
      cons->clear();
      uns->clear();
      return 0;
    }
  
  botha->clear();
  bothb->clear();
  onlya->clear();
  onlyb->clear();
  cons->clear();
  uns->clear();
  
  std::set<interval_t>::const_iterator aa = a.begin();
  std::set<interval_t>::const_iterator bb = b.begin();
  
  
  std::set<interval_t> current;
  
  
  // merge overlapping intervals
  std::map<interval_t,std::set<interval_t> > supera;
  std::map<interval_t,std::set<interval_t> > superb;

  uint64_t last = 0LL;
  uint64_t first = 0LL;
  std::set<interval_t> pool;
  while ( aa != a.end() )
    {
      interval_t wa = *aa;
      wa.expand( win );
      
      // if pool is empty, start adding
      if ( pool.size() == 0 ) 
	{
	  pool.insert( *aa );
	  first = wa.start;
	  last = wa.stop;	  
	}
      else
	{
	  // does this next one overlap?
	  if ( wa.start < last ) // yes
	    {
	      pool.insert( *aa );
	      if ( wa.stop > last ) last = wa.stop;
	    }
	  else // if no overlap, then add pool to super-list
	    {
	      supera[ interval_t( first , last ) ] = pool;	      
	      pool.clear();
	      // and make a new pool
	      pool.insert( *aa );
	      first = wa.start;
	      last = wa.stop;	  	      
	    }

	}
      ++aa;
    }
  
  // add final pool if not empty
  if ( pool.size() > 0 ) supera[ interval_t( first , last ) ] = pool;
  
  last = 0;
  first = 0;
  pool.clear();
  while ( bb != b.end() )
    {
      interval_t wb = *bb;
      wb.expand( win );
      
      // if pool is empty, start adding
      if ( pool.size() == 0 ) 
	{
	  pool.insert( *bb );
	  first = wb.start;
	  last = wb.stop;	  
	}
      else
	{
	  // does this next one overlap?
	  if ( wb.start < last )
	    {
	      pool.insert( *bb );
	      if ( wb.stop > last ) last = wb.stop;
	    }
	  else // if no overlap, then add pool to super-list
	    {
	      superb[ interval_t( first , last ) ] = pool;
	      pool.clear();
	      // and make a new pool
	      pool.insert( *bb );
	      first = wb.start;
	      last = wb.stop;	  	      

	    }
	}
      ++bb;
    }

  // add final pool if not empty
  if ( pool.size() > 0 ) superb[ interval_t( first , last ) ] = pool;

  // now we know that super-intervals do not overlap
  
  std::map<interval_t,std::set<interval_t> >::const_iterator saa = supera.begin();
  std::map<interval_t,std::set<interval_t> >::const_iterator sbb = superb.begin();
  
  while ( 1 ) 
    {

      if ( saa == supera.end() ) break;
      if ( sbb == superb.end() ) break;

      // std::cout << "SO " << saa->first.start << "-" << saa->first.stop << "(" << saa->second.size() << ")" 
      // 		<< " &  " 
      // 		<< sbb->first.start << "-" << sbb->first.stop << "(" << sbb->second.size() << ")"
      // 		<< "\n";
      
      // if using a 'window' around each event, expand here (i.e. based on the 'super' version )
      
      interval_t wa = saa->first;
      interval_t wb = sbb->first;
      // wa.expand( win );
      // wb.expand( win );
            
      // any overlap?
      if ( wa.overlaps( wb ) )
	{
	  
	  // consider the all-by-all possible overlaps now
	  const std::set<interval_t> & sa = saa->second;
	  const std::set<interval_t> & sb = sbb->second;
	  
	  std::set<interval_t>::const_iterator ii = sa.begin();
	  while ( ii != sa.end() )
	    {

	      interval_t wi = *ii;
	      wi.expand( win );

	      uint64_t astart = wi.start;
	      uint64_t astop  = wi.stop;
	      
	      std::set<interval_t>::const_iterator jj = sb.begin();
	      while ( jj != sb.end() )
		{
		  
		  interval_t wj = *jj;
		  wj.expand( win );
		    
		  // any overlap?
		  if ( wi.overlaps( wj ) )
		    {

		      // meets definition?
		      uint64_t bstart = wj.start;
		      uint64_t bstop = wj.stop;
		      
		      uint64_t constart = astart < bstart ? bstart : astart;
		      uint64_t constop  = astop  < bstop  ? astop  : bstop;
		      uint64_t conlen   = constop - constart + 1;
		      
		      uint64_t unionstart = astart < bstart ? astart : bstart;
		      uint64_t unionstop  = astop   < bstop ? bstop  : astop;
		      uint64_t unionlen   = unionstop - unionstart + 1;

		      double olap = (double)conlen / (double)unionlen;
		      if ( olap >= th ) 
			{
			  // std::cout << "adding " << astart << " " << astop << " ---- " 
			  // 	    << bstart << " "<< bstop << "\n";
			  cons->insert( interval_t( constart , constop ) );
			  uns->insert( interval_t( unionstart , unionstop ) );

			  // save unwindowed versions...
			  botha->insert( *ii );
			  bothb->insert( *jj );
			  
			}		      
		    }
		  ++jj;		  
		}
	      
	      ++ii;
	    }  
	} // end of overlap check at super-level

      // advance one of the blocks -- the one with the first stop
      if ( saa->first.stop < sbb->first.stop ) ++saa; else ++sbb;

    } // next set of super-sets

  //
  // now form the non-overlapping lists (i.e. any in 'a' that didn't make it into the alap or blap sets
  //

  aa = a.begin();
  while ( aa != a.end() )
    {
      if ( botha->find( *aa ) == botha->end() ) onlya->insert( *aa );
      ++aa;
    }

  bb = b.begin();
  while ( bb != b.end() )
    {
      if ( bothb->find( *bb ) == bothb->end() ) onlyb->insert( *bb );
      ++bb;
    }
  
  return cons->size();
}
