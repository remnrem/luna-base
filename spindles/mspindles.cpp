
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


#include "spindles.h"
#include "mspindles.h"
#include "edf/edf.h"
#include "db/db.h"

#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger;

void mspindles_t::add( const std::vector<spindle_t> & s , int fs , uint64_t len , int fc , int c , const std::string & label )
{
  S.push_back(s) ; 
  mins.push_back( ( len / (double)fs ) / 60.0 );
  frq.push_back( fc );
  ch.push_back( c );   
  run_label.push_back( label );
}


void mspindles_t::collate()
{
  
  // We have all spindles in S
  
  // from S (spindle_t), populate M (mspindle_t) by appropriately merging spindles
  
  M.clear();

  // list out all spindles, and merge
  std::set<sort_t> all;
  
  const int n = S.size();
  
  for (int i=0; i<n; i++) 
    {
      const double fc = frq[i];
      const int c  = ch[i];
      const std::vector<spindle_t> & sp = S[i];
      const std::string & l = run_label[i];

      std::vector<spindle_t>::const_iterator ii = sp.begin();
      
      while ( ii != sp.end() )
	{
	  sort_t s( ii->tp , fc , c , i, l , &(*ii) );	  
	  all.insert(s);
	  ++ii;
	}
     
      // next spindle set
    }
  
  //
  // now iterate over all time-sorted values
  //
  
  std::vector<sort_t> overlaps;
  uint64_t last_sp;
  uint64_t first_sp;

  std::set<sort_t>::const_iterator ss = all.begin();
  while ( ss != all.end() )
    {
      
      //std::cout << "check " << overlaps.size() << "\t" << ss->i.start << "\t" << ss->i.stop << "\n";

      // if empty, start an overlap bag
      if ( overlaps.empty() ) 
	{
	  overlaps.push_back( *ss );
	  first_sp = ss->i.start;
	  last_sp = ss->i.stop;	 
	  ++ss; // next spindle
	  //	  std::cout << " ... new group\n";
	  continue;
	}

      
      // first check for /any/ overlap
      bool has_overlap = ss->i.start <= last_sp;
      
      //std::cout << " OL CHK " << ss->i.start << "  " <<  last_sp << "\t" << has_overlap << "\n";
      
      // then refine
      if ( has_overlap )
	{
	  uint64_t max_start = first_sp;
	  uint64_t min_start = ss->i.start;  // i.e. because sorted
	  
	  uint64_t min_stop = last_sp < ss->i.stop ? last_sp : ss->i.stop;
	  uint64_t max_stop = last_sp < ss->i.stop ? ss->i.stop : last_sp ; 
	  
	  uint64_t o_intersection = min_stop - min_start + 1;
	  uint64_t o_union        = max_stop - max_start + 1;
	  
	  double metric = o_intersection / (double)o_union;	      
	  
	  has_overlap = metric > interval_th; 
	  
	  // or if one spindle is more or less completely subsumed by the other(s)
	  if ( ( last_sp - first_sp + 1 ) / o_intersection >= 0.8 ) has_overlap = true;
	  if ( ( ss->i.stop - ss->i.start + 1 ) / o_intersection >= 0.8 ) has_overlap = true;		  
	  
	}
      
      //
      // if we have overlap as defined, add to the group
      //

      if ( has_overlap )
	{		  
	  overlaps.push_back( *ss );
	  // extend range
	  if ( ss->i.stop > last_sp ) last_sp = ss->i.stop;		 
	}
      
          
      //
      // otherwise, add this group (which may be only a single event)
      //
      
      else
	{
	  //std::cout << " ... now proc'ing\n";
	  
	  proc_overlaps( overlaps );
	  
	  // clean up
	  overlaps.clear();
	  
	  // and add the next spindle (i.e. current one)
	  overlaps.push_back( *ss );
	  first_sp = ss->i.start;
	  last_sp = ss->i.stop;
	}
      
      ++ss;
    }
  
  
  //
  // empty the last bag here too..
  //
  
  proc_overlaps( overlaps );
  

  // 
  // Now, M should have been populated (from S) which contains the final mspindles_t set
  //  
	  
}
  

void mspindles_t::output( const signal_list_t & signals )
{
  
  //
  // Maximum number of minutes (for density)
  //

  double max_mins = 0;
  for (int i=0;i<mins.size();i++) if ( mins[i] > max_mins ) max_mins = mins[i]; 

  
  
  //
  // Report final list size
  //
  
  writer.var( "MSP_N" , "Number of merged spindles" );
  writer.var( "MSP_DENS" , "Merged spindle density");
  writer.var( "MSP_MINS" , "Denominator for merged spindle density" );
  
  writer.value( "MSP_N" , (int)M.size() );
  writer.value( "MSP_DENS" , M.size() / max_mins );
  writer.value( "MSP_MINS" , max_mins );
  
  //
  // frequency bins
  //
  
  const double f1 = 8;
  const double f2 = 16;
  const double fi = 0.5;
  
  std::map<double,int> bf;
  for (int i=0;i<M.size();i++)
    {
      for ( double f = f1 ; f <= f2 ; f+=fi )
	{
	  if ( M[i].frq >= f && M[i].frq < f+fi ) { bf[f]++; break; }
	}
    }

  for ( double f = f1 ; f <= f2 ; f+=fi )
    {
      writer.level( f , globals::freq_strat ); 
      writer.value( "MSP_FDENS" , bf[f] / max_mins );
    }
  writer.unlevel( globals::freq_strat ); 


  //
  // List by mspindle
  //

 
  clocktime_t starttime( edf->header.starttime );
  if ( ! starttime.valid ) 
    {
      logger << " ** could not find valid start-time in EDF header **\n";
      hms = false;
    }

  for (int i=0; i<M.size(); i++) 
    {
      
//       std::cout << "MSPINDLE " << i+1 << "\n";
//       const int nn = M[i].n();
//       for (int j = 0 ; j < nn ; j++ ) 
// 	std::cout << "\t" << M[i].spindles[j]->tp.start << " .. " << M[i].spindles[j]->tp.stop << "\n";
//       std::cout << "\n";

      writer.level( i+1 , "MSPINDLE" );
      writer.value( "MSP_F" , M[i].frq );
      writer.value( "MSP_SIZE" , M[i].n() );
      writer.value( "MSP_FL" , M[i].lwr_frq );
      writer.value( "MSP_FU" , M[i].upr_frq );
      writer.value( "MSP_DUR" , M[i].dur() );
      writer.value( "MSP_STAT" , M[i].stat );      
      writer.value( "MSP_START" , M[i].start * globals::tp_duration );
      writer.value( "MSP_STOP" , M[i].stop * globals::tp_duration );

      if ( hms )
	 {

	   double tp1_sec =  M[i].start / (double)globals::tp_1sec;
	   clocktime_t present1 = starttime;
	   present1.advance_seconds( tp1_sec );
	   // add down to 1/100th of a second
	   double tp1_extra = tp1_sec - (long)tp1_sec;
	   
	   double tp2_sec =  M[i].stop / (double)globals::tp_1sec;
	   clocktime_t present2 = starttime;
	   present2.advance_seconds( tp2_sec );
	   double tp2_extra = tp2_sec - (long)tp2_sec;
	   
	   writer.value( "MSP_START_HMS"  , present1.as_string() +  Helper::dbl2str_fixed( tp1_extra , 4  ).substr(1) );
	   writer.value( "MSP_STOP_HMS"   , present2.as_string() +  Helper::dbl2str_fixed( tp2_extra , 4  ).substr(1) );
	   
	 }

      
      if ( per_spindle_verbosity ) 
	{
	  // number of spindles in this group
	  const int n = M[i].n();

	  for (int s = 0 ; s < n ; s++ )
	    {
	      writer.level( s+1     , "SPINDLE" );
	      writer.value( "SCH"   , M[i].lab[s] );
	      writer.value( "START" , M[i].spindles[s]->tp.start * globals::tp_duration );
	      writer.value( "STOP"  , M[i].spindles[s]->tp.stop * globals::tp_duration );             
	      writer.value( "FFT"   , M[i].spindles[s]->fft ); 
	    }	  
	  writer.unlevel( "SPINDLE" );
	} // end of verbose output
           
      writer.unlevel( "MSPINDLE" );
    
    }

}


void mspindles_t::plot( const std::string & fname )
{	  
  
  //
  // Create plots for the merged spindle set?
  //
  
  // logger << " writing PDF of f-merged spindle traces to " << fname << "\n";	      
  
  // std::vector<spindle_t> mspindles_tp; // need to make a time-points
  
  // for (int i=0;i<M.size();i++)
  //   {
  //     spindle_t s( (*tp)[ mspindles[i].start ] , 
  // 		   (*tp)[ mspindles[i].stop ]  , 
  // 		   mspindles[i].start , 
  // 		   mspindles[i].stop ) ;
      
  //     mspindles_tp.push_back( s );
  //     // avg_map is already made above
  //   }	      
  
  // draw_mspindles( edf , param , fname , signals(s) , mspindles ) ; 
  
}



void mspindles_t::proc_overlaps( const std::vector<sort_t> & overlaps )
{
  
  // here we have a set of spindles with at least some overlap (based on interval_th)
  // need to determine which sets these go into, based on 
  //  - frequency differences (Fc)
  //  - interval for same Fc, difference channel 
  //  - interval for different Fc,
  
  const int ns = overlaps.size();

  //  std::cout << "ns size = " << ns << "\n";

  // // report

//   for (int i=0;i<ns;i++)
//     std::cout << "\t" << i << "\t" 
//    	      << overlaps[i].ch << " - " 
//    	      << overlaps[i].f << "\t" 
//    	      << overlaps[i].spindle->fft << "\t"
//    	      << overlaps[i].spindle->mean_stat << "\t"
//    	      << overlaps[i].i.as_string() << "\n";

  
  // take first unassigned spindle, assign to first mspindle set
  // consider all other spindles

  std::vector<std::vector<bool> > match( ns );
  for (int i=0;i<ns;i++) match[i].resize( ns , false );

  for (int i=0;i<ns;i++)
    {
      const sort_t & a = overlaps[i];

      for (int j=i+1;j<ns;j++)
      {

	const sort_t & b = overlaps[j];
	
	// base frequency comparison on estimated (FFT) frequency, rather than the 
	// target frequency

	bool same_frq = fabs( a.spindle->fft - b.spindle->fft ) <= frq_th ;
	bool same_ch  = a.ch == b.ch ; 
	bool overlap = a.i.overlaps( b.i ) ; 

	if ( overlap )
	  {
	    double o = a.i.prop_overlap( b.i );
	    
	    if ( ! same_frq ) overlap = false;
	    if ( same_ch ) 
	      { 
		if      ( o < within_ch_interval_th ) overlap = false;
		else if ( o < cross_ch_interval_th ) overlap = false;
	      }
	  }
	// set
	match[i][j] = match[j][i] = overlap;	
      }
    }
  
//     std::cout << "MATCH\n";
//     for (int i=0;i<ns;i++)
//       {
//         for (int j=0;j<ns;j++) std::cout << match[i][j] ; 
//         std::cout << "\n";
//       }


  //
  // We've created the adjacency matrix form, now find groups
  //

  // map each spindle to a mspindle set  
  std::vector<int> smap( ns , -1 ); 
  
  // assign first spindle to first mspindle group
  smap[0] = 0;
  
  int cnt = 1;

  for (int i=0;i<ns;i++)
    {
      
      if ( smap[i] == -1 ) smap[i] = cnt++;
      
      for(int j=0; j<ns; j++) 
	if ( match[i][j] ) 
	  {
	    if ( smap[j] == -1 ) smap[j] = smap[i];
	    else if ( smap[j] != smap[i] )
	      {
		int nsmap = smap[i] < smap[j] ? smap[i] : smap[j];
		for (int k=0;k<ns;k++) 
		  if ( smap[k] == smap[i] || smap[k] == smap[j] ) smap[k] = nsmap; 
	      }
	  }         
    }
  
  std::set<int> grps;
  for (int i=0;i<ns;i++) grps.insert( smap[i] );
  std::set<int>::const_iterator gg = grps.begin();
  while ( gg != grps.end() )
    {
      mspindle_t m;

      for (int i=0;i<ns;i++) 
	if ( smap[i] == *gg ) m.add( overlaps[i].spindle , overlaps[i].run , overlaps[i].label ); 
      
      // populate internal summary measures
      m.summarize();

      // add to list
      M.push_back( m ) ;
      
      ++gg;
    }

}


void mspindles_t::pairwise_statistics( int i , int j )
{
  
  // two sets, 'a' and 'b' of spindles, calculate intersection, etc
  // no merging per se
  
  std::set<interval_t> sa, sb, ba, bb, cons, uns, oa, ob;

  std::vector<spindle_t>::const_iterator ia = S[i].begin();
  while ( ia != S[i].end() )
    {
      sa.insert( ia->tp );
      ++ia;
    }

  std::vector<spindle_t>::const_iterator ib = S[j].begin();
  while ( ib != S[j].end() )
    {
      sb.insert( ib->tp );
      ++ib;
    }

  uint64_t win_msec = window <= 0 ? 0 : window * globals::tp_1sec;
  
  int olap = interval_t::intersect( sa, sb , &ba, &bb, &cons, &uns, &oa, &ob , interval_th , win_msec );
  
  double olapa = ba.size() / (double)S[i].size() ;

  double olapb = bb.size() / (double)S[j].size() ;

  // const std::string p1 = "SP_" + ch[i] + "_" + frq[i] + "_" + run_label[i] ; 
  // const std::string p2 = "SP_" + ch[j] + "_" + frq[j] + "_" + run_label[j] ; 

  const std::string p1 = "SP_" + run_label[i] ; 
  const std::string p2 = "SP_" + run_label[j] ; 

  writer.level( p1 + "x" + p2 , "PAIR" );

  writer.value( "OLAP" , ( olapa + olapb ) / 2.0 );
  writer.value( "A_IN_B" , olapa );
  writer.value( "B_IN_A" , olapb );
  
  
  //
  // Now list each spindle and flag whether or not it was consensus
  //
  
  if ( 0 ) 
    {      
      
      // std::vector<spindle_t>::const_iterator ii = a.begin();
      // while ( ii != a.end() )
      // 	{
      // 	  bool olap = oa.find( ii->tp ) == oa.end();
      // 	  std::cout << "SPINDLE_INTERSECTION" << "\t"
      // 		    << edf.id << "\t"
      // 		    << "[" << globals::current_tag << "]\t"
      // 		    << la << "\t"
      // 		    << olap << "\t"
      // 		    << ii->amp << "\t"
      // 		    << ii->dur << "\n";
      // 	  ++ii;
      // 	}
      
      // ii = b.begin();
      // while ( ii != b.end() )
      // 	{
      // 	  bool olap = ob.find( ii->tp ) == ob.end();
      // 	  std::cout << "SPINDLE_INTERSECTION" << "\t"
      // 		    << edf.id << "\t"
      // 		    << "[" << globals::current_tag << "]\t"
      // 		    << lb << "\t"
      // 		    << olap << "\t"
      // 		    << ii->amp << "\t"
      // 		    << ii->dur  << "\n";
      // 	  ++ii;
      // 	}
      
    }
  

  writer.unlevel( "PAIR" );


}
			 

void FFT_INTER()
{

  // //
  // // Spindle overlaps
  // //
  
  // if ( fft_intersections )
  //   {

  //     // get list of all types
  //     std::set<std::string> types;
      
  //     std::set<feature_t>::const_iterator ss = all_spindles.begin();
  //     while ( ss != all_spindles.end() )
  // 	{
  // 	  std::string t = ss->signal + "_" + ss->label ;
  // 	  types.insert( t );
  // 	  ++ss;
  // 	}
  
  //     // go over each pair of lists
  //     std::set<std::string>::const_iterator ii = types.begin();
  //     while ( ii != types.end() )
  // 	{
  // 	  std::set<std::string>::const_iterator jj = ii;
  // 	  while ( jj != types.end() )
  // 	    {
  // 	      if ( ii == jj ) { ++jj; continue; }
  // 	      // pick lists
  // 	      std::set<interval_t> sa, sb;
  // 	      std::set<feature_t>::const_iterator ss = all_spindles.begin();
  // 	      while ( ss != all_spindles.end() )
  // 		{
  // 		  if      ( ss->signal + "_" + ss->label == *ii ) sa.insert( ss->feature );
  // 		  else if ( ss->signal + "_" + ss->label == *jj ) sb.insert( ss->feature );
  // 		  ++ss;
  // 		}
	      
  // 	      // get intersection

  // 	      std::set<interval_t> botha, bothb;
  // 	      std::set<interval_t> cons, uns, onlya, onlyb;
	      
  // 	      double th = 0.5;
  // 	      double win = 1.0;
  // 	      uint64_t win_msec = win <= 0 ? 0 : win * globals::tp_1sec;

  // 	      int olap = interval_t::intersect( sa, sb, &botha, &bothb, &cons, &uns, &onlya, &onlyb, th, win_msec );
	        
  // 	      double olap1 = 2 * ( olap / (double)(sa.size() + sb.size() ) );
	      
  // 	      double olapa = botha.size() / (double)sa.size() ;
  // 	      double olapb = bothb.size() / (double)sb.size() ;
	      
  // 	      // average
  // 	      double olap2 = ( olapa + olapb ) / 2.0;
	      
	      
  // 	      // output  
	      
  // 	      writer.value( "PAIRED_PAIR" , *ii + "_" + *jj  ); 

  // 	      writer.value( "PAIRED_OLAP" , olap );
  // 	      writer.value( "PAIRED_SA" , sa.size() );
  // 	      writer.value( "PAIRED_SB" , sb.size() );	      
  // 	      writer.value( "PAIRED_OLAP1" , olap1 );
  // 	      writer.value( "PAIRED_OLAP2" , olap2 );
	      
  // 	      // next pair
  // 	      ++jj;
  // 	    }
  // 	  ++ii;
  // 	}
      
  //   }
  
}




void OLD_FFT_INTERSECTION()
{
  //
  // this code for reference more than anything else...
  //
  
  
  // was linked to 'fft-intersections' option
  

  // //
  // // 1. Merge all spindles
  // //
  
  // std::vector<interval_t> all_unique_spindles;
  
  // // and track inserts
  // std::vector<std::set<std::string> > flags;
  
  // std::set<feature_t> overlaps; // tmp
  // uint64_t last_sp;
  // uint64_t first_sp;
  
  // //      logger << "DET: have " << all_spindles.size() << " initial total spindles\n";

  // std::set<feature_t>::const_iterator ss = all_spindles.begin();
  // while ( ss != all_spindles.end() )
  // 	{
	  
  // 	  // std::cout << "considering " << ss->f << " : " << ss->i.start << " - " << ss->i.stop << "\n";
	  
  // 	  // if empty, start an overlap bag
  // 	  if ( overlaps.empty() ) 
  // 	    {
  // 	      overlaps.insert( *ss );
  // 	      first_sp = ss->feature.start;
  // 	      last_sp = ss->feature.stop;	      
  // 	      ++ss; // next spindle
  // 	      continue;
  // 	    }
	  
  // 	  // first check for /any/ overlap
  // 	  bool has_overlap = ss->feature.start <= last_sp;
	  
  // 	  // then refine
  // 	  if ( has_overlap )
  // 	    {
  // 	      int max_start = first_sp;
  // 	      int min_start = ss->feature.start;  // i.e. because sorted
	      
  // 	      int min_stop = last_sp < ss->feature.stop ? last_sp : ss->feature.stop;
  // 	      int max_stop = last_sp < ss->feature.stop ? ss->feature.stop : last_sp ; 
	      
  // 	      int o_intersection = min_stop - min_start + 1;
  // 	      int o_union        = max_stop - max_start + 1;
	      
  // 	      double metric = o_intersection / (double)o_union;
	      
  // 	      has_overlap = metric > fft_overlap_threshold;
	      
  // 	      // or if one spindle is more or less completely subsumed by the other(s)
  // 	      if ( ( last_sp - first_sp + 1 ) / o_intersection >= 0.8 ) has_overlap = true;
  // 	      if ( ( ss->feature.stop - ss->feature.start + 1 ) / o_intersection >= 0.8 ) has_overlap = true;		  
	      
  // 	    }
	  
  // 	  // if we have overlap as defined
  // 	  if ( has_overlap )
  // 	    {		  
  // 	      overlaps.insert( *ss );
  // 	      // extend range
  // 	      if ( ss->feature.stop > last_sp ) last_sp = ss->feature.stop;		 
  // 	    }
	  
	  
  // 	  //
  // 	  // Now merge
  // 	  //
	  
  // 	  else
  // 	    {

  // 	      // merged spindle goes from first_sp to last_sp;  
  // 	      interval_t mspindle( first_sp , last_sp );
	      
  // 	      // save interval...
  // 	      all_unique_spindles.push_back( mspindle );
	      
  // 	      // ...and info on contributing signals/frqs
  // 	      std::set<std::string> f;
  // 	      std::set<feature_t>::const_iterator oo = overlaps.begin();
  // 	      while ( oo != overlaps.end() )
  // 		{
  // 		  f.insert( oo->signal + "_" + oo->label );
  // 		  ++oo;
  // 		}
  // 	      flags.push_back(f);

  // 	      // clean up
  // 	      overlaps.clear();
	      
  // 	      // and add the next spindle (i.e. current one)
  // 	      overlaps.insert( *ss );
  // 	      first_sp = ss->feature.start;
  // 	      last_sp = ss->feature.stop;	      
  // 	    }
	  
  // 	  ++ss;
  // 	}
      
      
  //     //
  //     // empty the last bag here too..
  //     //
      
  //     // merged spindle goes from first_sp to last_sp;  
  //     interval_t mspindle( first_sp , last_sp );
  //     all_unique_spindles.push_back( mspindle );
  //     // ...and info on contributing signals/frqs
  //     std::set<std::string> f;
  //     std::set<feature_t>::const_iterator oo = overlaps.begin();
  //     while ( oo != overlaps.end() )
  // 	{
  // 	  f.insert( oo->signal + "_" + oo->label );
  // 	  ++oo;
  // 	}
  //     flags.push_back(f);
            
  //     //
  //     // Dump each
  //     //
      
  //     logger << " found " << all_unique_spindles.size() << " " << flags.size() << " uniq spins\n";

  //     std::set<std::string> all_flags;
  //     for (int i=0;i<flags.size();i++)
  // 	{
  // 	  std::set<std::string>::const_iterator ff=flags[i].begin();
  // 	  while ( ff != flags[i].end() )
  // 	    {
  // 	      all_flags.insert( *ff );
  // 	      ++ff;
  // 	    }
  // 	}

  //     //
  //     // header
  //     //

  //     std::cout << "MSPINDLE" << "\t" << "FLAGS" << "\t"
  // 		<< edf.id << "\t"
  // 		<< "[" << globals::current_tag << "]";
      
  //     std::set<std::string>::const_iterator ff = all_flags.begin();
  //     while ( ff != all_flags.end() )
  // 	{
  // 	  std::cout << "\t" << *ff;
  // 	  ++ff;
  // 	}
  //     std::cout << "\n";
      
      
  //     // also, for this individual, report flags for FFT
  //     // i.e. use this to check that all people line up...

  //     const int bl = 10; const int bu = 16;
  //     mslice_t mslice( edf , signals , interval_t(0,1) );

  //     std::string fft_sigs = "";
  //     for (int s=0;s<ns;s++) 
  // 	for (int b=bl;b<=bu;b++) 
  // 	  fft_sigs += "_" + mslice.label(s) + ":" + Helper::int2str( b );
      
  //     std::cout << "MSPINDLE" << "\t"
  // 		<< "FFT-FLAGS" << "\t"
  // 		<< edf.id << "\t"
  // 		<< "[" << globals::current_tag << "]\t"
  // 		<< fft_sigs << "\n";
      
      
      
  //     //
  //     // report for each merged spindle
  //     //

  //     for (int i=0;i<all_unique_spindles.size();i++)
  // 	{
  // 	  const interval_t spindle_interval = all_unique_spindles[i] ;
  // 	  double dur_sec = spindle_interval.duration() /(double)globals::tp_1sec; 
	  
  // 	  std::cout << "MSPINDLE" << "\t"
  // 		    << "FLAGGED" << "\t"
  // 		    << edf.id << "\t"
  // 		    << "[" << globals::current_tag << "]\t"
  // 		    << i+1 << "\t"
  // 		    << spindle_interval.start_sec() << "\t"
  // 		    << spindle_interval.stop_sec() << "\t"
  // 		    << dur_sec;
	  
  // 	  std::set<std::string>::const_iterator ff = all_flags.begin();
  // 	  while ( ff != all_flags.end() )
  // 	    {
  // 	      if ( flags[i].find( *ff ) != flags[i].end() ) std::cout << "\t1"; else std::cout << "\t0";
  // 	      ++ff;
  // 	    }

  // 	  std::cout << "\n";
	  

  // 	  //
  // 	  // perform FFT on segment
  // 	  //
	  
  // 	  mslice_t mslice( edf , signals , spindle_interval );
	  


  // 	  std::cout << "MSPINDLE" << "\t"
  // 		    << "FFT" << "\t"
  // 		    << edf.id << "\t"
  // 		    << "[" << globals::current_tag << "]" << "\t"
  // 		    << i+1;
	    
  // 	  // consider each signal
  // 	  for (int s=0;s<ns;s++)
  // 	    {
	  
  // 	      const int Fs = edf.header.sampling_freq( signals(0) );

  // 	      const std::vector<double> * data = mslice.channel[s]->pdata();
	      
  // 	      // use fixed 1024 points... note, need to change if v. different Fs...
  // 	      //const int npoints = data->size();
  // 	      const int NF = 1024;
  // 	      const int NF2 = NF < data->size() ? NF : data->size() ; 
	      
  // 	      std::vector<double> padded( NF , 0 );
  // 	      for(int j=0;j<NF2;j++) padded[j] = (*data)[j];
	      	      
  // 	      FFT fftseg( NF , Fs , FFT_FORWARD , WINDOW_HANN );	      
	     
  // 	      fftseg.apply( padded );

  // 	      std::map<int,double> bpow;
  // 	      double tpow = 0;
		
  // 	      int N = fftseg.cutoff;
  // 	      for (int f=0;f<N;f++)
  // 		{
  // 		  //if ( fftseg.frq[f] > (bu+1) ) break;
		  
  // 		  tpow += fftseg.X[f];
		  
  // 		  for (int b=bl;b<=bu;b++)
  // 		    {		      
  // 		      if ( fftseg.frq[f] >= b && fftseg.frq[f] < b+1 ) bpow[b] += fftseg.X[f]; 
  // 		    }
		  
  // 		}

	      
  // 	      // output (continue existing line)
  // 	      for (int b=bl;b<=bu;b++) std::cout << "\t" << log( bpow[b] / (double)tpow );   
	      
  // 	    } // next signal in this m-spindle
	  
  // 	  std::cout << "\n";
	  
  // 	} // next m-spindle


  //   } // end of fft-intersections
  
}
