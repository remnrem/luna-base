
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

#include "spindles/propag.h"

#include "db/db.h"
#include "helper/logger.h"

extern logger_t logger;
extern writer_t writer;

void sp_props_t::add_tp ( const std::vector<uint64_t> & tp ) 
{
  if ( tps.size() == 0 ) 
    tps = tp;
  else if ( tps.size() != tp.size() )
    Helper::halt( "internal error in prop(): must be similar intervals/sampling rates across signals" );
}

void sp_props_t::add( double f , const std::string & ch , const std::vector<spindle_t> & sp, const std::vector<double> & cwt )
{
  // store frequency as integer
  sp_idx_t idx( (uint64_t)( 1e9 * f ) , ch );
  sp_dat_t dat( sp , cwt );
  data[ idx ] = dat;
}

double sp_props_t::analyse( const std::set<double> & f , 
			    const std::set<std::string> & c , 
			    const std::string & seed , 			  
			    const double w ,
			    const bool verbose )
{

  
  // frequencies in data[] to consider
  // channels in data[] to consider
  // for now, assume a single frequency only
  
  if ( f.size() != 1 ) 
    Helper::halt( "expecting a single frequency for prop()" );

  uint64_t frq = *f.begin() * 1e9;
  
  //logger << "  spindle propagation analysis, seeding on " << seed << "\n";

  writer.level( seed , "SEED" );

  //
  // check we have everything as needed 
  //
  
  // TPs should match the coefficients for each individual
  const int n = tps.size();
  if ( n == 0 ) Helper::halt( "no time-points specified" );


  // originally, search within +/- w seconds
  //  but change to ignore 'w' and just search w/in the seed spindle window
    
  uint64_t w_tp = w * globals::tp_1sec;
  
  //
  // seed coefficients;
  //
  
  sp_idx_t idx( frq , seed );
  
  if ( data.find( idx ) == data.end() ) 
    Helper::halt( "could not seed on " + seed );
  
  //
  // Get map of all channels
  //
  
  std::map<std::string,double> cmap0;
  std::map<sp_idx_t,sp_dat_t>::const_iterator ii = data.begin();
  while ( ii != data.end() ) 
    {
      cmap0[ ii->first.ch ] = 0;
      ++ii;
    }
    
  //
  // for all other channels, record the nearest peak, and save the offset relative to midx
  // search within a +/- w window  (actually, now, within window of seed spindle)
  //
  
  std::map<std::string,double> cmap     = cmap0;
  std::map<std::string,double> cmap_amp = cmap0;
  std::map<std::string,double> cmap_n   = cmap0;

  // also take separate means for pre-seed and post-seed events
  std::map<std::string,double> pre_cmap     = cmap0;
  std::map<std::string,double> pre_cmap_amp = cmap0;
  std::map<std::string,double> pre_cmap_n   = cmap0;

  std::map<std::string,double> post_cmap     = cmap0;
  std::map<std::string,double> post_cmap_amp = cmap0;
  std::map<std::string,double> post_cmap_n   = cmap0;

  
  //
  // get all spindle peaks (based on max CWT) within spindle window
  //
  
  const sp_dat_t & dat = data[ idx ];
  
  const int np = dat.sp.size();

  if ( dat.coeff.size() != tps.size() )
    Helper::halt( "internal error in prop(): wrong TP/CWT size alignment" );
      
  for (int p=0; p<np; p++)
    {

      if ( verbose )
	writer.level( p+1 , "SPINDLE" );
      
      const spindle_t & sp = dat.sp[p];
      
      int midx = sp.start_sp;
      
      double mx = 0; // CWT always positive
      
      for (int i = sp.start_sp; i <= sp.stop_sp; i++)
	{
	  if ( dat.coeff[i] > mx ) 
	    {
	      mx = dat.coeff[i];
	      midx = i;
	    }
	}
      
      // midx is the peak

      std::map<std::string,double>::iterator cc = cmap.begin();
      while ( cc != cmap.end() ) 
	{
	  
	  // data for the other channel
	  const std::vector<double> & coeff = data.find( sp_idx_t( frq , cc->first ) )->second.coeff;
	  
	  // channel
	  const std::string ch = cc->first;

	  // skip self-seed
	  if ( seed == ch ) 
	    {
	      ++cc; continue;
	    }

	  if ( verbose )
	    writer.level( ch , globals::signal_strat );
	  
	  // mean for this channel:  for threshold
	  // as CWT are already baseline-adjusted when sent here
	  // so this means ignore if below the mean for the other spindle

	  // this implies 50% of spindle seed max amp 
	  const double th = 0.5 ;
	  
	  // startpoint for seach
	  
	  int nidx = midx;
	  double nx = coeff[midx];

	  if ( 0 ) // old search window +/- 1 second
	    {
	      int pidx = midx;
	      while ( 1 ) 
		{
		  if ( pidx == 0 ) break;
		  --pidx;
		  if ( tps[midx] - tps[pidx] > w_tp ) break;
		  if ( coeff[pidx] > nx )
		    {
		      nx = coeff[pidx];
		      nidx = pidx;
		    }
		}
	      
	      // other way
	      pidx = midx;
	      while ( 1 )
		{
		  ++pidx;
		  if ( pidx == np ) break;
		  if ( tps[pidx] - tps[midx] > w_tp ) break;
		  if ( coeff[pidx] > nx )
		    {
		      nx = coeff[pidx];
		      nidx = pidx;
		    }
		}
	    }

	  //
	  // new - use original seed window
	  //
	  
	  for (int i = sp.start_sp; i <= sp.stop_sp; i++)
	    {
	      if ( coeff[i] > nx )
		{
		  nx = coeff[i];
		  nidx = i;
		}	     
	    }

	  // express channel amplitude as a proportion of the seed peak amp.
	  //  and require that it is at least 50% of the detected spindle to include

	  nx /= mx;
	  	  
	  
	  // max (nearest to midx) for this channel will be 
	  
	  // all
	  // for (int i = sp.start_sp; i <= sp.stop_sp; i++)
	  //   {
	  //     std::cout << i << "\t" << dat.coeff[i] << "\t" << coeff[i] << "\n"; 
	  //   }
	  // std::cout << "  ---> seed ch " << seed << " " << ch << " " << midx << " " << nidx << "\n\n";
	  
	  
	  // only include if the other channel is above some threshold?
	  
	  if ( nx >= th )
	    {
	      
	      // this max point, relative to the index/seed spindle:
	      double offset_sec = globals::tp_duration * ( (double)tps[ nidx ] - (double)tps[ midx ] );

	      // is this channel prior to or after seed?
	      int pos = nidx < midx ? -1 : ( nidx == midx ? 0 : +1 ) ;
	      
	      if ( verbose )
		{
		  writer.value( "T" , offset_sec );
		  writer.value( "REL" , pos );
		}

	      // store overall 
	      ++cmap_n[ ch ];
	      cmap[ ch ] += offset_sec;
	      cmap_amp[ ch ] += nx ;

	      // also store specific means for pre-seed and post-seed events
	      if ( pos == -1 )
		{
		  ++pre_cmap_n[ ch ];
		  pre_cmap[ ch ] += offset_sec;
		  pre_cmap_amp[ ch ] += nx ;
		}
	      else if ( pos == +1 )
		{
		  ++post_cmap_n[ ch ];
		  post_cmap[ ch ] += offset_sec;
                  post_cmap_amp[ ch ] += nx ;
		}
	      
	    }

	  // next CH
	  ++cc;
	}
      
      if ( verbose ) 
	writer.unlevel( globals::signal_strat ); 
            
      // next spindle for the seed channel
      
    }

  if ( verbose )
    writer.unlevel( "SPINDLE" );
  

  
  //
  // get means for all cmaps
  //
  
  
  //      CH   ...  SEED   ...    CH 
  //  T        -ve   0     +ve    

  //  i.e. so get average of -1 value for the seed average
  //   that way a higher number means it occurs later (versus other channels)
  //   return this value for the seed from this function ( analyse() ) 
  //   and then take the set for all seeds (for a given spindle frequency) 
  //   and scale from 0..1 to indicate whether that channel tends to start 
  //   earlier (0.0) or later (1.0)
  //   

  std::map<std::string,double>::iterator cc = cmap.begin();
  
  double seed_avg = 0;
  int    seed_cnt = 0;

  while ( cc != cmap.end() )
    {

      if ( cc->first == seed )
	{
	  ++cc; continue;
	}
      
      writer.level( cc->first , globals::signal_strat );

      
      writer.value( "N" , cmap_n[ cc->first ] );
      writer.value( "P" , cmap_n[ cc->first ] / (double)np ); 
      
      if ( cmap_n[ cc->first ] > 0 ) 
	  {
	    
	    cmap[ cc->first ] /= cmap_n[ cc->first ];
	    cmap_amp[ cc->first ] /= cmap_n[ cc->first ];
	    
	    // mean time-offset of paired channel, relatve to seed spindle peak
	    writer.value( "T" , cmap[ cc->first ] );

	    // mean amplitude of paired-channel at seed spindle
	    writer.value( "A" , cmap_amp[ cc->first ] );	    
	    
	    // track for seed-level average (nb. take negative value
	    //  as we'll interpret this in terms of the seed being 
	    //  earlier or later ) 

	    seed_avg += -1 * cmap[ cc->first ] ; 
	    ++seed_cnt;
	  }

      //
      // pre/post stratified reports
      //

      writer.value( "N_PRESEED" , pre_cmap_n[ cc->first ] );
      writer.value( "P_PRESEED" , pre_cmap_n[ cc->first ] / (double)np );      
      
      writer.value( "N_POSTSEED" , post_cmap_n[ cc->first ] );
      writer.value( "P_POSTSEED" , post_cmap_n[ cc->first ] / (double)np );

      // scaled ( -1 to +1 ) pre/post metric: proportion of times CH is before versus after SEED
      const double pp = post_cmap_n[ cc->first ] / ( post_cmap_n[ cc->first ] + pre_cmap_n[ cc->first ] );

      writer.value( "PP" , 2 * ( pp - 0.5 ) ) ;

      if ( pre_cmap_n[ cc->first ] > 0 )
          {	    
            pre_cmap[ cc->first ] /= pre_cmap_n[ cc->first ];
            pre_cmap_amp[ cc->first ] /= pre_cmap_n[ cc->first ];
            writer.value( "T_PRESEED" , pre_cmap[ cc->first ] );
            writer.value( "A_PRESEED" , pre_cmap_amp[ cc->first ] );
          }

      if ( post_cmap_n[ cc->first ] > 0 )
	{	    
            post_cmap[ cc->first ] /= post_cmap_n[ cc->first ];
            post_cmap_amp[ cc->first ] /= post_cmap_n[ cc->first ];
            writer.value( "T_POSTSEED" , post_cmap[ cc->first ] );
            writer.value( "A_POSTSEED" , post_cmap_amp[ cc->first ] );
	}

      //
      // next channel
      //
      
      ++cc;
    }

  // report average over all other channels for this seed

  writer.unlevel( globals::signal_strat );
  writer.unlevel( "SEED" );
  
  return seed_cnt > 0 ? seed_avg / (double) seed_cnt : 0 ; 
  
}
