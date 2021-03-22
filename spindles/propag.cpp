
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

void sp_props_t::analyse( const std::set<double> & f , 
			  const std::set<std::string> & c , 
			  const std::string & seed , 			  
			  const double w )
{
  
  // frequencies in data[] to consider
  // channels in data[] to consider
  // for now, assume a single frequency only
  
  if ( f.size() != 1 ) 
    Helper::halt( "expecting a single frequency for prop()" );

  uint64_t frq = *f.begin() * 1e9;
  
  //
  
  logger << "  spindle propagation analysis, seeding on " << seed << "\n";
  
  // check we have everything as needed 
  
  // TPs should match the coefficients for each individual
  const int n = tps.size();
  if ( n == 0 ) Helper::halt( "no time-points specified" );
  
  uint64_t w_tp = w * globals::tp_1sec;
  
  // seed coefficients;
  
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
  // search within a +/- w window
  //
  
  std::map<std::string,double> cmap     = cmap0;
  std::map<std::string,double> cmap_amp = cmap0;
  std::map<std::string,double> cmap_n   = cmap0;
  

  //
  // get all spindle peaks (based on max CWT) within spindle window
  //
  
  const sp_dat_t & dat = data[ idx ];
  
  const int np = dat.sp.size();

  if ( dat.coeff.size() != tps.size() )
    Helper::halt( "internal error in prop(): wrong TP/CWT size alignment" );
      
  for (int p=0; p<np; p++)
    {
      
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
	  
	  // mean for this channel:  for threshold
	  // as CWT are already baseline-adjusted when sent here
	  // so this means ignore if below the mean for the other spindle
	  
	  const double th = 1.0 ;

	  // startpoint for seach
	  
	  int nidx = midx;
	  double nx = coeff[midx];
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
	  
	  // max (nearest to midx) for this channel will be 
	  
	  // only include if the other channel is above some threshold?
	  
	  if ( nx >= th )
	    {
	      
	      // this max point, relative to the index/seed spindle:
	      double offset_sec = globals::tp_duration * ( (double)tps[ nidx ] - (double)tps[ midx ] );
	      
	      ++cmap_n[ ch ];
	      cmap[ ch ] += offset_sec;
	      cmap_amp[ ch ] += nx ; 
	    }
	  
	  ++cc;
	}
      
      
      // next spindle for the seed channel
      
    }
  
  
  //
  // get means for all cmaps
  //
  
  writer.level( seed , "SEED" );
  writer.level( *f.begin() , globals::freq_strat );
  
  std::map<std::string,double>::iterator cc = cmap.begin();
  
  while ( cc != cmap.end() )
    {

      writer.level( cc->first , globals::signal_strat );
      writer.value( "N" , cmap_n[ cc->first ] );
	
	if ( cmap_n[ cc->first ] > 0 ) 
	  {
	  cmap[ cc->first ] /= cmap_n[ cc->first ];
	  cmap_amp[ cc->first ] /= cmap_n[ cc->first ];
	  writer.value( "T" , cmap[ cc->first ] );
	  writer.value( "P" , cmap_amp[ cc->first ] );	    
	}
      
      ++cc;
    }
  
  writer.unlevel( globals::signal_strat );
  writer.unlevel( globals::freq_strat );
  writer.unlevel( "SEED" );

}
