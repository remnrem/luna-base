
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


#include "correl.h"

#include "miscmath/miscmath.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "eval.h"

#include "stats/statistics.h"

#include "db/db.h"
#include "dsp/resample.h"

#include "helper/helper.h"
#include "helper/logger.h"

#include <string>


extern logger_t logger;

extern writer_t writer;
   
void dsptools::correlate_channels( edf_t & edf , param_t & param )
{
  
  std::string signal_label1 , signal_label2 ; 

  bool all_by_all = false;

  if ( param.has( "sig1" ) )
    {
      signal_label1 =  param.requires( "sig1" );
      signal_label2 =  param.requires( "sig2" );
    }
  else
    {
      all_by_all = true;
      signal_label1 = signal_label2 = param.requires( "sig" );
    }
  
  
  signal_list_t signals1 = edf.header.signal_list( signal_label1 );  

  signal_list_t signals2 = edf.header.signal_list( signal_label2 );  
  
  const int ns1 = signals1.size();
  const int ns2 = signals2.size();
  
  int sr = param.has( "sr" ) ? param.requires_int( "sr" ) : 0 ;

  
  writer.var( "R" , "Channel correlation (-1..1)" );


  //
  // compile signals
  //

  std::vector<int> sigs;
  std::map<int,std::string> sigset;

  for (int s=0;s<ns1;s++) 
    {
      if ( sigset.find( signals1(s) ) == sigset.end() )
	{
	  sigset[ signals1(s) ] =  signals1.label(s);
	  sigs.push_back( signals1(s) );
	}
    }


  for (int s=0;s<ns2;s++) 
    {
      if ( sigset.find( signals2(s) ) == sigset.end() )
	{
	  sigset[ signals2(s) ] =  signals2.label(s);
	  sigs.push_back( signals2(s) );
	}
    }

  const int ns = sigs.size();


  //
  // adjust all SRs now if needed
  //

  if ( sr )
    {
      for (int s=0;s<ns;s++)
	{
	  
	  if ( edf.header.is_annotation_channel( sigs[s] ) ) continue;
	  
	  if ( edf.header.sampling_freq( sigs[s] ) != sr ) 
	    {
	      // logger << "resampling channel " << sigset[ sigs[s] ]
	      // 	     << " from " << edf.header.sampling_freq( sigs[s] )
	      // 	     << " to " << sr << "\n";
	      resample_channel( edf, sigs[s] , sr );
	    }
	}
    }
  else
    {
      // check SR
      std::set<int> srs;
      for (int s=0;s<ns;s++)
	{
	  if ( edf.header.is_annotation_channel( sigs[s]) ) continue;
	  srs.insert( edf.header.sampling_freq( sigs[s] ) );
	}
      if ( srs.size() > 1 ) Helper::halt( "all sampling rates must be similar, use 'sr'" );
    }

  
  //
  // Epochs or whole signal?
  //
  
  bool epoched = param.has( "epoch" ) ;


  //
  // Number of pairwise comparisons
  //

  const int np = all_by_all ? ns*(ns-1) : ns1 * ns2;
  

  //
  // For channel-level summaries (not in per-epoch mode)
  //

  bool ch_summaries = param.has( "ch-high" ) || param.has( "ch-low" );
  double ch_over = param.has( "ch-high" ) ? param.requires_dbl( "ch-high" ) : 1;
  double ch_under = param.has( "ch-low" ) ? param.requires_dbl( "ch-low" ) : -1;
  std::map<int,double> summr_mean, summr_min, summr_max;
  std::map<int,int> summr_under, summr_over, summr_n;
  for (int s=0;s<ns;s++) { summr_min[ sigs[s] ] = 1; summr_max[ sigs[s] ] = -1; }
  if ( ch_summaries && ! all_by_all ) Helper::halt( "can only do ch-high/ch-low summaries with all-by-all CORREL (i.e. sig=X, not sig1=X sig2=Y)" );
  
  //
  // Start iterating over pairs
  //

  logger << "  calculating correlation for " << np << " channel pairs\n";

  //
  // Iterate over pairs
  //
  
  for (int i=0;i<ns1;i++)
    {
      
      if ( edf.header.is_annotation_channel( signals1(i) ) ) continue;

      writer.level( signals1.label(i) , "CH1" );
      
      for (int j=0;j<ns2;j++)
	{
	  
	  if ( edf.header.is_annotation_channel( signals2(j) ) ) continue;
	  
	  //
	  // if all-by-all, do not do redundant channels
	  //
	  
	  if ( all_by_all ) 
	    {
	      if ( j <= i ) continue;
	    }

	  //
	  // Store correlations
	  //
	  
	  std::vector<double> epoch_r;
	  double overall_r;
	  double mean_r;
	  double median_r;
	  
	  // stratify output by SIGNALS
	  //writer.level( signals1.label(i) + "x" + signals2.label(j) , "CHS" );
	  
	  writer.level( signals2.label(j) , "CH2" );
	  
	  if ( epoched ) 
	    {
	      
	      int epoch = edf.timeline.first_epoch();      
	      
	      while ( 1 ) 
		{
		  
		  int epoch = edf.timeline.next_epoch();      
		  if ( epoch == -1 ) break;

		  interval_t interval = edf.timeline.epoch( epoch );

 		  slice_t slice1( edf , signals1(i) , interval );
 		  slice_t slice2( edf , signals2(j) , interval );
		  
 		  const std::vector<double> * d1 = slice1.pdata();
 		  const std::vector<double> * d2 = slice2.pdata();

		  double r = Statistics::correlation( *d1 , *d2 );

		  //
		  // Output
		  //
		  
		  writer.epoch( edf.timeline.display_epoch( epoch ) );
		  
		  // -9 is NA code for a correlation
		  if ( r > -5 ) 
		    {
		      writer.value( "R" , r );
		      epoch_r.push_back( r );
		    }
		  
		} // next epoch
	      
	      
	      writer.unepoch();
	      
	      //
	      // Get mean/median correlation over epochs
	      //
	      
	      if ( epoch_r.size() > 0 ) 
		{
		  mean_r = MiscMath::mean( epoch_r );
		  
		  median_r = MiscMath::median( epoch_r );
		  
		  writer.value( "R_MEAN" , mean_r ); 
		  
		  writer.value( "R_MEDIAN" , median_r ); 
		}

	    }
	  else
	    {
	      
	      //
	      // Coherence for entire signal
	      //

	      interval_t interval = edf.timeline.wholetrace();
	      
	      slice_t slice1( edf , signals1(i) , interval );
	      slice_t slice2( edf , signals2(j) , interval );
	      
	      const std::vector<double> * d1 = slice1.pdata();
	      const std::vector<double> * d2 = slice2.pdata();
	      
	      
	      overall_r = Statistics::correlation( *d1 , *d2 );
	      
	      if ( overall_r > -5 ) 
		{
		  
		  writer.value( "R" , overall_r ); 
	      
		  // store channel-level summaries
		  if ( ch_summaries )
		    {
		      summr_mean[ signals1(i) ] += overall_r;
		      summr_mean[ signals2(j) ] += overall_r;
		      ++summr_n[ signals1(i) ];
		      ++summr_n[ signals2(j) ];
		      
		      if ( overall_r > summr_max[ signals1(i) ] ) summr_max[ signals1(i) ] = overall_r;
		      if ( overall_r > summr_max[ signals2(j) ] ) summr_max[ signals2(j) ] = overall_r;
		      
		      if ( overall_r < summr_min[ signals1(i) ] ) summr_min[ signals1(i) ] = overall_r;
		      if ( overall_r < summr_min[ signals2(j) ] ) summr_min[ signals2(j) ] = overall_r;
		      
		      if ( overall_r > ch_over ) { ++summr_over[ signals1(i) ]; ++summr_over[ signals2(j) ]; }
		      if ( overall_r < ch_under ) { ++summr_under[ signals1(i) ]; ++summr_under[ signals2(j) ]; }
		      
		    }	
		}
	    }	 
	}
    }
  
  writer.unlevel( "CH1" );
  writer.unlevel( "CH2" );

  
  //
  // Write channel-level summaries (min/max/mean/median/# channels over t1/# channels under t2)  correlations, if whole-signal (only)
  //

  if ( ch_summaries )
    {

      logger << "  writing channel-level summaries: mean min max";
      if ( ch_over < 1 ) logger << " #>" << ch_over ;
      if ( ch_under > -1 ) logger << " #<" << ch_under ;
      logger << "\n";
      
      for (int s=0;s<ns;s++)
	{
          if ( edf.header.is_annotation_channel( sigs[s]) ) continue;
	  writer.level( sigset[ sigs[s] ] , globals::signal_strat );
	  writer.value( "SUMM_MEAN" , summr_mean[ sigs[s] ] / (double)summr_n[ sigs[s] ] );
	  writer.value( "SUMM_MIN" , summr_min[ sigs[s] ] );
	  writer.value( "SUMM_MAX" , summr_max[ sigs[s] ] );
	  writer.value( "SUMM_HIGH" , summr_over[ sigs[s] ] );
	  writer.value( "SUMM_LOW" , summr_under[ sigs[s] ] );
	}
      writer.unlevel( globals::signal_strat );
      
    }
  
  
}

