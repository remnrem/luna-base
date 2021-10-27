
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

#include "artifacts/correct.h"

#include "helper/helper.h"

#include "dsp/emd.h"
#include "db/db.h"
#include "helper/logger.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "stats/glm.h"
#include <vector>

extern logger_t logger;
extern writer_t writer;

void dsptools::artifact_correction( edf_t & edf , param_t & param )
{

  //
  // segment-wise artifact correction
  //

  const bool no_annotations = true;

  // signals to-be-corrected
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , no_annotations );

  // template signals
  signal_list_t correctors = edf.header.signal_list( param.requires( "corr" ) );

  const int ns = signals.size();
  const int nc = correctors.size();
  
  //
  // EMD based model
  //

  const int emd_mode = param.has( "emd" ) ? param.requires_int( "emd" ) : 0; 
  const double emd_th = param.has( "th" ) ? param.requires_dbl( "th" ) : 0.9;
  const bool emd_corr = param.has( "emd-corr" ) ; 
  
  const int max_sift = 20;
  const int max_imf = 10;

  //
  // Regression-based model
  //
  
  const bool regression_mode = ! emd_mode ; 

  
  // regression mode: fix to 50% overlap
  // EMD mode: no overlap (e.g. expect 30s epochs)

  const double segment_size_sec = param.has( "segment-sec" )
    ? param.requires_dbl( "segment-sec" )
    : ( regression_mode ? 5 : 30 );
    
  const double segment_step_sec = regression_mode
    ? segment_size_sec / 2.0
    : segment_size_sec ;  
  
  const int ne = edf.timeline.first_epoch();

  //
  // requires similar SRs
  //

  std::vector<double> Fs  = edf.header.sampling_freq( signals );
  std::vector<double> FsC = edf.header.sampling_freq( correctors );

  int sr = Fs[0];
  for (int s=1;s<ns;s++)
    if ( Fs[s] != sr )
      Helper::halt( "all sampling rates must be similar for CORRECT" );
  for (int c=0;c<nc;c++)
    if ( FsC[c] != sr )
      Helper::halt( "all sampling rates must be similar for CORRECT" );



  //
  // Misc output 
  //

  
  logger << "  applying " << ( regression_mode ? "regression" : "EMD-based" )
	 << " correction for " << ns << " signals based on " << nc << " reference signals\n";
  
  logger << "  using a segment size of " << segment_size_sec << " seconds;";
  if ( regression_mode ) logger << " with 50% overlap";
  logger << "\n";

  
      
  //
  // Iterate over signals
  //

  for (int s=0; s<ns; s++)
    {


      //
      // Get data : signal
      //

      interval_t interval = edf.timeline.wholetrace();
      
      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();

      //
      // Get data : correctors 
      //

      matslice_t mslice( edf , correctors , interval );
      
      const Data::Matrix<double> & Z = mslice.data_ref();

      
      //
      // Process segment-wise
      //

      const int total_points = d->size();

      const int segment_points = segment_size_sec * sr;

      const int step_points  = segment_step_sec * sr;
           
      //
      // Iterate over segments
      //

      std::map<int,Data::Vector<double> > res;

      for (int p = 0; p < total_points ; p += step_points )
	{

	  // all done?

	  if ( p + segment_points > total_points ) break;

	  //
	  // actual signal
	  //

	  Data::Vector<double> y( segment_points );	  
	  for (int i=0; i<segment_points ; i++) y[i] = (*d)[ p + i ];

	  //
	  // correctors
	  //
	  
	  Data::Matrix<double> zz( segment_points , nc );
	  for (int i=0; i<segment_points ; i++)
	    for (int j=0; j<nc; j++)
	      zz( i , j ) = Z( p + i , j );
	  
	  
	  //	  
	  // regression-based correction?
	  //

	  if ( regression_mode )
	    {

	      // mean center
	      const double ymean = Statistics::mean( y );
	      for (int i=0; i<segment_points ; i++) y[i] -= ymean;

	      GLM glm( GLM::LINEAR );
	      
	      glm.set( y , zz );
	      
	      const bool okay = glm.fit();
	      
	      if ( ! okay )
		res[ p ] = y;
	      else
		res[ p ] = glm.get_residuals();
	    }

	  //
	  // Or EMD based correction?
	  //

	  if ( emd_mode )
	    {
	      emd_t emd;
	      emd.max_sift = max_sift;
	      emd.max_imf = max_imf;
	      const int nimf = emd.proc( y.data_pointer() );
	      
	      // do any segments highly correlate w/ any correctors?
	      //  OR.. any EMD of any correctors?
	      // if so, remove them; nb. 1-based encoding used for IMF, as residual == 0 
	      
	      std::set<int> remove;

	      //
	      // Iterate over correctors
	      //

	      for (int c=0; c<nc; c++)
		{

		  // use corrector 'as is', for do EMD of corrector?

		  if ( emd_corr )
		    {

		      //
		      // compare sig-EMD & corr-EMD
		      //

		      emd_t emdc;
		      emdc.max_sift = max_sift;
		      emdc.max_imf = max_imf;
		      const int nimfc = emdc.proc( zz.col(c).data_pointer() );

		      for (int i=0; i<nimf; i++)
			{			  
			  for (int j=0; j<nimfc; j++)
			    {
			      double r = Statistics::correlation( emdc.imf[j] , emd.imf[i] );
			      if ( r <= -emd_th || r >= emd_th ) { remove.insert( i+1 ) ; break; }
			    }
			  double r = Statistics::correlation( emdc.residual , emd.imf[i] );
			  if ( r <= -emd_th || r >= emd_th ) { remove.insert( i+1 ) ; break; }
			}
		      
		      // also consider sig-EMD residual
		      
		      for (int j=0; j<nimfc; j++)
			{
			  double r = Statistics::correlation( emdc.imf[j] , emd.residual );
			  if ( r <= -emd_th || r >= emd_th ) { remove.insert( 0 ) ; break; }
			}
		      double r = Statistics::correlation( emdc.residual , emd.residual );
		      if ( r <= -emd_th || r >= emd_th ) { remove.insert( 0 ) ; break; }
		      
		    }

		  else
		    {

		      //
		      // just look at correlations between sig-EMD and corrector-RAW
		      //
		      
		      for (int i=0; i<nimf; i++)
			{
			  double r = Statistics::correlation( *zz.col(c).data_pointer() , emd.imf[i] );
			  if ( r <= -emd_th || r >= emd_th ) { remove.insert( i+1 ) ; break; } 		    
			}
		      
		      // residual
                      double r = Statistics::correlation( *zz.col(c).data_pointer(), emd.residual );
                      if ( r <= -emd_th || r >= emd_th ) { remove.insert( 0 ) ; break; }
		    }
		  
		}

	      logger << " going to remove " << remove.size() << " components\n";
	      
	      //
	      // remove y components
	      //

	      std::set<int>::const_iterator rr = remove.begin();
	      while ( rr != remove.end() )
		{
		  logger << "  rr = " << *rr << "\n";
		  for (int i=0; i<segment_points ; i++)
		    y[i] -= *rr == 0 ? emd.residual[i] : emd.imf[*rr-1][i];
		  ++rr;
		}

	      //
	      // set 'y' as final answer
	      //

	      res[ p ] = y;
	      
	    }
	  
	  
	} // next segment
	  


      //
      // Create the final signal : as we have 50% overlap, always add weights 0..1
      //

      logger << " making the final signal\n";
      
      std::vector<int> cnt( total_points , 0 );
      std::vector<double> val( total_points , 0 );

      // weight to go in 50% overlap regions  0 w 2w 3w ... 1
      const double w = 1.0 / (double)( segment_points / 2.0 - 1 );

      // TODO: apply weights
      //       take note over discontinuities / first/last regions in a region
      //   i.e. could do epoch weise, but really should apply this to all contiguous segments.
      //   i.e. should run the same SEGMENTS type of logic, or use that command
      //   for now, ignore

      std::map<int,Data::Vector<double> >::const_iterator rr = res.begin();
      while ( rr != res.end() )
	{
	  
	  const int p = rr->first;
	  int c = 0;
	  for (int i=p; i<p+segment_points; i++)
	    {
	      ++cnt[i];
	      val[i] += rr->second[ c ];
	      ++c;
	    }
	  
	  ++rr;
	}
      
      for (int i=0; i<total_points; i++)
	{
	  if ( cnt[i] == 0 ) val[i] = (*d)[i];
	  else if ( cnt[i] > 1 ) val[i] /= (double)cnt[i];
	}

      logger << " updating signal\n";
      
      //
      // Update signal in EDF
      //

      edf.update_signal_retain_range( signals(s) , &val );

    }

}

