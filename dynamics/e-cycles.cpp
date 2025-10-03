
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

#include "dynamics/e-cycles.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "timeline/hypno.h"
#include "fftw/fftwrap.h"
#include "param.h"
#include "db/db.h"
#include "helper/logger.h"
#include "stats/Eigen/Dense"
#include "stats/eigen_ops.h"
#include "dsp/tv.h"

extern writer_t writer;
extern logger_t logger;

ecycles_t::ecycles_t( edf_t & edf , param_t & param )
{
  
  // get epoch stages
  // anchor on known W/R and NREM
  // allow arbitrary definition
  // apply weighted SVD
  // smooth outputs
  // find peaks

  //
  // Options
  //

  const std::string signal_label = param.requires( "sig" );
  
  const double maxf = param.has( "max" ) ? param.requires_dbl( "max" ) : 25 ;
  const double minf = param.has( "min" ) ? param.requires_dbl( "min" ) : 0.5 ; 

  // number of components to retain (and test)
  const int nc = param.has("nc" ) ? param.requires_int( "nc" ) : 10 ;

  // window size (in epochs)
  const int nw = param.has("w" ) ? param.requires_int( "w" ) : 19 ; 

  // winsorize original signals?
  const double winsor = param.has( "winsor" ) ? param.requires_dbl( "winsor" ) : 0 ; 
      
  // include/exclude epochs (based on stage annots)
  std::set<std::string> excludes, includes;
  if ( param.has( "exclude" ) ) excludes = param.strset( "exclude" );
  if ( param.has( "include" ) ) includes = param.strset( "include" );
  
  // by default, remove leading/trailing wake
  const bool remove_prepost_wake = param.has( "flanking-wake" ) ? param.yesno( "flanking-wake" ) : true;

  // include cycles (requries that HYPNO has been run)
  bool add_classical_cycles = param.has( "cycles" ) ? param.yesno( "cycles" ) : true;

  // use TV denoiser instead of median filter
  const bool denoiser = param.has( "lambda" );
  const double lambda = denoiser ? param.requires_dbl( "lambda" ) : 0 ; 
  
  //
  // Signals
  //

  const bool no_annotations = true;

  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );

  const int ns = signals.size();

  if ( ns == 0 ) return;


  //
  // Staging 
  //

  const int ne = edf.timeline.first_epoch();
  
  edf.annotations->make_sleep_stage( edf.timeline );

  bool has_staging = edf.timeline.hypnogram.construct( &edf.timeline , param , false ); // F - not verbose
  
  if ( has_staging && edf.timeline.hypnogram.empty() )
    has_staging = false;


  // sanity check
  if ( has_staging && edf.timeline.hypnogram.stages.size() != ne )
    Helper::halt( "internal error edf.timeline.hypnogram.stages.size() != ne" );

  std::vector<sleep_stage_t> stages = edf.timeline.hypnogram.stages;
  std::vector<bool> masked_epoch( ne , false );
  std::vector<sleep_stage_t> stages1;
  std::vector<int> cycles1;
  std::vector<int> epochs1;
  
  // tidy up staging so no surprises
  // only  LIGHTS_ON NREM1 NREM2 NREM3 REM WAKE UNKNOWN 

  for (int ss=0; ss < ne ; ss++ )
    {
      if ( stages[ss] == NREM4 ) stages[ss] = NREM3;
      
      if ( ! ( stages[ss] == NREM1
	       || stages[ss] == NREM2
	       || stages[ss] == NREM3
	       || stages[ss] == REM
	       || stages[ss] == WAKE
	       || stages[ss] == LIGHTS_ON ) )
	stages[ss] = UNKNOWN;

      // these will be spliced out
      masked_epoch[ss] = stages[ss] == LIGHTS_ON || stages[ss] == UNKNOWN ; 
      
    }


  //
  // NREM cycle information
  //

  std::vector<int> cycles;

  if ( add_classical_cycles )
    {

      // data must have already been epoched (by HYPNO, so that
      // we expected NREM cycle epoch-annotations)
      if ( ! edf.timeline.epoched() )
	Helper::halt( "data not epoched: run HYPNO first or set cycles=F" );      
      
      // do we have any cycles? 
      const bool has_cycles = edf.timeline.epoch_annotation( "_NREMC_1" ) 
	|| edf.timeline.epoch_annotation( "_NREMC_2" )
	|| edf.timeline.epoch_annotation( "_NREMC_3" )
	|| edf.timeline.epoch_annotation( "_NREMC_4" )
	|| edf.timeline.epoch_annotation( "_NREMC_5" )
	|| edf.timeline.epoch_annotation( "_NREMC_6" )
	|| edf.timeline.epoch_annotation( "_NREMC_7" )
	|| edf.timeline.epoch_annotation( "_NREMC_8" );
      
      if ( ! has_cycles )
	{
	  add_classical_cycles = false; 
	  logger << "  no valid classical NREM cycles found... skipping\n";
	}
      else
	{
	  cycles.resize( ne , 0 );
	  
	  int maxc = 0;
	  
	  for (int e=0;e<ne; e++)
	    {
	      
	      int c = 0;
	      // nb. uses legacy epoch-annotation encoding
	      // take up to 8 cycles
	      // nb. epoch_annot() takes 0..ne current epoch encoding (will map as needed)
	      if      ( edf.timeline.epoch_annotation( "_NREMC_1" , e ) ) c = 1;
	      else if ( edf.timeline.epoch_annotation( "_NREMC_2" , e ) ) c = 2;
	      else if ( edf.timeline.epoch_annotation( "_NREMC_3" , e ) ) c = 3;
	      else if ( edf.timeline.epoch_annotation( "_NREMC_4" , e ) ) c = 4;
	      else if ( edf.timeline.epoch_annotation( "_NREMC_5" , e ) ) c = 5;
	      else if ( edf.timeline.epoch_annotation( "_NREMC_6" , e ) ) c = 6;
	      else if ( edf.timeline.epoch_annotation( "_NREMC_7" , e ) ) c = 7;
	      else if ( edf.timeline.epoch_annotation( "_NREMC_8" , e ) ) c = 8;
	      
	      if ( c > maxc ) maxc = c;
	      
	      cycles[e] = c;
	      // next epoch
	    }
	  
	  logger << "  added flags for " << maxc << " classical NREM cycles\n";
	}
    }

  
  //
  // remove pre-/post-sleep wake
  //
  
   if ( remove_prepost_wake )
    {
      for (int ss=0; ss < ne ; ss++ )
	{
	  if ( stages[ss] == NREM1 || stages[ss] == NREM2 || stages[ss] == NREM3 || stages[ss] == REM )
	    break;
	  masked_epoch[ss] = true;
	}

      for (int ss= ne-1; ss >= 0 ; ss-- )
	{
          if ( stages[ss] == NREM1 || stages[ss] == NREM2 || stages[ss] == NREM3 || stages[ss] == REM )
            break;
          masked_epoch[ss] = true;
        }

    }

  
  //
  // Other masking?
  //

  if ( includes.size() )
    {
      std::set<sleep_stage_t> inc;
      std::set<std::string>::const_iterator ss = includes.begin();
      while ( ss != includes.end() )
	{
	  if      ( *ss == "N1" ) inc.insert( NREM1 );
	  else if ( *ss == "N2" ) inc.insert( NREM2 );
	  else if ( *ss == "N3" ) inc.insert( NREM3 );
	  else if ( *ss == "R" ) inc.insert( REM );
	  else if ( *ss == "W" ) inc.insert( WAKE ); 
	  ++ss;
	}

      for (int ss=0; ss < ne ; ss++ )
	if ( inc.find( stages[ss] ) == inc.end() ) 
	  masked_epoch[ss] = true;
    }

    if ( excludes.size() )
    {
      std::set<sleep_stage_t> exc;
      std::set<std::string>::const_iterator ss = excludes.begin();
      while ( ss != excludes.end() )
	{
	  if      ( *ss == "N1" ) exc.insert( NREM1 );
	  else if ( *ss == "N2" ) exc.insert( NREM2 );
	  else if ( *ss == "N3" ) exc.insert( NREM3 );
	  else if ( *ss == "R" ) exc.insert( REM );
	  else if ( *ss == "W" ) exc.insert( WAKE ); 
	  ++ss;
	}

      for (int ss=0; ss < ne ; ss++ )
	if ( exc.find( stages[ss] ) != exc.end() ) 
	  masked_epoch[ss] = true;
    }
    
  //
  // SR check
  //

  const int Fs = edf.header.sampling_freq( signals(0) );
  logger << "  sample rate = " << Fs << "\n";
  for (int s=1; s<ns; s++)
    if ( edf.header.sampling_freq( signals(s) ) != Fs )
      Helper::halt( "all sample rates must be similar" );
  
  //
  // Epoch-wise PSD
  //

  // figure out number of features per channel on the first epoch
  
  const int nf = ns * ( maxf - minf + 0.25 ) / 0.25 ;

  Eigen::MatrixXd X1 = Eigen::MatrixXd::Zero( ne , nf );

  logger << "  deriving " << nf << " spectral features across " << ns << " channels\n";
  
  //
  // Get size of spectra: fixed param here for now
  //
  
  const double overlap_sec = 2;
  const double segment_sec  = 4;  
  const int segment_points = segment_sec * Fs;
  const int noverlap_points  = overlap_sec * Fs;

  window_function_t window_function = WINDOW_TUKEY50;
  if      ( param.has( "no-window" ) ) window_function = WINDOW_NONE;
  else if ( param.has( "hann" ) ) window_function = WINDOW_HANN;
  else if ( param.has( "hamming" ) ) window_function = WINDOW_HAMMING;
  else if ( param.has( "tukey50" ) ) window_function = WINDOW_TUKEY50;

  const bool mean_centre_epoch = param.has( "center" ) || param.has( "centre" )
    || param.has( "mean-center" ) || param.has( "mean-centre" );
  
  // detrend signal first?
  const bool remove_linear_trend = param.has( "detrend" );

  if ( mean_centre_epoch && remove_linear_trend )
    Helper::halt( "cannot specify both mean-center and detrend" );

  
  const bool use_seg_median = true;
  const bool calc_seg_sd = false;
  const bool average_adj = false;
  const bool use_nextpow2 = false;

  logger << "  iterating over " << ne << " epochs\n";
  
  //
  // Populate spectral matrix, iterating over epochs
  //
      
  int row = 0;
  int ecnt = 0;
  
  edf.timeline.first_epoch();
  
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      
      
      if ( epoch == -1 ) break;

      // skip this epoch? (and do not advance row)
      if ( masked_epoch[ ecnt ] )
	{
	  ++ecnt;
	  continue; 
	}
      else
	{
	  stages1.push_back( stages[ ecnt ] );
	  if ( add_classical_cycles ) 
	    cycles1.push_back( cycles[ ecnt ] );
	  epochs1.push_back( edf.timeline.display_epoch( epoch ) ); // 1-based epoch count		    
	  ++ecnt;
	}
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      // Need to check segment length?
      
      if ( edf.timeline.generic_epochs() )
	{
	  // here, all epochs need to have the same segment length, so will skip here
	  // if the epoch is too short
	  if ( edf.timeline.epoch_length() < segment_sec )
	    Helper::halt( "cannot have epoch length shorter than segment size" );
	}	  


      //
      // Iterate over signals
      //

      int col = 0 ;

      for (int s=0; s<ns; s++)
	{
      
	  // Get data
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();
	  
	  //
	  // mean centre epoch?
	  //
	  
	  if ( mean_centre_epoch ) 
	    MiscMath::centre( d );
	  else if ( remove_linear_trend )
	    MiscMath::detrend( d );
	  
	   //
	   // pwelch() to obtain full PSD
	   //
	  
	   const int total_points = d->size();
	   
	   // implied number of segments
	   int noverlap_segments = floor( ( total_points - noverlap_points) 
					  / (double)( segment_points - noverlap_points ) );

	   PWELCH pwelch( *d , 
			  Fs , 
			  segment_sec , 
			  noverlap_segments , 
			  window_function , 
			  use_seg_median,
			  calc_seg_sd,
			  average_adj ,
			  use_nextpow2 );

	   //
	   // populate X1
	   //
	   
	   const double bin_fac = 1.0;
	   bin_t bin( minf , maxf , bin_fac );
	   bin.bin( pwelch.freq , pwelch.psd );
	   
	   for ( int i = 0 ; i < bin.bfa.size() ; i++ )
	     X1(row,col++) = 10*log10( bin.bspec[i] ) ;
	   	   
	   //
	   // next signal
	   //
	}

      ++row;

      //
      // next epoch
      //
    }


  //
  // if skipped
  //
  
  if ( row < ne )
    {
      logger << "  resizing X to " << row << " included rows from " << ne << " epochs\n";
      if ( row < nw ) Helper::halt( "too few valid epochs remaining" );	
      X1.conservativeResize( row , Eigen::NoChange );
    }



  //
  // Smoothing and norming
  //
  
  eigen_ops::robust_scale( X1 , true , true , -1 , true );
  //eigen_ops::scale( X1 , true , true );

  if ( denoiser ) // denoiser, assuming SD ~ 1 
    for (int v=0; v<nf; v++)
      dsptools::TV1D_denoise( X1.col(v) , lambda );
  else // median filter    
    for (int v=0; v<nf; v++)
      X1.col(v) = eigen_ops::median_filter( X1.col(v) , nw );

  //  std::cout << "X1 orig\n" << X1.col(0) << "\n";
  
  //
  // Weighting
  //
  
  Eigen::VectorXd W = Eigen::VectorXd::Ones( nf );

  for (int v=0; v<nf; v++)
    {
      //std::cout << " f = " << v << "\n";
      
      const Eigen::VectorXd & V = X1.col(v);
      //      std::cout << "V size = " << V.size() << " " << stages1.size() << " " << row << "\n";
      double t_nrem = 0, t_rem_wake = 0 , t_rem = 0 , t_wake = 0;
      int n_nrem = 0 , n_rem_wake = 0 , n_rem = 0 , n_wake = 0;
      
      for (int i=0; i<row; i++)
	{
	  if ( stages1[i] == NREM1 || stages1[i] == NREM2 || stages1[i] == NREM3 )
	    {
	      t_nrem += V[i];
	      n_nrem++;
	    }
	  else if ( stages1[i] == REM ) 
	    {
	      t_rem_wake += V[i];
	      n_rem_wake++;
	      t_rem += V[i];
	      n_rem++;
	    }
	  else if ( stages1[i] == WAKE )
	    {
	      t_rem_wake += V[i];
	      n_rem_wake++;
	      t_wake += V[i];
	      n_wake++;
	    }
	}
      
      
      // means
      if ( n_nrem ) t_nrem /= (double)n_nrem;
      if ( n_rem_wake ) t_rem_wake /= (double)n_rem_wake;
      if ( n_rem ) t_rem /= (double)n_rem;
      if ( n_wake ) t_wake /= (double)n_wake;
	    
      // by default compare NREM vs (REM+WAKE)      
      // make a weight?

      //      if ( n_nrem >= 10 && n_rem_wake >= 10 )
      //	W[v] = fabs( t_nrem - t_rem_wake );

      if ( n_nrem >= 10 && n_rem_wake >= 10 )
	W[v] = fabs( t_rem - t_wake );

            
    }

  // now W[] is populated (or all ones)
  // scale
  
  double Wmax = W.maxCoeff();
  W /= Wmax;
  
  for (int v=0; v<nf; v++)
    X1.col(v) *= W[v];
  
  
  //
  // SVD
  //

  logger << "  about to perform SVD...\n";

  Eigen::BDCSVD<Eigen::MatrixXd> svd( X1 , Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXd U = svd.matrixU();

  // components are already sorted in decreasing order, so so just take the first 'nc' components
  
  //  U[ ni x ni ]  -->   U[ ni x nc ]
  U.conservativeResize( Eigen::NoChange , nc );

  //
  // Smooth & norm U  
  //

  eigen_ops::robust_scale( U , true , true , winsor , true );
  
  if ( denoiser ) // denoiser, assuming SD ~ 1 
    for (int v=0; v<nc; v++)
      dsptools::TV1D_denoise( U.col(v) , lambda );
  else
    for (int v=0; v<nc; v++)
      U.col(v) = eigen_ops::median_filter( U.col(v) , nw );
     

  //
  // Orient all PCs to be R +ve, NR -ve
  //

  for (int v=0; v<nc; v++)
    {
      double t_nrem = 0, t_rem_wake = 0 ;
      int n_nrem = 0 , n_rem_wake = 0 ;
      
      for (int i=0; i<row; i++)
        {
          if ( stages1[i] == NREM1 || stages1[i] == NREM2 || stages1[i] == NREM3 )
            {
              t_nrem += U(i,v);
              n_nrem++;
            }
          else if ( stages1[i] == REM || stages1[i] == WAKE )
            {
              t_rem_wake += U(i,v);
              n_rem_wake++;
            }
    	}

      if ( n_nrem && n_rem_wake )
	{
	  t_nrem /= (double) n_nrem;
	  t_rem_wake /= (double) n_rem_wake;
	  // flip?
	  if ( t_nrem > t_rem_wake ) 
	    U.col(v) *= -1.0;
	}
      
    }

      
      

  
  //
  // output
  //
  
  for (int e=0; e<row; e++)
    {
      writer.epoch( epochs1[e] );
      
      for (int c=0; c<nc; c++)
	{
	  writer.level( c+1 , globals::comp_strat );	      
	  writer.value( "U" , U(e,c) );
	}
      writer.unlevel( globals::comp_strat );
      
      if ( add_classical_cycles ) 
	writer.value( "CYC" , cycles1[e] );
      
      writer.value( "SS" , globals::stage( stages1[e] ) );
      
    }
  writer.unepoch();

 
}
  

