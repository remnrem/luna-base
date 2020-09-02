
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

#include "dirent.h"

#include "stats/matrix.h"
#include "stats/statistics.h"
#include "stats/lda.h"
#include "miscmath/miscmath.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "dsp/resample.h"
#include "fftw/fftwrap.h"
#include "dsp/tv.h"

extern logger_t logger;

extern writer_t writer;


bool suds_t::verbose = false;
std::set<suds_indiv_t> suds_t::bank;
std::set<suds_indiv_t> suds_t::wbank;
int suds_t::nc;
int suds_t::ns;
std::vector<std::string> suds_t::siglab;
std::vector<double> suds_t::lwr;
std::vector<double> suds_t::upr;
std::vector<int> suds_t::fac;
std::vector<int> suds_t::sr;
double suds_t::wgt_percentile;
double suds_t::denoise_fac;
bool suds_t::standardize_u = true;
bool suds_t::use_best_guess = true;
bool suds_t::ignore_target_priors = false;
std::vector<double> suds_t::outlier_ths;
int suds_t::required_epoch_n = 5;


// 5 -> N1, N2, N3, R, W
// 3 -> NR, R, W (3-stage)
int suds_t::n_stages = 5; 
                            
bool suds_t::use_fixed_trainer_req;
std::vector<int> suds_t::fixed_trainer_req;
int suds_t::fake_ids;
std::string suds_t::fake_id_root;

std::string suds_t::eannot_file = "";
bool suds_t::eannot_ints = false;
std::string suds_t::eannot_prepend = "";
std::string suds_t::mat_dump_file = "";

std::vector<double> suds_t::lwr_h2, suds_t::upr_h2;
std::vector<double> suds_t::lwr_h3, suds_t::upr_h3;

// 1) Epoch-level PSD for 1+ channels (cs_EEG, etc) 
//    These can be given arbitrary labels for targets, but then
//    targets must match these

// 2) PSC to extract N components from the epoch-level

// 3) Outlier removal

// 4) If manual staging data present, fit LDA
//    -->  save in library

// 5) If manual staging data not present, predict using each of N trainers
//    (after projecting epoch-level PSD into the trainer space) 
//   --> evaluate how well this trainer worked by fitting LDA to predicted stages, as seeing how well that
//       manages to impute other trainers (with known staging)


void suds_indiv_t::add_trainer( edf_t & edf , param_t & param )
{

  // build a trainer
  int n_unique_stages = proc( edf , param , true );
  
  // only include recordings that have all five/three stages included, for now
  if ( n_unique_stages != suds_t::n_stages ) return;
  
  // save to disk
  write( edf , param ); 

}




int suds_indiv_t::proc( edf_t & edf , param_t & param , bool is_trainer )
{

  //
  // Is this individual a trainer (with known stages) or no?
  //

  trainer = is_trainer;
  
  const int nc = suds_t::nc;
  const int ns = suds_t::ns;

  
  //
  // Signals from this EDF
  //

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );

  if ( signals.size() != ns ) Helper::halt( "could not find specified signals" );
  
  //
  // Resample as needed
  //
  
  for (int s=0;s<ns;s++)
    {
      
      if ( edf.header.is_annotation_channel( signals(s) ) )
	Helper::halt( "cannot specificy annotation channel: " + signals.label(s) );
      
      if ( edf.header.sampling_freq( signals(s) ) != suds_t::sr[s] ) 
	dsptools::resample_channel( edf, signals(s) , suds_t::sr[s] );

    }


  //
  // Epoch 
  //

  const int ne = edf.timeline.first_epoch();

  // nb.  below:
  //   ne   total number of epochs
  //   nge  num of epochs with 'valid' staging (i.e. no UNKNOWN, etc)
  //   nve  of nge, number that are not statistical outliers for 1+ PSC [ stored in output ]

								      
  //
  // PSD
  //

  //  misc param
  double fft_segment_size = param.has( "segment-sec" ) 
    ? param.requires_dbl( "segment-sec" ) : 4 ;
  
  double fft_segment_overlap = param.has( "segment-overlap" ) 
    ? param.requires_dbl( "segment-overlap" ) : 2 ;
  
  if ( edf.timeline.epoch_length() <= ( fft_segment_size + fft_segment_overlap ) )
    {
      fft_segment_overlap = 0;
      fft_segment_size = edf.timeline.epoch_length();
    }

  window_function_t window_function = WINDOW_TUKEY50;	   
  if      ( param.has( "no-window" ) ) window_function = WINDOW_NONE;
  else if ( param.has( "hann" ) ) window_function = WINDOW_HANN;
  else if ( param.has( "hamming" ) ) window_function = WINDOW_HAMMING;
  else if ( param.has( "tukey50" ) ) window_function = WINDOW_TUKEY50;
  
  
  //
  // Get stage information (for trainers only)
  //

  std::vector<bool> retained( ne , true );

  bool has_prior_staging = false;
  
  if ( trainer )
    {
      edf.timeline.annotations.make_sleep_stage();
      
      if ( ! edf.timeline.hypnogram.construct( &edf.timeline , false ) )
	Helper::halt( "problem extracting stage information for trainer" );
      
      // total number of epochs does not match?
      if ( ne != edf.timeline.hypnogram.stages.size() )
	Helper::halt( "problem extracting stage information for trainer" );

      has_prior_staging = true;
      
    }
  else if ( ! suds_t::ignore_target_priors )
    {

      // for targets, manual/prior staging may exist, in which case we'll want to track it for comparison
      // unless we've been explicitly told to ignore it (ignore-prior --> suds_t::ignore_target_priors )

      edf.timeline.annotations.make_sleep_stage();

      has_prior_staging = edf.timeline.hypnogram.construct( &edf.timeline , false ) ;
      
      if ( has_prior_staging )
	{
	  // total number of epochs does not match?
	  if ( ne != edf.timeline.hypnogram.stages.size() )
	    Helper::halt( "problem extracting stage information for trainer" );
	}
      
    }

  
  // nb. overkill to have two staging classifications, but keep for now,
  // i.e. we may want to extend one type only

  // number of good (retained) epochs

  int nge = 0 ;
  
  if ( has_prior_staging )
    {
      
      obs_stage.resize( ne , SUDS_UNKNOWN );
      
      for (int ss=0; ss < ne ; ss++ )
	{
	  if ( edf.timeline.hypnogram.stages[ ss ] == UNSCORED
	       || edf.timeline.hypnogram.stages[ ss ] == LIGHTS_ON
	       || edf.timeline.hypnogram.stages[ ss ] == MOVEMENT
	       || edf.timeline.hypnogram.stages[ ss ] == UNKNOWN ) obs_stage[ss] = SUDS_UNKNOWN;
	  
	  else if ( edf.timeline.hypnogram.stages[ ss ] == WAKE ) obs_stage[ss] = SUDS_WAKE;
	  else if ( edf.timeline.hypnogram.stages[ ss ] == NREM1 ) obs_stage[ss] = SUDS_N1;
	  else if ( edf.timeline.hypnogram.stages[ ss ] == NREM2 ) obs_stage[ss] = SUDS_N2;
	  else if ( edf.timeline.hypnogram.stages[ ss ] == NREM3
		    || edf.timeline.hypnogram.stages[ ss ] == NREM4 )  obs_stage[ss] = SUDS_N3;
	  
	  else if ( edf.timeline.hypnogram.stages[ ss ] == REM ) obs_stage[ss] = SUDS_REM;

	  // expand retained class to exclude unknown staging info
	  // note: even for targets, this means we will not try to guess epochs
	  // that have an UNKNOWN value for the 
	  if ( obs_stage[ss] == SUDS_UNKNOWN ) { retained[ss] = false; } 
	  else ++nge;
	}
    }
  else 
    {

      //
      // for target individuals without staging, include all epochs
      //
      
      nge = ne;
      
    }


  //
  // for QC, estimate Hjorth parameters (only 2nd and 3rd used) over
  // epochs (for each signal) 
  //
  
  h2.resize( nge , ns );
  h3.resize( nge , ns );
    
  //
  // iterate over (retained) epochs
  //
    
  int en = 0 , en_good = 0;

  edf.timeline.first_epoch();
  
  epochs.clear();
  
  while ( 1 ) 
    {
      
      // select epoch
      int epoch = edf.timeline.next_epoch();      	  
      
      if ( epoch == -1 ) break;
      
      if ( en == ne ) Helper::halt( "internal error: over-counted epochs" );

      // retained? if not, skip
      if ( ! retained[ en ] ) 
	{
	  ++en;
	  continue;
	}

      // col counter for PSD aggregation matrix
      int col = 0;
      std::vector<double> firstrow;
      
      // iterate over signals      
      for (int s = 0 ; s < ns; s++ )
	{
	  // get data
	  interval_t interval = edf.timeline.epoch( epoch );	  
	  slice_t slice( edf , signals(s) , interval );
	  std::vector<double> * d = slice.nonconst_pdata();

	  // mean centre epoch
	  MiscMath::centre( d );
	  
	  // pwelch() to obtain full PSD
	  const double overlap_sec = fft_segment_overlap;
	  const double segment_sec  = fft_segment_size;
	  const int total_points = d->size();
	  const int segment_points = segment_sec * suds_t::sr[s];
	  const int noverlap_points  = overlap_sec * suds_t::sr[s];
	  
	  // implied number of segments
	  int noverlap_segments = floor( ( total_points - noverlap_points) 
					 / (double)( segment_points - noverlap_points ) );
	  
	  PWELCH pwelch( *d , 
			 suds_t::sr[s] , 
			 segment_sec , 
			 noverlap_segments , 
			 window_function );
	  	  	  
	  // using bin_t 	      
	  bin_t bin( suds_t::lwr[s] , suds_t::upr[s] , suds_t::fac[s] );	  
	  
	  bin.bin( pwelch.freq , pwelch.psd );

	  for ( int i = 0 ; i < bin.bfa.size() ; i++ )
	    {
	      
	      if ( bin.bfa[i] >= suds_t::lwr[s] && bin.bfb[i] <= suds_t::upr[s] )
		{
		  if ( en_good == 0 ) firstrow.push_back(  10*log10( bin.bspec[i] ) );
		  else PSD( en_good , col ) = 10*log10( bin.bspec[i] ) ; 		  
		  ++col;		  
		}	      
	    }
	  

	  // Hjorth parameters

	  double activity = 0 , mobility = 0 , complexity = 0;
	  MiscMath::hjorth( d , &activity , &mobility , &complexity );
	  h2[en_good][s] = mobility ;
	  h3[en_good][s] = complexity ;
	  
	  
	} // next signal

      
      // store/shape output if first go around
      nbins = col;
      if ( en_good == 0 )
	{
	  PSD.resize( nge , col );
	  for (int i=0;i<col;i++) PSD(0,i) = firstrow[i];
	}

      
      // increase epoch-number
      ++en;
      ++en_good;
      epochs.push_back( epoch );
    
    } // next epoch
  
  // all done: check

  if ( en_good != nge ) Helper::halt( "internal error: under-counted epochs" );


  //
  // Rescale PSD?
  //

  if ( suds_t::standardize_u )
    {
      logger << "  standardizing PSD\n";
      Statistics::standardize( PSD );
    }
  

  //
  // Get PSC
  //

  // mean-center columns

  U = PSD;
  
  Statistics::mean_center_cols( U );

  // SVD
  
  W.clear(); V.clear();
  W.resize( nbins ); 
  V.resize( nbins , nbins );

  bool okay = Statistics::svdcmp( U , W , V );
  if ( ! okay ) Helper::halt( "problem with SVD" );

  int rank = Statistics::orderSVD( U , W , V );
  if ( rank == 0 ) Helper::halt( "problem with input data, rank 0" );



  //
  // Outliers/smoothing
  //

  std::vector<bool> valid( nge , true );

  // track reasons for exclusion
  int nout_flat = 0;
  int nout_hjorth = 0;
  int nout_stat = 0;
  
  //
  // Exclusions based on H==0 parameters
  //
  
  for ( int s=0;s<ns;s++)
    {
      for (int i=0;i<nge;i++) 
	{
	  if ( h2( i, s ) < 1e-8 ) { valid[i] = false; nout_flat++; }
	  else if ( h3( i , s ) < 1e-8 ) { valid[i] = false; nout_flat++; } 
	}
    }

  //
  // For targets only, threshold epochs based on per-signal Hjorth
  // from trainers
  // 

  if ( ! trainer )
    {
      for ( int s=0;s<ns;s++)
	{
	  for (int i=0;i<nge;i++)
	    {
	      if ( h2( i, s ) <= suds_t::lwr_h2[s] || h2(i,s) >= suds_t::upr_h2[s] ||
		   h3( i, s ) <= suds_t::lwr_h3[s] || h3(i,s) >= suds_t::upr_h3[s] ) { valid[i] = false; nout_hjorth++; } 
	    }
	}
    }


  //
  // Spectra (after removing flat lines)
  //

  for (int o=0;o< suds_t::outlier_ths.size();o++)
    {
      for ( int j=0;j<nc;j++)
	{
	  std::vector<double> x;
	  for (int i=0;i<nge;i++) if ( valid[i] ) x.push_back( U(i,j) );
	  if ( x.size() < 2 ) Helper::halt( "no epochs left" );
	  double mean = MiscMath::mean( x );
	  double sd = MiscMath::sdev( x , mean );
	  double lwr = mean - suds_t::outlier_ths[o] * sd;
	  double upr = mean + suds_t::outlier_ths[o] * sd;
	  int c = 0;
	  for (int i=0;i<nge;i++)
	    {
	      if ( valid[i] )
		{
		  if ( x[c] < lwr || x[c] > upr ) { valid[i] = false; nout_stat++; } 
		  ++c;
		}
	    }
	}
    }





  int included = 0;
  for (int i=0;i<nge;i++)
    if ( valid[i] ) ++included;

  logger << "  of " << ne << " total epochs, valid staging for " << nge
         << ", and of those " << included << " passed outlier removal\n";

  logger << " (detailed outlier counts: " << nout_flat << ", " << nout_hjorth << ", " << nout_stat << ")\n";

  //
  // Remove bad epochs and repeat (SVD and smoothing)
  //

  // nve = number of valid epochs ( ne > nge > nve ) 

  
  nve = included;

  Data::Matrix<double> PSD2 = PSD;
  PSD.clear();
  PSD.resize( nve , nbins );
  std::vector<int> epochs2 = epochs;
  epochs.clear();

  if ( has_prior_staging )
    obs_stage_valid.clear();
  
  int r = 0;
  for (int i=0;i<PSD2.dim1() ; i++)
    {      
      if ( valid[i] )
	{
	  for (int j=0;j<nbins;j++)
	    PSD(r,j) = PSD2(i,j);
	  
	  epochs.push_back( epochs2[i] );
	  
	  if ( has_prior_staging )
	    obs_stage_valid.push_back( obs_stage[i] );
	  
	  ++r;
	}
    }

  
  // remove bad epochs for Hjorth parameters

  Data::Matrix<double> hh2 = h2;
  Data::Matrix<double> hh3 = h3;
  h2.clear(); h3.clear();
  h2.resize( nve , ns ); h3.resize( nve , ns );

  for (int s=0;s<ns;s++)
    {
      int r = 0;      
      for (int i=0; i < valid.size(); i++ )
	if ( valid[i] )
	  {
	    h2(r,s) = hh2(i,s);
	    h3(r,s) = hh3(i,s);
	    ++r;
	  }
    }


  //
  // Rescale PSD?
  //

  if ( suds_t::standardize_u )
    {
      //logger << "  standardizing PSD, round 2\n";
      Statistics::standardize( PSD );
    }
  else // just ensure we mean-center in any case
    {
      // mean-center columns (PSD)  
      Statistics::mean_center_cols( PSD );

    }


  //
  // Get PSC (post outlier removal)
  //

  
  // copy to U
  
  U  = PSD;
  
  // get PSCs

  W.clear(); V.clear();
  W.resize( nbins ); 
  V.resize( nbins , nbins );
  
  okay = Statistics::svdcmp( U , W , V );
  if ( ! okay ) Helper::halt( "problem with SVD" );

  rank = Statistics::orderSVD( U , W , V );
  if ( rank == 0 ) Helper::halt( "problem with input data, rank 0" );
  

  //
  // Smooth PSCs?
  //

  if ( suds_t::denoise_fac > 0 ) 
    for (int j=0;j<nc;j++)
      {
	std::vector<double> * col = U.col_nonconst_pointer(j)->data_nonconst_pointer();
	double sd = MiscMath::sdev( *col );
	double lambda = suds_t::denoise_fac * sd;
	dsptools::TV1D_denoise( *col , lambda );
      }
    

      
  //
  // Make variables for LDA: shrink down to 'nc'
  //

  Data::Matrix<double> U2 = U;
  U.clear();
  U.resize( nve , nc );
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      U(i,j) = U2(i,j);
  
  W.resize( nc );
  Data::Matrix<double> VV = V;
  V.clear();
  V.resize( VV.dim1() , nc );
  for (int i=0;i<VV.dim1();i++)
    for (int j=0;j<nc;j++)
      V(i,j) = VV(i,j);


  //
  // Summarize mean/SD for per-signal Hjorth parameters
  //

  mean_h2 = Statistics::mean( h2 );
  mean_h3 = Statistics::mean( h3 );

  sd_h2 = Statistics::sdev( h2 , mean_h2 ) ;
  sd_h3 = Statistics::sdev( h3 , mean_h3 ) ;


  // std::cout << "h2\n" << h2.print() << "\n";
  // std::cout << "h3\n" << h3.print() << "\n";

  //
  // make class labels ( trainer only )
  //
  
  if ( trainer )
    {
      
      y.clear();
      
      int c = 0;
      for ( int i = 0 ; i < ne ; i++ )
	{
	  if ( retained[i] )
	    {
	      if ( valid[c] )
		y.push_back( suds_t::str( obs_stage[i] ) );
	      ++c;
	    }
	}

      counts.clear();
      for (int i=0;i<y.size();i++) counts[y[i]]++;
      std::map<std::string,int>::const_iterator cc = counts.begin();
      while ( cc != counts.end() )
	{
	  logger << "   " << cc->second << " " << cc->first << " epochs\n";
	  ++cc;
	}
      
    }

  // for trainers, returns number of observed stages w/ at least suds_t::required_epoch_n 
  // -- i.e. should be suds_t::n_stages

  int nr = 0;
  std::map<std::string,int>::const_iterator cc = counts.begin();
  while ( cc != counts.end() )
    {
      if ( cc->second >= suds_t::required_epoch_n ) ++nr;      
      ++cc;
    }
  
  return nr;
  
}



void suds_indiv_t::write( edf_t & edf , param_t & param ) const
{

  const std::string folder = param.requires( "db" );

  const int nc = suds_t::nc;
  const int ns = suds_t::ns;
    
  //
  // Store as an epoch-level EDF
  //

  // Save: 
  //  ID
  //  U [ nve , nc ]  D [ nc ]  V [ nc , nc ] 
  //  stages [ nve ] 

  // create output folder if it does not exist
  std::string syscmd = "mkdir -p " + folder ;
  int retval = system( syscmd.c_str() );


  // for saving trainers: use EDF ID, or a fake ID?  (e.g. 'ids=suds')
  std::string suds_id = suds_t::fake_ids ? suds_t::fake_id_root + "_" + Helper::int2str( suds_t::fake_ids++ ) : edf.id;
  
  std::string filename = folder + globals::folder_delimiter + suds_id;

  logger << "  writing trainer data to " << filename << "\n";
  
  std::ofstream OUT1( filename.c_str() , std::ios::out );

  // file version code
  OUT1 << "SUDS\t1\n";
  
  OUT1 << "ID\t" << suds_id << "\n"
       << "N_VALID_EPOCHS\t" << nve << "\n"
       << "N_X\t" << nbins << "\n"
       << "N_SIGS\t" << ns << "\n"
       << "N_COMP\t" << nc << "\n";

  // channels , SR [ for comparability w/ other data ] 

  for (int s=0;s<ns;s++)
    {
      OUT1 << "\nCH\t" << suds_t::siglab[s] << "\n"
	   << "SR\t" << suds_t::sr[s] << "\n"
	   << "LWR\t" << suds_t::lwr[s] << "\n"
	   << "UPR\t" << suds_t::upr[s] << "\n"
	   << "FAC\t" << suds_t::fac[s] << "\n"
	   << "H2_MN\t" << mean_h2[s] << "\n"
	   << "H2_SD\t" << sd_h2[s] << "\n"
	   << "H3_MN\t" << mean_h3[s] << "\n"
	   << "H3_SD\t" << sd_h3[s] << "\n";
    }


  // stages
  OUT1 << "\nN_STAGES\t" << counts.size() << "\n";
  
  std::map<std::string,int>::const_iterator ss = counts.begin();
  while ( ss != counts.end() )
    {
      OUT1 << ss->first << "\t" << ss->second << "\n";
      ++ss;
    }
  
  // W
  OUT1 << "\nW[" << nc << "]";
  for (int j=0;j<nc;j++)
    OUT1 << " " << W[j];
  OUT1 << "\n";

  // V
  OUT1 << "\nV[" << nbins << "," << nc << "]";
  for (int i=0;i<nbins;i++)
    for (int j=0;j<nc;j++)
      OUT1 << " " << V(i,j);
  OUT1 << "\n";
  
  // stages
  OUT1 << "\nEPOCH_STAGE";
  for (int i=0;i<nve;i++)
    OUT1 << " " << epochs[i] << " " << y[i] ;
  OUT1 << "\n\n";

  // U (to reestimate LDA model upon loading, i.e
  //  to use lda.predict() on target   ; only needs to be nc rather than nbins
  OUT1 << "U[" << nve << "," << nc << "]";
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      OUT1 << " " << U(i,j);
  OUT1 << "\n\n";
  
  // X, RAW DATA (e.g. mean-centered PSD, but possibly other things)
  // i.e. if this trainer is being used as a 'weight trainer',
  // i.e. will project this individuals raw data into the target space
  OUT1 << "X[" << nve << "," << nbins << "]";
  for (int i=0;i<nve;i++)
    for (int j=0;j<nbins;j++)
      OUT1 << " " << PSD(i,j);
  OUT1 << "\n\n";
  
  OUT1.close();

  //
  // All done
  //

}


void suds_indiv_t::reload( const std::string & filename , bool load_rawx )
{
  
  //  logger << "  reloading trainer data from " << filename << "\n";
  
  std::ifstream IN1( filename.c_str() , std::ios::in );

  std::string dummy;
  std::string suds;
  int version;
  IN1 >> suds >> version;
  
  if ( suds != "SUDS" )
    Helper::halt( "bad file format for " + filename );
  
  int this_ns, this_nc;

  IN1 >> dummy >> id 
      >> dummy >> nve 
      >> dummy >> nbins 
      >> dummy >> this_ns 
      >> dummy >> this_nc ;
  
  if ( this_nc != suds_t::nc || this_ns != suds_t::ns )
    Helper::halt( "different trainer nc=" + Helper::int2str( this_nc ) + " in " + filename ); 
  
  const int nc = suds_t::nc;
  const int ns = suds_t::ns;
  
  mean_h2.clear(); mean_h3.clear();
  sd_h2.clear(); sd_h3.clear();

  for (int s=0;s<ns;s++)
    {
      std::string this_siglab;
      double this_lwr, this_upr;
      int this_sr, this_fac;
      double this_h2m, this_h2sd, this_h3m, this_h3sd;
      IN1 >> dummy >> this_siglab
	  >> dummy >> this_sr 
	  >> dummy >> this_lwr
	  >> dummy >> this_upr
	  >> dummy >> this_fac 
	  >> dummy >> this_h2m >> dummy >> this_h2sd 
	  >> dummy >> this_h3m >> dummy >> this_h3sd;

      mean_h2.push_back( this_h2m );
      mean_h3.push_back( this_h3m );
      sd_h2.push_back( this_h2sd );
      sd_h3.push_back( this_h3sd );

      if ( this_siglab != suds_t::siglab[s] ) Helper::halt( "different signals: " + this_siglab
							    + ", but expecting " + suds_t::siglab[s] );
      if ( this_sr != suds_t::sr[s] ) Helper::halt( "different SR: " + this_siglab
						    + ", but expecting " + suds_t::siglab[s] );
      if ( this_lwr != suds_t::lwr[s] ) Helper::halt( "different lower-freq: " + Helper::dbl2str( this_lwr )
						      + ", but expecting " + Helper::dbl2str( suds_t::lwr[s] ) );
      if ( this_upr != suds_t::upr[s] ) Helper::halt( "different upper-freq: " + Helper::dbl2str( this_upr )
						      + ", but expecting " + Helper::dbl2str( suds_t::upr[s] )) ;
      if ( this_fac != suds_t::fac[s] ) Helper::halt( "different fac: " + Helper::int2str( this_fac )
						      + ", but expecting " + Helper::int2str( suds_t::fac[s] ) ) ;
      
    }

  // stages
  int nstages;
  IN1 >> dummy >> nstages;
  for (int i=0;i<nstages;i++)
    {
      std::string sname;
      int scnt;
      IN1 >> sname >> scnt;
      counts[ sname ] = scnt;
    }
  
  
  // check that these equal suds_t values?
  
  // W [ only nc ]
  IN1 >> dummy;
  W.resize( nc );
  for (int j=0;j<nc;j++)
    IN1 >> W[j] ;

  // V [ only nc cols ] 
  IN1 >> dummy;
  V.resize( nbins, nc );
  for (int i=0;i<nbins;i++)
    for (int j=0;j<nc;j++)
      IN1 >> V(i,j);

  // stages
  IN1 >> dummy;
  y.resize( nve );
  epochs.resize( nve );
  for (int i=0;i<nve;i++)
    IN1 >> epochs[i] >> y[i] ;
  obs_stage = suds_t::type( y );
  
  // U (to reestimate LDA model upon loading, i.e
  //  to use lda.predict() on target 
  IN1 >> dummy;
  U.resize( nve , nc );
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      IN1 >> U(i,j) ;	

  // X values (e.g. mean-centered PSD)
  // i.e. if this trainer is being used as a 'weight trainer',
  // i.e. will project this individuals raw data into the target space
  
  if ( load_rawx )
    {      
      // PSD
      IN1 >> dummy;
      PSD.resize( nve , nbins );
      for (int i=0;i<nve;i++)
	for (int j=0;j<nbins;j++)
	  IN1 >> PSD(i,j) ;
    }

  IN1.close();
  
  //
  // All done
  //

}


void suds_t::attach_db( const std::string & folder , bool read_psd )
{

  std::set<suds_indiv_t> * b = read_psd ? &wbank : &bank ;
    
  // already done?
  if ( b->size() > 0 ) return;
  
  // find all files in this folder
  std::vector<std::string> trainer_ids;

  DIR *dir;
  struct dirent *ent;
  if ( (dir = opendir ( folder.c_str() ) )  != NULL )
    {
      /* print all the files and directories within directory */
      while ( (ent = readdir (dir)) != NULL) {
	std::string fname = ent->d_name;
	if (  fname != "." && fname != ".." ) 
	  trainer_ids.push_back( fname );
      }
      closedir (dir);
     }
  else
    {
      Helper::halt( "could not open directory " + folder );      
    }
  
  //
  // for primary trainers only (! read_psd ) track H2 and H3 distributions
  //

  Data::Matrix<double> h2m( trainer_ids.size() , suds_t::ns );
  Data::Matrix<double> h2sd( trainer_ids.size() , suds_t::ns );
  Data::Matrix<double> h3m( trainer_ids.size() , suds_t::ns );
  Data::Matrix<double> h3sd( trainer_ids.size() , suds_t::ns );

  
  //
  // load each
  //
  
  for ( int i=0; i<trainer_ids.size() ; i++)
    {
      suds_indiv_t trainer;
      
      trainer.reload( folder + globals::folder_delimiter + trainer_ids[i] , read_psd );      
      
      trainer.fit_lda();

      b->insert( trainer );

      if ( ! read_psd ) 
	{
	  for (int s=0;s<suds_t::ns;s++)
	    {
	      h2m(i,s) = trainer.mean_h2[s] ;
	      h2sd(i,s) = trainer.sd_h2[s] ;
	      h3m(i,s) = trainer.mean_h3[s] ;
	      h3sd(i,s) = trainer.sd_h3[s] ;
	    }
	}
          
    }
  
  logger << "  read " << b->size() << " trainers ("
	 << ( read_psd ? "with spectra" : "w/out spectra" )
	 << ") from " << folder << "\n";


  //
  // from primary trainers only, track signal-wise Hjorth limits
  //
  

  if ( ! read_psd )
    {
      suds_t::lwr_h2.resize( suds_t::ns );
      suds_t::upr_h2.resize( suds_t::ns );
      suds_t::lwr_h3.resize( suds_t::ns );
      suds_t::upr_h3.resize( suds_t::ns );
      
      Data::Vector<double> mean_h2m = Statistics::mean( h2m );
      Data::Vector<double> mean_h2sd = Statistics::mean( h2sd );
      Data::Vector<double> mean_h3m = Statistics::mean( h3m );
      Data::Vector<double> mean_h3sd = Statistics::mean( h3sd );

      for (int s=0;s<suds_t::ns;s++)
	{
	  suds_t::lwr_h2[s] = mean_h2m[s] - 3 * mean_h2sd[s];
	  suds_t::upr_h2[s] = mean_h2m[s] + 3 * mean_h2sd[s];

	  suds_t::lwr_h3[s] = mean_h3m[s] - 3 * mean_h3sd[s];
	  suds_t::upr_h3[s] = mean_h3m[s] + 3 * mean_h3sd[s];

	  if ( suds_t::lwr_h2[s] < 0 ) suds_t::lwr_h2[s] = 0;
	  if ( suds_t::lwr_h3[s] < 0 ) suds_t::lwr_h3[s] = 0;
	  
	  logger << "  thresholding " << suds_t::siglab[s] 
		 << " on H2: " << suds_t::lwr_h2[s] << " - " << suds_t::upr_h2[s] 
		 << " and H3: " << suds_t::lwr_h3[s] << " - " << suds_t::upr_h3[s] << "\n";
	}
    }

}



// fit LDA, io.e. after reloading U

void suds_indiv_t::fit_lda()
{

  lda_t lda( y , U );      

  model = lda.fit();

}


//
// make predictions given a different individual's signal data
//

lda_posteriors_t suds_indiv_t::predict( const suds_indiv_t & trainer )
{

  //
  // Project target (this) into trainer space:   U_targ = X_targ * V_trainer * D_trainer^{-1} 
  // subsetting to # of columns
  //

  if ( trainer.W.size() != suds_t::nc || trainer.V.dim2() != suds_t::nc )
    Helper::halt( "V of incorrect column dimension in suds_indiv_t::predict()");
  
  Data::Matrix<double> trainer_DW( suds_t::nc , suds_t::nc );  
  for (int i=0;i< suds_t::nc; i++)
    trainer_DW(i,i) = 1.0 / trainer.W[i];
  
  U_projected = PSD * trainer.V * trainer_DW;
  
  //
  // smooth U (projected)
  //
  
  if ( suds_t::denoise_fac > 0 ) 
    for (int j=0;j<suds_t::nc;j++)
      {
	std::vector<double> * col = U_projected.col_nonconst_pointer(j)->data_nonconst_pointer();
	double sd = MiscMath::sdev( *col );
	double lambda = suds_t::denoise_fac * sd;
	dsptools::TV1D_denoise( *col , lambda );
      }


  //
  // predict using trainer model
  //

  lda_posteriors_t pp = lda_t::predict( trainer.model , U_projected );

  return pp;
}




//
// Primary scoring routine
//

void suds_t::score( edf_t & edf , param_t & param ) {


  //
  // by this point, bank will be populated with N+ trainers
  //
  

  //
  // create a target 
  //

  suds_indiv_t target( edf.id ) ;

  target.proc( edf , param );


  //
  // Do we have prior staging available for this target?
  //

  bool prior_staging = target.obs_stage.size() != 0 ;

  //
  // for weight training, on use 'self' 
  //

  bool retrain_self = ! param.has( "retrain-all" );
  
  //
  // save weights for each trainer, based on re-predicting
  //

  Data::Vector<double> wgt_max( bank.size() ) ;
  Data::Vector<double> wgt_mean( bank.size() ) ;
  Data::Vector<double> wgt_n50( bank.size() ) ;

  //
  // Store Kappa3 for each trainer (valid w/ prior staging only)
  //
  
  Data::Vector<double> k3_prior( bank.size() ) ;
  
  //
  // Stats on trainers
  //

  Data::Vector<double> nr_trainer( bank.size() ) ; // number of unique imputed stages 

  //
  // Stats on weight trainers
  //

  std::map<std::string,double> wtrainer_mean_k3;
  

  //
  // iterate over trainers
  //
  
  std::set<suds_indiv_t>::const_iterator tt = bank.begin();

  int cntr = 0;

  while ( tt != bank.end() )
    {
      
      if ( (cntr+1) % 50 == 0 ) logger << "   ... " << (cntr+1) << "/" << bank.size() << " trainers\n";

      //
      // Extract this one trainer
      //

      const suds_indiv_t & trainer = *tt;

      
      //
      // Verbose output: TRAINER x WEIGHT TRAINER STATS: not yet
      //

      
      //
      // Predict target given trainer, after project target PSD into trainer-defined space 
      // ( i.e. this generates target.U_projected based on trainer, and then uses it to 
      //        predict target class given the trainer model )
      //
      
      lda_posteriors_t prediction = target.predict( trainer );

      
      //
      // Save predictions
      //

      target.prd_stage = suds_t::type( prediction.cl );   

      target.add( tt->id , prediction );
      

      //
      // Reweighting (using individuals specified in wbank, if any) 
      //
      // Consider that target's predicted stages (from this one particular trainer)
      // are in fact the real/observed stages for this target.    Now, the 'target'
      // stages and target model is used to predict other people (i.e. called 'weight trainers', 
      // and they are effectively targets in this context)
      //
      // Requires at least 2 predicted stages (of sufficient N) have been predicted by the trainer
      // before doing this step 
 
      //  trainer --> target                         : using trainer model (to define U)
      //              target ---> weight trainer1    : using target model (to define U)
      //              target ---> weight trainer2
      //              target ---> weight trainer3


      //
      // Now consider how well this predicts all the weight-trainers
      // i.e. where we also have true stage information
      //
      
      double max_kappa = 0;
      double mean_kappa = 0;
      int n_kappa50 = 0;
      int n_kappa_all = 0;
      
      std::map<std::string,int> counts;
      for (int i=0;i<prediction.cl.size();i++) counts[ prediction.cl[ i ] ]++;
      int nr= 0 ; 
      std::map<std::string,int>::const_iterator cc = counts.begin();
      while ( cc != counts.end() )
	{
	  if ( cc->second >= suds_t::required_epoch_n ) ++nr;
	  ++cc;
	}

      // save for output
      nr_trainer[ cntr ] = nr;


      //
      // If prior staging is available, report on kappa for this single trainer
      //

      double k3;
      if ( prior_staging )
	{
	  // obs_stage for valid epochs only
	  double kappa3 =  MiscMath::kappa( NRW( str( target.prd_stage ) ) , NRW( str( target.obs_stage_valid ) ) );	  
	  k3_prior[ cntr ] = kappa3;
	  
	  // tmp
	  k3 = kappa3;
	  //	  suds_t::tabulate( NRW( str( target.prd_stage ) ) , NRW( str( target.obs_stage_valid ) ) , true  );	  
	}
      
      
      
      bool okay_to_fit_model = nr > 1;

      if ( okay_to_fit_model )
	{
	  
	  // std::cout << "U\n" << target.U.print() << "\n";
	  // for (int i=0;i<prediction.cl.size();i++) std::cout << prediction.cl[i] << "\n";
	  // std::cout << "\n";
	  	  
	  //
	  // Generate model for prediction based on 'dummy' target (imputed) stages
	  // but U basd on the target's own SVD (i.e. not projected into trainer space);  
	  // This we use target.U, which is the original for the target, based on their own data
	  // (we ignore the U_projected which is based on the trainer model)
	  //

	  lda_t lda( prediction.cl , target.U );
      
	  target.model = lda.fit();
	  
	  std::set<suds_indiv_t>::iterator ww = wbank.begin();
	  while ( ww != wbank.end() )
	    {
	      
	      suds_indiv_t & weight_trainer = (suds_indiv_t&)(*ww);
	      
	      // only use self-training
	      if ( retrain_self )
		{
		  if ( trainer.id != weight_trainer.id ) { ++ww; continue; } 
		  std::cout << "WEIGHT TRAINER " << weight_trainer.id << "\n";
		}
	      
	      lda_posteriors_t reprediction = weight_trainer.predict( target );
	      
	      weight_trainer.prd_stage = suds_t::type( reprediction.cl );
	      
	      // obs_stage for predicted/valid epochs only
	      double kappa = MiscMath::kappa( NRW( reprediction.cl ) , NRW( str( weight_trainer.obs_stage ) ) );
	      
	      //logger << "trainer/wt kappa = " << trainer.id << " " << weight_trainer.id << " " << k3 << "\t" << kappa << "\n";
	      //suds_t::tabulate(  NRW( reprediction.cl ) , NRW(str( weight_trainer.obs_stage ) ) , true );
	      
	      ++n_kappa_all;
	      if ( kappa > 0.5 ) n_kappa50++;
	      if ( kappa > max_kappa ) max_kappa = kappa;
	      mean_kappa +=  kappa  ;
	      
	      if ( suds_t::verbose ) 
		wtrainer_mean_k3[ weight_trainer.id ] += kappa;
	      
	      ++ww;
	    }
	  
	}

          
          
     
      //
      // Trainer weights
      //
      
      if ( wbank.size() > 0 && okay_to_fit_model ) 
	{
	  wgt_max[ cntr ] = max_kappa;
	  wgt_mean[ cntr ] = ( mean_kappa ) / (double)n_kappa_all ;
	  wgt_n50[ cntr ] = n_kappa50;
	}
      
      //
      // Next trainer
      //

      ++cntr;
      ++tt;

    }



  //
  // Derive weights for each trainer based on KL divergence from trainer stage distribition to the mean
  // over all trainers
  //

  // normalized from 0..1 

  Data::Vector<double> wgt_kl = Statistics::unit_scale( target.wgt_kl() );
  
  //
  // Output all weights, and generate 'final' weigth
  //
  
  Data::Vector<double> wgt( bank.size() );

  tt = bank.begin();
  cntr = 0;

  while ( tt != bank.end() )
    {

      const suds_indiv_t & trainer = *tt;

      writer.level( trainer.id , "TRAINER" );

      writer.value( "NS" , nr_trainer[ cntr ] );

      if ( prior_staging )
	writer.value( "K3" , k3_prior[ cntr ] );
      
      writer.value( "WGT_KL"   , wgt_kl[ cntr ] );

      if ( wbank.size() > 0 ) 
	{
	  writer.value( "WGT_N50"  , wgt_n50[ cntr ] );
	  writer.value( "WGT_MAX"  , wgt_max[ cntr ] );
	  writer.value( "WGT_MEAN" , wgt_mean[ cntr ] ); // normalized	  
	}

      //
      // define 'final' weight: if weight trainers exist, 
      // using WGT_MEAN, otherwise WGT_KL
      //

      if ( wbank.size() > 0 )
	wgt[ cntr ] = wgt_mean[ cntr ] ;
      else
	wgt[ cntr ] = wgt_kl[ cntr ] ;
      
      ++tt;
      ++cntr;
    }
  writer.unlevel( "TRAINER" );
  

  //
  // Verbose output: mean weight trainer values
  //

  if ( suds_t::verbose && wbank.size() > 0 )
    {
      std::set<suds_indiv_t>::iterator ww = wbank.begin();
      while ( ww != wbank.end() )
	{	  
	  suds_indiv_t & weight_trainer = (suds_indiv_t&)(*ww);
	  writer.level( weight_trainer.id , "WTRAINER" );
	  double m = wtrainer_mean_k3[ weight_trainer.id ] / (double)bank.size() ;
	  writer.value( "K3" , m );
	  ++ww;
	}
      writer.unlevel( "WTRAINER" );
            
    }


  //
  // Normalize wgt / truncate at percentile?
  //

  if ( suds_t::wgt_percentile > 0 ) 
    {

      // get value X = top N% and set to 0/1 if below/above X
      double threshold = MiscMath::percentile( *wgt.data_pointer() , 1.0 - suds_t::wgt_percentile / 100.0 ) ;

      // binarize wgt 
      for (int i=0;i<wgt.size();i++)	
	wgt[i] = wgt[i] >= threshold ? 1 : 0 ;
	
    }
  

  //
  // Construct (weighted) posterior probabilities
  //    
  
  const int ne = target.prd_stage.size();

  target.prd_stage.clear();

  target.prd_stage.resize( SUDS_UNKNOWN );
  
  Data::Matrix<double> pp( ne , suds_t::n_stages );

  int ntrainers = 0;
  double tot_wgt = 0;

  std::map<std::string,Data::Matrix<double> >::iterator ii = target.target_posteriors.begin();
  while ( ii != target.target_posteriors.end() )
    {

      // get posteriors from this trainer 

      Data::Matrix<double> & m = ii->second;

      // force 0/1 encoding? i.e. 100% weight placed on most likely
      
      if ( suds_t::use_best_guess ) 
	{
	  suds_t::make01( m );
	}

      const suds_indiv_t & trainer = *bank.find( suds_indiv_t( ii->first ) );
            
      double w = wgt[ ntrainers ];

      tot_wgt += w;
            
      if ( pp.dim1() != m.dim1() || pp.dim2() != m.dim2() )
	Helper::halt( "internal error in compiling posteriors across trainers" );

      for (int i=0;i<ne;i++)
	for (int j=0;j<suds_t::n_stages;j++)
	  pp(i,j) += w * m(i,j);
      
      ++ntrainers;
      ++ii;
    }

  if ( suds_t::wgt_percentile > 0 ) 
    logger << "  constructed posteriors using " << (int)tot_wgt << " (of " << ntrainers << ") trainers\n";
  else
    logger << "  constructed posteriors using " << ntrainers << " trainers\n";


  //
  // Normalize (weighted) posteriors to sum to 1.0, and get MAP
  //

  double mean_maxpp = 0;

  for (int i=0;i<ne;i++)
    {
      // normalize
      for (int j=0;j<suds_t::n_stages;j++) // 5 or 3 stages
        pp(i,j) /= (double)tot_wgt;

      // track level of confidence for MAP
      mean_maxpp += suds_t::maxpp( pp.row(i) );
      
    }
  mean_maxpp /= (double)ne;



  //
  // Report, and calculate stage duraition too
  //

  

  std::map<int,int> e2e;
  for (int i=0;i<ne;i++) e2e[target.epochs[i]] = i ;  
  const int ne_all = edf.timeline.num_epochs();
    
  std::vector<std::string> final_prediction;

  const double epoch_sec = edf.timeline.epoch_length();

  std::map<std::string,double> prd_dur; // sum of PP
  std::map<std::string,double> prd2_dur; // based on mist likely
  std::map<std::string,double> obs_dur; // obserevd  (if present)... but based on same epochs as used in the staging (i.e. removing some outliers) 

  for (int i=0;i<ne_all;i++)
    {

      int e = -1;
      if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
      
      writer.epoch( edf.timeline.display_epoch( i ) );

      if ( e != -1 ) 
	{
	  
	  writer.value( "INC" , 1 );

	  if ( suds_t::n_stages == 5 )
	    {
	      writer.value( "PP_N1"  , pp(e,0) );
	      writer.value( "PP_N2"  , pp(e,1) );
	      writer.value( "PP_N3"  , pp(e,2) );
	      writer.value( "PP_REM" , pp(e,3) );
	      writer.value( "PP_W"   , pp(e,4) );
	      writer.value( "PP_NR"  , pp(e,0)+pp(e,1)+pp(e,2) );
	    }
	  else
	    {
	      writer.value( "PP_NR"  , pp(e,0) );
	      writer.value( "PP_REM" , pp(e,1) );
	      writer.value( "PP_W"   , pp(e,2) );
	    }
	  
	  // most likely value
	  std::string predss = max( pp.row(e) );
	  writer.value( "PRED" , predss );
	  final_prediction.push_back( predss );

	  // track stage duration (based on probabilistic calls)	  
	  if ( suds_t::n_stages == 5 )
	    {
	      prd_dur[ "N1" ]  += pp(e,0) * epoch_sec ;
	      prd_dur[ "N2" ]  += pp(e,1) * epoch_sec ;
	      prd_dur[ "N3" ]  += pp(e,2) * epoch_sec ;
	      prd_dur[ "REM" ] += pp(e,3) * epoch_sec ;
	      prd_dur[ "W" ]   += pp(e,4) * epoch_sec ;
	    }
	  else
	    {
	      prd_dur[ "NR" ]  += pp(e,0) * epoch_sec ;
	      prd_dur[ "REM" ] += pp(e,1) * epoch_sec ;
	      prd_dur[ "W" ]   += pp(e,2) * epoch_sec ;
	    }

	  // duration based on MAP estimate
	  prd2_dur[ predss ] += epoch_sec;

	  if ( prior_staging )
	    {
	      // discordance if prior/obs staging available	      


	      bool disc = predss !=  str( target.obs_stage[i] ) ;
	      writer.value( "DISC" , disc );

	      // collapse 5->3 ?
	      if ( suds_t::n_stages == 5 )
		{
		  bool disc3 = NRW( predss ) != NRW( str( target.obs_stage[i] ) ) ;
		  writer.value( "DISC3" , disc3 );
		}
	      
	      writer.value( "PRIOR" ,  str( target.obs_stage[i] ) );
	      
	      // comparable OBS duration
	      obs_dur[ str( target.obs_stage[i] ) ] += epoch_sec; 

	    }
	}
      else
	{
	  writer.value( "INC" , 0 );

	  // lookup from all stages
	  if ( prior_staging )
	    writer.value( "PRIOR" ,  str( target.obs_stage[i] ) );	  
	}
      
    }

  writer.unepoch();
  

  


  //
  // Baseline/individual level report
  // and cnfusion matrices to console (if prior staging available)
  //
  
  // stage durations (in minutes)

  if ( suds_t::n_stages == 5 )
    {
      writer.value( "DUR_PRD_N1" , prd_dur[ "N1" ] / 60.0 );
      writer.value( "DUR_PRD_N2" , prd_dur[ "N2" ] / 60.0 );
      writer.value( "DUR_PRD_N3" , prd_dur[ "N3" ] / 60.0 );
    }
  else
    {
      writer.value( "DUR_PRD_NR" , prd_dur[ "NR" ] / 60.0 );
    }
  writer.value( "DUR_PRD_REM" , prd_dur[ "REM" ] / 60.0 );
  writer.value( "DUR_PRD_W" , prd_dur[ "W" ] / 60.0 );

  // alternate estimates, based on most likely predicted epoch
  if ( suds_t::n_stages == 5 )
    {
      writer.value( "DUR_PRD2_N1" , prd2_dur[ "N1" ] / 60.0 );
      writer.value( "DUR_PRD2_N2" , prd2_dur[ "N2" ] / 60.0 );
      writer.value( "DUR_PRD2_N3" , prd2_dur[ "N3" ] / 60.0 );
    }
  else
    {
      writer.value( "DUR_PRD2_NR" , prd2_dur[ "NR" ] / 60.0 );
    }
  writer.value( "DUR_PRD2_REM" , prd2_dur[ "REM" ] / 60.0 );
  writer.value( "DUR_PRD2_W" , prd2_dur[ "W" ] / 60.0 );
  
  if ( prior_staging )
    {
      
      logger << std::fixed << std::setprecision(2);
      
      // original reporting (5 or 3 level)
      double kappa = MiscMath::kappa( final_prediction , str( target.obs_stage_valid ) );
      logger << "\n  Confusion matrix: " << suds_t::n_stages << "-level classification: kappa = " << kappa << "\n";
      suds_t::tabulate(  final_prediction , str( target.obs_stage_valid ) , true );
      writer.value( "K" , kappa );
      
      // collapse 5->3?
      if ( suds_t::n_stages == 5 )
	{
	  double kappa3 = MiscMath::kappa( NRW(final_prediction) , NRW(str( target.obs_stage_valid ) ) );
	  logger << "\n  Confusion matrix: 3-level classification: kappa = " << kappa3 << "\n";
	  suds_t::tabulate(  NRW(final_prediction) , NRW(str( target.obs_stage_valid ) ) , true );
	  writer.value( "K3" , kappa3 );
	}
      
      writer.value( "MAXPP" , mean_maxpp );

      //
      // estimates of observed stage duration (based on comparable set of epochs)
      //

      std::map<std::string,double>::const_iterator ss = obs_dur.begin();
      while ( ss != obs_dur.end() )
	{
	  writer.value( "DUR_OBS_" + ss->first , ss->second / 60.0 );
	  ++ss;
	}
      
      //
      // also, given correlations between weights and trainer kappas
      //

      writer.value( "R_K3_KL" , Statistics::correlation( wgt_kl , k3_prior) ); 

      if ( wbank.size() > 0 ) 
	{
	  writer.value( "R_K3_MAX" , Statistics::correlation( wgt_max  , k3_prior) ); 
	  writer.value( "R_K3_MEAN" , Statistics::correlation( wgt_mean , k3_prior) ); 
	  writer.value( "R_K3_N50" , Statistics::correlation( wgt_n50  , k3_prior) ); 
	}
					          
    }


  //
  // Verbose output?
  //

  if ( suds_t::mat_dump_file != "" ) 
    {
      
      std::string filename = Helper::expand( suds_t::mat_dump_file );
      std::ofstream OUT1( filename.c_str() , std::ios::out );
      
      logger << "  writing epoch-wise matrix to " << filename << "\n";

      OUT1 << "ID\tE";

      for (int i=0;i<target.PSD.dim2();i++)
	OUT1 << "\t" << "X" << (i+1);

      for (int i=0;i<target.U.dim2();i++)
	OUT1 << "\t" << "PSC" << (i+1);

      if ( suds_t::n_stages == 5 )
	{
	  OUT1 << "\tPP_N1"
	       << "\tPP_N2"
	       << "\tPP_N3"
	       << "\tPP_R"
	       << "\tPP_W";      
	}
      else
	{
	  OUT1 << "\tPP_NR"
	       << "\tPP_R"
	       << "\tPP_W";      
	}
      
      OUT1 << "\tPRD";

      if ( prior_staging ) OUT1 << "\tOBS";
      OUT1 << "\n";
      
      // each row/epoch
      for (int i=0;i<ne_all; i++)
	{

	  int e = -1;
	  if ( e2e.find( i ) != e2e.end() ) e = e2e[i];	  
	  
	  if ( e == -1 ) continue;

	  // only display good lines
	  OUT1 << target.id << "\t"
	       << edf.timeline.display_epoch( i ) ;
	  
	  for (int j=0;j<target.PSD.dim2();j++)
	    OUT1 << "\t" << target.PSD(e,j); 
	  
	  for (int j=0;j<target.U.dim2();j++)
	    OUT1 << "\t" << target.U(e,j); 
	  
	  for (int j=0;j<pp.dim2();j++)
	    OUT1 << "\t" << pp(e,j); 
	  
	  OUT1 << "\t" << final_prediction[e] ;

	  if ( prior_staging ) 
	    OUT1 << "\t" << str( target.obs_stage[i] ) ;

	  OUT1 << "\n";

	}

      OUT1.close();
    }




  //
  // Write .eannot file?
  //
  
  if ( suds_t::eannot_file != "" )
    {
      // expecting this will have individual wild-cards, but these will have 
      // been expanded already;  this is for home folder encoding ~ 
      std::string filename = Helper::expand( suds_t::eannot_file );

      logger << "\n  writing .eannot stage annotations "
	     << ( suds_t::eannot_ints ? "(as integeres) " : "" )
	     << " to " << filename << "\n";

      std::ofstream OUT1( filename.c_str() , std::ios::out );


      // make sure we output all epochs
      for (int i=0;i<ne_all;i++)
	{
	  
	  int e = -1;
	  if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
	  if ( e != -1 )
	    {
	      if ( suds_t::eannot_ints )
		OUT1 << suds_t::num( final_prediction[e] ) << "\n";
	      else
		OUT1 << final_prediction[e] << "\n";	      
	    }
	  else // could not score
	    {
	      if ( suds_t::eannot_ints ) OUT1 << suds_t::num( "?" ) << "\n";
	      else OUT1 << "?" << "\n";
	    }
	}
            
      OUT1.close();      
    }
  
}




void suds_indiv_t::add( const std::string & trainer_id , const lda_posteriors_t & prediction )
{
  
  target_posteriors[ trainer_id ] = prediction.pp;

  target_predictions[ trainer_id ] = suds_t::type( prediction.cl );
  
}



std::map<std::string,std::map<std::string,int> > suds_t::tabulate( const std::vector<std::string> & a , 
								   const std::vector<std::string> & b , 
								   const bool print  )
{
  std::map<std::string,std::map<std::string,int> > res;
  
  const int n = a.size();
  
  if ( n != b.size() ) 
    Helper::halt( "internal error: unequal vectors in tabulate()" );
  
  std::set<std::string> uniq;
  for (int i=0;i<n;i++)
    {
      res[ a[i] ][ b[i] ]++;
      uniq.insert( a[i] );
      uniq.insert( b[i] );
    }

  std::map<std::string,double> rows, cols;
  double tot = 0;
  std::set<std::string>::const_iterator uu = uniq.begin();
  while ( uu != uniq.end() )
    {
      std::set<std::string>::const_iterator jj = uniq.begin();
      while ( jj != uniq.end() )
	{
	  if ( res.find( *uu ) == res.end() )
	    res[ *uu ][ *jj ] = 0;
	  else
	    {
	      std::map<std::string,int> & rjj = res.find(*uu)->second;
	      if ( rjj.find( *jj ) == rjj.end() )
		res[ *uu ][ *jj ] = 0;
	    }	      

	  // col/row marginals
	  rows[ *uu ] += res[ *uu ][ *jj ] ;
	  cols[ *jj ] += res[ *uu ][ *jj ] ;
	  tot += res[ *uu ][ *jj ];
	  
	  ++jj;
	}
      ++uu;
    }
  
  
  if ( print )
    {

      logger << "\t  Obs:\n\t";
      std::set<std::string>::const_iterator uu = uniq.begin();
      while ( uu != uniq.end() )
	{
	  logger << "\t" << *uu;
	  ++uu;
	}
      logger << "\tTot\n";	
      
      logger << "  Pred:";
      uu = uniq.begin();
      while ( uu != uniq.end() )
	{
	  logger << "\t" << *uu;
	  std::set<std::string>::const_iterator jj = uniq.begin();
	  while ( jj != uniq.end() )
	    {
	      logger << "\t" << res[ *uu ][ *jj ];
	      ++jj;
	    }
	  // row sums
	  logger << "\t" << rows[ *uu ]/tot;
	  logger << "\n";
	  ++uu;
	}
      // col sums
      logger << "\tTot:";
      std::set<std::string>::const_iterator jj = uniq.begin();
      while ( jj != uniq.end() )
	{
	  logger << "\t" << cols[ *jj ]/tot;
	  ++jj;
	}
      logger << "\t1.00\n";
    }

  
  return res;
}



Data::Vector<double> suds_indiv_t::wgt_kl() const { 

  // returned weights
  const int nt = target_predictions.size();  

  Data::Vector<double> W( nt );

  if ( nt == 0 ) return W;

  Data::Matrix<double> Q( nt , suds_t::n_stages ) ;  

  int r = 0;
  std::map<std::string,std::vector<suds_stage_t> >::const_iterator ii = target_predictions.begin();
  while ( ii != target_predictions.end() ) 
    {

      const double ne = ii->second.size();

      if ( suds_t::n_stages == 5 ) 
	for (int e=0; e<ne; e++) 
	  {	    
	    if      ( ii->second[e] == SUDS_N1 ) Q(r,0)++;
	    else if ( ii->second[e] == SUDS_N2 ) Q(r,1)++;
	    else if ( ii->second[e] == SUDS_N3 ) Q(r,2)++;
	    else if ( ii->second[e] == SUDS_REM ) Q(r,3)++;
	    else if ( ii->second[e] == SUDS_WAKE ) Q(r,4)++;
	  }
      else
	for (int e=0; e<ne; e++) 
	  {	    
	    if      ( ii->second[e] == SUDS_NR ) Q(r,0)++;
	    else if ( ii->second[e] == SUDS_REM ) Q(r,1)++;
	    else if ( ii->second[e] == SUDS_WAKE ) Q(r,2)++;
	  }
	

      // normalize
      for (int s=0;s<suds_t::n_stages;s++) Q(r,s) /= ne;
      
      // next trainer
      ++ii;
      ++r;
    }

  
  // Means

  Data::Vector<double> P = Statistics::mean( Q );

  // divergence for each trainer from the mean
  r = 0;
  ii = target_predictions.begin();
  while ( ii != target_predictions.end() ) 
    {      
      // negative KL
      double ss = 0;
      for ( int s = 0 ; s < suds_t::n_stages ; s++ )
	if ( Q(r,s) > 0 ) ss += P[s] * log( P[s] / Q(r,s) );  
      W[r] = -ss;
      
      ++r;
      ++ii;
    }
    
  return W;
}


