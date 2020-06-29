
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

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

#include "dirent.h"

#include "stats/matrix.h"
#include "stats/statistics.h"
#include "stats/lda.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "dsp/resample.h"
#include "fftw/fftwrap.h"
#include "dsp/tv.h"

extern logger_t logger;

extern writer_t writer;

std::set<suds_indiv_t> suds_t::bank;
int suds_t::nc;
int suds_t::ns;
std::vector<std::string> suds_t::siglab;
std::vector<double> suds_t::lwr;
std::vector<double> suds_t::upr;
std::vector<double> suds_t::fac;
std::vector<int> suds_t::sr;
double suds_t::denoise_fac;
std::vector<double> suds_t::outlier_ths;


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
  proc( edf , param , true );
  
  // save disk
  write( edf , param ); 

}




void suds_indiv_t::proc( edf_t & edf , param_t & param , bool is_trainer )
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

  if ( signals.size() != ns ) Helper::halt( "could not find all channels" );
  
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
  else
    {
      // for targets, manual/prior staging may exist, in which case we'll want to track it for comparison
      edf.timeline.annotations.make_sleep_stage();

      has_prior_staging = edf.timeline.hypnogram.construct( &edf.timeline , false ) ;
      
      if ( has_prior_staging )
	{
	  // total number of epochs does not match?
	  if ( ne != edf.timeline.hypnogram.stages.size() )
	    Helper::halt( "problem extracting stage information for trainer" );
	}
      
    }
  
  // overkill to have two staging classifications, but keep for now,
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
	  if ( obs_stage[ss] == SUDS_UNKNOWN ) retained[ss] = false;
	  else ++nge;
	}
    }

  //
  // for target individuals without stag, include all epochs here (whether we have staging info or no)
  //

  if ( ! trainer )
    {
      nge = ne;
    }

  
  //
  // iterate over (retained) epochs
  //
    
  int en = 0;
  
  edf.timeline.first_epoch();
  
  while ( 1 ) 
    {
      
      // select epoch
      int epoch = edf.timeline.next_epoch();      	  
      if ( epoch == -1 ) break;

      if ( en == ne ) Helper::halt( "internal error: over-counted epochs" );

      // retained? if not, skip
      if ( ! retained[ en ] ) continue;
      
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
		  if ( en == 0 ) firstrow.push_back(  10*log10( bin.bspec[i] ) );
		  else PSD( en , col ) = 10*log10( bin.bspec[i] ) ; 		  
		  ++col;		  
		}	      
	    }
	  
	} // next signal

      // store/shape output if first go around
      nbins = col;
      if ( en == 0 )
	{
	  PSD.resize( nge , col );
	  for (int i=0;i<col;i++) PSD(0,i) = firstrow[i];
	}
      
      // increase epoch-number
      ++en;

    } // next epoch
  

  // all done: check

  if ( en != nge ) Helper::halt( "internal error: umder-counted epochs" );


  //
  // Get PSC
  //

  // mean-center columns

  U = PSD;
  
  Statistics::mean_center_cols( U );

  // get PSCs

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
		  if ( x[c] < lwr || x[c] > upr ) valid[i] = false;
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
  
  //
  // Remove bad epochs and repeat (SVD and smoothing)
  //

  // nve = number of valid epochs ( ne > nge > nve ) 

  nve = included;

  Data::Matrix<double> PSD2 = PSD;
  PSD.clear();
  PSD.resize( nve , nbins );
  int r = 0;
  for (int i=0;i<PSD2.dim1() ; i++)
    {
      if ( valid[i] )
	{
	  for (int j=0;j<nbins;j++)
	    PSD(r,j) = PSD2(i,j);
	  ++r;
	}
    }
  

  //
  // Get PSC (post outlier removal)
  //

  // mean-center columns (PSD)
  
  Statistics::mean_center_cols( PSD );
  
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
  // Rescale PSCs?  (e.g. unit variance?)
  //

  
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
  // Make variables for LDA (trainers only)
  //

  Data::Matrix<double> U2 = U;
  U.clear();
  U.resize( nve , nc );
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      U(i,j) = U2(i,j);
  

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

  std::string filename = folder + globals::folder_delimiter + edf.id;

  logger << "  writing trainer data to " << filename << "\n";
  
  std::ofstream OUT1( filename.c_str() , std::ios::out );
  
  OUT1 << edf.id << "\n"
       << nve << "\t"
       << nbins << "\t"
       << ns << "\t"
       << nc << "\n";

  // channels , SR [ for comparability w/ other data ] 

  for (int s=0;s<ns;s++)
    {
      OUT1 << suds_t::siglab[s] << "\t"
	   << suds_t::sr[s] << "\t"
	   << suds_t::lwr[s] << "\t"
	   << suds_t::upr[s] << "\t"
	   << suds_t::fac[s] << "\n";
    }

  // W
  for (int j=0;j<nc;j++)
    OUT1 << W[j] << "\n";

  // V
  for (int i=0;i<nbins;i++)
    for (int j=0;j<nc;j++)
      OUT1 << V(i,j) << "\n";

  // stages
  for (int i=0;i<nve;i++)
    OUT1 << y[i] << "\n";

  // U (to reestimate LDA model upon loading, i.e
  //  to use lda.predict() on target   ; only needs to be nc rather than nbins
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      OUT1 << U(i,j) << "\n";
  
  // PSD (mean-centered PSD)
  // i.e. if this trainer is being used as a 'weight trainer',
  // i.e. will project this individuals raw data into the target space
  for (int i=0;i<nve;i++)
    for (int j=0;j<nbins;j++)
      OUT1 << PSD(i,j) << "\n";
  
  OUT1.close();

  //
  // All done
  //

}



void suds_indiv_t::reload( const std::string & filename , bool load_psd )
{
  
  logger << "  reloading trainer data from " << filename << "\n";
  
  std::ifstream IN1( filename.c_str() , std::ios::in );

  int this_ns, this_nc;

  IN1 >> id 
      >> nve 
      >> nbins 
      >> this_ns 
      >> this_nc ;

  if ( this_nc != suds_t::nc || this_ns != suds_t::ns ) Helper::halt( "prob" );

  const int nc = suds_t::nc;
  const int ns = suds_t::ns;
  
  // siglab.resize( ns );
  // sr.resize( ns );
  // lwr.resize( ns );
  // upr.resize( ns );
  // fac.resize( ns );
  
  for (int s=0;s<ns;s++)
    {
      std::string this_siglab;
      double this_lwr, this_upr;
      int this_sr, this_fac;
      
      // IN1 >> siglab[s] 
      // 	  >> sr[s] 
      // 	  >> lwr[s] 
      // 	  >> upr[s] 
      // 	  >> fac[s] ;

      IN1 >> this_siglab
	  >> this_sr 
	  >> this_lwr
	  >> this_upr
	  >> this_fac ;

      if ( this_siglab != suds_t::siglab[s] ) Helper::halt( "prob" );
      if ( this_sr != suds_t::sr[s] ) Helper::halt( "prob" );
      if ( this_lwr != suds_t::lwr[s] ) Helper::halt( "prob" );
      if ( this_upr != suds_t::upr[s] ) Helper::halt( "prob" );
      if ( this_fac != suds_t::fac[s] ) Helper::halt( "prob" );
      
    }
  
  // check that these equal suds_t values?
  
  // W [ only nc ]
  W.resize( nc );
  for (int j=0;j<nc;j++)
    IN1 >> W[j] ;

  // V [ only nc cols ] 
  V.resize( nbins, nc );
  for (int i=0;i<nbins;i++)
    for (int j=0;j<nc;j++)
      IN1 >> V(i,j);

  // stages
  y.resize( nve );
  for (int i=0;i<nve;i++)
    IN1 >> y[i] ;

  // U (to reestimate LDA model upon loading, i.e
  //  to use lda.predict() on target 
  U.resize( nve , nc );
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      IN1 >> U(i,j) ;	

  
  // PSD (mean-centered PSD)
  // i.e. if this trainer is being used as a 'weight trainer',
  // i.e. will project this individuals raw data into the target space
  
  if ( load_psd )
    {      
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




void suds_t::attach_db( const std::string & folder )
{

  // already done?
  if ( bank.size() > 0 ) return;
  
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
	  {
	    trainer_ids.push_back( fname );
	    //std::cout << " jj [" << fname << "]\n";
	  }
      }
      closedir (dir);
     }
  else
    {
      Helper::halt( "could not open directory " + folder );      
    }

  logger << "  found " << trainer_ids.size() << " trainers in " << folder << "\n";
  
  //
  // load each
  //
  
  for ( int i=0; i<trainer_ids.size() ; i++)
    {
      suds_indiv_t trainer;
      
      trainer.reload( folder + globals::folder_delimiter + trainer_ids[i] , true );      
      
      trainer.fit_lda();

      bank.insert( trainer );
    }
  
  logger << "  read " << bank.size() << " trainers from " << folder << "\n";
      
}



// fit LDA, io.e. after reloading U

void suds_indiv_t::fit_lda()
{

  lda_t lda( y , U );      
  
  model = lda.fit();
  
}


//
// make predictions given a different individuals signal data
//

lda_posteriors_t suds_indiv_t::predict( const suds_indiv_t & trainer )
{

  //
  // Project target (this) into trainer space:   U_targ = X_targ * V_trainer * D_trainer^{-1} 
  //
  
  const int nc = trainer.W.size();  
  
  Data::Matrix<double> trainer_DW( nc , nc );

  for (int i=0;i< nc; i++)
    trainer_DW(i,i) = 1.0 / trainer.W[i];
  
  U = PSD * trainer.V * trainer_DW;
  
  //
  // smooth U
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
  // predict using trainer model
  //
  
  lda_posteriors_t pp = lda_t::predict( trainer.model , U );

  return pp;
}




//
// Primary scoring routine
//

void suds_t::score( edf_t & edf , param_t & param ) {

  // by this point, bank will be populated with N+ trainers

  // create a target 

  suds_indiv_t target;

  target.proc( edf , param );

  //
  // iterate over trainers
  //
  
  std::set<suds_indiv_t>::const_iterator tt = bank.begin();
  
  while ( tt != bank.end() )
    {

      const suds_indiv_t & trainer = *tt;
      
      lda_posteriors_t prediction = target.predict( trainer );

      // most likely stage assignment 
      
      target.prd_stage = suds_t::type( prediction.cl );      
            
      ++tt;

    }
  
}

