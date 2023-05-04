
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

#include "dsp/ged.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "db/db.h"
#include "helper/logger.h"
#include "dsp/ngaus.h"
#include "stats/eigen_ops.h"

extern writer_t writer;
extern logger_t logger;

void ged_wrapper( edf_t & edf , param_t & param )
{

  // signals

  const bool NO_ANNOTS = true;
  signal_list_t signals = edf.header.signal_list(  param.requires( "sig" ) , NO_ANNOTS );
  const int ns = signals.size();
  if ( ns < 2 ) return;  
  std::vector<double> Fs = edf.header.sampling_freq( signals );
  const int sr = Fs[0]; // these are checked to be the same by slice
   
  // get data
  // always run in whole-signal mode

  eigen_matslice_t mslice( edf , signals , edf.timeline.wholetrace() );
  
  Eigen::MatrixXd X = mslice.data_ref();
  
  // run-modes
  //  1 : whole signal covariaces, compare narrowband filter signals to broadband
  //    S - narrow-band filtered
  //    R - broad-band filtered
  //    --> generates a time-series (of top N) components
  //    --> for each, does a secondary GED of peak|trough|rise|all|broadband vs peak|trough|rise|all|broadband, with w1, and w1 seconds around each point

  // 2: point based GED but from cache
  //     cache = points  
  //     

  
  int run_mode = 1;
  
  if ( run_mode == 1 ) 
    {
      ged_runmode1( edf, param , X , sr ); 
      return;
    }


  if ( run_mode == 2 )
    {
      const std::vector<uint64_t> * tp = mslice.ptimepoints();
      
      ged_runmode2( edf , param, X , tp , sr );
    }
  
}



void ged_runmode2( edf_t & edf , param_t & param, Eigen::MatrixXd & Rd , const std::vector<uint64_t> * tp,  int sr )
{
  
  // options

  // S time-domain restriction:
  const std::string a1 = param.requires( "a1" );
  const double w1 = param.has( "w1" ) ? param.requires_dbl( "w1" ) : 0 ;
  const bool x1 = param.has( "x1" ) ;

  // R time-domain restriction (if any)
  const bool refall = ! param.has( "a2" );
  const std::string a2 = refall ? "" : param.value( "a2" );
  const double w2 = param.has( "w2" ) ? param.requires_dbl( "w2" ) : 0 ;
  const bool x2 = param.has( "x2" ) ;
  
  // get annotations
  annot_t * annot1 = edf.timeline.annotations.find( a1 );
  if ( annot1 == NULL ) Helper::halt( "could not find annotation " + a1 );
  
  annot_t * annot2 = NULL;
  if ( ! refall )
    {
      annot2 = edf.timeline.annotations.find( a2 );
      if ( annot2 == NULL ) Helper::halt( "could not find annotation " + a2 );
    }

  //
  // From original matrix Rd, get the two S and R
  //

  // S covariance 

  Eigen::MatrixXd Sd = eigen_ops::subset_rows( Rd , tp , annot1 , w1 , x1 );
  logger << "  reduced S matrix to " <<	Sd.rows() << " from " << Rd.rows() << "\n";
  Eigen::MatrixXd S = eigen_ops::covariance( Sd );

  // R covariance
  Eigen::MatrixXd R;
  if ( refall )
    {
      R = eigen_ops::covariance( Rd );
    }
  else
    {
      Eigen::MatrixXd Rd2 = eigen_ops::subset_rows( Rd , tp , annot2 , w2 , x2 );
      logger << "  reduced R matrix to " << Rd2.rows() << " from " << Rd.rows() << "\n";
      R = eigen_ops::covariance( Rd2 );
    }
  
  //
  // GED
  //

  ged_t ged;
  ged.covar( S, R );
  ged.calc();


  //
  // Spatial map for S
  //

  int mxch;
  
  Eigen::VectorXd map1 = ged.map( ged.largest_idx , S , &mxch );

  std::cout << "map\n" << map1 << "\n";
  
  //
  // New time series
  //

  const std::string new_ts = param.has( "ts" ) ? param.value( "ts" ) : "" ;
  
  
}


void ged_runmode1( edf_t & edf , param_t & param, Eigen::MatrixXd & Rd , int sr )
{

  // options

  const double f1 = param.requires_dbl( "f1" );

  const double fwhm1 = param.requires_dbl( "fwhm1" );

  const double filteR = param.has( "f2" );
  
  const double f2 = filteR ? param.requires_dbl( "f2" ) : 0 ;

  const double fwhm2  = filteR ? param.requires_dbl( "fwhm2" ) : 0 ;

  const std::string new_ts = param.has( "ts" ) ? param.value( "ts" ) : "" ;

  // compute
  
  Eigen::MatrixXd Sd = Rd;

  const int ns = Sd.cols();

  logger << "  creating narrowband S, " << f1 << " Hz (" << fwhm1 << " FWHM Gaussian)\n";

  for (int s=0; s<ns; s++)
    Sd.col(s) = narrow_gaussian_t::filter( Sd.col(s) , sr, f1, fwhm1 );

  // narrow-band covariance
  Eigen::MatrixXd S = eigen_ops::covariance( Sd );
  
  // reference (broad-band or a different narrow-band) covariance
  if ( filteR ) 
    {
      logger << "  creating narrowband R, " << f2 << " Hz (" << fwhm2 << " FWHM Gaussian)\n";
      for (int s=0; s<ns; s++)
	Rd.col(s) = narrow_gaussian_t::filter( Rd.col(s) , sr, f2, fwhm2 );
    }

  Eigen::MatrixXd R = eigen_ops::covariance( Rd );
  
  // GED
  ged_t ged;
  ged.covar( S, R );
  ged.calc();

  // ged.W = eigenvectors
  // ged.L = eigenvalues
  // ged.largest_idx = idx of largest eigenvalue

  
  // spatial map for largest component of narrow-band S
  int mxch; // channel w/ largest component
  Eigen::VectorXd spatialMap = ged.map( ged.largest_idx , S , &mxch ); 
  
  // narrow-band time series component
  // mxch used to flip polarity of time series (to +ve corr w/ that channel) 
  Eigen::VectorXd ts = ged.time_series( ged.largest_idx, Sd , mxch );
  
  // add new time signal?
  if ( new_ts != "" )
    {
      std::vector<double> copy( ts.size() );
      Eigen::VectorXd::Map( &copy[0], ts.size() ) = ts;
      logger << "  adding channel " << new_ts << " with the new narrow-band time-series\n";
      edf.add_signal( new_ts, sr, copy );
    }

  // output map
  const bool NO_ANNOTS = true;
  signal_list_t signals = edf.header.signal_list(  param.requires( "sig" ) , NO_ANNOTS );
  for (int s=0;s<ns; s++)
    {
      writer.level( signals.label(s) , globals::signal_strat );
      writer.value( "W" , spatialMap[s] );
    }
  writer.unlevel( globals::signal_strat );

  // TODO
  // save can >> components for the f1 frequency

  // then identify troughs/peaks in broad band signal and compare against broad band whole signal..
  //  i.e re-trun GED w/ peak-locked covariance etc.
  

// %% identify troughs and get surrounding covariance matrices

// nwin = ceil(EEG.srate/thetafreq/8); % window size is 1/4 cycle (1/8 of either side)

// % find troughs
// troughs = find(diff(sign(diff( thetacomp )))>0)+1;
// troughs(troughs<nwin+1) = [];
// troughs(troughs>EEG.pnts-nwin-1) = [];


// covT = zeros(EEG.nbchan);

// % trough-locked covariance
// for ti=1:length(troughs)
//     tmpdat = EEG.data(:,troughs(ti)-nwin:troughs(ti)+nwin);
//     tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
//     covT   = covT + (tmpdat*tmpdat')/nwin;
// end
// covT = covT./ti;

// %% GED to get gamma peak/trough networks

// [evecs,evals] = eig(covT,bbcov);
// [~,compidx]   = sort(diag(evals)); % max component

// maps    = covT*evecs; % forward model of filter
// gamnet1 = maps(:,compidx(end));

// % fix sign
// [~,idx] = max(abs(gamnet1));
// gamnet1 = gamnet1 * sign(gamnet1(idx));



// % plot component maps for comparison with dipole projections
// figure(1)
// subplot(234), topoplotIndie(thetamap, EEG.chanlocs,'numcontour',0,'electrodes','off','shading','interp');
// subplot(235), topoplotIndie(gamnet1,  EEG.chanlocs,'numcontour',0,'electrodes','off','shading','interp');

// %% get time course and reconstruct topography and sources

// frex    = linspace(10,190,70);
// mvarpac = zeros(size(frex));

// for fi=1:length(frex)

//     % bandpass filter trough component
//     gam1comp = filterFGx(EEG.data,EEG.srate,frex(fi),15)' * evecs(:,compidx(end));

//     % find peaks and troughs and get distances
//     troughs  = find(diff(sign(diff( gam1comp )))>0)+1;
//     peeks    = find(diff(sign(diff( gam1comp )))<0)+1;
//     mvarpac(fi) = mean(gam1comp(peeks)) - mean(gam1comp(troughs));
// end


 
  
  
}



// The generalized eigenvalues and eigenvectors of a matrix pair A and
// B are scalars λ and vectors v such that Av=λBv. If D is a diagonal
// matrix with the eigenvalues on the diagonal, and V is a matrix with
// the eigenvectors as its columns, then AV=BVD. The matrix V is
// almost always invertible, in which case we have A=BVDV−1. This is
// called the generalized eigen-decomposition.

// The generalized eigenvalues and eigenvectors of a matrix pair may
// be complex, even when the matrices are real. Moreover, the
// generalized eigenvalue might be infinite if the matrix B is
// singular. To workaround this difficulty, the eigenvalues are
// provided as a pair of complex α and real β such that: λi=αi/βi. If
// βi is (nearly) zero, then one can consider the well defined left
// eigenvalue μ=βi/αi such that: μiAvi=Bvi, or even μiuTiA=uTiB where
// ui is called the left eigenvector.

// Call the function compute() to compute the generalized eigenvalues
// and eigenvectors of a given matrix pair. Alternatively, you can use
// the GeneralizedEigenSolver(const MatrixType&, const MatrixType&,
// bool) constructor which computes the eigenvalues and eigenvectors
// at construction time. Once the eigenvalue and eigenvectors are
// computed, they can be retrieved with the eigenvalues() and
// eigenvectors() functions.

#include "helper/helper.h"
#include "db/db.h"
#include "edf/edf.h"
#include "eval.h"

extern writer_t writer;
extern logger_t logger;

// get covariance matrices from data 
void ged_t::data( const Eigen::MatrixXd & Sd, const Eigen::MatrixXd & Rd )
{
  if ( S.rows() < 2 || Rd.rows() < 2 )
    Helper::halt( "bad data for ged_t::data()" );
  
  // centered data matrices
  Eigen::MatrixXd Sc = Sd.rowwise() - Sd.colwise().mean();
  Eigen::MatrixXd Rc = Rd.rowwise() - Rd.colwise().mean();
  
  // covariance matrices
  S = ( Sc.adjoint() * Sc ) / double( Sc.rows() - 1 );
  R = ( Rc.adjoint() * Rc ) / double( Rc.rows() - 1 );
  
}
  

void ged_t::calc()
{

  if ( S.rows() == 0 || S.rows() != R.rows() )
    Helper::halt( "bad covar for ged_t::calc()" );

  es.compute( S, R );

  W = es.eigenvectors();
  L = es.eigenvalues();

  // get largest eigenvalue index
  L.maxCoeff( & largest_idx );
    
}


Eigen::VectorXd ged_t::time_series( const int e , const Eigen::MatrixXd & D , const int maxe )
{
  // narrow-band time series component
  Eigen::VectorXd ts = D * W.col(e);
  
  // fix sign of ts by correl w/ EEG : ts[]  and Sd[,maxe] 
  Eigen::MatrixXd D2( ts.size() , 2 );
  D2 << ts, D.col(maxe);
  Eigen::MatrixXd C1 = eigen_ops::covariance( D2 );
  if ( C1(0,1) < 0 ) ts *= -1;
  return ts;
}



Eigen::VectorXd ged_t::map( const int e , const Eigen::MatrixXd & C , int * maxch )
{
  // spatial map for covariance C
  // (sign-fixed to make the largest element positive)                                                                                                                                                                                               
  Eigen::VectorXd spatialMap = (C*W).col( e );  
  spatialMap.cwiseAbs().maxCoeff( maxch );
  if ( spatialMap[*maxch] < 0 ) spatialMap *= -1;
  return spatialMap;

}

