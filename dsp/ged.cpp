
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
#include "param.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "db/db.h"
#include "helper/logger.h"
#include "dsp/ngaus.h"
#include "dsp/hilbert.h"
#include "stats/eigen_ops.h"
#include "defs/defs.h"
#include "annot/annot.h"
#include "timeline/timeline.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <set>

extern writer_t writer;
extern logger_t logger;


// ---------------------------------------------------------------------------
// ged_input_mode_t helpers
// ---------------------------------------------------------------------------

std::string ged_input_mode_str( ged_input_mode_t m )
{
  if ( m == ged_input_mode_t::NB  ) return "nb";
  if ( m == ged_input_mode_t::ENV ) return "env";
  return "raw";
}

ged_input_mode_t ged_input_mode_from_str( const std::string & s )
{
  if ( s == "nb"  ) return ged_input_mode_t::NB;
  if ( s == "env" ) return ged_input_mode_t::ENV;
  return ged_input_mode_t::RAW;
}


// ---------------------------------------------------------------------------
// ged_t static member definitions
// ---------------------------------------------------------------------------

bool ged_t::has_group_solution = false;
std::string ged_t::group_solution_fname = "";
std::vector<std::string> ged_t::group_channels;
Eigen::MatrixXd ged_t::group_W;
Eigen::VectorXd ged_t::group_L;
std::vector<double> ged_t::group_ch_means;
std::vector<double> ged_t::group_ch_sds;
ged_input_mode_t ged_t::group_input_mode = ged_input_mode_t::RAW;
double ged_t::group_nb_f    = 0.0;
double ged_t::group_nb_fwhm = 0.0;
double ged_t::group_env_lwr = 0.0;
double ged_t::group_env_upr = 0.0;
bool ged_t::group_z_scored  = false;

void ged_t::clear_group_solution()
{
  has_group_solution = false;
  group_solution_fname = "";
  group_channels.clear();
  group_W.resize(0,0);
  group_L.resize(0);
  group_ch_means.clear();
  group_ch_sds.clear();
  group_nb_f = group_nb_fwhm = group_env_lwr = group_env_upr = 0.0;
  group_z_scored = false;
}


// ---------------------------------------------------------------------------
// ged_t::data() — fixed (was checking S.rows() before S assigned)
// ---------------------------------------------------------------------------

void ged_t::data( const Eigen::MatrixXd & Sd , const Eigen::MatrixXd & Rd )
{
  if ( Sd.rows() < 2 || Rd.rows() < 2 )
    Helper::halt( "ged_t::data(): fewer than 2 rows in data matrix" );

  Eigen::MatrixXd Sc = Sd.rowwise() - Sd.colwise().mean();
  Eigen::MatrixXd Rc = Rd.rowwise() - Rd.colwise().mean();

  S = ( Sc.adjoint() * Sc ) / double( Sc.rows() - 1 );
  R = ( Rc.adjoint() * Rc ) / double( Rc.rows() - 1 );

  n_S = (int)Sd.rows();
  n_R = (int)Rd.rows();
}


// ---------------------------------------------------------------------------
// ged_t::regularize_R()
// ---------------------------------------------------------------------------

void ged_t::regularize_R( const double alpha )
{
  if ( alpha <= 0.0 ) return;
  const int nc = (int)R.rows();
  if ( nc == 0 ) return;
  const double diag_val = alpha * R.trace() / double(nc);
  R += diag_val * Eigen::MatrixXd::Identity( nc , nc );
}


// ---------------------------------------------------------------------------
// ged_t::calc()
// ---------------------------------------------------------------------------

void ged_t::calc()
{
  solved = false;

  if ( S.rows() == 0 || S.rows() != R.rows() )
    Helper::halt( "ged_t::calc(): bad covariance matrices" );

  es.compute( S , R );

  if ( es.info() != Eigen::Success )
    {
      logger << "  ** WARNING: GED solver failed (R may be singular; try reg=0.01 or larger)\n";
      return;
    }

  W = es.eigenvectors();
  L = es.eigenvalues();

  L.maxCoeff( &largest_idx );

  solved = true;
}


// ---------------------------------------------------------------------------
// ged_t::sorted_indices() — descending eigenvalue order
// ---------------------------------------------------------------------------

std::vector<int> ged_t::sorted_indices() const
{
  const int nc = (int)L.size();
  std::vector<int> idx( nc );
  std::iota( idx.begin() , idx.end() , 0 );
  std::sort( idx.begin() , idx.end() ,
             [&]( int a , int b ) { return L[a] > L[b]; } );
  return idx;
}


// ---------------------------------------------------------------------------
// ged_t::eigenvalue_ratio()
// ---------------------------------------------------------------------------

double ged_t::eigenvalue_ratio( const int e ) const
{
  const double total = L.sum();
  if ( total <= 0.0 ) return 0.0;
  return L[e] / total;
}


// ---------------------------------------------------------------------------
// ged_t::map() — forward model
// ---------------------------------------------------------------------------

Eigen::VectorXd ged_t::map( const int e , const Eigen::MatrixXd & C , int * maxch )
{
  Eigen::VectorXd spatialMap = ( C * W ).col( e );
  spatialMap.cwiseAbs().maxCoeff( maxch );
  if ( spatialMap[ *maxch ] < 0 ) spatialMap *= -1;
  return spatialMap;
}


// ---------------------------------------------------------------------------
// ged_t::time_series() — projection
// ---------------------------------------------------------------------------

Eigen::VectorXd ged_t::time_series( const int e , const Eigen::MatrixXd & D , const int maxch )
{
  Eigen::VectorXd ts = D * W.col( e );

  // sign-fix: ensure positive correlation with dominant channel
  Eigen::MatrixXd D2( ts.size() , 2 );
  D2.col(0) = ts;
  D2.col(1) = D.col( maxch );
  Eigen::MatrixXd C1 = eigen_ops::covariance( D2 );
  if ( C1(0,1) < 0 ) ts *= -1;

  return ts;
}


// ---------------------------------------------------------------------------
// ged_t::focality()
// ---------------------------------------------------------------------------

double ged_t::focality( const Eigen::VectorXd & spatialMap )
{
  const int nc = (int)spatialMap.size();
  if ( nc < 2 ) return 1.0;

  Eigen::VectorXd w = spatialMap.cwiseAbs();
  const double total = w.sum();
  if ( total <= 0.0 ) return 0.0;
  w /= total;

  // entropy
  double H = 0.0;
  for ( int i = 0; i < nc; i++ )
    if ( w[i] > 0.0 ) H -= w[i] * std::log( w[i] );

  const double H_max = std::log( double(nc) );
  if ( H_max <= 0.0 ) return 1.0;
  return 1.0 - H / H_max;
}


// ---------------------------------------------------------------------------
// ged_t::ap_index() — weighted mean AP coordinate
// ---------------------------------------------------------------------------

double ged_t::ap_index( const Eigen::VectorXd & spatialMap ,
                        const std::vector<std::string> & ch_names ,
                        const std::map<std::string,double> & chan_ap )
{
  double num = 0.0 , den = 0.0;
  const int nc = (int)spatialMap.size();
  for ( int i = 0; i < nc; i++ )
    {
      std::string uc = Helper::toupper( ch_names[i] );
      auto it = chan_ap.find( uc );
      if ( it == chan_ap.end() ) continue;
      double w = std::abs( spatialMap[i] );
      num += w * it->second;
      den += w;
    }
  if ( den <= 0.0 ) return std::numeric_limits<double>::quiet_NaN();
  return num / den;
}


// ---------------------------------------------------------------------------
// ged_t::lat_index() — (R_wt - L_wt) / (R_wt + L_wt)
// ---------------------------------------------------------------------------

double ged_t::lat_index( const Eigen::VectorXd & spatialMap ,
                         const std::vector<std::string> & ch_names ,
                         const std::map<std::string,int> & chan_lat )
{
  double L_wt = 0.0 , R_wt = 0.0;
  const int nc = (int)spatialMap.size();
  for ( int i = 0; i < nc; i++ )
    {
      std::string uc = Helper::toupper( ch_names[i] );
      auto it = chan_lat.find( uc );
      if ( it == chan_lat.end() ) continue;
      double w = std::abs( spatialMap[i] );
      if ( it->second < 0 ) L_wt += w;
      else if ( it->second > 0 ) R_wt += w;
    }
  const double denom = R_wt + L_wt;
  if ( denom <= 0.0 ) return std::numeric_limits<double>::quiet_NaN();
  return ( R_wt - L_wt ) / denom;
}


// ---------------------------------------------------------------------------
// Static transform helpers
// ---------------------------------------------------------------------------

void ged_t::apply_nb_filter( Eigen::MatrixXd & X , int sr , double f , double fwhm )
{
  for ( int s = 0; s < (int)X.cols(); s++ )
    X.col(s) = narrow_gaussian_t::filter( X.col(s) , sr , f , fwhm );
}

void ged_t::apply_env_filter( Eigen::MatrixXd & X , int sr ,
                               double lwr , double upr ,
                               double ripple , double tw )
{
  for ( int s = 0; s < (int)X.cols(); s++ )
    {
      std::vector<double> col( X.rows() );
      Eigen::VectorXd::Map( &col[0] , X.rows() ) = X.col(s);
      hilbert_t h( col , sr , lwr , upr , ripple , tw );
      const std::vector<double> * mag = h.magnitude();
      if ( mag == nullptr || (int)mag->size() != X.rows() )
        Helper::halt( "GED: Hilbert envelope failed for a channel" );
      for ( int r = 0; r < (int)X.rows(); r++ )
        X(r,s) = (*mag)[r];
    }
}

void ged_t::zscore_columns( Eigen::MatrixXd & X ,
                             std::vector<double> * means ,
                             std::vector<double> * sds )
{
  const int nc = (int)X.cols();
  const int nr = (int)X.rows();
  if ( means ) { means->resize( nc ); }
  if ( sds   ) { sds->resize( nc ); }

  for ( int s = 0; s < nc; s++ )
    {
      double mu  = X.col(s).mean();
      double var = ( X.col(s).array() - mu ).square().sum() / double( nr - 1 );
      double sd  = ( var > 0.0 ) ? std::sqrt(var) : 1.0;

      X.col(s) = ( X.col(s).array() - mu ) / sd;

      if ( means ) (*means)[s] = mu;
      if ( sds   ) (*sds)[s]   = sd;
    }
}

void ged_t::apply_zscore_stored( Eigen::MatrixXd & X ,
                                  const std::vector<double> & means ,
                                  const std::vector<double> & sds )
{
  if ( means.empty() ) return;
  const int nc = (int)X.cols();
  for ( int s = 0; s < nc && s < (int)means.size(); s++ )
    {
      double sd = ( s < (int)sds.size() && sds[s] > 0.0 ) ? sds[s] : 1.0;
      X.col(s) = ( X.col(s).array() - means[s] ) / sd;
    }
}


// ---------------------------------------------------------------------------
// ged_apply_input_transform() — dispatcher called before covariance build
// ---------------------------------------------------------------------------

void ged_apply_input_transform( Eigen::MatrixXd & X ,
                                ged_input_mode_t mode , int sr ,
                                double nb_f , double nb_fwhm ,
                                double env_lwr , double env_upr ,
                                double env_ripple , double env_tw ,
                                bool z_score ,
                                std::vector<double> * col_means ,
                                std::vector<double> * col_sds )
{
  if ( mode == ged_input_mode_t::NB )
    ged_t::apply_nb_filter( X , sr , nb_f , nb_fwhm );
  else if ( mode == ged_input_mode_t::ENV )
    ged_t::apply_env_filter( X , sr , env_lwr , env_upr , env_ripple , env_tw );

  if ( z_score )
    ged_t::zscore_columns( X , col_means , col_sds );
}


// ---------------------------------------------------------------------------
// ged_read_clocs() — simple tab/space delimited: label ap lat
// ---------------------------------------------------------------------------

bool ged_read_clocs( const std::string & fname ,
                     std::map<std::string,double> * chan_ap ,
                     std::map<std::string,int>    * chan_lat )
{
  if ( fname.empty() ) return false;
  std::ifstream IN( fname.c_str() );
  if ( ! IN.good() )
    {
      logger << "  ** WARNING: could not open clocs file " << fname << "\n";
      return false;
    }

  int cnt = 0;
  while ( ! IN.eof() )
    {
      std::string line;
      std::getline( IN , line );
      if ( line.empty() || line[0] == '#' ) continue;
      std::istringstream ss( line );
      std::string label;
      double ap;
      int lat;
      if ( ! ( ss >> label >> ap >> lat ) ) continue;
      std::string uc = Helper::toupper( label );
      if ( chan_ap  ) (*chan_ap)[uc]  = ap;
      if ( chan_lat ) (*chan_lat)[uc] = lat;
      ++cnt;
    }
  IN.close();
  return cnt > 0;
}


// ---------------------------------------------------------------------------
// ged_build_row_stages() — map each timepoint row to a stage label
// ---------------------------------------------------------------------------

std::vector<std::string> ged_build_row_stages( edf_t & edf ,
                                               const std::vector<uint64_t> & tp )
{
  const int nr = (int)tp.size();
  std::vector<std::string> row_stage( nr , "?" );

  if ( ! edf.timeline.epoched() ) return row_stage;

  // Build sorted epoch boundary vector
  std::vector<uint64_t> epoch_starts;
  std::vector<std::string> epoch_stage_labels;

  const std::vector<std::string> stage_codes = { "W","N1","N2","N3","R" };

  int ne = edf.timeline.first_epoch();
  while ( true )
    {
      int e = edf.timeline.next_epoch();
      if ( e == -1 ) break;
      interval_t intv = edf.timeline.epoch( e );
      epoch_starts.push_back( intv.start );

      std::string stg = "?";
      for ( const auto & sc : stage_codes )
        if ( edf.timeline.epoch_annotation( sc , e ) ) { stg = sc; break; }
      epoch_stage_labels.push_back( stg );
    }

  if ( epoch_starts.empty() ) return row_stage;

  // For each timepoint row, binary-search for its epoch
  for ( int r = 0; r < nr; r++ )
    {
      auto it = std::upper_bound( epoch_starts.begin() , epoch_starts.end() , tp[r] );
      if ( it != epoch_starts.begin() )
        {
          --it;
          int e_idx = (int)( it - epoch_starts.begin() );
          if ( e_idx < (int)epoch_stage_labels.size() )
            row_stage[r] = epoch_stage_labels[ e_idx ];
        }
    }

  return row_stage;
}


// ---------------------------------------------------------------------------
// ged_write_component_outputs() — write eigenvalues, maps for nc_out components
// stage_label == "" → whole-recording (no SS stratifier)
// ---------------------------------------------------------------------------

void ged_write_component_outputs( const ged_t & ged ,
                                  const std::vector<std::string> & ch_names ,
                                  const Eigen::MatrixXd & S ,
                                  int nc_out ,
                                  const std::string & stage_label ,
                                  const std::map<std::string,double> * chan_ap ,
                                  const std::map<std::string,int>    * chan_lat )
{
  if ( ! ged.solved ) return;

  const int nc  = (int)ged.L.size();
  if ( nc_out < 0 || nc_out > nc ) nc_out = nc;

  std::vector<int> si = ged.sorted_indices();

  // Scalar summaries (no COMP stratifier, optionally SS-stratified)
  if ( stage_label != "" )
    writer.level( stage_label , globals::stage_strat );

  writer.value( "N_S" , ged.n_S );
  writer.value( "N_R" , ged.n_R );

  // Per-component outputs
  for ( int rank = 0; rank < nc_out; rank++ )
    {
      const int e = si[ rank ];

      writer.level( rank + 1 , "COMP" );

      if ( stage_label != "" )
        writer.level( stage_label , globals::stage_strat );

      writer.value( "LAMBDA"      , ged.L[e] );
      writer.value( "LAMBDA_R"    , ged.eigenvalue_ratio( e ) );
      writer.value( "LAMBDA_RANK" , rank + 1 );

      // forward model
      int mxch = 0;
      Eigen::VectorXd spatialMap = const_cast<ged_t&>(ged).map( e , S , &mxch );

      writer.value( "FOC" , ged_t::focality( spatialMap ) );

      if ( chan_ap && !chan_ap->empty() )
        {
          double ap = ged_t::ap_index( spatialMap , ch_names , *chan_ap );
          if ( ! std::isnan(ap) ) writer.value( "AP" , ap );
        }
      if ( chan_lat && !chan_lat->empty() )
        {
          double lat = ged_t::lat_index( spatialMap , ch_names , *chan_lat );
          if ( ! std::isnan(lat) ) writer.value( "LAT" , lat );
        }

      // per-channel W and MAP
      for ( int s = 0; s < (int)ch_names.size() && s < nc; s++ )
        {
          writer.level( ch_names[s] , globals::signal_strat );
          writer.value( "W"   , ged.W( s , e ) );
          writer.value( "MAP" , spatialMap[s] );
          writer.unlevel( globals::signal_strat );
        }

      if ( stage_label != "" )
        writer.unlevel( globals::stage_strat );

      writer.unlevel( "COMP" );
    }

  if ( stage_label != "" )
    writer.unlevel( globals::stage_strat );
}


// ---------------------------------------------------------------------------
// ged_epoch_power() — per-epoch discriminant power: w^T C_epoch w
// ---------------------------------------------------------------------------

std::vector<std::pair<int,double>> ged_epoch_power( edf_t & edf ,
                                                    const signal_list_t & signals ,
                                                    const Eigen::VectorXd & w ,
                                                    ged_input_mode_t mode , int sr ,
                                                    double nb_f , double nb_fwhm ,
                                                    double env_lwr , double env_upr ,
                                                    double env_ripple , double env_tw ,
                                                    const std::vector<double> & col_means ,
                                                    const std::vector<double> & col_sds )
{
  std::vector<std::pair<int,double>> result;

  if ( ! edf.timeline.epoched() ) return result;

  const int ns = (int)w.size();

  int ne = edf.timeline.first_epoch();
  while ( true )
    {
      int e = edf.timeline.next_epoch();
      if ( e == -1 ) break;

      int display_e = edf.timeline.display_epoch( e );
      interval_t intv = edf.timeline.epoch( e );

      eigen_matslice_t eslice( edf , signals , intv );
      Eigen::MatrixXd Xe = eslice.data_ref();

      if ( Xe.rows() < 2 || (int)Xe.cols() != ns ) continue;

      // apply same transform as whole-trace (using stored params)
      if ( mode == ged_input_mode_t::NB )
        ged_t::apply_nb_filter( Xe , sr , nb_f , nb_fwhm );
      else if ( mode == ged_input_mode_t::ENV )
        ged_t::apply_env_filter( Xe , sr , env_lwr , env_upr , env_ripple , env_tw );

      if ( ! col_means.empty() )
        ged_t::apply_zscore_stored( Xe , col_means , col_sds );

      if ( Xe.rows() < ns ) continue;  // rank-deficient skip

      Eigen::MatrixXd Ce = eigen_ops::covariance( Xe );
      double pwr = w.transpose() * Ce * w;

      result.push_back( { display_e , pwr } );
    }

  return result;
}


// ---------------------------------------------------------------------------
// Main wrapper — dispatches to appropriate run mode
// ---------------------------------------------------------------------------

void ged_wrapper( edf_t & edf , param_t & param )
{
  const bool NO_ANNOTS = true;
  signal_list_t signals = edf.header.signal_list( param.requires("sig") , NO_ANNOTS );
  const int ns = (int)signals.size();
  if ( ns < 2 )
    {
      logger << "  ** GED: need at least 2 channels; skipping\n";
      return;
    }

  std::vector<double> Fs = edf.header.sampling_freq( signals );

  // Validate uniform sampling rate
  for ( int s = 1; s < ns; s++ )
    if ( std::abs( Fs[s] - Fs[0] ) > 0.5 )
      Helper::halt( "GED: all channels must have the same sampling rate" );

  const int sr = (int)Fs[0];

  // Get whole-trace data
  eigen_matslice_t mslice( edf , signals , edf.timeline.wholetrace() );
  Eigen::MatrixXd X = mslice.data_ref();

  // Load group solution if requested
  if ( param.has("load") )
    {
      const std::string sol_file = param.value("load");
      if ( ! ged_t::has_group_solution || ged_t::group_solution_fname != sol_file )
        ged_group_t::load_solution( sol_file );

      const std::vector<uint64_t> * tp = mslice.ptimepoints();
      ged_apply_group( edf , param , X , tp , sr );
      return;
    }

  // Determine run mode from parameters
  const bool has_a1    = param.has("a1");
  const std::string input_str = param.has("input") ? param.value("input") : "raw";
  const ged_input_mode_t mode = ged_input_mode_from_str( input_str );

  if ( ! has_a1 && mode == ged_input_mode_t::RAW )
    Helper::halt( "GED: specify a1= (annotation-locked) or input=nb/env (filter-based)" );

  if ( has_a1 )
    {
      const std::vector<uint64_t> * tp = mslice.ptimepoints();
      ged_runmode2( edf , param , X , tp , sr );
    }
  else
    {
      ged_runmode1( edf , param , X , sr );
    }
}


// ---------------------------------------------------------------------------
// ged_runmode1 — narrowband/envelope vs. broadband (no annotation time-lock)
// ---------------------------------------------------------------------------

void ged_runmode1( edf_t & edf , param_t & param , Eigen::MatrixXd & Rd , int sr )
{
  // Input transform parameters
  const std::string input_str = param.has("input") ? param.value("input") : "nb";
  const ged_input_mode_t mode = ged_input_mode_from_str( input_str );

  double nb_f    = param.has("nb-f")    ? param.requires_dbl("nb-f")    :
                   param.has("f1")      ? param.requires_dbl("f1")      : 0.0;
  double nb_fwhm = param.has("nb-fwhm") ? param.requires_dbl("nb-fwhm") :
                   param.has("fwhm1")   ? param.requires_dbl("fwhm1")   : 0.0;
  double env_lwr = param.has("env-lwr") ? param.requires_dbl("env-lwr") : 0.0;
  double env_upr = param.has("env-upr") ? param.requires_dbl("env-upr") : 0.0;
  double env_ripple = param.has("env-ripple") ? param.requires_dbl("env-ripple") : 0.02;
  double env_tw     = param.has("env-tw")     ? param.requires_dbl("env-tw")     : 1.0;
  const bool z_score = param.has("z");
  const double reg   = param.has("reg") ? param.requires_dbl("reg") : 0.0;
  int nc_out = param.has("nc") ? param.requires_int("nc") : -1;
  const std::string new_ts = param.has("ts") ? param.value("ts") : "";
  const bool do_stages = param.has("stages");

  // Load channel locations if provided
  std::map<std::string,double> chan_ap;
  std::map<std::string,int>    chan_lat;
  bool has_clocs = false;
  if ( param.has("clocs") )
    has_clocs = ged_read_clocs( param.value("clocs") , &chan_ap , &chan_lat );

  // Channel names
  const bool NO_ANNOTS = true;
  signal_list_t signals = edf.header.signal_list( param.requires("sig") , NO_ANNOTS );
  const int ns = (int)signals.size();
  std::vector<std::string> ch_names( ns );
  for ( int s = 0; s < ns; s++ ) ch_names[s] = signals.label(s);

  // Validate mode
  if ( mode == ged_input_mode_t::NB && ( nb_f <= 0.0 || nb_fwhm <= 0.0 ) )
    Helper::halt( "GED mode 1 NB: specify nb-f and nb-fwhm (or f1/fwhm1)" );
  if ( mode == ged_input_mode_t::ENV && ( env_lwr <= 0.0 || env_upr <= 0.0 ) )
    Helper::halt( "GED mode 1 ENV: specify env-lwr and env-upr" );

  // S matrix: transform the data
  Eigen::MatrixXd Sd = Rd;
  std::vector<double> col_means , col_sds;

  logger << "  GED mode 1: building S from " << input_str << " transform\n";
  ged_apply_input_transform( Sd , mode , sr , nb_f , nb_fwhm ,
                             env_lwr , env_upr , env_ripple , env_tw ,
                             z_score , &col_means , &col_sds );

  Eigen::MatrixXd S = eigen_ops::covariance( Sd );

  // R matrix: broadband (raw), optionally z-scored same as S
  Eigen::MatrixXd Rd_r = Rd;
  if ( z_score && ! col_means.empty() )
    ged_t::apply_zscore_stored( Rd_r , col_means , col_sds );
  Eigen::MatrixXd R = eigen_ops::covariance( Rd_r );

  // GED
  ged_t ged;
  ged.covar( S , R );
  ged.n_S = (int)Sd.rows();
  ged.n_R = (int)Rd_r.rows();
  if ( reg > 0.0 ) ged.regularize_R( reg );
  ged.calc();

  if ( ! ged.solved ) return;

  if ( ns < nc_out || nc_out < 0 ) nc_out = ns;

  // Primary output
  ged_write_component_outputs( ged , ch_names , S , nc_out , "" ,
                               has_clocs ? &chan_ap  : nullptr ,
                               has_clocs ? &chan_lat : nullptr );

  // Time series for top component
  if ( new_ts != "" )
    {
      std::vector<int> si = ged.sorted_indices();
      int top_e = si[0];
      int mxch  = 0;
      ged.map( top_e , S , &mxch );
      Eigen::VectorXd ts = ged.time_series( top_e , Sd , mxch );
      std::vector<double> tsvec( ts.size() );
      Eigen::VectorXd::Map( &tsvec[0] , ts.size() ) = ts;
      // note: if new_ts already exists as a channel it will be duplicated; avoid running GED ts= twice
      logger << "  adding channel " << new_ts << "\n";
      edf.add_signal( new_ts , sr , tsvec );
    }

  // Epoch power
  if ( edf.timeline.epoched() )
    {
      std::vector<int> si = ged.sorted_indices();
      Eigen::VectorXd w = ged.W.col( si[0] );
      auto ep = ged_epoch_power( edf , signals , w , mode , sr ,
                                 nb_f , nb_fwhm , env_lwr , env_upr , env_ripple , env_tw ,
                                 col_means , col_sds );
      if ( ! ep.empty() )
        {
          double sum = 0.0 , sum2 = 0.0;
          for ( const auto & p : ep ) { sum += p.second; sum2 += p.second * p.second; }
          double mu  = sum / ep.size();
          double var = sum2 / ep.size() - mu*mu;
          double sd  = ( var > 0.0 ) ? std::sqrt(var) : 0.0;
          double cv  = ( mu  > 0.0 ) ? sd / mu        : 0.0;
          writer.value( "EPOW_MEAN" , mu );
          writer.value( "EPOW_SD"   , sd );
          writer.value( "EPOW_CV"   , cv );
          for ( const auto & p : ep )
            {
              writer.level( p.first , globals::epoch_strat );
              writer.value( "EPOW" , p.second );
              writer.unlevel( globals::epoch_strat );
            }
        }
    }

  // Stage-stratified GED
  if ( do_stages && edf.timeline.epoched() )
    {
      // Build timepoints vector from Rd row count using wholetrace mslice
      eigen_matslice_t mslice1( edf , signals , edf.timeline.wholetrace() );
      const std::vector<uint64_t> & tp = *mslice1.ptimepoints();
      std::vector<std::string> row_stages = ged_build_row_stages( edf , tp );
      const std::vector<std::string> stage_codes = { "W","N1","N2","N3","R" };

      for ( const auto & stg : stage_codes )
        {
          // Select rows belonging to this stage
          std::vector<int> rows;
          for ( int r = 0; r < (int)row_stages.size(); r++ )
            if ( row_stages[r] == stg ) rows.push_back(r);

          if ( (int)rows.size() < ns + 1 ) continue;

          Eigen::MatrixXd Xs( rows.size() , ns );
          for ( int i = 0; i < (int)rows.size(); i++ )
            Xs.row(i) = Rd.row( rows[i] );

          // S: apply spectral transform (no z-score here — use stored params)
          Eigen::MatrixXd Xs_t = Xs;
          if ( mode == ged_input_mode_t::NB )
            ged_t::apply_nb_filter( Xs_t , sr , nb_f , nb_fwhm );
          else if ( mode == ged_input_mode_t::ENV )
            ged_t::apply_env_filter( Xs_t , sr , env_lwr , env_upr , env_ripple , env_tw );
          if ( z_score && ! col_means.empty() )
            ged_t::apply_zscore_stored( Xs_t , col_means , col_sds );
          Eigen::MatrixXd Ss = eigen_ops::covariance( Xs_t );

          // R: broadband, z-score with stored params
          Eigen::MatrixXd Xs_r = Xs;
          if ( z_score && ! col_means.empty() )
            ged_t::apply_zscore_stored( Xs_r , col_means , col_sds );
          Eigen::MatrixXd Rs = eigen_ops::covariance( Xs_r );

          ged_t geds;
          geds.covar( Ss , Rs );
          geds.n_S = (int)rows.size();
          geds.n_R = (int)rows.size();
          if ( reg > 0.0 ) geds.regularize_R( reg );
          geds.calc();

          if ( geds.solved )
            ged_write_component_outputs( geds , ch_names , Ss , nc_out , stg ,
                                         has_clocs ? &chan_ap  : nullptr ,
                                         has_clocs ? &chan_lat : nullptr );
        }
    }

  // save-cov
  if ( param.has("save-cov") )
    {
      ged_group_t::save_covar( param.value("save-cov") ,
                               edf.id , ch_names , S , R ,
                               ged.n_S , ged.n_R ,
                               mode , nb_f , nb_fwhm , env_lwr , env_upr ,
                               z_score , col_means , col_sds );
    }
}


// ---------------------------------------------------------------------------
// ged_runmode2 — annotation-locked S vs. reference R
// ---------------------------------------------------------------------------

void ged_runmode2( edf_t & edf , param_t & param ,
                   Eigen::MatrixXd & Rd , const std::vector<uint64_t> * tp , int sr )
{
  // Parameters
  const std::string a1 = param.requires("a1");
  const double w1      = param.has("w1") ? param.requires_dbl("w1") : 0.0;
  const bool x1        = param.has("x1");
  const bool refall    = ! param.has("a2");
  const std::string a2 = refall ? "" : param.value("a2");
  const double w2      = param.has("w2") ? param.requires_dbl("w2") : 0.0;
  const bool x2        = param.has("x2");

  const std::string input_str = param.has("input") ? param.value("input") : "raw";
  const ged_input_mode_t mode = ged_input_mode_from_str( input_str );
  double nb_f    = param.has("nb-f")    ? param.requires_dbl("nb-f")    : 0.0;
  double nb_fwhm = param.has("nb-fwhm") ? param.requires_dbl("nb-fwhm") : 0.0;
  double env_lwr = param.has("env-lwr") ? param.requires_dbl("env-lwr") : 0.0;
  double env_upr = param.has("env-upr") ? param.requires_dbl("env-upr") : 0.0;
  double env_ripple = param.has("env-ripple") ? param.requires_dbl("env-ripple") : 0.02;
  double env_tw     = param.has("env-tw")     ? param.requires_dbl("env-tw")     : 1.0;
  const bool z_score   = param.has("z");
  const double reg     = param.has("reg") ? param.requires_dbl("reg") : 0.0;
  int nc_out = param.has("nc") ? param.requires_int("nc") : -1;
  const std::string new_ts  = param.has("ts")  ? param.value("ts")  : "";
  const bool do_stages = param.has("stages");
  const bool do_epoch  = param.has("win") || edf.timeline.epoched();

  // Channel names
  const bool NO_ANNOTS = true;
  signal_list_t signals = edf.header.signal_list( param.requires("sig") , NO_ANNOTS );
  const int ns = (int)signals.size();
  std::vector<std::string> ch_names( ns );
  for ( int s = 0; s < ns; s++ ) ch_names[s] = signals.label(s);

  // Clocs
  std::map<std::string,double> chan_ap;
  std::map<std::string,int>    chan_lat;
  bool has_clocs = false;
  if ( param.has("clocs") )
    has_clocs = ged_read_clocs( param.value("clocs") , &chan_ap , &chan_lat );

  // Annotations
  annot_t * annot1 = edf.annotations->find( a1 );
  if ( annot1 == nullptr )
    Helper::halt( "GED: could not find annotation " + a1 );
  annot_t * annot2 = nullptr;
  if ( ! refall )
    {
      annot2 = edf.annotations->find( a2 );
      if ( annot2 == nullptr )
        Helper::halt( "GED: could not find annotation " + a2 );
    }

  // Apply input transform to whole-trace data
  std::vector<double> col_means , col_sds;
  ged_apply_input_transform( Rd , mode , sr , nb_f , nb_fwhm ,
                             env_lwr , env_upr , env_ripple , env_tw ,
                             z_score , &col_means , &col_sds );

  // Build S and R covariance matrices
  Eigen::MatrixXd Sd = eigen_ops::subset_rows( Rd , tp , annot1 , w1 , x1 );
  logger << "  GED: S matrix from " << Sd.rows() << " of " << Rd.rows() << " rows\n";

  if ( Sd.rows() < ns + 1 )
    {
      logger << "  ** GED: too few timepoints for S matrix; skipping\n";
      return;
    }

  Eigen::MatrixXd S = eigen_ops::covariance( Sd );

  Eigen::MatrixXd R;
  int n_R_val = 0;
  if ( refall )
    {
      R = eigen_ops::covariance( Rd );
      n_R_val = (int)Rd.rows();
    }
  else
    {
      Eigen::MatrixXd Rd2 = eigen_ops::subset_rows( Rd , tp , annot2 , w2 , x2 );
      logger << "  GED: R matrix from " << Rd2.rows() << " of " << Rd.rows() << " rows\n";
      if ( Rd2.rows() < ns + 1 )
        {
          logger << "  ** GED: too few timepoints for R matrix; skipping\n";
          return;
        }
      R = eigen_ops::covariance( Rd2 );
      n_R_val = (int)Rd2.rows();
    }

  // GED
  ged_t ged;
  ged.covar( S , R );
  ged.n_S = (int)Sd.rows();
  ged.n_R = n_R_val;
  if ( reg > 0.0 ) ged.regularize_R( reg );
  ged.calc();

  if ( ! ged.solved ) return;

  if ( nc_out < 0 || nc_out > ns ) nc_out = ns;

  // Primary output
  ged_write_component_outputs( ged , ch_names , S , nc_out , "" ,
                               has_clocs ? &chan_ap  : nullptr ,
                               has_clocs ? &chan_lat : nullptr );

  // Time series
  if ( new_ts != "" )
    {
      int ts_comp = param.has("ts-comp") ? param.requires_int("ts-comp") - 1 : 0;
      std::vector<int> si = ged.sorted_indices();
      if ( ts_comp < 0 || ts_comp >= (int)si.size() ) ts_comp = 0;
      int top_e = si[ ts_comp ];
      int mxch  = 0;
      ged.map( top_e , S , &mxch );
      Eigen::VectorXd ts = ged.time_series( top_e , Rd , mxch );
      std::vector<double> tsvec( ts.size() );
      Eigen::VectorXd::Map( &tsvec[0] , ts.size() ) = ts;
      if ( edf.header.has_signal( new_ts ) )
        edf.drop_signal( edf.header.signal_list(new_ts,NO_ANNOTS)(0) );
      logger << "  adding channel " << new_ts << "\n";
      edf.add_signal( new_ts , sr , tsvec );
    }

  // Epoch power
  if ( do_epoch && edf.timeline.epoched() )
    {
      std::vector<int> si = ged.sorted_indices();
      Eigen::VectorXd w = ged.W.col( si[0] );
      auto ep = ged_epoch_power( edf , signals , w , mode , sr ,
                                 nb_f , nb_fwhm , env_lwr , env_upr , env_ripple , env_tw ,
                                 col_means , col_sds );
      if ( ! ep.empty() )
        {
          double sum = 0.0 , sum2 = 0.0;
          for ( const auto & p : ep ) { sum += p.second; sum2 += p.second * p.second; }
          double mu  = sum / ep.size();
          double var = sum2 / ep.size() - mu*mu;
          double sd  = ( var > 0.0 ) ? std::sqrt(var) : 0.0;
          double cv  = ( mu  > 0.0 ) ? sd / mu        : 0.0;
          writer.value( "EPOW_MEAN" , mu );
          writer.value( "EPOW_SD"   , sd );
          writer.value( "EPOW_CV"   , cv );

          // Get stage labels per epoch for optional stage-epoch cross
          std::vector<std::string> ep_stages;
          if ( do_stages )
            {
              int ne = edf.timeline.first_epoch();
              while (true)
                {
                  int e = edf.timeline.next_epoch();
                  if ( e == -1 ) break;
                  std::string stg = "?";
                  const std::vector<std::string> sc2 = {"W","N1","N2","N3","R"};
                  for ( const auto & c : sc2 )
                    if ( edf.timeline.epoch_annotation(c,e) ) { stg=c; break; }
                  ep_stages.push_back( stg );
                }
            }

          for ( int i = 0; i < (int)ep.size(); i++ )
            {
              writer.level( ep[i].first , globals::epoch_strat );
              writer.value( "EPOW" , ep[i].second );
              if ( do_stages && i < (int)ep_stages.size() )
                writer.value( "EPOW_SS" , ep_stages[i] );
              writer.unlevel( globals::epoch_strat );
            }
        }
    }

  // Stage-stratified GED
  if ( do_stages )
    {
      std::vector<std::string> row_stages = ged_build_row_stages( edf , *tp );
      const std::vector<std::string> stage_codes = { "W","N1","N2","N3","R" };

      for ( const auto & stg : stage_codes )
        {
          std::vector<int> rows;
          for ( int r = 0; r < (int)row_stages.size(); r++ )
            if ( row_stages[r] == stg ) rows.push_back(r);

          if ( (int)rows.size() < ns + 1 ) continue;

          // Build stage-subset of the TRANSFORMED data matrix
          Eigen::MatrixXd Xstg( rows.size() , ns );
          std::vector<uint64_t> tp_stg( rows.size() );
          for ( int i = 0; i < (int)rows.size(); i++ )
            {
              Xstg.row(i) = Rd.row( rows[i] );
              tp_stg[i] = (*tp)[ rows[i] ];
            }

          // S for this stage
          Eigen::MatrixXd Sd_s = eigen_ops::subset_rows( Xstg , &tp_stg , annot1 , w1 , x1 );
          if ( Sd_s.rows() < ns + 1 ) continue;
          Eigen::MatrixXd Ss = eigen_ops::covariance( Sd_s );

          // R for this stage
          Eigen::MatrixXd Rs;
          int n_R_s = 0;
          if ( refall )
            {
              Rs = eigen_ops::covariance( Xstg );
              n_R_s = (int)Xstg.rows();
            }
          else
            {
              Eigen::MatrixXd Rd2_s = eigen_ops::subset_rows( Xstg , &tp_stg , annot2 , w2 , x2 );
              if ( Rd2_s.rows() < ns + 1 ) continue;
              Rs = eigen_ops::covariance( Rd2_s );
              n_R_s = (int)Rd2_s.rows();
            }

          ged_t geds;
          geds.covar( Ss , Rs );
          geds.n_S = (int)Sd_s.rows();
          geds.n_R = n_R_s;
          if ( reg > 0.0 ) geds.regularize_R( reg );
          geds.calc();

          if ( geds.solved )
            ged_write_component_outputs( geds , ch_names , Ss , nc_out , stg ,
                                         has_clocs ? &chan_ap  : nullptr ,
                                         has_clocs ? &chan_lat : nullptr );
        }
    }

  // Save per-individual covariance to accumulation file
  if ( param.has("save-cov") )
    {
      ged_group_t::save_covar( param.value("save-cov") ,
                               edf.id , ch_names , S , R ,
                               ged.n_S , ged.n_R ,
                               mode , nb_f , nb_fwhm , env_lwr , env_upr ,
                               z_score , col_means , col_sds );
    }
}


// ---------------------------------------------------------------------------
// ged_apply_group() — project individual through pre-loaded group solution
// ---------------------------------------------------------------------------

void ged_apply_group( edf_t & edf , param_t & param ,
                      Eigen::MatrixXd & X ,
                      const std::vector<uint64_t> * tp , int sr )
{
  if ( ! ged_t::has_group_solution )
    Helper::halt( "GED: no group solution loaded" );

  const bool NO_ANNOTS = true;
  signal_list_t signals = edf.header.signal_list( param.requires("sig") , NO_ANNOTS );
  const int ns_indiv = (int)signals.size();
  const int nc_group = (int)ged_t::group_channels.size();

  // Build permutation: individual channels -> group channel order
  std::vector<int> perm( nc_group , -1 );
  for ( int g = 0; g < nc_group; g++ )
    {
      for ( int s = 0; s < ns_indiv; s++ )
        {
          if ( Helper::toupper( signals.label(s) ) ==
               Helper::toupper( ged_t::group_channels[g] ) )
            { perm[g] = s; break; }
        }
      if ( perm[g] == -1 )
        Helper::halt( "GED load: individual missing group channel: " + ged_t::group_channels[g] );
    }

  // Reorder X columns to match group channel order
  Eigen::MatrixXd Xg( X.rows() , nc_group );
  for ( int g = 0; g < nc_group; g++ )
    Xg.col(g) = X.col( perm[g] );

  // Apply same transform as group solution (using stored means/sds)
  const ged_input_mode_t mode = ged_t::group_input_mode;
  double nb_f    = ged_t::group_nb_f;
  double nb_fwhm = ged_t::group_nb_fwhm;
  double env_lwr = ged_t::group_env_lwr;
  double env_upr = ged_t::group_env_upr;
  const double env_ripple = 0.02;
  const double env_tw     = 1.0;

  if ( mode == ged_input_mode_t::NB )
    ged_t::apply_nb_filter( Xg , sr , nb_f , nb_fwhm );
  else if ( mode == ged_input_mode_t::ENV )
    ged_t::apply_env_filter( Xg , sr , env_lwr , env_upr , env_ripple , env_tw );

  if ( ged_t::group_z_scored )
    ged_t::apply_zscore_stored( Xg , ged_t::group_ch_means , ged_t::group_ch_sds );

  // Individual covariance
  Eigen::MatrixXd Ci = eigen_ops::covariance( Xg );

  // Per-component scores
  const int nc_sol = (int)ged_t::group_L.size();
  int nc_out = param.has("nc") ? param.requires_int("nc") : nc_sol;
  if ( nc_out < 0 || nc_out > nc_sol ) nc_out = nc_sol;

  for ( int k = 0; k < nc_out; k++ )
    {
      Eigen::VectorXd w = ged_t::group_W.col(k);
      double pwr   = w.transpose() * Ci * w;
      double lam_g = ged_t::group_L[k];
      double pwr_r = ( lam_g > 0.0 ) ? pwr / lam_g : 0.0;

      writer.level( k + 1 , "COMP" );
      writer.value( "POWER"   , pwr );
      writer.value( "POWER_R" , pwr_r );
      writer.unlevel( "COMP" );
    }

  // Time series for top (or specified) component
  const std::string new_ts = param.has("ts") ? param.value("ts") : "";
  if ( new_ts != "" )
    {
      int ts_comp = param.has("ts-comp") ? param.requires_int("ts-comp") - 1 : 0;
      if ( ts_comp < 0 || ts_comp >= nc_sol ) ts_comp = 0;
      Eigen::VectorXd w = ged_t::group_W.col( ts_comp );
      Eigen::VectorXd ts_vec = Xg * w;
      // sign-fix against dominant channel
      int mxch = 0;
      (Ci * w).cwiseAbs().maxCoeff( &mxch );
      if ( (Ci * w)[mxch] < 0 ) ts_vec *= -1;
      std::vector<double> tsvec( ts_vec.size() );
      Eigen::VectorXd::Map( &tsvec[0] , ts_vec.size() ) = ts_vec;
      if ( edf.header.has_signal( new_ts ) )
        edf.drop_signal( edf.header.signal_list(new_ts,NO_ANNOTS)(0) );
      edf.add_signal( new_ts , sr , tsvec );
      logger << "  GED load: added channel " << new_ts << "\n";
    }

  // Epoch power for top component
  if ( edf.timeline.epoched() )
    {
      Eigen::VectorXd w = ged_t::group_W.col(0);
      // Build signal_list for group channels
      signal_list_t grp_signals = edf.header.signal_list(
          param.requires("sig") , NO_ANNOTS );

      auto ep = ged_epoch_power( edf , grp_signals , w , mode , sr ,
                                 nb_f , nb_fwhm , env_lwr , env_upr ,
                                 env_ripple , env_tw ,
                                 ged_t::group_ch_means , ged_t::group_ch_sds );
      if ( ! ep.empty() )
        {
          double sum=0, sum2=0;
          for ( const auto & p : ep ) { sum+=p.second; sum2+=p.second*p.second; }
          double mu  = sum/ep.size();
          double var = sum2/ep.size() - mu*mu;
          double sd  = (var>0)?std::sqrt(var):0;
          double cv  = (mu>0)?sd/mu:0;
          writer.value( "EPOW_MEAN" , mu );
          writer.value( "EPOW_SD"   , sd );
          writer.value( "EPOW_CV"   , cv );
          for ( const auto & p : ep )
            {
              writer.level( p.first , globals::epoch_strat );
              writer.value( "EPOW" , p.second );
              writer.unlevel( globals::epoch_strat );
            }
        }
    }

  // Stage-stratified scores
  if ( param.has("stages") && edf.timeline.epoched() )
    {
      // rebuild with group-ordered tp pointer from original mslice
      // (X was the raw whole-trace data, tp is from original mslice)
      std::vector<std::string> row_stages = ged_build_row_stages( edf , *tp );
      const std::vector<std::string> stage_codes = { "W","N1","N2","N3","R" };

      for ( const auto & stg : stage_codes )
        {
          std::vector<int> rows;
          for ( int r=0; r<(int)row_stages.size(); r++ )
            if ( row_stages[r] == stg ) rows.push_back(r);
          if ( rows.empty() ) continue;

          Eigen::MatrixXd Xs( rows.size() , nc_group );
          for ( int i=0; i<(int)rows.size(); i++ )
            Xs.row(i) = Xg.row( rows[i] );

          Eigen::MatrixXd Cs = eigen_ops::covariance( Xs );

          writer.level( stg , globals::stage_strat );
          for ( int k=0; k<nc_out; k++ )
            {
              Eigen::VectorXd w = ged_t::group_W.col(k);
              double pwr   = w.transpose() * Cs * w;
              double lam_g = ged_t::group_L[k];
              double pwr_r = (lam_g>0)? pwr/lam_g : 0.0;
              writer.level( k+1 , "COMP" );
              writer.value( "POWER"   , pwr );
              writer.value( "POWER_R" , pwr_r );
              writer.unlevel( "COMP" );
            }
          writer.unlevel( globals::stage_strat );
        }
    }
}


// ===========================================================================
// ged_group_t — text accumulation and solution I/O
//
// Accumulation file: one self-contained 5-line record per individual.
// Files can be concatenated with 'cat' for parallel workflows.
//
//   GED_RECORD id=<id> nc=<nc> n_S=<n> n_R=<n> mode=<RAW|NB|ENV> \
//              nb_f=<f> nb_fwhm=<fwhm> env_lwr=<l> env_upr=<u> z=<0|1> \
//              ch=<ch1,ch2,...,chN>
//   S <nc*nc doubles row-major>
//   R <nc*nc doubles row-major>
//   MEANS <nc doubles>
//   SDS <nc doubles>
//
// Solution file: single 5-line text record.
//
//   GED_SOLUTION nc=<nc> nc_sol=<nc_sol> mode=... ch=<ch1,...,chN>
//   L <nc_sol doubles>
//   MEANS <nc doubles>
//   SDS <nc doubles>
//   W <nc*nc_sol doubles column-major>
//
// ===========================================================================

// Parse a key=value token list into a map
static std::map<std::string,std::string> ged_parse_kv( const std::string & line )
{
  std::map<std::string,std::string> kv;
  std::istringstream ss( line );
  std::string tok;
  while ( ss >> tok )
    {
      const size_t eq = tok.find('=');
      if ( eq != std::string::npos )
        kv[ tok.substr(0,eq) ] = tok.substr(eq+1);
    }
  return kv;
}


// ---------------------------------------------------------------------------
// ged_group_t::save_covar() — append per-individual S/R to accumulation file
// ---------------------------------------------------------------------------

void ged_group_t::save_covar( const std::string & fname ,
                              const std::string & id ,
                              const std::vector<std::string> & ch_names ,
                              const Eigen::MatrixXd & S ,
                              const Eigen::MatrixXd & R ,
                              int n_S , int n_R ,
                              ged_input_mode_t mode ,
                              double nb_f , double nb_fwhm ,
                              double env_lwr , double env_upr ,
                              bool z_scored ,
                              const std::vector<double> & col_means ,
                              const std::vector<double> & col_sds )
{
  const int nc = (int)ch_names.size();

  // If file exists, validate that first record's metadata matches
  std::ifstream CHKF( fname.c_str() );
  if ( CHKF.good() )
    {
      std::string line;
      while ( std::getline( CHKF , line ) )
        if ( line.substr(0,10) == "GED_RECORD" ) break;
      CHKF.close();

      if ( line.substr(0,10) == "GED_RECORD" )
        {
          auto kv = ged_parse_kv( line );
          if ( std::stoi( kv["nc"] ) != nc )
            Helper::halt( "GED save-cov: channel count mismatch in existing file" );
          if ( kv["mode"] != ged_input_mode_str(mode) )
            Helper::halt( "GED save-cov: input mode mismatch in existing file" );
          // check channel names
          std::string ch_str;
          for ( int s=0; s<nc; s++ ) ch_str += (s?",":"") + ch_names[s];
          if ( kv["ch"] != ch_str )
            Helper::halt( "GED save-cov: channel names mismatch in existing file" );
        }
    }

  // Append record (creates file if new)
  std::ofstream OUT( fname.c_str() , std::ios::app );
  if ( ! OUT.good() )
    Helper::halt( "GED save-cov: cannot write to " + fname );

  OUT << std::setprecision(15);

  // Header line
  std::string ch_str;
  for ( int s=0; s<nc; s++ ) ch_str += (s?",":"") + ch_names[s];
  OUT << "GED_RECORD"
      << " id="       << id
      << " nc="       << nc
      << " n_S="      << n_S
      << " n_R="      << n_R
      << " mode="     << ged_input_mode_str(mode)
      << " nb_f="     << nb_f
      << " nb_fwhm="  << nb_fwhm
      << " env_lwr="  << env_lwr
      << " env_upr="  << env_upr
      << " z="        << (z_scored ? 1 : 0)
      << " ch="       << ch_str
      << "\n";

  // S matrix row-major
  OUT << "S";
  for ( int r=0; r<nc; r++ )
    for ( int c=0; c<nc; c++ )
      OUT << " " << S(r,c);
  OUT << "\n";

  // R matrix row-major
  OUT << "R";
  for ( int r=0; r<nc; r++ )
    for ( int c=0; c<nc; c++ )
      OUT << " " << R(r,c);
  OUT << "\n";

  // Means (zeros if not z-scored)
  OUT << "MEANS";
  for ( int s=0; s<nc; s++ )
    OUT << " " << ( z_scored && s<(int)col_means.size() ? col_means[s] : 0.0 );
  OUT << "\n";

  // SDs (ones if not z-scored)
  OUT << "SDS";
  for ( int s=0; s<nc; s++ )
    OUT << " " << ( z_scored && s<(int)col_sds.size() ? col_sds[s] : 1.0 );
  OUT << "\n";

  OUT.close();
  logger << "  GED: appended covariance for " << id << " to " << fname << "\n";
}


// ---------------------------------------------------------------------------
// ged_group_t::run_group() — read accumulated text file(s), compute group GED, save solution
// ---------------------------------------------------------------------------

void ged_group_t::run_group( param_t & param )
{
  const std::string dat_file = param.requires("dat");
  const std::string sol_file = param.requires("sol");
  const bool trace_norm = param.has("trace-norm");
  const double reg      = param.has("reg") ? param.requires_dbl("reg") : 0.0;
  const int min_n       = param.has("min-n") ? param.requires_int("min-n") : 5;
  int nc_out = param.has("nc") ? param.requires_int("nc") : -1;

  std::ifstream IN( dat_file.c_str() );
  if ( ! IN.good() )
    Helper::halt( "ged --ged-group: cannot open " + dat_file );

  // State initialised from first GED_RECORD encountered
  int nc = -1;
  std::vector<std::string> ch_names;
  ged_input_mode_t mode = ged_input_mode_t::RAW;
  double nb_f = 0, nb_fwhm = 0, env_lwr = 0, env_upr = 0;
  bool z_scored = false;

  Eigen::MatrixXd S_sum;
  Eigen::MatrixXd R_sum;
  std::vector<double> means_sum;
  std::vector<double> sds_sum;
  int n_indiv = 0;

  std::string line;
  while ( std::getline( IN , line ) )
    {
      if ( line.substr(0,10) != "GED_RECORD" ) continue;

      auto kv = ged_parse_kv( line );

      const int nc_i   = std::stoi( kv["nc"] );
      const int n_S_i  = std::stoi( kv["n_S"] );
      const int n_R_i  = std::stoi( kv["n_R"] );
      const std::string id_i = kv["id"];

      // On first record: initialise group state
      if ( nc < 0 )
        {
          nc = nc_i;
          mode     = ged_input_mode_from_str( kv["mode"] );
          nb_f     = std::stod( kv["nb_f"] );
          nb_fwhm  = std::stod( kv["nb_fwhm"] );
          env_lwr  = std::stod( kv["env_lwr"] );
          env_upr  = std::stod( kv["env_upr"] );
          z_scored = ( kv["z"] == "1" );
          // parse comma-separated channel list
          std::istringstream chs( kv["ch"] );
          std::string ch;
          while ( std::getline(chs, ch, ',') ) ch_names.push_back(ch);
          S_sum    = Eigen::MatrixXd::Zero( nc , nc );
          R_sum    = Eigen::MatrixXd::Zero( nc , nc );
          means_sum.assign( nc , 0.0 );
          sds_sum.assign(   nc , 0.0 );
          logger << "  GED group: " << nc << " channels, mode=" << ged_input_mode_str(mode) << "\n";
        }
      else
        {
          // Validate consistency
          if ( nc_i != nc )
            Helper::halt( "GED group: channel count mismatch for " + id_i );
          if ( ged_input_mode_from_str(kv["mode"]) != mode )
            Helper::halt( "GED group: input mode mismatch for " + id_i );
        }

      // Read the next 4 data lines: S, R, MEANS, SDS
      Eigen::MatrixXd Si( nc , nc ) , Ri( nc , nc );
      std::vector<double> means_i( nc ) , sds_i( nc );

      // S
      std::string dline;
      std::getline( IN , dline );
      { std::istringstream ds(dline); std::string tag; ds >> tag;
        for (int r=0;r<nc;r++) for (int c=0;c<nc;c++) ds >> Si(r,c); }
      // R
      std::getline( IN , dline );
      { std::istringstream ds(dline); std::string tag; ds >> tag;
        for (int r=0;r<nc;r++) for (int c=0;c<nc;c++) ds >> Ri(r,c); }
      // MEANS
      std::getline( IN , dline );
      { std::istringstream ds(dline); std::string tag; ds >> tag;
        for (int s=0;s<nc;s++) ds >> means_i[s]; }
      // SDS
      std::getline( IN , dline );
      { std::istringstream ds(dline); std::string tag; ds >> tag;
        for (int s=0;s<nc;s++) ds >> sds_i[s]; }

      if ( trace_norm )
        {
          double ts = Si.trace(); if ( ts > 0 ) Si /= ts;
          double tr = Ri.trace(); if ( tr > 0 ) Ri /= tr;
        }

      S_sum += Si;
      R_sum += Ri;
      for ( int s=0; s<nc; s++ ) { means_sum[s] += means_i[s]; sds_sum[s] += sds_i[s]; }
      ++n_indiv;

      logger << "  GED group: read " << id_i << " (n_S=" << n_S_i << ", n_R=" << n_R_i << ")\n";
    }
  IN.close();

  if ( n_indiv < min_n )
    Helper::halt( "GED group: only " + Helper::int2str(n_indiv) + " individuals; min-n=" + Helper::int2str(min_n) );

  logger << "  GED group: " << n_indiv << " individuals\n";

  Eigen::MatrixXd S_grp = S_sum / double( n_indiv );
  Eigen::MatrixXd R_grp = R_sum / double( n_indiv );
  std::vector<double> means_grp( nc ) , sds_grp( nc );
  for ( int s=0; s<nc; s++ )
    { means_grp[s] = means_sum[s] / n_indiv; sds_grp[s] = sds_sum[s] / n_indiv; }

  // GED
  ged_t ged;
  ged.covar( S_grp , R_grp );
  ged.n_S = n_indiv;
  ged.n_R = n_indiv;
  if ( reg > 0.0 ) ged.regularize_R( reg );
  ged.calc();

  if ( ! ged.solved )
    Helper::halt( "GED group: solver failed; try reg=0.01 or larger" );

  // Write solution in descending eigenvalue order
  std::vector<int> si = ged.sorted_indices();
  if ( nc_out < 0 || nc_out > nc ) nc_out = nc;

  Eigen::VectorXd L_out( nc_out );
  Eigen::MatrixXd W_out( nc , nc_out );
  for ( int k=0; k<nc_out; k++ )
    {
      L_out[k] = ged.L[ si[k] ];
      W_out.col(k) = ged.W.col( si[k] );
    }

  ged_group_t::write_solution( sol_file , ch_names , W_out , L_out ,
                               mode , nb_f , nb_fwhm , env_lwr , env_upr ,
                               z_scored , means_grp , sds_grp );

  logger << "  GED group: solution written to " << sol_file << "\n";

  // Log eigenvalues
  for ( int k=0; k<nc_out; k++ )
    logger << "  GED group: component " << k+1 << " lambda=" << L_out[k] << "\n";
}


// ---------------------------------------------------------------------------
// ged_group_t::write_solution()
// ---------------------------------------------------------------------------

void ged_group_t::write_solution( const std::string & fname ,
                                  const std::vector<std::string> & ch_names ,
                                  const Eigen::MatrixXd & W ,
                                  const Eigen::VectorXd & L ,
                                  ged_input_mode_t mode ,
                                  double nb_f , double nb_fwhm ,
                                  double env_lwr , double env_upr ,
                                  bool z_scored ,
                                  const std::vector<double> & ch_means ,
                                  const std::vector<double> & ch_sds )
{
  const int nc     = (int)ch_names.size();
  const int nc_sol = (int)L.size();

  std::ofstream OUT( fname.c_str() );
  if ( ! OUT.good() )
    Helper::halt( "GED: cannot write solution to " + fname );

  OUT << std::setprecision(15);

  std::string ch_str;
  for ( int s=0; s<nc; s++ ) ch_str += (s?",":"") + ch_names[s];

  OUT << "GED_SOLUTION"
      << " nc="      << nc
      << " nc_sol="  << nc_sol
      << " mode="    << ged_input_mode_str(mode)
      << " nb_f="    << nb_f
      << " nb_fwhm=" << nb_fwhm
      << " env_lwr=" << env_lwr
      << " env_upr=" << env_upr
      << " z="       << (z_scored ? 1 : 0)
      << " ch="      << ch_str
      << "\n";

  // Eigenvalues (descending)
  OUT << "L";
  for ( int k=0; k<nc_sol; k++ ) OUT << " " << L[k];
  OUT << "\n";

  // Means
  OUT << "MEANS";
  for ( int s=0; s<nc; s++ )
    OUT << " " << ( z_scored && s<(int)ch_means.size() ? ch_means[s] : 0.0 );
  OUT << "\n";

  // SDs
  OUT << "SDS";
  for ( int s=0; s<nc; s++ )
    OUT << " " << ( z_scored && s<(int)ch_sds.size() ? ch_sds[s] : 1.0 );
  OUT << "\n";

  // Eigenvectors column-major (nc rows x nc_sol cols)
  OUT << "W";
  for ( int c=0; c<nc_sol; c++ )
    for ( int r=0; r<nc; r++ )
      OUT << " " << W(r,c);
  OUT << "\n";

  OUT.close();
}


// ---------------------------------------------------------------------------
// ged_group_t::load_solution()
// ---------------------------------------------------------------------------

void ged_group_t::load_solution( const std::string & fname )
{
  std::ifstream IN( fname.c_str() );
  if ( ! IN.good() )
    Helper::halt( "GED load: cannot open solution file " + fname );

  std::string line;

  // Header line
  if ( ! std::getline( IN , line ) || line.substr(0,12) != "GED_SOLUTION" )
    Helper::halt( "GED load: not a valid solution file: " + fname );

  auto kv = ged_parse_kv( line );
  const int nc     = std::stoi( kv["nc"] );
  const int nc_sol = std::stoi( kv["nc_sol"] );

  std::vector<std::string> ch_names;
  { std::istringstream chs( kv["ch"] ); std::string ch;
    while ( std::getline(chs,ch,',') ) ch_names.push_back(ch); }

  ged_input_mode_t mode = ged_input_mode_from_str( kv["mode"] );
  double nb_f    = std::stod( kv["nb_f"] );
  double nb_fwhm = std::stod( kv["nb_fwhm"] );
  double env_lwr = std::stod( kv["env_lwr"] );
  double env_upr = std::stod( kv["env_upr"] );
  bool z_scored  = ( kv["z"] == "1" );

  // L
  Eigen::VectorXd L( nc_sol );
  std::getline( IN , line );
  { std::istringstream ds(line); std::string tag; ds >> tag;
    for (int k=0;k<nc_sol;k++) ds >> L[k]; }

  // MEANS
  std::vector<double> ch_means(nc);
  std::getline( IN , line );
  { std::istringstream ds(line); std::string tag; ds >> tag;
    for (int s=0;s<nc;s++) ds >> ch_means[s]; }

  // SDS
  std::vector<double> ch_sds(nc);
  std::getline( IN , line );
  { std::istringstream ds(line); std::string tag; ds >> tag;
    for (int s=0;s<nc;s++) ds >> ch_sds[s]; }

  // W column-major
  Eigen::MatrixXd W( nc , nc_sol );
  std::getline( IN , line );
  { std::istringstream ds(line); std::string tag; ds >> tag;
    for (int c=0;c<nc_sol;c++) for (int r=0;r<nc;r++) ds >> W(r,c); }

  IN.close();

  // Store into ged_t statics
  ged_t::group_channels      = ch_names;
  ged_t::group_W             = W;
  ged_t::group_L             = L;
  ged_t::group_ch_means      = ch_means;
  ged_t::group_ch_sds        = ch_sds;
  ged_t::group_input_mode    = mode;
  ged_t::group_nb_f          = nb_f;
  ged_t::group_nb_fwhm       = nb_fwhm;
  ged_t::group_env_lwr       = env_lwr;
  ged_t::group_env_upr       = env_upr;
  ged_t::group_z_scored      = z_scored;
  ged_t::group_solution_fname = fname;
  ged_t::has_group_solution  = true;

  logger << "  GED: loaded group solution from " << fname
         << " (" << nc << " channels, " << nc_sol << " components)\n";
}
