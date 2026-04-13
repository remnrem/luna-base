
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

#ifndef __GED_H__
#define __GED_H__

struct edf_t;
struct param_t;
struct signal_list_t;

#include "stats/Eigen/Dense"
#include <vector>
#include <map>
#include <string>
#include <fstream>


// ---------------------------------------------------------------------------
// Input transform mode
// ---------------------------------------------------------------------------

enum class ged_input_mode_t { RAW = 0 , NB = 1 , ENV = 2 };

std::string ged_input_mode_str( ged_input_mode_t );
ged_input_mode_t ged_input_mode_from_str( const std::string & );


// ---------------------------------------------------------------------------
// Core GED struct
// ---------------------------------------------------------------------------

struct ged_t
{

  ged_t() : n_S(0), n_R(0), largest_idx(-1), solved(false) { }

  // Entry point 1: provide raw data matrices (centers + computes covariance internally)
  void data( const Eigen::MatrixXd & Sd , const Eigen::MatrixXd & Rd );

  // Entry point 2: provide pre-computed covariance matrices
  void covar( const Eigen::MatrixXd & S_ , const Eigen::MatrixXd & R_ )
  {
    S = S_;
    R = R_;
  }

  // Core solver
  void calc();

  // Tikhonov regularization of R:  R += alpha * (trace(R)/p) * I
  // Call between covar()/data() and calc()
  void regularize_R( const double alpha );

  // --- Signal transformation helpers (static, in-place on data matrix) ---

  // Apply narrow Gaussian bandpass to each column (uses existing ngaus)
  static void apply_nb_filter( Eigen::MatrixXd & X , int sr ,
                                double f , double fwhm );

  // Apply bandpass + Hilbert magnitude envelope to each column
  static void apply_env_filter( Eigen::MatrixXd & X , int sr ,
                                 double lwr , double upr ,
                                 double ripple , double tw );

  // Z-score each column in place, optionally storing per-column mean/sd
  static void zscore_columns( Eigen::MatrixXd & X ,
                               std::vector<double> * means = nullptr ,
                               std::vector<double> * sds   = nullptr );

  // Apply pre-computed z-score (for epochs / projection mode)
  static void apply_zscore_stored( Eigen::MatrixXd & X ,
                                   const std::vector<double> & means ,
                                   const std::vector<double> & sds );

  // --- Output helpers ---

  // Sorted component indices: sorted_indices()[0] = index of largest eigenvalue, etc.
  std::vector<int> sorted_indices() const;

  // lambda_e / sum(lambda)
  double eigenvalue_ratio( const int e ) const;

  // Forward model for component e using covariance C; *maxch = index of max |weight|
  Eigen::VectorXd map( const int e , const Eigen::MatrixXd & C , int * maxch );

  // Project data D (timepoints x channels) through eigenvector e, sign-fixed against maxch
  Eigen::VectorXd time_series( const int e , const Eigen::MatrixXd & D , const int maxch );

  // Focality of a spatial map: 1 - entropy of normalized |weights|
  static double focality( const Eigen::VectorXd & spatialMap );

  // Anterior-posterior index from channel locations (NaN if map empty)
  static double ap_index( const Eigen::VectorXd & spatialMap ,
                          const std::vector<std::string> & ch_names ,
                          const std::map<std::string,double> & chan_ap );

  // Lateralization index: (R_wt - L_wt) / (R_wt + L_wt)  (NaN if no lateralized channels)
  static double lat_index( const Eigen::VectorXd & spatialMap ,
                           const std::vector<std::string> & ch_names ,
                           const std::map<std::string,int> & chan_lat );

  // --- Members ---

  Eigen::MatrixXd S;   // (nc x nc) signal covariance
  Eigen::MatrixXd R;   // (nc x nc) reference covariance

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es;

  Eigen::MatrixXd W;   // (nc x nc) eigenvectors (columns), in solver order
  Eigen::VectorXd L;   // (nc)      eigenvalues,   in solver order

  int largest_idx;     // column index of max eigenvalue in W / L
  int n_S;             // timepoints used to build S
  int n_R;             // timepoints used to build R
  bool solved;         // true after calc() completes successfully


  // ---------------------------------------------------------------------------
  // Static group-solution (PSC pattern): loaded once, reused across individuals
  // ---------------------------------------------------------------------------

  static bool has_group_solution;
  static std::string group_solution_fname;
  static std::vector<std::string> group_channels;
  static Eigen::MatrixXd group_W;       // (nc x nc_sol)
  static Eigen::VectorXd group_L;       // (nc_sol) descending
  static std::vector<double> group_ch_means;
  static std::vector<double> group_ch_sds;
  static ged_input_mode_t group_input_mode;
  static double group_nb_f;
  static double group_nb_fwhm;
  static double group_env_lwr;
  static double group_env_upr;
  static bool group_z_scored;

  static void clear_group_solution();

};


// ---------------------------------------------------------------------------
// Group workflow struct
// ---------------------------------------------------------------------------

struct ged_group_t
{
  // Write per-individual S/R covariance matrices to binary accumulation file
  // (appends; creates with header if new)
  static void save_covar( const std::string & fname ,
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
                          const std::vector<double> & col_sds );

  // Read accumulated file, compute group-average S/R, run GED, write solution file
  static void run_group( param_t & param );

  // Write complete group solution to .ged file
  static void write_solution( const std::string & fname ,
                              const std::vector<std::string> & ch_names ,
                              const Eigen::MatrixXd & W ,
                              const Eigen::VectorXd & L ,
                              ged_input_mode_t mode ,
                              double nb_f , double nb_fwhm ,
                              double env_lwr , double env_upr ,
                              bool z_scored ,
                              const std::vector<double> & ch_means ,
                              const std::vector<double> & ch_sds );

  // Load solution file into ged_t static members
  static void load_solution( const std::string & fname );

};


// ---------------------------------------------------------------------------
// Top-level wrapper and helper functions (ged.cpp)
// ---------------------------------------------------------------------------

// Main entry point — dispatches to run mode or projection mode
void ged_wrapper( edf_t & , param_t & );

// Mode 1: narrowband vs. broadband (filter-based, no annotation time-lock)
void ged_runmode1( edf_t & , param_t & , Eigen::MatrixXd & X , int sr );

// Mode 2: annotation-locked S vs. reference R
void ged_runmode2( edf_t & , param_t & , Eigen::MatrixXd & X ,
                   const std::vector<uint64_t> * tp , int sr );

// Projection mode: apply pre-loaded group solution to this individual
void ged_apply_group( edf_t & , param_t & , Eigen::MatrixXd & X ,
                      const std::vector<uint64_t> * tp , int sr );

// Apply input transform in-place; stores col_means/col_sds if z-scoring
void ged_apply_input_transform( Eigen::MatrixXd & X ,
                                ged_input_mode_t mode , int sr ,
                                double nb_f , double nb_fwhm ,
                                double env_lwr , double env_upr ,
                                double env_ripple , double env_tw ,
                                bool z_score ,
                                std::vector<double> * col_means = nullptr ,
                                std::vector<double> * col_sds   = nullptr );

// Write all component outputs for nc_out components (stage = "" for whole-recording)
void ged_write_component_outputs( const ged_t & ged ,
                                  const std::vector<std::string> & ch_names ,
                                  const Eigen::MatrixXd & S ,
                                  int nc_out ,
                                  const std::string & stage_label ,
                                  const std::map<std::string,double> * chan_ap ,
                                  const std::map<std::string,int>    * chan_lat );

// Compute per-epoch discriminant power: w^T C_epoch w
// Returns vector of (display_epoch, power) pairs
std::vector<std::pair<int,double>> ged_epoch_power( edf_t & edf ,
                                                    const signal_list_t & signals ,
                                                    const Eigen::VectorXd & w ,
                                                    ged_input_mode_t mode , int sr ,
                                                    double nb_f , double nb_fwhm ,
                                                    double env_lwr , double env_upr ,
                                                    double env_ripple , double env_tw ,
                                                    const std::vector<double> & col_means ,
                                                    const std::vector<double> & col_sds );

// Read channel locations file: col1=label, col2=ap_coord, col3=lat_class (-1/0/1)
bool ged_read_clocs( const std::string & fname ,
                     std::map<std::string,double> * chan_ap ,
                     std::map<std::string,int>    * chan_lat );

// Build row-to-stage mapping from timepoints and epoch annotations
// Returns vector (one per row of X) of stage label strings ("N1","N2","N3","R","W","?")
std::vector<std::string> ged_build_row_stages( edf_t & edf ,
                                               const std::vector<uint64_t> & tp );


#endif
