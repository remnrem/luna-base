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

#ifndef __AROUSALS2_H__
#define __AROUSALS2_H__

#include "stats/Eigen/Dense"
#include "intervals/intervals.h"

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

struct edf_t;
struct param_t;
struct signal_list_t;
class FFT;

struct arousals2_params_t {
  bool do_nrem = true;
  bool do_rem = true;
  bool verbose = false;
  bool nrem_emg_only = false;

  double epoch_inc = 0.5;
  double pre_sec = 10.0;
  double pre_gap_sec = 2.0;
  double post_sec = 0;

  double eeg_tilt_rise_th = 4.0;
  double eeg_fast_rise_th = 2.5;
  double eeg_tilt_rise_th_nrem = 4.0;
  double eeg_tilt_rise_th_rem = 4.0;
  double eeg_fast_rise_th_nrem = 2.5;
  double eeg_fast_rise_th_rem = 2.5;
  double eeg_active_th = 0.20;
  double emg_rise_th = 1.0;
  double emg_active_th = 0.5;
  double active_gap_sec = 1.0;
  double delta_artifact_th = 5.0;
  double artifact_frac = 0.8;
  double eeg_artifact_th = 0.7;

  double emg_median_sec = 2.5;

  double min_dur = 2.0;
  double max_dur = 15.0;
  double long_dur = 30.0;
  double arousal_dur = 3.0;
  double duration_level_frac = 0.5;
  double merge_gap_sec = 2.5;
  double pre_sleep_sec = 10.0;
  double emg_rise_min_dur = 1.0;
  double emg_rise_buffer = 2.0;
};

struct arousals2_t {

  arousals2_t( edf_t & edf , param_t & param );

 private:

  enum {
    IDX_DELTA = 0,
    IDX_THETA = 1,
    IDX_ALPHA = 2,
    IDX_BETA = 3,
    IDX_EMG = 4,
    IDX_EEG_TILT_RISE = 5,
    IDX_EEG_FAST_RISE = 6,
    IDX_EMG_RISE = 7,
    IDX_DELTA_ART = 8,
    NFTR = 9
  };

  void build_ftr_matrix( edf_t & edf ,
                         const signal_list_t & eeg_signals ,
                         const signal_list_t & emg_signals ,
                         Eigen::MatrixXd & Xeeg ,
                         Eigen::MatrixXd & Xemg ,
                         std::vector<int> * state ,
                         std::vector<int> * sequence ,
                         std::vector<double> * sec ,
                         const int wake_bridge ,
                         const double epoch_inc );

  Eigen::VectorXd calc_eeg_ftrs( const Eigen::MatrixXd & X , FFT & fftseg );
  Eigen::VectorXd calc_emg_ftrs( const Eigen::MatrixXd & X , const std::vector<double> & thr );

  Eigen::MatrixXd process_ftr_matrix( Eigen::MatrixXd * Xeeg ,
                                      Eigen::MatrixXd * Xemg ,
                                      const std::vector<int> & st ,
                                      const std::vector<int> & seq ,
                                      const std::vector<double> & sec ,
                                      const arousals2_params_t & p );

  std::vector<std::vector<std::vector<Eigen::VectorXd> > >
  assemble( const Eigen::MatrixXd & Xftr ,
            const std::vector<int> & state ,
            const std::vector<int> & seq ,
            const std::vector<double> & sec ,
            std::vector<std::vector<std::vector<double> > > * tt );

  std::map<std::string,std::set<interval_t> >
  event_heuristic( const std::vector<std::vector<std::vector<Eigen::VectorXd> > > & X ,
                   const std::vector<std::vector<std::vector<double> > > & tt ,
                   const std::vector<std::vector<std::vector<double> > > & eeg_artifact ,
                   const arousals2_params_t & p );

  void add_channels( const std::vector<std::vector<std::vector<Eigen::VectorXd> > > & X ,
                     const std::vector<std::vector<std::vector<double> > > & tt ,
                     const std::string & ch_prefix ,
                     const signal_list_t & eeg_signals ,
                     const double epoch_win ,
                     std::vector<std::vector<std::vector<double> > > * eeg_artifact );

  Eigen::VectorXd robust_mad_norm( const Eigen::VectorXd & x ,
                                   const std::vector<int> & st ,
                                   const int target_state );

  std::vector<std::pair<int,int> >
  merge_events_with_gap_sorted( const std::vector<std::pair<int,int> > & events ,
                                const int max_gap );

  std::vector<std::pair<int,int> >
  mask_to_intervals( const std::vector<bool> & mask );

  double local_median( const std::vector<double> & x , int a , int b ) const;
  double local_mean( const std::vector<double> & x , int a , int b ) const;

  edf_t * parent;
  int sr = 0;
};

#endif
