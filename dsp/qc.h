
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

#ifndef __QC_H__
#define __QC_H__

struct edf_t;
struct param_t;
struct signal_list_t;

#include <string>

namespace dsptools
{
  
  struct qc_t {
    qc_t( edf_t & edf , param_t & param );

  private:
    
    edf_t & edf;

    //
    // general
    //
    
    bool by_epoch;
    
    //
    // resp
    //

    void do_resp( signal_list_t & signals );
    int resp_min_sr;
    double resp_th;
    double resp_prop_th;
    double resp_window_dur;
    double resp_window_inc;
    double resp_p1_lwr;
    double resp_p1_upr;
    double resp_p2_lwr;
    double resp_p2_upr;
    double resp_epsilon;   
    bool resp_add_annot;
    std::string resp_annot_label;
    bool resp_add_channel;    
    std::string resp_channel_label;

    //
    // eeg
    //
    
    void do_eeg( signal_list_t & signals );
    int eeg_min_sr;
    double eeg_window_dur;
    double eeg_window_inc;
    double eeg_min_amp_th;
    double eeg_max_amp_th;
    double eeg_spectral_peakedness_th;
    double eeg_spectral_skewness_th;
    double eeg_h1_min, eeg_h1_max;
    double eeg_h2_min, eeg_h2_max;
    double eeg_h3_min, eeg_h3_max;    
    double eeg_eps;
    double eeg_skew_th, eeg_kurt_th;
    double eeg_fft_seg_sec, eeg_fft_inc_sec;
    
  };
   
}

#endif

