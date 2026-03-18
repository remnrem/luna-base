
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

    // per-window artifact flags
    double resp_flatline_floor;
    double resp_flatline_prop;
    double resp_clip_prop;
    double resp_jump_mad_th;

    // channel-level bad rules (from epoch flags)
    double resp_flag_prop;
    double resp_flag_run;

    bool resp_add_channel;
    std::string resp_channel_label;

    //
    // spo2 (param/output domain: oxy/OXY)
    //

    void do_spo2( signal_list_t & signals );
    int    spo2_min_sr;
    double spo2_window_dur;
    double spo2_window_inc;

    // per-epoch thresholds
    double spo2_range_lo;
    double spo2_range_hi;
    double spo2_range_prop;
    double spo2_flatline_th;
    int    spo2_flatline_k;
    double spo2_jump_th;
    double spo2_invalid_floor;
    double spo2_missing_prop;

    // channel-level bad rules
    double spo2_flag_prop;
    double spo2_flag_run;

    //
    // eeg
    //

    void do_eeg( signal_list_t & signals );
    int    eeg_min_sr;
    double eeg_window_dur;
    double eeg_window_inc;

    // flatline
    double eeg_flat_std_th;      // SD threshold (µV), default 2.0
    double eeg_flat_deriv_prop;  // near-zero derivative proportion, default 0.8
    double eeg_flat_deriv_eps;   // near-zero derivative threshold, default 1e-4

    // clipping
    double eeg_clip_prop;        // proportion of clipped samples, default 0.01

    // extreme amplitude
    double eeg_amp_th;           // |EEG| threshold (µV), default 500
    double eeg_amp_prop;         // proportion of samples exceeding threshold, default 0.05

    // spectral: HF contamination (P[hf_lo,hf_hi] / P[sig_lo,sig_hi])
    double eeg_hf_lo;            // HF band lower (Hz), default 20
    double eeg_hf_hi;            // HF band upper (Hz), default 40
    double eeg_sig_lo;           // signal band lower (Hz), default 0.5
    double eeg_sig_hi;           // signal band upper (Hz), default 20
    double eeg_hf_ratio_th;      // HF/signal power ratio threshold, default 1.5

    // spectral: line noise
    double eeg_ln_bw;            // ± Hz around 50/60 Hz, default 2
    double eeg_ln_lo;            // broadband lower (Hz), default 0.5
    double eeg_ln_hi;            // broadband upper (Hz), default 40
    double eeg_ln_ratio_th;      // line/broadband ratio threshold, default 0.5

    // spectral peakedness (output only; thresholds TBD)
    double eeg_spectral_peakedness_th;   // SPK flag threshold (0 = no flagging yet)
    double eeg_spectral_skewness_th;     // KURT flag threshold (0 = no flagging yet)

    // Hjorth pre-filter (IIR Butterworth, applied only before Hjorth computation)
    double eeg_hjorth_lo;        // filter lower cutoff (Hz), default 0.5
    double eeg_hjorth_hi;        // filter upper cutoff (Hz), default 40
    int    eeg_hjorth_order;     // Butterworth order, default 2

    // Hjorth outlier (robust z vs channel-night baseline)
    double eeg_hjorth_z_th;      // |robust z| threshold for activity/complexity, default 5

    // channel-level BAD rules
    double eeg_flag_prop;   // channel BAD if > this fraction of epochs bad, default 0.5
    double eeg_flag_run;      // channel BAD if max contiguous bad run >= this (sec), default 1200

    //
    // emg
    //

    void do_emg( signal_list_t & signals );
    int    emg_min_sr;
    double emg_window_dur;
    double emg_window_inc;

    // flatline
    double emg_flat_std_th;      // SD threshold (µV), default 0.5
    double emg_flat_deriv_prop;  // proportion of near-zero derivatives, default 0.8
    double emg_flat_deriv_eps;   // near-zero derivative threshold, default 1e-4

    // clipping
    double emg_clip_prop;        // proportion of clipped samples, default 0.01

    // broadband RMS
    double emg_rms_th;           // RMS amplitude threshold (µV), default 150

    // impulse artifact
    int    emg_impulse_n;        // min count of |z| > emg_impulse_z, default 5
    double emg_impulse_z;        // z-score threshold, default 8.0

    // channel-level BAD rules
    double emg_flag_prop;   // channel BAD if > this fraction of epochs bad, default 0.5
    double emg_flag_run;      // channel BAD if max contiguous bad run >= this (sec), default 1200

    // line noise screening (shared helper, used by EMG and others)
    double emg_ln_bw;            // ± Hz around 50/60 Hz, default 2
    double emg_ln_broad_lo;      // broadband lower bound, default 10 Hz
    double emg_ln_broad_hi;      // broadband upper bound, default 100 Hz (Nyquist-capped)
    double emg_ln_ratio_th;      // epoch flag threshold P_line/P_broad, default 0.5

    // Returns P_line / P_broad for one epoch's data, or -1 if not assessable
    double line_noise_ratio( const std::vector<double> & data , double Fs ) const;

    // Ensure signal slot s is in target_unit (e.g. "uV" or "mV").
    // Auto-converts between uV / mV / V as needed; warns if unit is unrecognised.
    void ensure_units( int s , const std::string & target_unit );

    //
    // ecg
    //

    void do_ecg( signal_list_t & signals );
    int    ecg_min_sr;
    double ecg_window_dur;
    double ecg_window_inc;

    // flatline
    double ecg_flat_std_th;      // SD threshold (mV), default 0.02
    double ecg_flat_deriv_prop;  // near-zero derivative proportion, default 0.8
    double ecg_flat_deriv_eps;   // near-zero derivative threshold, default 1e-4

    // clipping
    double ecg_clip_prop;        // proportion of clipped samples, default 0.01

    // HR physiologic limits
    double ecg_hr_min;           // minimum plausible HR (bpm), default 25
    double ecg_hr_max;           // maximum plausible HR (bpm), default 220

    // RR interval plausibility
    double ecg_rr_min_ms;        // minimum plausible RR (ms), default 300
    double ecg_rr_max_ms;        // maximum plausible RR (ms), default 2000
    double ecg_rr_flag_prop;      // max fraction of RR outside limits, default 0.20

    // beat detection failure
    int    ecg_min_beats;        // min beats per epoch to be "detectable", default 5

    // line noise (separate from BAD, consistent with other domains)
    double ecg_ln_bw;            // ± Hz around 50/60 Hz, default 2
    double ecg_sig_lo;           // cardiac signal band lower (Hz), default 5
    double ecg_sig_hi;           // cardiac signal band upper (Hz), default 25
    double ecg_ln_ratio_th;      // P_line/P_signal ratio threshold, default 0.40

    // channel-level BAD rules
    double ecg_flag_prop;   // channel BAD if > this fraction of epochs bad, default 0.50
    double ecg_flag_run;      // channel BAD if max contiguous bad run >= this (sec), default 1200

    // R-peak annotations
    bool ecg_add_peaks;      // if true, add point annotations for each detected R-peak

    //
    // eog
    //

    void do_eog( signal_list_t & signals );
    int    eog_min_sr;
    double eog_window_dur;
    double eog_window_inc;

    // flatline
    double eog_flat_std_th;      // SD threshold (µV), default 3
    double eog_flat_deriv_prop;  // near-zero derivative proportion, default 0.8
    double eog_flat_deriv_eps;   // near-zero derivative threshold, default 1e-4

    // clipping
    double eog_clip_prop;        // proportion of clipped samples, default 0.01

    // extreme amplitude
    double eog_amp_th;           // |EOG| threshold (µV), default 700
    double eog_amp_prop;         // proportion of samples exceeding threshold, default 0.05

    // spectral: HF contamination (P[hf_lo,hf_hi] / P[sig_lo,sig_hi])
    double eog_hf_lo;            // HF band lower (Hz), default 20
    double eog_hf_hi;            // HF band upper (Hz), default 40
    double eog_sig_lo;           // signal band lower (Hz), default 0.5
    double eog_sig_hi;           // signal band upper (Hz), default 20
    double eog_hf_ratio_th;      // HF/signal power ratio threshold, default 1.5

    // spectral: line noise (same 50/60 Hz logic as EMG, different broadband)
    double eog_ln_bw;            // ± Hz around 50/60 Hz, default 2
    double eog_ln_lo;            // broadband lower (Hz), default 0.5
    double eog_ln_hi;            // broadband upper (Hz), default 40
    double eog_ln_ratio_th;      // line/broadband ratio threshold, default 0.5

    // Hjorth pre-filter (IIR Butterworth bandpass, applied only before Hjorth computation)
    double eog_hjorth_lo;        // filter lower cutoff (Hz), default 0.5
    double eog_hjorth_hi;        // filter upper cutoff (Hz), default 20
    int    eog_hjorth_order;     // Butterworth order, default 2

    // Hjorth outlier (robust z vs channel-night baseline)
    double eog_hjorth_z_th;      // |robust z| threshold for activity/complexity, default 5

    // channel-level BAD rules
    double eog_flag_prop;   // channel BAD if > this fraction of epochs bad, default 0.5
    double eog_flag_run;      // channel BAD if max contiguous bad run >= this (sec), default 1200

    //
    // unified annotation settings (all domains)
    //

    bool add_annots;
    std::string annot_prefix;
    bool annot_by_channel;
    bool annot_show_domain;

  };
   
}

#endif

