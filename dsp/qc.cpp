
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

#include "dsp/qc.h"
#include "param.h"

#include "dsp/fir.h"
#include "dsp/resample.h"
#include "dsp/iir.h"
#include "dsp/ecgsuppression.h"
#include "stats/eigen_ops.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "fftw/fftwrap.h"

extern logger_t logger;

extern writer_t writer;

dsptools::qc_t::qc_t( edf_t & edf1 , param_t & param ) : edf(edf1) 
{
  
  //
  // break-down by signal types
  //

  //
  // get signals
  //

  const bool no_annotations = true;
  
  const std::string resp_signal_labels = param.has( "resp" ) ? param.value( "resp" ) : "" ;
  const std::string spo2_signal_labels = param.has( "oxy"  ) ? param.value( "oxy"  ) : "" ;
  const std::string eeg_signal_labels  = param.has( "eeg"  ) ? param.value( "eeg"  ) : "" ;
  const std::string emg_signal_labels  = param.has( "emg"  ) ? param.value( "emg"  ) : "" ;
  const std::string ecg_signal_labels  = param.has( "ecg"  ) ? param.value( "ecg"  ) : "" ;
  const std::string eog_signal_labels  = param.has( "eog"  ) ? param.value( "eog"  ) : "" ;

  signal_list_t resp_signals = edf.header.signal_list( resp_signal_labels , no_annotations );
  signal_list_t spo2_signals = edf.header.signal_list( spo2_signal_labels , no_annotations );
  signal_list_t eeg_signals  = edf.header.signal_list( eeg_signal_labels  , no_annotations );
  signal_list_t emg_signals  = edf.header.signal_list( emg_signal_labels  , no_annotations );
  signal_list_t ecg_signals  = edf.header.signal_list( ecg_signal_labels  , no_annotations );
  signal_list_t eog_signals  = edf.header.signal_list( eog_signal_labels  , no_annotations );

  
  //
  // parameters
  //

  // SNR threshold for respiratory signals

  resp_th = param.has( "resp-snr-th" ) ? param.requires_dbl( "resp-snr-th" ) : 10 ; 
  resp_window_dur = param.has( "resp-win" ) ? param.requires_dbl( "resp-win" ) : 120 ;
  resp_window_inc = param.has( "resp-inc" ) ? param.requires_dbl( "resp-inc" ) : 10 ;

  // allow 10% of night to be 'bad' (noise-level 2 or worse) 
  
  resp_p1_lwr = param.has( "resp-p1-lwr" ) ? param.requires_dbl( "resp-p1-lwr" ) : 0.1 ;
  resp_p1_upr = param.has( "resp-p1-upr" ) ? param.requires_dbl( "resp-p1-upr" ) : 1 ;
  resp_p2_lwr = param.has( "resp-p2-lwr" ) ? param.requires_dbl( "resp-p2-lwr" ) : 1 ;
  resp_p2_upr = param.has( "resp-p2-upr" ) ? param.requires_dbl( "resp-p2-upr" ) : 10 ;
    
  resp_min_sr      = param.has( "resp-min-sr"      ) ? param.requires_int( "resp-min-sr"      ) : 32 ;
  resp_hard_min_sr = param.has( "resp-hard-min-sr" ) ? param.requires_int( "resp-hard-min-sr" ) : 4  ;
  resp_epsilon     = param.has( "resp-epsilon"      ) ? param.requires_dbl( "resp-epsilon"      ) : 1e-8 ;

  // per-window artifact thresholds
  resp_flatline_floor  = param.has( "resp-flat-floor" )     ? param.requires_dbl( "resp-flat-floor" )     : 1e-4 ;
  resp_flatline_prop   = param.has( "resp-flat-prop" )      ? param.requires_dbl( "resp-flat-prop" )      : 0.5  ;
  resp_clip_prop       = param.has( "resp-clip-prop" )      ? param.requires_dbl( "resp-clip-prop" )      : 0.01 ;
  resp_jump_mad_th     = param.has( "resp-jump-th" )        ? param.requires_dbl( "resp-jump-th" )        : 5.0  ;

  // channel-level bad rules
  resp_flag_prop  = param.has( "resp-flag-prop" ) ? param.requires_dbl( "resp-flag-prop" ) : 0.5   ;
  resp_flag_run     = param.has( "resp-flag-run" )        ? param.requires_dbl( "resp-flag-run" )        : 3600.0 ;

  //
  // eeg
  //

  eeg_min_sr         = param.has( "eeg-min-sr"     ) ? param.requires_int( "eeg-min-sr"     ) : 100   ;
  eeg_window_dur     = param.has( "eeg-win"         ) ? param.requires_dbl( "eeg-win"         ) : 30.0  ;
  eeg_window_inc     = param.has( "eeg-inc"         ) ? param.requires_dbl( "eeg-inc"         ) : 30.0  ;

  // flatline
  eeg_flat_std_th    = param.has( "eeg-flat-th"     ) ? param.requires_dbl( "eeg-flat-th"     ) : 2.0   ;
  eeg_flat_deriv_prop= param.has( "eeg-flat-prop"   ) ? param.requires_dbl( "eeg-flat-prop"   ) : 0.8   ;
  eeg_flat_deriv_eps = param.has( "eeg-flat-eps"    ) ? param.requires_dbl( "eeg-flat-eps"    ) : 1e-4  ;

  // clipping
  eeg_clip_prop      = param.has( "eeg-clip-prop"   ) ? param.requires_dbl( "eeg-clip-prop"   ) : 0.01  ;

  // extreme amplitude
  eeg_amp_th         = param.has( "eeg-amp-th"      ) ? param.requires_dbl( "eeg-amp-th"      ) : 500.0 ;
  eeg_amp_prop       = param.has( "eeg-amp-prop"    ) ? param.requires_dbl( "eeg-amp-prop"    ) : 0.05  ;

  // spectral: HF contamination
  eeg_hf_lo          = param.has( "eeg-hf-lo"       ) ? param.requires_dbl( "eeg-hf-lo"       ) : 20.0  ;
  eeg_hf_hi          = param.has( "eeg-hf-hi"       ) ? param.requires_dbl( "eeg-hf-hi"       ) : 40.0  ;
  eeg_sig_lo         = param.has( "eeg-sig-lo"      ) ? param.requires_dbl( "eeg-sig-lo"      ) : 0.5   ;
  eeg_sig_hi         = param.has( "eeg-sig-hi"      ) ? param.requires_dbl( "eeg-sig-hi"      ) : 20.0  ;
  eeg_hf_ratio_th    = param.has( "eeg-hf-th"       ) ? param.requires_dbl( "eeg-hf-th"       ) : 1.5   ;

  // spectral: line noise
  eeg_ln_bw          = param.has( "eeg-ln-bw"       ) ? param.requires_dbl( "eeg-ln-bw"       ) : 2.0   ;
  eeg_ln_lo          = param.has( "eeg-ln-lo"       ) ? param.requires_dbl( "eeg-ln-lo"       ) : 0.5   ;
  eeg_ln_hi          = param.has( "eeg-ln-hi"       ) ? param.requires_dbl( "eeg-ln-hi"       ) : 40.0  ;
  eeg_ln_ratio_th    = param.has( "eeg-ln-th"       ) ? param.requires_dbl( "eeg-ln-th"       ) : 0.3   ;

  // spectral peakedness
  eeg_spectral_peakedness_th = param.has( "eeg-peak-th" ) ? param.requires_dbl( "eeg-peak-th" ) : 15.0 ;
  eeg_spectral_skewness_th   = param.has( "eeg-kurt-th" ) ? param.requires_dbl( "eeg-kurt-th" ) : 5.0 ;

  // Hjorth pre-filter
  eeg_hjorth_lo      = param.has( "eeg-hjorth-lo"   ) ? param.requires_dbl( "eeg-hjorth-lo"   ) : 0.5  ;
  eeg_hjorth_hi      = param.has( "eeg-hjorth-hi"   ) ? param.requires_dbl( "eeg-hjorth-hi"   ) : 40.0 ;
  eeg_hjorth_order   = param.has( "eeg-hjorth-ord"  ) ? param.requires_int( "eeg-hjorth-ord"  ) : 2    ;
  eeg_hjorth_z_th    = param.has( "eeg-hjorth-z"    ) ? param.requires_dbl( "eeg-hjorth-z"    ) : 10.0 ;

  // channel-level BAD rules
  eeg_flag_prop = param.has( "eeg-flag-prop"    ) ? param.requires_dbl( "eeg-flag-prop"    ) : 0.5  ;
  eeg_flag_run    = param.has( "eeg-flag-run"     ) ? param.requires_dbl( "eeg-flag-run"     ) : 3600.0 ;
  
  //
  // spo2 (param/output domain: oxy/OXY)
  //

  spo2_min_sr      = param.has( "oxy-min-sr"          ) ? param.requires_int( "oxy-min-sr"          ) : 1     ;
  spo2_window_dur  = param.has( "oxy-win"              ) ? param.requires_dbl( "oxy-win"              ) : 30.0  ;
  spo2_window_inc  = param.has( "oxy-inc"              ) ? param.requires_dbl( "oxy-inc"              ) : 30.0  ;

  spo2_range_lo    = param.has( "oxy-range-lo"         ) ? param.requires_dbl( "oxy-range-lo"         ) : 50.0  ;
  spo2_range_hi    = param.has( "oxy-range-hi"         ) ? param.requires_dbl( "oxy-range-hi"         ) : 100.0 ;
  spo2_range_prop  = param.has( "oxy-range-prop"       ) ? param.requires_dbl( "oxy-range-prop"       ) : 0.1   ;
  spo2_flatline_th = param.has( "oxy-flat-th"          ) ? param.requires_dbl( "oxy-flat-th"          ) : 0.2   ;
  spo2_flatline_k  = param.has( "oxy-flat-k"           ) ? param.requires_int( "oxy-flat-k"           ) : 10    ;
  spo2_jump_th     = param.has( "oxy-jump-th"          ) ? param.requires_dbl( "oxy-jump-th"          ) : 20.0  ;
  spo2_invalid_floor = param.has( "oxy-invalid-floor"  ) ? param.requires_dbl( "oxy-invalid-floor"    ) : 0.0   ;
  spo2_missing_prop  = param.has( "oxy-missing-prop"   ) ? param.requires_dbl( "oxy-missing-prop"     ) : 0.1   ;

  spo2_flag_prop = param.has( "oxy-flag-prop" ) ? param.requires_dbl( "oxy-flag-prop"  ) : 0.5   ;
  spo2_flag_run    = param.has( "oxy-flag-run"        ) ? param.requires_dbl( "oxy-flag-run"         ) : 3600.0 ;

  //
  // emg
  //

  emg_min_sr         = param.has( "emg-min-sr"    ) ? param.requires_int( "emg-min-sr"    ) : 100  ;
  emg_window_dur     = param.has( "emg-win"       ) ? param.requires_dbl( "emg-win"       ) : 30.0 ;
  emg_window_inc     = param.has( "emg-inc"       ) ? param.requires_dbl( "emg-inc"       ) : 30.0 ;

  // flatline
  emg_flat_std_th    = param.has( "emg-flat-th"   ) ? param.requires_dbl( "emg-flat-th"   ) : 0.5  ;
  emg_flat_deriv_prop= param.has( "emg-flat-prop" ) ? param.requires_dbl( "emg-flat-prop" ) : 0.8  ;
  emg_flat_deriv_eps = param.has( "emg-flat-eps"  ) ? param.requires_dbl( "emg-flat-eps"  ) : 1e-4 ;

  // clipping
  emg_clip_prop      = param.has( "emg-clip-prop" ) ? param.requires_dbl( "emg-clip-prop" ) : 0.01 ;

  // broadband RMS
  emg_rms_th         = param.has( "emg-rms-th"    ) ? param.requires_dbl( "emg-rms-th"    ) : 150.0 ;

  // impulse artifact
  emg_impulse_n      = param.has( "emg-imp-n"     ) ? param.requires_int( "emg-imp-n"     ) : 5    ;
  emg_impulse_z      = param.has( "emg-imp-z"     ) ? param.requires_dbl( "emg-imp-z"     ) : 8.0  ;

  // channel-level BAD rules
  emg_flag_prop = param.has( "emg-flag-prop"  ) ? param.requires_dbl( "emg-flag-prop"  ) : 0.5  ;
  emg_flag_run    = param.has( "emg-flag-run"   ) ? param.requires_dbl( "emg-flag-run"   ) : 3600.0 ;

  // line-noise screening
  emg_ln_bw          = param.has( "emg-ln-bw"     ) ? param.requires_dbl( "emg-ln-bw"     ) : 2.0   ;
  emg_ln_broad_lo    = param.has( "emg-ln-lo"     ) ? param.requires_dbl( "emg-ln-lo"     ) : 10.0  ;
  emg_ln_broad_hi    = param.has( "emg-ln-hi"     ) ? param.requires_dbl( "emg-ln-hi"     ) : 100.0 ;
  emg_ln_ratio_th    = param.has( "emg-ln-th"     ) ? param.requires_dbl( "emg-ln-th"     ) : 0.3   ;

  //
  // ecg
  //

  ecg_min_sr         = param.has( "ecg-min-sr"    ) ? param.requires_int( "ecg-min-sr"    ) : 128  ;
  ecg_window_dur     = param.has( "ecg-win"       ) ? param.requires_dbl( "ecg-win"       ) : 30.0 ;
  ecg_window_inc     = param.has( "ecg-inc"       ) ? param.requires_dbl( "ecg-inc"       ) : 30.0 ;

  // flatline
  ecg_flat_std_th    = param.has( "ecg-flat-th"   ) ? param.requires_dbl( "ecg-flat-th"   ) : 0.02 ;
  ecg_flat_deriv_prop= param.has( "ecg-flat-prop" ) ? param.requires_dbl( "ecg-flat-prop" ) : 0.8  ;
  ecg_flat_deriv_eps = param.has( "ecg-flat-eps"  ) ? param.requires_dbl( "ecg-flat-eps"  ) : 1e-4 ;

  // clipping
  ecg_clip_prop      = param.has( "ecg-clip-prop" ) ? param.requires_dbl( "ecg-clip-prop" ) : 0.01 ;

  // HR / RR plausibility
  ecg_hr_min         = param.has( "ecg-hr-min"    ) ? param.requires_dbl( "ecg-hr-min"    ) : 40.0 ;
  ecg_hr_max         = param.has( "ecg-hr-max"    ) ? param.requires_dbl( "ecg-hr-max"    ) : 140.0;
  ecg_rr_min_ms      = param.has( "ecg-rr-min"    ) ? param.requires_dbl( "ecg-rr-min"    ) : 430.0;
  ecg_rr_max_ms      = param.has( "ecg-rr-max"    ) ? param.requires_dbl( "ecg-rr-max"    ) : 1500.0;
  ecg_rr_flag_prop   = param.has( "ecg-rr-prop"   ) ? param.requires_dbl( "ecg-rr-prop"   ) : 0.05 ;
  ecg_min_beats      = param.has( "ecg-min-beats" ) ? param.requires_int( "ecg-min-beats" ) : 10   ;

  // line noise
  ecg_ln_bw          = param.has( "ecg-ln-bw"     ) ? param.requires_dbl( "ecg-ln-bw"     ) : 2.0  ;
  ecg_sig_lo         = param.has( "ecg-sig-lo"    ) ? param.requires_dbl( "ecg-sig-lo"    ) : 5.0  ;
  ecg_sig_hi         = param.has( "ecg-sig-hi"    ) ? param.requires_dbl( "ecg-sig-hi"    ) : 25.0 ;
  ecg_ln_ratio_th    = param.has( "ecg-ln-th"     ) ? param.requires_dbl( "ecg-ln-th"     ) : 0.40 ;

  // channel BAD rules
  ecg_flag_prop = param.has( "ecg-flag-prop"  ) ? param.requires_dbl( "ecg-flag-prop"  ) : 0.50  ;
  ecg_flag_run    = param.has( "ecg-flag-run"   ) ? param.requires_dbl( "ecg-flag-run"   ) : 3600.0 ;

  // R-peak annotations
  ecg_add_peaks = param.has( "ecg-add-peaks" );

  //
  // eog
  //

  eog_min_sr         = param.has( "eog-min-sr"    ) ? param.requires_int( "eog-min-sr"    ) : 32    ;
  eog_window_dur     = param.has( "eog-win"       ) ? param.requires_dbl( "eog-win"       ) : 30.0  ;
  eog_window_inc     = param.has( "eog-inc"       ) ? param.requires_dbl( "eog-inc"       ) : 30.0  ;

  // flatline
  eog_flat_std_th    = param.has( "eog-flat-th"   ) ? param.requires_dbl( "eog-flat-th"   ) : 3.0   ;
  eog_flat_deriv_prop= param.has( "eog-flat-prop" ) ? param.requires_dbl( "eog-flat-prop" ) : 0.8   ;
  eog_flat_deriv_eps = param.has( "eog-flat-eps"  ) ? param.requires_dbl( "eog-flat-eps"  ) : 1e-4  ;

  // clipping
  eog_clip_prop      = param.has( "eog-clip-prop" ) ? param.requires_dbl( "eog-clip-prop" ) : 0.01  ;

  // extreme amplitude
  eog_amp_th         = param.has( "eog-amp-th"    ) ? param.requires_dbl( "eog-amp-th"    ) : 700.0 ;
  eog_amp_prop       = param.has( "eog-amp-prop"  ) ? param.requires_dbl( "eog-amp-prop"  ) : 0.05  ;

  // spectral: HF contamination
  eog_hf_lo          = param.has( "eog-hf-lo"     ) ? param.requires_dbl( "eog-hf-lo"     ) : 20.0  ;
  eog_hf_hi          = param.has( "eog-hf-hi"     ) ? param.requires_dbl( "eog-hf-hi"     ) : 40.0  ;
  eog_sig_lo         = param.has( "eog-sig-lo"    ) ? param.requires_dbl( "eog-sig-lo"    ) : 0.3   ;
  eog_sig_hi         = param.has( "eog-sig-hi"    ) ? param.requires_dbl( "eog-sig-hi"    ) : 20.0  ;
  eog_hf_ratio_th    = param.has( "eog-hf-th"     ) ? param.requires_dbl( "eog-hf-th"     ) : 1.5   ;

  // spectral: line noise
  eog_ln_bw          = param.has( "eog-ln-bw"     ) ? param.requires_dbl( "eog-ln-bw"     ) : 2.0   ;
  eog_ln_lo          = param.has( "eog-ln-lo"     ) ? param.requires_dbl( "eog-ln-lo"     ) : 0.3   ;
  eog_ln_hi          = param.has( "eog-ln-hi"     ) ? param.requires_dbl( "eog-ln-hi"     ) : 40.0  ;
  eog_ln_ratio_th    = param.has( "eog-ln-th"     ) ? param.requires_dbl( "eog-ln-th"     ) : 0.3   ;

  // Hjorth pre-filter (IIR Butterworth, applied only before Hjorth computation)
  eog_hjorth_lo      = param.has( "eog-hjorth-lo"  ) ? param.requires_dbl( "eog-hjorth-lo"  ) : 0.5 ;
  eog_hjorth_hi      = param.has( "eog-hjorth-hi"  ) ? param.requires_dbl( "eog-hjorth-hi"  ) : 20.0 ;
  eog_hjorth_order   = param.has( "eog-hjorth-ord" ) ? param.requires_int( "eog-hjorth-ord" ) : 2 ;

  // Hjorth outlier
  eog_hjorth_z_th    = param.has( "eog-hjorth-z"  ) ? param.requires_dbl( "eog-hjorth-z"  ) : 10.0  ;

  // channel-level BAD rules
  eog_flag_prop = param.has( "eog-flag-prop"  ) ? param.requires_dbl( "eog-flag-prop"  ) : 0.5  ;
  eog_flag_run    = param.has( "eog-flag-run"   ) ? param.requires_dbl( "eog-flag-run"   ) : 3600.0 ;

  //
  // Generic options
  //

  // verbose epoch-level output

  by_epoch = param.has( "epoch" );
  
  
  
  //
  // Annotations & channels
  //

  resp_add_channel = param.has( "resp-add-channel" );
  if ( resp_add_channel ) resp_channel_label = param.value( "resp-add-channel" );

  //
  // Unified annotation settings (all domains)
  //
  //  Default: by-channel, no domain suffix  --> QC_<CH>,  QC_LN_<CH>
  //  annot=X        --> by-channel, prefix X  (empty value → "QC")
  //  annot-domain=X --> by-domain,  prefix X  (empty value → "QC"); no channel suffix
  //  annot-show=F   --> suppress all annotations
  //

  if ( param.has( "annot-domain" ) )
    {
      annot_by_channel  = false;
      annot_show_domain = true;
      annot_prefix = param.empty( "annot-domain" ) ? "QC" : param.value( "annot-domain" );
    }
  else
    {
      annot_by_channel  = true;
      annot_show_domain = false;
      if ( param.has( "annot" ) )
	annot_prefix = param.empty( "annot" ) ? "QC" : param.value( "annot" );
      else
	annot_prefix = "QC";
    }

  add_annots = ! ( param.has( "annot-show" ) && param.value( "annot-show" ) == "F" );

  
  //
  // Respiratory signals
  //

  do_resp( resp_signals );

  //
  // SpO2
  //

  do_spo2( spo2_signals );

  //
  // EEG
  //

  do_eeg( eeg_signals );

  //
  // EMG
  //

  do_emg( emg_signals );

  //
  // ECG
  //

  do_ecg( ecg_signals );

  //
  // EOG
  //

  do_eog( eog_signals );

}





void dsptools::qc_t::do_resp( signal_list_t & signals )
{

  const int ns = signals.size();

  if ( ns == 0 ) return;

  logger << "  ----------------------------------------\n"
	 << "  checking " << ns << " respiratory channel(s):\n"
	 << "     window size (sec)          = " << resp_window_dur    << "  (resp-win)\n"
	 << "     window step (sec)          = " << resp_window_inc    << "  (resp-inc)\n"
	 << "     signal range (Hz)          = " << resp_p1_lwr << " - " << resp_p1_upr
	                                                                  << "  (resp-p1-lwr, resp-p1-upr)\n"
	 << "     noise range (Hz)           = " << resp_p2_lwr << " - " << resp_p2_upr
	                                                                  << "  (resp-p2-lwr, resp-p2-upr)\n"
	 << "     SNR threshold              = " << resp_th            << "  (resp-snr-th)\n"
	 << "     minimum SR (Hz)            = " << resp_min_sr        << "  (resp-min-sr)\n"
	 << "     flatline floor             = " << resp_flatline_floor << "  (resp-flat-floor)\n"
	 << "     flatline prop. threshold   = " << resp_flatline_prop << "  (resp-flat-prop)\n"
	 << "     clip prop. threshold       = " << resp_clip_prop     << "  (resp-clip-prop)\n"
	 << "     jump p95-p05 multiplier    = " << resp_jump_mad_th   << "  (resp-jump-th)\n"
	 << "     bad epoch proportion       = " << resp_flag_prop << "  (resp-flag-prop)\n"
	 << "     bad contiguous run (sec)   = " << resp_flag_run   << "  (resp-flag-run)\n"
	 << "  outputs:\n"
	 << "     add noise waveform channel = " << ( resp_add_channel ? resp_channel_label : "no" )
	                                                                  << "  (resp-add-channel)\n"
	 << "     annotations                = " << ( add_annots
	                                               ? annot_prefix
	                                                 + ( annot_show_domain ? "_RESP" : "" )
	                                                 + ( annot_by_channel  ? "_<CH>" : "" )
	                                               : "none" )
	                                                                  << "  (annot=X, annot-domain=X, annot-show=F)\n";

  logger << "  changing EPOCH settings (including 'splice-gaps' mode)...\n";
  edf.timeline.set_epochs_to_span_gaps( true );
  edf.timeline.set_epoch( resp_window_dur , resp_window_inc );
  

  //
  // iterate over all signals       
  //
  
  for (int s=0; s<ns; s++)
    {

      const double Fs = edf.header.sampling_freq( signals(s) );

      //
      // Hard floor: below this SR no meaningful RESP analysis is possible;
      // skip the channel and flag it in the output.
      //

      if ( Fs < resp_hard_min_sr )
	{
	  logger << "  ** skipping " << signals.label(s)
		 << " (SR=" << Fs << " Hz is below resp-hard-min-sr=" << resp_hard_min_sr
		 << "; channel omitted from RESP analysis)\n";
	  writer.level( signals.label(s) , globals::signal_strat );
	  writer.level( 1 , "RESP" );
	  writer.value( "FLAGGED" , 1 );
	  writer.value( "LOW_SR"  , 1 );
	  writer.unlevel( "RESP" );
	  writer.level( "RESP" , "DOMAIN" );
	  writer.value( "FLAGGED" , 1 );
	  writer.value( "LOW_SR"  , 1 );
	  writer.unlevel( "DOMAIN" );
	  continue;
	}

      //
      // Soft minimum: if SR is below resp_min_sr but above the hard floor,
      // upsample to resp_min_sr so that the full noise band (resp_p2_upr)
      // and analysis window sizes are fully supported.  Note: upsampling does
      // not recover spectral information above the original Nyquist (Fs/2),
      // so the noise-band estimate will be naturally truncated for low-SR
      // channels; the check remains meaningful for the respiratory signal
      // band and flatline/clip/jump criteria.
      //

      const bool   do_upsample  = ( Fs < resp_min_sr );
      const double effective_Fs = do_upsample ? (double)resp_min_sr : Fs;

      if ( do_upsample )
	logger << "  ** upsampling " << signals.label(s)
	       << " from " << Fs << " to " << resp_min_sr
	       << " Hz for RESP analysis (original Nyquist = " << Fs/2.0 << " Hz)\n";

      //
      // start processing this signal
      //
      
      writer.level( signals.label(s) , globals::signal_strat );

      writer.level( 1 , "RESP" );  // CH x RESP table (RESP==1 for all rows)

      // for each (overlapping) epoch, track SNR as well as the start/stop sample-points
      // (i.e. use the latter to remake the original signal)

      std::vector<double> snr;
      std::vector<double> P1;
      std::vector<bool> valid;
      std::vector<std::pair<int,int> > smps;

      // per-window artifact flags
      std::vector<bool> flag_flat;
      std::vector<bool> flag_clip;
      std::vector<bool> flag_jump;
      std::vector<bool> flag_bad_epoch;
      std::vector<interval_t> epoch_intervals;

      const int ne = edf.timeline.first_epoch();

      if ( s == 0 )
	logger << "  considering " << ne << " windows per channel\n";

      while ( 1 )
	{

	  int epoch = edf.timeline.next_epoch();

	  if ( epoch == -1 ) break;

	  interval_t interval = edf.timeline.epoch( epoch );
	  epoch_intervals.push_back( interval );

	  // note: final true returns sample points
	  slice_t slice( edf , signals(s) , interval, 1, false , true );

	  std::vector<double> * d = slice.nonconst_pdata();

	  const std::vector<int> * sp = slice.psmps();

	  // number of samples at the original SR (needed for sp indexing below)
	  const int n_orig = (int)d->size();

	  // upsample epoch to effective_Fs if channel SR is below resp_min_sr
	  std::vector<double> d_up;
	  if ( do_upsample )
	    {
	      d_up = dsptools::resample( d , Fs , effective_Fs );
	      d = &d_up;
	    }

	  const int n = (int)d->size();

	  //
	  // Per-window artifact checks (must run before FFT modifies d)
	  //

	  // 1. Flatline / dropout: proportion of flat (near-constant) samples
	  const double flat_p = MiscMath::flat( *d , resp_flatline_floor );
	  const bool flag_flat_e = (flat_p >= resp_flatline_prop);

	  // 2. Clipping / saturation: proportion of samples at ADC limits
	  const double clip_p = MiscMath::clipped( *d );
	  const bool flag_clip_e = (clip_p >= resp_clip_prop);

	  // 3. Step / jump artifact: any consecutive diff > jump_mad_th * epoch p95-p05 range
	  //    p95-p05 range is used instead of MAD/IQR to avoid false positives on
	  //    non-stationary resp/effort signals (e.g. nasal pressure flat during
	  //    apnea then recovering): MAD and IQR are suppressed when the flat/apneic
	  //    fraction exceeds 50% or 75% respectively, making the threshold
	  //    hypersensitive to normal physiological transitions.  The p95-p05 range
	  //    remains non-zero as long as at least ~5% of the epoch contains signal
	  //    (i.e. ~6s of breathing in a 120s window).
	  bool flag_jump_e = false;
	  double epoch_range = 0;
	  {
	    std::vector<double> tmp = *d;
	    const double p05 = MiscMath::percentile( tmp , 0.05 );
	    std::vector<double> tmp2 = *d;
	    const double p95 = MiscMath::percentile( tmp2 , 0.95 );
	    epoch_range = p95 - p05;
	    if ( epoch_range > 0 )
	      for (int i = 1; i < n; i++)
		if ( fabs( (*d)[i] - (*d)[i-1] ) > resp_jump_mad_th * epoch_range )
		  { flag_jump_e = true; break; }
	  }

	  const bool epoch_bad = flag_flat_e || flag_clip_e || flag_jump_e;

	  flag_flat.push_back( flag_flat_e );
	  flag_clip.push_back( flag_clip_e );
	  flag_jump.push_back( flag_jump_e );
	  flag_bad_epoch.push_back( epoch_bad );

	  FFT fftseg( n , n , effective_Fs , FFT_FORWARD , WINDOW_NONE );

	  fftseg.apply( &((*d)[0]) , n );
	  
	  double p1 = 0 , p2 = 0;

	  for (int f=0; f<fftseg.cutoff; f++)
	    {
	      if ( fftseg.frq[f] >= resp_p1_lwr && fftseg.frq[f] < resp_p1_upr )
		{ p1 += fftseg.X[f];  } 
	      else if ( fftseg.frq[f] >= resp_p2_lwr && fftseg.frq[f] < resp_p2_upr )
		{ p2 += fftseg.X[f];  } 
	      else if ( fftseg.frq[f] > resp_p2_upr ) break;
	    }
	  	  
	  // catch degenerate cases	  
	  const bool valid1 = ! ( fabs( p1 ) < resp_epsilon || fabs( p2 ) < resp_epsilon  );
	  
	  // get SNR (?how to handle degenerate cases)
	  const double snr1 = p1 / p2 ; 
	  
	  // sample span of this window (for reconstructing an original signal)
	  //  although note that epochs are not contiguous in sample-space, and
	  //  so below, we only extract non-missing points from the start/stop
	  //  ranges stored in smps[]

	  // sp is indexed at original SR; use n_orig (not n) when upsampled
	  std::pair<int,int> smps1( (*sp)[0] , (*sp)[n_orig-1] );
	  
	  // track for this epoch
	  snr.push_back( snr1 );
	  P1.push_back( p1 );
	  valid.push_back( valid1 );
	  smps.push_back( smps1 );

	} // next epoch
      

      //
      // 90th-percentile SNR must be > 10x to consider the signal
      //
      
      const double snr90 = MiscMath::percentile( snr , 0.9 );
      
      bool valid_signal = snr90 > resp_th ;  // defaults to 10x
      
      
      //
      // Get P1 scaling factor (based on log-mean of over valid regions)
      //

      // calculate:  P1scale = 10^mean(log10(P1(StoN_x>minStoN)))
      //   i.e. mean of P1 for SNR-positive regions (using log-transform for robustness)

      int n1 = 0;
      double sum = 0;
      
      const int total_epochs = snr.size();
      int valid_epochs = 0;
      double snr1 = 0;

      for (int e=0; e<total_epochs; e++)
	{
	  if ( valid[e] )
	    {
	      ++valid_epochs;
	      if ( snr[e] > resp_th )
		{
		  sum += log10( P1[e] );
		  ++n1;
		}
	    }
	}
      
      const double P1scale = n1 > 0 ? pow( 10.0 , sum / (double)n1 ) : 1 ; 

      //
      // Make criteria (0 for invalid regions, will be flagged below) 
      //
      
      std::vector<double> criteria( total_epochs );

      for (int e=0; e<total_epochs; e++)
	criteria[e] = valid[e] ? ( P1[e] / P1scale ) * snr[e] : 0 ;
      
      
      //
      // Verbose outputs
      //

      if ( by_epoch )
	{
	  for (int e=0; e<total_epochs; e++)
	    {
	      writer.level( e+1 , "WIN" ) ;
	      writer.value( "P1" , P1[e] );
	      writer.value( "SNR" , snr[e] );
	      writer.value( "CRIT" , criteria[e] );
	      writer.value( "FLAT" , (int)flag_flat[e] );
	      writer.value( "CLIP" , (int)flag_clip[e] );
	      writer.value( "JUMP" , (int)flag_jump[e] );
	      writer.value( "FLAG_EPOCH" , (int)flag_bad_epoch[e] );
	    }
	  writer.unlevel( "WIN" );
	}
            
      //
      // Derive % of signal that is noise
      //

      // noisewav below is created as a continuous signal and is
      // assigned values 1, 2 or 3. 

      // 1 : criteria < 1/4    = 0.25
      // 2 : criteria < 0.5/4  = 0.125 
      // 3 : criteria < 0.24/4 = 0.0625 

      // each sample spanned by multiple windows: assign the worst (3 worse than 2 worse than 1)

      
      // pull original signal for entire trace just to get N and time-points
      // but also need to pull the samples here (as above the sample-points saved
      // were w.r.t to original EDF, not any masked version)

      // as we're using gap-spliced epoching, some points in the start/stop interval
      // saved (in 'smps[]') may be gaps;  we can just ignore them here, as they
      // won't be in sp2slot

      // --> any un-epoched portions (e.g. <10 second parts at the end of the record)
      //     won't be assigned any noise values (defaults to 0) but that's probably
      //     fine - i.e. should be a smaller value
      
      // n.b. extra terms to pull psmps() [ last true ] 
      slice_t slice( edf , signals(s) , edf.timeline.wholetrace(),1,false,true);
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      const int n = tp->size();
      
      // make mapping of sample-points to slots
      const std::vector<int> * sp = slice.psmps();
      std::map<int,int> sp2slot;
      for (int i=0; i<n; i++)
	sp2slot[ (*sp)[i] ] = i;
      
      // new noise waveform (initialize at 0) 
      std::vector<double> noisewav( n , 0 );
      
      // iterate over each window:
      std::vector<double> fac = { 1 / 4.0 , 0.5 / 4.0 , 0.25 / 4.0 };
      
      for (int f=0; f<3; f++) // <1, <2, <3	
	{
	  // set for this window
	  for (int e=0; e<total_epochs; e++)
	    {	      
	      if ( criteria[e] < fac[f] )
		{		  
		  const int start = smps[e].first;
		  const int stop  = smps[e].second;
		  
		  for (int p=start; p<=stop; p++)		    
		    {
		      std::map<int,int>::const_iterator ii = sp2slot.find( p );
		      
                      if ( ii != sp2slot.end() )
                        {
			  const int slot = ii->second;		      
			  noisewav[ slot ] = f+1;		
			}
		    }
		}
	    }
	}

      
      //
      // Estimate fraction of recording that is noise (nb: allow for floating point
      // issues, though this should be fine) 
      //

      double f1 = 0 , f2 = 0 , f3 = 0;
      for (int p=0; p<n; p++)
	{
	  if ( noisewav[p] >= .99 )
	    {
	      ++f1;
	      if ( noisewav[p] >= 1.99 )
		{
		  ++f2;
		  if ( noisewav[p] >= 2.99 )
		    {
		      ++f3;
		    }
		}
	    }
	}
      
      // as a proportion of all samples
      f1 /= (double)n;
      f2 /= (double)n;
      f3 /= (double)n;
      
      //
      // Score 30s epochs: a 30s epoch is bad if bad 120s windows
      // (flat/clip/jump OR poor spectral SNR) cover >= half of it.
      // FLAGGED, N_FLAG_EPOCH and MAX_FLAG_RUN are derived from these counts.
      //

      std::vector<interval_t> bad_win_intervals;
      for (int e = 0; e < total_epochs; e++)
	if ( flag_bad_epoch[e] || criteria[e] < fac[1] )
	  bad_win_intervals.push_back( epoch_intervals[e] );

      int n_bad_30 = 0, total_30 = 0;
      int max_consec_30 = 0, cur_consec_30 = 0;
      const uint64_t half_epoch_tp = 15LLU * globals::tp_1sec;
      std::vector<interval_t> bad_30_epochs;

      edf.timeline.set_epochs_to_span_gaps( false );
      edf.timeline.set_epoch( 30.0 , 30.0 );
      edf.timeline.first_epoch();
      while ( 1 )
	{
	  int epoch30 = edf.timeline.next_epoch();
	  if ( epoch30 == -1 ) break;
	  ++total_30;
	  interval_t iv = edf.timeline.epoch( epoch30 );

	  uint64_t covered = 0LLU;
	  uint64_t union_start = 0LLU , union_stop = 0LLU;
	  bool in_union = false;
	  for ( const auto & bw : bad_win_intervals )
	    {
	      if ( bw.start >= iv.stop  ) break;
	      if ( bw.stop  <= iv.start ) continue;
	      uint64_t cs = bw.start > iv.start ? bw.start : iv.start;
	      uint64_t ce = bw.stop  < iv.stop  ? bw.stop  : iv.stop;
	      if ( !in_union )
		{ union_start = cs; union_stop = ce; in_union = true; }
	      else if ( cs <= union_stop )
		{ if ( ce > union_stop ) union_stop = ce; }
	      else
		{ covered += union_stop - union_start; union_start = cs; union_stop = ce; }
	    }
	  if ( in_union ) covered += union_stop - union_start;

	  if ( covered >= half_epoch_tp )
	    {
	      ++n_bad_30;
	      if ( ++cur_consec_30 > max_consec_30 ) max_consec_30 = cur_consec_30;
	      bad_30_epochs.push_back( iv );
	    }
	  else
	    cur_consec_30 = 0;
	}

      edf.timeline.set_epochs_to_span_gaps( true );
      edf.timeline.set_epoch( resp_window_dur , resp_window_inc );

      const double p_bad_30 = total_30 > 0 ? (double)n_bad_30 / total_30 : 0;
      const double max_bad_run_30_sec = max_consec_30 * 30.0;

      const bool bad = ( p_bad_30          >= resp_flag_prop ) ||
	               ( max_bad_run_30_sec >= resp_flag_run  );


      //
      // RESP-specific outputs  (CH x RESP table)
      //

      writer.value( "N_VALID_WIN" , valid_epochs );
      writer.value( "P_VALID_WIN" , valid_epochs /(double)total_epochs );
      writer.value( "P_NOISE1" , f1 );
      writer.value( "P_NOISE2" , f2 );
      writer.value( "P_NOISE3" , f3 );

      // per-reason epoch proportions (RESP-specific)
      int nf_r = 0 , nc_r = 0 , nj_r = 0;
      for (int e=0; e<total_epochs; e++)
	{
	  if ( flag_flat[e]  ) ++nf_r;
	  if ( flag_clip[e]  ) ++nc_r;
	  if ( flag_jump[e]  ) ++nj_r;
	}
      writer.value( "N_EPOCH"     , total_30 );
      writer.value( "P_FLAG_EPOCH" , p_bad_30 );
      writer.value( "P_FLAT"      , total_epochs > 0 ? (double)nf_r / total_epochs : 0 );
      writer.value( "P_CLIP"      , total_epochs > 0 ? (double)nc_r / total_epochs : 0 );
      writer.value( "P_JUMP"      , total_epochs > 0 ? (double)nj_r / total_epochs : 0 );

      writer.unlevel( "RESP" );

      //
      // Cross-domain summary outputs  (CH x DOMAIN table)
      //

      writer.level( "RESP" , "DOMAIN" );
      writer.value( "FLAGGED"      , (int)(bad) );
      writer.value( "N_FLAG_EPOCH" , n_bad_30 );
      writer.value( "MAX_FLAG_RUN" , max_bad_run_30_sec );
      writer.unlevel( "DOMAIN" );

      
      //
      // add as a channel?
      //
      
      if ( resp_add_channel )
	{
	  const std::string lab = resp_channel_label + "_" + signals.label(s) ;
	  edf.add_signal( lab , Fs , noisewav );
	  logger << "  adding new QC signal " << lab << ", " << Fs << " Hz\n"; 
	}

      //
      // Annotations: bad epochs (merged, gap-aware)
      //

      if ( add_annots )
	{
	  std::string alab = annot_prefix;
	  if ( annot_show_domain )  alab += "_RESP";
	  if ( annot_by_channel )   alab += "_" + signals.label(s);

	  annot_t * a = edf.annotations->add( alab );

	  // Emit merged runs of consecutive bad 30s epochs as annotations
	  bool in_run = false;
	  uint64_t run_start = 0LLU , run_stop = 0LLU;
	  for ( const auto & iv : bad_30_epochs )
	    {
	      if ( !in_run )
		{ run_start = iv.start; run_stop = iv.stop; in_run = true; }
	      else if ( iv.start <= run_stop )
		{ if ( iv.stop > run_stop ) run_stop = iv.stop; }
	      else
		{
		  a->add( "." , interval_t( run_start , run_stop ) , signals.label(s) );
		  run_start = iv.start; run_stop = iv.stop;
		}
	    }
	  if ( in_run )
	    a->add( "." , interval_t( run_start , run_stop ) , signals.label(s) );

	  logger << "  " << signals.label(s) << " : " << n_bad_30
		 << " bad 30s epochs flagged (" << (int)bad_win_intervals.size()
		 << " bad 120s windows)\n";
	}
      
    } // next signal

  writer.unlevel( globals::signal_strat );

}



void dsptools::qc_t::do_eeg( signal_list_t & signals )
{

  const int ns = signals.size();

  if ( ns == 0 ) return;

  logger << "  ----------------------------------------\n"
	 << "  checking " << ns << " EEG channel(s):\n"
	 << "     window size (sec)           = " << eeg_window_dur     << "  (eeg-win)\n"
	 << "     window step (sec)           = " << eeg_window_inc     << "  (eeg-inc)\n"
	 << "     minimum SR (Hz)             = " << eeg_min_sr         << "  (eeg-min-sr)\n"
	 << "     flatline SD threshold (µV)  = " << eeg_flat_std_th    << "  (eeg-flat-th)\n"
	 << "     flatline deriv. proportion  = " << eeg_flat_deriv_prop << "  (eeg-flat-prop)\n"
	 << "     flatline deriv. epsilon     = " << eeg_flat_deriv_eps  << "  (eeg-flat-eps)\n"
	 << "     clip proportion             = " << eeg_clip_prop      << "  (eeg-clip-prop)\n"
	 << "     amplitude threshold (µV)    = " << eeg_amp_th         << "  (eeg-amp-th)\n"
	 << "     amplitude proportion        = " << eeg_amp_prop       << "  (eeg-amp-prop)\n"
	 << "     HF band (Hz)               = " << eeg_hf_lo << " - " << eeg_hf_hi
                                                                          << "  (eeg-hf-lo, eeg-hf-hi)\n"
	 << "     signal band (Hz)            = " << eeg_sig_lo << " - " << eeg_sig_hi
                                                                          << "  (eeg-sig-lo, eeg-sig-hi)\n"
	 << "     HF/signal ratio threshold   = " << eeg_hf_ratio_th    << "  (eeg-hf-th)\n"
	 << "     line noise ± bw (Hz)        = " << eeg_ln_bw          << "  (eeg-ln-bw)\n"
	 << "     line noise broad range (Hz) = " << eeg_ln_lo << " - " << eeg_ln_hi
                                                                          << "  (eeg-ln-lo, eeg-ln-hi)\n"
	 << "     line noise ratio threshold  = " << eeg_ln_ratio_th    << "  (eeg-ln-th)\n"
	 << "     spectral peakedness (SPK)   = " << eeg_spectral_peakedness_th << "  (eeg-peak-th)\n"
	 << "     spectral kurtosis (KURT)    = " << eeg_spectral_skewness_th   << "  (eeg-kurt-th)\n"
	 << "     Hjorth pre-filter (Hz)      = " << eeg_hjorth_lo << " - " << eeg_hjorth_hi
                                                                          << "  (eeg-hjorth-lo, eeg-hjorth-hi)\n"
	 << "     Hjorth filter order         = " << eeg_hjorth_order   << "  (eeg-hjorth-ord)\n"
	 << "     Hjorth outlier |z| thresh.  = " << eeg_hjorth_z_th    << "  (eeg-hjorth-z)\n"
	 << "     channel bad epoch prop.     = " << eeg_flag_prop << "  (eeg-flag-prop)\n"
	 << "     channel bad run (sec)       = " << eeg_flag_run    << "  (eeg-flag-run)\n";

  // PWELCH parameters (4s segments, 50% overlap)
  const double seg_sec  = 4.0;
  const double step_sec = 2.0;

  // peakedness frequency range (output only; thresholds TBD)
  const double peak_minf          = 2.0;
  const double peak_maxf          = 28.0;
  const int    peak_median_filt_n = 11;

  edf.timeline.set_epochs_to_span_gaps( false );
  edf.timeline.set_epoch( eeg_window_dur , eeg_window_inc );

  for ( int s = 0; s < ns; s++ )
    {

      const double Fs      = edf.header.sampling_freq( signals(s) );
      const double nyquist = Fs / 2.0;

      if ( Fs < eeg_min_sr )
	{
	  logger << "  ** skipping " << signals.label(s)
		 << " (SR=" << Fs << " < eeg-min-sr=" << eeg_min_sr << ")\n";
	  continue;
	}

      // Ensure µV scale (auto-converts mV/V; warns if unrecognised)
      ensure_units( signals(s) , "uV" );

      // HF ratio check requires HF upper bound below Nyquist
      const bool do_hf_ratio = ( nyquist > eeg_hf_hi );
      if ( ! do_hf_ratio )
	logger << "  ** note: " << signals.label(s)
	       << " Nyquist (" << nyquist << " Hz) <= eeg-hf-hi (" << eeg_hf_hi
	       << " Hz); HF ratio check skipped\n";

      // IIR Butterworth bandpass for Hjorth pre-filtering (raw signal used for all other checks)
      const double hjorth_lo_safe = eeg_hjorth_lo < nyquist ? eeg_hjorth_lo : nyquist * 0.1;
      const double hjorth_hi_safe = eeg_hjorth_hi < nyquist ? eeg_hjorth_hi : nyquist * 0.9;
      iir_t hjorth_iir;
      bool  hjorth_filter_ok = false;
      if ( hjorth_lo_safe < hjorth_hi_safe )
	{
	  hjorth_iir.init( BUTTERWORTH_BANDPASS , eeg_hjorth_order , Fs , hjorth_lo_safe , hjorth_hi_safe );
	  hjorth_filter_ok = true;
	}
      else
	logger << "  ** note: " << signals.label(s)
	       << " could not set up Hjorth bandpass filter (Fs=" << Fs << "); Hjorth computed on raw signal\n";

      writer.level( signals.label(s) , globals::signal_strat );
      writer.level( 1 , "EEG" );

      const int ne = edf.timeline.first_epoch();
      if ( s == 0 ) logger << "  considering " << ne << " epochs per channel\n";

      // per-epoch metric storage (two-pass: Hjorth outlier needs full distribution)
      std::vector<bool>       flag_flat, flag_clip, flag_amp, flag_hf, flag_ln, flag_peak;
      std::vector<bool>       flag_hjorth;     // filled in second pass
      std::vector<double>     ep_sd, ep_deriv;
      std::vector<double>     ep_hf_ratio, ep_ln_ratio;
      std::vector<double>     ep_spk, ep_kurt; // spectral peakedness (output only)
      std::vector<double>     ep_act, ep_cplx; // Hjorth activity, complexity
      std::vector<interval_t> epoch_intervals;

      //
      // First pass: compute all per-epoch metrics
      //

      while ( 1 )
	{
	  int epoch = edf.timeline.next_epoch();
	  if ( epoch == -1 ) break;

	  interval_t interval = edf.timeline.epoch( epoch );
	  epoch_intervals.push_back( interval );

	  slice_t slice( edf , signals(s) , interval , 1 , false , false );
	  const std::vector<double> d = *slice.nonconst_pdata();
	  const int n = (int)d.size();

	  //
	  // 1. Flatline: SD < threshold OR near-zero derivative prop >= threshold
	  //

	  const double sd    = MiscMath::sdev( d );
	  const double deriv = MiscMath::flat( d , eeg_flat_deriv_eps );
	  const bool f_flat  = ( sd < eeg_flat_std_th ) || ( deriv >= eeg_flat_deriv_prop );

	  //
	  // 2. Clipping: >= eeg_clip_prop of samples at ADC limits
	  //

	  const bool f_clip = ( MiscMath::clipped( d ) >= eeg_clip_prop );

	  //
	  // 3. Extreme amplitude: >= eeg_amp_prop samples with |EEG| > eeg_amp_th
	  //

	  int n_amp = 0;
	  for ( int i = 0; i < n; i++ )
	    if ( std::fabs( d[i] ) > eeg_amp_th ) ++n_amp;
	  const bool f_amp = ( n > 0 ) && ( (double)n_amp / n >= eeg_amp_prop );

	  //
	  // 4. Spectral: HF contamination + line noise + peakedness (one PWELCH call)
	  //

	  double hf_ratio = -1.0;
	  double ln_ratio = -1.0;
	  double spk      = 0.0;
	  double kurt_s   = 0.0;
	  bool   f_hf     = false; // HF ratio only (stays in BAD)
	  bool   f_ln     = false; // line noise (tracked separately)

	  {
	    const int seg_pts  = (int)( seg_sec  * Fs );
	    const int ovlp_pts = (int)( step_sec * Fs );

	    if ( seg_pts <= n && seg_pts > ovlp_pts )
	      {
		const int nsegs = (int)floor( ( n - ovlp_pts )
					     / (double)( seg_pts - ovlp_pts ) );
		if ( nsegs >= 1 )
		  {
		    PWELCH pw( d , Fs , seg_sec , nsegs , WINDOW_TUKEY50 );

		    // 4a. HF ratio: P[20-40 Hz] / P[0.5-20 Hz]
		    if ( do_hf_ratio )
		      {
			const double p_hf  = pw.psdsum( eeg_hf_lo  , eeg_hf_hi  );
			const double p_sig = pw.psdsum( eeg_sig_lo , eeg_sig_hi );
			if ( p_sig > 0.0 )
			  {
			    hf_ratio = p_hf / p_sig;
			    if ( hf_ratio > eeg_hf_ratio_th ) f_hf = true;
			  }
		      }

		    // 4b. Line noise: max(P50, P60) / P[0.5-40 Hz]
		    const double broad_hi = eeg_ln_hi < nyquist - 0.5 ? eeg_ln_hi : nyquist - 0.5;
		    if ( broad_hi > eeg_ln_lo )
		      {
			const double p_broad = pw.psdsum( eeg_ln_lo , broad_hi );
			if ( p_broad > 0.0 )
			  {
			    const double f50_lo = 50.0 - eeg_ln_bw;
			    const double f50_hi = 50.0 + eeg_ln_bw;
			    const double f60_lo = 60.0 - eeg_ln_bw;
			    const double f60_hi = 60.0 + eeg_ln_bw;
			    const double p_50 = ( f50_hi < nyquist ) ? pw.psdsum( f50_lo , f50_hi ) : 0.0;
			    const double p_60 = ( f60_hi < nyquist ) ? pw.psdsum( f60_lo , f60_hi ) : 0.0;
			    const double p_line = p_50 > p_60 ? p_50 : p_60;
			    if ( p_50 > 0.0 || p_60 > 0.0 )
			      {
				ln_ratio = p_line / p_broad;
				if ( ln_ratio > eeg_ln_ratio_th ) f_ln = true;
			      }
			  }
		      }

		    // 4c. Spectral peakedness (output only, thresholds TBD)
		    {
		      std::vector<double> frq, logged;
		      for ( int i = 0; i < (int)pw.freq.size(); i++ )
			{
			  if ( pw.freq[i] < peak_minf ) continue;
			  if ( pw.freq[i] > peak_maxf ) break;
			  const double pv = pw.psd[i] > 0.0 ? pw.psd[i] : 1e-4;
			  frq.push_back( pw.freq[i] );
			  logged.push_back( 10.0 * log10( pv ) );
			}
		      if ( (int)frq.size() >= (int)( peak_median_filt_n * 1.5 ) )
			psd_shape_metrics( frq , logged , peak_median_filt_n ,
					   &spk , &kurt_s , NULL , NULL , NULL );
		    }
		  }
	      }
	  }

	  //
	  // 5. Hjorth parameters (activity, complexity — store for outlier pass)
	  //    Pre-filter with IIR Butterworth; raw signal d unchanged for all other checks.
	  //

	  double activity = 0 , mobility = 0 , complexity = 0;
	  if ( hjorth_filter_ok )
	    {
	      const std::vector<double> df = hjorth_iir.apply( d );
	      MiscMath::hjorth( &df , &activity , &mobility , &complexity , ! globals::legacy_hjorth );
	    }
	  else
	    {
	      MiscMath::hjorth( &d , &activity , &mobility , &complexity , ! globals::legacy_hjorth );
	    }

	  const bool f_peak = ( ( eeg_spectral_peakedness_th > 0.0 && spk    > eeg_spectral_peakedness_th ) ||
				( eeg_spectral_skewness_th   > 0.0 && kurt_s > eeg_spectral_skewness_th   ) );

	  flag_flat.push_back( f_flat );
	  flag_clip.push_back( f_clip );
	  flag_amp.push_back( f_amp );
	  flag_hf.push_back( f_hf );
	  flag_ln.push_back( f_ln );
	  flag_peak.push_back( f_peak );
	  ep_sd.push_back( sd );
	  ep_deriv.push_back( deriv );
	  ep_hf_ratio.push_back( hf_ratio );
	  ep_ln_ratio.push_back( ln_ratio );
	  ep_spk.push_back( spk );
	  ep_kurt.push_back( kurt_s );
	  ep_act.push_back( activity );
	  ep_cplx.push_back( complexity );

	} // next epoch

      //
      // Second pass: Hjorth robust z-score outlier detection
      // flag if |robust-z| > eeg_hjorth_z_th for activity OR complexity
      //

      const int total_epochs = (int)ep_act.size();
      flag_hjorth.assign( total_epochs , false );

      if ( total_epochs > 2 )
	{
	  auto flag_by_mad = [&]( const std::vector<double> & vals )
	  {
	    std::vector<double> tmp = vals;
	    const double med = MiscMath::median( tmp );
	    std::vector<double> abs_dev( total_epochs );
	    for ( int e = 0; e < total_epochs; e++ )
	      abs_dev[e] = std::fabs( vals[e] - med );
	    std::vector<double> tmp2 = abs_dev;
	    const double mad   = MiscMath::median( tmp2 );
	    const double scale = mad > 0.0 ? 1.4826 * mad : 1.0;
	    for ( int e = 0; e < total_epochs; e++ )
	      if ( std::fabs( vals[e] - med ) / scale > eeg_hjorth_z_th )
		flag_hjorth[e] = true;
	  };

	  flag_by_mad( ep_act  );
	  flag_by_mad( ep_cplx );
	}

      //
      // Combine flags
      //

      std::vector<bool> flag_bad( total_epochs );
      for ( int e = 0; e < total_epochs; e++ )
	flag_bad[e] = flag_flat[e] || flag_clip[e] || flag_amp[e]
	           || flag_hf[e]   || flag_peak[e]  || flag_hjorth[e];

      //
      // Epoch-level output (deferred so Hjorth flags are ready)
      //

      if ( by_epoch )
	{
	  for ( int e = 0; e < total_epochs; e++ )
	    {
	      writer.level( e + 1 , "WIN" );
	      writer.value( "SD"        , ep_sd[e]       );
	      writer.value( "DERIV"     , ep_deriv[e]    );
	      writer.value( "HF_RATIO"  , ep_hf_ratio[e] );
	      writer.value( "LN_RATIO"  , ep_ln_ratio[e] );
	      writer.value( "SPK"       , ep_spk[e]      );
	      writer.value( "KURT"      , ep_kurt[e]     );
	      writer.value( "ACT"       , ep_act[e]      );
	      writer.value( "CPLX"      , ep_cplx[e]     );
	      writer.value( "FLAT"      , (int)flag_flat[e]   );
	      writer.value( "CLIP"      , (int)flag_clip[e]   );
	      writer.value( "AMP"       , (int)flag_amp[e]    );
	      writer.value( "HF"        , (int)flag_hf[e]     );
	      writer.value( "LN"        , (int)flag_ln[e]     );
	      writer.value( "PEAK"      , (int)flag_peak[e]   );
	      writer.value( "HJORTH"    , (int)flag_hjorth[e] );
	      writer.value( "FLAG_EPOCH" , (int)flag_bad[e]    );
	    }
	  writer.unlevel( "WIN" );
	}

      //
      // Channel-level summary
      //

      int n_bad = 0 , n_flat = 0 , n_clip = 0 , n_amp = 0 , n_hf = 0 , n_ln = 0 , n_peak = 0 , n_hjorth = 0;
      int max_consec = 0 , cur_consec = 0;
      int max_ln_consec = 0 , cur_ln_consec = 0;

      for ( int e = 0; e < total_epochs; e++ )
	{
	  if ( flag_flat[e]   ) ++n_flat;
	  if ( flag_clip[e]   ) ++n_clip;
	  if ( flag_amp[e]    ) ++n_amp;
	  if ( flag_hf[e]     ) ++n_hf;
	  if ( flag_peak[e]   ) ++n_peak;
	  if ( flag_hjorth[e] ) ++n_hjorth;
	  if ( flag_ln[e] )
	    {
	      ++n_ln;
	      if ( ++cur_ln_consec > max_ln_consec ) max_ln_consec = cur_ln_consec;
	    }
	  else cur_ln_consec = 0;
	  if ( flag_bad[e] )
	    {
	      ++n_bad;
	      if ( ++cur_consec > max_consec ) max_consec = cur_consec;
	    }
	  else cur_consec = 0;
	}

      const double prop_bad    = total_epochs > 0 ? n_bad / (double)total_epochs : 0.0;
      const double max_bad_run = max_consec * eeg_window_inc;
      const bool   ch_bad      = ( prop_bad > eeg_flag_prop ) ||
	                          ( max_bad_run >= eeg_flag_run );

      const double prop_ln     = total_epochs > 0 ? n_ln / (double)total_epochs : 0.0;
      const double max_ln_run  = max_ln_consec * eeg_window_inc;
      const bool   ch_ln_bad   = ( prop_ln > eeg_flag_prop ) ||
	                          ( max_ln_run >= eeg_flag_run );

      writer.value( "N_FLAG_EPOCH" , n_bad      );
      writer.value( "PROP_FLAG"    , prop_bad   );
      writer.value( "MAX_FLAG_RUN" , max_bad_run );
      writer.value( "FLAGGED"         , (int)ch_bad );
      writer.value( "P_FLAT"      , total_epochs > 0 ? n_flat   / (double)total_epochs : 0.0 );
      writer.value( "P_CLIP"      , total_epochs > 0 ? n_clip   / (double)total_epochs : 0.0 );
      writer.value( "P_AMP"       , total_epochs > 0 ? n_amp    / (double)total_epochs : 0.0 );
      writer.value( "P_HF"        , total_epochs > 0 ? n_hf     / (double)total_epochs : 0.0 );
      writer.value( "P_PEAK"      , total_epochs > 0 ? n_peak   / (double)total_epochs : 0.0 );
      writer.value( "P_HJORTH"    , total_epochs > 0 ? n_hjorth / (double)total_epochs : 0.0 );
      writer.value( "P_LN"        , prop_ln    );
      writer.value( "MAX_LN_RUN"  , max_ln_run );
      writer.value( "LN_FLAG"      , (int)ch_ln_bad );

      writer.unlevel( "EEG" );

      // Cross-domain summary (CH x DOMAIN table)
      writer.level( "EEG" , "DOMAIN" );
      writer.value( "FLAGGED"         , (int)ch_bad    );
      writer.value( "N_FLAG_EPOCH" , n_bad          );
      writer.value( "MAX_FLAG_RUN" , max_bad_run    );
      writer.value( "LN_FLAG"      , (int)ch_ln_bad );
      writer.value( "MAX_LN_RUN"  , max_ln_run     );
      writer.unlevel( "DOMAIN" );

      //
      // Annotations (merged, gap-aware)
      //

      if ( add_annots )
	{
	  std::string alab = annot_prefix;
	  if ( annot_show_domain ) alab += "_EEG";
	  if ( annot_by_channel  ) alab += "_" + signals.label(s);

	  annot_t * a = edf.annotations->add( alab );

	  bool in_bad = false;
	  uint64_t bad_start = 0LLU , bad_stop = 0LLU;

	  for ( int e = 0; e < total_epochs; e++ )
	    {
	      if ( flag_bad[e] )
		{
		  if ( ! in_bad )
		    {
		      bad_start = epoch_intervals[e].start;
		      bad_stop  = epoch_intervals[e].stop;
		      in_bad = true;
		    }
		  else
		    {
		      if ( epoch_intervals[e].start > bad_stop )
			{
			  a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );
			  bad_start = epoch_intervals[e].start;
			}
		      if ( epoch_intervals[e].stop > bad_stop )
			bad_stop = epoch_intervals[e].stop;
		    }
		}
	      else if ( in_bad )
		{
		  a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );
		  in_bad = false;
		}
	    }

	  if ( in_bad )
	    a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );

	  // separate LN annotations
	  if ( n_ln > 0 )
	    {
	      std::string ln_lab = annot_prefix;
	      if ( annot_show_domain ) ln_lab += "_EEG";
	      ln_lab += "_LN";
	      if ( annot_by_channel ) ln_lab += "_" + signals.label(s);

	      annot_t * a_ln = edf.annotations->add( ln_lab );
	      bool in_ln = false;
	      uint64_t ln_start = 0LLU , ln_stop = 0LLU;

	      for ( int e = 0; e < total_epochs; e++ )
		{
		  if ( flag_ln[e] )
		    {
		      if ( ! in_ln ) { ln_start = epoch_intervals[e].start; ln_stop = epoch_intervals[e].stop; in_ln = true; }
		      else { if ( epoch_intervals[e].start > ln_stop ) { a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) ); ln_start = epoch_intervals[e].start; } if ( epoch_intervals[e].stop > ln_stop ) ln_stop = epoch_intervals[e].stop; }
		    }
		  else if ( in_ln ) { a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) ); in_ln = false; }
		}
	      if ( in_ln ) a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) );
	    }
	}

      logger << "  " << signals.label(s)
	     << " : " << n_bad << "/" << total_epochs << " flagged epochs"
	     << "  (flat=" << n_flat << " clip=" << n_clip << " amp=" << n_amp
	     << " hf=" << n_hf << " peak=" << n_peak << " hjorth=" << n_hjorth << ")"
	     << " | ln=" << n_ln
	     << ( ch_ln_bad ? " --> LN_FLAG" : "" )
	     << ( ch_bad    ? " --> FLAGGED channel" : "" ) << "\n";

    } // next channel

  writer.unlevel( globals::signal_strat );

}



void dsptools::qc_t::do_spo2( signal_list_t & signals )
{

  const int ns = signals.size();

  if ( ns == 0 ) return;

  logger << "  ----------------------------------------\n"
	 << "  checking " << ns << " SpO2/OXY channel(s):\n"
	 << "     window size (sec)          = " << spo2_window_dur     << "  (oxy-win)\n"
	 << "     window step (sec)          = " << spo2_window_inc     << "  (oxy-inc)\n"
	 << "     minimum SR (Hz)            = " << spo2_min_sr         << "  (oxy-min-sr)\n"
	 << "     valid range (%)            = " << spo2_range_lo << " - " << spo2_range_hi
	                                                                   << "  (oxy-range-lo, oxy-range-hi)\n"
	 << "     out-of-range prop.         = " << spo2_range_prop     << "  (oxy-range-prop)\n"
	 << "     flatline threshold (%)     = " << spo2_flatline_th    << "  (oxy-flat-th)\n"
	 << "     flatline run (epochs)      = " << spo2_flatline_k     << "  (oxy-flat-k)\n"
	 << "     jump threshold (%)         = " << spo2_jump_th        << "  (oxy-jump-th)\n"
	 << "     invalid floor (%)          = " << spo2_invalid_floor  << "  (oxy-invalid-floor)\n"
	 << "     missing prop.              = " << spo2_missing_prop   << "  (oxy-missing-prop)\n"
	 << "     bad epoch proportion       = " << spo2_flag_prop << "  (oxy-flag-prop)\n"
	 << "     bad contiguous run (sec)   = " << spo2_flag_run    << "  (oxy-flag-run)\n"
	 << "  scale detection: auto (0-1 → 0-100 if max <= 1.5)\n"
	 << "  outputs:\n"
	 << "     annotations                = " << ( add_annots
	                                               ? annot_prefix
	                                                 + ( annot_show_domain ? "_OXY" : "" )
	                                                 + ( annot_by_channel  ? "_<CH>" : "" )
	                                               : "none" )
	                                                                   << "  (annot=X, annot-domain=X, annot-show=F)\n";

  edf.timeline.set_epochs_to_span_gaps( false );
  edf.timeline.set_epoch( spo2_window_dur , spo2_window_inc );


  //
  // iterate over all signals
  //

  for (int s=0; s<ns; s++)
    {

      const double Fs = edf.header.sampling_freq( signals(s) );

      if ( Fs < spo2_min_sr )
	Helper::halt( signals.label(s) + " has a sample rate of "
		      + Helper::dbl2str( Fs ) + ", lower oxy-min-sr="
		      + Helper::dbl2str( spo2_min_sr ) );

      writer.level( signals.label(s) , globals::signal_strat );

      //
      // Scale detection: pull whole trace, check whether 0-1 or 0-100
      //

      slice_t wslice( edf , signals(s) , edf.timeline.wholetrace() , 1 , false , false );
      const std::vector<double> * wd = wslice.pdata();
      double sig_max = wd->empty() ? 0 : *std::max_element( wd->begin() , wd->end() );
      const bool rescaled = sig_max <= 1.5;

      if ( rescaled )
	logger << "  " << signals.label(s) << ": max = " << sig_max
	       << " — assuming 0-1 scale, rescaling to 0-100\n";

      //
      // Per-epoch artifact flags
      //

      std::vector<bool> flag_range;
      std::vector<bool> flag_flat;
      std::vector<bool> flag_jump;
      std::vector<bool> flag_missing;
      std::vector<bool> flag_bad_epoch;
      std::vector<interval_t> epoch_intervals;
      std::vector<double> epoch_emin, epoch_emax;

      const int ne = edf.timeline.first_epoch();

      if ( s == 0 )
	logger << "  considering " << ne << " windows per channel\n";

      while ( 1 )
	{
	  int epoch = edf.timeline.next_epoch();
	  if ( epoch == -1 ) break;

	  interval_t interval = edf.timeline.epoch( epoch );
	  epoch_intervals.push_back( interval );
	  slice_t slice( edf , signals(s) , interval , 1 , false , false );
	  std::vector<double> d = *slice.pdata();
	  const int n = d.size();

	  // apply scale
	  if ( rescaled )
	    for (int i=0; i<n; i++) d[i] *= 100.0;

	  // 1. Out-of-range: value < range_lo or > range_hi for >= range_prop of samples
	  int n_range = 0;
	  for (int i=0; i<n; i++)
	    if ( d[i] < spo2_range_lo || d[i] > spo2_range_hi ) ++n_range;
	  const bool flag_range_e = n > 0 && ( (double)n_range / n >= spo2_range_prop );

	  // 2. Flatline: computed post-loop via sliding window; store per-epoch min/max here
	  double emin = d[0] , emax = d[0];
	  for (int i=1; i<n; i++)
	    { if ( d[i] < emin ) emin = d[i]; if ( d[i] > emax ) emax = d[i]; }
	  epoch_emin.push_back( emin );
	  epoch_emax.push_back( emax );
	  const bool flag_flat_e = false; // set below via sliding window

	  // 3. Jump: |d[i] - d[i-1]| >= jump_th at least once
	  bool flag_jump_e = false;
	  for (int i=1; i<n; i++)
	    if ( fabs( d[i] - d[i-1] ) >= spo2_jump_th ) { flag_jump_e = true; break; }

	  // 4. Missing/invalid: value <= invalid_floor for >= missing_prop of samples
	  int n_missing = 0;
	  for (int i=0; i<n; i++)
	    if ( d[i] <= spo2_invalid_floor ) ++n_missing;
	  const bool flag_missing_e = n > 0 && ( (double)n_missing / n >= spo2_missing_prop );

	  const bool epoch_bad = flag_range_e || flag_flat_e || flag_jump_e || flag_missing_e;

	  flag_range.push_back( flag_range_e );
	  flag_flat.push_back( flag_flat_e );
	  flag_jump.push_back( flag_jump_e );
	  flag_missing.push_back( flag_missing_e );
	  flag_bad_epoch.push_back( epoch_bad );

	} // next epoch


      //
      // Sliding-window flatline check: flag all K epochs in any run where
      // max(emax) - min(emin) across K consecutive epochs < spo2_flatline_th
      //

      const int total_epochs = flag_flat.size();
      if ( spo2_flatline_k >= 1 && total_epochs >= spo2_flatline_k )
	{
	  for (int i=0; i <= total_epochs - spo2_flatline_k; i++)
	    {
	      double wmin = epoch_emin[i] , wmax = epoch_emax[i];
	      for (int j=1; j<spo2_flatline_k; j++)
		{
		  if ( epoch_emin[i+j] < wmin ) wmin = epoch_emin[i+j];
		  if ( epoch_emax[i+j] > wmax ) wmax = epoch_emax[i+j];
		}
	      if ( wmax - wmin < spo2_flatline_th )
		for (int j=0; j<spo2_flatline_k; j++)
		  flag_flat[i+j] = true;
	    }
	}

      // recompute flag_bad_epoch with updated flatline flags
      for (int e=0; e<total_epochs; e++)
	flag_bad_epoch[e] = flag_range[e] || flag_flat[e] || flag_jump[e] || flag_missing[e];


      //
      // Channel-level bad rules
      //

      int n_bad_epoch = 0;
      int max_consec = 0 , cur_consec = 0;
      for (int e=0; e<total_epochs; e++)
	{
	  if ( flag_bad_epoch[e] )
	    { ++n_bad_epoch; if ( ++cur_consec > max_consec ) max_consec = cur_consec; }
	  else cur_consec = 0;
	}

      const double p_bad_epoch    = total_epochs > 0 ? (double)n_bad_epoch / total_epochs : 0 ;
      const double max_bad_run_sec = max_consec * spo2_window_inc;
      const bool bad = ( p_bad_epoch >= spo2_flag_prop ) ||
	               ( max_bad_run_sec >= spo2_flag_run );


      //
      // OXY-specific outputs  (CH x OXY table)
      //

      writer.level( 1 , "OXY" );

      writer.value( "RESCALED" , (int)rescaled );

      // per-reason epoch proportions (OXY-specific)
      int nr_o = 0 , nf_o = 0 , nj_o = 0 , nm_o = 0;
      for (int e=0; e<total_epochs; e++)
	{
	  if ( flag_range[e]   ) ++nr_o;
	  if ( flag_flat[e]    ) ++nf_o;
	  if ( flag_jump[e]    ) ++nj_o;
	  if ( flag_missing[e] ) ++nm_o;
	}
      writer.value( "P_FLAG_EPOCH" , p_bad_epoch );
      writer.value( "P_RANGE"     , total_epochs > 0 ? (double)nr_o / total_epochs : 0 );
      writer.value( "P_FLAT"      , total_epochs > 0 ? (double)nf_o / total_epochs : 0 );
      writer.value( "P_JUMP"      , total_epochs > 0 ? (double)nj_o / total_epochs : 0 );
      writer.value( "P_MISSING"   , total_epochs > 0 ? (double)nm_o / total_epochs : 0 );

      if ( by_epoch )
	{
	  for (int e=0; e<total_epochs; e++)
	    {
	      writer.level( e+1 , "WIN" );
	      writer.value( "RANGE"     , (int)flag_range[e] );
	      writer.value( "FLAT"      , (int)flag_flat[e] );
	      writer.value( "JUMP"      , (int)flag_jump[e] );
	      writer.value( "MISSING"   , (int)flag_missing[e] );
	      writer.value( "FLAG_EPOCH" , (int)flag_bad_epoch[e] );
	    }
	  writer.unlevel( "WIN" );
	}

      writer.unlevel( "OXY" );

      //
      // Cross-domain summary outputs  (CH x DOMAIN table)
      //

      writer.level( "OXY" , "DOMAIN" );
      writer.value( "FLAGGED"         , (int)bad );
      writer.value( "N_FLAG_EPOCH" , n_bad_epoch );
      writer.value( "MAX_FLAG_RUN" , max_bad_run_sec );
      writer.unlevel( "DOMAIN" );

      //
      // Annotations: bad epochs (merged, gap-aware)
      //

      if ( add_annots )
	{
	  std::string alab = annot_prefix;
	  if ( annot_show_domain )  alab += "_OXY";
	  if ( annot_by_channel )   alab += "_" + signals.label(s);

	  annot_t * a = edf.annotations->add( alab );

	  bool in_bad = false;
	  uint64_t bad_start = 0LLU , bad_stop = 0LLU;
	  int ac = 0;

	  for (int e=0; e<total_epochs; e++)
	    {
	      if ( flag_bad_epoch[e] )
		{
		  if ( ! in_bad )
		    {
		      bad_start = epoch_intervals[e].start;
		      bad_stop  = epoch_intervals[e].stop;
		      in_bad = true;
		    }
		  else
		    {
		      if ( epoch_intervals[e].start > bad_stop )
			{
			  a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );
			  ++ac;
			  bad_start = epoch_intervals[e].start;
			}
		      if ( epoch_intervals[e].stop > bad_stop )
			bad_stop = epoch_intervals[e].stop;
		    }
		}
	      else if ( in_bad )
		{
		  a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );
		  ++ac;
		  in_bad = false;
		}
	    }

	  if ( in_bad )
	    {
	      a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );
	      ++ac;
	    }

	  logger << "  added " << ac << " " << alab << " annotations for "
		 << signals.label(s) << "\n";
	}

    } // next signal

  writer.unlevel( globals::signal_strat );

}



// ---------------------------------------------------------------------------
//  line_noise_ratio()
//
//  Computes the ratio P_line / P_broad for a single epoch's data using
//  Welch's PSD (4-second sub-windows, 50% overlap, Tukey50 window).
//
//  P_line  = max( integral PSD over [50-bw, 50+bw] Hz ,
//                 integral PSD over [60-bw, 60+bw] Hz )
//              (only bands reachable below Nyquist are tested)
//  P_broad = integral PSD over [broad_lo, min(broad_hi, Nyquist-0.5)] Hz
//
//  Returns -1 if the signal is too short or neither line-noise band is
//  below Nyquist.
//
//  NOTE: if data have been notch-filtered at 50/60 Hz prior to calling
//  this function the ratio will be artificially suppressed and the test
//  will always pass.  Run on raw (unfiltered) data for meaningful results.
// ---------------------------------------------------------------------------

double dsptools::qc_t::line_noise_ratio( const std::vector<double> & data , double Fs ) const
{
  const double nyquist = Fs / 2.0;

  // broadband upper limit is min(user-specified ceiling, Nyquist - 0.5 Hz)
  const double broad_hi = emg_ln_broad_hi < nyquist - 0.5 ? emg_ln_broad_hi : nyquist - 0.5;

  // need at least some broadband range
  if ( broad_hi <= emg_ln_broad_lo ) return -1.0;

  // need at least one 4-second sub-window
  const double seg_sec  = 4.0;
  const double step_sec = 2.0; // 50% overlap
  const int total_points   = (int)data.size();
  const int segment_points = (int)( seg_sec  * Fs );
  const int noverlap_pts   = (int)( step_sec * Fs );

  if ( segment_points > total_points || segment_points <= noverlap_pts ) return -1.0;

  // number of 50%-overlapping segments that fit in the epoch
  const int noverlap_segments = (int)floor( ( total_points - noverlap_pts )
					    / (double)( segment_points - noverlap_pts ) );
  if ( noverlap_segments < 1 ) return -1.0;

  PWELCH pwelch( data , Fs , seg_sec , noverlap_segments , WINDOW_TUKEY50 );

  const double p_broad = pwelch.psdsum( emg_ln_broad_lo , broad_hi );
  if ( p_broad <= 0.0 ) return -1.0;

  // 50 Hz band
  const double f50_lo = 50.0 - emg_ln_bw;
  const double f50_hi = 50.0 + emg_ln_bw;
  const double p_50   = ( f50_hi < nyquist ) ? pwelch.psdsum( f50_lo , f50_hi ) : 0.0;

  // 60 Hz band
  const double f60_lo = 60.0 - emg_ln_bw;
  const double f60_hi = 60.0 + emg_ln_bw;
  const double p_60   = ( f60_hi < nyquist ) ? pwelch.psdsum( f60_lo , f60_hi ) : 0.0;

  if ( p_50 <= 0.0 && p_60 <= 0.0 ) return -1.0; // neither band reachable

  const double p_line = p_50 > p_60 ? p_50 : p_60;

  return p_line / p_broad;
}


// ---------------------------------------------------------------------------
//  ensure_units()
//
//  Checks the physical dimension of signal slot s against target_unit.
//  If the signal is in a compatible voltage unit (uV / mV / V) but not
//  already in target_unit, it is rescaled in-place via edf_t::rescale().
//  If the unit is unrecognised a warning is logged and processing continues
//  (amplitude-based thresholds in the caller may not be meaningful).
//
//  Typical usage:
//    ensure_units( signals(s), "uV" );   // EEG / EMG / EOG
//    ensure_units( signals(s), "mV" );   // ECG (mV is the clinical standard)
// ---------------------------------------------------------------------------

void dsptools::qc_t::ensure_units( int s , const std::string & target_unit )
{
  const std::string pdim = Helper::trim( edf.header.phys_dimension[s] );

  // already correct — nothing to do
  if ( Helper::imatch( pdim , target_unit ) ) return;

  const bool is_uV = Helper::imatch( pdim , "uV" );
  const bool is_mV = Helper::imatch( pdim , "mV" );
  const bool is_V  = Helper::imatch( pdim , "V"  );

  if ( is_uV || is_mV || is_V )
    edf.rescale( s , target_unit ); // updates in-memory signal data + header
  else
    logger << "  ** warning: " << edf.header.label[s]
	   << " physical dimension is '" << ( pdim.empty() ? "(empty)" : pdim )
	   << "' -- not a recognised voltage unit (uV/mV/V);"
	   << " amplitude thresholds assuming " << target_unit
	   << " may not apply\n";
}


void dsptools::qc_t::do_emg( signal_list_t & signals )
{

  const int ns = signals.size();

  if ( ns == 0 ) return;

  logger << "  ----------------------------------------\n"
	 << "  checking " << ns << " EMG channel(s):\n"
	 << "     window size (sec)           = " << emg_window_dur     << "  (emg-win)\n"
	 << "     window step (sec)           = " << emg_window_inc     << "  (emg-inc)\n"
	 << "     minimum SR (Hz)             = " << emg_min_sr         << "  (emg-min-sr)\n"
	 << "     flatline SD threshold (µV)  = " << emg_flat_std_th    << "  (emg-flat-th)\n"
	 << "     flatline deriv. proportion  = " << emg_flat_deriv_prop<< "  (emg-flat-prop)\n"
	 << "     flatline deriv. epsilon     = " << emg_flat_deriv_eps << "  (emg-flat-eps)\n"
	 << "     clip proportion             = " << emg_clip_prop      << "  (emg-clip-prop)\n"
	 << "     RMS threshold (µV)          = " << emg_rms_th         << "  (emg-rms-th)\n"
	 << "     impulse |z| threshold       = " << emg_impulse_z      << "  (emg-imp-z)\n"
	 << "     impulse min count           = " << emg_impulse_n      << "  (emg-imp-n)\n"
	 << "     line noise ± bw (Hz)        = " << emg_ln_bw          << "  (emg-ln-bw)\n"
	 << "     line noise broadband (Hz)   = " << emg_ln_broad_lo << " - " << emg_ln_broad_hi
                                                                              << "  (emg-ln-lo, emg-ln-hi)\n"
	 << "     line noise ratio threshold  = " << emg_ln_ratio_th    << "  (emg-ln-th)\n"
	 << "     channel bad epoch prop.     = " << emg_flag_prop << "  (emg-flag-prop)\n"
	 << "     channel bad run (sec)       = " << emg_flag_run    << "  (emg-flag-run)\n";

  // set epoch grid
  edf.timeline.set_epochs_to_span_gaps( false );
  edf.timeline.set_epoch( emg_window_dur , emg_window_inc );

  //
  // iterate over channels
  //

  for (int s = 0; s < ns; s++)
    {

      const double Fs = edf.header.sampling_freq( signals(s) );

      if ( Fs < emg_min_sr )
	{
	  logger << "  ** skipping " << signals.label(s)
		 << " (SR=" << Fs << " < emg-min-sr=" << emg_min_sr << ")\n";
	  continue;
	}

      // Ensure µV scale (auto-converts mV/V; warns if unrecognised)
      ensure_units( signals(s) , "uV" );

      writer.level( signals.label(s) , globals::signal_strat );
      writer.level( 1 , "EMG" );

      const int ne = edf.timeline.first_epoch();

      if ( s == 0 )
	logger << "  considering " << ne << " epochs per channel\n";

      // per-epoch flag vectors
      std::vector<bool>       flag_flat;
      std::vector<bool>       flag_clip;
      std::vector<bool>       flag_rms;
      std::vector<bool>       flag_ln;
      std::vector<bool>       flag_imp;
      std::vector<bool>       flag_bad;
      std::vector<double>     ln_ratio_vec;
      std::vector<interval_t> epoch_intervals;

      while ( 1 )
	{
	  int epoch = edf.timeline.next_epoch();
	  if ( epoch == -1 ) break;

	  interval_t interval = edf.timeline.epoch( epoch );
	  epoch_intervals.push_back( interval );

	  slice_t slice( edf , signals(s) , interval , 1 , false , false );
	  const std::vector<double> d = *slice.nonconst_pdata();
	  const int n = (int)d.size();

	  //
	  // 1. Flatline: SD < threshold OR near-zero derivative prop >= threshold
	  //

	  const double ep_sd     = MiscMath::sdev( d );
	  const double ep_deriv  = MiscMath::flat( d , emg_flat_deriv_eps );
	  const bool   f_flat    = ( ep_sd < emg_flat_std_th ) ||
	                           ( ep_deriv >= emg_flat_deriv_prop );

	  //
	  // 2. Clipping: >= emg_clip_prop of samples at ADC limits
	  //

	  const bool f_clip = ( MiscMath::clipped( d ) >= emg_clip_prop );

	  //
	  // 3a. Broadband RMS: RMS > emg_rms_th
	  //

	  const bool f_rms = ( MiscMath::rms( d ) > emg_rms_th );

	  //
	  // 3b. Line noise: P_line/P_broad > emg_ln_ratio_th
	  //

	  const double ln_ratio  = line_noise_ratio( d , Fs );
	  const bool   f_ln      = ( ln_ratio >= 0.0 ) && ( ln_ratio > emg_ln_ratio_th );

	  //
	  // 4. Impulse artifact: >= emg_impulse_n samples with |z| > emg_impulse_z
	  //

	  bool f_imp = false;
	  {
	    const double mu = MiscMath::mean( d );
	    const double sd = ep_sd; // reuse value computed above
	    if ( sd > 0.0 )
	      {
		int n_imp = 0;
		for (int i = 0; i < n; i++)
		  if ( std::fabs( (d[i] - mu) / sd ) > emg_impulse_z )
		    if ( ++n_imp >= emg_impulse_n ) { f_imp = true; break; }
	      }
	  }

	  const bool f_any = f_flat || f_clip || f_rms || f_imp; // f_ln tracked separately

	  flag_flat.push_back( f_flat );
	  flag_clip.push_back( f_clip );
	  flag_rms.push_back( f_rms );
	  flag_ln.push_back( f_ln );
	  flag_imp.push_back( f_imp );
	  flag_bad.push_back( f_any );
	  ln_ratio_vec.push_back( ln_ratio );

	  //
	  // Epoch-level verbose output
	  //

	  if ( by_epoch )
	    {
	      writer.level( epoch + 1 , "WIN" );
	      writer.value( "SD"        , ep_sd    );
	      writer.value( "DERIV"     , ep_deriv );
	      writer.value( "RMS"       , MiscMath::rms( d ) );
	      writer.value( "LN_RATIO"  , ln_ratio );
	      writer.value( "FLAT"      , (int)f_flat );
	      writer.value( "CLIP"      , (int)f_clip );
	      writer.value( "HI_RMS"    , (int)f_rms  );
	      writer.value( "HI_LN"     , (int)f_ln   );
	      writer.value( "IMP"       , (int)f_imp  );
	      writer.value( "FLAG_EPOCH" , (int)f_any  );
	    }

	} // next epoch

      if ( by_epoch )
	writer.unlevel( "WIN" );

      //
      // Channel-level summary
      //

      const int total_epochs = (int)flag_bad.size();

      int n_bad = 0 , n_flat = 0 , n_clip = 0 , n_rms = 0 , n_ln = 0 , n_imp = 0;
      int max_consec = 0 , cur_consec = 0;
      int max_ln_consec = 0 , cur_ln_consec = 0;

      for (int e = 0; e < total_epochs; e++)
	{
	  if ( flag_flat[e] ) ++n_flat;
	  if ( flag_clip[e] ) ++n_clip;
	  if ( flag_rms[e]  ) ++n_rms;
	  if ( flag_imp[e]  ) ++n_imp;
	  if ( flag_ln[e] )
	    {
	      ++n_ln;
	      if ( ++cur_ln_consec > max_ln_consec ) max_ln_consec = cur_ln_consec;
	    }
	  else cur_ln_consec = 0;
	  if ( flag_bad[e]  )
	    {
	      ++n_bad;
	      if ( ++cur_consec > max_consec ) max_consec = cur_consec;
	    }
	  else cur_consec = 0;
	}

      const double prop_bad      = total_epochs > 0 ? n_bad / (double)total_epochs : 0.0;
      const double max_bad_run   = max_consec * emg_window_inc;
      const bool   ch_bad        = ( prop_bad > emg_flag_prop ) ||
	                            ( max_bad_run >= emg_flag_run );

      const double prop_ln       = total_epochs > 0 ? n_ln / (double)total_epochs : 0.0;
      const double max_ln_run    = max_ln_consec * emg_window_inc;
      const bool   ch_ln_bad     = ( prop_ln > emg_flag_prop ) ||
	                            ( max_ln_run >= emg_flag_run );

      writer.value( "N_FLAG_EPOCH" , n_bad );
      writer.value( "PROP_FLAG"    , prop_bad );
      writer.value( "MAX_FLAG_RUN" , max_bad_run );
      writer.value( "FLAGGED"         , (int)ch_bad );
      writer.value( "P_FLAT"      , total_epochs > 0 ? n_flat / (double)total_epochs : 0.0 );
      writer.value( "P_CLIP"      , total_epochs > 0 ? n_clip / (double)total_epochs : 0.0 );
      writer.value( "P_HI_RMS"    , total_epochs > 0 ? n_rms  / (double)total_epochs : 0.0 );
      writer.value( "P_IMP"       , total_epochs > 0 ? n_imp  / (double)total_epochs : 0.0 );
      writer.value( "P_LN"        , prop_ln );
      writer.value( "MAX_LN_RUN"  , max_ln_run );
      writer.value( "LN_FLAG"      , (int)ch_ln_bad );

      writer.unlevel( "EMG" );

      //
      // Cross-domain summary (CH x DOMAIN table)
      //

      writer.level( "EMG" , "DOMAIN" );
      writer.value( "FLAGGED"         , (int)ch_bad    );
      writer.value( "N_FLAG_EPOCH" , n_bad          );
      writer.value( "MAX_FLAG_RUN" , max_bad_run    );
      writer.value( "LN_FLAG"      , (int)ch_ln_bad );
      writer.value( "MAX_LN_RUN"  , max_ln_run     );
      writer.unlevel( "DOMAIN" );

      //
      // Annotations (merged, gap-aware — same pattern as RESP)
      //

      if ( add_annots )
	{
	  std::string alab = annot_prefix;
	  if ( annot_show_domain ) alab += "_EMG";
	  if ( annot_by_channel  ) alab += "_" + signals.label(s);

	  annot_t * a = edf.annotations->add( alab );

	  bool in_bad = false;
	  uint64_t bad_start = 0LLU , bad_stop = 0LLU;

	  for (int e = 0; e < total_epochs; e++)
	    {
	      if ( flag_bad[e] )
		{
		  if ( ! in_bad )
		    {
		      bad_start = epoch_intervals[e].start;
		      bad_stop  = epoch_intervals[e].stop;
		      in_bad = true;
		    }
		  else
		    {
		      if ( epoch_intervals[e].start > bad_stop )
			{
			  a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );
			  bad_start = epoch_intervals[e].start;
			}
		      if ( epoch_intervals[e].stop > bad_stop )
			bad_stop = epoch_intervals[e].stop;
		    }
		}
	      else if ( in_bad )
		{
		  a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );
		  in_bad = false;
		}
	    }

	  if ( in_bad )
	    a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );

	  // separate LN annotations
	  if ( n_ln > 0 )
	    {
	      std::string ln_lab = annot_prefix;
	      if ( annot_show_domain ) ln_lab += "_EMG";
	      ln_lab += "_LN";
	      if ( annot_by_channel ) ln_lab += "_" + signals.label(s);

	      annot_t * a_ln = edf.annotations->add( ln_lab );
	      bool in_ln = false;
	      uint64_t ln_start = 0LLU , ln_stop = 0LLU;

	      for ( int e = 0; e < total_epochs; e++ )
		{
		  if ( flag_ln[e] )
		    {
		      if ( ! in_ln ) { ln_start = epoch_intervals[e].start; ln_stop = epoch_intervals[e].stop; in_ln = true; }
		      else { if ( epoch_intervals[e].start > ln_stop ) { a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) ); ln_start = epoch_intervals[e].start; } if ( epoch_intervals[e].stop > ln_stop ) ln_stop = epoch_intervals[e].stop; }
		    }
		  else if ( in_ln ) { a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) ); in_ln = false; }
		}
	      if ( in_ln ) a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) );
	    }
	}

      logger << "  " << signals.label(s)
	     << " : " << n_bad << "/" << total_epochs << " flagged epochs"
	     << "  (flat=" << n_flat << " clip=" << n_clip
	     << " rms=" << n_rms << " imp=" << n_imp << ") | ln=" << n_ln
	     << ( ch_ln_bad ? " --> LN_FLAG" : "" )
	     << ( ch_bad ? " --> FLAGGED channel" : "" ) << "\n";

    } // next signal

  writer.unlevel( globals::signal_strat );

}



void dsptools::qc_t::do_ecg( signal_list_t & signals )
{

  const int ns = signals.size();

  if ( ns == 0 ) return;

  logger << "  ----------------------------------------\n"
	 << "  checking " << ns << " ECG channel(s):\n"
	 << "     window size (sec)           = " << ecg_window_dur      << "  (ecg-win)\n"
	 << "     window step (sec)           = " << ecg_window_inc      << "  (ecg-inc)\n"
	 << "     minimum SR (Hz)             = " << ecg_min_sr          << "  (ecg-min-sr)\n"
	 << "     flatline SD threshold (mV)  = " << ecg_flat_std_th     << "  (ecg-flat-th)\n"
	 << "     flatline deriv. proportion  = " << ecg_flat_deriv_prop << "  (ecg-flat-prop)\n"
	 << "     flatline deriv. epsilon     = " << ecg_flat_deriv_eps  << "  (ecg-flat-eps)\n"
	 << "     clip proportion             = " << ecg_clip_prop       << "  (ecg-clip-prop)\n"
	 << "     HR range (bpm)              = " << ecg_hr_min << " - " << ecg_hr_max
	                                                                   << "  (ecg-hr-min, ecg-hr-max)\n"
	 << "     RR range (ms)               = " << ecg_rr_min_ms << " - " << ecg_rr_max_ms
	                                                                   << "  (ecg-rr-min, ecg-rr-max)\n"
	 << "     bad-RR prop. threshold      = " << ecg_rr_flag_prop     << "  (ecg-rr-prop)\n"
	 << "     min beats per epoch         = " << ecg_min_beats       << "  (ecg-min-beats)\n"
	 << "     line noise ± bw (Hz)        = " << ecg_ln_bw           << "  (ecg-ln-bw)\n"
	 << "     cardiac signal band (Hz)    = " << ecg_sig_lo << " - " << ecg_sig_hi
	                                                                   << "  (ecg-sig-lo, ecg-sig-hi)\n"
	 << "     line noise ratio threshold  = " << ecg_ln_ratio_th     << "  (ecg-ln-th)\n"
	 << "     channel bad epoch prop.     = " << ecg_flag_prop  << "  (ecg-flag-prop)\n"
	 << "     channel bad run (sec)       = " << ecg_flag_run     << "  (ecg-flag-run)\n";

  edf.timeline.set_epochs_to_span_gaps( false );
  edf.timeline.set_epoch( ecg_window_dur , ecg_window_inc );

  for ( int s = 0; s < ns; s++ )
    {

      const double Fs      = edf.header.sampling_freq( signals(s) );
      const double nyquist = Fs / 2.0;

      if ( Fs < ecg_min_sr )
	{
	  logger << "  ** skipping " << signals.label(s)
		 << " (SR=" << Fs << " < ecg-min-sr=" << ecg_min_sr << ")\n";
	  continue;
	}

      // ECG: mV is the clinical standard; also accepts uV or V (auto-converts)
      ensure_units( signals(s) , "mV" );

      //
      // Whole-signal R-peak detection (done once, used across all epochs)
      //

      const interval_t whole = edf.timeline.wholetrace();
      slice_t wslice( edf , signals(s) , whole );
      const std::vector<double>   * wd  = wslice.pdata();
      const std::vector<uint64_t> * wtp = wslice.ptimepoints();

      rpeak_opt_t ropt;   // use defaults
      rpeaks_t whole_peaks = dsptools::mpeakdetect2( wd , wtp , (int)Fs , ropt );

      const int whole_np = (int)whole_peaks.R_t.size();
      logger << "  " << signals.label(s)
	     << ": detected " << whole_np << " R-peaks in full signal"
	     << ( whole_peaks.inverted ? " (inverted)" : "" ) << "\n";

      //
      // Optionally add one point annotation per R-peak (for visual QC)
      //

      if ( ecg_add_peaks )
	{
	  annot_t * pa = edf.annotations->add( "Rpk_" + signals.label(s) );
	  for ( int i = 0; i < whole_np; i++ )
	    pa->add( "." , interval_t( whole_peaks.R_t[i] , whole_peaks.R_t[i] ) ,
		    signals.label(s) );
	}

      //
      // Per-epoch pass
      //

      writer.level( signals.label(s) , globals::signal_strat );
      writer.level( 1 , "ECG" );

      const int ne = edf.timeline.first_epoch();
      if ( s == 0 ) logger << "  considering " << ne << " epochs per channel\n";

      std::vector<bool>      flag_flat, flag_clip, flag_rr, flag_ln;
      std::vector<double>    ep_np, ep_hr, ep_rr_mean, ep_rr_bad;
      std::vector<interval_t> epoch_intervals;

      while ( 1 )
	{
	  int epoch = edf.timeline.next_epoch();
	  if ( epoch == -1 ) break;

	  interval_t interval = edf.timeline.epoch( epoch );
	  epoch_intervals.push_back( interval );

	  slice_t slice( edf , signals(s) , interval , 1 , false , false );
	  const std::vector<double> d = *slice.nonconst_pdata();
	  const int n = (int)d.size();

	  bool   f_flat = false , f_clip = false , f_rr = false , f_ln = false;
	  double ep_hr_val = 0.0 , ep_rr_val = 0.0 , ep_rr_bad_val = 0.0;
	  int    np_ep = 0;

	  if ( n == 0 )
	    {
	      flag_flat.push_back( false );
	      flag_clip.push_back( false );
	      flag_rr.push_back( true );
	      flag_ln.push_back( false );
	      ep_np.push_back( 0 );
	      ep_hr.push_back( 0 );
	      ep_rr_mean.push_back( 0 );
	      ep_rr_bad.push_back( 1.0 );
	      continue;
	    }

	  //
	  // 1. Flatline: SD < threshold OR near-zero derivative prop >= threshold
	  //

	  const double sd    = MiscMath::sdev( d );
	  const double deriv = MiscMath::flat( d , ecg_flat_deriv_eps );
	  f_flat = ( sd < ecg_flat_std_th ) || ( deriv >= ecg_flat_deriv_prop );

	  //
	  // 2. Clipping
	  //

	  f_clip = ( MiscMath::clipped( d ) >= ecg_clip_prop );

	  //
	  // 3. Per-epoch RR from whole-signal peaks
	  //

	  std::vector<uint64_t> ep_beats = whole_peaks.beats( interval );
	  np_ep = (int)ep_beats.size();

	  if ( np_ep >= 2 )
	    {
	      int    n_bad_rr  = 0;
	      double sum_rr    = 0.0;   // all intervals (for output)
	      double sum_rr_ok = 0.0;   // clean intervals only (for HR estimate)
	      int    n_rr_ok   = 0;
	      for ( int i = 1 ; i < np_ep ; i++ )
		{
		  const double rr_ms = globals::tp_duration
		                     * (double)( ep_beats[i] - ep_beats[i-1] ) * 1000.0;
		  sum_rr += rr_ms;
		  if ( rr_ms < ecg_rr_min_ms || rr_ms > ecg_rr_max_ms )
		    ++n_bad_rr;
		  else
		    { sum_rr_ok += rr_ms; ++n_rr_ok; }
		}
	      ep_rr_val     = sum_rr / (double)( np_ep - 1 );  // mean of all (output)
	      ep_rr_bad_val = (double)n_bad_rr / (double)( np_ep - 1 );

	      // HR estimated from clean RR intervals only; if none are clean, use raw mean
	      const double rr_for_hr = n_rr_ok > 0 ? sum_rr_ok / n_rr_ok : ep_rr_val;
	      ep_hr_val = rr_for_hr > 0.0 ? 60000.0 / rr_for_hr : 0.0;

	      f_rr = ( ep_hr_val < ecg_hr_min  ||
		       ep_hr_val > ecg_hr_max  ||
		       ep_rr_bad_val > ecg_rr_flag_prop );
	    }
	  else
	    {
	      ep_hr_val    = 0.0;
	      ep_rr_val    = 0.0;
	      ep_rr_bad_val = 1.0;
	      f_rr = true;   // < 2 beats → always flag
	    }

	  // additionally flag if fewer than ecg_min_beats present
	  if ( np_ep < ecg_min_beats ) f_rr = true;

	  //
	  // 4. Line noise (separate from BAD, same pattern as other domains)
	  //

	  {
	    const double seg_sec  = 4.0;
	    const double step_sec = 2.0;
	    const int    seg_pts  = (int)( seg_sec  * Fs );
	    const int    ovlp_pts = (int)( step_sec * Fs );

	    const double sig_hi_safe = ecg_sig_hi < nyquist - 0.5 ? ecg_sig_hi : nyquist - 0.5;

	    if ( seg_pts <= n && seg_pts > ovlp_pts && sig_hi_safe > ecg_sig_lo )
	      {
		const int nsegs = (int)floor( (double)( n - ovlp_pts )
					     / (double)( seg_pts - ovlp_pts ) );
		if ( nsegs >= 1 )
		  {
		    PWELCH pw( d , Fs , seg_sec , nsegs , WINDOW_TUKEY50 );

		    const double p_sig = pw.psdsum( ecg_sig_lo , sig_hi_safe );
		    if ( p_sig > 0.0 )
		      {
			const double p_ln50 = nyquist > 50.0 + ecg_ln_bw
			  ? pw.psdsum( 50.0 - ecg_ln_bw , 50.0 + ecg_ln_bw ) : 0.0;
			const double p_ln60 = nyquist > 60.0 + ecg_ln_bw
			  ? pw.psdsum( 60.0 - ecg_ln_bw , 60.0 + ecg_ln_bw ) : 0.0;
			const double p_ln_max = std::max( p_ln50 , p_ln60 );
			f_ln = ( p_ln_max / p_sig > ecg_ln_ratio_th );
		      }
		  }
	      }
	  }

	  // BAD = flatline OR clip OR RR issues  (line noise tracked separately)
	  flag_flat.push_back( f_flat );
	  flag_clip.push_back( f_clip );
	  flag_rr.push_back( f_rr );
	  flag_ln.push_back( f_ln );

	  ep_np.push_back( (double)np_ep );
	  ep_hr.push_back( ep_hr_val );
	  ep_rr_mean.push_back( ep_rr_val );
	  ep_rr_bad.push_back( ep_rr_bad_val );

	} // end epoch loop

      const int total_epochs = (int)flag_flat.size();

      // Combined BAD flag
      std::vector<bool> flag_bad( total_epochs );
      for ( int e = 0; e < total_epochs; e++ )
	flag_bad[e] = flag_flat[e] || flag_clip[e] || flag_rr[e];

      //
      // Per-epoch output
      //

      if ( by_epoch )
	{
	  for ( int e = 0; e < total_epochs; e++ )
	    {
	      writer.level( e + 1 , "WIN" );
	      writer.value( "NP"        , ep_np[e]       );
	      writer.value( "HR"        , ep_hr[e]        );
	      writer.value( "RR"        , ep_rr_mean[e]   );
	      writer.value( "P_RR"      , ep_rr_bad[e]    );
	      writer.value( "FLAT"      , (int)flag_flat[e] );
	      writer.value( "CLIP"      , (int)flag_clip[e] );
	      writer.value( "RR_FLAG"    , (int)flag_rr[e]   );
	      writer.value( "LN"        , (int)flag_ln[e]   );
	      writer.value( "FLAG_EPOCH" , (int)flag_bad[e]  );
	    }
	  writer.unlevel( "WIN" );
	}

      //
      // Channel-level statistics
      //

      int n_bad = 0 , n_flat = 0 , n_clip = 0 , n_rr = 0 , n_ln = 0;
      int max_consec    = 0 , cur_consec    = 0;
      int max_ln_consec = 0 , cur_ln_consec = 0;

      for ( int e = 0; e < total_epochs; e++ )
	{
	  if ( flag_flat[e] ) ++n_flat;
	  if ( flag_clip[e] ) ++n_clip;
	  if ( flag_rr[e]   ) ++n_rr;
	  if ( flag_ln[e] )
	    {
	      ++n_ln;
	      if ( ++cur_ln_consec > max_ln_consec ) max_ln_consec = cur_ln_consec;
	    }
	  else cur_ln_consec = 0;
	  if ( flag_bad[e] )
	    {
	      ++n_bad;
	      if ( ++cur_consec > max_consec ) max_consec = cur_consec;
	    }
	  else cur_consec = 0;
	}

      const double prop_bad    = total_epochs > 0 ? n_bad / (double)total_epochs : 0.0;
      const double max_bad_run = max_consec    * ecg_window_inc;
      const bool   ch_bad      = ( prop_bad > ecg_flag_prop ) ||
	                          ( max_bad_run >= ecg_flag_run );

      const double prop_ln     = total_epochs > 0 ? n_ln / (double)total_epochs : 0.0;
      const double max_ln_run  = max_ln_consec * ecg_window_inc;
      const bool   ch_ln_bad   = ( prop_ln > ecg_flag_prop ) ||
	                          ( max_ln_run >= ecg_flag_run );

      writer.value( "N_PEAKS"     , whole_np    );
      writer.value( "N_FLAG_EPOCH" , n_bad       );
      writer.value( "PROP_FLAG"    , prop_bad    );
      writer.value( "MAX_FLAG_RUN" , max_bad_run );
      writer.value( "FLAGGED"         , (int)ch_bad );
      writer.value( "P_FLAT"      , total_epochs > 0 ? n_flat / (double)total_epochs : 0.0 );
      writer.value( "P_CLIP"      , total_epochs > 0 ? n_clip / (double)total_epochs : 0.0 );
      writer.value( "P_RR"        , total_epochs > 0 ? n_rr   / (double)total_epochs : 0.0 );
      writer.value( "P_LN"        , prop_ln     );
      writer.value( "MAX_LN_RUN"  , max_ln_run  );
      writer.value( "LN_FLAG"      , (int)ch_ln_bad );

      writer.unlevel( "ECG" );

      // Cross-domain summary (CH x DOMAIN table)
      writer.level( "ECG" , "DOMAIN" );
      writer.value( "FLAGGED"         , (int)ch_bad    );
      writer.value( "N_FLAG_EPOCH" , n_bad          );
      writer.value( "MAX_FLAG_RUN" , max_bad_run    );
      writer.value( "LN_FLAG"      , (int)ch_ln_bad );
      writer.value( "MAX_LN_RUN"  , max_ln_run     );
      writer.unlevel( "DOMAIN" );

      //
      // Annotations (merged, gap-aware)
      //

      if ( add_annots && total_epochs > 0 )
	{
	  // bad epochs (merged runs)
	  {
	    std::string alab = annot_prefix;
	    if ( annot_show_domain ) alab += "_ECG";
	    if ( annot_by_channel  ) alab += "_" + signals.label(s);

	    annot_t * a = edf.annotations->add( alab );

	    bool in_bad = false;
	    uint64_t bad_start = 0LLU , bad_stop = 0LLU;

	    for ( int e = 0; e < total_epochs; e++ )
	      {
		if ( flag_bad[e] )
		  {
		    if ( ! in_bad )
		      {
			bad_start = epoch_intervals[e].start;
			bad_stop  = epoch_intervals[e].stop;
			in_bad = true;
		      }
		    else
		      {
			if ( epoch_intervals[e].start > bad_stop )
			  {
			    a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );
			    bad_start = epoch_intervals[e].start;
			  }
			if ( epoch_intervals[e].stop > bad_stop )
			  bad_stop = epoch_intervals[e].stop;
		      }
		  }
		else if ( in_bad )
		  {
		    a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );
		    in_bad = false;
		  }
	      }
	    if ( in_bad )
	      a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );

	  }

	  // line-noise epochs (separate LN annotations)
	  if ( n_ln > 0 )
	    {
	      std::string ln_lab = annot_prefix;
	      if ( annot_show_domain ) ln_lab += "_ECG";
	      ln_lab += "_LN";
	      if ( annot_by_channel ) ln_lab += "_" + signals.label(s);

	      annot_t * a_ln = edf.annotations->add( ln_lab );
	      bool in_ln = false;
	      uint64_t ln_start = 0LLU , ln_stop = 0LLU;

	      for ( int e = 0; e < total_epochs; e++ )
		{
		  if ( flag_ln[e] )
		    {
		      if ( ! in_ln ) { ln_start = epoch_intervals[e].start; ln_stop = epoch_intervals[e].stop; in_ln = true; }
		      else { if ( epoch_intervals[e].start > ln_stop ) { a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) ); ln_start = epoch_intervals[e].start; } if ( epoch_intervals[e].stop > ln_stop ) ln_stop = epoch_intervals[e].stop; }
		    }
		  else if ( in_ln ) { a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) ); in_ln = false; }
		}
	      if ( in_ln ) a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) );
	    }
	}

      logger << "  " << signals.label(s)
	     << " : " << n_bad << "/" << total_epochs << " flagged epochs"
	     << "  (flat=" << n_flat << " clip=" << n_clip << " rr=" << n_rr << ")"
	     << "  peaks=" << whole_np
	     << " | ln=" << n_ln
	     << ( ch_ln_bad ? " --> LN_FLAG" : "" )
	     << ( ch_bad    ? " --> FLAGGED channel" : "" ) << "\n";

    } // next channel

  writer.unlevel( globals::signal_strat );

}



void dsptools::qc_t::do_eog( signal_list_t & signals )
{

  const int ns = signals.size();

  if ( ns == 0 ) return;

  logger << "  ----------------------------------------\n"
	 << "  checking " << ns << " EOG channel(s):\n"
	 << "     window size (sec)           = " << eog_window_dur    << "  (eog-win)\n"
	 << "     window step (sec)           = " << eog_window_inc    << "  (eog-inc)\n"
	 << "     minimum SR (Hz)             = " << eog_min_sr        << "  (eog-min-sr)\n"
	 << "     flatline SD threshold (µV)  = " << eog_flat_std_th   << "  (eog-flat-th)\n"
	 << "     flatline deriv. proportion  = " << eog_flat_deriv_prop << "  (eog-flat-prop)\n"
	 << "     flatline deriv. epsilon     = " << eog_flat_deriv_eps  << "  (eog-flat-eps)\n"
	 << "     clip proportion             = " << eog_clip_prop     << "  (eog-clip-prop)\n"
	 << "     amplitude threshold (µV)    = " << eog_amp_th        << "  (eog-amp-th)\n"
	 << "     amplitude proportion        = " << eog_amp_prop      << "  (eog-amp-prop)\n"
	 << "     HF band (Hz)                = " << eog_hf_lo << " - " << eog_hf_hi
                                                                         << "  (eog-hf-lo, eog-hf-hi)\n"
	 << "     signal band (Hz)            = " << eog_sig_lo << " - " << eog_sig_hi
                                                                         << "  (eog-sig-lo, eog-sig-hi)\n"
	 << "     HF/signal ratio threshold   = " << eog_hf_ratio_th   << "  (eog-hf-th)\n"
	 << "     line noise ± bw (Hz)        = " << eog_ln_bw         << "  (eog-ln-bw)\n"
	 << "     line noise broad range (Hz) = " << eog_ln_lo << " - " << eog_ln_hi
                                                                         << "  (eog-ln-lo, eog-ln-hi)\n"
	 << "     line noise ratio threshold  = " << eog_ln_ratio_th   << "  (eog-ln-th)\n"
	 << "     Hjorth pre-filter (Hz)       = " << eog_hjorth_lo << " - " << eog_hjorth_hi
                                                                         << "  (eog-hjorth-lo, eog-hjorth-hi)\n"
	 << "     Hjorth filter order         = " << eog_hjorth_order  << "  (eog-hjorth-ord)\n"
	 << "     Hjorth outlier |z| thresh.  = " << eog_hjorth_z_th   << "  (eog-hjorth-z)\n"
	 << "     channel bad epoch prop.     = " << eog_flag_prop << "  (eog-flag-prop)\n"
	 << "     channel bad run (sec)       = " << eog_flag_run    << "  (eog-flag-run)\n";

  edf.timeline.set_epochs_to_span_gaps( false );
  edf.timeline.set_epoch( eog_window_dur , eog_window_inc );

  for ( int s = 0; s < ns; s++ )
    {

      // Ensure µV scale (auto-converts mV/V; warns if unrecognised)
      ensure_units( signals(s) , "uV" );

      const double Fs      = edf.header.sampling_freq( signals(s) );
      const double nyquist = Fs / 2.0;

      if ( Fs < eog_min_sr )
	{
	  logger << "  ** skipping " << signals.label(s)
		 << " (SR=" << Fs << " < eog-min-sr=" << eog_min_sr << ")\n";
	  continue;
	}

      // HF ratio check requires HF upper bound to be below Nyquist
      const bool do_hf_ratio = ( nyquist > eog_hf_hi );
      if ( ! do_hf_ratio )
	logger << "  ** note: " << signals.label(s)
	       << " Nyquist (" << nyquist << " Hz) <= eog-hf-hi (" << eog_hf_hi
	       << " Hz); HF ratio check skipped\n";

      // IIR Butterworth bandpass for Hjorth pre-filtering (raw signal used for all other checks)
      const double hjorth_lo_safe = eog_hjorth_lo < nyquist ? eog_hjorth_lo : nyquist * 0.1;
      const double hjorth_hi_safe = eog_hjorth_hi < nyquist ? eog_hjorth_hi : nyquist * 0.9;
      iir_t hjorth_iir;
      bool  hjorth_filter_ok = false;
      if ( hjorth_lo_safe < hjorth_hi_safe )
	{
	  hjorth_iir.init( BUTTERWORTH_BANDPASS , eog_hjorth_order , Fs , hjorth_lo_safe , hjorth_hi_safe );
	  hjorth_filter_ok = true;
	}
      else
	logger << "  ** note: " << signals.label(s)
	       << " could not set up Hjorth bandpass filter (Fs=" << Fs << "); Hjorth computed on raw signal\n";

      writer.level( signals.label(s) , globals::signal_strat );
      writer.level( 1 , "EOG" );

      const int ne = edf.timeline.first_epoch();
      if ( s == 0 ) logger << "  considering " << ne << " epochs per channel\n";

      // per-epoch metric storage (two-pass: Hjorth outlier needs full distribution)
      std::vector<bool>       flag_flat, flag_clip, flag_amp, flag_hf, flag_ln;
      std::vector<bool>       flag_hjorth;        // filled in second pass
      std::vector<double>     ep_sd, ep_deriv;
      std::vector<double>     ep_hf_ratio, ep_ln_ratio;
      std::vector<double>     ep_act, ep_cplx;    // Hjorth activity, complexity
      std::vector<interval_t> epoch_intervals;

      //
      // First pass: compute all per-epoch metrics
      //

      while ( 1 )
	{
	  int epoch = edf.timeline.next_epoch();
	  if ( epoch == -1 ) break;

	  interval_t interval = edf.timeline.epoch( epoch );
	  epoch_intervals.push_back( interval );

	  slice_t slice( edf , signals(s) , interval , 1 , false , false );
	  const std::vector<double> d = *slice.nonconst_pdata();
	  const int n = (int)d.size();

	  //
	  // 1. Flatline: SD < threshold OR near-zero derivative prop >= threshold
	  //

	  const double sd    = MiscMath::sdev( d );
	  const double deriv = MiscMath::flat( d , eog_flat_deriv_eps );
	  const bool f_flat  = ( sd < eog_flat_std_th ) || ( deriv >= eog_flat_deriv_prop );

	  //
	  // 2. Clipping: >= eog_clip_prop of samples at ADC limits
	  //

	  const bool f_clip = ( MiscMath::clipped( d ) >= eog_clip_prop );

	  //
	  // 3. Extreme amplitude: >= eog_amp_prop samples with |EOG| > eog_amp_th
	  //

	  int n_amp = 0;
	  for ( int i = 0; i < n; i++ )
	    if ( std::fabs( d[i] ) > eog_amp_th ) ++n_amp;
	  const bool f_amp = ( n > 0 ) && ( (double)n_amp / n >= eog_amp_prop );

	  //
	  // 4. Spectral: HF contamination + line noise (one PWELCH call)
	  //

	  double hf_ratio = -1.0;
	  double ln_ratio = -1.0;
	  bool   f_hf     = false; // HF ratio only (stays in BAD)
	  bool   f_ln     = false; // line noise (tracked separately)

	  {
	    const double seg_sec = 4.0;
	    const double step_sec = 2.0;
	    const int seg_pts  = (int)( seg_sec  * Fs );
	    const int ovlp_pts = (int)( step_sec * Fs );

	    if ( seg_pts <= n && seg_pts > ovlp_pts )
	      {
		const int nsegs = (int)floor( ( n - ovlp_pts )
					     / (double)( seg_pts - ovlp_pts ) );
		if ( nsegs >= 1 )
		  {
		    PWELCH pw( d , Fs , seg_sec , nsegs , WINDOW_TUKEY50 );

		    // 4a. HF ratio: P[20-40 Hz] / P[0.3-20 Hz]
		    if ( do_hf_ratio )
		      {
			const double p_hf  = pw.psdsum( eog_hf_lo  , eog_hf_hi  );
			const double p_sig = pw.psdsum( eog_sig_lo , eog_sig_hi );
			if ( p_sig > 0.0 )
			  {
			    hf_ratio = p_hf / p_sig;
			    if ( hf_ratio > eog_hf_ratio_th ) f_hf = true;
			  }
		      }

		    // 4b. Line noise: max(P50, P60) / P[0.3-40 Hz]
		    const double broad_hi = eog_ln_hi < nyquist - 0.5 ? eog_ln_hi : nyquist - 0.5;
		    if ( broad_hi > eog_ln_lo )
		      {
			const double p_broad = pw.psdsum( eog_ln_lo , broad_hi );
			if ( p_broad > 0.0 )
			  {
			    const double f50_lo = 50.0 - eog_ln_bw;
			    const double f50_hi = 50.0 + eog_ln_bw;
			    const double f60_lo = 60.0 - eog_ln_bw;
			    const double f60_hi = 60.0 + eog_ln_bw;
			    const double p_50 = ( f50_hi < nyquist ) ? pw.psdsum( f50_lo , f50_hi ) : 0.0;
			    const double p_60 = ( f60_hi < nyquist ) ? pw.psdsum( f60_lo , f60_hi ) : 0.0;
			    const double p_line = p_50 > p_60 ? p_50 : p_60;
			    if ( p_50 > 0.0 || p_60 > 0.0 )
			      {
				ln_ratio = p_line / p_broad;
				if ( ln_ratio > eog_ln_ratio_th ) f_ln = true;
			      }
			  }
		      }
		  }
	      }
	  }

	  //
	  // 5. Hjorth parameters (activity, complexity — store for outlier pass)
	  //    Pre-filter with IIR Butterworth to reduce HF noise sensitivity;
	  //    raw signal d is unchanged for all other checks above.
	  //

	  double activity = 0 , mobility = 0 , complexity = 0;
	  if ( hjorth_filter_ok )
	    {
	      const std::vector<double> df = hjorth_iir.apply( d );
	      MiscMath::hjorth( &df , &activity , &mobility , &complexity , ! globals::legacy_hjorth );
	    }
	  else
	    {
	      MiscMath::hjorth( &d , &activity , &mobility , &complexity , ! globals::legacy_hjorth );
	    }

	  flag_flat.push_back( f_flat );
	  flag_clip.push_back( f_clip );
	  flag_amp.push_back( f_amp );
	  flag_hf.push_back( f_hf );
	  flag_ln.push_back( f_ln );
	  ep_sd.push_back( sd );
	  ep_deriv.push_back( deriv );
	  ep_hf_ratio.push_back( hf_ratio );
	  ep_ln_ratio.push_back( ln_ratio );
	  ep_act.push_back( activity );
	  ep_cplx.push_back( complexity );

	} // next epoch

      //
      // Second pass: Hjorth robust z-score outlier detection
      // flag if |robust-z| > eog_hjorth_z_th for activity OR complexity
      //

      const int total_epochs = (int)ep_act.size();
      flag_hjorth.assign( total_epochs , false );

      if ( total_epochs > 2 )
	{
	  // helper: flag outliers in a distribution via MAD robust-z
	  auto flag_by_mad = [&]( const std::vector<double> & vals )
	  {
	    std::vector<double> tmp = vals;
	    const double med = MiscMath::median( tmp );
	    std::vector<double> abs_dev( total_epochs );
	    for ( int e = 0; e < total_epochs; e++ )
	      abs_dev[e] = std::fabs( vals[e] - med );
	    std::vector<double> tmp2 = abs_dev;
	    const double mad   = MiscMath::median( tmp2 );
	    const double scale = mad > 0.0 ? 1.4826 * mad : 1.0;
	    for ( int e = 0; e < total_epochs; e++ )
	      if ( std::fabs( vals[e] - med ) / scale > eog_hjorth_z_th )
		flag_hjorth[e] = true;
	  };

	  flag_by_mad( ep_act  );   // activity outliers
	  flag_by_mad( ep_cplx );   // complexity outliers
	}

      //
      // Combine flags
      //

      std::vector<bool> flag_bad( total_epochs );
      for ( int e = 0; e < total_epochs; e++ )
	flag_bad[e] = flag_flat[e] || flag_clip[e] || flag_amp[e]
	           || flag_hf[e]   || flag_hjorth[e]; // f_ln tracked separately

      //
      // Epoch-level verbose output (deferred to here so Hjorth flags are ready)
      //

      if ( by_epoch )
	{
	  for ( int e = 0; e < total_epochs; e++ )
	    {
	      writer.level( e + 1 , "WIN" );
	      writer.value( "SD"        , ep_sd[e]       );
	      writer.value( "DERIV"     , ep_deriv[e]    );
	      writer.value( "HF_RATIO"  , ep_hf_ratio[e] );
	      writer.value( "LN_RATIO"  , ep_ln_ratio[e] );
	      writer.value( "ACT"       , ep_act[e]      );
	      writer.value( "CPLX"      , ep_cplx[e]     );
	      writer.value( "FLAT"      , (int)flag_flat[e]   );
	      writer.value( "CLIP"      , (int)flag_clip[e]   );
	      writer.value( "AMP"       , (int)flag_amp[e]    );
	      writer.value( "HF"        , (int)flag_hf[e]     );
	      writer.value( "LN"        , (int)flag_ln[e]     );
	      writer.value( "HJORTH"    , (int)flag_hjorth[e] );
	      writer.value( "FLAG_EPOCH" , (int)flag_bad[e]    );
	    }
	  writer.unlevel( "WIN" );
	}

      //
      // Channel-level summary
      //

      int n_bad = 0 , n_flat = 0 , n_clip = 0 , n_amp = 0 , n_hf = 0 , n_ln = 0 , n_hjorth = 0;
      int max_consec = 0 , cur_consec = 0;
      int max_ln_consec = 0 , cur_ln_consec = 0;

      for ( int e = 0; e < total_epochs; e++ )
	{
	  if ( flag_flat[e]   ) ++n_flat;
	  if ( flag_clip[e]   ) ++n_clip;
	  if ( flag_amp[e]    ) ++n_amp;
	  if ( flag_hf[e]     ) ++n_hf;
	  if ( flag_hjorth[e] ) ++n_hjorth;
	  if ( flag_ln[e] )
	    {
	      ++n_ln;
	      if ( ++cur_ln_consec > max_ln_consec ) max_ln_consec = cur_ln_consec;
	    }
	  else cur_ln_consec = 0;
	  if ( flag_bad[e] )
	    {
	      ++n_bad;
	      if ( ++cur_consec > max_consec ) max_consec = cur_consec;
	    }
	  else cur_consec = 0;
	}

      const double prop_bad    = total_epochs > 0 ? n_bad / (double)total_epochs : 0.0;
      const double max_bad_run = max_consec * eog_window_inc;
      const bool   ch_bad      = ( prop_bad > eog_flag_prop ) ||
	                          ( max_bad_run >= eog_flag_run );

      const double prop_ln     = total_epochs > 0 ? n_ln / (double)total_epochs : 0.0;
      const double max_ln_run  = max_ln_consec * eog_window_inc;
      const bool   ch_ln_bad   = ( prop_ln > eog_flag_prop ) ||
	                          ( max_ln_run >= eog_flag_run );

      writer.value( "N_FLAG_EPOCH" , n_bad      );
      writer.value( "PROP_FLAG"    , prop_bad   );
      writer.value( "MAX_FLAG_RUN" , max_bad_run );
      writer.value( "FLAGGED"         , (int)ch_bad );
      writer.value( "P_FLAT"      , total_epochs > 0 ? n_flat   / (double)total_epochs : 0.0 );
      writer.value( "P_CLIP"      , total_epochs > 0 ? n_clip   / (double)total_epochs : 0.0 );
      writer.value( "P_AMP"       , total_epochs > 0 ? n_amp    / (double)total_epochs : 0.0 );
      writer.value( "P_HF"        , total_epochs > 0 ? n_hf     / (double)total_epochs : 0.0 );
      writer.value( "P_HJORTH"    , total_epochs > 0 ? n_hjorth / (double)total_epochs : 0.0 );
      writer.value( "P_LN"        , prop_ln    );
      writer.value( "MAX_LN_RUN"  , max_ln_run );
      writer.value( "LN_FLAG"      , (int)ch_ln_bad );

      writer.unlevel( "EOG" );

      //
      // Cross-domain summary (CH x DOMAIN table)
      //

      writer.level( "EOG" , "DOMAIN" );
      writer.value( "FLAGGED"         , (int)ch_bad    );
      writer.value( "N_FLAG_EPOCH" , n_bad          );
      writer.value( "MAX_FLAG_RUN" , max_bad_run    );
      writer.value( "LN_FLAG"      , (int)ch_ln_bad );
      writer.value( "MAX_LN_RUN"  , max_ln_run     );
      writer.unlevel( "DOMAIN" );

      //
      // Annotations (merged, gap-aware)
      //

      if ( add_annots )
	{
	  std::string alab = annot_prefix;
	  if ( annot_show_domain ) alab += "_EOG";
	  if ( annot_by_channel  ) alab += "_" + signals.label(s);

	  annot_t * a = edf.annotations->add( alab );

	  bool in_bad = false;
	  uint64_t bad_start = 0LLU , bad_stop = 0LLU;

	  for ( int e = 0; e < total_epochs; e++ )
	    {
	      if ( flag_bad[e] )
		{
		  if ( ! in_bad )
		    {
		      bad_start = epoch_intervals[e].start;
		      bad_stop  = epoch_intervals[e].stop;
		      in_bad = true;
		    }
		  else
		    {
		      if ( epoch_intervals[e].start > bad_stop )
			{
			  a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );
			  bad_start = epoch_intervals[e].start;
			}
		      if ( epoch_intervals[e].stop > bad_stop )
			bad_stop = epoch_intervals[e].stop;
		    }
		}
	      else if ( in_bad )
		{
		  a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );
		  in_bad = false;
		}
	    }

	  if ( in_bad )
	    a->add( "." , interval_t( bad_start , bad_stop ) , signals.label(s) );

	  // separate LN annotations
	  if ( n_ln > 0 )
	    {
	      std::string ln_lab = annot_prefix;
	      if ( annot_show_domain ) ln_lab += "_EOG";
	      ln_lab += "_LN";
	      if ( annot_by_channel ) ln_lab += "_" + signals.label(s);

	      annot_t * a_ln = edf.annotations->add( ln_lab );
	      bool in_ln = false;
	      uint64_t ln_start = 0LLU , ln_stop = 0LLU;

	      for ( int e = 0; e < total_epochs; e++ )
		{
		  if ( flag_ln[e] )
		    {
		      if ( ! in_ln ) { ln_start = epoch_intervals[e].start; ln_stop = epoch_intervals[e].stop; in_ln = true; }
		      else { if ( epoch_intervals[e].start > ln_stop ) { a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) ); ln_start = epoch_intervals[e].start; } if ( epoch_intervals[e].stop > ln_stop ) ln_stop = epoch_intervals[e].stop; }
		    }
		  else if ( in_ln ) { a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) ); in_ln = false; }
		}
	      if ( in_ln ) a_ln->add( "." , interval_t( ln_start , ln_stop ) , signals.label(s) );
	    }
	}

      logger << "  " << signals.label(s)
	     << " : " << n_bad << "/" << total_epochs << " flagged epochs"
	     << "  (flat=" << n_flat << " clip=" << n_clip << " amp=" << n_amp
	     << " hf=" << n_hf << " hjorth=" << n_hjorth << ")"
	     << " | ln=" << n_ln
	     << ( ch_ln_bad ? " --> LN_FLAG" : "" )
	     << ( ch_bad    ? " --> FLAGGED channel" : "" ) << "\n";

    } // next channel

  writer.unlevel( globals::signal_strat );

}



//
// ECG QC
//

/*
  a. Signal-to-Noise Ratio (SNR)
    Estimate using power ratio:
    SNR=10⋅log⁡10(PsignalPnoise)
    SNR=10⋅log10​(Pnoise​Psignal​​)

    Often approximated by comparing total power to high-frequency noise (>40–50 Hz) or baseline wander power (<0.5 Hz)


  b. Flatline detection
    Long segments of zero or nearly constant values
    Threshold: e.g., std(signal[t:t+2s]) < 0.01 mV

  c. Amplitude clipping
    Clipped values (at ADC limit): e.g., repeated max/min values
    May indicate hardware or digitization problems


  a. RR interval plausibility
    Valid range: 300–2000 ms (30–200 bpm)
    Flag RR intervals outside this range
  

   🧪 In Practice (Example Thresholds)
QC Metric	Good Range	Flags
RR interval range	300–2000 ms	<300 or >2000 ms
Flatline duration	<1 sec	>1 sec constant
QRS width	70–120 ms	<60 or >150 ms
Line noise power	Low (<5% total)	High (>10%)
Power <0.5 Hz	Low	>20% total power

*/
