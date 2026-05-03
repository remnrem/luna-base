#ifndef __LUNAPI_WAVEFORMS_H__
#define __LUNAPI_WAVEFORMS_H__

#include <limits>
#include <map>
#include <string>
#include <vector>

struct edf_t;

struct waveform_feature_options_t {
  bool catch_enabled = true;
  bool catch24 = false;
  bool basic_stats = true;
};

struct waveform_block_t {
  std::string label;
  std::string unit;
  int sr = 0;
  double data_start_sec = std::numeric_limits<double>::quiet_NaN();
  double data_stop_sec = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> rel_time;
  std::vector<double> values;
  int feature_qc = -1;
  std::map<std::string,double> features;
};

struct waveform_event_t {
  std::string annot;
  std::string instance;
  std::string annot_ch;
  std::string meta;
  double annot_start_sec = std::numeric_limits<double>::quiet_NaN();
  double annot_stop_sec = std::numeric_limits<double>::quiet_NaN();
  double anchor_sec = std::numeric_limits<double>::quiet_NaN();
  double wave_start_sec = std::numeric_limits<double>::quiet_NaN();
  double wave_stop_sec = std::numeric_limits<double>::quiet_NaN();
  std::vector<waveform_block_t> blocks;
};

struct waveform_extract_result_t {
  std::vector<std::string> requested_annots;
  std::vector<std::string> requested_channels;
  std::vector<std::string> feature_names;
  std::vector<waveform_event_t> events;
  std::map<std::string,int> dropped;
  std::string align = "mid";
  std::string mode = "fixed";
  double pre_secs = 0.0;
  double post_secs = 0.0;
  double flank_left_secs = 0.0;
  double flank_right_secs = 0.0;
  int total_events = 0;
};

namespace waveform_core
{
  waveform_extract_result_t extract_fixed_window_waveforms(
    edf_t & edf ,
    const std::vector<std::string> & annots ,
    const std::vector<std::string> & channels ,
    const double pre_secs ,
    const double post_secs ,
    const std::string & align = "mid" );

  waveform_extract_result_t extract_annotation_window_waveforms(
    edf_t & edf ,
    const std::vector<std::string> & annots ,
    const std::vector<std::string> & channels ,
    const double flank_left_secs ,
    const double flank_right_secs ,
    const std::string & align = "mid" ,
    const std::string & require = "full" );

  std::vector<std::string> feature_names(
    const waveform_feature_options_t & options = waveform_feature_options_t() );

  void compute_features(
    waveform_extract_result_t * result ,
    const waveform_feature_options_t & options = waveform_feature_options_t() );
}

#endif
