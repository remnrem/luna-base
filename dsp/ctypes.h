
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

#ifdef HAS_LGBM

#ifndef __CTYPES_H__
#define __CTYPES_H__

struct edf_t;
struct param_t;
struct signal_list_t;

#include "lgbm/lgbm.h"
#include "stats/Eigen/Dense"

enum ctypes_ch_t
  {
    // neural
    CTYPE_EEG = 1 ,           // any EEG, any referece type, includes mastoids
    
    // ocular
    CTYPE_EOG = 2 , 
    
    // neuromuscular
    CTYPE_EMG = 3 , 

    // cardiovascular
    CTYPE_ECG = 4 ,           // EKG
    CTYPE_PLETH = 5 ,         // PPG waveform (oximetry)
    CTYPE_SPO2 = 6,           // PPG-derived 
    CTYPE_HR = 7,             // pulse rate / HR 
    
    // respiratory
    CTYPE_AIRFLOW = 8,        // any flow channel (thermistor, nasal canula)
    CTYPE_RESP = 9 ,         // any resp effort channel
    CTYPE_CO2ET = 10 ,       // capnography / peak value per breath (EtCO2)
    CTYPE_CO2WAVEFORM = 11 , // capnography / waveform
    CTYPE_SNORE = 12 ,        // audio/mic
    
    // Other common
    CTYPE_POSITION = 13 ,     // any type (e.g. 2, 5 or more classes) 
    
    // Auxilliary 
    CTYPE_MARKER = 14 ,       // generic light / event marker / status channel
    CTYPE_FLAT = 15 ,         // effectively empty channels    
    CTYPE_ARTIFACT = 0 ,      // generic high-artifact channel (non-physiologic signal)


    // special models; may be included, or may be from special models

    // generic EMG -> chin vs leg EMG
    CTYPE_CHIN = 16 ,      
    CTYPE_LEG = 17,        

    // generic AIRFLOW -> THERMISTOR vs NASAL CANNULA
    CTYPE_THERM = 18 , 
    CTYPE_NPRES = 19 

  };


struct ctypes_pred_t {
  bool valid;
  std::map<ctypes_ch_t,double> posteriors;  
};

struct ctypes_specific_ftrs_t {
  
  double flatline_frac;
  double clip_frac; 
  double unique_frac;  
  double most_common_value_frac;
  double transition_rate;
  
  double logH1;
  double H2;
  double H3;
  double line_length;
  double kurtosis;
  double skewness;  
  double ZCR;
  
  double spectral_centroid;
  double spectral_edge; // freq at 90%
  double spectral_bandwidth;
  double spectral_flatness;
  double spectral_entropy;
  double spectral_lowpower;
  double spectral_highpower;
  
  double acf1;
  double acf_1s;
  double acf_decay;
  double acf_peak; // maximun r[lag] between 0.5–5 seconds
  double acf_min; // minimum r[lag] between 0.1–1.0 seconds
  
};

struct ctypes_ftrs_t {

  // overall features 
  double n_epochs;
  double pct_rem;
  double pct_n3;

  // 1 Hz median 
  ctypes_specific_ftrs_t x1;         
  ctypes_specific_ftrs_t x1_p10;
  ctypes_specific_ftrs_t x1_p90;

  // 1 Hz first-diff
  ctypes_specific_ftrs_t x1_diff;    
  ctypes_specific_ftrs_t x1_diff_p10;    
  ctypes_specific_ftrs_t x1_diff_p90;    

  // 1 Hz envelope 
  ctypes_specific_ftrs_t a1;
  ctypes_specific_ftrs_t a1_p10;
  ctypes_specific_ftrs_t a1_p90;         
  
  // 128 Hz 
  ctypes_specific_ftrs_t x128;
  ctypes_specific_ftrs_t x128_p10;
  ctypes_specific_ftrs_t x128_p90;

  // 128 hz first diff
  ctypes_specific_ftrs_t x128_diff;  
  ctypes_specific_ftrs_t x128_diff_p10;
  ctypes_specific_ftrs_t x128_diff_p90; 

};




struct spectral_stats_t
{
  spectral_stats_t( const std::vector<double> & x , double sr )
  {
    compute( x , sr );
  }

  void compute(const std::vector<double>& x, double sr );    
  
  double centroid;
  double edge90;
  double bandwidth;
  double flatness;
  double entropy;
  double high;
  double low;
 
};


struct ctypes_t {
  
  ctypes_t( edf_t & , param_t & );

  void proc( edf_t & , param_t & );
  
  void calc_1Hz_stats( const std::vector<double> & x , double,
			 const std::vector<int> & , 
			 const std::vector<std::pair<int,int>> & ,
			 size_t , 
			 ctypes_ftrs_t * , bool );

  void calc_128Hz_stats( const std::vector<double> & x , double,
			 const std::vector<int> & , 
			 const std::vector<std::pair<int,int>> & ,
			 size_t , 
			 ctypes_ftrs_t * , bool );
  
  ctypes_specific_ftrs_t calc_specific_stats( const std::vector<double> & x , int Fs );

  void aggregate1( const std::vector<ctypes_ftrs_t> & aggr, ctypes_ftrs_t * ftr );

  void aggregate128( const std::vector<ctypes_ftrs_t> & aggr, ctypes_ftrs_t * ftr );
		    

  
private:
  
    
  //
  // feature map
  //

  using member_ptr = double ctypes_specific_ftrs_t::*;

  const std::map<std::string, member_ptr> feature_map = {
    { "flat", &ctypes_specific_ftrs_t::flatline_frac},
    { "clip", &ctypes_specific_ftrs_t::clip_frac}, 
    { "ufrac", &ctypes_specific_ftrs_t::unique_frac},  
    { "cfrac", &ctypes_specific_ftrs_t::most_common_value_frac},
    { "trans", &ctypes_specific_ftrs_t::transition_rate},    
    { "h1", &ctypes_specific_ftrs_t::logH1},
    { "h2", &ctypes_specific_ftrs_t::H2},
    { "h3", &ctypes_specific_ftrs_t::H3},
    { "line_length", &ctypes_specific_ftrs_t::line_length},
    { "kurtosis", &ctypes_specific_ftrs_t::kurtosis},
    { "skew", &ctypes_specific_ftrs_t::skewness},  
    { "zcr", &ctypes_specific_ftrs_t::ZCR},
    { "spec_centroid", &ctypes_specific_ftrs_t::spectral_centroid},
    { "spec_edge", &ctypes_specific_ftrs_t::spectral_edge},
    { "spec_bandwidth", &ctypes_specific_ftrs_t::spectral_bandwidth},
    { "spec_flatness", &ctypes_specific_ftrs_t::spectral_flatness},
    { "spec_entropy", &ctypes_specific_ftrs_t::spectral_entropy},
    { "spec_low", &ctypes_specific_ftrs_t::spectral_lowpower},
    { "spec_high", &ctypes_specific_ftrs_t::spectral_highpower},
    { "acf1", &ctypes_specific_ftrs_t::acf1},
    { "acf_1s", &ctypes_specific_ftrs_t::acf_1s},
    { "acf_delay", &ctypes_specific_ftrs_t::acf_decay},
    { "acf_peak", &ctypes_specific_ftrs_t::acf_peak}, 
    { "acf_min", &ctypes_specific_ftrs_t::acf_min }    
  };
  
  
  static double median_inplace(std::vector<double>& v);
  
  static inline bool is_finite(double v) { return std::isfinite(v); } 
  
  
  static inline bool is_nan(double x);

  static double percentile_inplace(std::vector<double> & v, double p);

  static inline void compute_p10_med_p90(std::vector<double> & vals,
					 double & out_p10,
					 double & out_med,
					 double & out_p90);
  
  // threshold for doing per-epoch stats (128Hz signal)
  double edge95_th;
  double edge95_prop;
  double edge95_mean;

  // threshold for doing whole-night stats (1Hz signal)
  double edge95_mean2;
  double edge95_th2;
  double edge95_prop2;

  // transformations

  void make_1s( const std::vector<double>& x, double sr ,
		std::vector<double> * x1, std::vector<double> * a1 );

  std::vector<double> make_diff(const std::vector<double>& x);

  // std::vector<double> make_corr(const std::vector<double>& x,
  // 				double Fs,               
  // 				double win_sec );

  // std::vector<double> make_corr_fast(const std::vector<double>& x,
  // 				     double Fs,
  // 				     double win_sec,
  // 				     int decim);
  
  
  void normalize( std::vector<double> * );
  
  // stats

  double mean(const std::vector<double>& x);

  double sd_sample(const std::vector<double>& x);

  double mean_ptr(const double* x, size_t n);

  double sd_sample_ptr(const double* x, size_t n);

  double quantile(std::vector<double> xs, double p);

  double mad(const std::vector<double>& x);

  int sgn(double v, double eps = 0.0);

  double line_length(const std::vector<double>& x);

  double line_length_norm_sd(const std::vector<double>& x);

  double zero_cross_rate(const std::vector<double>& x, double eps = 0.0);

  double flatline_fraction(const std::vector<double>& x);

  double clip_fraction(const std::vector<double>& x);

  long long qbin(double v, double eps = 1e-2 );

  double unique_frac_q(const std::vector<double>& x, double eps = 1e-2 );

  double most_common_value_frac_q(const std::vector<double>& x, double eps = 1e-2 );
  
  double transition_rate_q(const std::vector<double>& x, double eps = 1e-2 );


  
  // LGBM
  
  static lgbm_t lgbm;

  static std::string lgbm_model_loaded;

  static ctypes_pred_t predict( const ctypes_ftrs_t & );
			       
};

#endif

#endif
