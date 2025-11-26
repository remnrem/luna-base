
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

#include "stats/catch22/catch22.h"

#include "helper/logger.h"
#include "helper/helper.h"

extern logger_t logger;

#include <vector>
#include <map>

// a simple wrapper around 'CAnonical Time-series CHaracteristics' 
//   Lubba et al. catch22: CAnonical Time-series CHaracteristics, Data Min Knowl Disc 33, 1821 (2019).
//   https://github.com/DynamicsAndNeuralSystems/catch22?tab=readme-ov-file
//   https://time-series-features.gitbook.io/catch22/
//   version 0.4.0

std::vector<std::string> catch22_t::names = { 
  "DN_HistogramMode_5", 
  "DN_HistogramMode_10", 
  "CO_f1ecac", 
  "CO_FirstMin_ac", 
  "CO_HistogramAMI_even_2_5", 
  "CO_trev_1_num", 
  "MD_hrv_classic_pnn40", 
  "SB_BinaryStats_mean_longstretch1", 
  "SB_TransitionMatrix_3ac_sumdiagcov", 
  "PD_PeriodicityWang_th0_01", 
  "CO_Embed2_Dist_tau_d_expfit_meandiff", 
  "IN_AutoMutualInfoStats_40_gaussian_fmmi", 
  "FC_LocalSimple_mean1_tauresrat", 
  "DN_OutlierInclude_p_001_mdrmd", 
  "DN_OutlierInclude_n_001_mdrmd", 
  "SP_Summaries_welch_rect_area_5_1", 
  "SB_BinaryStats_diff_longstretch0", 
  "SB_MotifThree_quantile_hh", 
  "SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1", 
  "SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1", 
  "SP_Summaries_welch_rect_centroid", 
  "FC_LocalSimple_mean3_stderr", 
  "DN_Mean", 
  "DN_Spread_Std" };

std::vector<std::string> catch22_t::short_names = { 
  "mode_5",
  "mode_10",
  "acf_timescale",
  "acf_first_min",
  "ami2",
  "trev",
  "high_fluctuation",
  "stretch_high",
  "transition_matrix",
  "periodicity",
  "embedding_dist",
  "ami_timescale",
  "whiten_timescale",
  "outlier_timing_pos",
  "outlier_timing_neg",
  "centroid_freq",
  "stretch_decreasing",
  "entropy_pairs",
  "rs_range",
  "dfa",
  "low_freq_power",
  "forecast_error",
  "mean",
  "SD" } ;


catch22_t::catch22_t( bool catch24 )
{

  qc = -1;
  
  nstats = catch24 ? 24 : 22;
   
}

bool catch22_t::calc( const double y[], int size )
{
   
  qc = quality_check( y , size );

  // 0 good
  // 1 too short
  // 2 has Inf
  // 3 has NaN

  if ( qc == 0 )
    run_features( y , size , nstats == 24 );

  return qc == 0 ;
  
}

double catch22_t::stat( const int i ) const
{
  if ( i >= 0 && i < nstats )
    return res.find( names[i] )->second ;
  else
    return std::numeric_limits<double>::quiet_NaN() ;
}

std::string catch22_t::name( const int i ) 
{
  if ( i >= 0 && i < 24 ) 
    return names[i];
  else
    return ".";
}

std::string catch22_t::short_name( const int i )
{
  if ( i >= 0 && i < 24 ) 
    return short_names[i];
  else
    return ".";
}


// check if data qualifies to be caught22

int catch22_t::quality_check(const double y[], const int size)
{
  
  int minSize = 10;
  
  if(size < minSize)
    {
      return 1;
    }
  for(int i = 0; i < size; i++)
    {
      double val = y[i];
      if(val == INFINITY || -val == INFINITY)
        {
	  return 2;
        }
      if(isnan(val))
        {
	  return 3;
        }
    }
  return 0;
}



void catch22_t::run_features(const double y[], int size , bool catch24 )
{
  
  res.clear();
  
  // normalize
  std::vector<double> y_zscored( size ); 
  zscore_norm2(y, size, y_zscored.data());

  // GOOD
  res[ "DN_OutlierInclude_n_001_mdrmd" ] = DN_OutlierInclude_n_001_mdrmd(y_zscored.data(), size);
  
  // GOOD  
  res[ "DN_OutlierInclude_p_001_mdrmd" ] = DN_OutlierInclude_p_001_mdrmd(y_zscored.data(), size);
  
  // GOOD  
  res[ "DN_HistogramMode_5" ] = DN_HistogramMode_5(y_zscored.data(), size);
  
  // GOOD
  res[ "DN_HistogramMode_10" ] = DN_HistogramMode_10(y_zscored.data(), size);
  
  //GOOD
  res[ "CO_Embed2_Dist_tau_d_expfit_meandiff" ] = CO_Embed2_Dist_tau_d_expfit_meandiff(y_zscored.data(), size);
  
  //GOOD (memory leak?)
  res[ "CO_f1ecac" ] = CO_f1ecac(y_zscored.data(), size);
  
  //GOOD
  res[ "CO_FirstMin_ac" ] = CO_FirstMin_ac(y_zscored.data(), size);
  
  // GOOD (memory leak?)
  res[ "CO_HistogramAMI_even_2_5" ] = CO_HistogramAMI_even_2_5(y_zscored.data(), size);
  
  // GOOD
  res[ "CO_trev_1_num" ] = CO_trev_1_num(y_zscored.data(), size);
  
  //GOOD
  res[ "FC_LocalSimple_mean1_tauresrat" ] = FC_LocalSimple_mean1_tauresrat(y_zscored.data(), size);
  
  //GOOD
  res[ "FC_LocalSimple_mean3_stderr" ] = FC_LocalSimple_mean3_stderr(y_zscored.data(), size);
  
  //GOOD (memory leak?)
  res[ "IN_AutoMutualInfoStats_40_gaussian_fmmi" ] = IN_AutoMutualInfoStats_40_gaussian_fmmi(y_zscored.data(), size);
  
  //GOOD
  res[ "MD_hrv_classic_pnn40" ] = MD_hrv_classic_pnn40(y_zscored.data(), size);
  
  //GOOD
  res[ "SB_BinaryStats_diff_longstretch0" ] = SB_BinaryStats_diff_longstretch0(y_zscored.data(), size);

  //GOOD
  res[ "SB_BinaryStats_mean_longstretch1" ] = SB_BinaryStats_mean_longstretch1(y_zscored.data(), size);
  
  //GOOD (memory leak?)
  res[ "SB_MotifThree_quantile_hh" ] = SB_MotifThree_quantile_hh(y_zscored.data(), size);
  
  //GOOD (memory leak?)
  res[ "SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1" ] = SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1(y_zscored.data(), size);
  
  //GOOD
  res[ "SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1" ] = SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1(y_zscored.data(), size);

  //GOOD
  res[ "SP_Summaries_welch_rect_area_5_1" ] = SP_Summaries_welch_rect_area_5_1(y_zscored.data(), size);
  
  //GOOD
  res[ "SP_Summaries_welch_rect_centroid" ] = SP_Summaries_welch_rect_centroid(y_zscored.data(), size);

  //OK, BUT filt in Butterworth sometimes diverges, now removed alltogether, let's see results.
  res[ "SB_TransitionMatrix_3ac_sumdiagcov" ] = SB_TransitionMatrix_3ac_sumdiagcov(y_zscored.data(), size);

  // GOOD
  res[ "PD_PeriodicityWang_th0_01" ] = PD_PeriodicityWang_th0_01(y_zscored.data(), size);
  
  if (catch24) {
    
    // GOOD
    res[ "DN_Mean" ] = DN_Mean(y, size);    
    // GOOD
    res[ "DN_Spread_Std" ] = DN_Spread_Std(y, size);
  }

      
}
  

