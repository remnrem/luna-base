
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

#ifndef __HB_H__
#define __HB_H__

#include "stats/Eigen/Dense"
#include "intervals/intervals.h"

#include <vector>

struct edf_t;
struct param_t;

// find_burden() return value
struct hb_find_burden_t {
  bool valid;
  int ne;
  double HB;
  double BaselineSat;
  Eigen::MatrixXd SpO2MtxDiff;
  Eigen::ArrayXd BaselineSatAll;
  int searchWin_lwr, searchWin_upr;
};

// peakdet() return value
struct hb_peakdet_t {
  std::vector<double> maxV; 
  std::vector<double> maxX;
  std::vector<double> minV;
  std::vector<double> minX;
};


struct hb_event_results_t {
  interval_t interval; 
  std::string type;
  bool arousal;
  sleep_stage_t stage;
  double meanSatPerEvent;
  double minSatPerEvent;
  double SatSrchWinPre;
  double SatSrchWinPost;
  double HBperEvent;
  double BaselineSatperEvent;
  double DesatDur;
  double ResatDur;
  double DesatStartMag;
  double DesatNadirMag;
  double DesatEndMag;
  double DesatStartTime;
  double DesatMag;
  double dHR_meanBsln_PU;
  double dHR_minBsln_PU;
  double meanBslnHR_PU;
  double minBslnHR_PU;
 
};


struct hb_results_t {

  hb_results_t() {
    valid = false;
    
    TotalAHI = 0.0;
    srchWinStart = 0;
    srchWinEnd = 0;
    
    AHI_N1 = AHI_N2 = AHI_N3 = AHI_REM = AHI_NREM = 0;
    HB =  HB_N1 = HB_N2 = HB_N3 = HB_REM = HB_NREM = 0;
    BaselineSat = BaselineSat_N1 = BaselineSat_N2 = BaselineSat_N3 = BaselineSat_REM = BaselineSat_NREM = 0;
  
    HB4 = HB3pa = HBtot = HBtot3 = 0;
    NREM_HB4 = NREM_HB3pa = NREM_HBtot = NREM_HBtot3 = 0;
    REM_HB4 = REM_HB3pa = REM_HBtot = REM_HBtot3 = 0;

    HB3alt = NREM_HB3alt = REM_HB3alt = 0;
    NREM_BaselineSat = 0; REM_BaselineSat = 0;

    TST = TST_NREM = TST_REM = TST_N1 = TST_N2 = TST_N3 = TIB = 0;
    dHR_meanbsline_PU = dHR_minbsline_PU = NdHR_PU = 0;

  }

  bool valid;

  // event-level
  std::vector<hb_event_results_t> events;

  // individual-level
  double TotalAHI;
  int srchWinStart;
  int srchWinEnd;

  double AHI_N1, AHI_N2, AHI_N3, AHI_REM, AHI_NREM;
  double HB, HB_N1, HB_N2, HB_N3, HB_REM, HB_NREM;  
  double BaselineSat, BaselineSat_N1, BaselineSat_N2, BaselineSat_N3, BaselineSat_REM, BaselineSat_NREM;

  double HB4, HB3pa, HBtot, HBtot3;
  double NREM_HB4, NREM_HB3pa, NREM_HBtot, NREM_HBtot3;
  double REM_HB4, REM_HB3pa, REM_HBtot, REM_HBtot3;

  double HB3alt, NREM_HB3alt, REM_HB3alt;
  double NREM_BaselineSat, REM_BaselineSat;

  double TST, TST_NREM, TST_REM, TST_N1, TST_N2, TST_N3, TIB;
  double dHR_meanbsline_PU, dHR_minbsline_PU, NdHR_PU;
  
};

struct hb_find_desats_t
{
  //%% Fins all desaturations with avegrage desat/resat of at least 1.5%
  Eigen::ArrayXd  MagDown;
  Eigen::ArrayXd  MagUp;     
  Eigen::ArrayXXi dsatStEnd; // start, peak, end (sample-points)
};

struct delta_hr_t {
  double dHR_meanbsline;
  double dHR_minbsline;
  Eigen::ArrayXd IndHRRep_mean;
  Eigen::ArrayXd IndHRRep_min;
  int ne;
  Eigen::ArrayXd meanbsline;
  Eigen::ArrayXd minbsline;
  std::vector<bool> nan;
};

struct hb_t {
  
  hb_t( edf_t & edf , param_t & );
  
  static hb_peakdet_t peakdet( const Eigen::ArrayXd & v , 
			       double delta ,
			       const bool );

  static hb_peakdet_t peakdet( const Eigen::ArrayXd & v , 
			       double delta ,
			       const std::vector<double> & x ,
			       const bool );
  
  static hb_find_burden_t find_burden( const Eigen::MatrixXd & SpO2Mtx ,
				       const Eigen::ArrayXd & SpO2Mean ,
				       const std::vector<double> Time ,
				       const double TST , 
				       const int MaxWin ,
				       const std::vector<bool> * incl = NULL );
  

  static sleep_stage_t modal_stage( const Eigen::ArrayXi & );
  
  static std::vector<bool> which_events( const std::vector<sleep_stage_t> & , const std::string & s ,
					 const std::vector<bool> * orig = NULL );

  static bool enough( const std::vector<bool> & , int th = 1 );
  
  static hb_find_desats_t find_desats( const Eigen::ArrayXd & , int , double );

  static double min( const std::vector<double> & s , int * idx );

  static delta_hr_t SummarizeHR_AA( const Eigen::MatrixXd & HRVe ,
				    const std::vector<double> & Start ,
				    const std::vector<double> & End , 
				    const std::vector<double> & Time_new );
  
};




#endif
