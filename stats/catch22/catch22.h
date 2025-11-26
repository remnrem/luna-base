
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

#ifndef __CATCH22_H__
#define __CATCH22_H__

#include "stats/Eigen/Dense"
#include <vector>
#include <map>

// a simple wrapper around 'CAnonical Time-series CHaracteristics' 
//   Lubba et al. catch22: CAnonical Time-series CHaracteristics, Data Min Knowl Disc 33, 1821 (2019).
//   https://github.com/DynamicsAndNeuralSystems/catch22?tab=readme-ov-file
//   https://time-series-features.gitbook.io/catch22/
//   version 0.4.0

#include "stats/catch22/DN_HistogramMode_5.h"
#include "stats/catch22/DN_HistogramMode_10.h"
#include "stats/catch22/DN_Mean.h"
#include "stats/catch22/DN_Spread_Std.h"
#include "stats/catch22/CO_AutoCorr.h"
#include "stats/catch22/DN_OutlierInclude.h"
#include "stats/catch22/FC_LocalSimple.h"
#include "stats/catch22/IN_AutoMutualInfoStats.h"
#include "stats/catch22/MD_hrv.h"
#include "stats/catch22/SB_BinaryStats.h"
#include "stats/catch22/SB_MotifThree.h"
#include "stats/catch22/SC_FluctAnal.h"
#include "stats/catch22/SP_Summaries.h"
#include "stats/catch22/SB_TransitionMatrix.h"
#include "stats/catch22/PD_PeriodicityWang.h"
#include "stats/catch22/c22_stats.h"

struct catch22_t {
  
  catch22_t( bool catch24 = false );
  bool calc( const double y[], int size );
  int quality_check(const double y[], const int size);
  void run_features(const double y[], int size , bool catch24 );

  // result cache
  std::map<std::string,double> res;

  // qc return code (non-zero = fail)
  int qc;

  bool valid() const { return qc == 0 && res.size() >= 22 ; }
  
  // names
  static std::vector<std::string> names, short_names;
  
  double stat( const int i ) const ;
  static std::string name( const int i );
  static std::string short_name( const int i );

  // 22 or 24 stats
  int nstats;
  
};


#endif 
