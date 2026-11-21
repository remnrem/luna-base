
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

#ifndef LUNA_PSD_APERIODIC_H
#define LUNA_PSD_APERIODIC_H

#include <map>
#include <string>
#include <vector>
#include "spectral/aperiodic-fit.h"

struct param_t;
struct edf_t;

// Optional PSD sidecar. Does not modify the PSD, its bands, or existing tracks.
class psd_aperiodic_t {
public:
  psd_aperiodic_t( param_t &, double min, double max );
  void epoch( const std::vector<double> & f, const std::vector<double> & p, bool bands,
              const std::vector<double> & data, const std::string & channel,
              int epoch, double start, double stop );
  void summary( bool bands );
  void add( edf_t &, const std::string & label, int observed, int expected, int pad1, int pad2 );
private:
  bool enabled, db, silent, epoch_output, spectrum_output, epoch_spectrum,
    adding, adding_spectrum, any_variable = false;
  double lwr, upr, threshold, bin_fac;
  double pivot = 1, min_db = 1;
  std::vector<double> regions, evaluation;
  int iterations, count = 0, failed_count = 0, fallback_count = 0,
    flat_count = 0, gap = 0, max_gap = 0;
  std::string prefix, channel_label;
  std::vector<double> frequencies, sum, valid_windows, fallback_windows;
  std::map<std::string,std::vector<double> > tracks;
  std::vector<std::string> skip;
  void report( const std::vector<double> &, const std::vector<double> &,
               const psd_aperiodic::fit_t &, bool bands, bool spectrum );
  void values( double original, double background );
  bool band_in_range( int band ) const;
};
#endif
