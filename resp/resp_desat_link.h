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

#ifndef __RESP_DESAT_LINK_H__
#define __RESP_DESAT_LINK_H__

#include <string>
#include <vector>

struct edf_t;
struct param_t;

struct resp_desat_link_event_t {
  double resp_start_sec;
  double resp_stop_sec;
  double resp_dur_sec;
  std::string resp_type;
  bool in_sleep;
  bool linked;
  std::string reject_reason;
  bool arousal_linked;
  std::string arousal_reject_reason;
  std::string arousal_type;
  double arousal_start_sec;
  double arousal_stop_sec;
  double arousal_dur_sec;
  double arousal_onset_rel_start_sec;
  double arousal_onset_rel_stop_sec;
  double baseline_center_t;
  double nadir_center_t;
  double recovery_center_t;
  double prev_rel_stop_sec;
  double next_rel_start_sec;
  double base_a_sec;
  double base_b_sec;
  double nad_a_sec;
  double nad_b_sec;
  double baseline_sec;
  double baseline_val;
  double nadir_sec;
  double nadir_val;
  double recovery_sec;
  double recovery_val;
  double drop;
};

struct resp_desat_link_t {
  resp_desat_link_t( edf_t & edf , param_t & param );
};

#endif
